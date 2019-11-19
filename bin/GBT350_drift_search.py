#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from operator import attrgetter
import glob, os, os.path, shutil, socket, tarfile, stat
import numpy, sys, time
from presto import sigproc
from presto import sifting
from presto import presto
from presto import psr_utils as pu

institution = "NRAOCV" 
base_tmp_dir = "/dev/shm/"
base_output_dir = "/home/sransom/results/GBT/drift"

#-------------------------------------------------------------------
# Tunable parameters for searching and folding
# (you probably don't need to tune any of them)
orig_N                = 1728000 # Number of samples to analyze at a time (~141 sec)
raw_N                 = 1900000 # Number of samples to step through .fits files
overlap_factor        = 0.5  # Overlap each orig_N samples by this fraction
rfifind_chunk_time    = 25600 * 0.00008192  # ~2.1 sec
singlepulse_threshold = 5.0  # threshold SNR for candidate determination
singlepulse_plot_SNR  = 5.5  # threshold SNR for singlepulse plot
singlepulse_maxwidth  = 0.1  # max pulse width in seconds
to_prepfold_sigma     = 6.0  # incoherent sum significance to fold candidates
max_lo_cands_to_fold  = 20   # Never fold more than this many lo-accel candidates
max_hi_cands_to_fold  = 10   # Never fold more than this many hi-accel candidates
numhits_to_fold       = 2    # Number of DMs with a detection needed to fold
low_DM_cutoff         = 1.0  # Lowest DM to consider as a "real" pulsar
lo_accel_numharm      = 16   # max harmonics
lo_accel_sigma        = 2.0  # threshold gaussian significance
lo_accel_zmax         = 0    # bins
lo_accel_flo          = 2.0  # Hz
hi_accel_numharm      = 8    # max harmonics
hi_accel_sigma        = 3.0  # threshold gaussian significance
hi_accel_zmax         = 50   # bins
hi_accel_flo          = 1.0  # Hz
low_T_to_search       = 50.0 # sec
# Sifting specific parameters (don't touch without good reason!)
sifting.sigma_threshold = to_prepfold_sigma-1.0  # incoherent power threshold (sigma)
sifting.c_pow_threshold = 100.0                  # coherent power threshold
sifting.r_err           = 1.1    # Fourier bin tolerence for candidate equivalence
sifting.short_period    = 0.0005 # Shortest period candidates to consider (s)
sifting.long_period     = 15.0   # Longest period candidates to consider (s)
sifting.harm_pow_cutoff = 8.0    # Power required in at least one harmonic
#-------------------------------------------------------------------

def find_masked_fraction(obs):
    """
    find_masked_fraction(obs):
        Parse the output file from an rfifind run and return the
        fraction of the data that was suggested to be masked.
    """
    rfifind_out = obs.basefilenm + "_rfifind.out"
    for line in open(rfifind_out):
         if "Number of  bad   intervals" in line:
              return float(line.split("(")[1].split("%")[0])/100.0
    # If there is a problem reading the file, return 100%
    return 100.0

def timed_execute(cmd, run_cmd=1):
    """
    timed_execute(cmd):
        Execute the command 'cmd' after logging the command
            to STDOUT.  Return the wall-clock amount of time
            the command took to execute.
    """
    sys.stdout.write("\n'"+cmd+"'\n")
    sys.stdout.flush()
    start = time.time()
    if run_cmd:  os.system(cmd)
    end = time.time()
    return end - start

def get_folding_command(cand, obs, ddplans):
    """
    get_folding_command(cand, obs, ddplans):
        Return a command for prepfold for folding the subbands using
            an obs_info instance, a list of the ddplans, and a candidate 
            instance that describes the observations and searches.
    """
    # Folding rules are based on the facts that we want:
    #   1.  Between 24 and 200 bins in the profiles
    #   2.  For most candidates, we want to search length = 101 p/pd/DM cubes
    #       (The side of the cube is always 2*M*N+1 where M is the "factor",
    #       either -npfact (for p and pd) or -ndmfact, and N is the number of bins
    #       in the profile).  A search of 101^3 points is pretty fast.
    #   3.  For slow pulsars (where N=100 or 200), since we'll have to search
    #       many points, we'll use fewer intervals in time (-npart 30)
    #   4.  For the slowest pulsars, in order to avoid RFI, we'll
    #       not search in period-derivative.
    zmax = cand.filename.split("_")[-1]
    outfilenm = obs.basefilenm+"_DM%s_Z%s"%(cand.DMstr, zmax)
    hidms = [x.lodm for x in ddplans[1:]] + [2000]
    dfacts = [x.downsamp for x in ddplans]
    for hidm, dfact in zip(hidms, dfacts):
        if cand.DM < hidm:
            downsamp = dfact
            break
    if downsamp==1:
        filfile = obs.fil_filenm
    else:
        filfile = obs.basefilenm+"_DS%d.fil"%downsamp
    p = 1.0 / cand.f
    if (p < 0.002):
        Mp, Mdm, N = 2, 2, 24
        otheropts = "-npart 50 -ndmfact 3"
    elif p < 0.05:
        Mp, Mdm, N = 2, 1, 50
        otheropts = "-npart 40 -pstep 1 -pdstep 2 -dmstep 3"
    elif p < 0.5:
        Mp, Mdm, N = 1, 1, 100
        otheropts = "-npart 30 -pstep 1 -pdstep 2 -dmstep 1"
    else:
        Mp, Mdm, N = 1, 1, 200
        otheropts = "-npart 30 -nopdsearch -pstep 1 -pdstep 2 -dmstep 1"
    return "prepfold -noxwin -accelcand %d -accelfile %s.cand -dm %.2f -o %s %s -n %d -npfact %d -ndmfact %d %s" % \
           (cand.candnum, cand.filename, cand.DM, outfilenm,
            otheropts, N, Mp, Mdm, filfile)

class obs_info(object):
    """
    class obs_info(fil_filenm)
        A class describing the observation and the analysis.
    """
    def __init__(self, fil_filenm):
        self.fil_filenm = fil_filenm
        self.basefilenm = fil_filenm[:fil_filenm.find(".fil")]
        filhdr, hdrlen = sigproc.read_header(fil_filenm)
        self.MJD = filhdr['tstart']
        self.nchans = filhdr['nchans']
        self.ra_rad = sigproc.ra2radians(filhdr['src_raj'])
        self.ra_string = pu.coord_to_string(*pu.rad_to_hms(self.ra_rad))
        self.dec_rad = sigproc.dec2radians(filhdr['src_dej'])
        self.dec_string = pu.coord_to_string(*pu.rad_to_dms(self.dec_rad))
        self.str_coords = "J"+"".join(self.ra_string.split(":")[:2])
        if self.dec_rad >= 0.0:  self.str_coords += "+"
        self.str_coords += "".join(self.dec_string.split(":")[:2])
        self.az = filhdr['az_start']
        self.el = 90.0-filhdr['za_start']
        fillen = os.stat(fil_filenm)[6]
        self.raw_N = (fillen-hdrlen)/(filhdr['nbits']/8)/filhdr['nchans']
        self.dt = filhdr['tsamp']
        self.raw_T = self.raw_N * self.dt
        self.N = orig_N
        self.T = self.N * self.dt
        # Determine the average barycentric velocity of the observation
        self.baryv = presto.get_baryv(self.ra_string, self.dec_string,
                                      self.MJD, self.T, obs="GB")
        # Where to dump all the results
        # Directory structure is under the base_output_directory
        # according to base/MJD/filenmbase/beam
        self.outputdir = os.path.join(base_output_dir,
                                      str(int(self.MJD)),
                                      self.str_coords)
        # Figure out which host we are processing on
        self.hostname = socket.gethostname()
        # The fraction of the data recommended to be masked by rfifind
        self.masked_fraction = 0.0
        # Initialize our timers
        self.rfifind_time = 0.0
        self.downsample_time = 0.0
        self.dedispersing_time = 0.0
        self.FFT_time = 0.0
        self.lo_accelsearch_time = 0.0
        self.hi_accelsearch_time = 0.0
        self.singlepulse_time = 0.0
        self.sifting_time = 0.0
        self.folding_time = 0.0
        self.total_time = 0.0
        # Inialize some candidate counters
        self.num_sifted_cands = 0
        self.num_folded_cands = 0
        self.num_single_cands = 0
        
    def write_report(self, filenm):
        report_file = open(filenm, "w")
        report_file.write("---------------------------------------------------------\n")
        report_file.write("%s was processed on %s\n"%(self.fil_filenm, self.hostname))
        report_file.write("Ending UTC time:  %s\n"%(time.asctime(time.gmtime())))
        report_file.write("Total wall time:  %.1f s (%.2f hrs)\n"%\
                          (self.total_time, self.total_time/3600.0))
        report_file.write("Fraction of data masked:  %.2f%%\n"%\
                          (self.masked_fraction*100.0))
        report_file.write("---------------------------------------------------------\n")
        report_file.write("          rfifind time = %7.1f sec (%5.2f%%)\n"%\
                          (self.rfifind_time, self.rfifind_time/self.total_time*100.0))
        report_file.write("     dedispersing time = %7.1f sec (%5.2f%%)\n"%\
                          (self.dedispersing_time, self.dedispersing_time/self.total_time*100.0))
        report_file.write("     single-pulse time = %7.1f sec (%5.2f%%)\n"%\
                          (self.singlepulse_time, self.singlepulse_time/self.total_time*100.0))
        report_file.write("              FFT time = %7.1f sec (%5.2f%%)\n"%\
                          (self.FFT_time, self.FFT_time/self.total_time*100.0))
        report_file.write("   lo-accelsearch time = %7.1f sec (%5.2f%%)\n"%\
                          (self.lo_accelsearch_time, self.lo_accelsearch_time/self.total_time*100.0))
        report_file.write("   hi-accelsearch time = %7.1f sec (%5.2f%%)\n"%\
                          (self.hi_accelsearch_time, self.hi_accelsearch_time/self.total_time*100.0))
        report_file.write("          sifting time = %7.1f sec (%5.2f%%)\n"%\
                          (self.sifting_time, self.sifting_time/self.total_time*100.0))
        report_file.write("          folding time = %7.1f sec (%5.2f%%)\n"%\
                          (self.folding_time, self.folding_time/self.total_time*100.0))
        report_file.write("---------------------------------------------------------\n")
        report_file.close()
        
class dedisp_plan(object):
    """
    class dedisp_plan(lodm, dmstep, dmsperpass, numpasses, numsub, downsamp)
       A class describing a de-dispersion plan for prepsubband in detail.
    """
    def __init__(self, lodm, dmstep, dmsperpass, numpasses, numsub, downsamp):
        self.lodm = float(lodm)
        self.dmstep = float(dmstep)
        self.dmsperpass = int(dmsperpass)
        self.numpasses = int(numpasses)
        self.numsub = int(numsub)
        self.downsamp = int(downsamp)
        self.sub_dmstep = self.dmsperpass * self.dmstep
        self.dmlist = []  # These are strings for comparison with filenames
        self.subdmlist = []
        for ii in range(self.numpasses):
             self.subdmlist.append("%.2f"%(self.lodm + (ii+0.5)*self.sub_dmstep))
             lodm = self.lodm + ii * self.sub_dmstep
             dmlist = ["%.2f"%dm for dm in \
                       numpy.arange(self.dmsperpass)*self.dmstep + lodm]
             self.dmlist.append(dmlist)

def main(fil_filenm, workdir, ddplans):
   
    # Change to the specified working directory
    os.chdir(workdir)

    # Get information on the observation and the job
    job = obs_info(fil_filenm)
    if job.raw_T < low_T_to_search:
        print("The observation is too short (%.2f s) to search."%job.raw_T)
        sys.exit()
    job.total_time = time.time()
    ddplans = ddplans[job.nchans]
    
    # Use whatever .zaplist is found in the current directory
    default_zaplist = glob.glob("*.zaplist")[0]

    # Make sure the output directory (and parent directories) exist
    try:
        os.makedirs(job.outputdir)
        os.chmod(job.outputdir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
    except: pass

    # Make sure the tmp directory (in a tmpfs mount) exists
    tmpdir = os.path.join(base_tmp_dir, job.basefilenm)
    try:
        os.makedirs(tmpdir)
    except: pass

    print("\nBeginning GBT350 driftscan search of '%s'"%job.fil_filenm)
    print("UTC time is:  %s"%(time.asctime(time.gmtime())))

    # rfifind the filterbank file
    cmd = "rfifind -time %.17g -o %s %s > %s_rfifind.out"%\
          (rfifind_chunk_time, job.basefilenm,
           job.fil_filenm, job.basefilenm)
    job.rfifind_time += timed_execute(cmd)
    maskfilenm = job.basefilenm + "_rfifind.mask"
    # Find the fraction that was suggested to be masked
    # Note:  Should we stop processing if the fraction is
    #        above some large value?  Maybe 30%?
    job.masked_fraction = find_masked_fraction(job)
    
    # Iterate over the stages of the overall de-dispersion plan
    dmstrs = []
    for ddplan in ddplans:

        # Make a downsampled filterbank file
        if ddplan.downsamp > 1:
            cmd = "downsample_filterbank.py %d %s"%(ddplan.downsamp, job.fil_filenm)
            job.downsample_time += timed_execute(cmd)
            fil_filenm = job.fil_filenm[:job.fil_filenm.find(".fil")] + \
                         "_DS%d.fil"%ddplan.downsamp
        else:
            fil_filenm = job.fil_filenm
            
        # Iterate over the individual passes through the .fil file
        for passnum in range(ddplan.numpasses):
            subbasenm = "%s_DM%s"%(job.basefilenm, ddplan.subdmlist[passnum])

            # Now de-disperse 
            cmd = "prepsubband -mask %s -lodm %.2f -dmstep %.2f -nsub %d -numdms %d -numout %d -o %s/%s %s"%\
                  (maskfilenm, ddplan.lodm+passnum*ddplan.sub_dmstep, ddplan.dmstep,
                   ddplan.numsub, ddplan.dmsperpass, job.N/ddplan.downsamp,
                   tmpdir, job.basefilenm, fil_filenm)
            job.dedispersing_time += timed_execute(cmd)
            
            # Iterate over all the new DMs
            for dmstr in ddplan.dmlist[passnum]:
                dmstrs.append(dmstr)
                basenm = os.path.join(tmpdir, job.basefilenm+"_DM"+dmstr)
                datnm = basenm+".dat"
                fftnm = basenm+".fft"
                infnm = basenm+".inf"

                # Do the single-pulse search
                cmd = "single_pulse_search.py -p -m %f -t %f %s"%\
                      (singlepulse_maxwidth, singlepulse_threshold, datnm)
                job.singlepulse_time += timed_execute(cmd)
                try:
                    shutil.move(basenm+".singlepulse", workdir)
                except: pass

                # FFT, zap, and de-redden
                cmd = "realfft %s"%datnm
                job.FFT_time += timed_execute(cmd)
                cmd = "zapbirds -zap -zapfile %s -baryv %.6g %s"%\
                      (default_zaplist, job.baryv, fftnm)
                job.FFT_time += timed_execute(cmd)
                cmd = "rednoise %s"%fftnm
                job.FFT_time += timed_execute(cmd)
                try:
                    os.rename(basenm+"_red.fft", fftnm)
                except: pass
                
                # Do the low-acceleration search
                cmd = "accelsearch -numharm %d -sigma %f -zmax %d -flo %f %s"%\
                      (lo_accel_numharm, lo_accel_sigma, lo_accel_zmax, lo_accel_flo, fftnm)
                job.lo_accelsearch_time += timed_execute(cmd)
                try:
                    os.remove(basenm+"_ACCEL_%d.txtcand"%lo_accel_zmax)
                except: pass
                try:  # This prevents errors if there are no cand files to copy
                    shutil.move(basenm+"_ACCEL_%d.cand"%lo_accel_zmax, workdir)
                    shutil.move(basenm+"_ACCEL_%d"%lo_accel_zmax, workdir)
                except: pass
        
                # Do the high-acceleration search
                cmd = "accelsearch -numharm %d -sigma %f -zmax %d -flo %f %s"%\
                      (hi_accel_numharm, hi_accel_sigma, hi_accel_zmax, hi_accel_flo, fftnm)
                job.hi_accelsearch_time += timed_execute(cmd)
                try:
                    os.remove(basenm+"_ACCEL_%d.txtcand"%hi_accel_zmax)
                except: pass
                try:  # This prevents errors if there are no cand files to copy
                    shutil.move(basenm+"_ACCEL_%d.cand"%hi_accel_zmax, workdir)
                    shutil.move(basenm+"_ACCEL_%d"%hi_accel_zmax, workdir)
                except: pass

                # Move the .inf files
                try:
                    shutil.move(infnm, workdir)
                except: pass
                # Remove the .dat and .fft files
                try:
                    os.remove(datnm)
                except: pass
                try:
                    os.remove(fftnm)
                except: pass

    # Make the single-pulse plots
    basedmb = job.basefilenm+"_DM"
    basedme = ".singlepulse "
    # The following will make plots for DM ranges:
    #    0-30, 20-110, 100-310, 300-1000+
    dmglobs = [basedmb+"[0-9].[0-9][0-9]"+basedme +
               basedmb+"[012][0-9].[0-9][0-9]"+basedme,
               basedmb+"[2-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"10[0-9].[0-9][0-9]"+basedme,
               basedmb+"[12][0-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"30[0-9].[0-9][0-9]"+basedme,
               basedmb+"[3-9][0-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"1[0-9][0-9][0-9].[0-9][0-9]"+basedme]
    dmrangestrs = ["0-30", "20-110", "100-310", "300-1000+"]
    psname = job.basefilenm+"_singlepulse.ps"
    for dmglob, dmrangestr in zip(dmglobs, dmrangestrs):
        cmd = 'single_pulse_search.py -t %f -g "%s"' % \
              (singlepulse_plot_SNR, dmglob)
        job.singlepulse_time += timed_execute(cmd)
        try:
            os.rename(psname,
                      job.basefilenm+"_DMs%s_singlepulse.ps"%dmrangestr)
        except: pass

    # Sift through the candidates to choose the best to fold
    
    job.sifting_time = time.time()

    lo_accel_cands = sifting.read_candidates(glob.glob("*ACCEL_%d"%lo_accel_zmax))
    if len(lo_accel_cands):
        lo_accel_cands = sifting.remove_duplicate_candidates(lo_accel_cands)
    if len(lo_accel_cands):
        lo_accel_cands = sifting.remove_DM_problems(lo_accel_cands, numhits_to_fold,
                                                    dmstrs, low_DM_cutoff)
    if len(lo_accel_cands):
        lo_accel_cands.sort(key=attrgetter('sigma'), reverse=True)
        sifting.write_candlist(lo_accel_cands,
                               job.basefilenm+".accelcands_Z%d"%lo_accel_zmax)
        
    hi_accel_cands = sifting.read_candidates(glob.glob("*ACCEL_%d"%hi_accel_zmax))
    if len(hi_accel_cands):
        hi_accel_cands = sifting.remove_duplicate_candidates(hi_accel_cands)
    if len(hi_accel_cands):
        hi_accel_cands = sifting.remove_DM_problems(hi_accel_cands, numhits_to_fold,
                                                    dmstrs, low_DM_cutoff)
    if len(hi_accel_cands):
        hi_accel_cands.sort(key=attrgetter('sigma'), reverse=True)
        sifting.write_candlist(hi_accel_cands,
                               job.basefilenm+".accelcands_Z%d"%hi_accel_zmax)

    try:
        cmd = "mv *.accelcands* "+job.outputdir
        os.system(cmd)
    except: pass
    job.sifting_time = time.time() - job.sifting_time

    # Fold the best candidates

    cands_folded = 0
    for cand in lo_accel_cands:
        if cands_folded == max_lo_cands_to_fold:
            break
        elif cand.sigma > to_prepfold_sigma:
            job.folding_time += timed_execute(get_folding_command(cand, job, ddplans))
            cands_folded += 1
    cands_folded = 0
    for cand in hi_accel_cands:
        if cands_folded == max_hi_cands_to_fold:
            break
        elif cand.sigma > to_prepfold_sigma:
            job.folding_time += timed_execute(get_folding_command(cand, job, ddplans))
            cands_folded += 1
    # Remove the bestprof files
    bpfiles = glob.glob("*.pfd.bestprof")
    for bpfile in bpfiles:
        os.remove(bpfile)

    # Now step through the .ps files and convert them to .png and gzip them

    psfiles = glob.glob("*.ps")
    for psfile in psfiles:
        if "singlepulse" in psfile:
            # For some reason the singlepulse files don't transform nicely...
            epsfile = psfile.replace(".ps", ".eps")
            os.system("eps2eps "+psfile+" "+epsfile)
            os.system("pstoimg -density 100 -crop a "+epsfile)
            try:
                os.remove(epsfile)
            except: pass
        else:
            os.system("pstoimg -density 100 -flip cw "+psfile)
        os.system("gzip "+psfile)
    
    # Tar up the results files 

    tar_suffixes = ["_ACCEL_%d.tgz"%lo_accel_zmax,
                    "_ACCEL_%d.tgz"%hi_accel_zmax,
                    "_ACCEL_%d.cand.tgz"%lo_accel_zmax,
                    "_ACCEL_%d.cand.tgz"%hi_accel_zmax,
                    "_singlepulse.tgz",
                    "_inf.tgz",
                    "_pfd.tgz"]
    tar_globs = ["*_ACCEL_%d"%lo_accel_zmax,
                 "*_ACCEL_%d"%hi_accel_zmax,
                 "*_ACCEL_%d.cand"%lo_accel_zmax,
                 "*_ACCEL_%d.cand"%hi_accel_zmax,
                 "*.singlepulse",
                 "*_DM[0-9]*.inf",
                 "*.pfd"]
    for (tar_suffix, tar_glob) in zip(tar_suffixes, tar_globs):
        tf = tarfile.open(job.basefilenm+tar_suffix, "w:gz")
        for infile in glob.glob(tar_glob):
            tf.add(infile)
            os.remove(infile)
        tf.close()
            
    # Remove all the downsampled .fil files

    filfiles = glob.glob("*_DS?.fil") + glob.glob("*_DS??.fil")
    for filfile in filfiles:
        os.remove(filfile)

    # Remove the tmp directory (in a tmpfs mount)
    try:
        os.rmdir(tmpdir)
    except: pass

    # And finish up

    job.total_time = time.time() - job.total_time
    print("\nFinished")
    print("UTC time is:  %s"%(time.asctime(time.gmtime())))

    # Write the job report

    job.write_report(job.basefilenm+".report")
    job.write_report(os.path.join(job.outputdir, job.basefilenm+".report"))

    # Move all the important stuff to the output directory
    cmd = "mv *rfifind.[bimors]* *.tgz *.ps.gz *.png *.report "+\
          job.outputdir
    os.system(cmd)

if __name__ == "__main__":
    # Create our de-dispersion plans
    ddplans = {1024:[], 2048:[]}
    if (0):
        # The following are the near-optimal values for 1024 and 2048 lags.
        # They keeps the total dispersive smearing (i.e.
        # not counting scattering) <1 ms up to a DM of ~100 pc cm^-3 for 1024-lag
        # data and ~200 pc cm^-3 for 2048-lag data.
        # For 1024 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[1024].append(dedisp_plan(   0.0,  0.02,      20,    91,   32,       1))
        ddplans[1024].append(dedisp_plan(  36.4,  0.03,      24,    30,   32,       2))
        ddplans[1024].append(dedisp_plan(  58.0,  0.05,      24,    35,   32,       4))
        ddplans[1024].append(dedisp_plan( 100.0,  0.10,      24,    40,   32,       8))
        ddplans[1024].append(dedisp_plan( 196.0,  0.30,      22,    45,   32,      16))
        ddplans[1024].append(dedisp_plan( 493.0,  0.50,      24,    30,   32,      32))
        ddplans[1024].append(dedisp_plan( 853.0,  1.00,      24,     7,   32,      64))
        # For 2048 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[2048].append(dedisp_plan(   0.0,  0.02,      20,   177,   32,       1))
        ddplans[2048].append(dedisp_plan(  70.8,  0.03,      24,    60,   32,       2))
        ddplans[2048].append(dedisp_plan( 114.0,  0.05,      24,    65,   32,       4))
        ddplans[2048].append(dedisp_plan( 192.0,  0.10,      24,    80,   32,       8))
        ddplans[2048].append(dedisp_plan( 384.0,  0.30,      22,    80,   32,      16))
        ddplans[2048].append(dedisp_plan( 912.0,  0.50,      24,     8,   32,      32))
    elif (0):
        #
        # If there is <=1GB of RAM per node, the following are preferred
        #
        # DDplan.py -f 350.0 -b 50.0 -n 1024 -t 0.00008192 -s 64 -r 0.2
        # For 1024 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[1024].append(dedisp_plan(   0.0,  0.03,      50,    37,   64,       2))
        ddplans[1024].append(dedisp_plan(  55.5,  0.05,      50,    17,   64,       4))
        ddplans[1024].append(dedisp_plan(  98.0,  0.10,      50,    19,   64,       8))
        ddplans[1024].append(dedisp_plan( 193.0,  0.20,      50,    19,   64,      16))
        ddplans[1024].append(dedisp_plan( 383.0,  0.50,      50,    19,   64,      32))
        ddplans[1024].append(dedisp_plan( 858.0,  1.00,      50,     3,   64,      64))
        # DDplan.py -f 350.0 -b 50.0 -n 2048 -t 0.00008192 -s 64 -r 0.2
        # For 2048 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[2048].append(dedisp_plan(   0.0,  0.03,      50,    74,   64,       2))
        ddplans[2048].append(dedisp_plan( 111.0,  0.05,      50,    33,   64,       4))
        ddplans[2048].append(dedisp_plan( 193.5,  0.10,      50,    38,   64,       8))
        ddplans[2048].append(dedisp_plan( 383.5,  0.20,      50,    38,   64,      16))
        ddplans[2048].append(dedisp_plan( 763.5,  0.50,      50,    10,   64,      32))
    elif (1):
        #
        # If there is 2GB or more RAM per node, the following are probably faster
        #
        # DDplan.py -f 350.0 -b 50.0 -n 1024 -t 0.00008192 -s 128 -r 0.2
        # For 1024 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[1024].append(dedisp_plan(   0.0,  0.03,     100,    19,  128,       2))
        ddplans[1024].append(dedisp_plan(  57.0,  0.05,     100,     8,  128,       4))
        ddplans[1024].append(dedisp_plan(  97.0,  0.10,     100,    10,  128,       8))
        ddplans[1024].append(dedisp_plan( 197.0,  0.20,     100,    10,  128,      16))
        ddplans[1024].append(dedisp_plan( 397.0,  0.50,     100,    10,  128,      32))
        ddplans[1024].append(dedisp_plan( 897.0,  1.00,     100,     2,  128,      64))
        # DDplan.py -f 350.0 -b 50.0 -n 2048 -t 0.00008192 -s 128 -r 0.2
        # For 2048 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[2048].append(dedisp_plan(   0.0,  0.03,     100,    37,  128,       2))
        ddplans[2048].append(dedisp_plan( 111.0,  0.05,     100,    17,  128,       4))
        ddplans[2048].append(dedisp_plan( 196.0,  0.10,     100,    19,  128,       8))
        ddplans[2048].append(dedisp_plan( 386.0,  0.20,     100,    19,  128,      16))
        ddplans[2048].append(dedisp_plan( 766.0,  0.50,     100,     5,  128,      32))
    elif (0):
        #
        # This is for "quick" processing
        #
        # DDplan.py -f 350.0 -b 50.0 -n 1024 -t 0.00008192 -s 128 -r 1.5
        # For 1024 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[1024].append(dedisp_plan(   0.0,  0.20,     100,    20,  128,      16))
        ddplans[1024].append(dedisp_plan( 400.0,  0.50,     100,    10,  128,      32))
        ddplans[1024].append(dedisp_plan( 900.0,  1.00,     100,     2,  128,      64))
        # DDplan.py -f 350.0 -b 50.0 -n 2048 -t 0.00008192 -s 128 -r 1.5
        # For 2048 chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans[2048].append(dedisp_plan(   0.0,  0.20,     100,    39,  128,      16))
        ddplans[2048].append(dedisp_plan( 780.0,  0.50,     100,     5,  128,      32))

    # Arguments to the search program are
    # sys.argv[1] = filterbank file name
    # sys.argv[2] = working directory name
    if len(sys.argv) >= 3:
        workdir = sys.argv[2]
        fil_filenm = sys.argv[1]
        main(fil_filenm, workdir, ddplans)
    elif len(sys.argv) == 2:
        fil_filenm = sys.argv[1]
        main(fil_filenm, '.', ddplans)
    else:
        print("GBT350_drift_search.py fil_filenm [workdir]")
