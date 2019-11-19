#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from operator import attrgetter
import glob, os, os.path, shutil, socket, tarfile, stat
import numpy, sys, time
from presto import presto
from presto import sifting
import astropy.io.fits as pyfits

institution = "NRAO" 
base_tmp_dir = "."  # "/dev/shm/" is a good choice
# This is where the output will be archived.
base_output_dir = "."

#-------------------------------------------------------------------
# Tunable parameters for searching and folding
# (you probably don't need to tune any of them)
raw_N                 = 1440000 # Number of samples to analyze (~118 secs)
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
foldnsubs               = 128    # Number of subbands to use when folding
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

def get_folding_command(cand, obs, ddplans, maskfile):
    """
    get_folding_command(cand, obs, ddplans, maskfile):
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
        fitsfile = obs.fits_filenm
    else:
        fitsfile = obs.dsbasefilenm+"_DS%d%s"%(downsamp,obs.fits_filenm[obs.fits_filenm.rfind("_"):])
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
    return "prepfold -mask %s -noxwin -accelcand %d -accelfile %s.cand -dm %.2f -o %s %s -n %d -npfact %d -ndmfact %d -nsub %d %s" % \
           (maskfile, cand.candnum, cand.filename, cand.DM, outfilenm,
            otheropts, N, Mp, Mdm, foldnsubs, fitsfile)

class obs_info(object):
    """
    class obs_info(fits_filenm)
        A class describing the observation and the analysis.
    """
    def __init__(self, fits_filenm):
        self.fits_filenm = fits_filenm
        self.basefilenm = fits_filenm[:fits_filenm.find(".fits")]
        self.dsbasefilenm = fits_filenm[:fits_filenm.rfind("_")]
        fitshandle=pyfits.open(fits_filenm)
        self.MJD = fitshandle[0].header['STT_IMJD']+fitshandle[0].header['STT_SMJD']/86400.0+fitshandle[0].header['STT_OFFS']/86400.0
        self.nchans = fitshandle[0].header['OBSNCHAN']
        self.ra_string = fitshandle[0].header['RA']
        self.dec_string = fitshandle[0].header['DEC']
        self.str_coords = "J"+"".join(self.ra_string.split(":")[:2])
        self.str_coords += "".join(self.dec_string.split(":")[:2])
        self.nbits=fitshandle[0].header['BITPIX']
        
        self.raw_N=fitshandle[1].header['NAXIS2']*fitshandle[1].header['NSBLK']
        self.dt=fitshandle[1].header['TBIN']*1000000
        self.raw_T = self.raw_N * self.dt
        self.N = raw_N
        if self.dt == 163.84:
          self.N=self.N/2
        self.T = self.N * self.dt
        self.srcname=fitshandle[0].header['SRC_NAME']
        # Determine the average barycentric velocity of the observation
        self.baryv = presto.get_baryv(self.ra_string, self.dec_string,
                                      self.MJD, self.T, obs="GB")
        # Where to dump all the results
        # Directory structure is under the base_output_directory
        # according to base/MJD/filenmbase/beam
        self.outputdir = os.path.join(base_output_dir,
                                      str(int(self.MJD)),
                                      self.srcname)
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
        report_file.write("%s was processed on %s\n"%(self.fits_filenm, self.hostname))
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

def remove_crosslist_duplicate_candidates(candlist1,candlist2):
    n1 = len(candlist1)
    n2 = len(candlist2)
    removelist1 = []
    removelist2 = []
    candlist2.sort(key=attrgetter('r'))
    candlist1.sort(key=attrgetter('r'))
    print("  Searching for crosslist dupes...")
    ii = 0
    while ii < n1:
        jj=0
        while jj < n2:
            if numpy.fabs(candlist1[ii].r-candlist2[jj].r) < sifting.r_err:
                if candlist1[ii].sigma > candlist2[jj].sigma:
                    print("Crosslist remove from candlist 2, %f > %f, %d:%f~%f" % \
                          (candlist1[ii].sigma, candlist2[jj].sigma, jj,
                           candlist1[ii].r, candlist2[jj].r))
                    if jj not in removelist2:
                        removelist2.append(jj)
                else:
                    print("Crosslist remove from candlist 1, %f > %f, %d:%f~%f" % \
                          (candlist2[jj].sigma, candlist1[ii].sigma, ii,
                           candlist1[ii].r, candlist2[jj].r))
                    if ii not in removelist1:
                        removelist1.append(ii)
            jj += 1
        ii += 1
    for ii in range(len(removelist2)-1,-1,-1):
        print("Removing %d from candlist2" % removelist2[ii])
        del(candlist2[removelist2[ii]])
    for ii in range(len(removelist1)-1,-1,-1):
        print("Removing %d from candlist1" % removelist1[ii])
        del(candlist1[removelist1[ii]])
    print("Removed %d crosslist candidates\n" % (len(removelist1)+len(removelist2)))
    print("Found %d candidates.  Sorting them by significance...\n" % (len(candlist1)+len(candlist2)))
    candlist1.sort(key=attrgetter('sigma'), reverse=True)
    candlist2.sort(key=attrgetter('sigma'), reverse=True)
    return candlist1,candlist2


def main(fits_filenm, workdir, ddplans):
   
    # Change to the specified working directory
    os.chdir(workdir)

    # Get information on the observation and the job
    job = obs_info(fits_filenm)
    if job.raw_T < low_T_to_search:
        print("The observation is too short (%.2f s) to search."%job.raw_T)
        sys.exit()
    job.total_time = time.time()
    if job.dt == 163.84:
        ddplans = ddplans[str(job.nchans)+"slow"]
    else:
        ddplans = ddplans[str(job.nchans)+"fast"]
    
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

    print("\nBeginning GBNCC search of '%s'"%job.fits_filenm)
    print("UTC time is:  %s"%(time.asctime(time.gmtime())))

    rfifindout=job.basefilenm+"_rfifind.out"
    rfifindmask=job.basefilenm+"_rfifind.mask"

    if not os.path.exists(rfifindout) or not os.path.exists(rfifindmask):
        # rfifind the filterbank file
        cmd = "rfifind -time %.17g -o %s %s > %s_rfifind.out"%\
              (rfifind_chunk_time, job.basefilenm,
               job.fits_filenm, job.basefilenm)
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
            cmd = "psrfits_subband -dstime %d -nsub %d -o %s_DS%d %s"%\
                  (ddplan.downsamp, job.nchans, job.dsbasefilenm, ddplan.downsamp, job.dsbasefilenm )
            job.downsample_time += timed_execute(cmd)
            fits_filenm = job.dsbasefilenm + "_DS%d%s"%\
                          (ddplan.downsamp,job.fits_filenm[job.fits_filenm.rfind("_"):])
        else:
            fits_filenm = job.fits_filenm

        # Iterate over the individual passes through the .fil file
        for passnum in range(ddplan.numpasses):
            subbasenm = "%s_DM%s"%(job.basefilenm, ddplan.subdmlist[passnum])

            # Now de-disperse 
            cmd = "prepsubband -mask %s -lodm %.2f -dmstep %.2f -nsub %d -numdms %d -numout %d -o %s/%s %s"%\
                  (maskfilenm, ddplan.lodm+passnum*ddplan.sub_dmstep,
                   ddplan.dmstep, ddplan.numsub,
                   ddplan.dmsperpass, job.N/ddplan.downsamp,
                   tmpdir, job.basefilenm, fits_filenm)
            job.dedispersing_time += timed_execute(cmd)
            
            # Do the single-pulse search
            cmd = "single_pulse_search.py -p -m %f -t %f %s/*.dat"%\
                (singlepulse_maxwidth, singlepulse_threshold, tmpdir)
            job.singlepulse_time += timed_execute(cmd)
            spfiles = glob.glob("%s/*.singlepulse"%tmpdir)
            for spfile in spfiles:
                try:
                    shutil.move(spfile, workdir)
                except: pass

            # Iterate over all the new DMs
            for dmstr in ddplan.dmlist[passnum]:
                dmstrs.append(dmstr)
                basenm = os.path.join(tmpdir, job.basefilenm+"_DM"+dmstr)
                datnm = basenm+".dat"
                fftnm = basenm+".fft"
                infnm = basenm+".inf"

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
        
    hi_accel_cands = sifting.read_candidates(glob.glob("*ACCEL_%d"%hi_accel_zmax))
    if len(hi_accel_cands):
        hi_accel_cands = sifting.remove_duplicate_candidates(hi_accel_cands)
    if len(hi_accel_cands):
        hi_accel_cands = sifting.remove_DM_problems(hi_accel_cands, numhits_to_fold,
                                                    dmstrs, low_DM_cutoff)

    if len(lo_accel_cands) and len(hi_accel_cands):
        lo_accel_cands, hi_accel_cands = remove_crosslist_duplicate_candidates(lo_accel_cands, hi_accel_cands)

    if len(lo_accel_cands):
        lo_accel_cands.sort(key=attrgetter('sigma'), reverse=True)
        sifting.write_candlist(lo_accel_cands,
                               job.basefilenm+".accelcands_Z%d"%lo_accel_zmax)
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
            job.folding_time += timed_execute(get_folding_command(cand, job, ddplans, maskfilenm))
            cands_folded += 1
    cands_folded = 0
    for cand in hi_accel_cands:
        if cands_folded == max_hi_cands_to_fold:
            break
        elif cand.sigma > to_prepfold_sigma:
            job.folding_time += timed_execute(get_folding_command(cand, job, ddplans, maskfilenm))
            cands_folded += 1
    # Remove the bestprof files
    bpfiles = glob.glob("*.pfd.bestprof")
    for bpfile in bpfiles:
        os.remove(bpfile)

    # Now step through the .ps files and convert them to .png and gzip them

    psfiles = glob.glob("*.ps")
    for psfile in psfiles:
        if "singlepulse" in psfile:
            os.system("pstoimg -density 200 -antialias -crop a "+psfile)
            try:
                os.remove(psfile)
            except: pass
        else:
            os.system("pstoimg -density 200 -antialias -flip cw "+psfile)
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
            
    # Remove all the downsampled .fits files

    fitsfiles = glob.glob("*_DS?*.fits") + glob.glob("*_DS??*.fits")
    for fitsfile in fitsfiles:
        os.remove(fitsfile)

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
    # All GBNCC data have 4096 channels, but the earliest data is sampled
    # at 163.84us rather than 81.92 us...
    ddplans = {'4096slow':[], '4096fast':[]}
    if (1):
        #
        # If there is <=1GB of RAM per CPU core, the following are preferred
        #
        # For 4096slow chan data:               lodm    dmstep  dms/call #calls #subs downsamp
        ddplans['4096slow'].append(dedisp_plan(   0.0,  0.02,      86,    81,  128,       1))
        ddplans['4096slow'].append(dedisp_plan(139.32,  0.03,     102,    27,  128,       2))
        ddplans['4096slow'].append(dedisp_plan(221.94,  0.05,     102,    33,  128,       4))
        ddplans['4096slow'].append(dedisp_plan(390.24,  0.10,     102,    11,  128,       8))
        # For 4096fast chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans['4096fast'].append(dedisp_plan(   0.0,  0.01,      86,    81,  128,       1))
        ddplans['4096fast'].append(dedisp_plan( 69.66,  0.02,      86,    33,  128,       2))
        ddplans['4096fast'].append(dedisp_plan(126.42,  0.03,     102,    29,  128,       4))
        ddplans['4096fast'].append(dedisp_plan(215.16,  0.05,     102,    33,  128,       8))
        ddplans['4096fast'].append(dedisp_plan(383.46,  0.10,     102,    12,  128,      16))
    else:
        # If there is >2GB of RAM per CPU core, the following are preferred
        #
        # For 4096slow chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans['4096slow'].append(dedisp_plan(   0.0,  0.02,     172,    41,  256,       1))
        ddplans['4096slow'].append(dedisp_plan(141.04,  0.03,     204,    14,  256,       2))
        ddplans['4096slow'].append(dedisp_plan(226.72,  0.05,     204,    16,  256,       4))
        ddplans['4096slow'].append(dedisp_plan(389.92,  0.10,     204,     6,  256,       8))
        # For 4096fast chan data:               lodm dmstep dms/call #calls #subs downsamp
        ddplans['4096fast'].append(dedisp_plan(   0.0,  0.01,     172,    41,  256,       1))
        ddplans['4096fast'].append(dedisp_plan( 70.52,  0.02,     172,    16,  256,       2))
        ddplans['4096fast'].append(dedisp_plan(125.56,  0.03,     204,    15,  256,       4))
        ddplans['4096fast'].append(dedisp_plan(217.36,  0.05,     204,    17,  256,       8))
        ddplans['4096fast'].append(dedisp_plan(390.76,  0.10,     204,     6,  256,      16))

    # Arguments to the search program are
    # sys.argv[1] = PSRFITS file name
    # sys.argv[2] = working directory name
    if len(sys.argv) >= 3:
        workdir = sys.argv[2]
        fits_filenm = sys.argv[1]
        main(fits_filenm, workdir, ddplans)
    elif len(sys.argv) == 2:
        fits_filenm = sys.argv[1]
        main(fits_filenm, '.', ddplans)
    else:
        print("GBNCC_search.py fits_filenm [workdir]")
