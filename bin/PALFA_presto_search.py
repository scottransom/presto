#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from operator import attrgetter
import glob, os, os.path, socket, struct, sys, time, tarfile
import numpy
from presto import psr_utils
from presto import presto
from presto import sifting
from presto import sigproc

# Calling convention:
#
# PALFA_presto_search.py fil_file working_dir
#
#   fil_file is the filterbank file name
#   working_dir is the scratch directory where the work should be done
#       In general, there should be ~30GB of scratch disk per beam.
#       If you do not have that much scratch space, you will likely
#       need to to use set the use_subbands flag below.

# Basic parameters
# institution is one of: 'UBC', 'NRAOCV', 'McGill', 'Columbia', 'Cornell', 'UTB'
institution           = "NRAOCV" 
base_output_directory = "/home/sransom/results/ALFA"
db_pointing_file      = "/home/sransom/results/ALFA/PALFA_coords_table.txt"

# The following determines if we'll dedisperse and fold using subbands.
# In general, it is a very good idea to use them if there is enough scratch
# space on the machines that are processing (~30GB/beam processed)
use_subbands          = True

# Tunable parameters for searching and folding
# (you probably don't need to tune any of them)
rfifind_chunk_time    = 2**15 * 0.000064  # ~2.1 sec for dt = 64us 
singlepulse_threshold = 5.0  # threshold SNR for candidate determination
singlepulse_plot_SNR  = 6.0  # threshold SNR for singlepulse plot
singlepulse_maxwidth  = 0.1  # max pulse width in seconds
to_prepfold_sigma     = 6.0  # incoherent sum significance to fold candidates
max_cands_to_fold     = 150   # Never fold more than this many candidates
numhits_to_fold       = 2    # Number of DMs with a detection needed to fold
low_DM_cutoff         = 2.0  # Lowest DM to consider as a "real" pulsar
lo_accel_numharm      = 16   # max harmonics
lo_accel_sigma        = 2.0  # threshold gaussian significance
lo_accel_zmax         = 0    # bins
lo_accel_flo          = 2.0  # Hz
hi_accel_numharm      = 8    # max harmonics
hi_accel_sigma        = 3.0  # threshold gaussian significance
hi_accel_zmax         = 50   # bins
hi_accel_flo          = 1.0  # Hz
low_T_to_search       = 20.0 # sec

# Sifting specific parameters (don't touch without good reason!)
sifting.sigma_threshold = to_prepfold_sigma-1.0  # incoherent power threshold (sigma)
sifting.c_pow_threshold = 100.0                  # coherent power threshold
sifting.r_err           = 1.1    # Fourier bin tolerence for candidate equivalence
sifting.short_period    = 0.0005 # Shortest period candidates to consider (s)
sifting.long_period     = 15.0   # Longest period candidates to consider (s)
sifting.harm_pow_cutoff = 8.0    # Power required in at least one harmonic

def fix_fil_posn(fil_filenm, hdrlen, ra, dec):
    """
    fix_fil_posn(fil_filenm, hdrlen, ra, dec):
        Modify the filterbank header and update the RA and DEC
            fields using the values given as input.  ra and dec
            should be in 'HH:MM:SS.SSSS' and 'DD:MM:SS.SSSS' format.
            hdrlen is the length of the filterbank header as
            reported by PRESTO's 'readfile' or SIGPROC's 'header'.
    """
    newra = float(ra.replace(":", ""))
    newdec = float(dec.replace(":", ""))
    header = open(fil_filenm).read(hdrlen)
    ra_ptr = header.find("src_raj")+len("src_raj")
    dec_ptr = header.find("src_dej")+len("src_dej")
    filfile = open(fil_filenm, 'rb+')
    filfile.seek(ra_ptr)
    filfile.write(struct.pack('d', newra))
    filfile.seek(dec_ptr)
    filfile.write(struct.pack('d', newdec))
    filfile.close()

def read_db_posn(orig_filenm, beam):
    """
    read_db_posn(orig_filenm, beam):
        Find the original WAPP filename in the db_pointing_file
           and return the sexagesimal position strings for
           the choen beam in that file.  Return None if not found.
    """
    offset = beam % 2
    for line in open(db_pointing_file):
        sline = line.split()
        if sline[0].strip() == orig_filenm:
            ra_str = sline[2*offset+1].strip()
            dec_str = sline[2*offset+2].strip()
            return ra_str, dec_str
    return None

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

def get_all_subdms(ddplans):
    """
    get_all_subdms(ddplans):
        Return a sorted array of the subdms from the list of ddplans.
    """
    subdmlist = []
    for ddplan in ddplans:
        subdmlist += [float(x) for x in ddplan.subdmlist]
    subdmlist.sort()
    subdmlist = numpy.asarray(subdmlist)
    return subdmlist

def find_closest_subbands(obs, subdms, DM):
    """
    find_closest_subbands(obs, subdms, DM):
        Return the basename of the closest set of subbands to DM
        given an obs_info class and a sorted array of the subdms.
    """
    subdm = subdms[numpy.fabs(subdms - DM).argmin()]
    return "subbands/%s_DM%.2f.sub[0-6]*"%(obs.basefilenm, subdm)

def timed_execute(cmd): 
    """
    timed_execute(cmd):
        Execute the command 'cmd' after logging the command
            to STDOUT.  Return the wall-clock amount of time
            the command took to execute.
    """
    sys.stdout.write("\n'"+cmd+"'\n")
    sys.stdout.flush()
    start = time.time()
    os.system(cmd)
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

    # Note:  the following calculations should probably only be done once,
    #        but in general, these calculation are effectively instantaneous
    #        compared to the folding itself
    if use_subbands:  # Fold the subbands
        subdms = get_all_subdms(ddplans)
        subfiles = find_closest_subbands(obs, subdms, cand.DM)
        foldfiles = subfiles
    else:  # Folding the downsampled filterbank files instead
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
        foldfiles = filfile
    p = 1.0 / cand.f
    if p < 0.002:
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
            otheropts, N, Mp, Mdm, foldfiles)

class obs_info(object):
    """
    class obs_info(fil_filenm)
        A class describing the observation and the analysis.
    """
    def __init__(self, fil_filenm):
        self.fil_filenm = fil_filenm
        self.basefilenm = fil_filenm.rstrip(".fil")
        self.beam = int(self.basefilenm[-1])
        filhdr, self.hdrlen = sigproc.read_header(fil_filenm)
        self.orig_filenm = filhdr['rawdatafile']
        self.MJD = filhdr['tstart']
        self.nchans = filhdr['nchans']
        self.ra_rad = sigproc.ra2radians(filhdr['src_raj'])
        self.ra_string = psr_utils.coord_to_string(\
            *psr_utils.rad_to_hms(self.ra_rad))
        self.dec_rad = sigproc.dec2radians(filhdr['src_dej'])
        self.dec_string = psr_utils.coord_to_string(\
            *psr_utils.rad_to_dms(self.dec_rad))
        self.az = filhdr['az_start']
        self.el = 90.0-filhdr['za_start']
        self.BW = abs(filhdr['foff']) * filhdr['nchans']
        self.dt = filhdr['tsamp']
        self.orig_N = sigproc.samples_per_file(fil_filenm, filhdr, self.hdrlen)
        self.orig_T = self.orig_N * self.dt
        self.N = psr_utils.choose_N(self.orig_N)
        self.T = self.N * self.dt
        # Update the RA and DEC from the database file if required
        newposn = read_db_posn(self.orig_filenm, self.beam)
        if newposn is not None:
            self.ra_string, self.dec_string = newposn
            # ... and use them to update the filterbank file
            fix_fil_posn(fil_filenm, self.hdrlen,
                         self.ra_string, self.dec_string)
        # Determine the average barycentric velocity of the observation
        self.baryv = presto.get_baryv(self.ra_string, self.dec_string,
                                      self.MJD, self.T, obs="AO")
        # Where to dump all the results
        # Directory structure is under the base_output_directory
        # according to base/MJD/filenmbase/beam
        self.outputdir = os.path.join(base_output_directory,
                                      str(int(self.MJD)),
                                      self.basefilenm[:-2],
                                      str(self.beam))
        # Figure out which host we are processing on
        self.hostname = socket.gethostname()
        # The fraction of the data recommended to be masked by rfifind
        self.masked_fraction = 0.0
        # Initialize our timers
        self.rfifind_time = 0.0
        self.downsample_time = 0.0
        self.subbanding_time = 0.0
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
        if use_subbands:
            report_file.write("       subbanding time = %7.1f sec (%5.2f%%)\n"%\
                              (self.subbanding_time, self.subbanding_time/self.total_time*100.0))
        else:
            report_file.write("     downsampling time = %7.1f sec (%5.2f%%)\n"%\
                              (self.downsample_time, self.downsample_time/self.total_time*100.0))
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
        # Downsample less for the subbands so that folding
        # candidates is more accurate
        self.sub_downsamp = self.downsamp / 2
        if self.sub_downsamp==0: self.sub_downsamp = 1
        # The total downsampling is:
        #   self.downsamp = self.sub_downsamp * self.dd_downsamp
        if self.downsamp==1: self.dd_downsamp = 1
        else: self.dd_downsamp = 2
        self.sub_dmstep = self.dmsperpass * self.dmstep
        self.dmlist = []  # These are strings for comparison with filenames
        self.subdmlist = []
        for ii in range(self.numpasses):
            self.subdmlist.append("%.2f"%(self.lodm + (ii+0.5)*self.sub_dmstep))
            lodm = self.lodm + ii * self.sub_dmstep
            dmlist = ["%.2f"%dm for dm in \
                      numpy.arange(self.dmsperpass)*self.dmstep + lodm]
            self.dmlist.append(dmlist)

# Create our de-dispersion plans (for 100MHz WAPPs)
# The following are the "optimal" values for the 100MHz
# survey.  It keeps the total dispersive smearing (i.e.
# not counting scattering <1 ms up to a DM of ~600 pc cm^-3
ddplans = []
if (1):
    # The values here are:       lodm dmstep dms/call #calls #subbands downsamp
    ddplans.append(dedisp_plan(   0.0,   0.3,      24,    26,       32,       1))
    ddplans.append(dedisp_plan( 187.2,   0.5,      24,    10,       32,       2))
    ddplans.append(dedisp_plan( 307.2,   1.0,      24,    11,       32,       4))
    ddplans.append(dedisp_plan( 571.2,   3.0,      24,     6,       32,       8))
else: # faster option that sacrifices a small amount of time resolution at the lowest DMs
    # The values here are:       lodm dmstep dms/call #calls #subbands downsamp
    ddplans.append(dedisp_plan(   0.0,   0.5,      22,    21,       32,       1))
    ddplans.append(dedisp_plan( 231.0,   0.5,      24,     6,       32,       2))
    ddplans.append(dedisp_plan( 303.0,   1.0,      24,    11,       32,       4))
    ddplans.append(dedisp_plan( 567.0,   3.0,      24,     7,       32,       8))
    
def main(fil_filenm, workdir):

    # Change to the specified working directory
    os.chdir(workdir)

    # Get information on the observation and the jbo
    job = obs_info(fil_filenm)
    if job.T < low_T_to_search:
        print("The observation is too short (%.2f s) to search."%job.T)
        sys.exit()
    job.total_time = time.time()
    
    # Use whatever .zaplist is found in the current directory
    default_zaplist = glob.glob("*.zaplist")[0]

    # Make sure the output directory (and parent directories) exist
    try:
        os.makedirs(job.outputdir)
    except: pass

    # Create a directory to hold all the subbands
    if use_subbands:
        try:
            os.makedirs("subbands")
        except: pass
    
    print("\nBeginning PALFA search of '%s'"%job.fil_filenm)
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

        # Make a downsampled filterbank file if we are not using subbands
        if not use_subbands:
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

            if use_subbands:
                # Create a set of subbands
                cmd = "prepsubband -sub -subdm %s -downsamp %d -nsub %d -mask %s -o subbands/%s %s > %s.subout"%\
                      (ddplan.subdmlist[passnum], ddplan.sub_downsamp,
                       ddplan.numsub, maskfilenm, job.basefilenm,
                       job.fil_filenm, subbasenm)
                job.subbanding_time += timed_execute(cmd)
            
                # Now de-disperse using the subbands
                cmd = "prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -numout %d -o %s subbands/%s.sub[0-9]* > %s.prepout"%\
                      (ddplan.lodm+passnum*ddplan.sub_dmstep, ddplan.dmstep,
                       ddplan.dmsperpass, ddplan.dd_downsamp, job.N/ddplan.downsamp,
                       job.basefilenm, subbasenm, subbasenm)
                job.dedispersing_time += timed_execute(cmd)

            else:  # Not using subbands
                cmd = "prepsubband -mask %s -lodm %.2f -dmstep %.2f -numdms %d -numout %d -o %s %s"%\
                      (maskfilenm, ddplan.lodm+passnum*ddplan.sub_dmstep, ddplan.dmstep,
                       ddplan.dmsperpass, job.N/ddplan.downsamp,
                       job.basefilenm, fil_filenm)
                job.dedispersing_time += timed_execute(cmd)
            
            # Iterate over all the new DMs
            for dmstr in ddplan.dmlist[passnum]:
                dmstrs.append(dmstr)
                basenm = job.basefilenm+"_DM"+dmstr
                datnm = basenm+".dat"
                fftnm = basenm+".fft"
                infnm = basenm+".inf"

                # Do the single-pulse search
                cmd = "single_pulse_search.py -p -m %f -t %f %s"%\
                      (singlepulse_maxwidth, singlepulse_threshold, datnm)
                job.singlepulse_time += timed_execute(cmd)

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
                cmd = "accelsearch -locpow -harmpolish -numharm %d -sigma %f -zmax %d -flo %f %s"%\
                      (lo_accel_numharm, lo_accel_sigma, lo_accel_zmax, lo_accel_flo, fftnm)
                job.lo_accelsearch_time += timed_execute(cmd)
                try:
                    os.remove(basenm+"_ACCEL_%d.txtcand"%lo_accel_zmax)
                except: pass
        
                # Do the high-acceleration search
                cmd = "accelsearch -locpow -harmpolish -numharm %d -sigma %f -zmax %d -flo %f %s"%\
                      (hi_accel_numharm, hi_accel_sigma, hi_accel_zmax, hi_accel_flo, fftnm)
                job.hi_accelsearch_time += timed_execute(cmd)
                try:
                    os.remove(basenm+"_ACCEL_%d.txtcand"%hi_accel_zmax)
                except: pass

                # Remove the .dat and .fft files
                try:
                    os.remove(datnm)
                    os.remove(fftnm)
                except: pass

    # Make the single-pulse plots
    basedmb = job.basefilenm+"_DM"
    basedme = ".singlepulse "
    # The following will make plots for DM ranges:
    #    0-110, 100-310, 300-1000+
    dmglobs = [basedmb+"[0-9].[0-9][0-9]"+basedme +
               basedmb+"[0-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"10[0-9].[0-9][0-9]"+basedme,
               basedmb+"[12][0-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"30[0-9].[0-9][0-9]"+basedme,
               basedmb+"[3-9][0-9][0-9].[0-9][0-9]"+basedme +
               basedmb+"1[0-9][0-9][0-9].[0-9][0-9]"+basedme]
    dmrangestrs = ["0-110", "100-310", "300-1000+"]
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

    all_accel_cands = lo_accel_cands + hi_accel_cands
    if len(all_accel_cands):
        all_accel_cands = sifting.remove_harmonics(all_accel_cands)
        # Note:  the candidates will be sorted in _sigma_ order, not _SNR_!
        all_accel_cands.sort(key=attrgetter('sigma'), reverse=True)
        sifting.write_candlist(all_accel_cands, job.basefilenm+".accelcands")

    try:
        cmd = "cp *.accelcands "+job.outputdir
        os.system(cmd)
    except: pass

    job.sifting_time = time.time() - job.sifting_time

    # Fold the best candidates

    cands_folded = 0
    for cand in all_accel_cands:
        if cands_folded == max_cands_to_fold:
            break
        if cand.sigma > to_prepfold_sigma:
            job.folding_time += timed_execute(get_folding_command(cand, job, ddplans))
            cands_folded += 1

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
    
    # NOTE:  need to add database commands

    # Tar up the results files 

    tar_suffixes = ["_ACCEL_%d.tgz"%lo_accel_zmax,
                    "_ACCEL_%d.tgz"%hi_accel_zmax,
                    "_ACCEL_%d.cand.tgz"%lo_accel_zmax,
                    "_ACCEL_%d.cand.tgz"%hi_accel_zmax,
                    "_singlepulse.tgz",
                    "_inf.tgz",
                    "_pfd.tgz",
                    "_bestprof.tgz"]
    tar_globs = ["*_ACCEL_%d"%lo_accel_zmax,
                 "*_ACCEL_%d"%hi_accel_zmax,
                 "*_ACCEL_%d.cand"%lo_accel_zmax,
                 "*_ACCEL_%d.cand"%hi_accel_zmax,
                 "*.singlepulse",
                 "*_DM[0-9]*.inf",
                 "*.pfd",
                 "*.pfd.bestprof"]
    for (tar_suffix, tar_glob) in zip(tar_suffixes, tar_globs):
        tf = tarfile.open(job.basefilenm+tar_suffix, "w:gz")
        for infile in glob.glob(tar_glob):
            tf.add(infile)
            os.remove(infile)
        tf.close()
           

    # And finish up

    job.total_time = time.time() - job.total_time
    print("\nFinished")
    print("UTC time is:  %s"%(time.asctime(time.gmtime())))

    # Write the job report

    job.write_report(job.basefilenm+".report")
    job.write_report(os.path.join(job.outputdir, job.basefilenm+".report"))

    # Copy all the important stuff to the output directory
    try:
        cmd = "cp *rfifind.[bimors]* *.ps.gz *.tgz *.png "+job.outputdir
        os.system(cmd)
    except: pass
    
if __name__ == "__main__":
    # Arguments to the search program are
    # sys.argv[1] = filterbank file name
    # sys.argv[2] = working directory name
    fil_filenm = sys.argv[1]
    workdir = sys.argv[2]
    main(fil_filenm, workdir)
