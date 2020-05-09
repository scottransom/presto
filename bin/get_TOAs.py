#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import getopt, sys
from presto import fftfit
from presto import psr_utils
import numpy as Num
from presto.prepfold import pfd
from presto.polycos import polycos
from presto.psr_constants import *

scopes = {'GBT':'1',
          'Arecibo':'3',
          'Parkes':'7',
          'GMRT': 'r',
          'IRAM': 's',
          'LWA1': 'x',
          'LWA': 'x',
          'VLA': 'c',
          'FAST': 'k',
          'MeerKAT': 'm',
          'Geocenter': 'o'}

scopes2 = {'GBT':'gbt',
          'Arecibo':'ao',
          'Parkes':'pks',
          'GMRT': 'gmrt',
          'LWA1': 'lwa1',
          'LWA': 'lwa1',
          'VLA': 'vla',
          'FAST': 'fast',
          'MeerKAT': 'mk',
          'Geocenter': 'coe'}

def measure_phase(profile, template, rotate_prof=True):
    """
    measure_phase(profile, template):
        Call FFTFIT on the profile and template to determine the
            following parameters: shift,eshift,snr,esnr,b,errb,ngood
            (returned as a tuple).  These are defined as in Taylor's
            talk at the Royal Society.
    """
    c,amp,pha = fftfit.cprof(template)
    pha1 = pha[0]
    if (rotate_prof):
        pha = Num.fmod(pha-Num.arange(1,len(pha)+1)*pha1,TWOPI)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood

def usage():
    sys.stderr.write("""
usage:  get_TOAs.py [options which must include -t or -g] pfd_file
  [-h, --help]                       : Display this help
  [-s numsub, --subbands=numsub]     : Divide the fold into numsub subbands
  [-n numTOAs, --numtoas=numTOAs]    : Divide the fold into numTOAs parts
  [-d DM, --dm=DM]                   : Re-combine subbands at DM
  [-f, --FFTFITouts]                 : Print all FFTFIT outputs and errors
  [-g gausswidth, --gaussian=width]  : Use a Gaussian template of FWHM width
                                       or, if the arg is a string, read the file
                                       to get multiple-gaussian parameters
  [-t templateprof, --template=prof] : The template .bestprof file to use
  [-k subs_list, --kill=subs_list]   : List of subbands to ignore
  [-i ints_list, --kints=ints_list]  : List of intervals to ignore
  [-o seconds, --offset=seconds]     : Add the offset in seconds to any TOAs
  [-e, --event]                      : The .pfd file was made with events
  [-r, --norotate]                   : Do not rotate the template for FFTFIT
  [-2, --tempo2]                     : Write Tempo2 format TOAs
  pfd_file                           : The .pfd file containing the folds

  The program generates TOAs from a .pfd file using Joe Taylor's
  FFTFIT program.  The TOAs are output to STDOUT.  Typically, the .pfd
  file is created using prepfold with the "-timing" flag and an
  appropriate .par file on either a topocentric time series or raw
  telescope data.  But barycentric folds or folds of barycentered
  events are also acceptable.  The number of bins in the folded profile
  must be a power of two for FFTFIT to work.  The most important thing 
  about the fold, though, is that it must have been made using "-nosearch"! 
  (Note: "-timing" implies "-nosearch" and forces a power-of-two number 
  of bins.)
  
  A typical example would be something like:
      
      get_TOAs.py -n 30 -t myprof.bestprof -k 0,20-23 myprof.pfd | \\
          tail -28 >> good.tim
      
  which would extract 30 TOAs (the default number of slices or parts
  in time for "prepfold -timing" is 60) from a fold made from some raw
  radio telescope data.  The command would ignore (i.e. zero-out)
  subbands 0, 20, 21, 22, and 23 (e.g.  due to interference) and then
  ignore the first 2 TOAs with the tail command.
  
  If you don't specify "-n", the default number of parts in the fold
  is assumed, but if you don't specify "-s", all the subbands (if any
  are present) are integrated together.
  
  If you specify the "-f" flag, an additional line of output is
  displayed for each TOA that shows the "b +/- berr" and "SNR +/-
  SNRerr" params from FFTFIT.
  
""")

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "herf2s:n:d:g:t:o:k:i:",
                                   ["help", "event", "norotate", "FFTFITouts",
                                    "tempo2","subbands=", "numtoas=", "dm=", "gaussian=",
                                    "template=", "offset=", "kill=", "kints="])
                                    
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    lowfreq = None
    DM = 0.0
    gaussianwidth = 0.1
    gaussfitfile = None
    templatefilenm = None
    rotate_prof = True
    numsubbands = 1
    numtoas = 1
    otherouts = 0
    offset = 0.0
    events = 0
    t2format = False
    kill = []
    kints = []
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--FFTFITouts"):
            otherouts = 1
        if o in ("-r", "--norotate"):
            rotate_prof = False
        if o in ("-2", "--tempo2"):
            t2format = True
        if o in ("-e", "--event"):
            lowfreq = 0.0
            DM = 0.0
            events = 1
        if o in ("-s", "--subbands"):
            numsubbands = int(a)
        if o in ("-n", "--numtoas"):
            numtoas = int(a)
            if numtoas==0:
                sys.exit()
        if o in ("-d", "--dm"):
            DM = float(a)
        if o in ("-g", "--gaussian"):
            try:
                gaussianwidth = float(a)
            except ValueError:
                gaussfitfile = a
        if o in ("-t", "--template"):
            templatefilenm = a
        if o in ("-o", "--offset"):
            offset = float(a)
        if o in ("-k", "--kill"):
            for subs in a.split(','):
                if (subs.find("-") > 0):
                    lo, hi = subs.split("-")
                    kill.extend(list(range(int(lo), int(hi)+1)))
                else:
                    kill.append(int(subs))
        if o in ("-i", "--kints"):
            for ints in a.split(','):
                if (ints.find("-") > 0):
                    lo, hi = ints.split("-")
                    kints.extend(list(range(int(lo), int(hi)+1)))
                else:
                    kints.append(int(ints))

    # Read the prepfold output file and the binary profiles
    fold_pfd = pfd(sys.argv[-1])

    # Check to make sure we can use this .pfd for timing purposes
    if not fold_pfd.use_for_timing():
        sys.stderr.write(
            "Error: '%s' was made allowing prepfold to search!\n" % \
            sys.argv[-1])
        sys.exit(2)
    
    # Read key information from the bestprof file
    if fold_pfd.bestprof:
        fold = fold_pfd.bestprof
    else:
        sys.stderr.write(
            "Error:  Can't open '%s.bestrof'!  Regenerate with show_pfd.\n" % \
            sys.argv[-1])
        sys.exit(2)
    timestep_sec = fold.T / numtoas
    timestep_day = timestep_sec / SECPERDAY
    fold.epoch = fold.epochi+fold.epochf

    # If the requested number of TOAs doesn't divide into the
    # number of time intervals, then exit
    if fold_pfd.npart % numtoas:
        sys.stderr.write(
            "Error: # of TOAs (%d) doesn't divide # of time intervals (%d)!\n" % \
            (numtoas, fold_pfd.npart))
        sys.exit(2)

    # Over-ride the DM that was used during the fold
    if (DM!=0.0):
        fold_pfd.bestdm = DM
    if (fold_pfd.numchan==1 and DM==0.0 and events):
        fold_pfd.bestdm = 0.0
        fold_pfd.numchan = 1

    # Kill any required channels and/or subband
    fold_pfd.kill_subbands(kill)

    # Kill any required intervals
    fold_pfd.kill_intervals(kints)

    # De-disperse at the requested DM
    # Also save the pulse period used in dedispersion calc
    fold_pfd.dedisperse(interp=1)
    p_dedisp = fold_pfd.proflen / fold_pfd.binspersec
    
    # Combine the profiles as required
    profs = fold_pfd.combine_profs(numtoas, numsubbands)

    # PRESTO de-disperses at the high frequency channel so determine a
    # correction to the middle of the band
    if not events:
        subpersumsub = fold_pfd.nsub/numsubbands
        # Calculate the center of the summed subband freqs and delays
        sumsubfreqs = (Num.arange(numsubbands)+0.5)*subpersumsub*fold_pfd.subdeltafreq + \
                      (fold_pfd.lofreq-0.5*fold_pfd.chan_wid)
        # Note:  In the following, we cannot use fold_pfd.hifreqdelay since that
        #        is based on the _barycentric_ high frequency (if the barycentric 
        #        conversion was available).  For TOAs, we need a topocentric
        #        delay, which is based on the topocentric frequency fold_pfd.hifreq
        sumsubdelays = (psr_utils.delay_from_DM(fold_pfd.bestdm, sumsubfreqs) -
                        psr_utils.delay_from_DM(fold_pfd.bestdm, fold_pfd.hifreq))
        sumsubdelays_phs = Num.fmod(sumsubdelays / p_dedisp, 1.0)
        # Save the "higest channel within a subband" freqs/delays for use in
        # later DM/timing correction. PBD 2011/11/03
        sumsubfreqs_hi = sumsubfreqs + \
                fold_pfd.subdeltafreq/2.0 - fold_pfd.chan_wid/2.0
        subdelays2 = psr_utils.delay_from_DM(fold_pfd.bestdm, sumsubfreqs) - \
                psr_utils.delay_from_DM(fold_pfd.bestdm, sumsubfreqs_hi)

    else:
        fold_pfd.subfreqs = Num.asarray([0.0])
        sumsubfreqs = Num.asarray([0.0])
        sumsubdelays = Num.asarray([0.0])
        subdelays2 = Num.asarray([0.0])
        sumsubdelays_phs = Num.asarray([0.0])

    # Read the template profile
    if templatefilenm is not None:
        template = psr_utils.read_profile(templatefilenm, normalize=1)
    else:
        if (gaussfitfile):
            template = psr_utils.read_gaussfitfile(gaussfitfile, fold_pfd.proflen)
        else:
            template = psr_utils.gaussian_profile(fold_pfd.proflen, 0.0, gaussianwidth)
        template = template / max(template)
    #from Pgplot import *
    #plotxy(template)
    #closeplot()
    # Determine the Telescope used
    if (not fold.topo):
        obs = '@'  # Solarsystem Barycenter
    else:
        try: 
            if t2format:
                obs = scopes2[fold_pfd.telescope.split()[0]]
            else:
                obs = scopes[fold_pfd.telescope.split()[0]]
        except KeyError:  sys.stderr.write("Unknown telescope!!! : " + fold_pfd.telescope)

    # Read the polyco file (if required)
    if (fold.psr and fold.topo):
        if ("polycos" in fold_pfd.__dict__ and
            not fold_pfd.polycos==0):
            pcs = fold_pfd.polycos
        else:
            pcs = polycos(fold.psr, sys.argv[-1]+".polycos")
        (fold.phs0, fold.f0) = pcs.get_phs_and_freq(fold.epochi, fold.epochf)
        fold.f1 = fold.f2 = 0.0
    else:
        pcs = None
        fold.phs0 = 0.0
        (fold.f0, fold.f1, fold.f2) = psr_utils.p_to_f(fold.p0, fold.p1, fold.p2)

    #
    # Calculate the TOAs
    #

    if t2format:
        print("FORMAT 1")
        
    for ii in range(numtoas):

        # The .pfd file was generated using -nosearch and a specified
        # folding period, p-dot, and p-dotdot (or f, f-dot, and f-dotdot).
        if (pcs is None):
            # Time at the middle of the interval in question
            midtime = fold.epoch + (ii+0.5)*timestep_day
            p = 1.0/psr_utils.calc_freq(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
            t0 = psr_utils.calc_t0(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
            t0i= int(t0 + 1e-9)
            t0f = t0 - t0i
        # The .pfd file was folded using polycos
        else:
            # Time at the middle of the interval in question
            mjdf = fold.epochf + (ii+0.5)*timestep_day
            (phs, f0) = pcs.get_phs_and_freq(fold.epochi, mjdf)
            phs -= fold.phs0
            p = 1.0/f0
            if (phs < 0.0): phs += 1.0 # Consistent with pat
            t0f = mjdf - phs*p/SECPERDAY
            t0i = fold.epochi

        for jj in range(numsubbands):
            prof = profs[ii][jj]

            # If we have zapped intervals or subbands, or added padding
            # sometimes we can get folds with no signal at all.  Skip these.

            if Num.std(prof)==0.0:
                sys.stderr.write("Skipping TOA %d for subband %d due to lack of signal\n"%(ii+1, jj+1))
                continue

            # Make sure that the template and the data have the same number of bins
            if (not len(template)==fold_pfd.proflen):
                if (not ((len(template)%fold_pfd.proflen)==0 or
                         (fold_pfd.proflen%len(template))==0)):
                    if not ii and not jj:
                        sys.stderr.write("WARNING!: Lengths of template (%d) and data (%d) are incompatible!  Skipping '%s'!\n" % (len(template), fold_pfd.proflen, fold_pfd.filenm))
                    continue
                # Interpolate the data
                if (len(template) > fold_pfd.proflen):
                    prof = psr_utils.linear_interpolate(prof, len(template)//fold_pfd.proflen)
                    if not ii and not jj:
                        sys.stderr.write("Note: Interpolating the data for '%s'\n"%fold_pfd.filenm)
                # Interpolate the template
                elif (1):
                    template = psr_utils.linear_interpolate(template, fold_pfd.proflen//len(template))
                    if not ii and not jj:
                        sys.stderr.write("Note: Interpolating the template for '%s'\n"%fold_pfd.filenm)
                # Downsample the data (Probably not a good idea)
                else:
                    prof = psr_utils.downsample(prof, fold_pfd.proflen//len(template))
                    if not ii and not jj:
                        sys.stderr.write("Note:  Downsampling the data for '%s'\n"%fold_pfd.filenm)

            try:
                tau = None
                if len(prof) & 2*len(prof):
                    sys.stderr.write("Profile length %d is not a power of two; unable to use FFTFIT.\n" % len(prof))
                elif len(template) & 2*len(template):
                    sys.stderr.write("Template length %d is not a power of two; unable to use FFTFIT.\n" % len(template))
                else:
                    # Try using FFTFIT first
                    shift,eshift,snr,esnr,b,errb,ngood = measure_phase(prof, template, rotate_prof)
                    # tau and tau_err are the predicted phase of the pulse arrival
                    tau, tau_err = shift/len(prof), eshift/len(prof)
                    # Note: "error" flags are shift = 0.0 and eshift = 999.0

                    # If that failed, use a time-domain correlation
                    if (Num.fabs(shift) < 1e-7 and
                        Num.fabs(eshift-999.0) < 1e-7):
                        sys.stderr.write("Warning!  Bad return from FFTFIT. May be due to inadequate signal-to-noise.\n")
                        tau = None
                if tau is None:
                    sys.stderr.write("Warning: using PRESTO correlation - reported error is incorrect...\n")
                    # Not enough structure in the template profile for FFTFIT
                    # so use time-domain correlations instead
                    tau = psr_utils.measure_phase_corr(prof, template)
                    # This needs to be changed
                    tau_err = 0.1/len(prof)

                # Calculate correction for dedispersion to true channel
                # center freqs that used a slightly different pulse
                # period.
                dd_phs_2 = subdelays2[jj] * (1.0/p - 1.0/p_dedisp)

                # Sum up several phase shifts
                tau_tot = Num.fmod(tau+sumsubdelays_phs[jj]+dd_phs_2+3.0, 1.0)
                if (tau_tot > 0.5): tau_tot -= 1.0

                # Send the TOA to STDOUT
                toaf = t0f + (tau_tot*p + offset)/SECPERDAY
                newdays = int(Num.floor(toaf))
                if t2format:
                    psr_utils.write_tempo2_toa(t0i+newdays, toaf-newdays,
                                                  tau_err*p*1000000.0,
                                                  sumsubfreqs[jj], 0.0, name=fold_pfd.pfd_filename, obs=obs)
                else:
                    psr_utils.write_princeton_toa(t0i+newdays, toaf-newdays,
                                                  tau_err*p*1000000.0,
                                                  sumsubfreqs[jj], 0.0, obs=obs)
                if (otherouts):
                    sys.stderr.write("FFTFIT results:  b = %.4g +/- %.4g   SNR = %.4g +/- %.4g" %
                          (b, errb, snr, esnr))

            except ValueError as xxx_todo_changeme:
                fftfit.error = xxx_todo_changeme
                pass
