#!/usr/bin/env python
import struct, getopt, sys
from umath import *
from Numeric import *
from fftfit import *
from infodata import infodata
from bestprof import bestprof
from prepfold import pfd
from polycos import polycos
from psr_constants import *
from types import StringType, FloatType, IntType

scopes = {'GBT':'1', 'Arecibo':'3', 'Parkes':'7', 'GMRT': 'r'}

def measure_phase(profile, template):
    """
    measure_phase(profile, template):
        Call FFTFIT on the profile and template to determine the
            following parameters: shift,eshift,snr,esnr,b,errb,ngood
            (returned as a tuple).  These are defined as in Taylor's
            talk at the Royal Society.
    """
    c,amp,pha = cprof(template)
    pha.savespace()
    pha1 = pha[0]
    pha = fmod(pha-arange(1,len(pha)+1)*pha1,TWOPI)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood

def usage():
    print """
usage:  get_TOAs.py [options which must include -t or -g] pfd_file
  [-h, --help]                       : Display this help
  [-s numsub, --subbands=numsub]     : Divide the fold into numsub subbands
  [-n numTOAs, --numtoas=numTOAs]    : Divide the fold into numTOAs parts
  [-d DM, --dm=DM]                   : Re-combine subbands at DM
  [-f, --FFTFITouts]                 : Print all FFTFIT outputs and errors
  [-g gausswidth, --gaussian=width]  : Use a Gaussian template of FWHM width
  [-t templateprof, --template=prof] : The template .bestprof file to use
  [-k subs_list, --kill=subs_list]   : List of subbands to ignore
  [-e, --event]                      : The .pfd file was made with events
  pfd_file                           : The .pfd file containing the folds

  The program generates TOAs from a .pfd file using Joe Taylor's
  FFTFIT program. The TOAs are output to STDOUT.  Typically, the .pfd
  file is created using prepfold with the "-timing" flag and an
  appropriate .par file on either a topocentric time series or raw
  telescope data.  But barycentric folds or folds of barycentered
  events are also acceptable.  The most important thing about the
  fold, though, is that it must have been made using "-nosearch"! 
  (Note: "-timing" implies "-nosearch")
  
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
  
"""

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hefs:n:d:g:t:o:k:e:",
                                   ["help", "event", "FFTFITouts", "subbands=", 
				    "numtoas=", "dm=", "gaussian=", "template=",
                                    "offset=", "kill="])
                                    
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
    templatefilenm = None
    numchannels = 1
    numsubbands = 1
    numtoas = 1
    otherouts = 0
    offset = 0.0
    events = 0
    kill = []
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--FFTFITouts"):
	    otherouts = 1
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
            gaussianwidth = float(a)
        if o in ("-t", "--template"):
            templatefilenm = a
        if o in ("-o", "--offset"):
            offset = float(a)
        if o in ("-k", "--kill"):
            for subs in a.split(','):
                if (subs.find("-") > 0):
                    lo, hi = subs.split("-")
                    kill.extend(range(int(lo), int(hi)+1))
                else:
                    kill.append(int(subs))

    # Read key information from the bestprof file
    fold = bestprof(sys.argv[-1]+".bestprof")

    # Read the prepfold output file and the binary profiles
    fold_pfd = pfd(sys.argv[-1])
    # Over-ride the DM that was used during the fold
    if (DM!=0.0):
        fold_pfd.bestdm = DM
    if (fold_pfd.numchan==1):
        if (DM==0.0 and numchannels==1):
            idata = infodata(fold_pfd.filenm[:fold_pfd.filenm.rfind('.')]+".inf")
	    if (events):
		fold_pfd.bestdm = 0.0
		fold_pfd.numchan = 1
	    else:
		fold_pfd.bestdm = idata.DM
		fold_pfd.numchan = idata.numchan
	else:
            fold_pfd.bestdm = DM
            fold_pfd.numchan = numchannels

    # Read the template profile
    if templatefilenm is not None:
        template_fold = bestprof(templatefilenm)
        template = template_fold.normalize()
    else:
        template = gaussian_profile(fold_pfd.proflen, 0.0, gaussianwidth)
        template = template / max(template)

    timestep_sec = fold.T / numtoas
    timestep_day = timestep_sec / SECPERDAY

    # PRESTO de-disperses at the high frequency channel so determine a
    # correction to the middle of the band
    if not events:
	binspersec = fold_pfd.fold_p1*fold_pfd.proflen
	chanpersub = fold_pfd.numchan/fold_pfd.nsub
	subdeltafreq = fold_pfd.chan_wid*chanpersub
	losubfreq = fold_pfd.lofreq + subdeltafreq - fold_pfd.chan_wid
	subfreqs = arange(fold_pfd.nsub, typecode='d')*subdeltafreq + losubfreq
	subdelays = delay_from_DM(fold_pfd.bestdm, subfreqs)
	hifreqdelay = subdelays[-1]
	subdelays = subdelays-hifreqdelay
	subdelays_bins = floor(subdelays*binspersec+0.5)
	subpersumsub = fold_pfd.nsub/numsubbands
	# Calculate the center of the summed subband freqs and delays
	sumsubfreqs = (arange(numsubbands)+0.5)*subpersumsub*subdeltafreq + \
                      (fold_pfd.lofreq-0.5*fold_pfd.chan_wid)
	sumsubdelays = (delay_from_DM(fold_pfd.bestdm, sumsubfreqs)-hifreqdelay)/SECPERDAY
    else:
	subdelays = asarray([0.0])
	subdelays_bins = asarray([0.0])
	sumsubdelays = asarray([0.0])
	subfreqs = asarray([0.0])
	sumsubfreqs = asarray([0.0])
	
    # Shift the profiles by the correct amount required by the best DM
    for ii in range(fold_pfd.npart):
        for jj in range(fold_pfd.nsub):
            if jj in kill:
                fold_pfd.profs[ii][jj] *= 0.0
            else:
                fold_pfd.profs[ii][jj] = rotate(fold_pfd.profs[ii][jj],
                                                int(subdelays_bins[jj]))
    #fold.epochf += DM_delay/SECPERDAY
    if fold.epochf > 1.0:
        fold.epochf -= 1.0
        fold.epochi += 1
    fold.epoch = fold.epochi+fold.epochf

    if (not fold.topo):
        obs = '@'  # Solarsystem Barycenter
    else:
        try: obs = scopes[fold_pfd.telescope.split()[0]]
	except KeyError:  print "Unknown telescope!!!"

    # Read the polyco file
    if (fold.psr and fold.topo):
        pcs = polycos(fold.psr, sys.argv[-1]+".polycos")
        (fold.phs0, fold.f0) = pcs.get_phs_and_freq(fold.epochi, fold.epochf)
        fold.f1 = fold.f2 = 0.0
    else:
        pcs = None
        fold.phs0 = 0.0
        (fold.f0, fold.f1, fold.f2) = p_to_f(fold.p0, fold.p1, fold.p2)

    # Combine the sub-integration profiles
    profs = zeros((numtoas, numsubbands, fold_pfd.proflen), 'd')
    dp = fold_pfd.npart/numtoas
    ds = fold_pfd.nsub/numsubbands
    for ii in range(numtoas):
        # Combine the subbands if required
        if (fold_pfd.nsub > 1):
            for jj in range(numsubbands):
                subprofs = add.reduce(fold_pfd.profs[:,jj*ds:(jj+1)*ds], 1)
                # Combine the time intervals
                profs[ii][jj] = add.reduce(subprofs[ii*dp:(ii+1)*dp])
        else:
            profs[ii][0] = add.reduce(fold_pfd.profs[ii*dp:(ii+1)*dp,0])

    #if (0):
    #    for ii in range(numtoas):
    #        profs[ii][0] = (profs[ii][0]-min(profs[ii][0]))/max(profs[ii][0])
    #    plot2d(profs[:,0], arange(64.0)/64.0, arange(numtoas, typecode='d'), image='antigrey')
    #    closeplot()

    # Calculate the TOAs
    for ii in range(numtoas):
        if (pcs is None):
            midtime = fold.epoch + (ii+0.5)*timestep_day
            p = 1.0/calc_freq(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
            t0 = calc_t0(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
        else:
            mjdf = fold.epochf + (ii+0.5)*timestep_day
            (phs, f0) = pcs.get_phs_and_freq(fold.epochi, mjdf)
            phs -= fold.phs0
            p = 1.0/fold.f0
            t0 = fold.epochi+mjdf - phs*p/SECPERDAY
        for jj in range(numsubbands):
            prof = profs[ii][jj]
            try:
		if (len(template)==len(prof)):
		    shift,eshift,snr,esnr,b,errb,ngood = measure_phase(prof, template)
		    tau, tau_err = shift/fold_pfd.proflen, eshift/fold_pfd.proflen
		else:
		    shift = 0.0
		    eshift = 999.0
		if (fabs(shift) < 1e-7 and
		    fabs(eshift-999.0) < 1e-7):
		    # Not enough structure in the template profile for FFTFIT
		    # so use time-domain correlations instead
		    tau = measure_phase_corr(prof, template)
		    # This needs to be changed
		    tau_err = 0.1/len(prof)
		write_princeton_toa(t0+(tau*p+offset)/SECPERDAY+sumsubdelays[jj],
                                    tau_err*p*1000000.0,
                                    sumsubfreqs[jj], fold_pfd.bestdm, obs=obs)
		if (otherouts):
		    print "FFTFIT results:  b = %.4g +/- %.4g   SNR = %.4g +/- %.4g" % \
		        (b, errb, snr, esnr)
	    except ValueError, fftfit.error:
                pass
