#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
from builtins import range
import getopt, sys
import os.path
from presto import fftfit
from presto import psr_utils
from presto import sinc_interp
from presto import Pgplot
import numpy as Num
from presto.prepfold import pfd
from presto.psr_constants import *


scopes = {'GBT':'1', 'Arecibo':'3', 'Parkes':'7', 'GMRT': 'r'}


def measure_phase(profile, template):
    """
    measure_phase(profile, template):
        Call FFTFIT on the profile and template to determine the
            following parameters: shift,eshift,snr,esnr,b,errb,ngood
            (returned as a tuple).  These are defined as in Taylor's
            talk at the Royal Society.
    """
    c,amp,pha = fftfit.cprof(template)
    pha1 = pha[0]
    pha = Num.fmod(pha-Num.arange(1,len(pha)+1)*pha1,TWOPI)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood

def parse_vals(valstring):
    """
    parse_vals(valstring):
       Return a list of integers that corresponds to each of the numbers
          in a string representation where '-' gives an inclusive range
          and ',' separates individual values or ranges.  For example:

          > parse_vals('5,8,10-13,17')

          would return:  [5, 8, 10, 11, 12, 13, 17]
    """
    if (len(valstring)==0 or
        (len(valstring)==1 and not valstring.isdigit())):
        return None
    vals = set()
    for xx in valstring.split(','):
        if (xx.find("-") > 0):
            lo, hi = xx.split("-")
            vals = vals.union(set(range(int(lo), int(hi)+1)))
        else:
            vals.add(int(xx))
    vals = list(vals)
    vals.sort()
    return vals


def usage():
    print("""
usage:  sum_profiles.py [options which must include -t or -g] profs_file
  [-h, --help]                          : Display this help
  [-b bkgd_cutoff, --background=cutoff] : Fractional cutoff for the background level
                                          or, if the arg is a string (i.e. containing
                                          ',' and/or '-'), use the bins specified (as
                                          for parse_vals()) as the background values
  [-f, --fitbaseline]                   : Fit a 3-rd order polynomial to the specified
                                          background values before determining the RMS
  [-p pulsebins, --pulsebins=pulsebins] : A 'parse_vals' string that specifies the bins
                                          to include when integrating flux.  The minimum
                                          value (plus 1/2 of the std-dev) is substracted.
                                          This is not usually needed, but is used when there
                                          are pulse artifacts or if you want the flux of
                                          only part of a profile.
  [-d DM, --dm=DM]                      : Re-combine subbands at DM
  [-n N, --numbins=N]                   : The number of bins to use in the resulting profile
  [-g gausswidth, --gaussian=width]     : Use a Gaussian template of FWHM width
                                          or, if the arg is a string, read the file
                                          to get multiple-gaussian parameters
  [-t templateprof, --template=prof]    : The template .bestprof file to use
  [-o outputfilenm, --output=filenm]    : The output file to use for the summed profile
  [-s SEFD, --sefd=SEFD]                : For rough flux calcs, the SEFD (i.e. Tsys/G)

  This program reads in a list of *.pfd files from profs_file and then
  de-disperses each of these using the DM specified.  The de-dispersed
  profiles are fit against a template to determine an absolute offset,
  and are then co-added together to produce a 'master' profile.  Each profile
  is scaled so that the RMS level of the off-pulse region is equivalent.

  To-do:  -- add a .par option so that the profiles can be added based on
             an ephemeris.

""")

if __name__ == '__main__':
    # Import Psyco if available
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hfb:p:d:n:g:t:o:s:",
                                   ["help", "fitbaselibe", "background=", "pulsebins=",
                                    "dm=", "numbins=", "gaussian=", "template=",
                                    "outputfilenm=", "sefd="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    if len(sys.argv)==1:
        usage()
        sys.exit(2)

    fitbaseline = False
    lowfreq = None
    DM = 0.0
    bkgd_cutoff = 0.1
    bkgd_vals = None
    gaussianwidth = 0.1
    gaussfitfile = None
    templatefilenm = None
    pulsebins = None
    numbins = 128
    SEFD = 0.0
    outfilenm = "sum_profiles.bestprof"
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--fitbaseline"):
            fitbaseline = True
        if o in ("-b", "--background"):
            if '-' in a or ',' in a:
                bkgd_vals = Num.asarray(parse_vals(a))
            else:
                try:
                    bkgd_cutoff = float(a)
                except ValueError:
                    bkgd_vals = Num.asarray(parse_vals(a))
        if o in ("-p", "--pulsebins"):
            pulsebins = Num.asarray(parse_vals(a))
        if o in ("-d", "--dm"):
            DM = float(a)
        if o in ("-n", "--numbins"):
            numbins = int(a)
        if o in ("-g", "--gaussian"):
            try:
                gaussianwidth = float(a)
            except ValueError:
                gaussfitfile = a
        if o in ("-t", "--template"):
            templatefilenm = a
        if o in ("-o", "--output"):
            outfilenm = a
        if o in ("-s", "--sefd"):
            SEFD = float(a)

    print("Creating a summed profile of length %d bins using DM = %f"%(numbins, DM))

    # Read the template profile or create an appropriate Gaussian
    if templatefilenm is not None:
        template = psr_utils.read_profile(templatefilenm)
        # Resample the template profile to have the correct number of bins (if required)
        if not len(template)==numbins:
            oldlen = len(template)
            template = sinc_interp.periodic_interp(template, numbins)[::oldlen]
    else:
        if gaussfitfile is not None:
            template = psr_utils.read_gaussfitfile(gaussfitfile, numbins)
        else:
            template = psr_utils.gaussian_profile(numbins, 0.0, gaussianwidth)
    # Normalize it
    template -= min(template)
    template /= max(template)
    # Rotate it so that it becomes a "true" template according to FFTFIT
    shift,eshift,snr,esnr,b,errb,ngood = measure_phase(template, template)
    template = psr_utils.fft_rotate(template, shift)
        
    # Determine the off-pulse bins
    if bkgd_vals is not None:
        Pgplot.plotxy(template, labx="Phase bins")
        Pgplot.plotxy(template[bkgd_vals], Num.arange(numbins)[bkgd_vals],
                      line=None, symbol=2, color='red')
        Pgplot.closeplot()
        offpulse_inds = bkgd_vals
        onpulse_inds = set(Num.arange(numbins)) - set(bkgd_vals)
    else:
        offpulse_inds = Num.compress(template<=bkgd_cutoff, Num.arange(numbins))
        onpulse_inds = Num.compress(template>bkgd_cutoff, Num.arange(numbins))
        Pgplot.plotxy(template)
        Pgplot.plotxy([bkgd_cutoff, bkgd_cutoff], [0.0, numbins], color='red')
        Pgplot.closeplot()
    # If the number of bins in the offpulse section is < 10% of the total
    # use the statistics in the .pfd file to set the RMS
    if (len(offpulse_inds) < 0.1*numbins):
        print("Number of off-pulse bins to use for RMS is too low.  Using .pfd stats.")
        usestats = 1
    else:
        usestats = 0    

    # Read the list of *.pfd files to process
    pfdfilenms = []
    killsubss = []
    killintss = []
    for line in open(sys.argv[-1]):
        if not line.startswith("#"):
            sline = line.split()
            pfdfilenm = sline[0]
            if len(sline)==1:
                killsubs, killints = None, None
            elif len(sline)==2:
                killsubs = parse_vals(sline[1])
                killints = None
            elif len(sline)>=3:
                killsubs = parse_vals(sline[1])
                killints = parse_vals(sline[2])
            if os.path.exists(pfdfilenm):
                pfdfilenms.append(pfdfilenm)
                killsubss.append(killsubs)
                killintss.append(killints)
            else:
                print("Can't find '%s'.  Skipping it."%pfdfilenm)

    sumprof = Num.zeros(numbins, dtype='d')

    base_T = None
    base_BW = None
    orig_fctr = None
    Tprerfi = 0.0
    Tpostrfi = 0.0
    avg_S = 0.0

    # Step through the profiles and determine the offsets
    for pfdfilenm, killsubs, killints in zip(pfdfilenms, killsubss, killintss):

        print("\n  Processing '%s'..."%pfdfilenm)

        # Read the fold data and de-disperse at the requested DM
        current_pfd = pfd(pfdfilenm)
        current_pfd.dedisperse(DM)
        # This corrects for any searching that prepfold did to find the peak
        current_pfd.adjust_period()
        T = current_pfd.T
        Tprerfi += T
        BW = current_pfd.nsub*current_pfd.subdeltafreq
        fctr = current_pfd.lofreq + 0.5*BW

        # If there are subbands to kill, kill em'
        if killsubs is not None:
            print("    killing subbands:  ", killsubs)
            current_pfd.kill_subbands(killsubs)
            BW *= (current_pfd.nsub-len(killsubs))/float(current_pfd.nsub)
        # If there are intervals to kill, kill em'
        if killints is not None:
            print("    killing intervals: ", killints)
            current_pfd.kill_intervals(killints)
            T *= (current_pfd.npart-len(killints))/float(current_pfd.npart)

        if base_T is None:
            base_T = T
        if base_BW is None:
            base_BW = BW
        if orig_fctr is None:
            orig_fctr = fctr
        else:
            if fctr != orig_fctr:
                print("Warning!:  fctr = %f, but original f_ctr = %f!" % (fctr, orig_fctr))
        Tpostrfi += T

        prof = current_pfd.profs.sum(0).sum(0)

        # Resample the current profile to have the correct number of bins
        if not len(prof)==numbins:
            oldlen = len(prof)
            prof = sinc_interp.periodic_interp(prof, numbins)[::oldlen]

        # Determine the amount to rotate the profile using FFTFIT
        shift,eshift,snr,esnr,b,errb,ngood = measure_phase(prof, template)

        # Rotate the profile to match the template
        newprof = psr_utils.fft_rotate(prof, shift)
        offpulse = newprof[offpulse_inds]

        # Remove a polynomial fit from the off-pulse region if required
        if fitbaseline:
            pfit = Num.poly1d(Num.polyfit(offpulse_inds, offpulse, 3))
            offpulse -= pfit(offpulse_inds)
            if 0:
                Pgplot.plotxy(offpulse)
                Pgplot.closeplot()

        # Determine the off-pulse RMS
        if usestats:
            print("Using raw data statistics instead of off-pulse region")
            offpulse_rms = Num.sqrt(current_pfd.varprof)
        else:
            offpulse_rms = offpulse.std()
            
        Ppsr = 1.0 / current_pfd.fold_p1      # Pulsar period
        tau_bin = Ppsr / current_pfd.proflen  # Duration of profile bin
        dt_per_bin = tau_bin / current_pfd.dt
        corr_rms = offpulse_rms / Num.sqrt(current_pfd.DOF_corr())
        print("samples/bin = ", current_pfd.dt_per_bin)
        print("RMSs  (uncorr, corr)  = ", offpulse_rms, corr_rms)

        # Now attempt to shift and scale the profile so that it has
        # an off-pulse mean of ~0 and an off-pulse RMS of ~1
        offset = Num.median(newprof[offpulse_inds])
        newprof -= offset
        newprof /= corr_rms
        if 0:
            Pgplot.plotxy(newprof, labx="Phase bins")
            if fitbaseline:
                Pgplot.plotxy((pfit(offpulse_inds)-offset)/corr_rms, offpulse_inds,
                              color='yellow')
            Pgplot.plotxy(newprof[offpulse_inds], offpulse_inds,
                          line=None, symbol=2, color='red')
            if pulsebins is not None:
                Pgplot.plotxy(newprof[pulsebins], pulsebins,
                              line=None, symbol=2, color='blue')
            Pgplot.closeplot()

        if pulsebins is None:
            SNR = newprof.sum()  # integrate everything
        else:
            SNR = newprof[pulsebins].sum()
        print("    Approx SNR = %.3f" % SNR)
        if SEFD:
            S = SEFD * SNR / Num.sqrt(2.0 * BW * T / numbins) / numbins
            avg_S += S
            print("    Approx flux density = %.3f mJy" % S)
        
        # Now weight the profile based on the observation duration
        # and BW as compared to the first profile 
        newprof *= Num.sqrt(T/base_T * BW/base_BW)
        
        if (0):
            Pgplot.plotxy(newprof)
            Pgplot.closeplot()

        # Add it to the summed profile
        sumprof += newprof

    # Now normalize, plot, and write the summed profile
    offpulse = sumprof[offpulse_inds]
    # Remove a polynomial fit from the off-pulse region if required
    if fitbaseline:
        pfit = Num.poly1d(Num.polyfit(offpulse_inds, offpulse, 3))
        offpulse -= pfit(offpulse_inds)
        Pgplot.plotxy(offpulse)
        Pgplot.closeplot()
    # Now attempt to shift and scale the profile so that it has
    # an off-pulse mean of ~0 and an off-pulse RMS of ~1
    sumprof -= Num.median(offpulse)
    sumprof *= Num.sqrt(current_pfd.DOF_corr()) / offpulse.std()
    print("\nSummed profile approx SNR = %.3f" % sum(sumprof))
    if SEFD:
        avg_S /= len(pfdfilenms)
        if pulsebins is None:
            SNR = sumprof.sum()  # integrate everything
        else:
            SNR = sumprof[pulsebins].sum()
        S = SEFD * SNR / Num.sqrt(2.0 * BW * Tpostrfi / numbins) / numbins
        print("     Approx sum profile flux density = %.3f mJy" % S)
        print("    Avg of individual flux densities = %.3f mJy" % avg_S)
        print("     Total (RFI cleaned) integration = %.0f s (%.2f hrs)" % \
              (Tpostrfi, Tpostrfi/3600.0))

    # Rotate the summed profile so that the max value is at the phase ~ 0.25 mark
    sumprof = psr_utils.rotate(sumprof, -len(sumprof)/4)
    Pgplot.plotxy(sumprof, Num.arange(numbins),
                  labx="Pulse Phase", laby="Relative Flux")
    Pgplot.closeplot()

    print("\n Writing profile to '%s'..."%(outfilenm), end=' ')

    outfile = open(outfilenm, "w")
    for ii, val in enumerate(sumprof):
        outfile.write("%04d  %20.15g\n"%(ii, val))
    outfile.close()
    
    print("Done\n")
    
