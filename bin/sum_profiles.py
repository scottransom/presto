#!/usr/bin/env python
import struct, getopt, sys, umath, fftfit, psr_utils, os.path, sinc_interp, Pgplot
import Numeric as Num
from Scientific.Statistics import standardDeviation, average
from infodata import infodata
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
    c,amp,pha = fftfit.cprof(template)
    pha.savespace()
    pha1 = pha[0]
    pha = umath.fmod(pha-Num.arange(1,len(pha)+1)*pha1,TWOPI)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood

def read_profile(filenm):
    prof = []
    for line in file(filenm):
        if line.startswith("#"): continue
        else: prof.append(float(line.split()[-1]))
    return Num.asarray(prof)

def usage():
    print """
usage:  sum_profiles.py [options which must include -t or -g] profs_file
  [-h, --help]                          : Display this help
  [-b bkgd_cutoff, --background=cutoff] : Fractional cutoff for the background level
  [-d DM, --dm=DM]                      : Re-combine subbands at DM
  [-n N, --numbins=N]                   : The number of bins to use in the resulting profile
  [-g gausswidth, --gaussian=width]     : Use a Gaussian template of FWHM width
  [-t templateprof, --template=prof]    : The template .bestprof file to use
  [-o outputfilenm, --output=filenm]    : The output file to use for the summed profile

  This program reads in a list of *.pfd files from profs_file and then
  de-disperses each of these using the DM specified.  The de-dispersed
  profiles are fit against a template to determine an absolute offset,
  and are then co-added together to produce a 'master' profile.  Each profile
  is scaled so that the RMS level of the off-pulse region is equivalent.

  To-do:  -- add a .par option so that the profiles can be added based on
             an ephemeris.
          -- add ability to kill subbands or intervals in the profs_file?

"""

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:d:n:g:t:o:",
                                   ["help", "background=", "dm=",
                                    "numbins=", "gaussian=", "template="])
                                    
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    lowfreq = None
    DM = 0.0
    bkgd_cutoff = 0.1
    gaussianwidth = 0.1
    templatefilenm = None
    numbins = 128
    outfilenm = "sum_profiles.bestprof"
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-b", "--background"):
            bkgd_cutoff = float(a)
        if o in ("-d", "--dm"):
            DM = float(a)
        if o in ("-n", "--numbins"):
            numbins = int(a)
        if o in ("-g", "--gaussian"):
            gaussianwidth = float(a)
        if o in ("-t", "--template"):
            templatefilenm = a
        if o in ("-o", "--output"):
            outfilenm = a

    print "Creating a summed profile of length %d bins using DM = %f"%(numbins, DM)

    # Read the template profile or create an appropriate Gaussian
    if templatefilenm is not None:
        template = read_profile(templatefilenm)
        # Resample the template profile to have the correct number of bins (if required)
        if not len(template)==numbins:
            oldlen = len(template)
            template = sinc_interp.periodic_interp(template, numbins)[::oldlen]
    else:
        template = psr_utils.gaussian_profile(numbins, 0.0, gaussianwidth)
    # Normalize it
    template -= min(template)
    template /= max(template)

    if (0):
        Pgplot.plotxy(template)
        Pgplot.closeplot()

    # Determine the off-pulse bins
    offpulse_inds = Num.compress(template<=bkgd_cutoff, Num.arange(numbins))
    onpulse_inds = Num.compress(template>bkgd_cutoff, Num.arange(numbins))
    if (1):
        Pgplot.plotxy(template)
        Pgplot.plotxy([bkgd_cutoff, bkgd_cutoff], [0.0, numbins], color='red')
        Pgplot.closeplot()

    # Read the list of *.pfd files to process
    pfdfilenms = []
    for line in file(sys.argv[-1]):
        pfdfilenm = line.split()[0]
        if not pfdfilenm.startswith("#"):
            if os.path.exists(pfdfilenm):
                pfdfilenms.append(pfdfilenm)
            else:
                print "Can't find '%s'.  Skipping it."

    sumprof = Num.zeros(numbins, typecode='d')

    # Step through the profiles and determine the offsets
    for pfdfilenm in pfdfilenms:

        print "\n  Processing '%s'..."%pfdfilenm

        # Read the fold data and de-disperse at the requested DM
        current_pfd = pfd(pfdfilenm)
        current_pfd.dedisperse(DM)
        prof = current_pfd.sumprof

        # Resample the current profile to have the correct number of bins
        if not len(prof)==numbins:
            oldlen = len(prof)
            prof = sinc_interp.periodic_interp(prof, numbins)[::oldlen]

        # Determine the amount to rotate the profile using FFTFIT
        shift,eshift,snr,esnr,b,errb,ngood = measure_phase(prof, template)

        # Rotate the profile to match the template
        newprof = psr_utils.interp_rotate(prof, shift)

        # Now shift and scale the profile
        offpulse = Num.take(newprof, offpulse_inds)
        newprof -= average(offpulse)
        newprof /= standardDeviation(offpulse)

        if (0):
            Pgplot.plotxy(newprof)
            Pgplot.closeplot()

        # Add it to the summed profile
        sumprof += newprof

    # Now normalize, plot, and write the summed profile
    offpulse = Num.take(sumprof, offpulse_inds)
    sumprof -= average(offpulse)
    sumprof /= standardDeviation(offpulse)
    Pgplot.plotxy(sumprof, Num.arange(numbins),
                  labx="Pulse Phase", laby="Relative Flux")
    Pgplot.closeplot()

    print "\n Writing profile to '%s'..."%(outfilenm),

    outfile = file(outfilenm, "w")
    for ii, val in enumerate(sumprof):
        outfile.write("%04d  %20.15g\n"%(ii, val))
    outfile.close()
    
    print "Done\n"
    
