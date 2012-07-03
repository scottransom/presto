#/usr/bin/env python

"""Inject a fake pulsar into real data, creating
a filterbank file.

Patrick Lazarus, June 26, 2012
"""
import sys
import optparse
import warnings
import copy

import numpy as np
import matplotlib.pyplot as plt

import filterbank
import psr_utils

NUMSECS = 1.0 # Number of seconds of data to use to determine global scale
              # when repacking floating-point data into integers


def parse_gaussians(gaussian_strs):
    if not gaussian_strs:
        warnings.warn("Using default Gaussian profile (Amplitude=1.0 " \
                        "FWHM=0.05, and phase=0.0)")
        gaussians = [(1.0, 0.05, 0.0)]
    else:
        gaussians = []
        for gaussian_str in gaussian_strs:
            split = gaussian_str.split()
            if len(split) != 3:
                raise ValueError("Bad number of gaussian components " \
                        "should be 3, got %d" % len(split))
            amp = float(split[0])
            fwhm = float(split[1])
            phs = float(split[2])

            if fwhm <= 0 or fwhm > 1:
                raise ValueError("FWHM provided (%g) is bad, should be " \
                        "between 0 and 1" % fwhm)
            if phs < 0 or fwhm > 1:
                raise ValueError("Phase provided (%g) is bad, should be " \
                        "between 0 and 1" % phs)
            gaussians.append((amp,fwhm,phs))
    return gaussians


def main():
    print "%d input files provided" % len(args)
    for fn in args:
        print "Injecting pulsar signal into: %s" % fn
        fil = filterbank.FilterbankFile(fn, delayed_read=True)
        nbin = int(np.round(options.period/fil.tsamp)) # Number of bins 
                                                       # across the profile
        delays = psr_utils.delay_from_DM(options.dm, fil.frequencies)
        delay_bins = np.round(delays/fil.tsamp).astype(int)
    
        gaussians = parse_gaussians(options.gaussians)
        prof = np.zeros(nbin)
        for amp, fwhm, phs in gaussians:
            prof += amp*psr_utils.gaussian_profile(nbin, phs, fwhm)

        if options.dryrun:
            plt.figure()
            plt.plot(prof, 'k-')
            plt.show()
            sys.exit()

        profs = np.empty((fil.nchans, nbin))
        for ichan, bindelay in enumerate(delay_bins):
            profs[ichan,:] = psr_utils.rotate(prof, -bindelay)

        # Read the first second of data to get the global scaling to use
        fil.seek_to_data_start()
        onesec_nspec = int(np.round(1.0/fil.tsamp))
        onesec_data = fil.view_as_spectra(fil.read_Nsamples(onesec_nspec))
        nprofs = onesec_nspec/nbin
        remainder = onesec_nspec % nbin
        for iprof in np.arange(nprofs):
            onesec_data[:, iprof*nbin:(iprof+1)*nbin] += profs
        onesec_data[:, -remainder:] += profs[:, :remainder]
        minimum = np.min(onesec_data)
        median = np.median(onesec_data)
        # Set median to 1/3 of dynamic range
        global_scale = (256.0/3.0) / median
        del onesec_data
        
        # Start an output file
        outfil = filterbank.Filterbank(copy.deepcopy(fil.header), np.array([]))
        outfil.dtype = 'uint8'
        outfil.header['source_name'] += "_fake"
        outfn = options.outname % fil.header 
        print "Creating out file: %s" % outfn
        outfil.write_filterbank_file(outfn) 
        outfil = filterbank.FilterbankFile(outfn)

        fil.seek_to_data_start()
        spectra = fil.view_as_spectra(fil.read_Nsamples(nbin))
        while spectra.size:
            spectra += profs[:, :spectra.shape[1]]
            spectra = np.clip((spectra-minimum)*global_scale, 0, 256)
            outfil.append_spectra(spectra)
            spectra = fil.view_as_spectra(fil.read_Nsamples(nbin))


if __name__ == '__main__':
    parser = optparse.OptionParser(prog='injectpsr.py', \
                    version="v0.1 Patrick Lazarus (June 26, 2012)")
    parser.add_option("--dm", dest='dm', action='store', type='float', \
                    help="The DM of the (fake) injected pulsar signal. " \
                        "(This argument is required.", \
                    default=None)
    parser.add_option("-p", "--period", dest='period', action='store', \
                    default=None, type='float', \
                    help="The period (in seconds) of the (fake) injected " \
                        "pulsar signal. (This argument is required.)")
    parser.add_option("-g", "--gaussian", dest='gaussians', action='append', \
                    help="A string of 3 parameters defining a gaussian " \
                        "component to be injected. Be sure to quote the " \
                        "3 parameters together. The params are: 'amplitude " \
                        "FWHM phase'. Amplitude is not related to SNR in " \
                        "any way, FWHM and phase should be between 0 and 1. " \
                        "(Default: if no compoments are provided " \
                        "a Gaussian with amplitude=1.0, FWHM=0.05, and " \
                        "phase=0.0 will be used.)", \
                    default=[])
    parser.add_option("-n", "--dryrun", dest="dryrun", action="store_true", \
                    help="Show the pulse profile to be injected and exit. " \
                        "(Default: do not show profile, inject it)")
    parser.add_option("-o", "--outname", dest='outname', action='store', \
                    default="injected.fil", \
                    help="The name of the output file.")
    (options, args) = parser.parse_args()
    if options.period is None or options.dm is None:
        raise ValueError("Both a period and a DM _must_ be provided!")
    main()
