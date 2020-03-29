#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
import numpy as np
import matplotlib.pyplot as plt
import sys
from presto import psrfits
from optparse import OptionParser

full_usage = """
"""

usage = "Calculate the average and stdev bandpass of PSRFITS search data"

def write_bandpass(filenm, freqs, means, stdevs):
    of = open(filenm, 'w')
    of.write("# Chan   Freq(MHz)     Mean       StDev\n")
    for ii, (freq, mean, stdev) in enumerate(zip(freqs, means, stdevs)):
        of.write("%6d  %9.3f  %9.3f  %9.3f\n" % (ii, freq, mean, stdev))
    of.close()

def plot_bandpass(freqs, means, stdevs, outfile=None):
    plt.plot(freqs, means, '-k',
             freqs, means+stdevs, '-r',
             freqs, means-stdevs, '-r')
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Relative power or Counts")
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)

def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Show the bandpass in an x-window as well")
    parser.add_option("-p", "--plot", action="store_true", dest="plot",
                      default=False, help="Show the bandpass in a .png plot as well")
    parser.add_option("-n", "--nomods", action="store_true", dest="nomods",
                      default=False, help="Do not apply offsets/scales (default applies)")
    parser.add_option("-w", "--weights", action="store_true", dest="weights",
                      default=False, help="Apply weights (default doesn't apply_")
    parser.add_option("-f", "--first", type="int", dest="subfirst", default=0,
                      help="First subint to compute stats for")
    parser.add_option("-s", "--skip", type="int", dest="subskip", default=10,
                      help="Number of subints to skip during stats calculations")
    parser.add_option("-o", "--outfile", type="string", dest="outfile", default=None,
                      help="Output filename (default will be INFILE.bandpass")
    (opts, args) = parser.parse_args()
    if len(args)==0:
        print(full_usage)
        sys.exit(0)

    for infile in args:
        print("Processing '%s'" % (infile))
        pf = psrfits.PsrfitsFile(infile)
        if opts.nomods:
            # for a bandpass histogram of raw bits
            htot = np.zeros(1<<pf.nbits)
        subints = np.arange(opts.subfirst, pf.specinfo.num_subint,
                            opts.subskip).astype(np.int)
        means = np.zeros((len(subints), pf.nchan))
        stdevs = np.zeros((len(subints), pf.nchan))
        for ii, subint in enumerate(subints):
            print("%.0f%%.." % (100.0 * float(subint) / pf.specinfo.num_subint), end=' ')
            sys.stdout.flush()
            specs = pf.read_subint(subint, apply_weights=opts.weights,
                                   apply_scales=not opts.nomods,
                                   apply_offsets=not opts.nomods,
                                   apply_zero_off=not opts.nomods)
            if opts.nomods:
                h, b = np.histogram(specs.flatten(),
                                    bins=np.arange((1<<pf.nbits)+1))
                htot += h
            means[ii] = specs.mean(axis=0)
            stdevs[ii] = specs.std(axis=0)
        print("%.0f%%" % (100.0))
        med_mean = np.median(means, axis=0)
        med_stdev = np.median(stdevs, axis=0)
        outfilenm = infile+".bandpass" if opts.outfile is None else opts.outfile
        plotfilenm = outfilenm+".png" if opts.plot else None
        if opts.xwin or opts.plot:
            plot_bandpass(pf.freqs, med_mean, med_stdev, outfile=plotfilenm)
        if opts.nomods:
            htot = htot / htot.sum()
            print("# Bits histogram")
            print("# val  fract")
            print("#---------------")
            for b, h in zip(b, htot):
                print("%3d  %6.4f" % (b, h))
        write_bandpass(outfilenm, pf.freqs, med_mean, med_stdev)

if __name__=='__main__':
    main()
