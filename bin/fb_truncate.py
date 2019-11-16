#!/usr/bin/env python

"""
A script to truncate a filterbank file in time/frequency.

Patrick Lazarus, Aug 27, 2012
"""
from __future__ import print_function

import sys
import copy
from argparse import ArgumentParser
import numpy as np
from presto import filterbank

BLOCKSIZE = 1e5 # Number of spectra to manipulate at once



def main(args):
    infn = args[0]
    print("Reading filterbank file (%s)" % infn)
    fil = filterbank.FilterbankFile(infn)
    if options.start_time is None:
        startbin = 0
    else:
        startbin = int(np.round(options.start_time/fil.tsamp))
    if options.end_time is None:
        endbin = fil.nspec
    else:
        endbin = int(np.round(options.end_time/fil.tsamp))+1

    new_nspec = endbin-startbin
    if new_nspec <= 0:
        raise ValueError("Bad number of spectra to be written (%d). " \
                            "Check start/end times." % new_nspec)
   
    # Determine lo/hi channels to write to file
    # If high frequencies come first in spectra 'hichan' refers to 
    # the lo-freq cutoff and 'lochan' refers to the hi-freq cutoff.
    if options.lo_freq is None:
        if fil.foff > 0:
            lochan = 0
        else:
            hichan = fil.nchans
    else:
        ichan = int(np.round((options.lo_freq-fil.fch1)/fil.foff))
        if fil.foff > 0:
            lochan = ichan
        else:
            hichan = ichan+1
    if options.hi_freq is None:
        if fil.foff > 0:
            hichan = fil.nchans
        else:
            lochan = 0
    else:
        ichan = int(np.round((options.hi_freq-fil.fch1)/fil.foff))
        if fil.foff > 0:
            hichan = ichan+1
        else:
            lochan = ichan

    new_nchans = hichan-lochan
    if new_nchans <= 0:
        raise ValueError("Bad number of channels to be written (%d). " \
                            "Check lo/hi frequencies." % new_nchans)

    print("Will extract")
    print("    %d bins (%d to %d incl.)" % (new_nspec, startbin, endbin-1))
    print("    (Original num bins: %d)" % fil.nspec)
    print("    %d channels (%d to %d incl.)" % (new_nchans, lochan, hichan-1))
    print("    (Original num chans: %d)" % fil.nchans)

    # Create output file
    outfn = options.outname % fil.header
    print("Creating out file: %s" % outfn)
    outhdr = copy.deepcopy(fil.header)
    outhdr['nchans'] = new_nchans
    outhdr['fch1'] = fil.frequencies[lochan]
    filterbank.create_filterbank_file(outfn, outhdr, nbits=fil.nbits)
    outfil = filterbank.FilterbankFile(outfn, mode='write')

    # Write data
    sys.stdout.write(" %3.0f %%\r" % 0)
    sys.stdout.flush()
    nblocks = int(new_nspec/options.block_size)
    remainder = new_nspec % options.block_size
    oldprogress = -1
    for iblock in np.arange(nblocks):
        lobin = iblock*options.block_size + startbin
        hibin = lobin+options.block_size
        spectra = fil.get_spectra(lobin, hibin)
        spectra = spectra[:,lochan:hichan] # restrict channels
        outfil.append_spectra(spectra)
        progress = int(100.0*((hibin-startbin)/new_nspec))
        if progress > oldprogress: 
            sys.stdout.write(" %3.0f %%\r" % progress)
            sys.stdout.flush()
            oldprogress = progress
    # Read all remaining spectra
    if remainder:
        spectra = fil.get_spectra(endbin-remainder, endbin)
        spectra = spectra[:,lochan:hichan] # restrict channels
        outfil.append_spectra(spectra)
    sys.stdout.write("Done   \n")
    sys.stdout.flush()


if __name__ == '__main__':
    parser = ArgumentParser(description="v0.1 Patrick Lazarus (Aug. 28, 2012)")
    parser.add_argument("-L", "--lo-freq", dest="lo_freq", type=float,
                    help="Desired low frequency for output file. Note: "
                        "actual low frequency will be rounded to the nearest"
                        "channel (Default: Don't truncate low-freq channels)",
                    default=None)
    parser.add_argument("-H", "--hi-freq", dest="hi_freq", type=float,
                    help="Desired high frequency for output file. Note: "
                        "actual high frequency will be rounded to the nearest"
                        "channel (Default: Don't truncate high-freq channels)",
                    default=None)
    parser.add_argument("-s", "--start-time", dest="start_time", type=float,
                    help="Start of desired range of input file to write "
                        "to output file. Note: The actual start time will "
                        "be rounded to the nearest sample.(Default: Don't "
                        "truncate from start of file.)", default=None)
    parser.add_argument("-e", "--end-time", dest="end_time", type=float,
                    help="End of desired range of input file to write "
                        "to output file. Note: The actual end time will "
                        "be rounded to the nearest sample. (Default: "
                        "Don't truncate from end of file.)", default=None)
    parser.add_argument("--block-size", dest='block_size', default=BLOCKSIZE,
                    type=float,
                    help="Number of spectra per block. This is the amount "
                        "of data manipulated/written at a time. (Default: "
                        " %d spectra)" % BLOCKSIZE)
    parser.add_argument("-o", "--outname", dest='outname', action='store', required=True,
                    help="The name of the output file.")
    (options, args) = parser.parse_args()

    main(args)
