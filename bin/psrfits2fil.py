#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import numpy as np
from presto import psrfits
from presto import filterbank
from presto import sigproc
import optparse
import sys
import os
import time

def translate_header(psrfits_file):
    fits_hdr = psrfits_file.header
    subint_hdr = psrfits_file.fits['SUBINT'].header 
    subint_data = psrfits_file.fits['SUBINT'].data
    fil_header = {} 

    if fits_hdr['TELESCOP'] in sigproc.telescope_ids:
        fil_header["telescope_id"] = \
                    sigproc.telescope_ids[fits_hdr['TELESCOP']]
    else:
        fil_header["telescope_id"] = -1
    if fits_hdr['BACKEND'] in sigproc.machine_ids:
        fil_header["machine_id"] = \
                    sigproc.machine_ids[fits_hdr['BACKEND']]
    else:
        fil_header["machine_id"] = -1

    fil_header["data_type"] = 1 # filterbank
    fn = psrfits_file.filename
    fil_header["rawdatafile"] = os.path.basename(fn)
    fil_header["source_name"] = fits_hdr['SRC_NAME']
    fil_header["barycentric"] = 0 # always not barycentered?
    fil_header["pulsarcentric"] = 0 # whats pulsarcentric?
    fil_header["az_start"] = subint_data[0]['TEL_AZ']
    fil_header["za_start"] = subint_data[0]['TEL_ZEN']
    fil_header["src_raj"] = float(fits_hdr['RA'].replace(':',''))
    fil_header["src_dej"] = float(fits_hdr['DEC'].replace(':',''))
    fil_header["tstart"] = fits_hdr['STT_IMJD'] + \
                           fits_hdr['STT_SMJD']/86400.0 + \
                           fits_hdr['STT_OFFS']/86400.0
    fil_header["tsamp"] = subint_hdr['TBIN']
    fil_header["nbits"] = None # set by user. Input should always be 4-bit.

    # first channel (fch1) in sigproc is the highest freq
    # foff is negative to signify this
    fil_header["fch1"] = fits_hdr['OBSFREQ'] + \
                         np.abs(fits_hdr['OBSBW'])/2.0 - \
                         np.abs(subint_hdr['CHAN_BW'])/2.0
    fil_header["foff"] = -1.0*np.abs(subint_hdr['CHAN_BW'])
    fil_header["nchans"] = subint_hdr['NCHAN']
    fil_header["nifs"] = subint_hdr['NPOL']

    return fil_header

def main(fits_fn, outfn, nbits, \
            apply_weights, apply_scales, apply_offsets):
    start = time.time()
    psrfits_file = psrfits.PsrfitsFile(fits_fn)

    fil_header = translate_header(psrfits_file)
    fil_header['nbits'] = nbits
    outfil = filterbank.create_filterbank_file(outfn, fil_header, \
                                        nbits=nbits)

    # if frequency channels are in ascending order
    # band will need to be flipped
    if psrfits_file.fits['SUBINT'].header['CHAN_BW'] > 0:
        flip_band=True
        print("\nFits file frequencies in ascending order.")
        print("\tFlipping frequency band.\n")
    else:
        flip_band=False

    # check nbits for input
    if psrfits_file.nbits < 4:
        raise ValueError('Does not support %d-bit data' % \
                        psrfits_file.nbits)

    if nbits != 32:
        print("\nCalculating statistics on first subintegration...")
        subint0 = psrfits_file.read_subint(0, \
                        apply_weights, apply_scales, apply_offsets)
        #new_max = np.mean(subint0) + 3*np.std(subint0)
        new_max = 3 * np.median(subint0)
        print("\t3*median =",new_max)
        if new_max > 2.0**nbits:
            scale = True
            scale_fac = new_max / ( 2.0**nbits )
            print("\tScaling data by",1/scale_fac)
            print("\tValues larger than",new_max,"(pre-scaling) "\
                  "will be set to",2.0**nbits - 1,"\n")

        else:
            scale = False
            scale_fac = 1
            print("\tNo scaling necessary")
            print("\tValues larger than",2.0**nbits-1,"(2^nbits) will "\
                  "be set to ",2.0**nbits-1,"\n")
    else:
        scale_fac = 1
        print("\nNo scaling necessary for 32-bit float output file.")

    print("Writing data...")
    sys.stdout.flush()
    oldpcnt = ""
    for isub in range(int(psrfits_file.nsubints)):
        subint = psrfits_file.read_subint(isub, apply_weights, apply_scales, apply_offsets)
        if flip_band:
            subint = np.fliplr(subint)
        subint /= scale_fac
        outfil.append_spectra(subint)
        pcnt = "%d" % (isub*100.0/psrfits_file.nsubints)
        if pcnt != oldpcnt:
            sys.stdout.write("% 4s%% complete\r" % pcnt)
            sys.stdout.flush()

    print("Done               ")
    outfil.close()

    print("Runtime:",time.time() - start)


if __name__=='__main__':
    parser = optparse.OptionParser(prog='psrfits2fil.py',
                    version="v0.2 Paul Scholz, Patrick Lazarus (Sept 2012)",
                    usage = "usage: %prog [options] input_fits")
    parser.add_option("-n",dest='nbits', action='store',
                      default=8, type='int',
                      help="The number of bits in the output .fil file. " +
                           "Default=8")
    parser.add_option("-o",dest='outfn',action='store',
                      default=None, type='string',
                      help="The filename of the output filterbank file. " +
                           "Default: same as .fits input but with .fil extn")
    parser.add_option("--noweights", dest='apply_weights',
                      default=True, action="store_false",
                      help="Do not apply weights when converting data.")
    parser.add_option("--noscales", dest='apply_scales',
                      default=True, action="store_false",
                      help="Do not apply scales when converting data.")
    parser.add_option("--nooffsets", dest='apply_offsets',
                      default=True, action="store_false",
                      help="Do not apply offsets when converting data.")
    (options, args) = parser.parse_args()

    fits_fn = args[0]

    if options.outfn:
        outfn = options.outfn
    else:
        outfn = '.'.join(fits_fn.split('.')[:-1]) + '.fil'

    main(fits_fn, outfn, options.nbits, options.apply_weights,
            options.apply_scales, options.apply_offsets)
