#!/usr/bin/env python

#Begun on 2/13/2014 from a copy of psrfits2fil.py //NEG-D//
from __future__ import print_function
from builtins import range
import numpy as np
from astropy.io import fits as pyfits
from presto import filterbank
import optparse
import sys
import os
import time
from math import *
# 2/13/2014 //NEG-D///
import pyslalib.slalib as sla
from presto import presto
fil_header_keys = [
    "telescope_id",
    "machine_id",
    "data_type", 
    "rawdatafile",
    "source_name", 
    "barycentric", 
    "pulsarcentric", 
    "az_start",  
    "za_start",  
    "src_raj",  
    "src_dej",  
    "tstart",  
    "tsamp",  
    "nbits", 
    "fch1",  
    "foff",
    "nchans", 
    "nifs" ]

telescope_ids = {"Fake": 0, "Arecibo": 1, "ARECIBO 305m": 1, "Ooty": 2, "Nancay": 3,
                 "Parkes": 4, "Jodrell": 5, "GBT": 6, "GMRT": 7,
                 "Effelsberg": 8, "ATA": 9, "UTR-2": 10, "LOFAR": 11}

machine_ids = {"FAKE": 0, "PSPM": 1, "Wapp": 2, "WAPP": 2, "AOFTM": 3,
               "BCPM1": 4, "OOTY": 5, "SCAMP": 6, "GBT Pulsar Spigot": 7, 
               "SPIGOT": 7, "BG/P": 11, "pdev": 11}

def read_4bit(data):                       
    """ 
    Unpack 4-bit PSRFITS data that has been read in as bytes
     by pyfits.

         Input: array of unsigned 8-bit ints
         Output: unpacked array
    """
    first_piece = np.bitwise_and(15,data)
    #second_piece = np.bitwise_and(240,data) / 16
    second_piece = data >> 4
    #return np.array([first_piece,second_piece]).transpose().flatten()
    return np.dstack([first_piece,second_piece]).flatten()

def read_subint(fits,i_subint,nchan,nsamps, apply_weights=True, \
                apply_scales=True, apply_offsets=True, input_nbits=4):
    """
    Read a 4-bitized PSRFITS subint from a open pyfits file object.
     Applys scales, weights, and offsets to the data.

         Input: fits - open pyfits file object
                i_subint - index of subint (first subint is 0)
                nchan - number of frequency channels
                nsamps - number of time samples per subint
                apply_weights - If True, apply weights. 
                        (Default: apply weights)
                apply_scales - If True, apply scales. 
                        (Default: apply scales)
                apply_offsets - If True, apply offsets. 
                        (Default: apply offsets)
         Output: subint data with scales, weights, and offsets
                 applied in float32 dtype with shape (nsamps,nchan).
    """ 

    if input_nbits == 4:
        data = read_4bit(fits['SUBINT'].data[i_subint]['DATA'])
    elif input_nbits > 4:
        data = fits['SUBINT'].data[i_subint]['DATA']
    if apply_weights:
        offsets = fits['SUBINT'].data[i_subint]['DAT_OFFS']
    else:
        offsets = 0
    if apply_scales:
        scales = fits['SUBINT'].data[i_subint]['DAT_SCL']
    else:
        scales = 1
    if apply_weights:
        weights = fits['SUBINT'].data[i_subint]['DAT_WTS']
    else:
        weights = 1
    data = data.reshape((nsamps,nchan))
    data_wso = ((data * scales) + offsets) * weights
    return data_wso

def translate_header(fits_file,skip,output_subints):
    fits_hdr = fits_file['PRIMARY'].header
    subint_hdr = fits_file['SUBINT'].header 
    fil_header = dict.fromkeys(fil_header_keys,None)

    if fits_hdr['TELESCOP'] in telescope_ids:
        fil_header["telescope_id"] = telescope_ids[fits_hdr['TELESCOP']]
    else:
        fil_header["telescope_id"] = -1
    if fits_hdr['BACKEND'] in machine_ids:
        fil_header["machine_id"] = machine_ids[fits_hdr['BACKEND']]
    else:
        fil_header["machine_id"] = -1

    fil_header["data_type"] = 1 # filterbank
    # Get filename in a way that is safe for old versions of pyfits
    # (i.e. using private attribute)

    #Let's get the center of our intended output file.  
    #Right now I'm assuming we're in search mode and that we'll take
    #The subint offset from the first file and add skip + 
    #output_subints / 2 subints to get the position of the telescope.
    time_subint = subint_hdr['TBIN']*subint_hdr['NSBLK']
    time_offset = subint_hdr['NSUBOFFS'] * time_subint + \
                  skip * time_subint 

    #Let's set some variables that we'll use to update our position
    #via pyslalib
    type = "A"         #We need this for SLA
    dtmp = 0.0         #Delta UT
    atm  = 1010.0      #Local atmospheric pressure in mB
    temp = 283.0       #Local temperature in DegK
    humid = 0.5        #Local relative humidity in range 0.0-1.0
    tlr = 0.0065       #Tropospheric lapse rate (DegL/meter)
    eq = 2000.0        #Julian epoch of mean place
    microns = 3e8 / fits_hdr['OBSFREQ']*1e6*1e6
    az = radians(fits_file['SUBINT'].data[0]['TEL_AZ']) #azimuth
    za = radians(fits_file['SUBINT'].data[0]['TEL_ZEN']) #zenith
    tstart = fits_hdr['STT_IMJD'] + \
                           fits_hdr['STT_SMJD']/86400.0 + \
                           fits_hdr['STT_OFFS']/86400.0 + \
                           time_offset/86400.0  #tstart of new file
    #Now let's get the MJD of the center of the new file 
    #so that we can update the RA/DEC for the new file

    MJD = tstart + (time_subint * output_subints / 2)/86400.0

    #Ok, I think we're set to update our position.  Let's use slalib
    #to get the position of the telescope 
  
    telecope,telescope_name,tel_lon,tel_lat,tel_hgt = \
      sla.sla_obs(0,fits_hdr['TELESCOP'])
    #Now we have tel_lon,tel_lat and tel_ght.  We need to flip the lon

    tel_lon = -tel_lon
    
    #Let's get the geocentric apparent RA and DEC
    rap,dap = sla.sla_oap(type,az,za,MJD,dtmp,tel_lon,tel_lat,tel_hgt,\
                      dtmp,dtmp,temp,atm,humid,microns,tlr)

    #OK, for reals we are almost there  Let's convert to FK5 (mean place)

    rmn,dmn = sla.sla_amp(rap,dap,MJD,eq)
    hours,hmin,hsec = presto.hours2hms(degrees(rmn)/15.0)
    deg,dmin,dsec = presto.deg2dms(degrees(dmn))


    fil_header["src_raj"] = float(hours*10000.0 + hmin*100.0 +  hsec)
    fil_header["src_dej"] = float(deg*10000.0 + dmin*100.0 + dsec)

    fn = fits_file.filename()
    fil_header["rawdatafile"] = os.path.basename(fn) 
    fil_header["source_name"] = fits_hdr['SRC_NAME']
    fil_header["barycentric"] = 0 # always not barycentered?
    fil_header["pulsarcentric"] = 0 # whats pulsarcentric?
    fil_header["az_start"] = fits_file['SUBINT'].data[0]['TEL_AZ']
    fil_header["za_start"] = fits_file['SUBINT'].data[0]['TEL_ZEN']
#   fil_header["src_raj"] = float(fits_hdr['RA'].replace(':',''))
#   fil_header["src_dej"] = float(fits_hdr['DEC'].replace(':',''))
    fil_header["tstart"] = fits_hdr['STT_IMJD'] + \
                           fits_hdr['STT_SMJD']/86400.0 + \
                           fits_hdr['STT_OFFS']/86400.0 + \
                           time_offset/86400.0
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
            apply_weights, apply_scales, apply_offsets,
            output_subints, skip):
    start = time.time()


    firstfits = pyfits.open(fits_fn[0],memmap=True)

    nchan = firstfits['SUBINT'].header['NCHAN']
    nsamps = firstfits['SUBINT'].header['NSBLK']
    nsubints = firstfits['SUBINT'].header['NAXIS2']




    fil_header = translate_header(firstfits,skip,output_subints) 
    fil_header['nbits'] = nbits
    outfil = filterbank.create_filterbank_file(outfn, fil_header, \
                                        nbits=nbits)


    # if frequency channels are in ascending order
    # band will need to be flipped
    if firstfits['SUBINT'].header['CHAN_BW'] > 0:
        flip_band=True
        print("\nFits file frequencies in ascending order.")
        print("\tFlipping frequency band.\n")
    else:
        flip_band=False

    # check nbits for input
    input_nbits = firstfits['SUBINT'].header['NBITS']
    if input_nbits < 4:
        raise ValueError('Does not support %d-bit data' % input_nbits)

    if nbits != 32:
        print("\nCalculating statistics on first subintegration...")
        subint0 = read_subint(firstfits,0,nchan,nsamps, \
                        apply_weights, apply_scales, apply_offsets, \
                        input_nbits=input_nbits)
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

    firstfits.close()

    fits = [pyfits.open(filename,memmap=True) for filename in fits_fn]

    print("Writing data...")
    sys.stdout.flush()
    oldpcnt = ""
    for i in range(skip+1,output_subints+skip+1):
        index = (i-1)//320
        subint_in_file = i-1-(index * 320)
        subint = read_subint(fits[index],subint_in_file,nchan,nsamps, \
                    apply_weights, apply_scales, apply_offsets, \
                    input_nbits=input_nbits)
        if flip_band:
            subint = np.fliplr(subint)
        subint /= scale_fac
        outfil.append_spectra(subint)
        pcnt = "%d" % (i*100.0/output_subints)
        if pcnt != oldpcnt:
           sys.stdout.write("% 4s%% complete\r" % pcnt)
           sys.stdout.flush()

    print("Done               ")
    outfil.close()

    print("Runtime:",time.time() - start)

if __name__=='__main__':
    parser = optparse.OptionParser(prog='guppidrift2fil.py', \
                    version="v0.1 Nate Garver-Daniels (Mar 2014)", \
                    usage = "usage: %prog [options] input_fits files to read")
    parser.add_option("-n",dest='nbits', action='store',
                      default=8, type='int',
                      help="The number of bits in the output .fil file. " +\
                           "Default=8")
    parser.add_option("-o",dest='outfn',action='store',
                      default=None, type='string',  
                      help="The filename of the output filterbank file. " +\
                           "Default: same as .fits input but with .fil extn")
    parser.add_option("--noweights", dest='apply_weights', \
                      default=True, action="store_false", \
                      help="Do not apply weights when converting data.")
    parser.add_option("--noscales", dest='apply_scales', \
                      default=True, action="store_false", \
                      help="Do not apply scales when converting data.")
    parser.add_option("--nooffsets", dest='apply_offsets', \
                      default=True, action="store_false", \
                      help="Do not apply offsets when converting data.")
    parser.add_option("--nsubint", dest='number_subints',\
                      default="320", type='int',\
                      help="Specify the number of subints to write.")
    parser.add_option("--skip", dest='skip',action='store', \
                      default="0", type='int',\
                      help="Specify number of seconds to skip.")
    (options, args) = parser.parse_args()
    
    fits_fn = args

    print(fits_fn)

    if options.outfn:
        outfn = options.outfn
    else:
        parser.print_help()
        parser.error("You must supply a valid output file with -o")
 
    main(fits_fn, outfn, options.nbits, options.apply_weights, \
            options.apply_scales, options.apply_offsets, \
            options.number_subints, options.skip)
