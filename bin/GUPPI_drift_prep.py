#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import sys, os, random
from astropy.io import fits
from presto import sigproc
from presto import psr_utils as pu
import math as math


def guppi_subint_per_file(guppi_filenm):
    """
    guppi_samples_per_file(spigot_filenm,):
        Return the number of subints present in the GUPPI FITs file.
    """
    fitsfile = fits.open(guppi_filenm,memmap=True)
    nsubint = fitsfile['SUBINT'].header['NAXIS2']
    fitsfile.close()
    return nsubint
 
def guppi_time_per_subint(guppi_filenm):
    fitsfile = fits.open(guppi_filenm,memmap=True)
    subint_hdr = fitsfile['SUBINT'].header
    time_subint = subint_hdr['TBIN']*subint_hdr['NSBLK']
    fitsfile.close()
    return time_subint

debug = 0 

if __name__=="__main__":
    if (len(sys.argv) < 3):
        print("usage:  GUPPI_drift_prep.py NUM guppi_fits_files")
        print("    NUM is the 'beam' number in the scan.  It starts ")
        print("        with 0 and goes to NMAX.  If NUM is < 0, NMAX")
        print("        is sent to STDOUT by the program.")
        sys.exit()

    new_obs_length = 141 #approximate number of seconds we want the new 
                         #observations to last

    infilenms = sys.argv[2:]
    print(len(infilenms))
    time_subint = guppi_time_per_subint(infilenms[0])


    orig_N  = int(math.floor(new_obs_length / time_subint)) # Number of subints to analyze at a time (~141 sec)
    print(orig_N)
    raw_N   = orig_N     # Number of subints to step through .fits files
    overlap_factor = 0.5  # Overlap each orig_N samples by this fraction
    overlap_subints = int(orig_N * overlap_factor)
    print(overlap_subints)
    nom_subint_per_file = 320

    # Now see how much data we have to work with
    subints_per_file = []
    numinfiles = len(infilenms)
    print(numinfiles)
    for ii in range(numinfiles):
        subints = guppi_subint_per_file(infilenms[ii])
        if ((subints < nom_subint_per_file) and (ii < numinfiles-1)):
            print("Warning!  '%s' only has %d samples!"%\
                  (infilenms[ii], subints))
            print("    You need to fix that file!")
            sys.exit(-1)
        subints_per_file.append(subints)
    total_subints = sum(subints_per_file)
    print(total_subints)
    num = int(sys.argv[1])
    nmax = total_subints//overlap_subints-1
    if num < 0:
        print(nmax)
        sys.exit(0)
    if num > nmax:
        print("NUM > NMAX (%d)!  Exiting!"%nmax)
        sys.exit(-1)

    # Now figure out which file is the first
    first_subint = num * overlap_subints
    accum_subints = 0
    for ii in range(len(subints_per_file)):
        next_accum_subints = accum_subints + subints_per_file[ii]
        if next_accum_subints > first_subint:
            first_filenm = infilenms[ii]
            # How much data to skip in the first file
            skip = first_subint - accum_subints
            # How many total files we need
            first_file_subints = subints_per_file[ii]-skip
            numfiles = (raw_N - first_file_subints) // nom_subint_per_file + 1
            if ((raw_N - first_file_subints) % nom_subint_per_file):
                numfiles += 1
            if debug:
                print("first_filenum  = ", ii)
                print("1st subint     = ", first_subint)
                print("1st filenam    = ", infilenms[ii])
                print("skip           = ", skip)
                print("1st_file_samps = ", first_file_subints)
                print("numfiles       = ", numfiles)
            break
        else:
            accum_subints += subints_per_file[ii]
    
    # Now make a command line option for guppidrift2fil.py
    tmpfilenm = "tmp%d.fil"%random.randint(0,2**30)
    cmd = "guppidrift2fil.py --skip=%d --nsubint=%d -o %s " % \
          (skip, raw_N, tmpfilenm)
    for goodfile in infilenms[ii:ii+numfiles]:
        cmd += "%s "%goodfile
        print(cmd)
    os.system(cmd)

    # Now read the header to determine what the correct filename
    # should be.  Use that to rename the fil file.

    filhdr, hdrlen = sigproc.read_header(tmpfilenm)
    MJDi = int(filhdr['tstart'])
    ra_rad = sigproc.ra2radians(filhdr['src_raj'])
    ra_string = pu.coord_to_string(*pu.rad_to_hms(ra_rad))
    dec_rad = sigproc.dec2radians(filhdr['src_dej'])
    dec_string = pu.coord_to_string(*pu.rad_to_dms(dec_rad))
    str_coords = "".join(ra_string.split(":")[:2])
    if dec_rad >= 0.0:  str_coords += "+"
    str_coords += "".join(dec_string.split(":")[:2])
    filfilenm = "GBT350drift_%d_%s.fil" % (MJDi, str_coords)
    os.rename(tmpfilenm, filfilenm)
    print("Renamed '%s' to '%s'." % (tmpfilenm, filfilenm))
    

    
