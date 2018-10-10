#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import sys, os, random
from presto import sigproc
from presto import psr_utils as pu


def spigot_samples_per_file(spigot_filenm):
    """
    spigot_samples_per_file(spigot_filenm,):
        Return the number of samples present in the Spigot FITs file.
    """
    hdrlen = 184320
    bytes_per_sample = 2048
    filelen = os.stat(spigot_filenm)[6]
    return int((filelen-hdrlen)/bytes_per_sample)

debug = 1

if __name__=="__main__":
    if (len(sys.argv) < 3):
        print("usage:  GBT350_drift_prep.py NUM spigot_fits_files")
        print("    NUM is the 'beam' number in the scan.  It starts ")
        print("        with 0 and goes to NMAX.  If NUM is < 0, NMAX")
        print("        is sent to STDOUT by the program.")
        sys.exit()

    orig_N  = 1728000     # Number of samples to analyze at a time (~141 sec)
    raw_N   = 1900000     # Number of samples to step through .fits files
    overlap_factor = 0.5  # Overlap each orig_N samples by this fraction
    overlap_samples = int(orig_N * overlap_factor)
    nom_samps_per_file = 976896

    # Now see how much data we have to work with
    samples_per_file = []
    infilenms = sys.argv[2:]
    numinfiles = len(infilenms)
    for ii in range(numinfiles):
        samps = spigot_samples_per_file(infilenms[ii])
        if ((samps < nom_samps_per_file) and (ii < numinfiles-1)):
            print("Warning!  '%s' only has %d samples!"%\
                  (infilenms[ii], samps))
            print("    You need to fix that file!")
            sys.exit(-1)
        samples_per_file.append(samps)
    total_samples = sum(samples_per_file)
    num = int(sys.argv[1])
    nmax = total_samples/overlap_samples-1
    if num < 0:
        print(nmax)
        sys.exit(0)
    if num > nmax:
        print("NUM > NMAX (%d)!  Exiting!"%nmax)
        sys.exit(-1)

    # Now figure out which file is the first
    first_sample = num * overlap_samples
    accum_samples = 0
    for ii in range(len(samples_per_file)):
        next_accum_samples = accum_samples + samples_per_file[ii]
        if next_accum_samples > first_sample:
            first_filenm = infilenms[ii]
            # How much data to skip in the first file
            skip = first_sample - accum_samples
            # How many total files we need
            first_file_samples = samples_per_file[ii]-skip
            numfiles = (raw_N - first_file_samples) / nom_samps_per_file + 1
            if ((raw_N - first_file_samples) % nom_samps_per_file):
                numfiles += 1
            if debug:
                print("first_filenum  = ", ii)
                print("1st sample     = ", first_sample)
                print("1st filenam    = ", infilenms[ii])
                print("skip           = ", skip)
                print("1st_file_samps = ", first_file_samples)
                print("numfiles       = ", numfiles)
            break
        else:
            accum_samples += samples_per_file[ii]
    
    # Now make a command line option for spigot2filterbank
    tmpfilenm = "tmp%d.fil"%random.randint(0,2**30)
    cmd = "spigot2filterbank -skip %d -numout %d -o %s " % \
          (skip, raw_N, tmpfilenm)
    for goodfile in infilenms[ii:ii+numfiles]:
        cmd += "%s "%goodfile
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
    

    
