#!/usr/bin/env python
from __future__ import print_function
import astropy.coordinates as coords
from builtins import str
import os
import os.path
import argparse
import sys
import getpass
import numpy as np
from presto import sigproc

BLOCKSIZE = 10000  # Amount of data to copy at a time
                   # from input file to output file (in samples)

def write_inf_file(datfn, hdr, hdrlen):
    """Write a PRESTO .inf file given a .dat file and
        a dictionary of SIGPROC-style header values.

        Inputs:
            datfn: The PRESTO .dat file to write a .inf file for.
            hdr: A dictionary of SIGPROC header values, as produced
                by PRESTO's sigproc.read_header.
            hdrlen: Length (in bytes) of SIGPROC file's header.

        Output:
            inffn: The corresponding .inf file that is created.
    """
    if not datfn.endswith(".dat"):
        raise ValueError("Was expecting a file name ending with '.dat'. "
                         "Got: %s" % datfn)
    size = os.path.getsize(datfn)
    if size % 4:
        raise ValueError("Bad size (%d bytes) for PRESTO .dat file (%s)"
                         "Should be multiple of 4 because samples are "
                         "32-bit floats." % (size, datfn))
    N = size / 4  # Number of samples
    pos = coords.SkyCoord(sigproc.ra2radians(hdr['src_raj']),
                          sigproc.dec2radians(hdr['src_dej']),
                          frame='icrs', unit='rad')
    rastr, decstr = pos.to_string('hmsdms', sep=':',
                                  precision=4, pad=True).split()
    inffn = datfn[:-4]+".inf"
    with open(inffn, 'w') as ff:
        ff.write(" Data file name without suffix          =  %s\n" %
                 os.path.basename(datfn))
        ff.write(" Telescope used                         =  %s\n" %
                 sigproc.ids_to_telescope[hdr['telescope_id']])
        ff.write(" Instrument used                        =  %s\n" %
                 sigproc.ids_to_machine.get('machine_id', 'UNKNOWN'))
        ff.write(" Object being observed                  =  %s\n" %
                 hdr['source_name'])
        ff.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" %
                 rastr)
        ff.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" %
                 decstr)
        ff.write(" Data observed by                       =  UNKNOWN\n")
        ff.write(" Epoch of observation (MJD)             =  %05.15f\n" %
                 hdr['tstart'])
        ff.write(" Barycentered?           (1=yes, 0=no)  =  %d\n" %
                 hdr['barycentric'])
        ff.write(" Number of bins in the time series      =  %d\n" % N)
        ff.write(" Width of each time series bin (sec)    =  %.15g\n" %
                 hdr['tsamp'])
        ff.write(" Any breaks in the data? (1 yes, 0 no)  =  0\n")
        if hdr.has_key('pulsarcentric'):
            ff.write(" Orbit removed?          (1=yes, 0=no)  =  %d\n" %
                     hdr['pulsarcentric'])
        ff.write(" Dispersion measure (cm-3 pc)           =  %f\n" %
                 hdr['refdm'])
        ff.write(" Central freq of low channel (Mhz)      =  %f\n" %
                 hdr['fch1'])
        if hdr.has_key('foff'):
            ff.write(" Total bandwidth (Mhz)                  =  %f\n" %
                     (hdr['nchans']*hdr['foff']))
        else: # what else can we do?
            ff.write(" Total bandwidth (Mhz)                  =  %f\n" %
                     100.0)
        ff.write(" Number of channels                     =  %d\n" %
                 hdr['nchans'])
        if hdr.has_key('foff'):
            ff.write(" Channel bandwidth (Mhz)                =  %d\n" %
                     hdr['foff'])
        else: # what else can we do?
            ff.write(" Channel bandwidth (Mhz)                =  %d\n" %
                     100.0)
        ff.write(" Data analyzed by                       =  %s\n" %
                 getpass.getuser())
        ff.write(" Any additional notes:\n"
                 "    File converted from SIGPROC .tim time series\n"
                 "    with PRESTO's tim2dat.py, written by Patrick Lazarus\n")
    return inffn


def convert_tim_to_dat(tim):
    """Convert a SIGPROC time series .tim file to a
        PRESTO .dat time series

        Input:
            tim: The SIGPROC .tim time series file to convert.

        Output:
            datfn: The PRESTO .dat time series file
    """
    if not tim.endswith(".tim"):
        raise ValueError("Was expecting a file name ending with '.tim'. "
                         "Got: %s" % tim)
    path, fn = os.path.split(tim)
    basenm = fn[:-4]
    outfn = os.path.join(path, basenm+".dat")
    hdr, hdrlen = sigproc.read_header(tim)
    N = sigproc.samples_per_file(tim, hdr, hdrlen)
    Ndone = 0
    status = -1
    with open(tim, 'rb') as inff, open(outfn, 'wb') as outff:
        inff.seek(hdrlen)
        data = np.fromfile(inff, dtype='float32', count=BLOCKSIZE)
        while data.size:
            data.tofile(outff)
            Ndone += data.size
            data = np.fromfile(inff, dtype='float32', count=BLOCKSIZE)
            newstatus = int(100.0*Ndone/N)
            if newstatus > status:
                sys.stdout.write(" %d %%\r" % newstatus)
                sys.stdout.flush()
                status = newstatus
    return outfn


def main():
    for tim in args.timfiles:
        print("Working on %s" % tim)

        if args.write_dat:
            try:
                datfn = convert_tim_to_dat(tim)
                print("    Wrote PRESTO time series: %s" % datfn)
            except ValueError as e:
                sys.stderr.write("Error encountered when converting on %s" % tim)
                sys.stderr.write(str(e))
        else:
            datfn = tim[-3:]+"dat"
        hdr, hdrlen = sigproc.read_header(tim)
        inffn = write_inf_file(datfn, hdr, hdrlen)
        print("    Wrote info data: %s" % inffn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("timfiles", nargs="+",
                        help="SIGPROC .tim time series files to convert to "
                             "PRESTO *.dat time series files")
    parser.add_argument("--inf-only", dest='write_dat', action='store_false',
                        help="Only produce the .inf file, not the .dat file."
                             "(Default: create both .dat and .inf files)")
    args = parser.parse_args()
    main()
