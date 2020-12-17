#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import sys
import numpy as num
from presto import sigproc


if __name__ == "__main__":
    if len(sys.argv)==1:
        print("\nusage:  downsample_filterbank_hdr.py DS_fact infile.fil\n")
        sys.exit()
    DS_fact = int(sys.argv[1])
    basefilenm = sys.argv[2][:sys.argv[2].find(".fil")]
    
    filhdr = {}
    newhdr = "".encode('utf-8')
    infile = open(sys.argv[2], 'rb')

    # Determine the full size of the file
    infile.seek(0, 2)
    infilelen = infile.tell()
    infile.seek(0, 0)
    
    outfile = open(basefilenm+"_DS%d.fil"%DS_fact, 'wb')

    # Loop over the values in the .fil file
    while 1:
        param, val = sigproc.read_hdr_val(infile, stdout=False)
        filhdr[param] = val

        if param=="tsamp":
            val *= DS_fact

        # Append to the new hdr string
        newhdr += sigproc.addto_hdr(param, val)

        # Break out of the loop if the header is over
        if param=="HEADER_END":
            break

    # Write the new header to the output file
    outfile.write(newhdr)

    nchans = filhdr['nchans']

    # Remove the header length from infilelen and then
    # determine how many spectra are in the file
    infilelen -= infile.tell()
    numspec = infilelen // nchans
    if infilelen % nchans:
        print("Whoops!  File length calculation is not right...")

    # Now loop over the spectra
    for ii in range(numspec // DS_fact):
        try:
            x = num.fromfile(infile, dtype=num.ubyte, count=DS_fact*nchans)
            x.shape = (DS_fact, nchans)
            dsx = (x.mean(0)+0.5).astype(num.ubyte)
            dsx.tofile(outfile)
        except:
            break

    infile.close()
    outfile.close()
