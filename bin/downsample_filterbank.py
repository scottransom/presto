#!/usr/bin/env python
import sys, sigproc
import numpy as num

if __name__ == "__main__":
    if len(sys.argv)==1:
        print "\nusage:  downsample_filterbank_hdr.py DS_fact infile.fil\n"
        sys.exit()
    DS_fact = int(sys.argv[1])
    basefilenm = sys.argv[2][:sys.argv[2].find(".fil")]
    
    filhdr = {}
    newhdr = ""
    infile = open(sys.argv[2], 'rb')
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

    # Now loop over the spectra
    nchans = filhdr['nchans']
    while (1):
        try:
            x = num.fromfile(infile, dtype=num.ubyte, count=DS_fact*nchans)
            x.shape = (DS_fact, nchans)
            dsx = (x.mean(0)+0.5).astype(num.ubyte)
            dsx.tofile(outfile)
        except:
            break

    infile.close()
    outfile.close()
