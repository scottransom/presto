#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
import os
import struct
import sys
import math
import warnings
from presto.psr_constants import ARCSECTORAD

telescope_ids = {"Fake": 0, "Arecibo": 1, "ARECIBO 305m": 1, 
                 "Ooty": 2, "Nancay": 3, "Parkes": 4, "Jodrell": 5,
                 "GBT": 6, "GMRT": 7, "Effelsberg": 8, "ATA": 9,
                 "SRT": 10, "LOFAR": 11, "VLA": 12, "CHIME": 20,
                 "FAST": 21, "MeerKAT": 64, "KAT-7": 65}
ids_to_telescope = dict(list(zip(list(telescope_ids.values()), list(telescope_ids.keys()))))

machine_ids = {"FAKE": 0, "PSPM": 1, "Wapp": 2, "WAPP": 2, "AOFTM": 3,
               "BCPM1": 4, "BPP": 4, "OOTY": 5, "SCAMP": 6,
               "GBT Pulsar Spigot": 7, "SPIGOT": 7, "BG/P": 11,
               "PDEV": 12, "CHIME+PSR": 20, "KAT": 64, "KAT-DC2": 65}
ids_to_machine = dict(list(zip(list(machine_ids.values()), list(machine_ids.keys()))))

header_params = {
    "HEADER_START": 'flag',
    "telescope_id": 'i',
    "machine_id": 'i',
    "data_type": 'i', 
    "rawdatafile": 'str',
    "source_name": 'str', 
    "barycentric": 'i', 
    "pulsarcentric": 'i', 
    "az_start": 'd',  
    "za_start": 'd',  
    "src_raj": 'd',  
    "src_dej": 'd',  
    "tstart": 'd',  
    "tsamp": 'd',  
    "nbits": 'i', 
    "signed": 'b', 
    "nsamples": 'i', 
    "nbeams": "i",
    "ibeam": "i",
    "fch1": 'd',  
    "foff": 'd',
    "FREQUENCY_START": 'flag',
    "fchannel": 'd',  
    "FREQUENCY_END": 'flag',
    "nchans": 'i', 
    "nifs": 'i', 
    "refdm": 'd',  
    "period": 'd',  
    "npuls": 'q',
    "nbins": 'i', 
    "HEADER_END": 'flag'}

def dec2radians(src_dej):
    """
    dec2radians(src_dej):
       Convert the SIGPROC-style DDMMSS.SSSS declination to radians
    """
    sign = 1.0
    if (src_dej < 0): sign = -1.0;
    xx = math.fabs(src_dej)
    dd = int(math.floor(xx / 10000.0))
    mm = int(math.floor((xx - dd * 10000.0) / 100.0))
    ss = xx - dd * 10000.0 - mm * 100.0
    return sign * ARCSECTORAD * (60.0 * (60.0 * dd + mm) + ss)

def ra2radians(src_raj):
    """
    ra2radians(src_raj):
       Convert the SIGPROC-style HHMMSS.SSSS right ascension to radians
    """
    return 15.0 * dec2radians(src_raj)

def read_doubleval(filfile, stdout=False):
    dblval = struct.unpack('d', filfile.read(8))[0]
    if stdout:
        print("  double value = '%20.15f'"%dblval)
    return dblval

def read_intval(filfile, stdout=False):
    intval = struct.unpack('i', filfile.read(4))[0]
    if stdout:
        print("  int value = '%d'"%intval)
    return intval

def read_charval(filfile, stdout=False):
    charval = struct.unpack('b', filfile.read(1))[0]
    if stdout:
        print(" char value = '%d'"%charval)
    return charval

def read_longintval(filfile, stdout=False):
    longintval = struct.unpack('q', filfile.read(8))[0]
    if stdout:
        print("  long int value = '%d'"%longintval)
    return longintval

def read_string(filfile, stdout=False):
    strlen = struct.unpack('i', filfile.read(4))[0]
    strval = filfile.read(strlen)
    if stdout:
        print("  string = '%s'"%strval)
    return strval.decode('utf-8')

def read_paramname(filfile, stdout=False):
    paramname = read_string(filfile, stdout=False)
    if stdout:
        print("Read '%s'"%paramname)
    return paramname

def read_hdr_val(filfile, stdout=False):
    paramname = read_paramname(filfile, stdout)
    try:
        if header_params[paramname] == 'd':
            return paramname, read_doubleval(filfile, stdout)
        elif header_params[paramname] == 'i':
            return paramname, read_intval(filfile, stdout)
        elif header_params[paramname] == 'q':
            return paramname, read_longintval(filfile, stdout)
        elif header_params[paramname] == 'b':
            return paramname, read_charval(filfile, stdout)
        elif header_params[paramname] == 'str':
            return paramname, read_string(filfile, stdout)
        elif header_params[paramname] == 'flag':
            return paramname, None
    except KeyError:
        warnings.warn("key '%s' is unknown!" % paramname)
        return None, None

def prep_string(string):
    return struct.pack('i', len(string))+string.encode('utf-8')

def prep_double(name, value):
    return prep_string(name)+struct.pack('d', float(value))

def prep_int(name, value):
    return prep_string(name)+struct.pack('i', int(value))

def prep_char(name, value):
    return prep_string(name)+struct.pack('b', int(value))

def addto_hdr(paramname, value):
    try:
        if header_params[paramname] == 'd':
            return prep_double(paramname, value)
        elif header_params[paramname] == 'i':
            return prep_int(paramname, value)
        elif header_params[paramname] == 'b':
            return prep_char(paramname, value)
        elif header_params[paramname] == 'str':
            return prep_string(paramname) + prep_string(value)
        elif header_params[paramname] == 'flag':
            return prep_string(paramname)
    except KeyError:
        warnings.warn("key '%s' is unknown!" % paramname)
        return ""

def read_header(infile):
    """
    read_header(infile):
       Read a SIGPROC-style header and return the keys/values in a dictionary,
          as well as the length of the header: (hdrdict, hdrlen)
    """
    hdrdict = {}
    if type(infile) == type("abc"):
        infile = open(infile,'rb')
    param = ""
    while (param != "HEADER_END"):
        param, val = read_hdr_val(infile, stdout=False)
        hdrdict[param] = val
    hdrlen = infile.tell()
    infile.close()
    return hdrdict, hdrlen

def samples_per_file(infile, hdrdict, hdrlen):
    """
    samples_per_file(infile, hdrdict, hdrlen):
       Given an input SIGPROC-style filterbank file and a header
           dictionary and length (as returned by read_header()),
           return the number of (time-domain) samples in the file.
    """
    numbytes = os.stat(infile)[6] - hdrlen
    bytes_per_sample = hdrdict['nchans'] * (hdrdict['nbits']/8)
    if numbytes % bytes_per_sample:
        print("Warning!:  File does not appear to be of the correct length!")
    numsamples = numbytes / bytes_per_sample
    return numsamples

if __name__ == "__main__":
    if len(sys.argv)==1:
        print("\nusage:  mod_filterbank_hdr.py infile.fil [outfile.fil]\n")
        sys.exit()
    filhdr = {}
    newhdr = ""
    infile = open(sys.argv[1], 'rb')

    # Loop over the values in the .fil file
    while 1:
        param, val = read_hdr_val(infile, stdout=True)
        filhdr[param] = val

        # Add lines here to correct stuff
        #if param=="nchans":  val = 768

        # Append to the new hdr string
        # newhdr += addto_hdr(param, val)

        # Break out of the loop if the header is over
        if param=="HEADER_END":  break
            
    if len(sys.argv) > 2:
        print("Writing new header to '%s'"%sys.argv[2])
        outfile = open(sys.argv[2], 'wb')
        outfile.write(newhdr)
        outfile.close()
    else:
        print(filhdr)
