"""
A module for reading filterbank files.

Patrick Lazarus, June 26, 2012
(Minor modification from file originally from June 6th, 2009)
"""

import sys
import warnings
import os.path
import numpy as np
import sigproc


DEBUG = False

class Filterbank(object):
    def __init__(self, header, data):
        self.header = header
        self.data = data
        self.number_of_samples = self.data.size / self.nchans
        self.spectra = self.view_as_spectra(self.data)
        self.frequencies = self.fch1 + self.foff*np.arange(self.nchans)
        self.is_hifreq_first = (self.foff < 0)

    def view_as_spectra(self, data):
        """Return a view with data as spectra.

            Input:
                data: A 1-D numpy array of filterbank data.

            Output:
                spectra: A view of the input data arranged as spectra.
        """
        spectra = data.view()
        spectra.shape = (data.size/self.nchans, self.nchans)
        spectra = spectra.transpose()
        return spectra


    def __getattr__(self, name):
        if DEBUG:
            print "Fetching header param (%s)" % name
        return self.header[name]

    def print_header(self):
        """Print header parameters and values.
        """
        for param in sorted(self.header.keys()):
            if param in ("HEADER_START", "HEADER_END"):
                continue
            print "%s: %s" % (param, self.header[param])

    def write_filterbank_file(self, outfn):
        """Write filterbank header and data to file.
 
            Input:
                outfn: The outfile filterbank file's name.
 
            Output:
                None
        """
        outfile = open(outfn, 'wb')
        outfile.write(sigproc.addto_hdr("HEADER_START", None))
        for paramname in self.header.keys():
            if DEBUG:
                print "Writing header param (%s)" % paramname
            value = self.header[paramname]
            outfile.write(sigproc.addto_hdr(paramname, value))
        outfile.write(sigproc.addto_hdr("HEADER_END", None))
        self.data.astype(self.dtype).tofile(outfile)
        outfile.close()
 

class FilterbankFile(Filterbank):
    def __init__(self, filfn, delayed_read=False):
        if not os.path.isfile(filfn):
            raise ValueError("ERROR: File does not exist!\n\t(%s)" % filfn)
        else:
            self.filename = filfn
            self.filfile = open(filfn, 'r+b')
            header = self.read_header()
            self.header_size = self.filfile.tell()
            self.data_size = os.stat(self.filename)[6] - self.header_size
            bytes_per_sample = header['nchans']*header['nbits'] / 8
            # Calculate additional information
            # Such as: datatype, numsamps, datasize, hdrsize
            if header['nbits'] not in [32, 16, 8]:
                raise NotImplementedError("'filterbank.py' only supports " \
                                            "files with 8- or 16-bit " \
                                            "integers, or 32-bit floats " \
                                            "(nbits provided: %g)!" % \
                                            header['nbits'])
            if header['nbits'] == 32:
                self.dtype = 'float32'
            else:
                self.dtype = 'uint%d' % header['nbits']
            if self.data_size % bytes_per_sample:
                warnings.warn("Not an integer number of samples in file.")
            if delayed_read:
                data = np.array([])
            else:
                data = self.read_all_samples()
            super(FilterbankFile, self).__init__(header, data)

    def close(self):
        self.filfile.close()

    def __del__(self):
        self.close()

    def read_header(self):
        header = {}
        self.filfile.seek(0)
        paramname = ""
        while (paramname != 'HEADER_END'):
            if DEBUG:
                print "File location: %d" % self.filfile.tell()
            paramname, val = sigproc.read_hdr_val(self.filfile, stdout=DEBUG)
            if DEBUG:
                print "Read param %s (value: %s)" % (paramname, val)
            if paramname not in ["HEADER_START", "HEADER_END"]:
                header[paramname] = val
        return header

    def read_sample(self):
        """
        Read one sample starting at current
        position. Return as numpy array.
        
        NOTE: No checks are made to see if 
              current position is start of
              a sample.
        """
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.nchans)

    def read_all_samples(self):
        """
        Read all samples from file.
        Return as numpy array.
        """
        self.seek_to_data_start()
        return np.fromfile(self.filfile, dtype=self.dtype)

    def read_Nsamples(self, N):
        """
        Read N samples starting at current
        position. Return as numpy array.
        
        NOTE: No checks are made to see if 
              current position is start of
              a sample.
        """
        return np.fromfile(self.filfile, dtype=self.dtype, count=self.nchans*N)
    

    def seek_to_data_start(self):
        self.filfile.seek(self.header_size)

    def seek_to_sample(self, sampnum):
        """
        Seek to sample 'sampnum'. First sample is
        sampnum=0.
        """
        self.filfile.seek(self.header_size + self.nbits/8*self.nchans*sampnum)
        
    def seek_to_position(self, posn):
        """
        See to position 'posn' relative to
        beginning of file.
        """
        self.filfile.seek(posn)

    def append_spectra(self, spectra):
        """Append spectra to the filterbank file.
            
            Input:
                spectra: An array of spectra.

            Outputs:
                None
        """
        if spectra.shape[0] != self.nchans:
            raise ValueError("Bad number of channels (%d), expected (%d)" % \
                            (spectra.shape[0], self.nchans))
        spectra.T.flatten().astype(self.dtype).tofile(self.filfile)
        self.data_size += int(spectra.size*self.nbits/8.0)
        self.number_of_samples += spectra.shape[1]

def main():
    fil = FilterbankFile(sys.argv[1])
    fil.print_header()

if __name__ == '__main__':
    main()
