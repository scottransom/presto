"""
A module for reading filterbank files.

Patrick Lazarus, June 26, 2012
(Minor modification from file originally from June 6th, 2009)
"""

import sys
import warnings
import os
import os.path
import numpy as np
import sigproc


DEBUG = False

def create_filterbank_file(outfn, header, spectra=None, nbits=8, verbose=False):
    """Write filterbank header and spectra to file.

        Input:
            outfn: The outfile filterbank file's name.
            header: A dictionary of header paramters and values.
            spectra: Spectra to write to file. (Default: don't write
                any spectra - i.e. write out header only)
            nbits: The number of bits per sample of the filterbank file.
                (Default: 8 - i.e. each sample is an 8-bit integer)
            verbose: If True, be verbose (Default: be quiet)

        Output:
            fbfile: The resulting FilterbankFile object opened
                in read-write mode.
    """
    dtype = get_dtype(nbits) # Get dtype. This will check to ensure
                             # 'nbits' is valid.
    header['nbits'] = nbits
    outfile = open(outfn, 'wb')
    outfile.write(sigproc.addto_hdr("HEADER_START", None))
    for paramname in header.keys():
        if verbose:
            print "Writing header param (%s)" % paramname
        value = header[paramname]
        outfile.write(sigproc.addto_hdr(paramname, value))
    outfile.write(sigproc.addto_hdr("HEADER_END", None))
    if spectra:
        spectra.flatten().astype(dtype).tofile(outfile)
    outfile.close()
    return FilterbankFile(outfn, read_only=False)

def is_float(nbits):
    """For a given number of bits per sample return
        true if it corresponds to floating-point samples
        in filterbank files.

        Input:
            nbits: Number of bits per sample, as recorded in the filterbank
                file's header.

        Output:
            isfloat: True, if 'nbits' indicates the data in the file
                are encoded as floats.
    """
    check_nbits(nbits)
    if nbits == 32:
        return True
    else:
        return False


def check_nbits(nbits):
    """Given a number of bits per sample check to make
        sure 'filterbank.py' can cope with it.

        An exception is raise if 'filterbank.py' cannot cope.

        Input:
            nbits: Number of bits per sample, as recorded in the filterbank
                file's header.

        Output:
            None
    """
    if nbits not in [32, 16, 8]:
        raise ValueError("'filterbank.py' only supports " \
                                    "files with 8- or 16-bit " \
                                    "integers, or 32-bit floats " \
                                    "(nbits provided: %g)!" % nbits)


def get_dtype(nbits):
    """For a given number of bits per sample return
        a numpy-recognized dtype.

        Input:
            nbits: Number of bits per sample, as recorded in the filterbank
                file's header.

        Output:
            dtype: A numpy-recognized dtype string.
    """
    check_nbits(nbits)
    if is_float(nbits):
        dtype = 'float%d' % nbits
    else:
        dtype = 'uint%d' % nbits
    return dtype


def read_header(filename, verbose=False):
    """Read the header of a filterbank file, and return
        a dictionary of header paramters and the header's
        size in bytes.

        Inputs:
            filename: Name of the filterbank file.
            verbose: If True, be verbose. (Default: be quiet)

        Outputs:
            header: A dictionary of header paramters.
            header_size: The size of the header in bytes.
    """
    header = {}
    filfile = open(filename, 'rb')
    filfile.seek(0)
    paramname = ""
    while (paramname != 'HEADER_END'):
        if verbose:
            print "File location: %d" % filfile.tell()
        paramname, val = sigproc.read_hdr_val(filfile, stdout=verbose)
        if verbose:
            print "Read param %s (value: %s)" % (paramname, val)
        if paramname not in ["HEADER_START", "HEADER_END"]:
            header[paramname] = val
    header_size = filfile.tell()
    filfile.close()
    return header, header_size


class FilterbankFile(object):
    def __init__(self, filfn, read_only=True):
        if not os.path.isfile(filfn):
            raise ValueError("ERROR: File does not exist!\n\t(%s)" % filfn)
        self.filename = filfn
        self.read_only = read_only
        self.header, self.header_size = read_header(self.filename)
        self.frequencies = self.fch1 + self.foff*np.arange(self.nchans)
        self.is_hifreq_first = (self.foff < 0)
        self.bytes_per_spectrum = self.nchans*self.nbits / 8
        
        # Get info about dtype
        self.dtype = get_dtype(self.nbits)
        if is_float(self.nbits):
            tinfo = np.finfo(self.dtype)
        else:
            tinfo = np.iinfo(self.dtype)
        self.dtype_min = tinfo.min
        self.dtype_max = tinfo.max

        self.needs_sync = True
        self.sync_spectra()

        if not self.read_only:
            self.file_append_mode = open(self.filename, 'ab')

    def __del__(self):
        self.close()
        
    def close(self):
        if not self.read_only:
            self.file_append_mode.close()

    def sync_spectra(self):
        self.file_size = os.stat(self.filename)[6]
        self.data_size = self.file_size - self.header_size
        self.nspec = self.data_size / self.bytes_per_spectrum

        if self.data_size % self.bytes_per_spectrum:
            warnings.warn("Not an integer number of spectra in file.")

        mode = (self.read_only and 'r') or 'r+'
        self.spectra = np.memmap(self.filename, dtype=self.dtype, \
                                    mode=mode, offset=self.header_size, \
                                    shape=(self.nspec, self.nchans))
        self.needs_sync = False

    def get_timeslice(self, start, stop):
        startbins = int(np.round(start/self.tsamp))
        stopbins = int(np.round(stop/self.tsamp))
        return self.get_spectra(startbins, stopbins)

    def get_spectra(self, start, stop):
        if self.needs_sync:
            warnings.warn("Spectra haven't been sync'ed recently.")
        return self.spectra[start:stop]

    def append_spectra(self, spectra):
        """Append spectra to the file if is not read-only.
            
            Input:
                spectra: The spectra to append. The new spectra
                    must have the correct number of channels (ie
                    dimension of axis=1.

            Outputs:
                None
        """
        if self.read_only:
            raise ValueError("FilterbankFile object for '%s' is read-only." % \
                        self.filename)
        nspec, nchans = spectra.shape
        if nchans != self.nchans:
            raise ValueError("Cannot append spectra. Incorrect shape. " \
                        "Number of channels in file: %d; Number of " \
                        "channels in spectra to append: %d" % \
                        (self.nchans, nchans))
        data = spectra.flatten()
        np.clip(data, self.dtype_min, self.dtype_max, out=data)
        self.file_append_mode.write(data.astype(self.dtype))
        self.file_append_mode.flush()
        os.fsync(self.file_append_mode)
        self.needs_sync = True

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


def main():
    fil = FilterbankFile(sys.argv[1])
    fil.print_header()


if __name__ == '__main__':
    main()
