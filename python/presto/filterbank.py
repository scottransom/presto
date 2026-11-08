"""
A module for reading filterbank files.

Patrick Lazarus, June 26, 2012
(Minor modification from file originally from June 6th, 2009)
"""

from __future__ import annotations

import sys
import os
import os.path

import numpy as np

from presto import sigproc
from presto import spectra


DEBUG = False


def create_filterbank_file(
    outfn: str,
    header: dict,
    spectra: np.ndarray | None = None,
    nbits: int = 8,
    verbose: bool = False,
    mode: str = "append",
) -> FilterbankFile:
    """
    Write a filterbank header and spectra to a file.

    Parameters
    ----------
    outfn : str
        The output filterbank file's name.
    header : dict
        A dictionary of header parameters and values.
    spectra : numpy.ndarray, optional
        Spectra to write to the file. If None (default), only the header is
        written.
    nbits : int, optional
        The number of bits per sample of the filterbank file (default 8, so
        each sample is an 8-bit integer). This value always overrides the
        value in the `header` dictionary.
    verbose : bool, optional
        If True, be verbose (default False).
    mode : str, optional
        The mode for writing, either "append" or "write" (default "append").

    Returns
    -------
    FilterbankFile
        The resulting FilterbankFile object, opened in read-write mode.
    """
    # Get dtype. This also checks that 'nbits' is valid.
    dtype = get_dtype(nbits, signed=bool(header.get("signed", 0)))
    header["nbits"] = nbits
    outfile = open(outfn, "wb")
    outfile.write(sigproc.addto_hdr("HEADER_START", None))
    for paramname in list(header.keys()):
        if paramname not in sigproc.header_params:
            # Only add recognized parameters
            continue
        if verbose:
            print("Writing header param (%s)" % paramname)
        value = header[paramname]
        outfile.write(sigproc.addto_hdr(paramname, value))
    outfile.write(sigproc.addto_hdr("HEADER_END", None))
    if spectra is not None:
        if nbits < 8:
            spectra = np.clip(spectra, 0, 2**nbits - 1)
        _encode_spectra(spectra, nbits, header["nchans"], dtype).tofile(outfile)
    outfile.close()
    return FilterbankFile(outfn, mode=mode)


def is_float(nbits: int) -> bool:
    """
    Return whether a number of bits per sample corresponds to floats.

    Parameters
    ----------
    nbits : int
        The number of bits per sample, as recorded in the filterbank file's
        header.

    Returns
    -------
    bool
        True if `nbits` indicates the data are encoded as floats.
    """
    check_nbits(nbits)
    if nbits == 32:
        return True
    else:
        return False


def check_nbits(nbits: int) -> None:
    """
    Check that a number of bits per sample is supported.

    Parameters
    ----------
    nbits : int
        The number of bits per sample, as recorded in the filterbank file's
        header.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``filterbank.py`` cannot cope with the given `nbits` (1-, 2-, 4-,
        8- and 16-bit integers and 32-bit floats are supported).
    """
    if nbits not in [32, 16, 8, 4, 2, 1]:
        raise ValueError(
            "'filterbank.py' only supports "
            "files with 1-, 2-, 4-, 8- or 16-bit "
            "integers, or 32-bit floats "
            "(nbits provided: %g)!" % nbits
        )


def get_dtype(nbits: int, signed: bool = False) -> str:
    """
    Return a numpy dtype string for a number of bits per sample.

    Parameters
    ----------
    nbits : int
        The number of bits per sample, as recorded in the filterbank file's
        header.
    signed : bool, optional
        Whether the header's ``signed`` key is set (default False). As in
        PRESTO's C reader and in SIGPROC, it applies only to 8-bit data.

    Returns
    -------
    str
        A numpy-recognized dtype string (e.g. ``"uint8"``, ``"int8"`` or
        ``"float32"``). 1-, 2- and 4-bit samples are unpacked into ``"uint8"``.
    """
    check_nbits(nbits)
    if is_float(nbits):
        dtype = "float%d" % nbits
    elif signed and nbits == 8:
        dtype = "int8"
    elif nbits < 8:
        dtype = "uint8"
    else:
        dtype = "uint%d" % nbits
    return dtype


def pack_bits(samples: np.ndarray, nbits: int) -> np.ndarray:
    """
    Pack 1-, 2- or 4-bit samples into bytes, in SIGPROC order.

    Each byte holds ``8 // nbits`` consecutive channels, with the first of them
    in the least significant bits. This matches PRESTO's C reader and SIGPROC.

    Parameters
    ----------
    samples : numpy.ndarray
        Integer samples of shape ``(nspec, nchans)``, each in
        ``0 .. 2**nbits - 1``. ``nchans * nbits`` must be a multiple of 8.
    nbits : int
        Bits per sample: 1, 2 or 4.

    Returns
    -------
    numpy.ndarray
        The packed bytes, shape ``(nspec, nchans * nbits // 8)``, dtype uint8.
    """
    per_byte = 8 // nbits
    nspec, nchans = samples.shape
    if nchans % per_byte:
        raise ValueError(
            "%d channels at %d bits do not fill whole bytes" % (nchans, nbits)
        )
    vals = samples.astype(np.uint8).reshape(nspec, nchans // per_byte, per_byte)
    packed = np.zeros((nspec, nchans // per_byte), dtype=np.uint8)
    for ii in range(per_byte):
        packed |= vals[:, :, ii] << np.uint8(ii * nbits)
    return packed


def unpack_bits(packed: np.ndarray, nbits: int, nchans: int) -> np.ndarray:
    """
    Unpack bytes of 1-, 2- or 4-bit samples (the inverse of `pack_bits`).

    Parameters
    ----------
    packed : numpy.ndarray
        Packed uint8 data holding a whole number of spectra.
    nbits : int
        Bits per sample: 1, 2 or 4.
    nchans : int
        Number of channels per spectrum.

    Returns
    -------
    numpy.ndarray
        The samples, shape ``(nspec, nchans)``, dtype uint8.
    """
    per_byte = 8 // nbits
    mask = np.uint8(2**nbits - 1)
    packed = np.asarray(packed, dtype=np.uint8).reshape(-1, nchans // per_byte)
    vals = np.empty(packed.shape + (per_byte,), dtype=np.uint8)
    for ii in range(per_byte):
        vals[:, :, ii] = (packed >> np.uint8(ii * nbits)) & mask
    return vals.reshape(-1, nchans)


def _encode_spectra(data: np.ndarray, nbits: int, nchans: int, dtype: str) -> np.ndarray:
    """Return spectra as the array of bytes to write: cast, and packed if nbits < 8."""
    data = np.asarray(data).reshape(-1, nchans)
    if nbits < 8:
        return pack_bits(data, nbits)
    return data.astype(dtype)


def read_header(filename: str, verbose: bool = False) -> tuple[dict, int]:
    """
    Read the header of a filterbank file.

    Parameters
    ----------
    filename : str
        The name of the filterbank file.
    verbose : bool, optional
        If True, be verbose (default False).

    Returns
    -------
    header : dict
        A dictionary of header parameters.
    header_size : int
        The size of the header in bytes.
    """
    header = {}
    filfile = open(filename, "rb")
    filfile.seek(0)
    paramname = ""
    while paramname != "HEADER_END":
        if verbose:
            print("File location: %d" % filfile.tell())
        paramname, val = sigproc.read_hdr_val(filfile, stdout=verbose)
        if verbose:
            print("Read param %s (value: %s)" % (paramname, val))
        if paramname not in ["HEADER_START", "HEADER_END"]:
            header[paramname] = val
    header_size = filfile.tell()
    filfile.close()
    return header, header_size


class FilterbankFile(object):
    """
    A SIGPROC filterbank file.

    Header parameters are accessible as attributes via ``__getattr__``
    (e.g. ``fch1``, ``foff``, ``nchans``, ``nbits``, ``tsamp``), in addition
    to the derived attributes documented below.

    Parameters
    ----------
    filfn : str
        The path to the filterbank file.
    mode : str, optional
        How to open the file: "read"/"readonly", "write"/"readwrite", or
        "append" (default "readonly").

    Attributes
    ----------
    filename : str
        The path to the filterbank file.
    header : dict
        The parsed header parameters.
    header_size : int
        The size of the header in bytes.
    frequencies : numpy.ndarray
        The channel center frequencies (MHz).
    is_hifreq_first : bool
        True if the highest frequency channel is first (``foff < 0``).
    bytes_per_spectrum : int
        The number of bytes in a single spectrum.
    nspec : int
        The number of spectra in the file.
    isfold : bool
        True if the file is a folded-filterbank file.
    dt : float
        The time per spectrum in sec (the sample time, or period/nbins for a
        folded file).
    signed : bool
        Whether the header's ``signed`` key is set (meaningful for 8-bit data).
    dtype : str
        The numpy dtype string of the samples (``int8`` for signed 8-bit data).
    dtype_min, dtype_max : int or float
        The minimum and maximum representable sample values.

    Raises
    ------
    ValueError
        If the file does not exist, or `mode` is unrecognized.
    """

    def __init__(self, filfn: str, mode: str = "readonly"):
        self.filename = filfn
        self.filfile = None
        if not os.path.isfile(filfn):
            raise ValueError("ERROR: File does not exist!\n\t(%s)" % filfn)
        self.header, self.header_size = read_header(self.filename)
        self.frequencies = self.fch1 + self.foff * np.arange(self.nchans)
        self.is_hifreq_first = self.foff < 0
        self.bytes_per_spectrum = self.nchans * self.nbits // 8
        data_size = os.path.getsize(self.filename) - self.header_size
        self.nspec = data_size // self.bytes_per_spectrum

        # Check if this file is a folded-filterbank file
        if (
            "npuls" in self.header
            and "period" in self.header
            and "nbins" in self.header
            and "tsamp" not in self.header
        ):
            # Foleded file
            self.isfold = True
            self.dt = self.period / self.nbins
        else:
            self.isfold = False
            self.dt = self.tsamp

        # Get info about dtype
        self.signed = bool(self.header.get("signed", 0))
        self.dtype = get_dtype(self.nbits, signed=self.signed)
        if is_float(self.nbits):
            tinfo = np.finfo(self.dtype)
        else:
            tinfo = np.iinfo(self.dtype)
        self.dtype_min = tinfo.min
        self.dtype_max = tinfo.max
        if self.nbits < 8:
            # Samples are unpacked into uint8 but only span nbits
            self.dtype_max = 2**self.nbits - 1
            if (self.nchans * self.nbits) % 8:
                raise ValueError(
                    "%d channels at %d bits do not fill whole bytes"
                    % (self.nchans, self.nbits)
                )

        if mode.lower() in ("read", "readonly"):
            self.filfile = open(self.filename, "rb")
        elif mode.lower() in ("write", "readwrite"):
            self.filfile = open(self.filename, "r+b")
        elif mode.lower() == "append":
            self.filfile = open(self.filename, "a+b")
        else:
            raise ValueError("Unrecognized mode (%s)!" % mode)

    @property
    def freqs(self) -> np.ndarray:
        # Alias for frequencies
        return self.frequencies

    @property
    def nchan(self) -> int:
        # more aliases..
        return self.nchans

    def close(self) -> None:
        """Close the underlying file if it is open."""
        if self.filfile is not None:
            self.filfile.close()

    def get_timeslice(self, start: float, stop: float) -> spectra.Spectra:
        """
        Return the spectra between two times.

        Parameters
        ----------
        start : float
            The start time in sec.
        stop : float
            The stop time in sec.

        Returns
        -------
        spectra.Spectra
            The spectra covering the requested time slice.
        """
        startspec = int(np.round(start / self.tsamp))
        stopspec = int(np.round(stop / self.tsamp))
        return self.get_spectra(startspec, stopspec - startspec)

    def get_spectra(self, start: int, nspec: int) -> spectra.Spectra:
        """
        Return a block of spectra from the file.

        Parameters
        ----------
        start : int
            The index of the first spectrum to read.
        nspec : int
            The number of spectra to read (clipped to the end of the file).

        Returns
        -------
        spectra.Spectra
            The requested spectra.
        """
        stop = min(start + nspec, self.nspec)
        pos = self.header_size + start * self.bytes_per_spectrum
        # Compute number of elements to read
        nspec = int(stop) - int(start)
        nspec = max(0, nspec)
        self.filfile.seek(pos, os.SEEK_SET)
        if self.nbits < 8:
            packed = np.fromfile(
                self.filfile, dtype=np.uint8, count=nspec * self.bytes_per_spectrum
            )
            spectra_dat = unpack_bits(packed, self.nbits, self.nchans)
        else:
            spectra_dat = np.fromfile(
                self.filfile, dtype=self.dtype, count=nspec * self.nchans
            )
        spectra_dat = np.reshape(spectra_dat, (nspec, self.nchans))
        spec = spectra.Spectra(
            self.freqs,
            self.tsamp,
            spectra_dat.T,
            starttime=start * self.tsamp,
            dm=0.0,
        )
        return spec

    def append_spectra(self, spectra: np.ndarray) -> None:
        """
        Append spectra to the file if it is not read-only.

        Parameters
        ----------
        spectra : numpy.ndarray
            The spectra to append. The new spectra must have the correct
            number of channels (i.e. the dimension of axis=1).

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the file is read-only, or the spectra have the wrong number
            of channels.
        """
        if self.filfile.mode.lower() in ("r", "rb"):
            raise ValueError(
                "FilterbankFile object for '%s' is read-only." % self.filename
            )
        nspec, nchans = spectra.shape
        if nchans != self.nchans:
            raise ValueError(
                "Cannot append spectra. Incorrect shape. "
                "Number of channels in file: %d; Number of "
                "channels in spectra to append: %d" % (self.nchans, nchans)
            )
        data = spectra.flatten()
        np.clip(data, self.dtype_min, self.dtype_max, out=data)
        # Move to end of file
        self.filfile.seek(0, os.SEEK_END)
        self.filfile.write(_encode_spectra(data, self.nbits, self.nchans, self.dtype))
        self.nspec += nspec
        # self.filfile.flush()
        # os.fsync(self.filfile)

    def write_spectra(self, spectra: np.ndarray, ispec: int) -> None:
        """
        Write spectra to the file if it is writable.

        Parameters
        ----------
        spectra : numpy.ndarray
            The spectra to write. The new spectra must have the correct
            number of channels (i.e. the dimension of axis=1).
        ispec : int
            The index of the spectrum at which to start writing.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the file is not writable, the spectra have the wrong number
            of channels, or the write would go past the end of the file.
        """
        if "r+" not in self.filfile.mode.lower():
            raise ValueError(
                "FilterbankFile object for '%s' is not writable." % self.filename
            )
        nspec, nchans = spectra.shape
        if nchans != self.nchans:
            raise ValueError(
                "Cannot write spectra. Incorrect shape. "
                "Number of channels in file: %d; Number of "
                "channels in spectra to write: %d" % (self.nchans, nchans)
            )
        if ispec > self.nspec:
            raise ValueError(
                "Cannot write past end of file! "
                "Present number of spectra: %d; "
                "Requested index of write: %d" % (self.nspec, ispec)
            )
        data = spectra.flatten()
        np.clip(data, self.dtype_min, self.dtype_max, out=data)
        # Move to requested position
        pos = self.header_size + ispec * self.bytes_per_spectrum
        self.filfile.seek(pos, os.SEEK_SET)
        self.filfile.write(_encode_spectra(data, self.nbits, self.nchans, self.dtype))
        if nspec + ispec > self.nspec:
            self.nspec = nspec + ispec

    def __getattr__(self, name: str):
        if name in self.header:
            if DEBUG:
                print("Fetching header param (%s)" % name)
            val = self.header[name]
        else:
            raise ValueError("No FilterbankFile attribute called '%s'" % name)
        return val

    def print_header(self) -> None:
        """Print the header parameters and values."""
        for param in sorted(self.header.keys()):
            if param in ("HEADER_START", "HEADER_END"):
                continue
            print("%s: %s" % (param, self.header[param]))


def main() -> None:
    fil = FilterbankFile(sys.argv[1])
    fil.print_header()


if __name__ == "__main__":
    main()
