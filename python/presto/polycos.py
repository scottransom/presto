from __future__ import annotations

import os
import re
import sys
import subprocess
from typing import TextIO

import numpy as np

from presto import parfile

# Constants
NUMCOEFFS_DEFAULT = 12
SPAN_DEFAULT = 60  # span of each polyco in minutes

# Telescope name to TEMPO observatory code conversion (kept for backwards
# compatibility -- polyco generation now uses tempo2 observatory names).
# The codes are TEMPO's, from its obsys.dat.  Never hand them to tempo2: many
# of its one-character aliases mean other sites (e.g. "y" is Jodrell Bank, not
# CHIME, "6" is DSS-43, not the VLA, and "k" is Tabley, not FAST), so
# create_polycos() goes from code to telescope name to tempo2 site name.
telescope_to_id = {
    "GBT": "1",
    "Arecibo": "3",
    "VLA": "6",
    "Parkes": "7",
    "Jodrell": "8",
    "GB43m": "a",
    "GB 140FT": "a",
    "Nancay": "f",
    "Effelsberg": "g",
    "WSRT": "i",
    "FAST": "k",
    "GMRT": "r",
    "CHIME": "y",
    "MeerKAT": "m",
    "Hobart": "4",
    "Princeton": "5",
    "GB 300FT": "9",
    "GB 85-3": "b",
    "MOST": "e",
    "Jodrell MkII": "h",
    "SHAO": "s",
    "LOFAR": "t",
    "MWA": "u",
    "Pico Veleta": "v",
    "SRT": "z",
    "Geocenter": "0",
    "Barycenter": "@",
}

# TEMPO observatory code to Telescope name conversion
id_to_telescope = {
    "1": "GBT",
    "3": "Arecibo",
    "6": "VLA",
    "7": "Parkes",
    "8": "Jodrell",
    "a": "GB 140FT",
    "f": "Nancay",
    "g": "Effelsberg",
    "i": "WSRT",
    "k": "FAST",
    "r": "GMRT",
    "y": "CHIME",
    "m": "MeerKAT",
    "c": "VLA",  # TEMPO's "VLA SITE"
    "4": "Hobart",
    "5": "Princeton",
    "9": "GB 300FT",
    "b": "GB 85-3",
    "e": "MOST",
    "h": "Jodrell MkII",
    "s": "SHAO",
    "t": "LOFAR",
    "u": "MWA",
    "v": "Pico Veleta",
    "z": "SRT",
    "0": "Geocenter",
    "@": "Barycenter",
}

# Telescope name to tempo2 observatory name conversion
telescope_to_tempo2_name = {
    "GBT": "gbt",
    "Arecibo": "ao",
    "VLA": "vla",
    "Parkes": "pks",
    "Jodrell": "jb",
    "GB43m": "gb140",
    "GB 140FT": "gb140",
    "Nancay": "ncy",
    "Effelsberg": "eff",
    "ATA": "ata",
    "LOFAR": "lofar",
    "WSRT": "wsrt",
    "FAST": "fast",
    "GMRT": "gmrt",
    "CHIME": "chime",
    "MWA": "mwa",
    "LWA": "lwa1",
    "SRT": "srt",
    # tempo2 writes the site name with "%5s", so one of 5 or more characters
    # (e.g. "meerkat") runs into F0 on the polyco's second line.  polyco()
    # copes, but other readers may not.  "m" is tempo2's alias for MeerKAT.
    "MeerKAT": "m",
    "KAT-7": "k7",
    "Hobart": "hob",
    "Princeton": "princeton",
    "GB 300FT": "gb300",
    "GB 85-3": "gb853",
    "MOST": "mo",
    "Jodrell MkII": "jbmk2",
    "SHAO": "shao",
    "Pico Veleta": "pv",
    "Geocenter": "coe",
    "Barycenter": "@",
}

# Telescope name to track length (max hour angle) conversion
telescope_to_maxha = {
    "GBT": 12,
    "Arecibo": 3,
    "FAST": 5,
    "VLA": 6,
    "Parkes": 12,
    "Jodrell": 12,
    "GB43m": 12,
    "GB 140FT": 12,
    "Nancay": 4,
    "Effelsberg": 12,
    "ATA": 12,
    "LOFAR": 12,
    "WSRT": 12,
    "GMRT": 12,
    "CHIME": 1,
    "MWA": 12,
    "LWA": 12,
    "SRT": 12,
    "MeerKAT": 12,
    "KAT-7": 12,
    "Hobart": 12,
    "Princeton": 12,
    "GB 300FT": 12,
    "GB 85-3": 12,
    "MOST": 12,
    "Jodrell MkII": 12,
    "SHAO": 12,
    "Pico Veleta": 12,
    "Geocenter": 12,
    "Barycenter": 12,
}


class polyco(object):
    """
    A single TEMPO1-format polyco entry.

    One entry is parsed from an open polyco file, giving the reference epoch
    and the phase polynomial coefficients used to predict pulsar phase and
    spin frequency near that epoch.

    Parameters
    ----------
    fileptr : file object
        An open polyco file, positioned at the start of an entry. If the
        end of the file is reached, `psr` is set to None.

    Attributes
    ----------
    psr : str or None
        The pulsar name, or None if the end of the file was reached.
    date : str
        The date of the reference epoch.
    UTC : str
        The UTC of the reference epoch.
    TMIDi, TMIDf, TMID : float
        The integer, fractional, and full MJD of the reference epoch.
    DM : float
        The dispersion measure (pc cm^-3).
    doppler : float
        The observatory-frame Doppler shift (v/c).
    log10rms : float
        The log10 of the fit RMS in periods.
    RPHASE : float
        The reference phase.
    F0 : float
        The reference rotation frequency (Hz).
    obs : str
        The observatory code.
    dataspan : int
        The span of the entry in minutes.
    numcoeff : int
        The number of polynomial coefficients.
    obsfreq : float
        The observing frequency (MHz).
    binphase : float
        The binary/orbital phase. Only set when the polyco's second line has
        exactly seven whitespace-separated fields; tempo2's ``-tempo1`` output
        uses an eight-field layout, so in practice this is usually absent. Not
        currently used elsewhere in PRESTO.
    coeffs : numpy.ndarray
        The phase polynomial coefficients.
    phasepoly : numpy.polynomial.polynomial.Polynomial
        The phase polynomial built from `coeffs`.
    """

    def __init__(self, fileptr: TextIO):
        line = fileptr.readline()
        if line == "":
            self.psr = None
        else:
            sl = line.split()
            self.psr = sl[0]
            self.date = sl[1]
            self.UTC = sl[2]
            self.TMIDi = float(sl[3].split(".")[0])
            self.TMIDf = float("0." + sl[3].split(".")[1])
            self.TMID = self.TMIDi + self.TMIDf
            self.DM = float(sl[4])
            if len(sl) == 7:
                self.doppler = float(sl[5]) * 1e-4
                self.log10rms = float(sl[6])
            else:
                self.log10rms = "-" + sl[-1].split("-")[-1]
                self.doppler = float(sl[-1][: sl[-1].find(self.log10rms)]) * 1e-4
                self.log10rms = float(self.log10rms)
            sl = fileptr.readline().split()
            # tempo2 prints the site name with "%5s", so one of 5 or more
            # characters runs straight into F0: "11.188185028329meerkat".
            fused = re.fullmatch(r"([-+]?\d*\.\d*)([^\d.]\S*)", sl[1])
            if fused:
                sl[1:2] = fused.groups()
            self.RPHASE = float(sl[0])
            self.F0 = float(sl[1])
            self.obs = sl[2]
            self.dataspan = int(sl[3])
            self.numcoeff = int(sl[4])
            self.obsfreq = float(sl[5])
            if len(sl) == 7:
                self.binphase = float(sl[6])
            self.coeffs = np.zeros(self.numcoeff, "d")
            for linenum in range(self.numcoeff // 3):
                sl = fileptr.readline().split()
                self.coeffs[linenum * 3 + 0] = float(sl[0].replace("D", "E"))
                self.coeffs[linenum * 3 + 1] = float(sl[1].replace("D", "E"))
                self.coeffs[linenum * 3 + 2] = float(sl[2].replace("D", "E"))
            if self.numcoeff % 3 != 0:  # get remaining terms if needed
                sl = fileptr.readline().split()
                nlines = self.numcoeff // 3
                for coeffnum in range(len(sl)):
                    self.coeffs[nlines * 3 + coeffnum] = float(
                        sl[coeffnum].replace("D", "E")
                    )
            self.phasepoly = np.polynomial.polynomial.Polynomial(self.coeffs)

    def phase(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar phase (0-1) at a given MJD.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted pulsar phase, in the range [0, 1).
        """
        return self.rotation(mjdi, mjdf) % 1

    def rotation(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar rotation at a given MJD.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted (fractional) pulsar rotation number.
        """
        DT = ((mjdi - self.TMIDi) + (mjdf - self.TMIDf)) * 1440.0
        phase = self.phasepoly(DT)
        # phase = self.coeffs[self.numcoeff-1]
        # for ii in range(self.numcoeff-1, 0, -1):
        #    phase = DT*phase + self.coeffs[ii-1]
        phase += self.RPHASE + DT * 60.0 * self.F0
        return phase

    def freq(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar spin frequency at a given MJD.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted pulsar spin frequency (Hz).
        """
        DT = ((mjdi - self.TMIDi) + (mjdf - self.TMIDf)) * 1440.0
        psrfreq = 0.0
        for ii in range(self.numcoeff - 1, 0, -1):
            psrfreq = DT * psrfreq + ii * self.coeffs[ii]
        return self.F0 + psrfreq / 60.0


class polycos(object):
    """
    A collection of polyco entries for a single pulsar.

    Reads all the entries for `psrname` from a polyco file and provides
    methods to predict phase, rotation, spin frequency, and v/c at a given
    MJD by selecting the appropriate entry.

    Parameters
    ----------
    psrname : str
        The pulsar name to select entries for.
    filenm : str, optional
        The polyco file to read (default "polyco.dat").

    Attributes
    ----------
    psr : str
        The pulsar name.
    file : str
        The polyco filename that was read.
    polycos : list of polyco
        The individual polyco entries for this pulsar.
    TMIDs : numpy.ndarray
        The reference-epoch MJDs of the entries.
    dataspan : int
        The span of each entry in minutes.
    validrange : float
        Half the data span, in days (the +/- validity window of an entry).
    """

    def __init__(self, psrname: str, filenm: str = "polyco.dat"):
        self.psr = psrname
        self.file = filenm
        self.polycos = []
        self.TMIDs = []
        infile = open(filenm, "r")
        tmppoly = polyco(infile)
        while tmppoly.psr:
            if len(self.polycos):
                if tmppoly.dataspan != self.dataspan:
                    sys.stderr.write("Data span is changing!\n")
            else:
                self.dataspan = tmppoly.dataspan
            if tmppoly.psr == psrname:
                self.polycos.append(tmppoly)
                self.TMIDs.append(tmppoly.TMID)
            tmppoly = polyco(infile)
        sys.stderr.write("Read %d polycos for PSR %s\n" % (len(self.polycos), psrname))
        self.TMIDs = np.asarray(self.TMIDs)
        infile.close()
        self.validrange = 0.5 * self.dataspan / 1440.0

    def select_polyco(self, mjdi: int, mjdf: float) -> int:
        """
        Return the index of the polyco entry valid for the specified time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        int
            The index into `self.polycos` of the closest entry. A warning is
            written to stderr if the closest entry is outside its validity
            range.
        """
        goodpoly = np.argmin(np.fabs(self.TMIDs - (mjdi + mjdf)))
        if np.fabs(self.TMIDs[goodpoly] - (mjdi + mjdf)) > self.validrange:
            sys.stderr.write("Cannot find a valid polyco at %f!\n" % (mjdi + mjdf))
        return goodpoly

    def get_phase(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar phase for the specified time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted pulsar phase, in the range [0, 1).
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].phase(mjdi, mjdf)

    def get_rotation(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar rotation number for the specified time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted (fractional) pulsar rotation number.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].rotation(mjdi, mjdf)

    def get_freq(self, mjdi: int, mjdf: float) -> float:
        """
        Return the predicted pulsar spin frequency for the specified time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The predicted pulsar spin frequency (Hz).
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].freq(mjdi, mjdf)

    def get_phs_and_freq(self, mjdi: int, mjdf: float) -> tuple[float, float]:
        """
        Return the predicted pulsar phase and spin frequency for a time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        phase : float
            The predicted pulsar phase, in the range [0, 1).
        freq : float
            The predicted pulsar spin frequency (Hz).
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return (
            self.polycos[goodpoly].phase(mjdi, mjdf),
            self.polycos[goodpoly].freq(mjdi, mjdf),
        )

    def get_voverc(self, mjdi: int, mjdf: float) -> float:
        """
        Return the (approximate) topocentric v/c for the specified time.

        Parameters
        ----------
        mjdi : int
            The integer part of the MJD.
        mjdf : float
            The fractional part of the MJD.

        Returns
        -------
        float
            The approximate topocentric v/c (from the selected entry's
            Doppler term).
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].doppler


def create_polycos(
    parfn: str | parfile.psr_par,
    telescope_id: str,
    center_freq: float,
    start_mjd: int,
    end_mjd: int,
    max_hour_angle: int | None = None,
    span: int = SPAN_DEFAULT,
    numcoeffs: int = NUMCOEFFS_DEFAULT,
    keep_file: bool = False,
) -> polycos:
    """Create a polycos object from a parfile using tempo2.

    This runs ``tempo2 -tempo1 -polyco``, which generates TEMPO1-format
    polycos, so the external ``tempo2`` executable (with its ``$TEMPO2``
    runtime directory, both available from conda-forge) is required.
    TEMPO itself is no longer used.

    Parameters
    ----------
    parfn : str or parfile.psr_par
        The parfile's filename, or a parfile object.
    telescope_id : str
        The TEMPO 1-character telescope identifier (e.g. '1' for the
        GBT), kept for backwards compatibility.
    center_freq : float
        The observation's center frequency in MHz.
    start_mjd : int
        MJD on which the polycos should start.
    end_mjd : int
        MJD until which the polycos should extend.
    max_hour_angle : int, optional
        The maximum hour angle.  Default: the value appropriate for the
        given telescope (see `telescope_to_maxha`).
    span : int, optional
        Span of each set of polycos in minutes.  Default: 60.
    numcoeffs : int, optional
        Number of polynomial coefficients.  Default: 12.
    keep_file : bool, optional
        If True, keep the generated file as 'polyco.dat'.  Default:
        delete it.

    Returns
    -------
    polycos
        A polycos object built from the generated file.
    """
    if isinstance(parfn, (str, bytes)):
        # assume parfn is a filename
        par = parfile.psr_par(parfn)
    else:
        # assume par is already a parfile.psr_par object
        par = parfn

    telescope_name = id_to_telescope[telescope_id.strip()]
    if max_hour_angle is None:
        max_hour_angle = telescope_to_maxha[telescope_name]
    tempo2_name = telescope_to_tempo2_name[telescope_name]

    if hasattr(par, "PSR"):
        psrname = par.PSR
    else:
        psrname = par.PSRJ
    command = 'tempo2 -tempo1 -f %s -polyco "%d %d %d %d %d %s %0.5f"' % (
        par.FILE,
        start_mjd,
        end_mjd,
        span,
        numcoeffs,
        max_hour_angle,
        tempo2_name,
        center_freq,
    )
    tempo2 = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Note: tempo2 routinely prints warnings to stderr, so only check
    # the return code and that the output file actually appeared
    if tempo2.returncode != 0 or not os.path.exists("polyco_new.dat"):
        raise Tempo2Error(
            "The following was encountered when running "
            "tempo2 to generate polycos from the input "
            "parfile (%s):\n\n%s\n%s\n" % (parfn, tempo2.stdout, tempo2.stderr)
        )
    os.rename("polyco_new.dat", "polyco.dat")
    new_polycos = polycos(psrname, filenm="polyco.dat")
    # Remove other files created by tempo2
    for fn in ("newpolyco.dat", "polyco.tim"):
        if os.path.exists(fn):
            os.remove(fn)
    if not keep_file:
        os.remove("polyco.dat")
    return new_polycos


class Tempo2Error(Exception):
    pass


# Backwards-compatible alias
TempoError = Tempo2Error
