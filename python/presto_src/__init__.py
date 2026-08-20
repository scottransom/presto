from __future__ import annotations

import os.path

import numpy as np

from .prestoswig import *  # noqa: F401,F403  (SWIG-generated C symbols)
from presto import Pgplot
from presto import psr_utils
from presto.observatories import obscode


def val_with_err(
    value: float, error: float, length: int = 0, digits: int = 2, latex: int = 0
) -> str:
    """
    Format a value and its error as in scientific journals.

    Returns a string with `value` rounded to the appropriate decimal place
    and `error` in parentheses (with `digits` decimal places).

    Parameters
    ----------
    value : float
        The value to format.
    error : float
        The uncertainty on `value`.
    length : int, optional
        The length of the returned string (auto if 0). Should be ~20 to show
        full double precision if the base-10 exponent of the error needs to
        be shown. If 0, a left-justified minimum-length string is returned;
        if > 0, the string is right-justified; if < 0, it is left-justified.
    digits : int, optional
        The number of decimal places on the error (default 2).
    latex : int, optional
        If nonzero, convert the string to LaTeX markup (default 0).

    Returns
    -------
    str
        The formatted ``value(error)`` string.
    """
    slen = 40
    outstr = " " * slen
    if abs(length) > slen:
        slen = abs(length)
    if digits == 2:
        slen = nice_output_2(outstr, value, error, length)
    else:
        slen = nice_output_1(outstr, value, error, length)
    outstr = outstr[:slen].strip()  # remove null termination and any space
    if length < 0:
        outstr = outstr + (20 - len(outstr)) * " "
    if length > 0:
        outstr = (20 - len(outstr)) * " " + outstr
    if latex:
        if outstr.find("x10") > 0:
            outstr = outstr.replace("x10^", r"$\times$10$^{") + "}$"
    return outstr


def read_inffile(filename: str, verbose: bool = True) -> "infodata":
    """
    Read a PRESTO ``.inf`` file into an ``infodata`` C structure.

    Parameters
    ----------
    filename : str
        The ``.inf`` filename (a trailing ``.inf`` is optional).
    verbose : bool, optional
        If True, print what is being read (default True).

    Returns
    -------
    infodata
        The C ``infodata`` structure holding the file's contents.
    """
    fname = filename[:-4] if (filename[-4:] == ".inf") else filename
    id = infodata()
    if verbose:
        print("Reading information from", '"' + fname + '.inf"')
    readinf(id, fname)
    return id


def write_inffile(infodata, verbose: bool = True) -> None:
    """
    Write a PRESTO ``.inf`` file from an ``infodata`` C structure.

    Parameters
    ----------
    infodata : infodata
        The C ``infodata`` structure to write out.
    verbose : bool, optional
        If True, print the file being written (default True).

    Returns
    -------
    None
    """
    if verbose:
        print("Writing .inf file to '%s.inf'" % infodata.name)
    writeinf(infodata)


def psrepoch(psrname: str, epoch: float, verbose: bool = True) -> "psrparams":
    """
    Return a ``psrparams`` C structure for a pulsar at a given epoch.

    Parameters
    ----------
    psrname : str
        The B1950 or J2000 name of the pulsar (without a leading PSR, J, or
        B), which is looked up in the ATNF database; or the path to a
        parfile, from which the pulsar information is read instead.
    epoch : float
        The epoch in MJD format.
    verbose : bool, optional
        If True, print what was retrieved (default True).

    Returns
    -------
    psrparams
        The C ``psrparams`` structure for the pulsar at `epoch`.
    """
    pp = psrparams()
    if os.path.isfile(psrname):
        get_psr_from_parfile(psrname, epoch, pp)
        if verbose:
            print('Retrieved data at MJD %f from "%s"' % (epoch, psrname))
    else:
        num = get_psr_at_epoch(psrname, epoch, pp)
        if verbose:
            print("Retrieved data at MJD %f for %s" % (epoch, pp.jname))
            print("The pulsar was #%d in the database." % num)
    return pp


def read_rzwcands(filename: str) -> list:
    """
    Read all of the rzw search candidates from a file.

    Parameters
    ----------
    filename : str
        The candidate file to read.

    Returns
    -------
    list of fourierprops
        The rzw search candidates.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = fourierprops()
    while read_rzw_cand(infile, nextcand):
        cands.append(nextcand)
        nextcand = fourierprops()
    infile.close()
    return cands


def read_rawbincands(filename: str) -> list:
    """
    Read all of the raw binary search candidates from a file.

    Parameters
    ----------
    filename : str
        The candidate file to read.

    Returns
    -------
    list of rawbincand
        The raw binary search candidates.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = rawbincand()
    while read_rawbin_cand(infile, nextcand):
        cands.append(nextcand)
        nextcand = rawbincand()
    infile.close()
    return cands


def next2_to_n(x: float) -> int:
    """
    Return the smallest power of two that is >= `x`.

    Parameters
    ----------
    x : float
        The value to round up to a power of two.

    Returns
    -------
    int
        The first ``2**n`` that is >= `x`.
    """
    i = 1
    while i < x:
        i = i << 1
    return i


def rfft(data: np.ndarray, sign: int = -1) -> np.ndarray:
    """
    Return the FFT of real-valued, 32-bit floating point data.

    Only the positive-frequency half of the FFT is returned, since the other
    half is symmetric. The Nyquist frequency is stored in the imaginary part
    of frequency 0 as per Numerical Recipes.

    Parameters
    ----------
    data : numpy.ndarray
        The real-valued input data (forward) or complex data (inverse).
    sign : int, optional
        -1 for a forward transform (default), +1 for an inverse transform.

    Returns
    -------
    numpy.ndarray
        The (complex, for a forward transform) FFT of `data`.
    """
    # Default to sign = -1 if the user gives a bad value
    if sign == -1 or sign != 1:
        tmp = np.array(data, copy=1).astype(np.float32)
        realfft(tmp, -1)
        return tmp.view(np.complex64)
    else:
        tmp = np.array(data.view(np.float32), copy=1).astype(np.float32)
        realfft(tmp, 1)
        return tmp.view(np.float32)


def spectralpower(fftarray: np.ndarray) -> np.ndarray:
    """
    Return the power spectrum of a complex FFT.

    Parameters
    ----------
    fftarray : numpy.ndarray
        A complex FFT.

    Returns
    -------
    numpy.ndarray
        The power spectrum.
    """
    return power_arr(np.asarray(fftarray).astype(np.complex64))


def spectralphase(fftarray: np.ndarray) -> np.ndarray:
    """
    Return the spectral phase (in degrees) of a complex FFT.

    Parameters
    ----------
    fftarray : numpy.ndarray
        A complex FFT.

    Returns
    -------
    numpy.ndarray
        The spectral phase in degrees.
    """
    return phase_arr(np.asarray(fftarray).astype(np.complex64))


def rzw_response(
    roffset: float, z: float, w: float, numbetween: int = 1, numkern: int | None = None
) -> np.ndarray:
    """
    Return the response of a signal offset from a Fourier bin.

    Parameters
    ----------
    roffset : float
        The signal's offset from a Fourier bin, in bins.
    z : float
        The Fourier f-dot.
    w : float
        The Fourier f-dotdot.
    numbetween : int, optional
        The Fourier interpolation factor (default 1).
    numkern : int, optional
        The length of the resulting kernel. If None (default), it is
        auto-determined.

    Returns
    -------
    numpy.ndarray
        The complex response kernel.
    """
    if numkern is None:
        numkern = w_resp_halfwidth(z, w, LOWACC)
    return gen_w_response(roffset, numbetween, numkern, z, w)


def maximize_r(data: np.ndarray, r: float, norm: float | None = None) -> list:
    """
    Optimize the detection of a signal at Fourier frequency `r` in an FFT.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at which to optimize.
    norm : float, optional
        The normalization to divide the power by. If None (default), the
        local power (``rderivs.locpow``) is used.

    Returns
    -------
    list
        ``[maxpow, rmax, rderivs]``: the maximum normalized power, the
        optimized frequency, and the ``rderivs`` structure for the peak.
    """
    rd = rderivs()
    (rmax, maxpow) = max_r_arr(data, r, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, rd]


def maximize_rz(
    data: np.ndarray, r: float, z: float, norm: float | None = None
) -> list:
    """
    Optimize the detection of a signal at `r`, `z` in the F-Fdot plane.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at which to optimize.
    z : float
        The Fourier f-dot at which to optimize.
    norm : float, optional
        The normalization to divide the power by. If None (default), the
        local power (``rderivs.locpow``) is used.

    Returns
    -------
    list
        ``[maxpow, rmax, zmax, rderivs]``: the maximum normalized power, the
        optimized frequency and f-dot, and the ``rderivs`` structure for the
        peak.
    """
    rd = rderivs()
    (rmax, zmax, maxpow) = max_rz_arr(data, r, z, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, zmax, rd]


def maximize_rz_harmonics(
    data: np.ndarray, r: float, z: float, numharm: int, norm: float | None = None
) -> list:
    """
    Optimize a signal at `r`, `z` in the F-Fdot plane, summing harmonics.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at which to optimize.
    z : float
        The Fourier f-dot at which to optimize.
    numharm : int
        The number of harmonics to sum.
    norm : float, optional
        The normalization to divide each harmonic's power by. If None
        (default), each harmonic's local power is used.

    Returns
    -------
    list
        ``[maxpow, rmax, zmax, rds]``: the maximum summed normalized power,
        the optimized frequency and f-dot, and a list of ``rderivs``
        structures (one per harmonic).
    """
    rds = [rderivs() for ii in range(numharm)]
    derivdata = np.zeros(7 * numharm, dtype=np.float64)
    rmax, zmax = max_rz_arr_harmonics(data, r, z, derivdata)
    maxpow = 0.0
    for ii in range(numharm):
        rds[ii].pow = derivdata[ii * 7 + 0]
        rds[ii].phs = derivdata[ii * 7 + 1]
        rds[ii].dpow = derivdata[ii * 7 + 2]
        rds[ii].dphs = derivdata[ii * 7 + 3]
        rds[ii].d2pow = derivdata[ii * 7 + 4]
        rds[ii].d2phs = derivdata[ii * 7 + 5]
        rds[ii].locpow = derivdata[ii * 7 + 6]
        maxpow += rds[ii].pow / rds[ii].locpow if norm is None else rds[ii].pow / norm
    return [maxpow, rmax, zmax, rds]


def maximize_rzw(
    data: np.ndarray, r: float, z: float, w: float, norm: float | None = None
) -> list:
    """
    Optimize a signal at `r`, `z`, `w` in the F-Fdot-Fdotdot plane.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at which to optimize.
    z : float
        The Fourier f-dot at which to optimize.
    w : float
        The Fourier f-dotdot at which to optimize.
    norm : float, optional
        The normalization to divide the power by. If None (default), the
        local power (``rderivs.locpow``) is used.

    Returns
    -------
    list
        ``[maxpow, rmax, zmax, wmax, rderivs]``: the maximum normalized
        power, the optimized frequency, f-dot, and f-dotdot, and the
        ``rderivs`` structure for the peak.
    """
    rd = rderivs()
    (rmax, zmax, wmax, maxpow) = max_rzw_arr(data, r, z, w, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, zmax, wmax, rd]


def maximize_rzw_harmonics(
    data: np.ndarray,
    r: float,
    z: float,
    w: float,
    numharm: int,
    norm: float | None = None,
) -> list:
    """
    Optimize a signal at `r`, `z`, `w` in the F-Fd-Fdd volume, summing harmonics.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at which to optimize.
    z : float
        The Fourier f-dot at which to optimize.
    w : float
        The Fourier f-dotdot at which to optimize.
    numharm : int
        The number of harmonics to sum.
    norm : float, optional
        The normalization to divide each harmonic's power by. If None
        (default), each harmonic's local power is used.

    Returns
    -------
    list
        ``[maxpow, rmax, zmax, wmax, rds]``: the maximum summed normalized
        power, the optimized frequency, f-dot, and f-dotdot, and a list of
        ``rderivs`` structures (one per harmonic).
    """
    rds = [rderivs() for ii in range(numharm)]
    derivdata = np.zeros(7 * numharm, dtype=np.float64)
    rmax, zmax, wmax = max_rzw_arr_harmonics(data, r, z, w, derivdata)
    maxpow = 0.0
    for ii in range(numharm):
        rds[ii].pow = derivdata[ii * 7 + 0]
        rds[ii].phs = derivdata[ii * 7 + 1]
        rds[ii].dpow = derivdata[ii * 7 + 2]
        rds[ii].dphs = derivdata[ii * 7 + 3]
        rds[ii].d2pow = derivdata[ii * 7 + 4]
        rds[ii].d2phs = derivdata[ii * 7 + 5]
        rds[ii].locpow = derivdata[ii * 7 + 6]
        maxpow += rds[ii].pow / rds[ii].locpow if norm is None else rds[ii].pow / norm
    return [maxpow, rmax, zmax, wmax, rds]


def search_fft(data: np.ndarray, numcands: int, norm: float | str = "default") -> list:
    """
    Search a short FFT for its highest candidates.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    numcands : int
        The number of highest candidates to return.
    norm : float or str, optional
        The value to multiply each power by to normalize the power spectrum.
        Defaults to ``1.0 / (frequency-0 value)``.

    Returns
    -------
    list
        A list of ``[power, frequency]`` pairs for the `numcands` highest
        candidates.
    """
    if norm == "default":
        norm = 1.0 / data[0].real
    hp = np.zeros(numcands, "f")
    hf = np.zeros(numcands, "f")
    search_minifft(data, len(data), norm, numcands, hp, hf)
    cands = []
    for i in range(numcands):
        cands.append([hp[i], hf[i]])
    return cands


def ffdot_plane(
    data: np.ndarray,
    lor: int,
    dr: float,
    numr: int,
    loz: float,
    dz: float,
    numz: int,
) -> np.ndarray:
    """
    Generate an F-Fdot plane with its lower-left corner at `lor`, `loz`.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    lor : int
        The lowest Fourier frequency (bins) of the plane.
    dr : float
        The frequency step, which must be the reciprocal of an integer
        (``1 / numbetween``).
    numr : int
        The number of frequency bins.
    loz : float
        The lowest Fourier f-dot of the plane.
    dz : float
        The f-dot step.
    numz : int
        The number of f-dot slices.

    Returns
    -------
    numpy.ndarray
        The complex F-Fdot plane.

    Notes
    -----
    ``r`` is considered to be the average frequency (``r = r0 + z / 2``).
    """
    lor = int(lor)
    numr = int(numr)
    numz = int(numz)
    numbetween = int(1.0 / dr)
    hiz = loz + (numz - 1) * dz
    maxabsz = max(abs(loz), abs(hiz))
    kern_half_width = z_resp_halfwidth(maxabsz, LOWACC)
    fftlen = next2_to_n(numr + 2 * numbetween * kern_half_width)
    ffd = corr_rz_plane(data, numbetween, lor, loz, hiz, numz, fftlen, LOWACC)
    return np.array(ffd[:, 0:numr], copy=1)


def fdotdot_vol(
    data: np.ndarray,
    lor: int,
    dr: float,
    numr: int,
    loz: float,
    dz: float,
    numz: int,
    low: float,
    dw: float,
    numw: int,
) -> np.ndarray:
    """
    Generate an F-Fdot-Fdotdot volume with its corner at `lor`, `loz`, `low`.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    lor : int
        The lowest Fourier frequency (bins) of the volume.
    dr : float
        The frequency step, which must be the reciprocal of an integer
        (``1 / numbetween``).
    numr : int
        The number of frequency bins.
    loz : float
        The lowest Fourier f-dot of the volume.
    dz : float
        The f-dot step.
    numz : int
        The number of f-dot slices.
    low : float
        The lowest Fourier f-dotdot of the volume.
    dw : float
        The f-dotdot step.
    numw : int
        The number of f-dotdot slices.

    Returns
    -------
    numpy.ndarray
        The complex F-Fdot-Fdotdot volume.

    Notes
    -----
    ``r`` is the average frequency (``r = r0 + w/6 + z0/2``) and ``z`` is the
    average f-dot (``z = z0 + w / 2``).
    """
    lor = int(lor)
    numr, numz, numw = int(numr), int(numz), int(numw)
    numbetween = int(1.0 / dr)
    hiz = loz + (numz - 1) * dz
    maxabsz = max(abs(loz), abs(hiz))
    hiw = low + (numw - 1) * dw
    maxabsw = max(abs(low), abs(hiw))
    kern_half_width = w_resp_halfwidth(maxabsz, maxabsw, LOWACC)
    fftlen = next2_to_n(numr + 2 * numbetween * kern_half_width)
    ffd = corr_rzw_vol(
        data, numbetween, lor, loz, hiz, numz, low, hiw, numw, fftlen, LOWACC
    )
    return np.array(ffd[:, :, 0:numr], copy=1)


def estimate_rz(psr, T: float, show: int = 0, device: str = "/XWIN") -> tuple:
    """
    Estimate a pulsar's average Fourier frequency and f-dot.

    Parameters
    ----------
    psr : psrparams
        The ``psrparams`` structure describing the pulsar.
    T : float
        The length of the observation in sec.
    show : int, optional
        If nonzero, display plots of ``r`` and ``z`` (default 0).
    device : str, optional
        The PGPLOT device to plot to if `show` is nonzero (default "/XWIN").

    Returns
    -------
    tuple of float
        ``(r, z)``: the average Fourier frequency (relative to the nominal
        Fourier frequency) and the average Fourier f-dot, both in bins.
    """
    startE = keplers_eqn(psr.orb.t, psr.orb.p, psr.orb.e, 1.0e-15)
    numorbpts = int(T / psr.orb.p + 1.0) * 1024 + 1
    dt = T / (numorbpts - 1)
    E = dorbint(startE, numorbpts, dt, psr.orb)
    z = z_from_e(E, psr, T)
    r = T / p_from_e(E, psr) - T / psr.p
    if show:
        times = np.arange(numorbpts) * dt
        Pgplot.plotxy(
            r, times, labx="Time", laby="Fourier Frequency (r)", device=device
        )
        if device == "/XWIN":
            print("Press enter to continue:")
            input()
        Pgplot.nextplotpage()
        Pgplot.plotxy(
            z,
            times,
            labx="Time",
            laby="Fourier Frequency Derivative (z)",
            device=device,
        )
        Pgplot.closeplot()
    return r.mean(), z.mean()


def alias(r: float, rny: float) -> float:
    """
    Convert between an aliased and the 'true' Fourier frequency of a signal.

    The transformation is symmetric about the Nyquist frequency, so it
    converts either direction.

    Parameters
    ----------
    r : float
        The signal's Fourier frequency to convert.
    rny : float
        The Nyquist frequency in bins. For an FFT of real data, this is the
        number of data points FFT'd divided by 2.

    Returns
    -------
    float
        The converted Fourier frequency.
    """
    return 2.0 * rny - r


def show_ffdot_plane(
    data,
    r,
    z,
    dr=0.125,
    dz=0.5,
    numr=300,
    numz=300,
    T=None,
    contours=None,
    title=None,
    image="astro",
    device="/XWIN",
    norm=1.0,
):
    """
    Show a color plot of the F-Fdot plane centered on `r`, `z`.

    Parameters
    ----------
    data : numpy.ndarray
        The input FFT.
    r : float
        The Fourier frequency at the center of the plot.
    z : float
        The Fourier f-dot at the center of the plot.
    dr, dz : float, optional
        The frequency and f-dot steps (defaults 0.125 and 0.5).
    numr, numz : int, optional
        The number of frequency bins and f-dot slices (defaults 300, 300).
    T : float, optional
        The observation length in sec. If given, axes are labeled in Hz.
    contours : optional
        Contour levels to overplot.
    title : str, optional
        The plot title.
    image : str, optional
        The image color map (default "astro").
    device : str, optional
        The PGPLOT device (default "/XWIN").
    norm : float, optional
        The power normalization (default 1.0).

    Returns
    -------
    None
    """
    ffdp = ffdot_plane(data, r, dr, numr, z, dz, numz)
    ffdpow = spectralpower(ffdp.ravel())
    ffdpow = np.reshape(ffdpow, (numz, numr))
    startbin = int(r - (numr * dr) / 2)
    startz = int(z - (numz * dz) / 2)
    x = np.arange(numr, dtype="d") * dr + startbin
    y = np.arange(numz, dtype="d") * dz + startz
    highpt = np.argmax(ffdpow.ravel())
    hir = highpt % numr
    hiz = highpt // numr
    print("")
    print("Fourier Freqs from ", min(x), "to", max(x), ".")
    print("Fourier Fdots from ", min(y), "to", max(y), ".")
    print("Maximum normalized power is ", ffdpow[hiz][hir])
    print(
        "The max value is located at:  r =",
        startbin + hir * dr,
        "  z =",
        startz + hiz * dz,
    )
    print("")
    if not T:
        Pgplot.plot2d(
            ffdpow,
            x,
            y,
            labx="Fourier Frequency (bins)",
            laby="Fourier Frequency Derivative",
            title=title,
            image=image,
            contours=contours,
            device=device,
        )
    else:
        Pgplot.plot2d(
            ffdpow,
            x / T,
            y / (T**2.0),
            labx="Frequency (hz)",
            laby="Frequency Derivative (Hz/sec)",
            rangex2=[x[0], x[-1]],
            rangey2=[y[0], y[-1]],
            labx2="Fourier Frequency",
            laby2="Fourier Frequency Derivative",
            title=title,
            image=image,
            contours=contours,
            device=device,
        )


def v_from_e(e: np.ndarray, psr) -> np.ndarray:
    """
    Return velocities from a vector of eccentric anomalies.

    Parameters
    ----------
    e : numpy.ndarray
        The eccentric anomalies.
    psr : psrparams
        The ``psrparams`` instance describing the pulsar.

    Returns
    -------
    numpy.ndarray
        The velocities in km/s.
    """
    oldw = psr.orb.w
    v = np.array(e, copy=1)
    E_to_v(v, psr.orb)
    psr.orb.w = oldw
    return v


def d_from_e(e: np.ndarray, psr) -> np.ndarray:
    """
    Return time delays from a vector of eccentric anomalies.

    Parameters
    ----------
    e : numpy.ndarray
        The eccentric anomalies.
    psr : psrparams
        The ``psrparams`` instance describing the pulsar.

    Returns
    -------
    numpy.ndarray
        The time delays in seconds.
    """
    oldw = psr.orb.w
    d = np.array(e, copy=1)
    E_to_phib(d, psr.orb)
    psr.orb.w = oldw
    return d


def p_from_e(e: np.ndarray, psr) -> np.ndarray:
    """
    Return pulsar periods from a vector of eccentric anomalies.

    Parameters
    ----------
    e : numpy.ndarray
        The eccentric anomalies.
    psr : psrparams
        The ``psrparams`` instance describing the pulsar.

    Returns
    -------
    numpy.ndarray
        The pulsar periods in seconds.
    """
    oldw = psr.orb.w
    psr.orb.w = psr.orb.w * DEGTORAD
    p = np.array(e, copy=1)
    E_to_p(p, psr.p, psr.orb)
    psr.orb.w = oldw
    return p


def z_from_e(e: np.ndarray, psr, T: float) -> np.ndarray:
    """
    Return Fourier f-dots from a vector of eccentric anomalies.

    Parameters
    ----------
    e : numpy.ndarray
        The eccentric anomalies.
    psr : psrparams
        The ``psrparams`` instance describing the pulsar.
    T : float
        The total length of the observation in seconds.

    Returns
    -------
    numpy.ndarray
        The Fourier f-dots in bins.
    """
    oldw = psr.orb.w
    psr.orb.w = psr.orb.w * DEGTORAD
    z = np.array(e, copy=1)
    E_to_z(z, psr.p, T, psr.orb)
    psr.orb.w = oldw
    return z


def pcorr(
    data: np.ndarray, kernel: np.ndarray, numbetween: int, lo: int, hi: int
) -> np.ndarray:
    """
    Correlate the raw complex vectors `data` and `kernel`.

    Parameters
    ----------
    data : numpy.ndarray
        The raw complex data vector.
    kernel : numpy.ndarray
        The raw complex kernel vector.
    numbetween : int
        The Fourier interpolation factor.
    lo : int
        The starting frequency of the result (inclusive).
    hi : int
        The ending frequency of the result (exclusive).

    Returns
    -------
    numpy.ndarray
        The correlation result, from frequency `lo` up to (not including)
        `hi`.
    """
    kern_half_width = len(kernel) // (2 * numbetween)
    result = np.zeros((hi - lo) * numbetween, "F")
    corr_complex(
        data,
        len(data),
        RAW,
        kernel,
        len(kernel),
        RAW,
        result,
        len(result),
        lo,
        numbetween,
        kern_half_width,
        CORR,
    )
    return result


def p_to_f(p: float, pd: float, pdd: float) -> list:
    """
    Convert period and its derivatives to frequency and derivatives.

    The conversion is symmetric, so it also converts from f to p.

    Parameters
    ----------
    p : float
        The period (or frequency).
    pd : float
        The period derivative (or frequency derivative).
    pdd : float
        The period second derivative (or frequency second derivative).

    Returns
    -------
    list of float
        ``[f, fd, fdd]``, the equivalent frequency counterparts.
    """
    f = 1.0 / p
    fd = -pd / (p * p)
    if pdd == 0.0:
        fdd = 0.0
    else:
        fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
    return [f, fd, fdd]


def bary_to_topo(pb, pbd, pbdd, infofilenm, ephem="DE200"):
    """
    Calculate topocentric spin parameters from barycentric ones.

    Uses least squares to compute the topocentric period, period
    derivative, and period second derivative corresponding to the given
    barycentric values. The observation data must be found in the info file.

    Parameters
    ----------
    pb : float
        The barycentric period.
    pbd : float
        The barycentric period derivative.
    pbdd : float
        The barycentric period second derivative.
    infofilenm : str
        The ``.inf`` file describing the observation.
    ephem : str, optional
        The solar-system ephemeris to use (default "DE200").

    Returns
    -------
    list of float
        ``[pt, ptd, ptdd]``, the topocentric period and its derivatives.

    Warnings
    --------
    This routine is stale: it imports ``numpy.linalg.old`` (removed from
    NumPy long ago), so it does not currently run. It is kept pending a
    rewrite.
    """
    from numpy.linalg.old import linear_least_squares

    if infofilenm[-4:] == ".inf":
        infofilenm = infofilenm[:-4]
    obs = read_inffile(infofilenm)
    T = obs.N * obs.dt
    dt = 10.0
    tto = obs.mjd_i + obs.mjd_f
    tts = np.arange(tto, tto + (T + dt) / SECPERDAY, dt / SECPERDAY)
    nn = len(tts)
    bts = np.zeros(nn, "d")
    vel = np.zeros(nn, "d")
    ra = psr_utils.coord_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
    dec = psr_utils.coord_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
    try:
        tel = obscode(obs.telescope)
    except KeyError:
        print("Telescope not recognized.")
        return 0
    barycenter(tts, bts, vel, ra, dec, tel, ephem)
    print("Topocentric start time = %17.11f" % tts[0])
    print("Barycentric start time = %17.11f" % bts[0])
    avgvel = np.add.reduce(vel) / nn
    print("Average Earth velocity = %10.5e c" % (avgvel))
    tts = np.arange(nn, dtype="d") * dt
    bts = (bts - bts[0]) * SECPERDAY
    [fb, fbd, fbdd] = p_to_f(pb, pbd, pbdd)
    b = fb * bts + fbd * bts**2.0 / 2.0 + fbdd * bts**3.0 / 6.0
    a = np.transpose(np.asarray([tts, tts**2.0, tts**3.0]))
    [ft, ftd, ftdd], residuals, rank, sv = linear_least_squares(a, b)
    [pt, ptd, ptdd] = p_to_f(ft, ftd, ftdd)
    print("    Topocentric period = %15.12f" % pt)
    print("     Topocentric p-dot = %15.9e" % ptd)
    print("  Topocentric p-dotdot = %15.9e" % ptdd)
    print("     Quick Topo period = %15.12f" % (pb * (1.0 + avgvel)))
    print("      Quick Topo p-dot = %15.9e" % (pbd * (1.0 + avgvel)))
    print("   Quick Topo p-dotdot = %15.9e" % (pbdd * (1.0 + avgvel)))
    return [pt, ptd, ptdd]


def measure_phase(
    profile: np.ndarray, template: np.ndarray, sigma: float, fwhm: float
) -> tuple:
    """
    Measure the phase offset of a profile from a template.

    Uses the TOA measurement technique from J. H. Taylor's talk
    *Pulsar Timing and Relativistic Gravity*. The phase returned is cyclic
    (0-1).

    Parameters
    ----------
    profile : numpy.ndarray
        The measured profile.
    template : numpy.ndarray
        A high-S/N template. It must have the same number of points as
        `profile`.
    sigma : float
        The RMS noise level of `profile`.
    fwhm : float
        The approximate width of the template pulse (0-1).

    Returns
    -------
    tuple
        ``(tau, tau_err, b, b_err, a)``: the phase, its 1-sigma error, the
        scaling factor and its 1-sigma error, and the DC offset.
    """
    from presto.simple_roots import newton_raphson

    N = len(profile)
    if not (N == len(template)):
        print("Lengths of 'profile' and 'template' must")
        print("  be equal in measure_phase().")
        return 0.0
    ft = rfft(profile)
    p0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    P_k = abs(ft)
    frotate(P_k, len(ft), 1)
    Theta_k = np.arctan2(-ft.imag, ft.real)
    frotate(Theta_k, len(ft), 1)
    ft = rfft(template)
    s0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    S_k = abs(ft)
    frotate(S_k, len(ft), 1)
    Phi_k = np.arctan2(-ft.imag, ft.real)
    frotate(Phi_k, len(ft), 1)
    # Estimate of the noise sigma (This needs to be checked)
    # Note:  Checked 10 Jul 2000.  Looks OK.
    sig = sigma * np.sqrt(N)
    k = np.arange(len(ft), dtype="d") + 1.0

    def fn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
        # Since Nyquist freq always has phase = 0.0
        k[-1] = 0.0
        return np.add.reduce(k * p * s * np.sin(phi - theta + k * tau))

    def dfn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
        # Since Nyquist freq always has phase = 0.0
        k[-1] = 0.0
        return np.add.reduce(k * k * p * s * np.cos(phi - theta + k * tau))

    numphases = 200
    ddchidt = np.zeros(numphases, "d")
    phases = np.arange(numphases, dtype="d") / float(numphases - 1) * TWOPI - PI
    for i in np.arange(numphases):
        ddchidt[i] = dfn(phases[i])
    maxdphase = phases[np.argmax(ddchidt)] + 0.5 * TWOPI / (numphases - 1.0)
    # Solve for tau
    tau = newton_raphson(
        fn, dfn, maxdphase - 0.5 * fwhm * TWOPI, maxdphase + 0.5 * fwhm * TWOPI
    )
    # Solve for b
    c = P_k * S_k * np.cos(Phi_k - Theta_k + k * tau)
    d = np.add.reduce(S_k**2.0)
    b = np.add.reduce(c) / d
    # tau sigma
    tau_err = sig * np.sqrt(1.0 / (2.0 * b * np.add.reduce(k**2.0 * c)))
    # b sigma  (Note:  This seems to be an underestimate...)
    b_err = sig * np.sqrt(1.0 / (2.0 * d))
    # Solve for a
    a = (p0 - b * s0) / float(N)
    return (tau / TWOPI, tau_err / TWOPI, b, b_err, a)


def get_baryv(ra: str, dec: str, mjd: float, T: float, obs: str = "PK") -> float:
    """
    Determine the average barycentric velocity towards a sky position.

    Parameters
    ----------
    ra : str
        The right ascension in ``hh:mm:ss.ssss`` format.
    dec : str
        The declination in ``dd:mm:ss.ssss`` format.
    mjd : float
        The start of the observation in MJD.
    T : float
        The observation length in seconds.
    obs : str, optional
        The standard two-character ITOA observatory code (PK, GB, AO, GM,
        JB, ...; see ``src/observatories.c``). Default "PK".

    Returns
    -------
    float
        The average barycentric velocity towards `ra`, `dec` during the
        observation, in units of v/c.
    """
    tts = np.linspace(mjd, mjd + T / 86400.0, 100)
    nn = len(tts)
    bts = np.zeros(nn, dtype=np.float64)
    vel = np.zeros(nn, dtype=np.float64)
    barycenter(tts, bts, vel, ra, dec, obs, "DE421")
    return vel.mean()


def fold(
    indata: np.ndarray,
    dt: float,
    nbins: int,
    f: float,
    fd: float = 0.0,
    fdd: float = 0.0,
    startphs: float = 0.0,
    tlo: float = 0.0,
    standard: bool = True,
) -> tuple:
    """
    Fold data using PRESTO's ``fold()`` code (as ``prepfold`` does).

    Parameters
    ----------
    indata : numpy.ndarray
        The array of floats to fold.
    dt : float
        The duration in sec of each `indata` bin.
    nbins : int
        The number of bins in the output profile.
    f : float
        The frequency to fold at (Hz).
    fd : float, optional
        The frequency derivative to fold at (default 0.0).
    fdd : float, optional
        The frequency second derivative to fold at (default 0.0).
    startphs : float, optional
        The phase (0-1) for the beginning of the first bin (default 0.0).
    tlo : float, optional
        The time (sec) of the start of the first bin, with respect to the
        reference time of `f`, `fd`, and `fdd` (default 0.0).
    standard : bool, optional
        If True (default), use traditional prepfold "drizzling"; otherwise
        put each input data point fully into a single profile bin.

    Returns
    -------
    prof : numpy.ndarray
        The double-precision folded profile of length `nbins`.
    phs : float
        The ending phase (0-1) of the fold.
    """
    prof = np.zeros(nbins, dtype=np.float64)
    data = indata.astype(np.float32)
    phs = simplefold(data, dt, tlo, prof, startphs, f, fd, fdd, 1 if standard else 0)
    return (prof, phs)


def compute_chi2(data: np.ndarray, avg: float, var: float) -> float:
    """Compute chi^2 as a pulsation test for a folded pulse profile 'data'

    To get the reduced-chi^2, you would typically divide the result by
    the number of profile bins minus 1 (but beware of prepfold's inter-bin
    correlations!  See DOF_corr() in prepfold.py for details.)

    See Leahy et al. 1983 for details:
    https://ui.adsabs.harvard.edu/abs/1983ApJ...266..160L/abstract

    Parameters
    ----------
    data : [double precision numpy array]
        A folded pulse profile on which to compute Z^2_N
    avg : [double]
        The average level of the data (should be the background average).
    var : [double]
        The variance of the data (should be the background variance).
        Beware prepfold's bin correlations!
    """
    return chisqr(data, avg, var)


def compute_Z2N(data: np.ndarray, N: int, var: float) -> float:
    """Compute Z^2_N statistic for a folded pulse profile 'data'

    See Bachetti et al. 2021 for details:
    https://ui.adsabs.harvard.edu/abs/2021ApJ...909...33B/abstract

    Parameters
    ----------
    data : [double precision numpy array]
        A folded pulse profile on which to compute Z^2_N
    N : [integer]
        The number of harmonics to include in the Z^2_N calculation
    var : [double]
        The variance of the data (should be the background variance).
        Beware prepfold's bin correlations!
    """
    return z2n(data, var, N)
