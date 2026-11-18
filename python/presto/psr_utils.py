from __future__ import annotations

import bisect
from collections.abc import Callable

import numpy as np
import numpy.fft as FFT
from scipy.special import ndtr, ndtri, chdtrc, chdtri, fdtrc, i0, kolmogorov
from scipy.optimize import leastsq
from scipy.optimize import bisect as bisect_solve

from presto import Pgplot, ppgplot, sinc_interp
import presto.psr_constants as pc


def span(Min: float, Max: float, Number: int) -> np.ndarray:
    """
    Create a range of evenly spaced floats.

    Parameters
    ----------
    Min : float
        Inclusive lower end of the range.
    Max : float
        Inclusive upper end of the range.
    Number : int
        Number of points to generate.

    Returns
    -------
    numpy.ndarray
        `Number` floats spanning [`Min`, `Max`] inclusive.
    """
    return np.linspace(Min, Max, Number)


def distance(width: int) -> np.ndarray:
    """
    Build a square array of distances from its center.

    Parameters
    ----------
    width : int
        Side length of the returned square array.

    Returns
    -------
    numpy.ndarray
        A `width` x `width` array with each point set to the geometric
        distance from the array's center.
    """
    x = np.arange(-width / 2.0 + 0.5, width / 2.0 + 0.5, 1.0) ** 2
    x = np.resize(x, (width, width))
    return np.sqrt(x + np.transpose(x))


def is_power_of_10(n: float) -> bool:
    """
    Test whether a number is a power of 10.

    Parameters
    ----------
    n : float
        The number to test (converted to an int internally).

    Returns
    -------
    bool
        True if `n` is a power of 10, False otherwise.
    """
    N = int(n)
    while N > 9 and N % 10 == 0:
        N //= 10
    return N == 1


def choose_N(orig_N: int) -> int:
    """
    Choose a highly factorable time series length.

    The returned value is larger than `orig_N` but highly factorable.
    Note that it must be divisible by at least the maximum downsample
    factor * 2 (currently 8 * 2 = 16).

    Parameters
    ----------
    orig_N : int
        The original (minimum desired) time series length.

    Returns
    -------
    int
        A highly factorable length >= `orig_N`, or 0 if `orig_N` < 10000.
    """
    # A list of 4-dgit numbers that are highly factorable by small primes
    # fmt: off
    goodfactors = [1000, 1008, 1024, 1056, 1120, 1152, 1200, 1232, 1280,
                   1296, 1344, 1408, 1440, 1536, 1568, 1584, 1600, 1680,
                   1728, 1760, 1792, 1920, 1936, 2000, 2016, 2048, 2112,
                   2160, 2240, 2304, 2352, 2400, 2464, 2560, 2592, 2640,
                   2688, 2800, 2816, 2880, 3024, 3072, 3136, 3168, 3200,
                   3360, 3456, 3520, 3584, 3600, 3696, 3840, 3872, 3888,
                   3920, 4000, 4032, 4096, 4224, 4320, 4400, 4480, 4608,
                   4704, 4752, 4800, 4928, 5040, 5120, 5184, 5280, 5376,
                   5488, 5600, 5632, 5760, 5808, 6000, 6048, 6144, 6160,
                   6272, 6336, 6400, 6480, 6720, 6912, 7040, 7056, 7168,
                   7200, 7392, 7680, 7744, 7776, 7840, 7920, 8000, 8064,
                   8192, 8400, 8448, 8624, 8640, 8800, 8960, 9072, 9216,
                   9408, 9504, 9600, 9680, 9856, 10000]
    # fmt: on
    if orig_N < 10000:
        return 0
    # Get the number represented by the first 4 digits of orig_N
    first4 = int(str(orig_N)[:4])
    # Now get the number that is just bigger than orig_N
    # that has its first 4 digits equal to "factor"
    for factor in goodfactors:
        if (
            factor == first4
            and orig_N % factor == 0
            and is_power_of_10(orig_N // factor)
        ):
            break
        if factor > first4:
            break
    new_N = factor
    while new_N < orig_N:
        new_N *= 10
    if new_N == orig_N:
        return orig_N
    # Finally, compare new_N to the closest power_of_two
    # greater than orig_N.  Take the closest.
    two_N = 2
    while two_N < orig_N:
        two_N *= 2
    return min(two_N, new_N)


def running_avg(arr: np.ndarray, navg: int) -> np.ndarray:
    """
    Compute a non-overlapping running average.

    Parameters
    ----------
    arr : array_like
        The input array. Its length must be divisible by `navg`.
    navg : int
        The number of consecutive bins to average together.

    Returns
    -------
    numpy.ndarray
        The running average of `navg` bins from `arr`.
    """
    a = np.asarray(arr, "d")
    a = np.reshape(a, (len(a) // navg, navg))
    return np.add.reduce(np.transpose(a)) / navg


def hist(
    data, bins, range=None, laby: str = "Number", **kwargs
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return and plot a histogram in one variable.

    Parameters
    ----------
    data : array_like
        A sequence of data points.
    bins : int
        The number of bins into which the data is to be sorted.
    range : tuple of float, optional
        A tuple ``(lo, hi)`` specifying the lower and upper end of the
        interval spanned by the bins. Data points outside this interval
        are ignored. If None, the smallest and largest data values are
        used to define the interval.
    laby : str, optional
        The y-axis label for the plot (default "Number").
    **kwargs
        Additional keyword arguments passed to ``Pgplot.plotbinned``.

    Returns
    -------
    xs : numpy.ndarray
        The bin centers.
    ys : numpy.ndarray
        The number of data points in each bin.
    """
    ys, bin_edges = np.histogram(data, bins, range)
    dx = bin_edges[1] - bin_edges[0]
    xs = bin_edges[:-1] + 0.5 * dx
    maxy = int(1.1 * max(ys))
    if maxy < max(ys):
        maxy = max(ys) + 1.0
    if "rangey" not in list(kwargs.keys()):
        kwargs["rangey"] = [0, maxy]
    Pgplot.plotbinned(ys, xs, laby=laby, **kwargs)
    return (xs, ys)


def KS_test(
    data: np.ndarray, cumdist: Callable[[np.ndarray], np.ndarray], output: int = 0
) -> tuple[float, float]:
    """
    Perform a Kolmogorov-Smirnov test.

    Parameters
    ----------
    data : array_like
        The data points to test.
    cumdist : callable
        The cumulative-distribution function to compare the data against.
        It is called with a sorted copy of `data`.
    output : int, optional
        If nonzero, print the D and P values (default 0).

    Returns
    -------
    D : float
        The maximum distance between the cumulative distributions.
    P : float
        The probability that the data is drawn from the specified
        distribution.
    """
    nn = len(data)
    sdata = np.sort(np.asarray(data))
    D1 = np.maximum.reduce(np.absolute(cumdist(sdata) - np.arange(nn, dtype="d") / nn))
    D2 = np.maximum.reduce(
        np.absolute(cumdist(sdata) - np.arange(1, nn + 1, dtype="d") / nn)
    )
    D = max((D1, D2))
    P = kolmogorov(np.sqrt(nn) * D)
    if output:
        print("Max distance between the cumulative distributions (D) = %.5g" % D)
        print("Prob the data is from the specified distrbution   (P) = %.3g" % P)
    return (D, P)


def weighted_mean(
    arrin,
    weights_in,
    inputmean: float | None = None,
    calcerr: bool = False,
    sdev: bool = False,
) -> tuple[float, ...]:
    """
    Calculate the weighted mean, error, and optionally standard deviation.

    By default the error is calculated assuming the weights are 1/err^2,
    but if `calcerr` is True this assumption is dropped and the error is
    determined from the weighted scatter.

    Parameters
    ----------
    arrin : array_like
        A numpy array or a sequence that can be converted.
    weights_in : array_like
        A set of weights for each element in `arrin`.
    inputmean : float, optional
        An input mean value, around which the mean is calculated. If None,
        the weighted mean of `arrin` is used.
    calcerr : bool, optional
        If True, calculate the error as
        ``sqrt((w**2 * (arr - mean)**2).sum()) / weights.sum()``. By default
        (False) the error is ``1 / sqrt(weights.sum())``.
    sdev : bool, optional
        If True, also return the weighted standard deviation.

    Returns
    -------
    wmean : float
        The weighted mean.
    werr : float
        The weighted error.
    wsdev : float
        The weighted standard deviation. Only returned if `sdev` is True.

    Notes
    -----
    Converted from IDL: 2006-10-23, Erin Sheldon, NYU.
    """
    # no copy made if they are already arrays
    arr = np.array(arrin, ndmin=1)
    # Weights is forced to be type double. All resulting calculations
    # will also be double
    weights = np.array(weights_in, ndmin=1, dtype="f8")
    wtot = weights.sum()
    # user has input a mean value
    if inputmean is None:
        wmean = (weights * arr).sum() / wtot
    else:
        wmean = float(inputmean)
    # how should error be calculated?
    if calcerr:
        werr2 = (weights**2 * (arr - wmean) ** 2).sum()
        werr = np.sqrt(werr2) / wtot
    else:
        werr = 1.0 / np.sqrt(wtot)
    # should output include the weighted standard deviation?
    if sdev:
        wvar = (weights * (arr - wmean) ** 2).sum() / wtot
        wsdev = np.sqrt(wvar)
        return wmean, werr, wsdev
    else:
        return wmean, werr


def MJD_to_JD(MJD: float) -> float:
    """
    Convert Modified Julian Date (MJD) to Julian Date (JD).

    Parameters
    ----------
    MJD : float
        The Modified Julian Date.

    Returns
    -------
    float
        The equivalent Julian Date.
    """
    return MJD + 2400000.5


def JD_to_MJD(JD: float) -> float:
    """
    Convert Julian Date (JD) to Modified Julian Date (MJD).

    Parameters
    ----------
    JD : float
        The Julian Date.

    Returns
    -------
    float
        The equivalent Modified Julian Date.
    """
    return JD - 2400000.5


def MJD_to_Julian_Epoch(MJD: float) -> float:
    """
    Convert Modified Julian Date (MJD) to Julian Epoch.

    Parameters
    ----------
    MJD : float
        The Modified Julian Date.

    Returns
    -------
    float
        The equivalent Julian Epoch.
    """
    return 2000.0 + (MJD - 51544.5) / 365.25


def Julian_Epoch_to_MJD(jepoch: float) -> float:
    """
    Convert Julian Epoch to Modified Julian Date (MJD).

    Parameters
    ----------
    jepoch : float
        The Julian Epoch.

    Returns
    -------
    float
        The equivalent Modified Julian Date.
    """
    return 51544.5 + (jepoch - 2000.0) * 365.25


def MJD_to_Besselian_Epoch(MJD: float) -> float:
    """
    Convert Modified Julian Date (MJD) to Besselian Epoch.

    Parameters
    ----------
    MJD : float
        The Modified Julian Date.

    Returns
    -------
    float
        The equivalent Besselian Epoch.
    """
    return 1900.0 + (MJD - 15019.81352) / 365.242198781


def Besselian_Epoch_to_MJD(bepoch: float) -> float:
    """
    Convert Besselian Epoch to Modified Julian Date (MJD).

    Parameters
    ----------
    bepoch : float
        The Besselian Epoch.

    Returns
    -------
    float
        The equivalent Modified Julian Date.
    """
    return 15019.81352 + (bepoch - 1900.0) * 365.242198781


def rad_to_dms(rad: float) -> tuple[int, int, float]:
    """
    Convert radians to degrees, minutes, and seconds of arc.

    Parameters
    ----------
    rad : float
        The angle in radians.

    Returns
    -------
    tuple of (int, int, float)
        The angle as (degrees, minutes, seconds) of arc. The sign is
        carried on whichever leading component is nonzero.
    """
    if rad < 0.0:
        sign = -1
    else:
        sign = 1
    arc = pc.RADTODEG * np.fmod(np.fabs(rad), pc.PI)
    d = int(arc)
    arc = (arc - d) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    if sign == -1 and d == 0:
        return (sign * d, sign * m, sign * s)
    else:
        return (sign * d, m, s)


def dms_to_rad(deg: float, min: float, sec: float) -> float:
    """
    Convert degrees, minutes, and seconds of arc to radians.

    Parameters
    ----------
    deg : float
        Degrees of arc.
    min : float
        Minutes of arc.
    sec : float
        Seconds of arc.

    Returns
    -------
    float
        The equivalent angle in radians.
    """
    if deg < 0.0:
        sign = -1
    elif deg == 0.0 and (min < 0.0 or sec < 0.0):
        sign = -1
    else:
        sign = 1
    return (
        sign
        * pc.ARCSECTORAD
        * (60.0 * (60.0 * np.fabs(deg) + np.fabs(min)) + np.fabs(sec))
    )


def dms_to_deg(deg: float, min: float, sec: float) -> float:
    """
    Convert degrees, minutes, and seconds of arc to decimal degrees.

    Parameters
    ----------
    deg : float
        Degrees of arc.
    min : float
        Minutes of arc.
    sec : float
        Seconds of arc.

    Returns
    -------
    float
        The equivalent angle in decimal degrees.
    """
    return pc.RADTODEG * dms_to_rad(deg, min, sec)


def rad_to_hms(rad: float) -> tuple[int, int, float]:
    """
    Convert radians to hours, minutes, and seconds of arc.

    Parameters
    ----------
    rad : float
        The angle in radians.

    Returns
    -------
    tuple of (int, int, float)
        The angle as (hours, minutes, seconds).
    """
    rad = np.fmod(rad, pc.TWOPI)
    if rad < 0.0:
        rad = rad + pc.TWOPI
    arc = pc.RADTOHRS * rad
    h = int(arc)
    arc = (arc - h) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (h, m, s)


def hms_to_rad(hour: float, min: float, sec: float) -> float:
    """
    Convert hours, minutes, and seconds of arc to radians.

    Parameters
    ----------
    hour : float
        Hours.
    min : float
        Minutes.
    sec : float
        Seconds.

    Returns
    -------
    float
        The equivalent angle in radians.
    """
    if hour < 0.0:
        sign = -1
    else:
        sign = 1
    return (
        sign
        * pc.SECTORAD
        * (60.0 * (60.0 * np.fabs(hour) + np.fabs(min)) + np.fabs(sec))
    )


def hms_to_hrs(hour: float, min: float, sec: float) -> float:
    """
    Convert hours, minutes, and seconds of arc to decimal hours.

    Parameters
    ----------
    hour : float
        Hours.
    min : float
        Minutes.
    sec : float
        Seconds.

    Returns
    -------
    float
        The equivalent value in decimal hours.
    """
    return pc.RADTOHRS * hms_to_rad(hour, min, sec)


def coord_to_string(h_or_d: float, m: float, s: float) -> str:
    """
    Format an RA or DEC coordinate as a string.

    Parameters
    ----------
    h_or_d : float
        Hours (for RA) or degrees (for DEC).
    m : float
        Minutes of arc.
    s : float
        Seconds of arc.

    Returns
    -------
    str
        The coordinate as ``'hh:mm:ss.ssss'`` (RA) or ``'dd:mm:ss.ssss'``
        (DEC).
    """
    # The sign goes on the most-significant nonzero field, so values between 0 and
    # -1 deg arrive as (0, -m, -s) or (0, 0, -s).
    neg = h_or_d < 0 or (h_or_d == 0 and (m < 0.0 or s < 0.0))
    h_or_d, m, s = int(abs(h_or_d)), int(abs(m)), round(abs(s), 4)
    # Carry after rounding, so that e.g. 59.99996 s is the next minute, not 60.0000
    if s >= 60.0:
        s -= 60.0
        m += 1
    if m >= 60:
        m -= 60
        h_or_d += 1
    return "%s%.2d:%.2d:%07.4f" % ("-" if neg else "", h_or_d, m, s)


def ra_to_rad(ra_string: str) -> float:
    """
    Convert an RA string to radians.

    Parameters
    ----------
    ra_string : str
        RA information as ``'hh:mm:ss.ssss'``.

    Returns
    -------
    float
        The equivalent decimal radians.
    """
    h, m, s = ra_string.split(":")
    return hms_to_rad(int(h), int(m), float(s))


def dec_to_rad(dec_string: str) -> float:
    """
    Convert a DEC string to radians.

    Parameters
    ----------
    dec_string : str
        DEC information as ``'dd:mm:ss.ssss'``.

    Returns
    -------
    float
        The equivalent decimal radians.
    """
    d, m, s = dec_string.split(":")
    if "-" in d and int(d) == 0:
        m, s = "-" + m, "-" + s
    return dms_to_rad(int(d), int(m), float(s))


def delta_m(flux_factor: float) -> float:
    """
    Return the change in magnitudes caused by a change in flux.

    Parameters
    ----------
    flux_factor : float
        The multiplicative change in flux.

    Returns
    -------
    float
        The corresponding change in magnitudes.
    """
    return -2.5 * np.log10(flux_factor)


def flux_factor(delta_m: float) -> float:
    """
    Return the change in flux caused by a change in magnitude.

    Parameters
    ----------
    delta_m : float
        The change in magnitude.

    Returns
    -------
    float
        The corresponding multiplicative change in flux.
    """
    return 10.0 ** (delta_m / -2.5)


def distance_modulus_to_distance(dm: float, absorption: float = 0.0) -> float:
    """
    Convert a distance modulus to a distance.

    Parameters
    ----------
    dm : float
        The distance modulus.
    absorption : float, optional
        The absorption in magnitudes (default 0.0).

    Returns
    -------
    float
        The distance in kpc.
    """
    return 10.0 ** (((dm - absorption) + 5.0) / 5.0) / 1000.0


def distance_to_distance_modulus(d: float, absorption: float = 0.0) -> float:
    """
    Convert a distance to a distance modulus.

    Parameters
    ----------
    d : float
        The distance in kpc.
    absorption : float, optional
        The absorption in magnitudes (default 0.0).

    Returns
    -------
    float
        The distance modulus.
    """
    return 5.0 * np.log10(d * 1000.0) - 5.0 + absorption


def true_anomaly(E: float, ecc: float) -> float:
    """
    Return the true anomaly given the eccentric anomaly.

    Parameters
    ----------
    E : float
        The eccentric anomaly in radians.
    ecc : float
        The orbital eccentricity.

    Returns
    -------
    float
        The true anomaly in radians.
    """
    return 2.0 * np.arctan(np.sqrt((1.0 + ecc) / (1.0 - ecc)) * np.tan(E / 2.0))


def mass_funct(pb: float, x: float) -> float:
    """
    Return the mass function of an orbit.

    Parameters
    ----------
    pb : float
        The binary period in days.
    x : float
        The projected semi-major axis in lt-sec.

    Returns
    -------
    float
        The mass function in solar masses.
    """
    pbs = pb * pc.SECPERDAY
    return 8015123.37129 * x**3.0 / (pbs * pbs)


def mass_funct2(mp: float, mc: float, i: float) -> float:
    """
    Return the mass function of an orbit given the masses and inclination.

    Parameters
    ----------
    mp : float
        The mass of the primary in solar masses.
    mc : float
        The mass of the companion in solar masses.
    i : float
        The orbital inclination in radians.

    Returns
    -------
    float
        The mass function in solar masses.

    Notes
    -----
    An 'average' orbit has cos(i) = 0.5, or i = 60 deg.
    """
    return (mc * np.sin(i)) ** 3.0 / (mc + mp) ** 2.0


def asini_c(pb: float, mf: float) -> float:
    """
    Return the orbital projected semi-major axis.

    Parameters
    ----------
    pb : float
        The binary period in sec.
    mf : float
        The mass function of the orbit.

    Returns
    -------
    float
        The projected semi-major axis in lt-sec.
    """
    return (mf * pb * pb / 8015123.37129) ** (1.0 / 3.0)


def TS99_WDmass(pb: float, pop: str = "I+II") -> float | None:
    """
    Predict the WD companion mass for an MSP-He WD system.

    From Tauris & Savonije, 1999, ApJ.

    Parameters
    ----------
    pb : float
        The orbital period in days.
    pop : {"I", "I+II", "II"}, optional
        The population of the stars that formed the system (default
        "I+II"). Population II stars are older and more metal poor.

    Returns
    -------
    float or None
        The predicted WD companion mass in solar masses, or None if `pop`
        is not one of the valid options.
    """
    vals = {
        "I": (4.50, 1.2e5, 0.120),
        "I+II": (4.75, 1.1e5, 0.115),
        "II": (5.00, 1.0e5, 0.110),
    }
    if pop not in vals.keys():
        print("Not a valid stellar pop: should be 'I', 'I+II', or 'II'")
        return None
    else:
        a, b, c = vals[pop]
        return (pb / b) ** (1.0 / a) + c


def ELL1_check(
    A1: float, E: float, TRES: float, NTOA: int, output: bool = False
) -> bool:
    """
    Check whether the ELL1 binary model can be safely used.

    To work properly, we should have
    ``asini/c * ecc**2 << timing precision / sqrt(# TOAs)``,
    i.e. ``A1 * E**2 << TRES / sqrt(NTOA)``.

    Parameters
    ----------
    A1 : float
        The projected semi-major axis (asini/c) in lt-sec.
    E : float
        The orbital eccentricity.
    TRES : float
        The timing residual (timing precision) in microseconds.
    NTOA : int
        The number of TOAs.
    output : bool, optional
        If True, print diagnostic information (default False).

    Returns
    -------
    bool
        True if ELL1 should be fine or acceptable, False if BT or DD is
        recommended instead.
    """
    lhs = A1 * E**2.0 * 1e6
    rhs = TRES / np.sqrt(NTOA)
    if output:
        print(
            "Condition is asini/c * ecc**2 << timing precision / sqrt(# TOAs) to use ELL1:"
        )
        print("     asini/c * ecc**2 = %8.3g us" % lhs)
        print("  TRES / sqrt(# TOAs) = %8.3g us" % rhs)
    if lhs * 50.0 < rhs:
        if output:
            print("Should be fine.")
        return True
    elif lhs * 5.0 < rhs:
        if output:
            print("Should be OK, but not optimal.")
        return True
    else:
        if output:
            print("Should probably use BT or DD instead.")
        return False


def accel_to_z(accel: float, T: float, reffreq: float, harm: int = 1) -> float:
    """
    Convert an acceleration to the accelsearch 'z'.

    Parameters
    ----------
    accel : float
        The acceleration in m/s/s.
    T : float
        The observation duration in seconds.
    reffreq : float
        The reference frequency in Hz.
    harm : int, optional
        The harmonic number (default 1).

    Returns
    -------
    float
        The accelsearch 'z' (i.e. number of bins drifted).
    """
    return accel * harm * reffreq * T * T / pc.SOL


def z_to_accel(z: float, T: float, reffreq: float, harm: int = 1) -> float:
    """
    Convert the accelsearch 'z' to an acceleration.

    Parameters
    ----------
    z : float
        The accelsearch 'z' (i.e. number of bins drifted).
    T : float
        The observation duration in seconds.
    reffreq : float
        The reference frequency in Hz.
    harm : int, optional
        The harmonic number (default 1).

    Returns
    -------
    float
        The acceleration in m/s/s.
    """
    return z * pc.SOL / (harm * reffreq * T * T)


def bins_to_accel(
    z: float, T: float, f: list[float] = [1.0, 1000.0], device: str = "/XWIN"
) -> np.ndarray | None:
    """
    Plot the acceleration corresponding to a number of Fourier bins drifted.

    Parameters
    ----------
    z : float
        The number of Fourier bins drifted during the observation.
    T : float
        The observation length in seconds.
    f : list of float, optional
        The frequency range ``[fmin, fmax]`` in Hz (default [1.0, 1000.0]).
    device : str, optional
        The PGPLOT device (default "/XWIN"). If falsy, the accelerations
        are returned instead of plotted.

    Returns
    -------
    numpy.ndarray or None
        The accelerations if `device` is falsy, otherwise None (a plot is
        produced).
    """
    fs = span(np.log10(f[0]), np.log10(f[1]), 1000)
    accels = z_to_accel(z, T, 10.0**fs)
    if device:
        Pgplot.plotxy(
            np.log10(accels),
            fs,
            logx=1,
            logy=1,
            labx="Frequency (Hz)",
            laby=r"Acceleration (m/s\u2\d)",
            device=device,
        )
        ppgplot.pgmtxt("T", -2.0, 0.75, 0.0, "T = %.0f sec" % T)
        ppgplot.pgmtxt("T", -3.5, 0.75, 0.0, r"r\B\u\.\d = %.1f bins" % z)
        if device != "/XWIN":
            Pgplot.closeplot()
    else:
        return accels


def pulsar_mass(pb: float, x: float, mc: float, inc: float) -> float:
    """
    Return the pulsar mass for a binary system.

    Parameters
    ----------
    pb : float
        The binary period in days.
    x : float
        The projected semi-major axis in lt-sec.
    mc : float
        The mass of the companion in solar mass units.
    inc : float
        The orbital inclination in degrees.

    Returns
    -------
    float
        The pulsar mass in solar mass units.
    """
    massfunct = mass_funct(pb, x)

    def localmf(mp, mc=mc, mf=massfunct, i=inc * pc.DEGTORAD):
        return mass_funct2(mp, mc, i) - mf

    return bisect_solve(localmf, 0.0, 1000.0)


def companion_mass(pb: float, x: float, inc: float = 60.0, mpsr: float = 1.4) -> float:
    """
    Return the companion mass for a binary system.

    Parameters
    ----------
    pb : float
        The binary period in days.
    x : float
        The projected semi-major axis in lt-sec.
    inc : float, optional
        The orbital inclination in degrees (default 60.0).
    mpsr : float, optional
        The mass of the pulsar in solar mass units (default 1.4).

    Returns
    -------
    float
        The companion mass in solar mass units.
    """
    massfunct = mass_funct(pb, x)

    def localmf(mc, mp=mpsr, mf=massfunct, i=inc * pc.DEGTORAD):
        return mass_funct2(mp, mc, i) - mf

    return bisect_solve(localmf, 0.0, 1000.0)


def companion_mass_limit(pb: float, x: float, mpsr: float = 1.4) -> float:
    """
    Return the lower limit of the companion mass in a binary system.

    The lower limit corresponds to an inclination of i = 90 degrees.

    Parameters
    ----------
    pb : float
        The binary period in days.
    x : float
        The projected semi-major axis in lt-sec.
    mpsr : float, optional
        The mass of the pulsar in solar mass units (default 1.4).

    Returns
    -------
    float
        The minimum companion mass in solar mass units.
    """
    return companion_mass(pb, x, inc=90.0, mpsr=mpsr)


def OMDOT(porb: float, e: float, Mp: float, Mc: float) -> float:
    """
    Return the predicted advance of periastron.

    Parameters
    ----------
    porb : float
        The orbital period in days.
    e : float
        The orbital eccentricity.
    Mp : float
        The pulsar mass in solar masses.
    Mc : float
        The companion mass in solar masses.

    Returns
    -------
    float
        The predicted advance of periastron in deg/yr.
    """
    return (
        3.0
        * (porb * pc.SECPERDAY / pc.TWOPI) ** (-5.0 / 3.0)
        * (pc.Tsun * (Mp + Mc)) ** (2.0 / 3.0)
        / (1.0 - e**2.0)
        * pc.RADTODEG
        * pc.SECPERJULYR
    )


def GAMMA(porb: float, e: float, Mp: float, Mc: float) -> float:
    """
    Return the predicted value of relativistic gamma.

    Parameters
    ----------
    porb : float
        The orbital period in days.
    e : float
        The orbital eccentricity.
    Mp : float
        The pulsar mass in solar masses.
    Mc : float
        The companion mass in solar masses.

    Returns
    -------
    float
        The predicted relativistic gamma in sec.
    """
    return (
        e
        * (porb * pc.SECPERDAY / pc.TWOPI) ** (1.0 / 3.0)
        * pc.Tsun ** (2.0 / 3.0)
        * (Mp + Mc) ** (-4.0 / 3.0)
        * Mc
        * (Mp + 2.0 * Mc)
    )


def PBDOT(porb: float, e: float, Mp: float, Mc: float) -> float:
    """
    Return the predicted orbital period derivative.

    Parameters
    ----------
    porb : float
        The orbital period in days.
    e : float
        The orbital eccentricity.
    Mp : float
        The pulsar mass in solar masses.
    Mc : float
        The companion mass in solar masses.

    Returns
    -------
    float
        The predicted orbital period derivative in s/s.
    """
    return (
        -192.0
        * pc.PI
        / 5.0
        * (porb * pc.SECPERDAY / pc.TWOPI) ** (-5.0 / 3.0)
        * (1.0 + 73.0 / 24.0 * e**2.0 + 37.0 / 96.0 * e**4.0)
        * (1.0 - e**2.0) ** (-7.0 / 2.0)
        * pc.Tsun ** (5.0 / 3.0)
        * Mp
        * Mc
        * (Mp + Mc) ** (-1.0 / 3.0)
    )


def OMDOT_to_Mtot(OMDOT: float, porb: float, e: float) -> float:
    """
    Return the total system mass given an advance of periastron.

    Parameters
    ----------
    OMDOT : float
        The advance of periastron in deg/yr.
    porb : float
        The orbital period in days.
    e : float
        The orbital eccentricity.

    Returns
    -------
    float
        The total mass of the system in solar units.
    """
    wd = OMDOT / pc.SECPERJULYR * pc.DEGTORAD  # rad/s
    return (
        wd / 3.0 * (1.0 - e * e) * (porb * pc.SECPERDAY / pc.TWOPI) ** (5.0 / 3.0)
    ) ** (3.0 / 2.0) / pc.Tsun


def GAMMA_to_Mc(gamma: float, porb: float, e: float, Mp: float) -> float:
    """
    Return the predicted companion mass given relativistic gamma.

    Parameters
    ----------
    gamma : float
        The relativistic gamma in sec.
    porb : float
        The orbital period in days.
    e : float
        The orbital eccentricity.
    Mp : float
        The pulsar mass in solar units.

    Returns
    -------
    float
        The predicted companion mass in solar units.
    """

    def funct(mc, mp=Mp, porb=porb, e=e, gamma=gamma):
        return GAMMA(porb, e, mp, mc) - gamma

    return bisect_solve(funct, 0.01, 20.0)


def shklovskii_effect(pm: float, D: float) -> float:
    """
    Return the apparent acceleration due to the Shklovskii effect.

    Parameters
    ----------
    pm : float
        The proper motion in mas/yr.
    D : float
        The distance in kpc.

    Returns
    -------
    float
        The 'acceleration' due to the transverse Doppler effect, a_pm/C,
        equivalently Pdot_pm/P.
    """
    return (
        (pm / 1000.0 * pc.ARCSECTORAD / pc.SECPERJULYR) ** 2.0
        * pc.KMPERKPC
        * D
        / (pc.C / 1000.0)
    )


def galactic_accel_simple(
    l: float, b: float, D: float, v_o: float = 240.0, R_o: float = 8.34
) -> float:
    """
    Return the approximate projected Galactic acceleration/c (simple model).

    This is ``(a_p - a_ssb) dot n / c``, where a_p and a_ssb are
    acceleration vectors and n is the line-of-sight vector. It assumes a
    simple spherically symmetric isothermal sphere. This is eqn 2.4 of
    Phinney 1992.

    Parameters
    ----------
    l : float
        The Galactic longitude in degrees.
    b : float
        The Galactic latitude in degrees.
    D : float
        The distance in kpc.
    v_o : float, optional
        The circular velocity in km/s (default 240.0, from Reid et al 2014).
    R_o : float, optional
        The distance to the Galactic center in kpc (default 8.34, from Reid
        et al 2014).

    Returns
    -------
    float
        The projected acceleration/c in s^-1.
    """
    A_sun = v_o * v_o / (pc.C / 1000.0 * R_o * pc.KMPERKPC)
    d = D / R_o
    cbcl = np.cos(b * pc.DEGTORAD) * np.cos(l * pc.DEGTORAD)
    return -A_sun * (cbcl + (d - cbcl) / (1.0 + d * d - 2.0 * d * cbcl))


def galactic_accel(
    l: float, b: float, D: float, v_o: float = 240.0, R_o: float = 8.34
) -> float:
    """
    Return the approximate projected Galactic acceleration/c.

    This is ``(a_p - a_ssb) dot n / c``, where a_p and a_ssb are
    acceleration vectors and n is the line-of-sight vector. This is eqn 5
    of Nice & Taylor 1995.

    Parameters
    ----------
    l : float
        The Galactic longitude in degrees.
    b : float
        The Galactic latitude in degrees.
    D : float
        The distance in kpc.
    v_o : float, optional
        The circular velocity in km/s (default 240.0, from Reid et al 2014).
    R_o : float, optional
        The distance to the Galactic center in kpc (default 8.34, from Reid
        et al 2014).

    Returns
    -------
    float
        The projected acceleration/c in s^-1.
    """
    A_sun = v_o * v_o / (pc.C / 1000.0 * R_o * pc.KMPERKPC)
    cb = np.cos(b * pc.DEGTORAD)
    cl = np.cos(l * pc.DEGTORAD)
    sl = np.sin(l * pc.DEGTORAD)
    beta = D / R_o * cb - cl
    return -A_sun * cb * (cl + beta / (sl**2 + beta**2))


def gal_z_accel(l: float, b: float, D: float) -> float:
    """
    Return the approximate projected acceleration/c towards the Galactic plane.

    This is ``(a_p - a_ssb) dot n / c`` caused by the acceleration of the
    pulsar towards the plane of the Galaxy. This is eqn 3+4 of Nice &
    Taylor 1995.

    Parameters
    ----------
    l : float
        The Galactic longitude in degrees.
    b : float
        The Galactic latitude in degrees.
    D : float
        The distance in kpc.

    Returns
    -------
    float
        The projected acceleration/c in s^-1.
    """
    sb = np.sin(b * pc.DEGTORAD)
    z = D * sb
    az = 1.08e-19 * (1.25 * z / np.sqrt(z**2 + 0.0324) + 0.58 * z)
    return az * sb


def beam_halfwidth(obs_freq: float, dish_diam: float) -> float:
    """
    Return the telescope beam halfwidth.

    Parameters
    ----------
    obs_freq : float
        The observing frequency in MHz.
    dish_diam : float
        The telescope diameter in m.

    Returns
    -------
    float
        The beam halfwidth in arcmin.
    """
    return 1.2 * pc.SOL / (obs_freq * 10.0**6) / dish_diam * pc.RADTODEG * 60 / 2


def limiting_flux_dens(
    Ttot: float,
    G: float,
    BW: float,
    T: float,
    P: float = 0.01,
    W: float = 0.05,
    polar: int = 2,
    factor: float = 15.0,
) -> float:
    """
    Return the approximate limiting flux density for a pulsar survey.

    This is a *very* approximate calculation. For a better calculation,
    see Cordes and Chernoff, ApJ, 482, p971, App. A.

    Parameters
    ----------
    Ttot : float
        The sky + system temperature in K.
    G : float
        The forward gain of the antenna in K/Jy.
    BW : float
        The observing bandwidth in MHz.
    T : float
        The integration time in sec.
    P : float, optional
        The pulsar period in sec (default 0.01).
    W : float, optional
        The duty cycle of the pulsar, 0-1 (default 0.05).
    polar : int, optional
        The number of polarizations (default 2).
    factor : float, optional
        Normalization factor accounting for limiting SNR, hardware
        limitations, etc. (default 15.0).

    Returns
    -------
    float
        The approximate limiting flux density in mJy.

    Notes
    -----
    Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy.
    """
    w = W * P
    return np.sqrt(w / ((P - w) * polar * BW * T)) * factor * Ttot / G


def dm_info(
    dm: float | None = None,
    dmstep: float = 1.0,
    freq: float = 1390.0,
    numchan: int = 512,
    chanwidth: float = 0.5,
) -> None:
    """
    Print info about potential DM smearing during an observation.

    Parameters
    ----------
    dm : float, optional
        The dispersion measure in pc cm^-3. If given, the per-channel
        smearing is also printed.
    dmstep : float, optional
        The DM step size in pc cm^-3 (default 1.0).
    freq : float, optional
        The center frequency in MHz (default 1390.0).
    numchan : int, optional
        The number of channels (default 512).
    chanwidth : float, optional
        The channel width in MHz (default 0.5).

    Returns
    -------
    None
    """
    BW = chanwidth * numchan
    print("      Center freq (MHz) = %.3f" % (freq))
    print("     Number of channels = %d" % (numchan))
    print("    Channel width (MHz) = %.3g" % (chanwidth))
    print("  Total bandwidth (MHz) = %.3g" % (BW))
    print("   DM offset (0.5*step) = %.3g" % (0.5 * dmstep))
    print(
        "  Smearing over BW (ms) = %.3g" % (1000.0 * dm_smear(0.5 * dmstep, BW, freq))
    )
    if dm:
        print(
            " Smearing per chan (ms) = %.3g" % (1000.0 * dm_smear(dm, chanwidth, freq))
        )


def best_dm_step(
    maxsmear: float = 0.1,
    dt: float = 0.00080,
    dm: float = 0.0,
    freq: float = 1390.0,
    numchan: int = 512,
    chanwidth: float = 0.5,
) -> float:
    """
    Return the DM step needed to keep total smearing below a threshold.

    Parameters
    ----------
    maxsmear : float, optional
        The maximum total smearing in ms (default 0.1).
    dt : float, optional
        The sample time in sec (default 0.00080).
    dm : float, optional
        The dispersion measure in pc cm^-3 (default 0.0).
    freq : float, optional
        The center frequency in MHz (default 1390.0).
    numchan : int, optional
        The number of channels (default 512).
    chanwidth : float, optional
        The channel width in MHz (default 0.5).

    Returns
    -------
    float
        The required DM step in pc cm^-3, or 0.0 if the requested total
        smearing is smaller than the fixed smearing components.
    """
    BW = chanwidth * numchan
    tau_tot = maxsmear / 1000.0
    tau_chan = dm_smear(dm, chanwidth, freq)
    tau_samp = dt
    if tau_tot**2.0 < (tau_chan**2.0 + tau_samp**2.0):
        print(
            "The requested total smearing is smaller than one or more of the components."
        )
        return 0.0
    else:
        return (
            0.0001205
            * freq**3.0
            * 2.0
            / BW
            * np.sqrt(tau_tot**2.0 - tau_chan**2.0 - tau_samp**2.0)
        )


def dm_smear_approx(dm: float, BW: float, center_freq: float) -> float:
    """
    Return the DM smearing over a bandwidth (small-BW approximation).

    This version assumes the bandwidth is small compared to `center_freq`.

    Parameters
    ----------
    dm : float
        The dispersion measure in pc cm^-3.
    BW : float
        The bandwidth in MHz.
    center_freq : float
        The center frequency in MHz.

    Returns
    -------
    float
        The smearing in sec.
    """
    return dm * BW / (0.0001205 * center_freq * center_freq * center_freq)


def dm_smear(dm: float, BW: float, center_freq: float) -> float:
    """
    Return the DM smearing over a bandwidth.

    Parameters
    ----------
    dm : float
        The dispersion measure in pc cm^-3.
    BW : float
        The bandwidth in MHz.
    center_freq : float
        The center frequency in MHz.

    Returns
    -------
    float
        The smearing in sec.
    """
    return delay_from_DM(dm, center_freq - 0.5 * BW) - delay_from_DM(
        dm, center_freq + 0.5 * BW
    )


def diagonal_DM(dt: float, chanBW: float, center_freq: float) -> float:
    """
    Return the "diagonal DM".

    The diagonal DM is the DM for which the smearing across one channel is
    equal to the sample time.

    Parameters
    ----------
    dt : float
        The sample time in sec.
    chanBW : float
        The channel bandwidth in MHz.
    center_freq : float
        The center frequency in MHz.

    Returns
    -------
    float
        The diagonal DM in pc cm^-3.
    """
    return (0.0001205 * center_freq * center_freq * center_freq) * dt / chanBW


def pulse_broadening(DM: float, f_ctr: float) -> float:
    """
    Return the approximate scattering pulse broadening.

    Based on the rough relation in Cordes' 'Pulsar Observations I' paper.
    The approximate error is 0.65 in log(tau).

    Parameters
    ----------
    DM : float
        The dispersion measure in pc cm^-3.
    f_ctr : float
        The center frequency in MHz.

    Returns
    -------
    float
        The approximate pulse broadening (tau) in ms.
    """
    logDM = np.log10(DM)
    return (
        10.0
        ** (-3.59 + 0.129 * logDM + 1.02 * logDM**2.0 - 4.4 * np.log10(f_ctr / 1000.0))
        / 1000.0
    )


def rrat_period(times: np.ndarray, numperiods: int = 20, output: bool = True) -> float:
    """
    Determine a RRAT pulse period by brute-force search.

    Parameters
    ----------
    times : array_like
        The (real!) single-pulse arrival times.
    numperiods : int, optional
        The number of integer pulses to try between the first two pulses
        (default 20).
    output : bool, optional
        If True, print diagnostic information (default True).

    Returns
    -------
    float
        The most likely (refined) RRAT period, in the same units as `times`.
    """
    ts = np.asarray(sorted(times))
    ps = (ts[1] - ts[0]) / np.arange(1, numperiods + 1)
    dts = np.diff(ts)
    xs = dts / ps[:, np.newaxis]
    metric = np.sum(np.fabs((xs - xs.round())), axis=1)
    pnum = metric.argmin()
    numrots = xs.round()[pnum].sum()
    p = (ts[-1] - ts[0]) / numrots
    if output:
        print(
            "Min, avg, std metric values are %.4f, %.4f, %.4f"
            % (metric.min(), metric.mean(), metric.std())
        )
        print(" Approx period is likely:", ps[pnum])
        print("Refined period is likely:", p)
        print("Rotations between pulses are:")
        print(dts / p)
    return p


def rrat_period_multiday(
    days_times: list, numperiods: int = 20, output: bool = True
) -> float:
    """
    Determine a RRAT pulse period from multi-day data by brute-force search.

    Parameters
    ----------
    days_times : list of array_like
        A list where each entry is the list of single-pulse arrival times
        for a single day/observation, e.g.
        ``[[times, from, one, day], [times, from, another, day], ...]``.
    numperiods : int, optional
        The maximum number of periods to try in the smallest interval
        between pulses (default 20).
    output : bool, optional
        If True, print diagnostic information (default True).

    Returns
    -------
    float
        The most likely (refined) RRAT period, in the same units as the
        input times.
    """
    all_dt = []
    for times in days_times:
        daily_dt = np.diff(sorted(times))
        all_dt.extend(daily_dt.tolist())

    dts = np.asarray(sorted(all_dt))
    ps = dts[0] / np.arange(1, numperiods + 1)
    xs = dts / ps[:, np.newaxis]
    metric = np.sum(np.fabs((xs - xs.round())), axis=1)
    pnum = metric.argmin()

    numrots = xs.round()[pnum].sum()
    p = dts.sum() / numrots

    if output:
        print(
            "Min, avg, std metric values are %.4f, %.4f, %.4f"
            % (metric.min(), metric.mean(), metric.std())
        )
        print(" Approx period is likely:", ps[pnum])
        print("Refined period is likely:", p)
        print("Rotations between pulses are:")
        print(dts / p)
    return p


def guess_DMstep(DM: float, dt: float, BW: float, f_ctr: float) -> float:
    """
    Choose a reasonable DM step.

    The step is set so that the maximum smearing across the bandwidth
    equals the sampling time.

    Parameters
    ----------
    DM : float
        The dispersion measure in pc cm^-3.
    dt : float
        The sampling time in sec.
    BW : float
        The bandwidth in MHz.
    f_ctr : float
        The center frequency in MHz.

    Returns
    -------
    float
        A reasonable DM step in pc cm^-3.
    """
    return dt * 0.0001205 * f_ctr**3.0 / (0.5 * BW)


def delay_from_DM(DM: float, freq_emitted: float | np.ndarray) -> float | np.ndarray:
    """
    Return the dispersive delay.

    Parameters
    ----------
    DM : float
        The dispersion measure in cm^-3 pc.
    freq_emitted : float or numpy.ndarray
        The emitted frequency (or frequencies) of the pulsar in MHz.

    Returns
    -------
    float or numpy.ndarray
        The delay in seconds caused by dispersion. Zero where
        `freq_emitted` is not positive.
    """
    if type(freq_emitted) is type(0.0):
        if freq_emitted > 0.0:
            return DM / (0.000241 * freq_emitted * freq_emitted)
        else:
            return 0.0
    else:
        return np.where(
            freq_emitted > 0.0, DM / (0.000241 * freq_emitted * freq_emitted), 0.0
        )


def delay_from_foffsets(
    df: float, dfd: float, dfdd: float, times: np.ndarray
) -> np.ndarray:
    """
    Return the phase delays caused by frequency offsets.

    Parameters
    ----------
    df : float
        The offset in frequency.
    dfd : float
        The offset in the first frequency derivative.
    dfdd : float
        The offset in the second frequency derivative.
    times : numpy.ndarray
        The times in seconds at which to evaluate the delays.

    Returns
    -------
    numpy.ndarray
        The delays in phase.
    """
    f_delays = df * times
    fd_delays = dfd * times**2 / 2.0
    fdd_delays = dfdd * times**3 / 6.0
    return f_delays + fd_delays + fdd_delays


def smear_plot(
    dm: list[float] = [1.0, 1000.0],
    dmstep: float = 1.0,
    subdmstep: float = 10.0,
    freq: float = 1390.0,
    numchan: int = 512,
    numsub: int = 32,
    chanwidth: float = 0.5,
    dt: float = 0.000125,
    device: str = "/xwin",
) -> None:
    """
    Plot the expected smearing from various effects in a radio pulsar search.

    Parameters
    ----------
    dm : list of float, optional
        The DM range ``[dm_min, dm_max]`` in pc cm^-3 (default [1.0, 1000.0]).
    dmstep : float, optional
        The DM step size in pc cm^-3 (default 1.0).
    subdmstep : float, optional
        The subband DM step size in pc cm^-3 (default 10.0).
    freq : float, optional
        The center frequency in MHz (default 1390.0).
    numchan : int, optional
        The number of channels (default 512).
    numsub : int, optional
        The number of subbands (default 32).
    chanwidth : float, optional
        The channel width in MHz (default 0.5).
    dt : float, optional
        The sample time in sec (default 0.000125).
    device : str, optional
        The PGPLOT device (default "/xwin").

    Returns
    -------
    None
    """
    numpts = 500
    BW = numchan * chanwidth
    subBW = numchan / numsub * chanwidth
    maxDMerror = 0.5 * dmstep
    maxsubDMerror = 0.5 * subdmstep
    ldms = span(np.log10(dm[0]), np.log10(dm[1]), numpts)
    dms = 10.0**ldms
    # Smearing from sample rate
    dts = np.zeros(numpts) + 1000.0 * dt
    # Smearing due to the intrinsic channel width
    chan_smear = 1000.0 * dm_smear(dms, chanwidth, freq)
    # Smearing across the full BW due to max DM mismatch
    BW_smear = np.zeros(numpts) + 1000.0 * dm_smear(maxDMerror, BW, freq)
    # Smearing in each subband due to max DM mismatch
    subband_smear = np.zeros(numpts) + 1000.0 * dm_smear(maxsubDMerror, subBW, freq)
    total_smear = np.sqrt(
        dts**2.0 + chan_smear**2.0 + subband_smear**2.0 + BW_smear**2.0
    )
    maxval = np.log10(2.0 * max(total_smear))
    minval = np.log10(
        0.5 * min([min(dts), min(chan_smear), min(BW_smear), min(subband_smear)])
    )
    Pgplot.plotxy(
        np.log10(total_smear),
        ldms,
        rangey=[minval, maxval],
        logx=1,
        logy=1,
        labx="Dispersion Measure",
        laby="Smearing (ms)",
        device=device,
    )
    ppgplot.pgsch(0.8)
    ppgplot.pgmtxt("t", 1.5, 1.0 / 12.0, 0.5, r"\(2156)\dcenter\u = %gMHz" % freq)
    ppgplot.pgmtxt("t", 1.5, 3.0 / 12.0, 0.5, r"N\dchan\u = %d" % numchan)
    ppgplot.pgmtxt("t", 1.5, 5.0 / 12.0, 0.5, r"N\dsub\u = %d" % numsub)
    ppgplot.pgmtxt("t", 1.5, 7.0 / 12.0, 0.5, r"BW\dchan\u = %gMHz" % chanwidth)
    ppgplot.pgmtxt("t", 1.5, 9.0 / 12.0, 0.5, r"\gDDM = %g" % dmstep)
    ppgplot.pgmtxt("t", 1.5, 11.0 / 12.0, 0.5, r"\gDDM\dsub\u = %g" % subdmstep)
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt("b", -7.5, 0.95, 1.0, "Total")
    Pgplot.plotxy(np.log10(dts), ldms, color="green", logx=1, logy=1)
    ppgplot.pgmtxt("b", -6.0, 0.95, 1.0, "Sample Rate")
    Pgplot.plotxy(np.log10(chan_smear), ldms, color="purple", logx=1, logy=1)
    ppgplot.pgmtxt("b", -4.5, 0.95, 1.0, "Channel")
    Pgplot.plotxy(np.log10(BW_smear), ldms, color="red", logx=1, logy=1)
    ppgplot.pgmtxt("b", -3.0, 0.95, 1.0, "Full BW")
    Pgplot.plotxy(np.log10(subband_smear), ldms, color="blue", logx=1, logy=1)
    ppgplot.pgmtxt("b", -1.5, 0.95, 1.0, "Subband")
    ppgplot.pgsci(1)


def search_sensitivity(
    Ttot: float,
    G: float,
    BW: float,
    chan: int,
    freq: float,
    T: float,
    dm: float,
    ddm: float,
    dt: float,
    Pmin: float = 0.001,
    Pmax: float = 1.0,
    W: float = 0.1,
    polar: int = 2,
    factor: float = 15.0,
    pts: int = 1000,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the approximate limiting flux density for a pulsar survey.

    This is a *very* approximate calculation. For a better calculation,
    see Cordes and Chernoff, ApJ, 482, p971, App. A.

    Parameters
    ----------
    Ttot : float
        The sky + system temperature in K.
    G : float
        The forward gain of the antenna in K/Jy.
    BW : float
        The observing bandwidth in MHz.
    chan : int
        The number of channels in the filterbank.
    freq : float
        The central observing frequency in MHz.
    T : float
        The integration time in sec.
    dm : float
        The dispersion measure in pc cm^-3.
    ddm : float
        The dispersion measure step size in pc cm^-3.
    dt : float
        The sample time for each data point in sec.
    Pmin : float, optional
        The minimum pulsar period in sec (default 0.001).
    Pmax : float, optional
        The maximum pulsar period in sec (default 1.0).
    W : float, optional
        The duty cycle of the pulsar, 0-1 (default 0.1).
    polar : int, optional
        The number of polarizations (default 2).
    factor : float, optional
        Normalization factor accounting for limiting SNR, hardware
        limitations, etc. (default 15.0).
    pts : int, optional
        The number of points to calculate (default 1000).

    Returns
    -------
    periods : numpy.ndarray
        The pulsar periods in sec.
    S_min : numpy.ndarray
        The approximate limiting flux density in mJy at each period.

    Notes
    -----
    Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy.
    """
    periods = span(Pmin, Pmax, pts)
    widths = (
        np.sqrt(
            (W * periods) ** 2.0
            + dm_smear(dm, BW / chan, freq) ** 2.0
            + dm_smear(ddm / 2.0, BW, freq) ** 2.0
            + dt**2.0
        )
        / periods
    )
    return (
        periods,
        limiting_flux_dens(Ttot, G, BW, T, periods, widths, polar=polar, factor=factor),
    )


def smin_noise(Ttot: float, G: float, BW: float, dt: float) -> float:
    """
    Return the 1-sigma Gaussian noise level per time series bin.

    Default is for a sinusoidal pulse (i.e. W = P / 2) with freq << Nyquist
    frequency.

    Parameters
    ----------
    Ttot : float
        The sky + system temperature in K.
    G : float
        The forward gain of the antenna in K/Jy.
    BW : float
        The observing bandwidth in MHz.
    dt : float
        The time per time series bin in sec.

    Returns
    -------
    float
        The 1-sigma noise level in mJy.

    Notes
    -----
    Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy.
    """
    return Ttot / (G * np.sqrt(2 * BW * dt))


def read_profile(filenm: str, normalize: int = 0) -> np.ndarray:
    """
    Read a simple ASCII pulse profile.

    The file has one bin per line. Comments are allowed if they begin
    with '#'.

    Parameters
    ----------
    filenm : str
        The name of the ASCII profile file.
    normalize : int, optional
        If nonzero, pseudo-normalize the profile so its minimum is 0 and
        its maximum is 1 (default 0).

    Returns
    -------
    numpy.ndarray
        The pulse profile.
    """
    prof = []
    for line in open(filenm):
        if line.startswith("#"):
            continue
        else:
            prof.append(float(line.split()[-1]))
    prof = np.asarray(prof)
    if normalize:
        prof -= min(prof)
        prof /= max(prof)
    return prof


def calc_phs(
    MJD: float | np.ndarray, refMJD: float, *args: float
) -> float | np.ndarray:
    """
    Return the rotational phase (0-1) at a given MJD.

    Parameters
    ----------
    MJD : float or numpy.ndarray
        The MJD (or array of MJDs) at which to evaluate the phase.
    refMJD : float
        The reference MJD.
    *args : float
        The rotational frequency f0 and optional frequency derivatives
        (f1, f2, ...), in order.

    Returns
    -------
    float or numpy.ndarray
        The rotational phase (0-1).
    """
    t = (MJD - refMJD) * pc.SECPERDAY
    n = len(args)  # polynomial order
    nargs = np.concatenate(([0.0], args))
    taylor_coeffs = np.concatenate(
        ([0.0], np.cumprod(1.0 / (np.arange(float(n)) + 1.0)))
    )
    p = np.poly1d((taylor_coeffs * nargs)[::-1])
    return np.fmod(p(t), 1.0)


def calc_freq(
    MJD: float | np.ndarray, refMJD: float, *args: float
) -> float | np.ndarray:
    """
    Return the instantaneous frequency at a given MJD.

    Parameters
    ----------
    MJD : float or numpy.ndarray
        The MJD (or array of MJDs) at which to evaluate the frequency.
    refMJD : float
        The reference MJD.
    *args : float
        The rotational frequency f0 and optional frequency derivatives
        (f1, f2, ...), in order.

    Returns
    -------
    float or numpy.ndarray
        The instantaneous frequency.
    """
    t = (MJD - refMJD) * pc.SECPERDAY
    n = len(args)  # polynomial order
    taylor_coeffs = np.concatenate(
        ([1.0], np.cumprod(1.0 / (np.arange(float(n - 1)) + 1.0)))
    )
    p = np.poly1d((taylor_coeffs * args)[::-1])
    return p(t)


def calc_t0(MJD: float, refMJD: float, *args: float) -> float:
    """
    Return the closest previous MJD corresponding to phase=0.

    Parameters
    ----------
    MJD : float
        The MJD near which to find phase=0.
    refMJD : float
        The reference MJD.
    *args : float
        The spin frequency f0 and optional frequency derivatives
        (f1, f2, ...), in order.

    Returns
    -------
    float
        The closest previous MJD at which the pulse phase is zero.
    """
    phs = calc_phs(MJD, refMJD, *args)
    p = 1.0 / calc_freq(MJD, refMJD, *args)
    return MJD - phs * p / pc.SECPERDAY


def write_princeton_toa(
    toa_MJDi: int,
    toa_MJDf: float,
    toaerr: float,
    freq: float,
    dm: float,
    obs: str = "@",
    name: str = " " * 13,
) -> None:
    """
    Print a TOA in Princeton format.

    Parameters
    ----------
    toa_MJDi : int
        The integer part of the TOA MJD.
    toa_MJDf : float
        The fractional part of the TOA MJD.
    toaerr : float
        The TOA uncertainty in microseconds.
    freq : float
        The observing frequency in MHz.
    dm : float
        The DM correction in pc cm^-3. Only written if nonzero.
    obs : str, optional
        The one-character observatory code, '@' is barycenter (default "@").
    name : str, optional
        A 13-character name field (default 13 spaces).

    Returns
    -------
    None

    Notes
    -----
    Princeton format columns::

        1-1     Observatory (one-character code) '@' is barycenter
        2-2     must be blank
        16-24   Observing frequency (MHz)
        25-44   TOA (decimal point must be in column 30 or column 31)
        45-53   TOA uncertainty (microseconds)
        69-78   DM correction (pc cm^-3)
    """
    # Splice together the fractional and integer MJDs
    toa = "%5d" % int(toa_MJDi) + ("%.13f" % toa_MJDf)[1:]
    if dm != 0.0:
        print(
            obs
            + " %13s %8.3f %s %8.2f              %9.4f" % (name, freq, toa, toaerr, dm)
        )
    else:
        print(obs + " %13s %8.3f %s %8.2f" % (name, freq, toa, toaerr))


def write_tempo2_toa(
    toa_MJDi: int,
    toa_MJDf: float,
    toaerr: float,
    freq: float,
    dm: float,
    obs: str = "@",
    name: str = "unk",
    flags: str = "",
) -> None:
    """
    Print a TOA in Tempo2 format.

    The TOA format is ``file freq sat satErr siteID <flags>``. Note that
    the first line of the file should be "FORMAT 1".

    Parameters
    ----------
    toa_MJDi : int
        The integer part of the TOA MJD.
    toa_MJDf : float
        The fractional part of the TOA MJD.
    toaerr : float
        The TOA uncertainty in microseconds.
    freq : float
        The observing frequency in MHz.
    dm : float
        The DM correction in pc cm^-3. If nonzero, added as a ``-dm`` flag.
    obs : str, optional
        The observatory (site) code (default "@").
    name : str, optional
        The file/name field (default "unk").
    flags : str, optional
        Additional flags to append (default "").

    Returns
    -------
    None
    """
    toa = "%5d" % int(toa_MJDi) + ("%.13f" % toa_MJDf)[1:]
    if dm != 0.0:
        flags += "-dm %.4f" % (dm,)
    print("%s %f %s %.2f %s %s" % (name, freq, toa, toaerr, obs, flags))


def rotate(arr: np.ndarray, bins: int) -> np.ndarray:
    """
    Rotate an array to the left.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to rotate.
    bins : int
        The number of places to rotate to the left.

    Returns
    -------
    numpy.ndarray
        The rotated array.
    """
    bins = int(bins) % len(arr)
    if bins == 0:
        return arr
    else:
        return np.concatenate((arr[bins:], arr[:bins]))


def interp_rotate(arr: np.ndarray, bins: float, zoomfact: int = 10) -> np.ndarray:
    """
    Sinc-interpolate and rotate an array to the left.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to rotate.
    bins : float
        The number of places to rotate to the left. May be fractional; it
        is rounded to the closest whole number of interpolated bins.
    zoomfact : int, optional
        The interpolation zoom factor (default 10).

    Returns
    -------
    numpy.ndarray
        The rotated array, the same length as the original.
    """
    newlen = len(arr) * zoomfact
    rotbins = int(np.floor(bins * zoomfact + 0.5)) % newlen
    newarr = sinc_interp.periodic_interp(arr, zoomfact)
    return rotate(newarr, rotbins)[::zoomfact]


def fft_rotate(arr: np.ndarray, bins: float) -> np.ndarray:
    """
    Rotate an array to the left using the FFT Shift Theorem.

    Parameters
    ----------
    arr : array_like
        The array to rotate.
    bins : float
        The number of places to rotate to the left. May be fractional.

    Returns
    -------
    numpy.ndarray
        The rotated array, the same length as the original.
    """
    arr = np.asarray(arr)
    freqs = np.arange(arr.size / 2 + 1, dtype=float)
    phasor = np.exp(complex(0.0, pc.TWOPI) * freqs * bins / float(arr.size))
    return np.fft.irfft(phasor * np.fft.rfft(arr), arr.size)


def corr(profile: np.ndarray, template: np.ndarray) -> np.ndarray:
    """
    Cross-correlate a profile and a template using FFTs.

    Parameters
    ----------
    profile : numpy.ndarray
        The pulse profile.
    template : numpy.ndarray
        The template to correlate against.

    Returns
    -------
    numpy.ndarray
        The cross-correlation of `profile` and `template`.
    """
    return FFT.irfft(FFT.rfft(template) * np.conjugate(FFT.rfft(profile)), profile.size)


def autocorr(x: np.ndarray) -> np.ndarray:
    """
    Circular normalized auto-correlation of a real function.

    Parameters
    ----------
    x : numpy.ndarray
        The (real) function to auto-correlate.

    Returns
    -------
    numpy.ndarray
        The normalized auto-correlation. Only the first N/2+1 points are
        returned, as the remaining N/2-1 points are symmetric
        (corresponding to negative lags).
    """
    fftx = FFT.rfft(x)
    acf = FFT.irfft(fftx * np.conjugate(fftx), x.size)[: len(x) // 2 + 1]
    return acf / acf[0]


def maxphase(profile: np.ndarray, template: np.ndarray) -> float:
    """
    Return the phase offset that best matches a profile to a template.

    Parameters
    ----------
    profile : numpy.ndarray
        The pulse profile.
    template : numpy.ndarray
        The template.

    Returns
    -------
    float
        The phase offset (0-1) required to best match `profile` to
        `template`.
    """
    return float(np.argmax(corr(profile, template))) / len(profile)


def linear_interpolate(vector: np.ndarray, zoom: int = 10) -> np.ndarray:
    """
    Linearly interpolate a vector.

    Parameters
    ----------
    vector : numpy.ndarray
        The vector to interpolate.
    zoom : int, optional
        The interpolation factor (default 10).

    Returns
    -------
    numpy.ndarray
        The interpolated vector, `zoom` times longer than `vector`.
    """
    n = len(vector)
    ivect = np.zeros(zoom * n, dtype="d")
    nvect = np.concatenate((vector, vector[:1]))
    ivals = np.arange(zoom, dtype="d") / zoom
    loy = nvect[0]
    for ii in range(n):
        hiy = nvect[ii + 1]
        ivect[ii * zoom : (ii + 1) * zoom] = ivals * (hiy - loy) + loy
        loy = hiy
    return ivect


def downsample(vector: np.ndarray, factor: int) -> np.ndarray | int:
    """
    Downsample a vector by co-adding consecutive samples.

    Parameters
    ----------
    vector : numpy.ndarray
        The vector to downsample. Its length must be divisible by `factor`.
    factor : int
        The integer downsample factor.

    Returns
    -------
    numpy.ndarray or int
        The downsampled vector, or 0 if the length of `vector` is not
        divisible by `factor`.
    """
    if len(vector) % factor:
        print("Length of 'vector' is not divisible by 'factor'=%d!" % factor)
        return 0
    newvector = np.reshape(vector, (len(vector) // factor, factor))
    return np.add.reduce(newvector, 1)


def measure_phase_corr(
    profile: np.ndarray, template: np.ndarray, zoom: int = 10
) -> float:
    """
    Return the best-match phase offset between an interpolated profile and template.

    Each of `profile` and `template` is interpolated by a factor of `zoom`
    before the correlation is measured.

    Parameters
    ----------
    profile : numpy.ndarray
        The pulse profile.
    template : numpy.ndarray
        The template.
    zoom : int, optional
        The interpolation factor (default 10).

    Returns
    -------
    float
        The phase offset (0-1) required to best match `profile` to
        `template`.
    """
    zoomprof = zoomtemp = zoom
    if len(template) != len(profile):
        if len(template) % len(profile) == 0:
            zoomprof = zoom * len(template) // len(profile)
        else:
            print(
                "Warning!:  The lengths of the template (%d) and profile (%d)"
                % (len(template), len(profile))
            )
            print("           are not the same!")
    # itemp = linear_interpolate(rotate(template, np.argmax(template)), zoomtemp)
    itemp = linear_interpolate(template, zoomtemp)
    iprof = linear_interpolate(profile, zoomprof)
    return maxphase(iprof, itemp)


def harm_to_sum(fwhm: float) -> int:
    """
    Return the optimal number of harmonics to sum incoherently.

    For an MVMD (Modified Von Mises Distribution) profile.

    Parameters
    ----------
    fwhm : float
        The pulse full width at half-max (0-1).

    Returns
    -------
    int
        The optimal number of harmonics to sum incoherently.
    """
    # fmt: off
    fwhms = [0.0108, 0.0110, 0.0113, 0.0117, 0.0119, 0.0124, 0.0127, 0.0132,
             0.0134, 0.0140, 0.0145, 0.0151, 0.0154, 0.0160, 0.0167, 0.0173,
             0.0180, 0.0191, 0.0199, 0.0207, 0.0220, 0.0228, 0.0242, 0.0257,
             0.0273, 0.0295, 0.0313, 0.0338, 0.0366, 0.0396, 0.0437, 0.0482,
             0.0542, 0.0622, 0.0714, 0.0836, 0.1037, 0.1313, 0.1799, 0.2883]
    # fmt: on
    return len(fwhms) - bisect.bisect(fwhms, fwhm) + 1


def expcos_profile(N: int, phase: float, fwhm: float) -> np.ndarray:
    """
    Return an 'Exponentiated Sinusoid' pulse profile.

    The profile has an integrated 'flux' of 1 unit.

    Parameters
    ----------
    N : int
        The number of points in the profile.
    phase : float
        The pulse phase (0-1).
    fwhm : float
        The pulse full width at half-max (0.0 < fwhm <= 0.5).

    Returns
    -------
    numpy.ndarray
        The pulse profile with `N` bins.
    """
    from presto.simple_roots import secant

    def fwhm_func(k, fwhm=fwhm):
        if fwhm < 0.02:
            return np.arccos(1.0 - np.log(2.0) / k) / pc.PI - fwhm
        else:
            return np.arccos(np.log(0.5 * (np.exp(k) + np.exp(-k))) / k) / pc.PI - fwhm

    phsval = pc.TWOPI * np.arange(N, dtype="d") / float(N)
    phi = -phase * pc.TWOPI
    if fwhm >= 0.5:
        return np.cos(phsval + phi) + 1.0
    elif fwhm < 0.02:
        # The following is from expanding of iO(x) as x->Infinity.
        k = np.log(2.0) / (1.0 - np.cos(pc.PI * fwhm))
        # print("Expansion:  k = %f  FWHM = %f" % (k, fwhm_func(k, 0.0)))
        phsval = np.fmod(phsval + phi, pc.TWOPI)
        phsval = np.where(np.greater(phsval, pc.PI), phsval - pc.TWOPI, phsval)
        denom = (
            1
            + 1 / (8 * k)
            + 9 / (128 * k * k)
            + 75 / (1024 * k**3)
            + 3675 / (32768 * k**4)
            + 59535 / (262144 * k**5)
        ) / np.sqrt(pc.TWOPI * k)
        return np.where(
            np.greater(np.fabs(phsval / pc.TWOPI), 3.0 * fwhm),
            0.0,
            np.exp(k * (np.cos(phsval) - 1.0)) / denom,
        )
    else:
        k = secant(fwhm_func, 1e-8, 0.5)
        norm = 1.0 / (i0(k) - np.exp(-k))
        # print("Full Calc:  k = %f  FWHM = %f" % (k, fwhm_func(k, 0.0)))
    if k < 0.05:
        tmp = np.cos(phsval + phi)
        tmp2 = tmp * tmp
        return norm * (
            k * (tmp + 1)
            + k * k * (tmp2 - 1.0) / 2.0
            + k * k * k * (tmp2 * tmp + 1.0) / 6.0
        )
    else:
        return norm * (np.exp(k * np.cos(phsval + phi)) - np.exp(-k))


def read_gaussfitfile(
    gaussfitfile: str,
    proflen: int,
    rotate: bool | float = True,
    normalize: bool = False,
) -> np.ndarray | float:
    """
    Read a Gaussian-fit file as created by pygaussfit.py.

    Parameters
    ----------
    gaussfitfile : str
        The name of the Gaussian-fit file.
    proflen : int
        The number of bins to include in the resulting template.
    rotate : bool or float, optional
        If True, the resulting profile is rotated so that the Gaussian with
        the largest amplitude is placed at phase=0. If a number, all
        Gaussians are rotated leftward by that amount of pulse phase (0-1).
        If False, no rotation is performed (default True).
    normalize : bool, optional
        If True, the sum of all Gaussians will be 1.0 and the minimum value
        will be 0.0 (default False).

    Returns
    -------
    numpy.ndarray or float
        A `proflen`-length template array, or 0.0 if the numbers of phases,
        amplitudes, and FWHMs in the file are not the same.
    """
    phass = []
    ampls = []
    fwhms = []
    for line in open(gaussfitfile):
        if line.lstrip().startswith("phas"):
            phass.append(float(line.split()[2]))
        if line.lstrip().startswith("ampl"):
            ampls.append(float(line.split()[2]))
        if line.lstrip().startswith("fwhm"):
            fwhms.append(float(line.split()[2]))
    if not (len(phass) == len(ampls) == len(fwhms)):
        print(
            "Number of phases, amplitudes, and FWHMs are not the same in '%s'!"
            % gaussfitfile
        )
        return 0.0
    phass = np.asarray(phass)
    ampls = np.asarray(ampls)
    fwhms = np.asarray(fwhms)
    # Now sort them all according to decreasing amplitude
    new_order = np.argsort(ampls)
    new_order = new_order[::-1]
    ampls = np.take(ampls, new_order)
    phass = np.take(phass, new_order)
    fwhms = np.take(fwhms, new_order)
    if rotate is not False:
        if rotate is True:
            # Put the biggest gaussian at phase = 0.0
            phass = phass - phass[0]
        elif isinstance(rotate, float):
            phass = phass - rotate
        phass = phass % 1.0
    template = np.zeros(proflen, dtype="d")
    for ii in range(len(ampls)):
        template += ampls[ii] * gaussian_profile(proflen, phass[ii], fwhms[ii])
    if normalize:
        template -= template.min()
        template = template / template.sum()
    return template


def gaussian_profile(N: int, phase: float, fwhm: float) -> np.ndarray:
    """
    Return a Gaussian pulse profile.

    The profile has an integrated 'flux' of 1 unit.

    Parameters
    ----------
    N : int
        The number of points in the profile.
    phase : float
        The pulse phase (0-1).
    fwhm : float
        The Gaussian pulse's full width at half-max.

    Returns
    -------
    numpy.ndarray
        The Gaussian pulse profile with `N` bins.

    Notes
    -----
    The FWHM of a Gaussian is approx 2.35482 sigma.
    """
    sigma = fwhm / 2.35482
    mean = phase % 1.0  # Ensures between 0-1
    phss = np.arange(N, dtype=np.float64) / N - mean
    # Following two lines allow the Gaussian to wrap in phase
    phss[phss > 0.5] -= 1.0
    phss[phss < -0.5] += 1.0
    zs = np.fabs(phss) / sigma
    # The following avoids overflow by truncating the Gaussian at 20 sigma
    return np.where(
        zs < 20.0, np.exp(-0.5 * zs**2.0) / (sigma * np.sqrt(2 * np.pi)), 0.0
    )


def gauss_profile_params(profile: np.ndarray, output: int = 0) -> tuple[float, ...]:
    """
    Return the parameters of a best-fit Gaussian to a profile.

    Parameters
    ----------
    profile : array_like
        The pulse profile to fit.
    output : int, optional
        If nonzero, the fit is plotted and the return values are printed
        (default 0).

    Returns
    -------
    tuple of float
        A 6-tuple ``(flux, fwhm, phase, baseline, resid_avg, resid_std)``:

        - flux : best-fit Gaussian integrated 'flux'.
        - fwhm : best-fit Gaussian FWHM.
        - phase : best-fit Gaussian phase (0.0-1.0).
        - baseline : baseline (i.e. noise) average value.
        - resid_avg : residuals average value.
        - resid_std : residuals standard deviation.
    """
    profile = np.asarray(profile)

    def funct(afpo, profile):
        return (
            afpo[0] * gaussian_profile(len(profile), afpo[2], afpo[1])
            + afpo[3]
            - profile
        )

    ret = leastsq(
        funct,
        [
            profile.max() - profile.min(),
            0.25,
            profile.argmax() / float(len(profile)),
            profile.min(),
        ],
        args=(profile),
    )
    if output:
        phases = np.arange(0.0, 1.0, 1.0 / len(profile)) + 0.5 / len(profile)
        Pgplot.plotxy(
            profile,
            phases,
            rangex=[0.0, 1.0],
            labx="Pulse Phase",
            laby="Pulse Intensity",
        )
    bestfit = (
        ret[0][0] * gaussian_profile(len(profile), ret[0][2], ret[0][1]) + ret[0][3]
    )
    if output:
        Pgplot.plotxy(bestfit, phases, color="red")
        Pgplot.closeplot()
    residuals = bestfit - profile
    resid_avg = residuals.mean()
    resid_std = residuals.std()
    if output:
        Pgplot.plotxy(
            residuals,
            phases,
            rangex=[0.0, 1.0],
            rangey=[min(residuals) - 2 * resid_std, max(residuals) + 2 * resid_std],
            labx="Pulse Phase",
            laby="Residuals",
            line=None,
            symbol=3,
        )
        ppgplot.pgerrb(
            6, phases, residuals, np.zeros(len(residuals), "d") + resid_std, 2
        )
        Pgplot.plotxy([resid_avg, resid_avg], [0.0, 1.0], line=2)
        Pgplot.closeplot()
        print("")
        print("  Best-fit gaussian integrated 'flux'  = ", ret[0][0])
        print("               Best-fit gaussian FWHM  = ", ret[0][1])
        print("    Best-fit gaussian phase (0.0-1.0)  = ", ret[0][2])
        print("        Baseline (i.e. noise) average  = ", ret[0][3])
        print("                    Residuals average  = ", resid_avg)
        print("         Residuals standard deviation  = ", resid_std)
        print("")
    return (ret[0][0], ret[0][1], ret[0][2], ret[0][3], resid_avg, resid_std)


def twogauss_profile_params(profile: np.ndarray, output: int = 0) -> tuple[float, ...]:
    """
    Return the parameters of two best-fit Gaussians to a profile.

    Parameters
    ----------
    profile : array_like
        The pulse profile to fit.
    output : int, optional
        If nonzero, the fit is plotted and the return values are printed
        (default 0).

    Returns
    -------
    tuple of float
        A 9-tuple ``(flux1, fwhm1, phase1, flux2, fwhm2, phase2, baseline,
        resid_avg, resid_std)``:

        - flux1, fwhm1, phase1 : first best-fit Gaussian integrated 'flux',
          FWHM, and phase (0.0-1.0).
        - flux2, fwhm2, phase2 : second best-fit Gaussian integrated 'flux',
          FWHM, and phase (0.0-1.0).
        - baseline : baseline (i.e. noise) average value.
        - resid_avg : residuals average value.
        - resid_std : residuals standard deviation.
    """

    def yfunct(afpo, n):
        return (
            afpo[0] * gaussian_profile(n, afpo[2], afpo[1])
            + afpo[3] * gaussian_profile(n, afpo[5], afpo[4])
            + afpo[6]
        )

    def min_funct(afpo, profile):
        return yfunct(afpo, len(profile)) - profile

    ret = leastsq(
        min_funct,
        [
            max(profile) - min(profile),
            0.05,
            np.argmax(profile) / float(len(profile)),
            0.2 * max(profile) - min(profile),
            0.1,
            np.fmod(np.argmax(profile) / float(len(profile)) + 0.5, 1.0),
            min(profile),
        ],
        args=(profile),
    )
    if output:
        phases = np.arange(0.0, 1.0, 1.0 / len(profile)) + 0.5 / len(profile)
        Pgplot.plotxy(
            profile,
            phases,
            rangex=[0.0, 1.0],
            labx="Pulse Phase",
            laby="Pulse Intensity",
        )
    bestfit = yfunct(ret[0], len(profile))
    if output:
        Pgplot.plotxy(bestfit, phases, color="red")
        Pgplot.closeplot()
    residuals = bestfit - profile
    resid_avg = residuals.mean()
    resid_std = residuals.std()
    if output:
        Pgplot.plotxy(
            residuals,
            phases,
            rangex=[0.0, 1.0],
            rangey=[min(residuals) - 2 * resid_std, max(residuals) + 2 * resid_std],
            labx="Pulse Phase",
            laby="Residuals",
            line=None,
            symbol=3,
        )
        ppgplot.pgerrb(
            6, phases, residuals, np.zeros(len(residuals), "d") + resid_std, 2
        )
        Pgplot.plotxy([resid_avg, resid_avg], [0.0, 1.0], line=2)
        Pgplot.closeplot()
        print("")
        print("  Best-fit gaussian integrated 'flux'  = ", ret[0][0])
        print("               Best-fit gaussian FWHM  = ", ret[0][1])
        print("    Best-fit gaussian phase (0.0-1.0)  = ", ret[0][2])
        print("  Best-fit gaussian integrated 'flux'  = ", ret[0][3])
        print("               Best-fit gaussian FWHM  = ", ret[0][4])
        print("    Best-fit gaussian phase (0.0-1.0)  = ", ret[0][5])
        print("        Baseline (i.e. noise) average  = ", ret[0][6])
        print("                    Residuals average  = ", resid_avg)
        print("         Residuals standard deviation  = ", resid_std)
        print("")
    return (
        ret[0][0],
        ret[0][1],
        ret[0][2],
        ret[0][3],
        ret[0][4],
        ret[0][5],
        ret[0][6],
        resid_avg,
        resid_std,
    )


def estimate_flux_density(
    profile: np.ndarray,
    N: int,
    dt: float,
    Ttot: float,
    G: float,
    BW: float,
    prof_stdev: float,
    display: int = 0,
) -> float:
    """
    Return an estimate of the flux density of a pulsar.

    Parameters
    ----------
    profile : numpy.ndarray
        The pulse profile.
    N : int
        The number of time series bins folded.
    dt : float
        The time per time series bin in sec.
    Ttot : float
        The sky + system temperature in K.
    G : float
        The forward gain of the antenna in K/Jy.
    BW : float
        The observing bandwidth in MHz.
    prof_stdev : float
        The profile standard deviation.
    display : int, optional
        If nonzero, the Gaussian fit plots are shown (default 0).

    Returns
    -------
    float
        The estimated flux density in mJy.

    Notes
    -----
    Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy.
    """
    (amp, fwhm, phase, offset, resid_avg, resid_std) = gauss_profile_params(
        profile, display
    )
    T = N * dt
    norm_fact = (prof_stdev * len(profile)) / smin_noise(Ttot, G, BW, T / len(profile))
    return np.add.reduce(profile - offset) / norm_fact


def max_spike_power(FWHM: float) -> float:
    """
    Return the approximate max power ratio of a spike to a sinusoidal profile.

    This is the ratio of the highest power from a triangular spike pulse
    profile to the power from a perfect sinusoidal pulse profile. Both the
    spike and the sine are assumed to have an area under one full pulse of
    1 unit. A Gaussian profile gives almost identical powers as a spike
    profile of the same width. This expression was determined using a
    least-squares fit (max abs error ~ 0.016).

    Parameters
    ----------
    FWHM : float
        The full width at half-max of the spike (0.0 < FWHM <= 0.5).

    Returns
    -------
    float
        The approximate ratio of spike power to sinusoidal power.
    """
    return (
        (36.4165309504 * FWHM - 32.0107844537) * FWHM + 0.239948319674
    ) * FWHM + 4.00277916584


def num_spike_powers(FWHM: float) -> float:
    """
    Return the approximate number of high powers from a spike profile.

    This is the number of powers from a triangular spike pulse profile
    which are greater than one half the power of a perfect sinusoidal pulse
    profile. Both the spike and the sine are assumed to have an area under
    one full pulse of 1 unit. A Gaussian profile gives almost identical
    numbers of high powers as a spike profile of the same width. This
    expression was determined using a least-squares fit (errors get large
    as FWHM -> 0).

    Parameters
    ----------
    FWHM : float
        The full width at half-max of the spike (0.0 < FWHM <= 0.5).

    Returns
    -------
    float
        The approximate number of high powers.
    """
    return -3.95499721563e-05 / FWHM**2 + 0.562069634689 / FWHM - 0.683604041138


def incoherent_sum(amps: np.ndarray) -> np.ndarray:
    """
    Return the accumulated incoherently-summed powers.

    Parameters
    ----------
    amps : numpy.ndarray
        A series of complex Fourier amplitudes.

    Returns
    -------
    numpy.ndarray
        A vector showing the accumulated incoherently-summed powers.
    """
    return np.add.accumulate(np.absolute(amps) ** 2.0)


def coherent_sum(amps: np.ndarray) -> np.ndarray:
    """
    Return the accumulated coherently-summed powers.

    Parameters
    ----------
    amps : numpy.ndarray
        A series of complex Fourier amplitudes.

    Returns
    -------
    numpy.ndarray
        A vector showing the accumulated coherently-summed powers.
    """
    phss = np.arctan2(amps.imag, amps.real)
    phs0 = phss[0]
    phscorr = phs0 - np.fmod((np.arange(len(amps), dtype="d") + 1.0) * phs0, pc.TWOPI)
    sumamps = np.add.accumulate(amps * np.exp(complex(0.0, 1.0) * phscorr))
    return np.absolute(sumamps) ** 2.0


def dft_vector_response(
    roff: float, z: float = 0.0, w: float = 0.0, phs: float = 0.0, N: int = 1000
) -> np.ndarray:
    """
    Return the complex DFT response as a vector sum for a noise-less signal.

    Parameters
    ----------
    roff : float
        The Fourier frequency offset (roff=0 means we are exactly at the
        signal frequency).
    z : float, optional
        The average Fourier f-dot (default 0.0).
    w : float, optional
        The Fourier 2nd derivative (default 0.0).
    phs : float, optional
        An optional phase in radians (default 0.0).
    N : int, optional
        The number of vectors to sum (default 1000).

    Returns
    -------
    numpy.ndarray
        A complex vector addition of `N` vectors showing the DFT response.
    """
    r0 = roff - 0.5 * z + w / 12.0  # Make symmetric for all z and w
    z0 = z - 0.5 * w
    us = np.linspace(0.0, 1.0, N)
    phss = 2.0 * np.pi * (us * (us * (us * w / 6.0 + z0 / 2.0) + r0) + phs)
    return np.cumsum(np.exp(complex(0.0, 1.0) * phss)) / N


def prob_power(power: float | np.ndarray) -> float | np.ndarray:
    """
    Return the probability for noise to exceed a normalized power level.

    Parameters
    ----------
    power : float or numpy.ndarray
        The normalized power level(s).

    Returns
    -------
    float or numpy.ndarray
        The probability for noise to exceed `power` in a power spectrum.
    """
    return np.exp(-power)


def Ftest(chi2_1: float, dof_1: int, chi2_2: float, dof_2: int) -> float:
    """
    Compute an F-test comparing two models.

    Tests whether a model with extra parameters is significant compared to
    a simpler model. The probability is computed exactly like Sherpa's
    F-test routine (in Ciao) and is also described in the Wikipedia article
    on the F-test: http://en.wikipedia.org/wiki/F-test.

    Parameters
    ----------
    chi2_1 : float
        The (non-reduced) chi^2 of the original model.
    dof_1 : int
        The number of degrees of freedom of the original model.
    chi2_2 : float
        The (non-reduced) chi^2 of the new model (with more fit params).
    dof_2 : int
        The number of degrees of freedom of the new model.

    Returns
    -------
    float
        The probability that the improvement in chi2 is due to chance. A
        low probability means the new fit is quantitatively better, while a
        value near 1 means the new model should likely be rejected.
    """
    delta_chi2 = chi2_1 - chi2_2
    delta_dof = dof_1 - dof_2
    new_redchi2 = chi2_2 / dof_2
    F = (delta_chi2 / delta_dof) / new_redchi2
    return fdtrc(delta_dof, dof_2, F)


def equivalent_gaussian_sigma(p: float | np.ndarray) -> float | np.ndarray:
    """
    Return the equivalent Gaussian sigma for a cumulative probability.

    Return x such that Q(x) = p, where Q(x) is the cumulative normal
    distribution. For very small p, an extended-range approximation is
    used.

    Parameters
    ----------
    p : float or numpy.ndarray
        The cumulative Gaussian probability (or array of probabilities).

    Returns
    -------
    float or numpy.ndarray
        The equivalent Gaussian sigma.
    """
    if np.isscalar(p):
        logp = np.log(p)
        return ndtri(1.0 - p) if logp > -30.0 else extended_equiv_gaussian_sigma(logp)
    else:  # logp is an array
        return _vec_equivalent_gaussian_sigma(p)


_vec_equivalent_gaussian_sigma = np.vectorize(
    equivalent_gaussian_sigma, doc="Vectorized `equivalent_gaussian_sigma` over p"
)


def extended_equiv_gaussian_sigma(logp: float) -> float:
    """
    Return the equivalent Gaussian sigma for the log of a cumulative probability.

    Return x such that Q(x) = p, where Q(x) is the cumulative normal
    distribution. This version uses the rational approximation from
    Abramowitz and Stegun, eqn 26.2.23. Using log(P) as input gives a much
    extended range.

    Parameters
    ----------
    logp : float
        The natural log of the cumulative Gaussian probability.

    Returns
    -------
    float
        The equivalent Gaussian sigma.
    """
    t = np.sqrt(-2.0 * logp)
    num = 2.515517 + t * (0.802853 + t * 0.010328)
    denom = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308))
    return t - num / denom


def log_asymtotic_incomplete_gamma(a: float, z: float) -> float:
    """
    Return the log of the incomplete gamma function in its asymptotic limit.

    This is the limit as z -> infinity, from Abramowitz and Stegun eqn
    6.5.32.

    Parameters
    ----------
    a : float
        The first argument of the incomplete gamma function.
    z : float
        The second argument (assumed large).

    Returns
    -------
    float
        The log of the incomplete gamma function.
    """
    x = 1.0
    newxpart = 1.0
    term = 1.0
    ii = 1
    while np.fabs(newxpart) > 1e-15:
        term *= a - ii
        newxpart = term / z**ii
        x += newxpart
        ii += 1
    return (a - 1.0) * np.log(z) - z + np.log(x)


def log_asymtotic_gamma(z: float) -> float:
    """
    Return the log of the gamma function in its asymptotic limit.

    This is the limit as z -> infinity, from Abramowitz and Stegun eqn
    6.1.41.

    Parameters
    ----------
    z : float
        The argument (assumed large).

    Returns
    -------
    float
        The log of the gamma function.
    """
    x = (z - 0.5) * np.log(z) - z + 0.91893853320467267
    y = 1.0 / (z * z)
    x += (
        (
            (-5.9523809523809529e-4 * y + 7.9365079365079365079365e-4) * y
            - 2.7777777777777777777778e-3
        )
        * y
        + 8.3333333333333333333333e-2
    ) / z
    return x


def prob_sum_powers(power: float, nsum: int) -> float:
    """
    Return the probability for noise to exceed a sum of normalized powers.

    Parameters
    ----------
    power : float
        The summed normalized power level.
    nsum : int
        The number of normalized powers that were summed.

    Returns
    -------
    float
        The probability for noise to exceed `power` in the sum of `nsum`
        normalized powers from a power spectrum.
    """
    # Notes:
    # prob_sum_powers(power, nsum)
    # = scipy.special.gammaincc(nsum, power)
    # = statdists.chi_prob(power*2, nsum*2)
    # = scipy.special.chdtrc(nsum*2, power*2)
    # = Q(power*2|nsum*2)  (from A&S 26.4.19)
    # = Gamma(nsum,power)/Gamma(nsum)
    # = [Gamma(nsum) - gamma(nsum,power)]/Gamma(nsum)
    return chdtrc(2 * nsum, 2.0 * power)


def log_prob_sum_powers(power: float | np.ndarray, nsum: int) -> float | np.ndarray:
    """
    Return the log probability for noise to exceed a sum of normalized powers.

    This version allows the use of very large powers by using asymptotic
    expansions from Abramowitz and Stegun Chap 6.

    Parameters
    ----------
    power : float or numpy.ndarray
        The summed normalized power level(s).
    nsum : int
        The number of normalized powers that were summed.

    Returns
    -------
    float or numpy.ndarray
        The log of the probability for noise to exceed `power` in the sum
        of `nsum` normalized powers from a power spectrum.
    """
    # Notes:
    # prob_sum_powers(power, nsum)
    # = scipy.special.gammaincc(nsum, power)
    # = statdists.chi_prob(power*2, nsum*2)
    # = scipy.special.chdtrc(nsum*2, power*2)
    # = Q(power*2|nsum*2)  (from A&S 26.4.19)
    # = Gamma(nsum,power)/Gamma(nsum)
    # = [Gamma(nsum) - gamma(nsum,power)]/Gamma(nsum)
    #
    # For chi^2 dist with dof=2*nsum, mean=dof and var=2*dof
    # And our powers are 1/2 what they should be in chi^2 dist
    # Set our cutoff above ~10 sigma
    thresh = 0.5 * (2 * nsum + 10 * np.sqrt(4 * nsum))  # (mean + 10*std) / 2
    if np.isscalar(power):
        return (
            np.log(prob_sum_powers(power, nsum))
            if power < thresh
            else log_asymtotic_incomplete_gamma(nsum, power) - log_asymtotic_gamma(nsum)
        )
    else:  # power is an array
        return _vec_log_prob_sum_powers(power, nsum)


_vec_log_prob_sum_powers = np.vectorize(
    log_prob_sum_powers, doc="Vectorized `log_prob_sum_powers` over powers"
)


def sigma_power(power: float | np.ndarray) -> float | np.ndarray:
    """
    Return the approximate equivalent Gaussian sigma for a normalized power.

    Parameters
    ----------
    power : float or numpy.ndarray
        The normalized power level(s).

    Returns
    -------
    float or numpy.ndarray
        The approximate equivalent Gaussian sigma for noise to exceed
        `power` in a power spectrum.
    """
    if np.isscalar(power):
        return (
            np.sqrt(2.0 * power - np.log(pc.PI * power))
            if power > 36.0
            else equivalent_gaussian_sigma(prob_power(power))
        )
    else:  # power is an array
        return _vec_sigma_power(power)


_vec_sigma_power = np.vectorize(sigma_power, doc="Vectorized `sigma_power` over powers")


def sigma_sum_powers(power: float | np.ndarray, nsum: int) -> float | np.ndarray:
    """
    Return the approximate equivalent Gaussian sigma for a sum of powers.

    Parameters
    ----------
    power : float or numpy.ndarray
        The summed normalized power level(s).
    nsum : int
        The number of normalized powers that were summed.

    Returns
    -------
    float or numpy.ndarray
        The approximate equivalent Gaussian sigma for noise to exceed a sum
        of `nsum` normalized powers given by `power` in a power spectrum.
    """
    # For chi^2 dist with dof=2*nsum, mean=dof and var=2*dof
    # And our powers are 1/2 what they should be in chi^2 dist
    # Set our cutoff above ~10 sigma
    thresh = 0.5 * (2 * nsum + 10 * np.sqrt(4 * nsum))  # (mean + 10*std) / 2
    if np.isscalar(power):
        return (
            equivalent_gaussian_sigma(prob_sum_powers(power, nsum))
            if power < thresh
            else extended_equiv_gaussian_sigma(log_prob_sum_powers(power, nsum))
        )
    else:  # power is an array
        return _vec_sigma_sum_powers(power, nsum)


_vec_sigma_sum_powers = np.vectorize(
    sigma_sum_powers, doc="Vectorized `sigma_sum_powers` over powers"
)


def power_at_sigma(sigma: float) -> float:
    """
    Return the normalized power level equivalent to a detection significance.

    Parameters
    ----------
    sigma : float
        The detection significance in sigma.

    Returns
    -------
    float
        The approximate normalized power level equivalent to a detection of
        significance `sigma`.
    """
    return sigma**2 / 2.0 + np.log(np.sqrt(pc.PIBYTWO) * sigma)


def powersum_at_sigma(sigma: float, nsum: int) -> float:
    """
    Return the summed power level equivalent to a detection significance.

    Parameters
    ----------
    sigma : float
        The detection significance in sigma.
    nsum : int
        The number of normalized powers that were summed.

    Returns
    -------
    float
        The approximate sum of `nsum` normalized powers equivalent to a
        detection of significance `sigma`.
    """
    return 0.5 * chdtri(2.0 * nsum, 1.0 - ndtr(sigma))


def cand_sigma(N: int, power: float) -> float:
    """
    Return the sigma of a candidate found in a power spectrum.

    Parameters
    ----------
    N : int
        The number of bins searched in the power spectrum.
    power : float
        The normalized power level of the candidate.

    Returns
    -------
    float
        The sigma of the candidate, accounting for the number of bins
        searched.
    """
    return ndtri(1.0 - N * prob_power(power))


def fft_max_pulsed_frac(N: int, numphot: int, sigma: float = 3.0) -> float:
    """
    Return the approximate maximum pulsed fraction not found in an FFT search.

    For a sinusoidal signal that *wasn't* found in an FFT-based search.

    Parameters
    ----------
    N : int
        The number of bins searched in the FFT.
    numphot : int
        The number of photons present in the data.
    sigma : float, optional
        The confidence (in sigma) of the limit (default 3.0).

    Returns
    -------
    float
        The approximate maximum pulsed fraction.
    """
    # The following is the power level required to get a
    # noise spike that would appear somewhere in N bins
    # at the 'sigma' level
    power_required = -np.log((1.0 - ndtr(sigma)) / N)
    return np.sqrt(4.0 * numphot * power_required) / N


def p_to_f(p: float, pd: float, pdd: float | None = None) -> list[float]:
    """
    Convert period and derivatives to frequency and derivatives (or vice versa).

    Parameters
    ----------
    p : float
        The period (or frequency).
    pd : float
        The period derivative (or frequency derivative).
    pdd : float, optional
        The period second derivative (or frequency second derivative). If
        None, only the first two terms are returned.

    Returns
    -------
    list of float
        ``[f, fd]`` if `pdd` is None, otherwise ``[f, fd, fdd]``, the
        equivalent frequency counterparts (the conversion is symmetric, so
        it also converts f to p).
    """
    f = 1.0 / p
    fd = -pd / (p * p)
    if pdd is None:
        return [f, fd]
    else:
        if pdd == 0.0:
            fdd = 0.0
        else:
            fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
        return [f, fd, fdd]


def pferrs(
    porf: float,
    porferr: float,
    pdorfd: float | None = None,
    pdorfderr: float | None = None,
) -> list[float]:
    """
    Convert period/frequency errors and their derivative errors.

    Parameters
    ----------
    porf : float
        The period or frequency.
    porferr : float
        The error on `porf`.
    pdorfd : float, optional
        The pdot or fdot. If None, only `porf` and its error are converted.
    pdorfderr : float, optional
        The error on `pdorfd`.

    Returns
    -------
    list of float
        ``[1/porf, porferr/porf**2]`` if `pdorfd` is None, otherwise
        ``[forp, forperr, fdorpd, fdorpderr]``.
    """
    if pdorfd is None:
        return [1.0 / porf, porferr / porf**2.0]
    else:
        forperr = porferr / porf**2.0
        fdorpderr = np.sqrt(
            (4.0 * pdorfd**2.0 * porferr**2.0) / porf**6.0 + pdorfderr**2.0 / porf**4.0
        )
        [forp, fdorpd] = p_to_f(porf, pdorfd)
        return [forp, forperr, fdorpd, fdorpderr]


def pdot_from_B(p: float, B: float) -> float:
    """
    Return the pdot implied by a magnetic field strength.

    Parameters
    ----------
    p : float
        The spin period (or pdot) in sec.
    B : float
        The magnetic field strength in gauss.

    Returns
    -------
    float
        The pdot (or p) that the pulsar would experience.
    """
    return (B / 3.2e19) ** 2.0 / p


def pdot_from_age(p: float, age: float) -> float:
    """
    Return the pdot implied by a characteristic age.

    Parameters
    ----------
    p : float
        The spin period in sec.
    age : float
        The characteristic age in yrs.

    Returns
    -------
    float
        The pdot that the pulsar would experience.
    """
    return p / (2.0 * age * pc.SECPERJULYR)


def pdot_from_edot(p: float, edot: float, I: float = 1.0e45) -> float:
    """
    Return the pdot implied by a spin-down luminosity.

    Parameters
    ----------
    p : float
        The spin period in sec.
    edot : float
        The spin-down luminosity Edot in ergs/s.
    I : float, optional
        The moment of inertia in g cm^2 (default 1.0e45).

    Returns
    -------
    float
        The pdot that the pulsar would experience.
    """
    return (p**3.0 * edot) / (4.0 * pc.PI * pc.PI * I)


def pulsar_age(f: float, fdot: float, n: float = 3, fo: float = 1e99) -> float:
    """
    Return the age of a pulsar.

    By default the characteristic age is returned (assuming a braking index
    n=3 and an initial spin frequency fo >> f).

    Parameters
    ----------
    f : float
        The spin frequency in Hz.
    fdot : float
        The frequency derivative in Hz/s.
    n : float, optional
        The braking index (default 3).
    fo : float, optional
        The initial spin frequency in Hz (default 1e99).

    Returns
    -------
    float
        The age of the pulsar in years.
    """
    return -f / ((n - 1.0) * fdot) * (1.0 - (f / fo) ** (n - 1.0)) / pc.SECPERJULYR


def pulsar_edot(f: float, fdot: float, I: float = 1.0e45) -> float:
    """
    Return the pulsar spin-down luminosity (Edot).

    Parameters
    ----------
    f : float
        The spin frequency in Hz.
    fdot : float
        The frequency derivative in Hz/s.
    I : float, optional
        The NS moment of inertia in g cm^2 (default 1.0e45).

    Returns
    -------
    float
        The pulsar Edot in erg/s.
    """
    return -4.0 * pc.PI * pc.PI * I * f * fdot


def pulsar_B(f: float, fdot: float) -> float:
    """
    Return the estimated pulsar surface magnetic field strength.

    Parameters
    ----------
    f : float
        The spin frequency in Hz.
    fdot : float
        The frequency derivative in Hz/s.

    Returns
    -------
    float
        The estimated surface magnetic field strength in Gauss.
    """
    return 3.2e19 * np.sqrt(-fdot / f**3.0)


def pulsar_B_lightcyl(f: float, fdot: float) -> float:
    """
    Return the estimated pulsar magnetic field strength at the light cylinder.

    Parameters
    ----------
    f : float
        The spin frequency in Hz.
    fdot : float
        The frequency derivative in Hz/s.

    Returns
    -------
    float
        The estimated magnetic field strength at the light cylinder in
        Gauss.
    """
    p, pd = p_to_f(f, fdot)
    return 2.9e8 * p ** (-5.0 / 2.0) * np.sqrt(pd)


def psr_info(
    porf: float,
    pdorfd: float,
    time: float | None = None,
    input: str | None = None,
    I: float = 1e45,
) -> None:
    """
    Print a list of standard derived pulsar parameters.

    Derived from the period (or frequency) and its first derivative. The
    routine automatically assumes periods if `porf` <= 1.0 and frequencies
    otherwise. This can be overridden with `input`.

    Parameters
    ----------
    porf : float
        The period (in sec) or frequency (in Hz).
    pdorfd : float
        The period derivative or frequency derivative.
    time : float, optional
        The duration of an observation in sec. If given, the Fourier
        frequency 'r' and Fourier fdot 'z' are also printed.
    input : {"p", "f"}, optional
        Force interpretation of `porf` as a period ("p") or frequency ("f").
    I : float, optional
        The NS moment of inertia in g cm^2 (default 1e45).

    Returns
    -------
    None
    """
    if (input is None and porf > 1.0) or (input == "f" or input == "F"):
        pdorfd = -pdorfd / (porf * porf)
        porf = 1.0 / porf
    [f, fd] = p_to_f(porf, pdorfd)
    print("")
    print("             Period = %f s" % porf)
    print("              P-dot = %g s/s" % pdorfd)
    print("          Frequency = %f Hz" % f)
    print("              F-dot = %g Hz/s" % fd)
    if time:
        print("       Fourier Freq = %g bins" % (f * time))
        print("      Fourier F-dot = %g bins" % (fd * time * time))
    print("              E-dot = %g ergs/s" % pulsar_edot(f, fd, I))
    print("    Surface B Field = %g gauss" % pulsar_B(f, fd))
    print(" Characteristic Age = %g years" % pulsar_age(f, fd))
    print("          Assumed I = %g g cm^2" % I)
    print("")


def doppler(freq_observed: float, voverc: float) -> float:
    """
    Return the frequency emitted by a pulsar given the observed frequency.

    Parameters
    ----------
    freq_observed : float
        The observed frequency in MHz.
    voverc : float
        The radial velocity in units of v/c, with respect to the pulsar.

    Returns
    -------
    float
        The frequency emitted by the pulsar in MHz.
    """
    return freq_observed * (1.0 + voverc)
