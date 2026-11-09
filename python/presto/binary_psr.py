from __future__ import annotations

import numpy as np

from presto import parfile, psr_utils
from presto.psr_constants import (
    DEGTORAD,
    SECPERDAY,
    SECPERJULYR,
    SOL,
    TWOPI,
    Tsun,
)


def myasarray(a) -> np.ndarray:
    """
    Return `a` as a 1-D array, wrapping a scalar in a length-1 array.

    Parameters
    ----------
    a : scalar or array_like
        The value(s) to convert.

    Returns
    -------
    numpy.ndarray
        `a` as an array (a length-1 array if `a` was a scalar or empty).
    """
    a = np.asarray(a)
    if a.ndim == 0:  # Python or NumPy scalar
        a = a.reshape(1)
    if len(a) == 0:
        a = np.asarray([a])
    return a


def shapR(m2: float) -> float:
    """
    Return the Shapiro 'R' parameter.

    Parameters
    ----------
    m2 : float
        The companion mass in solar units.

    Returns
    -------
    float
        The Shapiro 'R' parameter in sec.
    """
    return Tsun * m2


def shapS(m1: float, m2: float, x: float, pb: float) -> float:
    """
    Return the Shapiro 'S' parameter (equal to sin(i)).

    Parameters
    ----------
    m1 : float
        The pulsar mass in solar units.
    m2 : float
        The companion mass in solar units.
    x : float
        The projected semi-major axis (asini/c) in sec.
    pb : float
        The orbital period in days.

    Returns
    -------
    float
        The Shapiro 'S' parameter, which is also equal to sin(i).
    """
    return (
        x
        * (pb * SECPERDAY / TWOPI) ** (-2.0 / 3.0)
        * Tsun ** (-1.0 / 3.0)
        * (m1 + m2) ** (2.0 / 3.0)
        * 1.0
        / m2
    )


# Note:  S is also equal to sin(i)


class binary_psr(object):
    """
    A binary pulsar, built from a parfile.

    Reads in a parfile (the only option for instantiation) of a binary
    pulsar and provides the calculation of the mean, eccentric, and true
    anomalies, orbital position, radial velocity, and predicted spin period
    as a function of time.

    Parameters
    ----------
    parfilenm : str
        The path to the pulsar's ``.par`` file. It must contain binary
        parameters (a ``BINARY`` line).
    """

    def __init__(self, parfilenm: str):
        self.par = parfile.psr_par(parfilenm)
        if not hasattr(self.par, "BINARY"):
            print(f"'{parfilenm}' doesn't contain parameters for a binary pulsar!")
            return None
        self.PBsec = self.par.PB * SECPERDAY
        self.T0 = self.par.T0

    def calc_anoms(
        self, MJD: float | np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return the mean, eccentric, and true anomalies at the given MJD(s).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The barycentric epoch(s).

        Returns
        -------
        mean_anom : numpy.ndarray
            The mean anomaly in radians.
        ecc_anom : numpy.ndarray
            The eccentric anomaly in radians.
        true_anom : numpy.ndarray
            The true anomaly in radians.
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0) * SECPERDAY
        sec_since_peri = np.fmod(difft, self.PBsec)
        sec_since_peri[sec_since_peri < 0.0] += self.PBsec
        mean_anom = sec_since_peri / self.PBsec * TWOPI
        ecc_anom = self.eccentric_anomaly(mean_anom)
        true_anom = psr_utils.true_anomaly(ecc_anom, self.par.E)
        return (mean_anom, ecc_anom, true_anom)

    def most_recent_peri(self, MJD: float | np.ndarray) -> np.ndarray:
        """
        Return the MJD(s) of the most recent periastron before the input MJD(s).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s) before which to find the most recent periastron.

        Returns
        -------
        numpy.ndarray
            The MJD(s) of the most recent periastron passages.
        """
        MJD = myasarray(MJD)
        difft = MJD - self.T0
        days_since_peri = np.fmod(difft, self.par.PB)
        days_since_peri[days_since_peri < 0.0] += self.par.PB
        return MJD - days_since_peri

    def eccentric_anomaly(self, mean_anomaly: np.ndarray) -> np.ndarray:
        """
        Return the eccentric anomaly for a set of mean anomalies.

        Solves Kepler's equation by simple iteration.

        Parameters
        ----------
        mean_anomaly : numpy.ndarray
            The mean anomalies in radians.

        Returns
        -------
        numpy.ndarray
            The eccentric anomalies in radians.
        """
        ma = np.fmod(mean_anomaly, TWOPI)
        ma = np.where(ma < 0.0, ma + TWOPI, ma)
        eccentricity = self.par.E
        ecc_anom_old = ma
        ecc_anom = ma + eccentricity * np.sin(ecc_anom_old)
        # This is a simple iteration to solve Kepler's Equation
        while np.maximum.reduce(np.fabs(ecc_anom - ecc_anom_old)) > 5e-15:
            ecc_anom_old = ecc_anom[:]
            ecc_anom = ma + eccentricity * np.sin(ecc_anom_old)
        return ecc_anom

    def calc_omega(self, MJD: float | np.ndarray) -> float | np.ndarray:
        """
        Return the argument of periastron at the given MJD(s).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s).

        Returns
        -------
        float or numpy.ndarray
            The argument of periastron (omega) in radians. If ``OMDOT`` is
            present in the parfile, this is an array; otherwise it is the
            constant ``OM``.
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0) * SECPERDAY
        if hasattr(self.par, "OMDOT"):
            # Note:  This is an array
            return (self.par.OM + difft / SECPERJULYR * self.par.OMDOT) * DEGTORAD
        else:
            return self.par.OM * DEGTORAD

    def radial_velocity(self, MJD: float | np.ndarray) -> np.ndarray:
        """
        Return the radial velocity of the pulsar at the given MJD(s).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s).

        Returns
        -------
        numpy.ndarray
            The radial velocity in km/s.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        c1 = TWOPI * self.par.A1 / self.PBsec
        c2 = np.cos(ws) * np.sqrt(1 - self.par.E * self.par.E)
        sws = np.sin(ws)
        cea = np.cos(ea)
        return (
            SOL / 1000.0 * c1 * (c2 * cea - sws * np.sin(ea)) / (1.0 - self.par.E * cea)
        )

    def doppler_period(self, MJD: float | np.ndarray) -> np.ndarray:
        """
        Return the observed (Doppler-shifted) spin period at the given MJD(s).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s).

        Returns
        -------
        numpy.ndarray
            The observed pulse spin period in sec.
        """
        vs = self.radial_velocity(MJD) * 1000.0  # m/s
        return self.par.P0 * (1.0 + vs / SOL)

    def position(
        self, MJD: float | np.ndarray, inc: float = 60.0, returnz: bool = False
    ) -> tuple[np.ndarray, ...]:
        """
        Return the pulsar's position relative to the center of mass.

        The 'x' coordinate is along the line of sight (+ towards us) and the
        'y' coordinate is in the plane of the sky (+ away from the line of
        nodes, - towards it), both in lt-sec. This places the observer at
        ``(+inf, 0.0)`` with the line of nodes extending towards
        ``(0.0, -inf)`` and the pulsar orbiting ``(0.0, 0.0)`` clockwise.

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s).
        inc : float, optional
            The orbital inclination in degrees (default 60.0).
        returnz : bool, optional
            If True, also return the 'z' (other in-the-sky) coordinate
            (default False).

        Returns
        -------
        tuple of numpy.ndarray
            ``(xs, ys)``, or ``(xs, ys, zs)`` if `returnz` is True. These
            correspond to the I, J, K vectors of Damour & Taylor (1992) via
            ``x = -K``, ``y = -I``, ``z = -J``.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        orb_phs = ta + ws
        sini = np.sin(inc * DEGTORAD)
        x = self.par.A1 / sini  # This is a since A1 is asini
        r = x * (1.0 - self.par.E * self.par.E) / (1.0 + self.par.E * np.cos(ta))
        if returnz:
            return (
                -r * np.sin(orb_phs) * sini,
                -r * np.cos(orb_phs),
                -r * np.sin(orb_phs) * np.cos(inc * DEGTORAD),
            )
        else:
            return -r * np.sin(orb_phs) * sini, -r * np.cos(orb_phs)

    def reflex_motion(
        self, MJD: float | np.ndarray, inc: float, Omega: float, dist: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the projected on-sky orbital reflex motion.

        The motion is referenced to `Omega`, the line of nodes, measured
        clockwise from East towards North. This is the definition of Omega
        used by e.g. Damour & Taylor (1992) and Kopeikin (1996), but note
        that it differs from most non-pulsar applications (in which Omega is
        measured counter-clockwise from North to East).

        Parameters
        ----------
        MJD : float or numpy.ndarray
            The epoch(s).
        inc : float
            The orbital inclination in degrees.
        Omega : float
            The longitude of the ascending node in degrees.
        dist : float
            The distance to the pulsar in kpc.

        Returns
        -------
        dRA : numpy.ndarray
            The reflex motion in right ascension (corrected by cos(dec)), in
            mas.
        dDEC : numpy.ndarray
            The reflex motion in declination, in mas.
        """
        xs, ys, zs = self.position(MJD, inc, returnz=True)
        ys = -ys / dist * 2.003988804115705e-03  # in mas, (i.e. DT92 "I")
        zs = -zs / dist * 2.003988804115705e-03  # in mas, (i.e. DT92 "J")
        sino, coso = np.sin(Omega * DEGTORAD), np.cos(Omega * DEGTORAD)
        # Convert from DT92 I, J to I_0, J_0 (= RA, Dec)
        dRA = (coso * ys - sino * zs) / np.cos(self.par.DEC_RAD)
        dDEC = sino * ys + coso * zs
        return dRA, dDEC

    def demodulate_TOAs(self, MJD: np.ndarray) -> np.ndarray:
        """
        Return orbitally de-modulated arrival times.

        Uses the iterative procedure described in Deeter, Boynton, and
        Pravdo (1981ApJ...247.1003D). This corrects for the fact that the
        emitted times are what you want when you only have the arrival
        times.

        Parameters
        ----------
        MJD : numpy.ndarray
            The arrival times (MJDs).

        Returns
        -------
        numpy.ndarray
            The de-modulated (emitted) times, in MJD.
        """
        ts = MJD[:]  # start of iteration
        dts = np.ones_like(MJD)
        # This is a simple Newton's Method iteration based on
        # the code orbdelay.c written by Deepto Chakrabarty
        while np.maximum.reduce(np.fabs(dts)) > 1e-10:
            # radial position in lt-days
            xs = -self.position(ts, inc=90.0)[0] / 86400.0
            # radial velocity in units of C
            dxs = self.radial_velocity(ts) * 1000.0 / SOL
            dts = (ts + xs - MJD) / (1.0 + dxs)
            ts = ts - dts
        return ts

    def shapiro_delays(self, R: float, S: float, MJD: float | np.ndarray) -> np.ndarray:
        """
        Return the predicted Shapiro delay at the given MJD(s).

        Parameters
        ----------
        R : float
            The Shapiro 'R' parameter (see :func:`shapR`).
        S : float
            The Shapiro 'S' parameter (see :func:`shapS`); equal to sin(i).
        MJD : float or numpy.ndarray
            The barycentric epoch(s).

        Returns
        -------
        numpy.ndarray
            The predicted Shapiro delay in microseconds.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        canoms = np.cos(ea)
        sanoms = np.sin(ea)
        ecc = self.par.E
        cw = np.cos(ws)
        sw = np.sin(ws)
        delay = (
            -2.0e6
            * R
            * np.log(
                1.0
                - ecc * canoms
                - S * (sw * (canoms - ecc) + np.sqrt((1.0 - ecc * ecc)) * cw * sanoms)
            )
        )
        return delay

    def shapiro_measurable(
        self, R: float, S: float, MJD: float | np.ndarray
    ) -> np.ndarray:
        """
        Return the predicted *measurable* Shapiro delay at the given MJD(s).

        This is eqn 28 in Freire & Wex 2010 and is only valid in the
        low-eccentricity limit.

        Parameters
        ----------
        R : float
            The Shapiro 'R' parameter (see :func:`shapR`).
        S : float
            The Shapiro 'S' parameter (see :func:`shapS`); equal to sin(i).
        MJD : float or numpy.ndarray
            The barycentric epoch(s).

        Returns
        -------
        numpy.ndarray
            The predicted measurable Shapiro delay in microseconds.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        Phi = ma + ws
        cbar = np.sqrt(1.0 - S**2.0)
        zeta = S / (1.0 + cbar)
        h3 = R * zeta**3.0
        sPhi = np.sin(Phi)
        delay = (
            -2.0e6
            * h3
            * (
                np.log(1.0 + zeta * zeta - 2.0 * zeta * sPhi) / zeta**3.0
                + 2.0 * sPhi / zeta**2.0
                - np.cos(2.0 * Phi) / zeta
            )
        )
        return delay


if __name__ == "__main__":
    import presto.Pgplot as pg

    # The following reproduces the RV plot in Hulse & Taylor, 1975
    psrA = binary_psr("B1913+16.par")
    T0 = 42320.933  # From Hulse & Taylor, 1975
    times = psr_utils.span(0.0, psrA.par.PB, 1000) + T0
    rv = psrA.radial_velocity(times)
    pg.plotxy(
        rv,
        (times - T0) * 24,
        labx="Hours since Periastron",
        laby="Radial Velocity (km.s)",
    )
    pg.closeplot()
