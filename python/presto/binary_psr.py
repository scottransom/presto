from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import numpy as Num
from presto import parfile, psr_utils
from presto.psr_constants import *

def myasarray(a):
    if type(a) in [type(1.0),type(1),type(1),type(1j)]:
        a = Num.asarray([a])
    if len(a) == 0:
        a = Num.asarray([a])
    return a

def shapR(m2):
    """
    shapR(m2):
        Return the Shapiro 'R' parameter (in sec) with m2 in
            solar units.
    """
    return Tsun * m2

def shapS(m1, m2, x, pb):
    """
    shapS(m1, m2, x, pb):
        Return the Shapiro 'S' parameter with m1 and m2 in
            solar units, x (asini/c) in sec, and pb in days.
            The Shapiro S param is also equal to sin(i).
    """
    return x * (pb*SECPERDAY/TWOPI)**(-2.0/3.0) * \
           Tsun**(-1.0/3.0) * (m1 + m2)**(2.0/3.0) * 1.0/m2

# Note:  S is also equal to sin(i)

class binary_psr(object):
    """
    class binary_psr

        This class reads in a parfile (the only option for instantiation) of
            a binary pulsar.  It allows access the calculation of the mean,
            eccentric, and true anomalies, orbital position, radial velocity,
            and predicted spin period as a function of time.
            
    """

    def __init__(self, parfilenm):
        self.par = parfile.psr_par(parfilenm)
        if not hasattr(self.par, 'BINARY'):
            print("'%s' doesn't contain parameters for a binary pulsar!")
            return None
        self.PBsec = self.par.PB*SECPERDAY
        self.T0 = self.par.T0

    def calc_anoms(self, MJD):
        """
        calc_anoms(MJD):
            Return a tuple of the mean, eccentric, and true anomalies (all
                in radians) at the barycentric epoch MJD(s).
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0)*SECPERDAY
        sec_since_peri = Num.fmod(difft, self.PBsec)
        sec_since_peri[sec_since_peri < 0.0] += self.PBsec
        mean_anom = sec_since_peri/self.PBsec*TWOPI
        ecc_anom = self.eccentric_anomaly(mean_anom)
        true_anom = psr_utils.true_anomaly(ecc_anom, self.par.E)
        return (mean_anom, ecc_anom, true_anom)

    def most_recent_peri(self, MJD):
        """
        most_recent_peri(MJD):
            Return the MJD(s) of the most recent periastrons that occurred
                before the input MJD(s).
        """
        MJD = myasarray(MJD)
        difft = MJD - self.T0
        days_since_peri = Num.fmod(difft, self.par.PB)
        if (days_since_peri < 0.0): days_since_peri += self.par.PB
        return MJD - days_since_peri

    def eccentric_anomaly(self, mean_anomaly):
        """
        eccentric_anomaly(mean_anomaly):
            Return the eccentric anomaly in radians, given a set of mean_anomalies
                in radians.
        """
        ma = Num.fmod(mean_anomaly, TWOPI)
        ma = Num.where(ma < 0.0, ma+TWOPI, ma)
        eccentricity = self.par.E
        ecc_anom_old = ma
        ecc_anom = ma + eccentricity*Num.sin(ecc_anom_old)
        # This is a simple iteration to solve Kepler's Equation
        while (Num.maximum.reduce(Num.fabs(ecc_anom-ecc_anom_old)) > 5e-15):
            ecc_anom_old = ecc_anom[:]
            ecc_anom = ma + eccentricity*Num.sin(ecc_anom_old)
        return ecc_anom

    def calc_omega(self, MJD):
        """
        calc_omega(MJD):
            Return the argument of periastron (omega in radians) at
            time (or times) MJD(s).
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0)*SECPERDAY
        if hasattr(self.par, 'OMDOT'):
            # Note:  This is an array
            return (self.par.OM + difft/SECPERJULYR*self.par.OMDOT)*DEGTORAD
        else:
            return self.par.OM*DEGTORAD

    def radial_velocity(self, MJD):
        """
        radial_velocity(MJD):
            Return the radial velocity of the pulsar (km/s) at the given MJD(s).
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        c1 = TWOPI*self.par.A1/self.PBsec
        c2 = Num.cos(ws)*Num.sqrt(1-self.par.E*self.par.E)
        sws = Num.sin(ws)
        cea = Num.cos(ea)
        return SOL/1000.0*c1*(c2*cea - sws*Num.sin(ea)) / (1.0 - self.par.E*cea)

    def doppler_period(self, MJD):
        """
        doppler_period(MJD):
            Return the observed pulse spin period in sec at the given MJD(s).
        """
        vs = self.radial_velocity(MJD)*1000.0 # m/s
        return self.par.P0*(1.0+vs/SOL)

    def position(self, MJD, inc=60.0, returnz=False):
        """
        position(MJD, inc=60.0, returnz=False):
            Return the 'x' (along the LOS with + being towards us) and 'y' (in the
                plane of the sky with + being away from the line of nodes and -
                being in the direction of the line of nodes) positions of the
                pulsar with respect to the center of mass in units of lt-sec.
                (Note:  This places the observer at (+inf,0.0) and the line of nodes
                extending towards (0.0,-inf) with the pulsar orbiting (0.0,0.0)
                clockwise).  'inc' is the inclination of the orbit in degrees.
                MJD can be an array.  The return value is (xs, ys).  If returnz
                is True, return (xs, ys, zs), where 'z' is the other in-the-sky
                direction.  These coordinates correspond to the I, J, and K vectors
                in Damour & Taylor (1992) in the following way:
                x = -K
                y = -I
                z = -J
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        orb_phs = ta + ws
        sini = Num.sin(inc*DEGTORAD)
        x = self.par.A1/sini # This is a since A1 is asini
        r = x*(1.0-self.par.E*self.par.E)/(1.0+self.par.E*Num.cos(ta))
        if returnz:
            return -r*Num.sin(orb_phs)*sini, -r*Num.cos(orb_phs), \
                -r*Num.sin(orb_phs)*Num.cos(inc*DEGTORAD)
        else:
            return -r*Num.sin(orb_phs)*sini, -r*Num.cos(orb_phs)

    def reflex_motion(self, MJD, inc, Omega, dist):
        """
        reflex_motion(MJD, inc, Omega, dist):
            Return the projected on-sky orbital reflex motion in mas referenced
                to Omega, which is the line-of-nodes, clockwise from East towards
                North.  This is the definition of Omega used by e.g. Damour & 
                Taylor (1992) and Kopeikin (1996), but note that it differs from
                most non-pulsar applications (in which Omega is measured counter-
                clockwise from North to East). The distance to the pulsar is in 
                kpc. The returned values are dRA (corrected by cos(dec)), dDEC.
        """
        xs, ys, zs = self.position(MJD, inc, returnz=True)
        ys = -ys / dist * 2.003988804115705e-03 # in mas, (i.e. DT92 "I")
        zs = -zs / dist * 2.003988804115705e-03 # in mas, (i.e. DT92 "J")
        sino, coso = Num.sin(omega*DEGTORAD), Num.cos(omega*DEGTORAD)
        # Convert from DT92 I, J to I_0, J_0 (= RA, Dec)
        dRA  = (coso * ys - sino * zs) / Num.cos(self.par.DEC_RAD)
        dDEC = (sino * ys + coso * zs)
        return dRA, dDEC

    def demodulate_TOAs(self, MJD):
        """
        demodulate_TOAs(MJD):
            Return arrival times correctly orbitally de-modulated using
                the iterative procedure described in Deeter, Boynton, and Pravdo
                (1981ApJ...247.1003D, thanks, Deepto!).  This corrects for the
                fact that the emitted times are what you want when you only
                have the arrival times.  MJD can be an array.  The returned
                values are in MJD as well.
        """
        ts = MJD[:]  # start of iteration
        dts = Num.ones_like(MJD)
        # This is a simple Newton's Method iteration based on
        # the code orbdelay.c written by Deepto Chakrabarty
        while (Num.maximum.reduce(Num.fabs(dts)) > 1e-10):
            # radial position in lt-days
            xs = -self.position(ts, inc=90.0)[0]/86400.0
            # radial velocity in units of C
            dxs = self.radial_velocity(ts)*1000.0/SOL
            dts = (ts + xs - MJD) / (1.0 + dxs)
            ts = ts - dts
        return ts

    def shapiro_delays(self, R, S, MJD):
        """
        shapiro_delays(R, S, MJD):
            Return the predicted Shapiro delay (in us) for a variety of
                barycentric MJDs, given the R and S parameters.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        canoms = Num.cos(ea)
        sanoms = Num.sin(ea)
        ecc = self.par.E
        cw = Num.cos(ws)
        sw = Num.sin(ws)
        delay = -2.0e6*R*Num.log(1.0 - ecc*canoms -
                                 S*(sw*(canoms-ecc) +
                                    Num.sqrt((1.0 - ecc*ecc)) * cw * sanoms))
        return delay


    def shapiro_measurable(self, R, S, MJD):
        """
        shapiro_measurable(R, S, MJD):
            Return the predicted _measurable_ Shapiro delay (in us) for a
                variety of barycentric MJDs, given the R and S parameters.  
                This is eqn 28 in Freire & Wex 2010 and is only valid in 
                the low eccentricity limit.
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        Phi = ma + ws
        cbar = Num.sqrt(1.0 - S**2.0)
        zeta = S / (1.0 + cbar)
        h3 = R * zeta**3.0
        sPhi = Num.sin(Phi)
        delay = -2.0e6 * h3 * (
            Num.log(1.0 + zeta*zeta - 2.0 * zeta * sPhi) / zeta**3.0 +
            2.0 * sPhi / zeta**2.0 -
            Num.cos(2.0 * Phi) / zeta)
        return delay


if __name__=='__main__':
    import presto.Pgplot as pg
    
    # The following reproduces the RV plot in Hulse & Taylor, 1975
    psrA = binary_psr("B1913+16.par")
    T0 = 42320.933 # From Hulse & Taylor, 1975
    times = psr_utils.span(0.0, psrA.par.PB, 1000) + T0
    rv = psrA.radial_velocity(times)
    pg.plotxy(rv, (times-T0)*24, \
        labx="Hours since Periastron", laby="Radial Velocity (km.s)")
    pg.closeplot()
    
