## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

from types import StringType, FloatType
from psr_utils import ra_to_rad, dec_to_rad, pferrs
#
# The following are the parameters that are accepted in a
# par file when trying to determine a pulsar ephemeris.
#
#  PEPOCH   Epoch of period/frequency parameters and position (MJD)
#  F0       Pulsar rotation frequency (s-2)
#  F        Alternative for F0
#  F1       Pulsar rotation frequency derivative (s^-2)
#  F2       Pulsar rotation frequency second derivative
#  P0       Pulsar period (s).
#  P        Alternative for P0
#  P1       Pulsar period derivative (10^-15).
#  DM       Dispersion measure (pc cm^-3)
#  A1       Projected pulsar semi-major axis of 1st orbit
#  E        Eccentricity of 1st orbit
#  T0       Epoch of periastron passage of 1st orbit (MJD)
#  TASC     Epoch of ascending node passage (MJD)
#  PB       Period of 1st orbit (days)
#  OM       Longitude of periastron passage, 2st orbit (deg)
#  EPS1     First Laplace parameter [eccentricity times sin(omega)]
#  EPS2     Second Laplace parameter [eccentricity times cos(omega)]
#  EPS1DOT  Time derivative of EPS1
#  EPS2DOT  Time derivative of EPS2
#  OMDOT    Rate of periastron advance (deg/yr) 
#  PBDOT    Rate of change of orbital period (10^-12) 
#  XDOT     Rate of change of projected semi-major axis (-12)
#  EDOT     Rate of change of eccentricity (-12)
#
#  The following are _not_ currently implemented:
#  F3, F4, F5,...  Higher order frequency derivative terms
#  OM2DOT   Second time derivative of angle of periastron (rad/s^2)
#  X2DOT    Second time derivative of projected semi-major axis (1/s)
#

float_keys = ["F", "F0", "F1", "F2", "F3", "F4", "F5", "F6",
              "P", "P0", "P1", "P2", "P3", "P4", "P5", "P6",
              "PEPOCH", "DM", "START", "FINISH", "NTOA",
              "TRES", "TZRMJD", "TZRFRQ", "TZRSITE", "NITS",
              "A1", "XDOT", "E", "EDOT", "T0", "PB", "PBDOT", "OM", "OMDOT",
              "EPS1", "EPS2", "EPS1DOT", "EPS2DOT", "TASC", "LAMBDA", "BETA",
              "RA_RAD", "DEC_RAD", "GAMMA", "SINI", "M2", "MTOT", "FB0", "FB1", "FB2"]
str_keys = ["FILE", "PSR", "RAJ", "DECJ", "EPHEM", "CLK", "BINARY"]

class psr_par:
    def __init__(self, parfilenm):
        self.FILE = parfilenm
        pf = open(parfilenm)
        for line in pf.readlines():
            splitline = line.split()
            key = splitline[0]
            if key in str_keys:
                self.__dict__[key] = splitline[1]
            elif key in float_keys:
                try:
                    self.__dict__[key] = float(splitline[1])
                except ValueError:
                    pass
            if len(splitline)==4:
                self.__dict__[key+'_ERR'] = float(splitline[3])
        if self.__dict__.has_key('RAJ'):
            self.__dict__['RA_RAD'] = ra_to_rad(self.RAJ)
        if self.__dict__.has_key('DECJ'):
            self.__dict__['DEC_RAD'] = dec_to_rad(self.DECJ)
        if self.__dict__.has_key('P'):
            self.__dict__['P0'] = self.P
        if self.__dict__.has_key('P0'):
            self.__dict__['F0'] = 1.0/self.P0
        if self.__dict__.has_key('F0'):
            self.__dict__['P0'] = 1.0/self.F0
        if self.__dict__.has_key('FB0'):
            self.__dict__['PB'] = (1.0/self.FB0)/86400.0
        if self.__dict__.has_key('P0_ERR'):
            if self.__dict__.has_key('P1_ERR'):
                (f, self.__dict__['F0_ERR'], self.__dict__['F1'], self.__dict__['F1_ERR']) = \
                    pferrs(self.P0, self.P0_ERR, self.P1, self.P1_ERR)
            else:
                self.__dict__['F0_ERR'] = self.P0_ERR/(self.P0*self.P0)
        if self.__dict__.has_key('F0_ERR'):
	    if self.__dict__.has_key('F1_ERR'):
                (p, self.__dict__['P0_ERR'], self.__dict__['P1'], self.__dict__['P1_ERR']) = \
                    pferrs(self.F0, self.F0_ERR, self.F1, self.F1_ERR)
            else:
                self.__dict__['P0_ERR'] = self.F0_ERR/(self.F0*self.F0)
        pf.close()
    def __str__(self):
        out = ""
        for k, v in self.__dict__.items():
            if k[:2]!="__":
                if type(self.__dict__[k]) is StringType:
                    out += "%10s = '%s'\n" % (k, v)
                else:
                    out += "%10s = %-20.15g\n" % (k, v)
        return out

if __name__ == '__main__':
    a = psr_par("2140-2310A.par")
    print a
