from types import StringType, FloatType
import math, re
import psr_utils as pu
try:
    from slalib import sla_ecleq, sla_eqecl, sla_eqgal
    slalib = True
except ImportError:
    slalib = False
    
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

float_keys = ["PEPOCH", "POSEPOCH", "DM", "START", "FINISH", "NTOA",
              "TRES", "TZRMJD", "TZRFRQ", "TZRSITE", "NITS",
              "A1", "XDOT", "E", "ECC", "EDOT", "T0", "PB", "PBDOT", "OM", "OMDOT",
              "EPS1", "EPS2", "EPS1DOT", "EPS2DOT", "TASC", "LAMBDA", "BETA",
              "RA_RAD", "DEC_RAD", "GAMMA", "SINI", "M2", "MTOT",
              "ELAT", "ELONG", "PMLAMBDA", "PMBETA", "PX",
              "PMRA", "PMDEC", "PB_2", "A1_2", "E_2", "T0_2", "OM_2"]
floatn_keys = ["F", "P", "FB", "FD", "DMX_", "DMXEP_", "DMXR1_",
               "DMXR2_", "DMXF1_", "DMXF2_"]
str_keys = ["FILE", "PSR", "PSRJ", "RAJ", "DECJ", "EPHEM", "CLK", "BINARY"]

class psr_par:
    def __init__(self, parfilenm):
        self.FILE = parfilenm
        pf = open(parfilenm)
        for line in pf.readlines():
            # Skip comments
            if line[0]=='#':
                continue
            # Convert any 'D-' or 'D+' to 'E-' or 'E+'
            line = line.replace("D-", "E-")
            line = line.replace("D+", "E+")
            splitline = line.split()
            # Skip blank lines
            if len(splitline)==0:
                continue
            key = splitline[0]
            # Regex checks for non-digit chars, followed by digit chars
            m1 = re.search(r'(\D+)(\d+)$', key)
            # This one looks for the DMX[RF][12]_* params
            m2 = re.search(r'(\D+\d+_)(\d+)$', key)
            if key == "JUMP":
                if splitline[3] not in ['0', '1']:
                    setattr(self, key+'_%s'%splitline[2], float(splitline[3]))
                if len(splitline)==5:
                    if splitline[4] not in ['0', '1']:
                        setattr(self, key+'_%s'%splitline[2]+'_ERR',
                                float(splitline[4]))
                elif len(splitline)==6:
                    setattr(self, key+'_%s'%splitline[2]+'_ERR',
                            float(splitline[5]))
            if key in str_keys:
                setattr(self, key, splitline[1])
            elif key in float_keys:
                try:
                    setattr(self, key, float(splitline[1]))
                except ValueError:
                    pass
            elif m1 is not None:
                m = m1
                if m2 is not None:
                    m = m2
                if m.group(1) in floatn_keys:
                    try:
                        setattr(self, key, float(splitline[1]))
                    except ValueError:
                        pass
            if len(splitline)==3:  # Some parfiles don't have flags, but do have errors
                if splitline[2] not in ['0', '1']:
                    setattr(self, key+'_ERR', float(splitline[2]))
            if len(splitline)==4:
                setattr(self, key+'_ERR', float(splitline[3]))
        # Deal with Ecliptic coords
        if (hasattr(self, 'BETA') and hasattr(self, 'LAMBDA')):
            setattr(self, 'ELAT', self.BETA)
            setattr(self, 'ELONG', self.LAMBDA)
        if (slalib and hasattr(self, 'ELAT') and hasattr(self, 'ELONG')):
            if hasattr(self, 'POSEPOCH'):
                epoch = self.POSEPOCH
            else:
                epoch = self.PEPOCH
            ra_rad, dec_rad = sla_ecleq(self.ELONG*pu.DEGTORAD,
                                        self.ELAT*pu.DEGTORAD, epoch)
            rstr = pu.coord_to_string(*pu.rad_to_hms(ra_rad))
            dstr = pu.coord_to_string(*pu.rad_to_dms(dec_rad))
            setattr(self, 'RAJ', rstr)
            setattr(self, 'DECJ', dstr)
        if hasattr(self, 'RAJ'):
            setattr(self, 'RA_RAD', pu.ra_to_rad(self.RAJ))
        if hasattr(self, 'DECJ'):
            setattr(self, 'DEC_RAD', pu.dec_to_rad(self.DECJ))
        # Compute the Galactic coords
        if (slalib and hasattr(self, 'RA_RAD') and hasattr(self, 'DEC_RAD')):
            l, b = sla_eqgal(self.RA_RAD, self.DEC_RAD)
            setattr(self, 'GLONG', l*pu.RADTODEG)
            setattr(self, 'GLAT', b*pu.RADTODEG)
        # Compute the Ecliptic coords
        if (slalib and hasattr(self, 'RA_RAD') and hasattr(self, 'DEC_RAD')):
            if hasattr(self, 'POSEPOCH'):
                epoch = self.POSEPOCH
            else:
                epoch = self.PEPOCH
            elon, elat = sla_eqecl(self.RA_RAD, self.DEC_RAD, epoch)
            setattr(self, 'ELONG', elon*pu.RADTODEG)
            setattr(self, 'ELAT', elat*pu.RADTODEG)
        if hasattr(self, 'P'):
            setattr(self, 'P0', self.P)
        if hasattr(self, 'P0'):
            setattr(self, 'F0', 1.0/self.P0)
        if hasattr(self, 'F0'):
            setattr(self, 'P0', 1.0/self.F0)
        if hasattr(self, 'FB0'):
            setattr(self, 'PB', (1.0/self.FB0)/86400.0)
        if hasattr(self, 'P0_ERR'):
            if hasattr(self, 'P1_ERR'):
                f, ferr, fd, fderr = pu.pferrs(self.P0, self.P0_ERR,
                                               self.P1, self.P1_ERR)
                setattr(self, 'F0_ERR', ferr) 
                setattr(self, 'F1', fd) 
                setattr(self, 'F1_ERR', fderr) 
            else:
                f, fd, = pu.p_to_f(self.P0, self.P1)
                setattr(self, 'F0_ERR', self.P0_ERR/(self.P0*self.P0))
                setattr(self, 'F1', fd)
        else:
            if hasattr(self, 'P1'):
                f, fd, = pu.p_to_f(self.P0, self.P1)
                setattr(self, 'F1', fd)
            elif hasattr(self, 'F1'):
                p, pd, = pu.p_to_f(self.F0, self.F1)
                setattr(self, 'P1', pd)
        if (hasattr(self, 'F0_ERR') and hasattr(self, 'F1_ERR')):
            p, perr, pd, pderr = pu.pferrs(self.F0, self.F0_ERR, 
                                           self.F1, self.F1_ERR)
            setattr(self, 'P0_ERR', perr) 
            setattr(self, 'P1', pd) 
            setattr(self, 'P1_ERR', pderr) 
        elif (hasattr(self, 'F0') and hasattr(self, 'F0_ERR')):
            setattr(self, 'P0_ERR', self.F0_ERR/(self.F0*self.F0))
        if hasattr(self, 'EPS1') and hasattr(self, 'EPS2'):
            ecc = math.sqrt(self.EPS1 * self.EPS1 + self.EPS2 * self.EPS2)
            omega = math.atan2(self.EPS1, self.EPS2)
            setattr(self, 'E', ecc)
            setattr(self, 'OM', omega * pu.RADTODEG)
            setattr(self, 'T0', self.TASC + self.PB * omega/pu.TWOPI)
        if hasattr(self, 'PB') and hasattr(self, 'A1') and not \
               (hasattr(self, 'E') or hasattr(self, 'ECC')):
            setattr(self, 'E', 0.0)
        if hasattr(self, 'T0') and not hasattr(self, 'TASC'):
            setattr(self, 'TASC', self.T0 - self.PB * self.OM/360.0)
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
