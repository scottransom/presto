## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

from math import *
import struct, os, os.path

## The most recent catalogs are available here:
## 
## http://www.atnf.csiro.au/research/pulsar/catalogue/

## For the radio pulsars, listed parameters in order are: pulsar J2000
## name, J2000 right ascension (h m s) and its error, J2000 declination
## (deg. ' ") and its error, Galactic longitude and latitude (deg.),
## pulse period (s) and its error, pulse period derivative (10^-15) and
## its error, epoch of the period (MJD), dispersion measure (cm^-3 pc)
## and its error, rotation measure (rad m^-2) and its error, pulsar
## distance (kpc), in most cases estimated from the dispersion measure
## using the Taylor & Cordes (1993) model for the Galactic free electron
## density, 400 MHz mean flux density (mJy) and its error, 1400 MHz mean
## flux density (mJy) and its error, a binary coded integer representing
## surveys in which the pulsar has been detected, and codes for surveys
## in which the pulsar has been detected. Errors are in the last quoted
## digit of the value.

## Binary parameters in order are: J2000 pulsar name, binary orbital
## period (d) and its error, projected semi-major axis (s) and its error,
## orbit eccentricity and its error, longitude of periastron (deg.) and
## its error, epoch of periastron (MJD) and its error. For some
## near-circular orbit pulsars the parameters eps1 = ecc * sin(omega),
## eps2 = ecc * cos(omega) and the epoch of ascending node (MJD) are also
## given.

## For the high-energy-only pulsars, the parameters are pulsar J2000
## name, J2000 right ascension (h m s) and its error, J2000 declination
## (deg. ' ") and its error, Galactic longitude and latitude (deg.),
## pulse period (s) and its error, pulse period derivative (10^-15) and
## its error, epoch of the period (MJD), estimated pulsar distance (kpc),
## a binary coded integer representing surveys in which the pulsar has
## been detected, and codes for surveys in which the pulsar has been
## detected.

ARCSECTORAD = float('4.8481368110953599358991410235794797595635330237270e-6')
RADTOARCSEC = float('206264.80624709635515647335733077861319665970087963')
SECTORAD    = float('7.2722052166430399038487115353692196393452995355905e-5')
RADTOSEC    = float('13750.987083139757010431557155385240879777313391975')
RADTODEG    = float('57.295779513082320876798154814105170332405472466564')
DEGTORAD    = float('1.7453292519943295769236907684886127134428718885417e-2')
RADTOHRS    = float('3.8197186342054880584532103209403446888270314977710')
HRSTORAD    = float('2.6179938779914943653855361527329190701643078328126e-1')
PI          = float('3.1415926535897932384626433832795028841971693993751')
TWOPI       = float('6.2831853071795864769252867665590057683943387987502')

def coord_to_string(h_or_d, m, s):
    """
    coord_to_string(h_or_d, m, s):
       Return a formatted string of RA or DEC values as
       'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
    """
    if (s >= 10.0):
        return "%.2d:%.2d:%.4f" % (h_or_d, m, s)
    else:
        return "%.2d:%.2d:0%.4f" % (h_or_d, m, s)

def rad_to_dms(rad):
    """
    rad_to_dms(rad):
       Convert radians to degrees, minutes, and seconds of arc.
    """
    if (rad < 0.0): sign = -1
    else: sign = 1
    arc = RADTODEG * fmod(fabs(rad), PI)
    d = int(arc)
    arc = (arc - d) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (sign * d, m, s)

def rad_to_hms(rad):
    """
    rad_to_hms(rad):
       Convert radians to hours, minutes, and seconds of arc.
    """
    rad = fmod(rad, TWOPI)
    if (rad < 0.0): rad = rad + TWOPI
    arc = RADTOHRS * rad
    h = int(arc)
    arc = (arc - h) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (h, m, s)

def hms_to_rad(hmslist):
    hour = float(hmslist[0])
    min = float(hmslist[1])
    sec = 0.0
    if (len(hmslist)==3):
        sec = float(hmslist[2])
    if (hour < 0.0): sign = -1
    else: sign = 1
    return sign*SECTORAD*(60.0*(60.0*fabs(hour)+fabs(min))+fabs(sec))

def val_and_err(valstr, errstr):
    exponent = 0
    coord = 0
    if (valstr=='*'):
        return (None, None)
    elif valstr.count(":"): # This is a coordinate
        coord = 1
        val = hms_to_rad(valstr.split(':'))
    else:
        val = float(valstr)
        eplace = valstr.find('E')
        if (eplace > 0):
            exponent = int(valstr[eplace+1:])
            valstr = valstr[:eplace]
    if (errstr=='*'):
        err = 0.0
    else:
        err = float(errstr)
        decimal = valstr.find('.')
        if (decimal > 0): # Decimal was found
            err *= 10.0**-(len(valstr)-decimal-1-exponent)
    if coord:
        err *= SECTORAD
    return (val, err)

class psr:
    def __init__(self, line):
        parts = line.split()
        self.jname = parts[0][1:]
        (self.ra, self.raerr) = val_and_err(parts[1], parts[2])
        (self.dec, self.decerr) = val_and_err(parts[3], parts[4])
        self.dec /= 15.0
        self.decerr /= 15.0
	self.l = float(parts[5])
        self.b = float(parts[6])
        (self.p, self.perr) = val_and_err(parts[7], parts[8])
        (self.pd, self.pderr) = val_and_err(parts[9], parts[10])
        if (self.pd):
            self.pd *= 1.0e-15
            self.pderr *= 1.0e-15
        else:
            self.pd = 0.0
            self.pderr = 0.0
        if (parts[11]=='*'):
            self.pepoch = 51000.0 # Just to pick a reasonable value
        else:
            self.pepoch = float(parts[11])
        if (len(parts)==15):  # The high-energy/AXP tables
            if (parts[12]=='*'):
                self.dist = None
            else:
                self.dist = float(parts[12])
            self.dm = 0.0
            self.dmerr = 0.0
            self.s400 = 0.0
            self.s400err = 0.0
            self.s1400 = 0.0
            self.s1400err = 0.0
        else:  # The radio pulsar table
            (self.dm, self.dmerr) = val_and_err(parts[12], parts[13])
            if (self.dm is None):
                self.dm = 0.0
                self.dmerr = 0.0
            (self.rm, self.rmerr) = val_and_err(parts[14], parts[15])
            if (parts[16]=='*'):
                self.dist = None
            else:
                self.dist = float(parts[16])
            (self.s400, self.s400err) = val_and_err(parts[17], parts[18])
            (self.s1400, self.s1400err) = val_and_err(parts[19], parts[20])
        self.bname = ""
        self.alias = ""
        self.binary = 0
    def __cmp__(self, other):
        return cmp(self.jname, other.jname)
    def __str__(self):
        out = ''
        if (self.bname):
            out = out + "\nPulsar  B%s  (J%s)\n" % \
                  (self.bname, self.jname)
        else:
            out = out + "\nPulsar J%s\n" % (self.jname)
        if (self.alias):
            out = out + "                 Alias = %s\n" % self.alias
        (h, m, s) = rad_to_hms(self.ra)
        serr = RADTOSEC * self.raerr
        out = out + "            RA (J2000) = %s +/- %.4fs\n" % \
              (coord_to_string(h, m, s), serr)
        (d, m, s) = rad_to_dms(self.dec)
        serr = RADTOARCSEC * self.decerr
        out = out + "           DEC (J2000) = %s +/- %.4f\"\n" % \
              (coord_to_string(d, m, s), serr)
        out = out + "                (l, b) = (%.2f, %.2f)\n" % \
              (self.l, self.b)
        out = out + "          DM (cm-3 pc) = %.5g +/- %.5g\n" % \
              (self.dm, self.dmerr)
        if (self.s400 is not None):
            out = out + "        S_400MHz (mJy) = %.2g +/- %.2g\n" % \
            (self.s400, self.s400err)
        if (self.s1400 is not None):
            out = out + "       S_1400MHz (mJy) = %.2g +/- %.2g\n" % \
            (self.s1400, self.s1400err)
        if (self.dist is not None):
            out = out + "        Distance (kpc) = %.3g\n" % self.dist
        out = out + "            Period (s) = %.15g +/- %.15g\n" % \
              (self.p, self.perr)
        out = out + "           P-dot (s/s) = %.8g +/- %.8g\n" % \
              (self.pd, self.pderr)
        out = out + "           Epoch (MJD) = %.10g\n" % self.pepoch
        if (self.binary):
            out = out + "          P_binary (s) = %.10g +/- %.10g\n" % \
                  (self.pb*86400.0, self.pberr*86400.0)
            out = out + "          P_binary (d) = %.10g +/- %.10g\n" % \
                  (self.pb, self.pberr)
            out = out + "        a*sin(i)/c (s) = %.8g +/- %.8g\n" % \
                  (self.x, self.xerr)
            out = out + "          Eccentricity = %.8g +/- %.8g\n" % \
                  (self.e, self.eerr)
            if (self.e > 0.0):
                out = out + "    Long of Peri (deg) = %.10g +/- %.10g\n" % \
                      (self.w, self.werr)
                out = out + "    Time of Peri (MJD) = %.12g +/- %.12g\n" % \
                      (self.To, self.Toerr)
            else:
                out = out + "  T of Ascd Node (MJD) = %.12g +/- %.12g\n" % \
                      (self.To, self.Toerr)
        return out
    def add_bin_params(self, parts):
        self.binary = 1
        (self.pb, self.pberr) = val_and_err(parts[1], parts[2])
        (self.x, self.xerr) = val_and_err(parts[3], parts[4])
        if (self.x is None):
            self.x = 0.0
            self.xerr = 0.0
        (self.e, self.eerr) = val_and_err(parts[5], parts[6])
        if (self.e is None):
            self.e = 0.0
            self.eerr = 0.0
        (self.w, self.werr) = val_and_err(parts[7], parts[8])
        if (self.w is None):
            self.w = 0.0
            self.werr = 0.0
        (self.To, self.Toerr) = val_and_err(parts[9], parts[10])
        if (self.To is None):
            (self.To, self.Toerr) = val_and_err(parts[15], parts[16])
            if (self.To is None):
                self.To = 0.0
                self.Toerr = 0.0
        if (parts[11]!='*'):  # Use the 
            (eps1, eps2) = (float(parts[11]), float(parts[13]))
            self.e = sqrt(eps1*eps1 + eps2*eps2)
            self.eerr = 0.0001 # This needs fixing...
            self.w = RADTODEG*atan2(eps1, eps2)
            if (self.w < 0.0):
                self.w += 360.0
            self.werr = 1.0 # This needs fixing...
            (self.To, self.Toerr) = val_and_err(parts[15], parts[16])
    def pack_structs(self):
        out = struct.pack("13s9s10s12d", \
                          self.jname, self.bname, self.alias.lower(),
                          self.ra, self.raerr, self.dec, self.decerr,
                          self.p, self.perr, self.pd, self.pderr,
                          self.dm, self.dmerr, self.pepoch, self.binary)
        if self.binary:
            out = out + struct.pack("10d",
                                    self.pb, self.pberr, self.x, self.xerr,
                                    self.e, self.eerr, self.w, self.werr,
                                    self.To, self.Toerr)
        return out

pulsars = {}
num_binaries = 0

binflag = 0
presto_path = os.getenv("PRESTO")
infile = open(os.path.join(presto_path, "lib", "psr_export.dat"))
for line in infile.readlines()[1:]:
    if (line[0]=='J'):
        if (binflag==0):
            currentpulsar = psr(line)
            pulsars[currentpulsar.jname] = currentpulsar
        else:
            binprops = line.split()
            pulsars[binprops[0][1:]].add_bin_params(binprops)
            num_binaries += 1
    else:
        binflag = 1
infile.close()

if (0):
    binflag = 0
    infile = open("hep_export.dat")
    for line in infile.readlines()[1:]:
        if (line[0]=='J'):
            if (binflag==0):
                currentpulsar = psr(line)
                pulsars[currentpulsar.jname] = currentpulsar
            else:
                binprops = line.split()
                pulsars[binprops[0][1:]].add_bin_params(binprops)
                num_binaries += 1
        else:
            binflag = 1
    infile.close()

    binflag = 0
    infile = open("axp_export.dat")
    for line in infile.readlines()[1:]:
        if (line[0]=='J'):
            if (binflag==0):
                currentpulsar = psr(line)
                pulsars[currentpulsar.jname] = currentpulsar
            else:
                binprops = line.split()
                pulsars[binprops[0][1:]].add_bin_params(binprops)
                num_binaries += 1
        else:
            binflag = 1
    infile.close()

infile = open(os.path.join(presto_path, "lib", "aliases.txt"))
for line in infile.readlines()[1:]:
    if (line[0]=='J'):
        vals = line.split()
        if (len(vals)>=2):
            if (vals[1][0]=='B'):
                pulsars[vals[0][1:]].bname = vals[1][1:]
            else:
                pulsars[vals[0][1:]].alias = vals[1]
        if (len(vals)==3):
            pulsars[vals[0][1:]].alias = vals[2]
infile.close()

psrs = pulsars.values()
psrs.sort()

if __name__ == '__main__' :
    outfilename = "pulsars.cat"
    outfile = open(outfilename, "w")
    print "Writing %d pulsars (%d binaries) to %s" % \
          (len(psrs), num_binaries, outfilename)
    for ii, psr in enumerate(psrs):
        outfile.write(psr.pack_structs())
    outfile.close()

