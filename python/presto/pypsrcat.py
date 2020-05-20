from __future__ import print_function
from __future__ import absolute_import
from builtins import object
from operator import attrgetter
import struct
import os.path
import math
import csv
import astropy.coordinates as c
import astropy.units as u
from presto import presto
import presto.psr_utils as pu
import presto.psr_constants as pc

## The most recent catalogs are available here:
## 
## http://www.atnf.csiro.au/research/pulsar/psrcat/

version = 'v1.63'

## And here is the command used to get the data:
# Note version number now!
# http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=1.63&Name=Name&JName=JName&RaJ=RaJ&DecJ=DecJ&PMRA=PMRA&PMDec=PMDec&PX=PX&PosEpoch=PosEpoch&GL=GL&GB=GB&F0=F0&F1=F1&F2=F2&F3=F3&PEpoch=PEpoch&DM=DM&DM1=DM1&S400=S400&S1400=S1400&Binary=Binary&T0=T0&PB=PB&A1=A1&OM=OM&Ecc=Ecc&Tasc=Tasc&Eps1=Eps1&Eps2=Eps2&Dist=Dist&Assoc=Assoc&Survey=Survey&Type=Type&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+csv+with+errors&no_value=*&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=40&table_bottom.y=0

params = ["NAME", "PSRJ", "RAJ", "DECJ", "PMRA", "PMDEC", "PX", "POSEPOCH",
          "Gl", "Gb", "F0", "F1", "F2", "F3", "PEPOCH", "DM", "DM1",
          "S400", "S1400", "BINARY", "T0", "PB", "A1", "OM", "ECC",
          "TASC", "EPS1", "EPS2", "DIST", "ASSOC", "SURVEY", "PSR"]
params_with_errs = ["RAJ", "DECJ", "PMRA", "PMDEC", "PX", "P0", "P1", "F2", "F3",
                    "DM", "DM1", "S400", "S1400", "T0", "PB", "A1", "OM", "ECC",
                    "TASC", "EPS1", "EPS2"]
digits = '0123456789'

class psr(object):
    def __init__(self, parts, indices):
        # Do RAJ and DECJ first
        posn = c.SkyCoord(parts[indices['RAJ']]+" "+parts[indices['DECJ']],
                          frame=c.ICRS, unit=(u.hourangle, u.deg))
        for param in params:
            part_index = indices[param]
            if param=="NAME":
                if not parts[part_index]=='*':
                    self.name = parts[part_index][1:]
                else:
                    self.name = ""
            elif param=="PSRJ":
                if not parts[part_index]=='*':
                    self.jname = parts[part_index][1:]
                    if self.name == self.jname:
                        self.name = ""
            elif param=="RAJ":
                if not parts[part_index]=='*':
                    self.rajstr = parts[part_index]
                    self.ra = posn.ra.to(u.rad).value
                    self.raerr = float(parts[part_index+1]) * pc.SECTORAD
            elif param=="DECJ":
                if not parts[part_index]=='*':
                    self.decjstr = parts[part_index]
                    self.dec = posn.dec.to(u.rad).value
                    self.decerr = float(parts[part_index+1]) * pc.ARCSECTORAD
            elif param=="PMRA":
                if not parts[part_index]=='*':
                    self.pmra, self.pmraerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="PMDEC":
                if not parts[part_index]=='*':
                    self.pmdec, self.pmdecerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="PX":
                if not parts[part_index]=='*':
                    self.px, self.pxerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="POSEPOCH":
                if not parts[part_index]=='*':
                    self.posepoch = float(parts[part_index])
            elif param=="Gl":
                if not parts[part_index]=='*':
                    self.l = float(parts[part_index])
            elif param=="Gb":
                if not parts[part_index]=='*':
                    self.b = float(parts[part_index])
            elif param=="F0":
                if not parts[part_index]=='*':
                    self.f, self.ferr = float(parts[part_index]), float(parts[part_index+1])
                    self.p, self.perr = pu.pferrs(self.f, self.ferr)
                else:
                    self.f = self.ferr = self.p = self.perr = 0.0
                self.fd = self.fdd = self.fddd = 0.0
                self.pd = self.pdd = self.pddd = 0.0
                self.fderr = self.fdderr = self.fddderr = 0.0
                self.pderr = self.pdderr = self.pddderr = 0.0
            elif param=="F1":
                if not parts[part_index]=='*':
                    self.fd, self.fderr = float(parts[part_index]), float(parts[part_index+1])
                    self.p, self.perr, self.pd, self.pderr = pu.pferrs(self.f, self.ferr, self.fd, self.fderr)
            elif param=="F2":
                if not parts[part_index]=='*':
                    self.fdd, self.fdderr = float(parts[part_index]), float(parts[part_index+1])
                    self.p, self.pd, self.pdd = presto.p_to_f(self.f, self.fd, self.fdd)
            elif param=="F3":
                if not parts[part_index]=='*':
                    self.fddd, self.fddderr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="PEPOCH":
                if parts[part_index]=='*': 
                    self.pepoch = 51000.0 # Just to pick a reasonable value
                else: 
                    self.pepoch = float(parts[part_index])
            elif param=="DM":
                if not parts[part_index]=='*':
                    self.dm, self.dmerr = float(parts[part_index]), float(parts[part_index+1])
                else:
                    self.dm = self.dmerr = 0.0
            elif param=="DM1":
                if not parts[part_index]=='*':
                    self.ddm, self.ddmerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="S400":
                if not parts[part_index]=='*':
                    self.s400, self.s400err = float(parts[part_index]), float(parts[part_index+1])
                else:
                    self.s400 = None
            elif param=="S1400":
                if not parts[part_index]=='*':
                    self.s1400, self.s1400err = float(parts[part_index]), float(parts[part_index+1])
                else:
                    self.s1400 = None
            elif param=="BINARY":
                if not parts[part_index]=='*':
                    self.binary_model = parts[part_index]
                    self.binary = 1
                    self.pb = self.x = self.w = self.To = self.e = None
                    self.pberr = self.xerr = self.werr = self.Toerr =self.eerr = None
                else:
                    self.binary = 0
            elif param=="T0":
                if self.binary and not parts[part_index]=='*':
                    self.To, self.Toerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="PB":
                if self.binary and not parts[part_index]=='*':
                    self.pb, self.pberr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="A1":
                if self.binary and not parts[part_index]=='*':
                    self.x, self.xerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="OM":
                if self.binary and not parts[part_index]=='*':
                    self.w, self.werr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="ECC":
                if self.binary and not parts[part_index]=='*':
                    self.e, self.eerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="TASC":
                if self.binary and self.binary_model=="ELL1" and not parts[part_index]=='*':
                    self.To, self.Toerr = float(parts[part_index]), float(parts[part_index+1])
            elif param=="EPS1":
                if self.binary and self.binary_model=="ELL1" and not parts[part_index]=='*':
                    self.eps1, self.eps1err = float(parts[part_index]), float(parts[part_index+1])
            elif param=="EPS2":
                if self.binary and self.binary_model=="ELL1" and not parts[part_index]=='*':
                    self.eps2, self.eps2err = float(parts[part_index]), float(parts[part_index+1])
                    if not hasattr(self, 'eps1'): self.eps1 = 0.0
                    self.e = math.sqrt(self.eps1*self.eps1 + self.eps2*self.eps2)
                    self.eerr = 0.0001 # This needs fixing...
                    self.w = pc.RADTODEG*math.atan2(self.eps1, self.eps2)
                    if (self.w < 0.0): self.w += 360.0
                    self.werr = 1.0 # This needs fixing...
            elif param=="DIST":
                if not parts[part_index]=='*':
                    self.dist = float(parts[part_index])
                else:
                    self.dist = None
            elif param=="ASSOC":
                if not parts[part_index]=='*':
                    self.assoc = parts[part_index]
                else:
                    self.assoc = None
            elif param=="SURVEY":
                if not parts[part_index]=='*':
                    self.survey = parts[part_index]
                else:
                    self.survey = None
            elif param=="PSR":
                if not parts[part_index]=='*':
                    self.type = parts[part_index]
                else:
                    self.type = None
        self.alias = ""
    def __str__(self):
        out = ''
        if (self.name):
            out = out + "\nPulsar  B%s  (J%s)\n" % \
                  (self.name, self.jname)
        else:
            out = out + "\nPulsar J%s\n" % (self.jname)
        if (self.alias):
            out = out + "                 Alias = %s\n" % self.alias
        if (self.assoc is not None):
            out = out + "           Association = %s\n" % self.assoc
        if (self.survey is not None):
            out = out + "     Survey Detections = %s\n" % self.survey
            out = out + "    (Discoverer first)\n"
        if (self.type is not None):
            out = out + "                  Type = %s\n" % self.type
        (h, m, s) = pu.rad_to_hms(self.ra)
        serr = pc.RADTOSEC * self.raerr
        out = out + "            RA (J2000) = %s +/- %.4fs\n" % \
              (pu.coord_to_string(h, m, s), serr)
        (d, m, s) = pu.rad_to_dms(self.dec)
        serr = pc.RADTOARCSEC * self.decerr
        out = out + "           DEC (J2000) = %s +/- %.4f\"\n" % \
              (pu.coord_to_string(d, m, s), serr)
        out = out + "                (l, b) = (%.2f, %.2f)\n" % \
              (self.l, self.b)
        out = out + "          DM (cm-3 pc) = %.8g +/- %.5g\n" % \
              (self.dm, self.dmerr)
        if (self.s400 is not None):
            out = out + "        S_400MHz (mJy) = %.3g +/- %.2g\n" % \
                  (self.s400, self.s400err)
        if (self.s1400 is not None):
            out = out + "       S_1400MHz (mJy) = %.3g +/- %.2g\n" % \
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
            if self.x is not None:
                out = out + "        a*sin(i)/c (s) = %.8g +/- %.8g\n" % \
                      (self.x, self.xerr)
            if self.e is not None:
                out = out + "          Eccentricity = %.8g +/- %.8g\n" % \
                      (self.e, self.eerr)
                if (self.e > 0.0):
                    if self.w is not None:
                        out = out + "    Long of Peri (deg) = %.10g +/- %.10g\n" % \
                              (self.w, self.werr)
                    if self.To is not None:
                        out = out + "    Time of Peri (MJD) = %.12g +/- %.12g\n" % \
                              (self.To, self.Toerr)
            else:
                if self.To is not None:
                    out = out + "  T of Ascd Node (MJD) = %.12g +/- %.12g\n" % \
                          (self.To, self.Toerr)
        return out
    def pack_structs(self):
        out = struct.Struct("13s9s10s12d")
        packed = out.pack(self.jname.encode('utf-8'),
                          self.name.encode('utf-8'),
                          self.alias.lower().encode('utf-8'),
                          self.ra, self.raerr,
                          self.dec, self.decerr,
                          self.p, self.perr, self.pd, self.pderr,
                          self.dm, self.dmerr, self.pepoch, self.binary)
        if self.binary:
            if self.pb is None: self.pb = 0.0
            if self.pberr is None: self.pberr = 0.0
            if self.x is None: self.x = 0.0
            if self.xerr is None: self.xerr = 0.0
            if self.e is None: self.e = 0.0
            if self.eerr is None: self.eerr = 0.0
            if self.w is None: self.w = 0.0
            if self.werr is None: self.werr = 0.0
            if self.To is None: self.To = 0.0
            if self.Toerr is None: self.Toerr = 0.0
            packed += struct.pack("10d",
                                  self.pb, self.pberr, self.x, self.xerr,
                                  self.e, self.eerr, self.w, self.werr,
                                  self.To, self.Toerr)
        return packed

pulsars = {}
num_binaries = 0

# Read the file that was taken from the ATNF database
presto_path = os.getenv("PRESTO")
with open(os.path.join(presto_path, "lib", "psr_catalog.txt")) as csvfile:
    reader = csv.reader(csvfile, delimiter=';')
    first = next(reader)
    indices = {}
    for ii in range(len(first)):
        if first[ii]!='' and first!="#":
            indices[first[ii]] = ii
    units = next(reader)
    for row in reader:
        currentpulsar = psr(row, indices)
        pulsars[currentpulsar.jname] = currentpulsar
        if currentpulsar.binary: num_binaries += 1

# Now add the aliases to the pulsars
infile = open(os.path.join(presto_path, "lib", "aliases.txt"))
for line in infile.readlines()[1:]:
    if line[0]=='J':
        vals = line.split()
        jname = vals[0][1:]
        if jname in pulsars:
            pulsars[jname].alias = vals[2]
infile.close()

psrs = list(pulsars.values())
psrs.sort(key=attrgetter('jname'))

# Now create a new dictionary of pulsars with aliases
psr_aliases = {}
for psr in psrs:
    if psr.alias:
        psr_aliases[psr.alias] = psr

# No create a new master dictionary with all pulsar names and aliases
allpsrs = {}
for psr in psrs:
    allpsrs[psr.jname] = psr
    allpsrs["j%s" % psr.jname] = psr
    allpsrs["J%s" % psr.jname] = psr

    if psr.alias:
        allpsrs[psr.alias] = psr
    if psr.name:
        allpsrs[psr.name] = psr
        allpsrs["b%s" % psr.name] = psr
        allpsrs["B%s" % psr.name] = psr

# Add a couple important pulsars
for psr in psrs:
    if psr.jname=="1614-23":
        psr.jname=="1614-2318"
        psr.f = 29.8475387364133766
        psr.fd = -4.683105034721e-17
        psr.p, psr.pd = pu.p_to_f(psr.f, psr.fd)
        psr.x  = 1.327490
        psr.e  = 0.0
        psr.To = 52819.878171
        psr.pb = 3.15238573
        psr.w  = 0.0
        psr.dm = 52.43
        psr.l  = 351.91856
        psr.b  = 19.74496
        psr.dist = 1.80
    if psr.jname=="2204+27":
        psr.x  = 0.1
        psr.xerr  = 1.0
        psr.e  = 0.129
        psr.eerr  = 0.05
        psr.To = 57000.0
        psr.Toerr = 16.0
        psr.w  = 180.0
        psr.werr  = 180.0
        psr.pb = 32.0*24.0
        psr.pberr = 1.0

# If calling this as a main program, then write out the new pulsars.cat file
if __name__ == '__main__' :
    presto_path = os.getenv("PRESTO")
    outfilename = os.path.join(presto_path, "lib", "pulsars.cat")
    outfile = open(outfilename, "wb")
    print("Writing %d pulsars (%d binaries) to %s" % \
          (len(psrs), num_binaries, outfilename))
    for ii, psr in enumerate(psrs):
        try:
            outfile.write(psr.pack_structs())
        except:
            print(ii, psr.jname)
    outfile.close()

