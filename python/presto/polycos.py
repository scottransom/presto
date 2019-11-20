from __future__ import absolute_import
from builtins import range
from builtins import object
import os
import sys
import subprocess
from presto import parfile

import numpy as Num

# Constants
NUMCOEFFS_DEFAULT = 12
SPAN_DEFAULT = 60 # span of each polyco in minutes

# Telescope name to TEMPO observatory code conversion
telescope_to_id = {"GBT": '1', \
                   "Arecibo":' 3', \
                   "VLA": '6', \
                   "Parkes": '7', \
                   "Jodrell": '8', \
                   "GB43m": 'a', \
                   "GB 140FT": 'a', \
                   "Nancay": 'f', \
                   "Effelsberg": 'g', \
                   "WSRT": 'i', \
                   "FAST": 'k', \
                   "GMRT": 'r', \
                   "CHIME": 'y', \
                   "Geocenter": '0', \
                   "Barycenter": '@'}

# TEMPO observatory code to Telescope name conversion
id_to_telescope = {'1': "GBT", \
                   '3': "Arecibo", \
                   '6': "VLA", \
                   '7': "Parkes", \
                   '8': "Jodrell", \
                   'a': "GB43m", \
                   'a': "GB 140FT", \
                   'f': "Nancay", \
                   'g': "Effelsberg", \
                   'i': "WSRT", \
                   'k': "FAST", \
                   'r': "GMRT", \
                   'y': "CHIME", \
                   '0': "Geocenter", \
                   '@': "Barycenter"}

# Telescope name to track length (max hour angle) conversion
telescope_to_maxha = {"GBT": 12, \
                   "Arecibo": 3, \
                   "FAST": 5, \
                   "VLA": 6, \
                   "Parkes": 12, \
                   "Jodrell": 12, \
                   "GB43m": 12, \
                   "GB 140FT": 12, \
                   "Nancay": 4, \
                   "Effelsberg": 12, \
                   "WSRT": 12, \
                   "GMRT": 12, \
                   "CHIME": 1, \
                   "Geocenter": 12, \
                   "Barycenter": 12}

class polyco(object):
    def __init__(self, fileptr):
        line = fileptr.readline()
        if (line==""):
            self.psr = None
        else:
            sl = line.split()
            self.psr = sl[0]
            self.date = sl[1]
            self.UTC = sl[2]
            self.TMIDi = float(sl[3].split(".")[0])
            self.TMIDf = float("0."+sl[3].split(".")[1])
            self.TMID = self.TMIDi+self.TMIDf
            self.DM = float(sl[4])
            if (len(sl)==7):
                self.doppler = float(sl[5])*1e-4
                self.log10rms = float(sl[6])
            else:
                self.log10rms = "-"+sl[-1].split("-")[-1]
                self.doppler = float(sl[-1][:sl[-1].find(self.log10rms)])*1e-4
                self.log10rms = float(self.log10rms)
            sl = fileptr.readline().split()
            self.RPHASE = float(sl[0])
            self.F0 = float(sl[1])
            self.obs = sl[2]
            self.dataspan = int(sl[3])
            self.numcoeff = int(sl[4])
            self.obsfreq = float(sl[5])
            if (len(sl)==7):
                self.binphase = float(sl[6])
            self.coeffs = Num.zeros(self.numcoeff, 'd')
            for linenum in range(self.numcoeff//3):
                sl = fileptr.readline().split()
                self.coeffs[linenum*3+0] = float(sl[0].replace('D', 'E'))
                self.coeffs[linenum*3+1] = float(sl[1].replace('D', 'E'))
                self.coeffs[linenum*3+2] = float(sl[2].replace('D', 'E'))
            if self.numcoeff%3 != 0:  # get remaining terms if needed
                sl = fileptr.readline().split()
                nlines = self.numcoeff//3
                for coeffnum in range(len(sl)):
                    self.coeffs[nlines*3+coeffnum] = float(sl[coeffnum].replace('D', 'E'))
            self.phasepoly = Num.polynomial.polynomial.Polynomial(self.coeffs)

    def phase(self, mjdi, mjdf):
        """
        self.phase(mjdi, mjdf):
            Return the predicted pulsar phase at a given integer and frational MJD.
        """
        return self.rotation(mjdi, mjdf) % 1

    def rotation(self, mjdi, mjdf):
        """
        self.rotation(mjdi, mjdf):
            Return the predicted pulsar (fractional) rotation at a 
            given integer and fractional MJD.
        """
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        phase = self.phasepoly(DT)
        #phase = self.coeffs[self.numcoeff-1]
        #for ii in range(self.numcoeff-1, 0, -1):
        #    phase = DT*phase + self.coeffs[ii-1]
        phase += self.RPHASE + DT*60.0*self.F0
        return phase 

    def freq(self, mjdi, mjdf):
        """
        self.freq(mjdi, mjdf):
            Return the predicted pulsar spin frequency at a given integer and frational MJD.
        """
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        psrfreq = 0.0
        for ii in range(self.numcoeff-1, 0, -1):
            psrfreq = DT*psrfreq + ii*self.coeffs[ii]
        return self.F0 + psrfreq/60.0
        
class polycos(object):
    def __init__(self, psrname, filenm="polyco.dat"):
        self.psr = psrname
        self.file = filenm
        self.polycos = []
        self.TMIDs = []
        infile = open(filenm, "r")
        tmppoly = polyco(infile)
        while(tmppoly.psr):
            if (len(self.polycos)):
                if (tmppoly.dataspan != self.dataspan):
                    sys.stderr.write("Data span is changing!\n")
            else:
                self.dataspan = tmppoly.dataspan
            if (tmppoly.psr==psrname):
                self.polycos.append(tmppoly)
                self.TMIDs.append(tmppoly.TMID)
            tmppoly = polyco(infile)
        sys.stderr.write("Read %d polycos for PSR %s\n" % (len(self.polycos), psrname))
        self.TMIDs = Num.asarray(self.TMIDs)
        infile.close()
        self.validrange = 0.5*self.dataspan/1440.0

    def select_polyco(self, mjdi, mjdf):
        """
        self.select_polyco(mjdi, mjdf):
            Return the polyco number that is valid for the specified time.
        """
        goodpoly = Num.argmin(Num.fabs(self.TMIDs-(mjdi+mjdf)))
        if (Num.fabs(self.TMIDs[goodpoly]-(mjdi+mjdf)) > self.validrange):
            sys.stderr.write("Cannot find a valid polyco at %f!\n" % (mjdi+mjdf))
        return goodpoly

    def get_phase(self, mjdi, mjdf):
        """
        self.get_phase(mjdi, mjdf):
            Return the predicted pulsar phase for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].phase(mjdi, mjdf) 

    def get_rotation(self, mjdi, mjdf):
        """
        self.get_rotation(mjdi, mjdf):
            Return the predicted pulsar (fractional) rotation 
            number for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].rotation(mjdi, mjdf) 

    def get_freq(self, mjdi, mjdf):
        """
        self.get_freq(mjdi, mjdf):
            Return the predicted pulsar spin frquency for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].freq(mjdi, mjdf) 

    def get_phs_and_freq(self, mjdi, mjdf):
        """
        self.get_voverc(mjdi, mjdf):
            Return the predicted pulsar phase and spin frquency for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return (self.polycos[goodpoly].phase(mjdi, mjdf), 
                self.polycos[goodpoly].freq(mjdi, mjdf))

    def get_voverc(self, mjdi, mjdf):
        """
        self.get_voverc(mjdi, mjdf):
            Return the (approximate) topocentric v/c for the specified time.
        """
        goodpoly = self.select_polyco(mjdi, mjdf)
        return self.polycos[goodpoly].doppler

def create_polycos(parfn, telescope_id, center_freq, start_mjd, end_mjd, \
                    max_hour_angle=None, span=SPAN_DEFAULT, \
                    numcoeffs=NUMCOEFFS_DEFAULT, keep_file=False):
    """Create polycos object from a parfile.
        Inputs:
            parfn: parfile's filename, or a parfile object.
            telescope_id: The TEMPO 1-character telescope identifier.
            center_freq: The observation's center frequencies in MHz.
            start_mjd: MJD on which the polycos should start.
            end_mjd: MJD until the polycos should extend.
            max_hour_angle: The maximum hour angle as expected by tempo.
                (Default: Use default value chosen for given telescope).
            span: Span of each set of polycos in min.
                (Default: 60 min).
            numcoeffs: Number of coefficients to use.
                (Default: 12).
            keep_file: If true do not delete polyco.dat file.
                (Default: delete polyco.dat file).

        Output:
            new_polycos: a polycos object.
    """
    if type(parfn)==bytes:
        # assume parfn is a filename
        par = parfile.psr_par(parfn)
    else:
        # assume par is already a parfile.psr_par object
        par = parfn
   
    if max_hour_angle is None:
        telescope_name = id_to_telescope[telescope_id]
        max_hour_angle = telescope_to_maxha[telescope_name]
    
    tzfile = open("tz.in", "w")
    # Default parameters for prediction mode
    tzfile.write("%s %d %d %d %0.5f\n" % (telescope_id, max_hour_angle, \
                            SPAN_DEFAULT, NUMCOEFFS_DEFAULT, center_freq))
    # TEMPO ignores lines 2 and 3 in tz.in file
    tzfile.write("\n\n")
    if hasattr(par, "PSR"):
        psrname = par.PSR
    else:
        psrname = par.PSRJ
    tzfile.write("%s %d %d %d %0.5f\n" % (psrname, SPAN_DEFAULT, \
                        NUMCOEFFS_DEFAULT, max_hour_angle, center_freq)) 
    tzfile.close()
    tempo = subprocess.Popen("tempo -z -f %s" % par.FILE, shell=True, \
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                        stderr=subprocess.PIPE)
    (out, err) = tempo.communicate("%d %d\n" % (start_mjd, end_mjd))
    if err:
        raise TempoError("The following error message was encountered " \
                        "when running TEMPO to generate polycos from " \
                        "the input parfile (%s):\n\n%s\n" % (parfn, err))
    new_polycos = polycos(psrname, filenm='polyco.dat')
    # Remove files created by us and by TEMPO
    os.remove("tz.in")
    if not keep_file:
        os.remove("polyco.dat")
    return new_polycos


class TempoError(Exception):
    pass
