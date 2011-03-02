import numpy as Num
import sys

class polyco:
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
    def phase(self, mjdi, mjdf):
        """
        self.phase(mjdi, mjdf):
            Return the predicted pulsar phase at a given integer and frational MJD.
        """
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        phase = self.coeffs[self.numcoeff-1]
        for ii in range(self.numcoeff-1, 0, -1):
            phase = DT*phase + self.coeffs[ii-1]
        phase += self.RPHASE + DT*60.0*self.F0
        return phase - Num.floor(phase)
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
        
class polycos:

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
