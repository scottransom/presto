import Numeric
import umath
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
	    self.doppler = float(sl[4])
	    self.log10rms = float(sl[5])
	    sl = fileptr.readline().split()
	    self.RPHASE = float(sl[0])
	    self.F0 = float(sl[1])
	    self.obs = sl[2]
	    self.dataspan = int(sl[3])
	    self.numcoeff = int(sl[4])
	    self.obsfreq = float(sl[5])
            if (len(sl)==7):
                self.binphase = float(sl[6])
	    self.coeffs = Numeric.zeros(self.numcoeff, 'd')
	    for linenum in range(self.numcoeff/3):
		sl = fileptr.readline().split()
		self.coeffs[linenum*3+0] = float(sl[0])
		self.coeffs[linenum*3+1] = float(sl[1])
		self.coeffs[linenum*3+2] = float(sl[2])
    def phase(self, mjdi, mjdf):
        DT = ((mjdi-self.TMIDi)+(mjdf-self.TMIDf))*1440.0
        phase = self.coeffs[self.numcoeff-1]
        for ii in range(self.numcoeff-1, 0, -1):
            phase = DT*phase + self.coeffs[ii-1]
        phase += self.RPHASE + DT*60.0*self.F0
        return phase - umath.floor(phase)
    def freq(self, mjdi, mjdf):
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
        self.TMIDs = Numeric.asarray(self.TMIDs)
        infile.close()
        self.validrange = 0.5*self.dataspan/1440.0
    def get_phs_and_freq(self, mjdi, mjdf):
        goodpoly = Numeric.argmin(umath.fabs(self.TMIDs-(mjdi+mjdf)))
        if (umath.fabs(self.TMIDs[goodpoly]-(mjdi+mjdf)) > self.validrange):
            sys.stderr.write("Cannot find a valid polyco at %f!\n" % (mjdi+mjdf))
        return (self.polycos[goodpoly].phase(mjdi, mjdf), 
                self.polycos[goodpoly].freq(mjdi, mjdf))
