#!/usr/bin/env python
import struct, getopt, sys
from umath import *
from Numeric import *
from FFT import *
from fftfit import *
from types import StringType, FloatType, IntType
from Pgplot import *

PI = float('3.1415926535897932384626433832795028841971693993751')
TWOPI = PI*2.0
SECPERDAY = 86400.0
scopes = {'GBT':'1', 'Arecibo':'3', 'Parkes':'7', 'GMRT': 'r'}

class infodata:
    def __init__(self, filenm):
        infile = open(filenm)
        for line in infile.readlines():
            if line.startswith(" Data file name"):
                self.basenm = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                self.telescope = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                self.instrument = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                self.object = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                self.RA = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                self.DEC = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                self.observer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                self.epoch = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                self.bary = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                self.N = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                self.dt = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                self.breaks = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Type of observation"):
                self.waveband = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                self.beam_diam = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                self.DM = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                self.lofreq = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                self.BW = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                self.numchan = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                self.chan_width = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                self.analyzer = line.split("=")[-1].strip()
                continue
        infile.close()
        
class foldstats:
    def __init__(self, intuple):
        (self.numdata, self.data_avg, self.data_var, self.numprof, self.prof_avg, self.prof_var, self.redchi) = intuple
    def __str__(self):
        out = ""
        for k, v in self.__dict__.items():
            if k[:2]!="__":
                out += "%10s = '%s' " % (k, v)
        out += '\n'
        return out

class pfd:
    def __init__(self, filename):
        infile = open(filename, "rb")
        swapchar = '<' # this is little-endian
        data = infile.read(5*4)
        testswap = struct.unpack(swapchar+"i"*5, data)
        if min(abs(asarray(testswap))) > 100000:
            swapchar = '>' # this is big-endian
        (self.numdms, self.numperiods, self.numpdots, self.nsub, self.npart) = \
                      struct.unpack(swapchar+"i"*5, data)
        (self.proflen, self.numchan, self.pstep, self.pdstep, self.dmstep, self.ndmfact, self.npfact) = \
                       struct.unpack(swapchar+"i"*7, infile.read(7*4))
        self.filenm = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        self.candnm = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        self.telescope = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        self.pgdev = infile.read(struct.unpack(swapchar+"i", infile.read(4))[0])
        test = infile.read(16)
        has_posn = 1
        for ii in range(16):
            if test[ii] not in '0123456789:.-\0':
                has_posn = 0
                break
        if has_posn:
            self.rastr = test[:test.find('\0')]
            test = infile.read(16)
            self.decstr = test[:test.find('\0')]
            (self.dt, self.startT) = struct.unpack(swapchar+"dd", infile.read(2*8))
        else:
            self.rastr = "Unknown"
            self.decstr = "Unknown"
            (self.dt, self.startT) = struct.unpack(swapchar+"dd", test)
        (self.endT, self.tepoch, self.bepoch, self.avgvoverc, self.lofreq, self.chan_wid, self.bestdm) = \
                    struct.unpack(swapchar+"d"*7, infile.read(7*8))
        (self.topo_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.topo_p1, self.topo_p2, self.topo_p3) = struct.unpack(swapchar+"d"*3, infile.read(3*8))
        (self.bary_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.bary_p1, self.bary_p2, self.bary_p3) = struct.unpack(swapchar+"d"*3, infile.read(3*8))
        (self.fold_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.fold_p1, self.fold_p2, self.fold_p3) = struct.unpack(swapchar+"d"*3, infile.read(3*8))
        (self.orb_p, self.orb_e, self.orb_x, self.orb_w, self.orb_t, self.orb_pd, self.orb_wd) = \
                     struct.unpack(swapchar+"d"*7, infile.read(7*8))
        self.dms = asarray(struct.unpack(swapchar+"d"*self.numdms, infile.read(self.numdms*8)))
        if self.numdms==1:
            self.dms = self.dms[0]
        self.periods = asarray(struct.unpack(swapchar+"d"*self.numperiods, infile.read(self.numperiods*8)))
        self.pdots = asarray(struct.unpack(swapchar+"d"*self.numpdots, infile.read(self.numpdots*8)))
        self.numprofs = self.nsub*self.npart
        self.profs = asarray(struct.unpack(swapchar+"d"*self.numprofs*self.proflen, infile.read(self.numprofs*self.proflen*8)))
        self.profs = reshape(self.profs, (self.npart, self.nsub, self.proflen))
        self.stats = []
        for ii in range(self.npart):
            self.stats.append([])
            for jj in range(self.nsub):
                self.stats[ii].append(foldstats(struct.unpack(swapchar+"d"*7, infile.read(7*8))))
        infile.close()
    def __str__(self):
        out = ""
        for k, v in self.__dict__.items():
            if k[:2]!="__":
                if type(self.__dict__[k]) is StringType:
                    out += "%10s = '%s'\n" % (k, v)
                elif type(self.__dict__[k]) is IntType:
                    out += "%10s = %d\n" % (k, v)
                elif type(self.__dict__[k]) is FloatType:
                    out += "%10s = %-20.15g\n" % (k, v)
        return out

class bestprof:
    def __init__(self, filenm):
        infile = open(filenm)
        self.topo = 0
        self.profile = []
        for line in infile.readlines():
            if line[0]=="#":
                if line.startswith("# Input file"):
                    self.datnm = line.split("=")[-1][:-1]
                    continue
		if line.startswith("# Candidate"):
		    if line.startswith("# Candidate        =  PSR_"):
			self.psr = line.split("=")[-1].split("_")[1][:-1]
			continue
		    else:
			self.psr = None
                if line.startswith("# T_sample"):
                    self.dt = float(line.split("=")[-1])
                    continue
                if line.startswith("# Data Folded"):
                    self.N = float(line.split("=")[-1])
                    continue 
                if line.startswith("# Epoch_topo"):
                    try:
                        self.epochi = float(line.split("=")[-1].split(".")[0])
                        self.epochf = float("0."+line.split("=")[-1].split(".")[1])
                        self.topo = 1
                    except ValueError:
                        pass
                    continue
                if (not self.topo and line.startswith("# Epoch_bary")):
                    try:
                        self.epochi = float(line.split("=")[-1].split(".")[0])
                        self.epochf = float("0."+line.split("=")[-1].split(".")[1])
                    except ValueError:
                        pass
                if ((not self.topo and line.startswith("# P_bary")) or
                     (self.topo and line.startswith("# P_topo"))):
                    self.p0 = float(line.split("=")[-1].split("+")[0])/1000.0
                    continue
                if ((not self.topo and line.startswith("# P'_bary")) or
                     (self.topo and line.startswith("# P'_topo"))):
                    self.p1 = float(line.split("=")[-1].split("+")[0])
                    continue
                if ((not self.topo and line.startswith("# P''_bary")) or
                     (self.topo and line.startswith("# P''_topo"))):
                    self.p2 = float(line.split("=")[-1].split("+")[0])
                    continue
            else:
                self.profile.append(float(line.split()[-1]))
        infile.close()
        self.T = self.dt*self.N
        self.proflen = len(self.profile)

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
	    self.coeffs = zeros(self.numcoeff, 'd')
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
        return phase - floor(phase)
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
        self.TMIDs = asarray(self.TMIDs)
        infile.close()
        self.validrange = 0.5*self.dataspan/1440.0
    def get_phs_and_freq(self, mjdi, mjdf):
        goodpoly = argmin(fabs(self.TMIDs-(mjdi+mjdf)))
        if (fabs(self.TMIDs[goodpoly]-(mjdi+mjdf)) > self.validrange):
            sys.stderr.write("Cannot find a valid polyco at %f!\n" % (mjdi+mjdf))
        return (self.polycos[goodpoly].phase(mjdi, mjdf), 
                self.polycos[goodpoly].freq(mjdi, mjdf))

def delay_from_DM(DM, freq_emitted):
    """
    Return the delay in seconds caused by dispersion, given
    a Dispersion Measure (DM) in cm-3 pc, and the emitted
    frequency (freq_emitted) of the pulsar in MHz.
    """
    if (freq_emitted > 0.0):
        return DM/(0.000241*freq_emitted*freq_emitted)
    else:
        return 0.0

def indexed(sequence):
    """
    indexed(sequence):
        Return a sequence woven with its indices as (ii, sequence[ii]).
    """
    return [(ii,sequence[ii]) for ii in range(len(sequence))]

def p_to_f(p, pd, pdd=None):
   """
   p_to_f(p, pd, pdd=None):
      Convert period, period derivative and period second
      derivative to the equivalent frequency counterparts.
      Will also convert from f to p.
   """
   f = 1.0 / p
   fd = -pd / (p * p)
   if (pdd==None):
       return [f, fd]
   else:
       if (pdd==0.0):
           fdd = 0.0
       else:
           fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
       return [f, fd, fdd]

def measure_phase(profile, template):
    c,amp,pha = cprof(template)
    pha.savespace()
    pha1 = pha[0]
    pha = fmod(pha-arange(1,len(pha)+1)*pha1,TWOPI)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood

def calc_phs(closeMJD, epoch, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
    t = (closeMJD-epoch)*SECPERDAY
    return fmod(t*(f0 +
                   t*(f1/2.0 +
                      t*(f2/6.0 +
                         t*(f3/24.0 +
                            t*(f4/120.0 +
                               t*f5/720.0))))), 1.0)

def calc_freq(closeMJD, epoch, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
    t = (closeMJD-epoch)*SECPERDAY
    return f0 + t*(f1 +
                   t*(f2/2.0 +
                      t*(f3/6.0 +
                         t*(f4/24.0 +
                            t*f5/120.0))))

def calc_t0(closeMJD, epoch, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
    phs = calc_phs(closeMJD, epoch, f0, f1, f2, f3, f4, f5)
    p = 1.0/calc_freq(closeMJD, epoch, f0, f1, f2, f3, f4, f5)
    return closeMJD-phs*p/SECPERDAY

def gaussian_profile(N, phase, fwhm):
    """
    gaussian_profile(N, phase, fwhm):
        Return a gaussian pulse profile with 'N' bins and
        an integrated 'flux' of 1 unit.
            'N' = the number of points in the profile
            'phase' = the pulse phase (0-1)
            'fwhm' = the gaussian pulses full width at half-max
        Note:  The FWHM of a gaussian is approx 2.35482 sigma
    """
    sigma = fwhm/2.35482
    mean = phase%1.0
    phsval = arange(N, typecode='d')/float(N)
    if (mean < 0.5):
        phsval = where(phsval>mean+0.5, phsval-1.0, phsval)
    else:
        phsval = where(phsval<mean-0.5, phsval+1.0, phsval)
    return exp(-0.5*((phsval-mean)/sigma)**2)/(sigma*sqrt(TWOPI))

def write_princeton_toa(toa, toaerr, freq, dm, obs='@', name=' '*13):
    """
    Princeton Format

    columns     item
    1-1     Observatory (one-character code) '@' is barycenter
    2-2     must be blank
    16-24   Observing frequency (MHz)
    25-44   TOA (decimal point must be in column 30 or column 31)
    45-53   TOA uncertainty (microseconds)
    69-78   DM correction (pc cm^-3)
    """
    if dm!=0.0:
        print obs+" %13s %8.3f %19.13f %8.2f              %9.4f" % \
              (name, freq, toa, toaerr, dm)
    else:
        print obs+" %13s %8.3f %19.13f %8.2f" % \
              (name, freq, toa, toaerr)
        
def rotate(arr, bins):
    """
    rotate(arr, bins):
	Return an array rotated by 'bins' places to the left
    """
    bins = bins % len(arr)
    if bins==0:
        return arr
    else:
	return concatenate((arr[bins:], arr[:bins]))

def corr(profile, template):
    """
    corr(profile, template):
	Cross-correlate (using FFTs) a 'profile' and a 'template'.
    """
    return inverse_real_fft(real_fft(template) * real_fft(profile))

def maxphase(profile, template):
    """
    maxphase(profile, template):
	Return the phase offset required to get the 'profile' to best
            match the 'template'.
    """
    return float(argmax(corr(profile, template))) / len(profile)

def linear_interpolate(vector, zoom=10):
    """
    linear_interpolate(vector, zoom=10):
	Linearly interpolate 'vector' by a factor of 'zoom'.
    """
    n = len(vector)
    ivect = zeros(zoom*n, typecode='d')
    nvect = concatenate((vector, vector[:1]))
    ivals = arange(zoom, typecode='d')/zoom
    loy = nvect[0]
    for ii in range(n):
	hiy = nvect[ii+1]
	ivect[ii*zoom:(ii+1)*zoom] = ivals*(hiy-loy) + loy
	loy = hiy
    return ivect

def measure_phase_corr(profile, template, zoom=10):
    """
    measure_phase_corr(profile, template, zoom=10):
	Return the phase offset required to get the 'profile' to best
            match the 'template', each of which has been interpolated
	    by a factor of 'zoom'.
    """
    zoomprof = zoomtemp = zoom
    if (len(template) != len(profile)):
	if (len(template)%len(profile) == 0):
	    zoomprof = zoom*len(template)/len(profile)
	else:
	    print "Warning!:  The lengths of the template (%d) and profile (%d)" % \
	          (len(template), len(profile))
	    print "           are not the same!"
    itemp = linear_interpolate(rotate(template, argmax(template)), zoomtemp)
    iprof = linear_interpolate(profile, zoomprof)
    return maxphase(iprof, itemp)

def usage():
    print """
usage:  get_TOAs.py [options which must include -t or -g] pfd_file
  [-h, --help]                       : Display this help
  [-s numsub, --subbands=numsub]     : Divide the fold into numsub subbands
  [-n numTOAs, --numtoas=numTOAs]    : Divide the fold into numTOAs parts
  [-d DM, --dm=DM]                   : Re-combine subbands at DM
  [-f, --FFTFITouts]                 : Print all FFTFIT outputs and errors
  [-g gausswidth, --gaussian=width]  : Use a Gaussian template of FWHM width
  [-t templateprof, --template=prof] : The template .bestprof file to use
  [-k subs_list, --kill=subs_list]   : List of subbands to ignore
  [-e, --event]                      : The .pfd file was made with events
  pfd_file                           : The .pfd file containing the folds

  The program generates TOAs from a .pfd file using Joe Taylor's
  FFTFIT program. The TOAs are output to STDOUT.  Typically, the .pfd
  file is created using prepfold with the "-timing" flag and an
  appropriate .par file on either a topocentric time series or raw
  telescope data.  But barycentric folds or folds of barycentered
  events are also acceptable.  The most important thing about the
  fold, though, is that it must have been made using "-nosearch"! 
  (Note: "-timing" implies "-nosearch")
  
  A typical example would be something like:
      
      get_TOAs.py -n 30 -t myprof.bestprof -k 0,20-23 myprof.pfd | \\
          tail -28 >> good.tim
      
  which would extract 30 TOAs (the default number of slices or parts
  in time for "prepfold -timing" is 60) from a fold made from some raw
  radio telescope data.  The command would ignore (i.e. zero-out)
  subbands 0, 20, 21, 22, and 23 (e.g.  due to interference) and then
  ignore the first 2 TOAs with the tail command.
  
  If you don't specify "-n", the default number of parts in the fold
  is assumed, but if you don't specify "-s", all the subbands (if any
  are present) are integrated together.
  
  If you specify the "-f" flag, an additional line of output is
  displayed for each TOA that shows the "b +/- berr" and "SNR +/-
  SNRerr" params from FFTFIT.
  
"""

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hefs:n:d:g:t:o:k:e:",
                                   ["help", "event", "FFTFITouts", "subbands=", 
				    "numtoas=", "dm=", "gaussian=", "template=",
                                    "offset=", "kill="])
                                    
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    lowfreq = None
    DM = 0.0
    gaussianwidth = 0.1
    templatefilenm = None
    numchannels = 1
    numsubbands = 1
    numtoas = 1
    otherouts = 0
    offset = 0.0
    events = 0
    kill = []
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--FFTFITouts"):
	    otherouts = 1
	if o in ("-e", "--event"):
	    lowfreq = 0.0
	    DM = 0.0
	    events = 1
	if o in ("-s", "--subbands"):
            numsubbands = int(a)
        if o in ("-n", "--numtoas"):
            numtoas = int(a)
            if numtoas==0:
                sys.exit()
        if o in ("-d", "--dm"):
            DM = float(a)
        if o in ("-g", "--gaussian"):
            gaussianwidth = float(a)
        if o in ("-t", "--template"):
            templatefilenm = a
        if o in ("-o", "--offset"):
            offset = float(a)
        if o in ("-k", "--kill"):
            for subs in a.split(','):
                if (subs.find("-") > 0):
                    lo, hi = subs.split("-")
                    kill.extend(range(int(lo), int(hi)+1))
                else:
                    kill.append(int(subs))

    # Read key information from the bestprof file
    fold = bestprof(sys.argv[-1]+".bestprof")

    # Read the prepfold output file and the binary profiles
    fold_pfd = pfd(sys.argv[-1])
    # Over-ride the DM that was used during the fold
    if (DM!=0.0):
        fold_pfd.bestdm = DM
    if (fold_pfd.numchan==1):
        if (DM==0.0 and numchannels==1):
            idata = infodata(fold_pfd.filenm[:fold_pfd.filenm.rfind('.')]+".inf")
	    if (events):
		fold_pfd.bestdm = 0.0
		fold_pfd.numchan = 1
	    else:
		fold_pfd.bestdm = idata.DM
		fold_pfd.numchan = idata.numchan
	else:
            fold_pfd.bestdm = DM
            fold_pfd.numchan = numchannels

    # Read the template profile
    if templatefilenm is not None:
        template_fold = bestprof(templatefilenm)
        template = asarray(template_fold.profile)
        (mn, mx) = (min(template), max(template))
        template = (template-mn)/(mx-mn)
    else:
        template = gaussian_profile(fold_pfd.proflen, 0.0, gaussianwidth)
        template = template / max(template)

    timestep_sec = fold.T / numtoas
    timestep_day = timestep_sec / SECPERDAY

    # PRESTO de-disperses at the high frequency channel so determine a
    # correction to the middle of the band
    if not events:
	binspersec = fold_pfd.fold_p1*fold_pfd.proflen
	chanpersub = fold_pfd.numchan/fold_pfd.nsub
	subdeltafreq = fold_pfd.chan_wid*chanpersub
	losubfreq = fold_pfd.lofreq + subdeltafreq - fold_pfd.chan_wid
	subfreqs = arange(fold_pfd.nsub, typecode='d')*subdeltafreq + losubfreq
	subdelays = delay_from_DM(fold_pfd.bestdm, subfreqs)
	hifreqdelay = subdelays[-1]
	subdelays = subdelays-hifreqdelay
	subdelays_bins = floor(subdelays*binspersec+0.5)
	subpersumsub = fold_pfd.nsub/numsubbands
	# Calculate the center of the summed subband freqs and delays
	sumsubfreqs = (arange(numsubbands)+0.5)*subpersumsub*subdeltafreq + \
                      (fold_pfd.lofreq-0.5*fold_pfd.chan_wid)
	sumsubdelays = (delay_from_DM(fold_pfd.bestdm, sumsubfreqs)-hifreqdelay)/SECPERDAY
    else:
	subdelays = asarray([0.0])
	subdelays_bins = asarray([0.0])
	sumsubdelays = asarray([0.0])
	subfreqs = asarray([0.0])
	sumsubfreqs = asarray([0.0])
	
    # Shift the profiles by the correct amount required by the best DM
    for ii in range(fold_pfd.npart):
        for jj in range(fold_pfd.nsub):
            if jj in kill:
                fold_pfd.profs[ii][jj] *= 0.0
            else:
                fold_pfd.profs[ii][jj] = rotate(fold_pfd.profs[ii][jj],
                                                int(subdelays_bins[jj]))
    #fold.epochf += DM_delay/SECPERDAY
    if fold.epochf > 1.0:
        fold.epochf -= 1.0
        fold.epochi += 1
    fold.epoch = fold.epochi+fold.epochf

    if (not fold.topo):
        obs = '@'  # Solarsystem Barycenter
    else:
        try: obs = scopes[fold_pfd.telescope.split()[0]]
	except KeyError:  print "Unknown telescope!!!"

    # Read the polyco file
    if (fold.psr and fold.topo):
        pcs = polycos(fold.psr, sys.argv[-1]+".polycos")
        (fold.phs0, fold.f0) = pcs.get_phs_and_freq(fold.epochi, fold.epochf)
        fold.f1 = fold.f2 = 0.0
    else:
        pcs = None
        fold.phs0 = 0.0
        (fold.f0, fold.f1, fold.f2) = p_to_f(fold.p0, fold.p1, fold.p2)

    # Combine the sub-integration profiles
    profs = zeros((numtoas, numsubbands, fold_pfd.proflen), 'd')
    dp = fold_pfd.npart/numtoas
    ds = fold_pfd.nsub/numsubbands
    for ii in range(numtoas):
        # Combine the subbands if required
        if (fold_pfd.nsub > 1):
            for jj in range(numsubbands):
                subprofs = add.reduce(fold_pfd.profs[:,jj*ds:(jj+1)*ds], 1)
                # Combine the time intervals
                profs[ii][jj] = add.reduce(subprofs[ii*dp:(ii+1)*dp])
        else:
            profs[ii][0] = add.reduce(fold_pfd.profs[ii*dp:(ii+1)*dp,0])

    if (0):
        for ii in range(numtoas):
            profs[ii][0] = (profs[ii][0]-min(profs[ii][0]))/max(profs[ii][0])
        plot2d(profs[:,0], arange(64.0)/64.0, arange(numtoas, typecode='d'), image='antigrey')
        closeplot()

    # Calculate the TOAs
    for ii in range(numtoas):
        if (pcs is None):
            midtime = fold.epoch + (ii+0.5)*timestep_day
            p = 1.0/calc_freq(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
            t0 = calc_t0(midtime, fold.epoch, fold.f0, fold.f1, fold.f2)
        else:
            mjdf = fold.epochf + (ii+0.5)*timestep_day
            (phs, f0) = pcs.get_phs_and_freq(fold.epochi, mjdf)
            phs -= fold.phs0
            p = 1.0/fold.f0
            t0 = fold.epochi+mjdf - phs*p/SECPERDAY
        for jj in range(numsubbands):
            prof = profs[ii][jj]
            try:
		if (len(template)==len(prof)):
		    shift,eshift,snr,esnr,b,errb,ngood = measure_phase(prof, template)
		    tau, tau_err = shift/fold_pfd.proflen, eshift/fold_pfd.proflen
		else:
		    shift = 0.0
		    eshift = 999.0
		if (fabs(shift) < 1e-7 and
		    fabs(eshift-999.0) < 1e-7):
		    # Not enough structure in the template profile for FFTFIT
		    # so use time-domain correlations instead
		    tau = measure_phase_corr(prof, template)
		    # This needs to be changed
		    tau_err = 0.1/len(prof)
		write_princeton_toa(t0+(tau*p+offset)/SECPERDAY+sumsubdelays[jj],
                                    tau_err*p*1000000.0,
                                    sumsubfreqs[jj], fold_pfd.bestdm, obs=obs)
		if (otherouts):
		    print "FFTFIT results:  b = %.4g +/- %.4g   SNR = %.4g +/- %.4g" % \
		        (b, errb, snr, esnr)
	    except ValueError, fftfit.error:
                pass
