import umath, struct
import Numeric as num

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
        if min(umath.fabs(num.asarray(testswap))) > 100000:
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
        self.dms = num.asarray(struct.unpack(swapchar+"d"*self.numdms, infile.read(self.numdms*8)))
        if self.numdms==1:
            self.dms = self.dms[0]
        self.periods = num.asarray(struct.unpack(swapchar+"d"*self.numperiods, infile.read(self.numperiods*8)))
        self.pdots = num.asarray(struct.unpack(swapchar+"d"*self.numpdots, infile.read(self.numpdots*8)))
        self.numprofs = self.nsub*self.npart
        self.profs = num.asarray(struct.unpack(swapchar+"d"*self.numprofs*self.proflen, infile.read(self.numprofs*self.proflen*8)))
        self.profs = num.reshape(self.profs, (self.npart, self.nsub, self.proflen))
	self.binspersec = self.fold_p1*self.proflen
	self.chanpersub = self.numchan/self.nsub
	self.subdeltafreq = self.chan_wid*self.chanpersub
	self.losubfreq = self.lofreq + subdeltafreq - self.chan_wid
	self.subfreqs = num.arange(self.nsub, typecode='d')*self.subdeltafreq + self.losubfreq
	self.subdelays = delay_from_DM(self.bestdm, self.subfreqs)
	self.hifreqdelay = self.subdelays[-1]
	self.subdelays = self.subdelays-self.hifreqdelay
	self.subdelays_bins = floor(self.subdelays*self.binspersec+0.5)
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

    def dedisperse(self, DM=self.bestdm):
        
