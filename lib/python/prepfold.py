import umath
import struct
import psr_utils
import Numeric
import Pgplot
from types import StringType, FloatType, IntType
from bestprof import bestprof

class foldstats:
    def __init__(self, intuple):
        (self.numdata, self.data_avg, self.data_var, self.numprof, \
         self.prof_avg, self.prof_var, self.redchi) = intuple
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
        # See if the .bestprof file is around
        try:
            self.bestprof = bestprof(filename+".bestprof")
        except IOError:
            self.bestprof = 0
        swapchar = '<' # this is little-endian
        data = infile.read(5*4)
        testswap = struct.unpack(swapchar+"i"*5, data)
        if min(umath.fabs(Numeric.asarray(testswap))) > 100000:
            swapchar = '>' # this is big-endian
        (self.numdms, self.numperiods, self.numpdots, self.nsub, self.npart) = \
                      struct.unpack(swapchar+"i"*5, data)
        (self.proflen, self.numchan, self.pstep, self.pdstep, self.dmstep, \
         self.ndmfact, self.npfact) = struct.unpack(swapchar+"i"*7, infile.read(7*4))
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
        (self.endT, self.tepoch, self.bepoch, self.avgvoverc, self.lofreq, \
         self.chan_wid, self.bestdm) = struct.unpack(swapchar+"d"*7, infile.read(7*8))
        (self.topo_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.topo_p1, self.topo_p2, self.topo_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        (self.bary_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.bary_p1, self.bary_p2, self.bary_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        (self.fold_pow, tmp) = struct.unpack(swapchar+"f"*2, infile.read(2*4))
        (self.fold_p1, self.fold_p2, self.fold_p3) = struct.unpack(swapchar+"d"*3, \
                                                                   infile.read(3*8))
        (self.orb_p, self.orb_e, self.orb_x, self.orb_w, self.orb_t, self.orb_pd, \
         self.orb_wd) = struct.unpack(swapchar+"d"*7, infile.read(7*8))
        self.dms = Numeric.asarray(struct.unpack(swapchar+"d"*self.numdms, \
                                                 infile.read(self.numdms*8)))
        if self.numdms==1:
            self.dms = self.dms[0]
        self.periods = Numeric.asarray(struct.unpack(swapchar+"d"*self.numperiods, \
                                                     infile.read(self.numperiods*8)))
        self.pdots = Numeric.asarray(struct.unpack(swapchar+"d"*self.numpdots, \
                                                   infile.read(self.numpdots*8)))
        self.numprofs = self.nsub*self.npart
        self.profs = Numeric.asarray(struct.unpack(swapchar+"d"*self.numprofs*self.proflen, \
                                                   infile.read(self.numprofs*self.proflen*8)))
        self.profs = Numeric.reshape(self.profs, (self.npart, self.nsub, self.proflen))
	self.binspersec = self.fold_p1*self.proflen
	self.chanpersub = self.numchan/self.nsub
	self.subdeltafreq = self.chan_wid*self.chanpersub
	self.losubfreq = self.lofreq + self.subdeltafreq - self.chan_wid
	self.subfreqs = Numeric.arange(self.nsub, typecode='d')*self.subdeltafreq + \
                        self.losubfreq
        self.subdelays_bins = Numeric.zeros(self.nsub, typecode='d')
        self.killed_subbands = []
        self.killed_intervals = []
        self.stats = []
        for ii in range(self.npart):
            self.stats.append([])
            for jj in range(self.nsub):
                self.stats[ii].append(foldstats(struct.unpack(swapchar+"d"*7, \
                                                              infile.read(7*8))))
        self.avgprof = sum(sum(sum(self.profs)))/self.proflen
        self.varprof = self.calc_varprof()
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
    def dedisperse(self, DM=None):
        """
        pfd.dedisperse(DM=self.bestdm):
            Rotate (internally) the profiles so that they are de-dispersed
                at a dispersion measure of DM.
        """
        if DM is None:
            DM = self.bestdm
        self.subdelays = psr_utils.delay_from_DM(DM, self.subfreqs)
	self.hifreqdelay = self.subdelays[-1]
	self.subdelays = self.subdelays-self.hifreqdelay
        # The following allows multiple DM trials
        delaybins = self.subdelays*self.binspersec - self.subdelays_bins
	new_subdelays_bins = umath.floor(delaybins+0.5)
        for ii in range(self.npart):
            for jj in range(self.nsub):
                self.profs[ii][jj] = psr_utils.rotate(self.profs[ii][jj],
                                                      int(new_subdelays_bins[jj]))
        self.subdelays_bins += new_subdelays_bins
        self.sumprof = sum(sum(self.profs))
    def kill_intervals(self, killarr):
        """
        pfd.kill_intervals(killarr):
            Set all the subintervals (internally) from the list of
                subintervals (killarr) to all zeros, effectively 'killing' them.
        """
        killchan = Numeric.zeros((self.nsub, self.proflen), typecode='d')
        for part in killarr:
            self.profs[part,:,:] *= killchan
            self.killed_intervals.append(part)
        # Update the stats
        self.avgprof = sum(sum(sum(self.profs)))/self.proflen
        self.varprof = self.calc_varprof()
    def kill_subbands(self, killarr):
        """
        pfd.kill_subbands(killarr):
            Set all the profiles (internally) from the list of
                subbands (killarr) to all zeros, effectively 'killing' them.
        """
        killchan = Numeric.zeros((self.npart, self.proflen), typecode='d')
        for sub in killarr:
            self.profs[:,sub,:] *= killchan
            self.killed_subbands.append(sub)
        # Update the stats
        self.avgprof = sum(sum(sum(self.profs)))/self.proflen
        self.varprof = self.calc_varprof()
    def plot_sumprof(self, device='/xwin'):
        """
        plot_sumprof(self, device='/xwin'):
            Plot the dedispersed and summed profile.
        """
        if not self.__dict__.has_key('subdelays'):
            print "Dedispersing first..."
            self.dedisperse()
        normprof = self.sumprof - min(self.sumprof)
        normprof /= max(normprof)
        Pgplot.plotxy(normprof, labx="Phase Bins", laby="Normalized Flux",
                      device=device)
        Pgplot.closeplot()
    def plot_intervals(self, device='/xwin'):
        """
        pfd.plot_intervals(self, device='/xwin'):
            Plot the subband-summed profiles vs time.
        """
        if not self.__dict__.has_key('subdelays'):
            print "Dedispersing first..."
            self.dedisperse()
        profs = Numeric.sum(self.profs, 1)
        # Use the same scaling as in prepfold_plot.c
        global_max = Numeric.maximum.reduce(Numeric.maximum.reduce(profs))
        min_parts = Numeric.minimum.reduce(profs, 1)
        profs = (profs-min_parts[:,Numeric.NewAxis])/global_max
        Pgplot.plot2d(profs, rangex=[0.0,self.proflen], rangey=[0.0, self.npart],
                      labx="Phase Bins", laby="Time Intervals",
                      image='antigrey', device=device)
        Pgplot.closeplot()
    def plot_subbands(self, device='/xwin'):
        """
        pfd.plot_subbands(self, device='/xwin'):
            Plot the interval-summed profiles vs subband.
        """
        if not self.__dict__.has_key('subdelays'):
            print "Dedispersing first..."
            self.dedisperse()
        profs = Numeric.sum(self.profs)
        # Use the same scaling as in prepfold_plot.c
        global_max = Numeric.maximum.reduce(Numeric.maximum.reduce(profs))
        min_subs = Numeric.minimum.reduce(profs, 1)
        profs = (profs-min_subs[:,Numeric.NewAxis])/global_max
        lof = self.lofreq - 0.5*self.chan_wid
        hif = lof + self.chan_wid*self.numchan
        Pgplot.plot2d(profs, rangex=[0.0,self.proflen], rangey=[0.0, self.nsub],
                      labx="Phase Bins", laby="Subbands",
                      laby2="Frequency (MHz)", rangey2=[lof, hif],
                      image='antigrey', device=device)
        Pgplot.closeplot()
    def calc_varprof(self):
        """
        pfd.calc_varprof(self):
            This function calculates the summed profile variance of the
                current pfd file.  Killed profiles are ignored.
        """
        varprof = 0.0
        for part in range(self.npart):
            if part in self.killed_intervals: continue
            for sub in range(self.nsub):
                if sub in self.killed_subbands: continue
                varprof += self.stats[part][sub].prof_var
        return varprof
    def calc_redchi2(self):
        """
        pfd.calc_redchi2(self):
            Return the calculated reduced-chi^2 of the current summed profile.
        """
        if not self.__dict__.has_key('subdelays'):
            print "Dedispersing first..."
            self.dedisperse()
        return sum((self.sumprof-self.avgprof)**2.0/self.varprof)/(self.proflen-1.0)
    def plot_chi2_vs_DM(self, loDM, hiDM, N=100):
        """
        pfd.plot_chi2_vs_DM(self, loDM, hiDM, N=100):
            Plot (and return) an array showing the reduced-chi^2 versus
                DM (N DMs spanning loDM-hiDM).
        """
        # Sum the profiles in time
        profs = Numeric.sum(self.profs)
        DMs = psr_utils.span(loDM, hiDM, N)
        chis = Numeric.zeros(N, typecode='f')
        subdelays_bins = self.subdelays_bins.copy()
        for ii, DM in enumerate(DMs):
            subdelays = psr_utils.delay_from_DM(DM, self.subfreqs)
            hifreqdelay = subdelays[-1]
            subdelays = subdelays - hifreqdelay
            delaybins = subdelays*self.binspersec - subdelays_bins
            last_subdelays_bins = subdelays_bins
            new_subdelays_bins = umath.floor(delaybins+0.5)
            for jj in range(self.nsub):
                profs[jj] = psr_utils.rotate(profs[jj], int(new_subdelays_bins[jj]))
            subdelays_bins += new_subdelays_bins
            sumprof = sum(profs)
            chis[ii] = sum((sumprof-self.avgprof)**2.0/self.varprof)/(self.proflen-1.0)
        # Now plot it
        Pgplot.plotxy(chis, DMs, labx="DM", laby="Reduced-\gx\u2\d")
        Pgplot.closeplot()
        return (chis, DMs)
    def plot_chi2_vs_sub(self):
        """
        pfd.plot_chi2_vs_sub(self):
            Plot (and return) an array showing the reduced-chi^2 versus
                the subband number.
        """
        # Sum the profiles in each subband
        profs = Numeric.sum(self.profs)
        # Compute the averages and variances for the subbands
        avgs = Numeric.add.reduce(profs, 1)/self.proflen
        vars = []
        for sub in range(self.nsub):
            var = 0.0
            if sub in self.killed_subbands:
                vars.append(var)
                continue
            for part in range(self.npart):
                if part in self.killed_intervals:
                    continue
                var += self.stats[part][sub].prof_var
            vars.append(var)
        chis = Numeric.zeros(self.nsub, typecode='f')
        for ii in range(self.nsub):
            chis[ii] = sum((profs[ii]-avgs[ii])**2.0/vars[ii])/(self.proflen-1.0)
        # Now plot it
        Pgplot.plotxy(chis, labx="Subband Number", laby="Reduced-\gx\u2\d")
        Pgplot.closeplot()
        return chis

if __name__ == "__main__":
    
    #testpfd = pfd("/home/ransom/tmp_pfd/M13_52724_W234_PSR_1641+3627C.pfd")
    testpfd = pfd("/home/ransom/tmp_pfd/M5_52725_W234_PSR_1518+0204A.pfd")
    print testpfd

    testpfd.kill_subbands([6,7,8,9,30,31,32,33])

    testpfd.kill_intervals([2,3,4,5,6])
        
    testpfd.plot_chi2_vs_sub()

    testpfd.plot_chi2_vs_DM(0.0, 50.0, 1000)
    
    testpfd.dedisperse()
    testpfd.plot_subbands()
    testpfd.plot_sumprof()
    print testpfd.calc_redchi2()

    testpfd.dedisperse(33.0)
    testpfd.plot_subbands()
    testpfd.plot_sumprof()
    print testpfd.calc_redchi2()

    testpfd.plot_intervals()
