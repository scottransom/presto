from math import sqrt, sin, pi
from Numeric import *
from RandomArray import random, randint, normal
from presto import rfft, TWOPI, spectralpower, maximize_r
from Pgplot import plotxy, plotbinned, closeplot
from Statistics import average, standardDeviation, histogram

# This code tests the statistics of photon based power
# spectra and signals.  It is useful for checking
# x-ray statistics equations...

Ntrials = 50                     # Number of trials per set of inputs
Nbinss = 2**arange(12, 16)       # Number of points in FFT
Nphots = 200 * arange(1, 6)      # Number of photons in time series
pfracs = arange(0.1, 1.1, 0.2)   # Pulsed cts / background cts
pulsetype = 'Sine'               # Can be either 'Sine' or 'Gaussian'
width = 0.1                      # +/-1 sigma width of Gaussian pulse in phase

def sinc(x):
    if (x==0.0):
        return 1.0
    else:
        return sin(x)/x

def secant(func, oldx, x, TOL=1e-6): # f(x)=func(x)
    """
    Similar to Newton's method, but the derivative is estimated
    by divided difference using only function calls.  A root is
    estimated by x = x - f(x) (x - oldx)/(f(x) - f(oldx))
    where oldx = x[i-1] and x = x[i].
    """
    oldf, f = func(oldx), func(x)
    if (abs(f) > abs(oldf)):         # swap so that f(x) is closer to 0
        oldx, x = x, oldx
        oldf, f = f, oldf
    count = 0
    while 1:
        dx = f * (x - oldx) / float(f - oldf)
        if abs(dx) < TOL * (1 + abs(x)): return x - dx
        oldx, x = x, x - dx
        oldf, f = f, func(x)
        count = count + 1
        print "secant(%d): x=%s, f(x)=%s" % (count, x, f)

class trials:
    def __init__(self, Nbins, Nphot, pfrac):
        self.Nphot = Nphot
        self.Nbins = Nbins
        self.pfrac = pfrac
        self.freq = Nbins / 8
        self.pulsedphot = int(pfrac * Nphot)
        self.backgdphot = Nphot - self.pulsedphot
        self.avg_theo_pows = 1.0 + pfrac**2.0 * Nphot / 4.0 * \
                             sinc(pi * self.freq / Nbins)**2.0
        self.sig_theo_pows = sqrt(2 * self.avg_theo_pows)
        self.meas_pows = []
    def calc_time_series(self):
        data = zeros(self.Nbins, 'f')
        data.savespace()
        # Place the non-pulsed photons
        points = randint(0, self.Nbins, self.backgdphot)
        for point in points:
            data[point] = data[point] + 1.0
        # Place the pulsed photons
        pulses = randint(0, self.freq, self.pulsedphot)
        if pulsetype=='Sine':
            x = arange(10001, typecode='d') * TWOPI/10000
            coscdf = (x + sin(x))/TWOPI
            uvar = random(self.pulsedphot)
            phases = take((x + TWOPI * random()) % TWOPI,
                          searchsorted(coscdf, uvar)) / TWOPI
            # plotxy(coscdf, x)
            # closeplot()
            # hist = histogram(phases, 100, [0.0, TWOPI])
            # plotbinned(hist[:,1], hist[:,0])
            # closeplot()
        elif pulsetype=='Gaussian':
            phases = normal(0.0, width, self.pulsedphot) % 1.0 + random()
        points = ((pulses + phases) / self.freq * self.Nbins).astype('i')
        for point in points:
            data[point] = data[point] + 1.0
            # tmp = average(data)
            # print 'Average Data = ', tmp
            # print '   sqrt(avg) = ', sqrt(tmp)
            # print ' StdDev Data = ', standardDeviation(data)
        return data
    def do_trials(self, Ntrials):
        for trial in xrange(Ntrials):
            data = self.calc_time_series()
            ft = rfft(data)
            # pows = spectralpower(ft)/ft[0].real
            # print 'Number of counts =', sum(data)
            # print 'Max power =', max(pows)
            # print 'Max  freq =', argmax(pows)
            # plotxy(pows)
            # closeplot()
            [hipow, hifreq, rderivs] = \
                    maximize_r(ft, self.freq, norm = ft[0].real)
            self.meas_pows.append(hipow)
        self.calc_props()
    def calc_props(self):
        self.avg_meas_pows = average(self.meas_pows)
        self.sig_meas_pows = standardDeviation(self.meas_pows)
        self.avg_diff_pows = self.avg_meas_pows - self.avg_theo_pows
        self.sig_diff_pows = self.sig_meas_pows
    def print_props(self):
        print 'Nbins = %-7d  Nphot = %-7d  pfrac = %-.2f' % \
              (self.Nbins, self.Nphot, self.pfrac)
        print '----------------------------------------------'
        print ' Theo Power = %-6.2f +/- %-5.2f' % (self.avg_theo_pows,
                                                   self.sig_theo_pows)
        print ' Meas Power = %-6.2f +/- %-5.2f' % (self.avg_meas_pows,
                                                   self.sig_meas_pows)
        print ' Difference = %-6.2f +/- %-5.2f' % (self.avg_diff_pows,
                                                   self.sig_diff_pows)
        print ''

# The main loop

results = []
for Nbins in Nbinss:
    for Nphot in Nphots:
        for pfrac in pfracs:
            current_trial = trials(Nbins, Nphot, pfrac)
            current_trial.do_trials(Ntrials)
            current_trial.print_props()
            results.append(current_trial)
