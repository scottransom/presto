#!/usr/bin/env python
from __future__ import print_function
import numpy as num
import sys
import presto.events as evts
from presto import kuiper
from presto.Pgplot import *


if len(sys.argv) != 2:
    print("\nusage: {} file\n".format(sys.argv[0]))
    sys.exit(1)


def calc_phases(events, f, fd):
    return num.fmod(events*(f+(0.5*fd*events)), 1.0)


events = num.loadtxt(sys.argv[1])
events.sort()
print("Read %d events from '%s'." % (events.size, sys.argv[1]))
minT, maxT = events.min(), events.max()
events -= minT
T = maxT - minT
if T > 100:
    print("Assuming that the events are in seconds (T = %.1f sec)" % T)
else:
    events *= 86400.0
    print("Assuming that the events are in days  (T = %.3f days)" % T)
    T *= 86400.0

fctr = float(sys.argv[2])
fdctr = float(sys.argv[3])
osamp = 10
df = 1.0 / (osamp * T)
dfd = 4.0 / (osamp * T * T)
nn = 101 # number of f and fd trials

print("osamp = %d, nn = %d" % (osamp, nn))
print(" fd = %g" % df)
print("dfd = %g" % dfd)

n = (nn-1)/2
fs = num.linspace(fctr-n*df, fctr+n*df, nn)
fds = num.linspace(fdctr-n*dfd, fdctr+n*dfd, nn)

kuipers = num.zeros((nn, nn), dtype=float)
htests = num.zeros((nn, nn), dtype=float)

minPk = minPh = 1.0
for ii, fd in enumerate(fds):
    print(ii)
    for jj, f in enumerate(fs):
        phases = calc_phases(events, f, fd)
        D, Pk = kuiper.kuiper_uniform_test(phases)
        h, harm = evts.Htest_exact(phases)
        Ph = evts.Hstat_prob(h)
        kuipers[ii, jj] = Pk
        htests[ii, jj] = Ph
        #print D, Pk, h, harm, Ph
        if Pk < minPk:
            minPk, fk, fdk = Pk, f, fd
        if Ph < minPh:
            minPh, fh, fdh, bestharm = Ph, f, fd, harm

print()
print("Min P(kuiper) = %.2e at f = %g, fd = %g" % (minPk, fk, fdk))
print("Min P(h-test) = %.2e at f = %g, fd = %g, (%d harmonics)" % \
      (minPh, fh, fdh, bestharm))

sigmas = num.asarray([3.0, 5.0, 7.0])
contours = num.log10(1.0-evts.gauss_sigma_to_prob(sigmas))[::-1]

print("Kuiper")
plot2d(num.log10(kuipers), fs, fds,
       labx="Frequency (Hz)", laby="F-dot (Hz/s)")
       #contours=contours, color='black', width=6)
closeplot()

print("H-test")
plot2d(num.log10(htests), fs, fds,
       labx="Frequency (Hz)", laby="F-dot (Hz/s)")
       #contours=contours, color='black')
closeplot()
