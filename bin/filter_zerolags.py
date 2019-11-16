#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import numpy as N
import sys, scipy.io, scipy.signal


if len(sys.argv) != 2:
    print("\nusage: {} file\n".format(sys.argv[0]))
    sys.exit(1)


plot=0

infilenm = sys.argv[1]
basename = infilenm[:infilenm.find(".zerolags")]

dt = 0.00008192
flo      = 2.0   # cutoff freq in Hz
passband = 0.8   # Fractional freq where passband ends
stopband = 1.2   # Fractional freq where stopband starts
max_pass_atten = 3.0  # dB
min_stop_atten = 30.0 # dB

zls = N.fromfile(infilenm, 'f')
numpts = len(zls)

if (plot):
    from presto.Pgplot import *
    plotxy(zls)

T = numpts*dt
fnyq = numpts/(T*2)
cutoff = flo/fnyq

# Determine an "average" standard deviation
stds = []
for ii in range(100):
    loind = int(N.random.rand() * (numpts-1001))
    hiind = loind + 1000
    stds.append(N.std(zls[loind:hiind]))
goodstd = N.median(stds)

# First, detrend the data in a piecewise linear fashion
# where the pieces are defined by jumps in the data

num_bps = 0
max_num_break_points = 100
break_points = N.zeros(max_num_break_points)

dzls = N.fabs(zls[1:] - zls[:-1])
argsort_jumps = dzls.argsort()[::-1]

ii = 0
dind = 200
while (ii < numpts-1):
    index = argsort_jumps[ii]
    # Don't allow the breakpoints to be within 100 of each other
    if num_bps and min(N.fabs(break_points[:num_bps] - index)) < 100:
        ii += 1
        continue
    if index > 100:
        lomean = N.mean(zls[index-dind:index-20])
    else:
        lomean = N.mean(zls[:index])
    if index < numpts-1-dind:
        himean = N.mean(zls[index+20:index+dind])
    else:
        himean = N.mean(zls[index:])
    if N.fabs(himean - lomean) > goodstd:
        break_points[num_bps] = index
        num_bps += 1
        if num_bps == max_num_break_points: break
    ii += 1
    if dzls[index] < 3.0 * goodstd:  break
if (num_bps):
    break_points = break_points[:num_bps]
    break_points.sort()
    detrend_zls = scipy.signal.detrend(zls, bp=break_points)
    print("%s: Found %d breakpoints for detrending "%(basename, num_bps), break_points)
else:
    detrend_zls = scipy.signal.detrend(zls)
    print("%s: Found 0 breakpoints for detrending"%basename)
    
# Now high-pass filter the data to get rid of the not-so-drastic
# power fluctuations

logood = numpts/5
newzls = N.zeros(numpts+logood, 'f')
newzls[:logood] += detrend_zls[:logood][::-1]
newzls[logood:] += detrend_zls

b, a = scipy.signal.iirdesign(passband*cutoff, stopband*cutoff,
                              max_pass_atten, min_stop_atten,
                              analog=0, ftype='ellip', output='ba')
filtered_zls = scipy.signal.lfilter(b, a, newzls)[logood:]

if (plot):
    plotxy(filtered_zls, color='red')
    closeplot()

# This is the total offset that we will apply to the data
offset = (zls - detrend_zls) + filtered_zls
outfilenm = basename+".offset"
outfile = open(outfilenm, 'wb')
scipy.io.fwrite(outfile, numpts, offset, 'f')
outfile.close()
