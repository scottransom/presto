from __future__ import print_function
from builtins import range
from Numeric import *
from presto import *
from LeastSquares import leastSquaresFit
from orbitstuff import *

# Observation parameters
dt = 0.000125       # The duration of each data sample
N = 2**28           # The number of points in the observation
T = N*dt            # The total observation time
ctype = 'NS'        # The type of binary companion: 'WD', 'NS', or 'BH'

mpsr = mc = 0.0
psr = fake_mspsr(companion = ctype)
psr.orb.e = 0.0
psr.orb.t = 0.0
psr.orb.w = 0.0
z = 2*pi*psr.orb.x/psr.p
print('')
print('   PSR mass              =', mpsr)
print('   Companion mass        =', mc)
print('   PSR period (s)        =', psr.p)
print('   PSR frequency (hz)    =', 1.0/psr.p)
print('   Orbit period (s)      =', psr.orb.p)
print('   Orbit asini/c (lt-s)  =', psr.orb.x)
print('   Orbit eccentricity    =', psr.orb.e)
print('   Orbit angle (deg)     =', psr.orb.w)
print('   Orbit time (s)        =', psr.orb.t)
print('   Orbit Fourier Freq    =', T/psr.orb.p)
print('   Orbit z               =', z)
print('')
m = 0
kernel = presto.gen_bin_response(0.0, 1, psr.p, T, psr.orb , 
                               presto.LOWACC, m)
fftlen = next2_to_n(len(kernel))
comb= zeros(fftlen, 'F')
comb[0:len(kernel)] = kernel
minifft = rfft(spectralpower(comb))
minipow = spectralpower(minifft)
miniphs = spectralphase(minifft)
miniphs = array(miniphs, typecode='d', copy=1)
for i in range(len(minifft)-1):
    if miniphs[i+1]-miniphs[i] < 0.0:
        miniphs[i+1:] = miniphs[i+1:] + 360.0

