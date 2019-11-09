import numpy as np
from presto import presto
from numpy.random import standard_normal as norm
from numpy.random import uniform
import time

N = 2**17
noiseamp = 1.0
numharm = 1
numtrials = 100

us = np.arange(N, dtype=np.float64) / N # normalized time coordinate

rztime = 0.0
rzwtime = 0.0
rzerrs = np.zeros((numtrials, 3))
rzwerrs = np.zeros((numtrials, 4))
theo_max_pow = N**2.0/4.0

for n in range(numtrials): 
    r = N/(4*numharm) + uniform(0.0, 1.0, 1)[0] # average freq over "observation"
    z = uniform(-100, 100, 1)[0] # average fourier f-dot
    w = uniform(-600, 600, 1)[0] # fourier freq double deriv
    w = 0.0 # fourier freq double deriv
    data = np.zeros_like(us)
    for ii in range(numharm):
        rh = r * (ii + 1)
        zh = z * (ii + 1)
        wh = w * (ii + 1)
        r0 = rh - 0.5 * zh + wh / 12.0 # Make symmetric for all z and w
        z0 = zh - 0.5 * wh
        phss = 2.0 * np.pi * (us * (us * (us * wh/6.0 + z0/2.0) + r0))
        data += np.cos(phss)
    data += noiseamp * norm(N)
    ft = presto.rfft(data)

    offset = uniform(-1.0, 1.0, 3) * np.array([0.5, 2.0, 20.0]) / (0.5 * numharm)

    a = time.clock()
    if (numharm > 1):
        [maxpow, rmax, zmax, rds] = presto.maximize_rz_harmonics(ft, r+offset[0],
                                                                 z+offset[1], numharm,
                                                                 norm=1.0)
    else:
        [maxpow, rmax, zmax, rd] = presto.maximize_rz(ft, r+offset[0],
                                                      z+offset[1],
                                                      norm=1.0)
    rztime += time.clock() - a
    rzerrs[n] = (maxpow/numharm - theo_max_pow) / theo_max_pow, rmax - r, zmax - z

    a = time.clock()
    if (numharm > 1):
        [maxpow, rmax, zmax, wmax, rds] = presto.maximize_rzw_harmonics(ft, r+offset[0],
                                                                        z+offset[1],
                                                                        w+offset[2], numharm,
                                                                        norm=1.0)
    else:
        [maxpow, rmax, zmax, wmax, rd] = presto.maximize_rzw(ft, r+offset[0],
                                                             z+offset[1],
                                                             w+offset[2],
                                                             norm=1.0)
    rzwtime += time.clock() - a
    rzwerrs[n] = (maxpow/numharm - theo_max_pow) / theo_max_pow, rmax - r, zmax - z, wmax - w

print("Time for  rz: %s".format(rztime / numtrials))
print("Time for rzw: %s".format(rzwtime / numtrials))

print("rzerrs:")
print("  avg: %6.3f %6.3f %6.3f" % tuple(rzerrs.mean(axis=0)))
print("  std: %6.3f %6.3f %6.3f" % tuple(rzerrs.std(axis=0)))

print("rzwerrs:")
print("  avg: %6.3f %6.3f %6.3f %6.3f" % tuple(rzwerrs.mean(axis=0)))
print("  std: %6.3f %6.3f %6.3f %6.3f" % tuple(rzwerrs.std(axis=0)))
