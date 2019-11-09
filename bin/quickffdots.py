#!/usr/bin/env python
from __future__ import print_function
from builtins import range
from presto.infodata import *
from presto.presto import ffdot_plane, spectralpower
from pylab import *
import numpy as N
import sys

numharm = 4

# Contour values as a fraction of max in "window"
rel_convals = N.asarray([0.5, 0.8, 0.95, 0.995])
rel_alphas  = N.asarray([0.3, 0.5,  0.7,   1.0])
# Absolute power values for contours
abs_convals = N.asarray([5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 1e6])
abs_alphas  = N.asarray([0.1,  0.3,  0.4,  0.5, 0.65,   0.8, 1.0])

if len(sys.argv) < 3:
    print("usage:  quickffdots.py fftfile freq(Hz)")
    sys.exit(0)
    
if 0:
    convals = rel_convals
    alphas  = rel_alphas
else:
    convals = abs_convals
    alphas  = abs_alphas
    
if (numharm > 6):
    print("Numharm must be < 6!")
    sys.exit(0)

# Each of the harmonics must have an even number of bins interpolated
# (i.e. integer number of interpolated pointed between the integer bins)
drs = {1: 1.0/4, 2: 1.0/8, 3: 1.0/12, 4: 1.0/12, 5: 1.0/60, 6: 1.0/60}

# The following are for the lowest harmonic
zmax = 20
dr = drs[numharm]
dz = dr * 4.0
numrs = int(round(4*zmax/dr))
numzs = int(round(2*zmax/dz)) + 1

infilenm = sys.argv[1]
infile = open(infilenm, 'rb')
idata = infodata(infilenm[:-4]+".inf")
idata.T = idata.N * idata.dt

ctrr = int(round(float(sys.argv[2]) * idata.T))
startr = int(ctrr - numrs/2 * dr)
ctrfreq = ctrr / idata.T

ffdps = []
maxvals = []
maxargs = []

for harmnum in range(1, numharm+1):
    print("Computing harmonic", harmnum)
    ldr = dr * harmnum
    ldz = dz * harmnum
    lor = startr * harmnum
    loz = 0.0 - (numzs-1)/2 * ldz
    rs = N.arange(numrs) * ldr + lor
    zs = N.arange(numzs) * ldz + loz
    if harmnum==1:
        rs0, zs0 = rs[:], zs[:]
    lo_file_r = int(rs.min()) - 1000
    hi_file_r = int(rs.max()) + 1000

    # Read and normalize the raw spectrum
    infile.seek(lo_file_r * 8, 0)
    fftamps = N.fromfile(infile, 'F', hi_file_r-lo_file_r+1)
    fftpows = spectralpower(fftamps)
    pownorm = 1.0 / (1.442695 * N.median(fftpows))
    fftamps *= sqrt(pownorm)

    ffd = ffdot_plane(fftamps, lor-lo_file_r, ldr, numrs, loz, ldz, numzs)
    ffd_pows = (ffd * ffd.conj()).real
    ffdps.append(ffd_pows)
    if (harmnum==1):
        sumpows = N.zeros_like(ffd_pows)
    sumpows += ffd_pows

    argmax = ffd_pows.argmax()
    maxvals.append(ffd_pows.max())
    maxargs.append((argmax // numrs, argmax % numrs))
    print("  Maximum power for harmonic %d = %.2f"%(harmnum, maxvals[-1]))

    if (convals.max() < 1.5): # Using relative contours
        print("Using relative contours..")
        pow_contours = convals * maxvals[-1]
    else:
        # Only choose the contours with values < the max power
        highcut = N.compress(abs_convals > maxvals[-1],
                             N.arange(len(abs_convals))).min()
        pow_contours = N.empty(highcut+1, dtype=float)
        pow_contours[:highcut] = abs_convals[:highcut]
        pow_contours[highcut] = 1e6
        alphas = abs_alphas[:highcut+1]

    if harmnum==1:   # 'Red'
        colorvals = [(alpha, 0.0, 0.0) for alpha in alphas]
    elif harmnum==2: # 'Green'
        colorvals = [(0.0, alpha, 0.0) for alpha in alphas]
    elif harmnum==3: # 'Blue'
        colorvals = [(0.0, 0.0, alpha) for alpha in alphas]
    elif harmnum==4: # 'Cyan'
        colorvals = [(0.0, alpha, alpha) for alpha in alphas]
    elif harmnum==5: # 'Magenta'
        colorvals = [(alpha, 0.0, alpha) for alpha in alphas]
    elif harmnum==6: # 'Yellow'
        colorvals = [(alpha, alpha, 0.0) for alpha in alphas]
    limits = [rs0.min(), rs0.max(), zs0.min(), zs0.max()]

    cstr = "".join(["%.2f "%x for x in pow_contours[:-1]])
    print("  Contour levels at powers = "+cstr)
    
    # Plot the contours
    contour(ffd_pows, pow_contours, origin='lower',
            alpha=0.3, colors=colorvals, extent=limits)
    # Plot fill between the last two contours
    contourf(ffd_pows, [pow_contours[-2], pow_contours[-1]],
             origin='lower', colors=[colorvals[-1]], alpha=0.3,
             extent=limits)
    xlabel("Average Fourier Frequency (bins)")
    ylabel("Fourier Frequency Derivative z (bins)")
    if harmnum==1:
        axhline(0.0, linewidth=1, color='black', alpha=0.3)

print("\nMax summed power = %.2f"%(sumpows.max()))
argmax = sumpows.argmax()
maxr = rs0[argmax % numrs]
maxz = zs0[argmax // numrs]
maxf = maxr/idata.T
maxfd = maxz/(idata.T*idata.T)
initf = (maxr - 0.5*maxz)/idata.T
print("  at r =", maxr, " (%.10f Hz)"%maxf)
print("  at z =", maxz, " (%.10g Hz/s)"%maxfd)
print("Folding command would be: ")
print("  prepfold -f %.10f -fd %.6g ..." %(initf, maxfd))

infile.close()
show()



