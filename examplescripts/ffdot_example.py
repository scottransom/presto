from __future__ import print_function
import numpy as num
from presto import presto
import presto.ppgplot as ppgplot
from presto.Pgplot import pgpalette
from numpy.random import standard_normal as norm
import time

N = 2**14
r = N/4.0 # average freq over "observation"
#r = N/4.0 + 0.5 # average freq over "observation"
rint = num.floor(r)
dr = 1.0/32.0
dz = 0.18
np = 512 # number of pixels across for f-fdot image
z = 10.0 # average fourier f-dot
w = -40.0 # fourier freq double deriv
noise = 0.0
noise = 1.0*norm(N)

us = num.arange(N, dtype=num.float64) / N # normalized time coordinate
r0 = r - 0.5 * z + w / 12.0 # Make symmetric for all z and w
z0 = z - 0.5 * w
phss = 2.0 * num.pi * (us * (us * (us * w/6.0 + z0/2.0) + r0))
ft = presto.rfft(num.cos(phss)+noise)
ffdot = presto.ffdot_plane(ft, rint-np/2*dr, dr, np, 0.0-np/2*dz, dz, np)
pffdot = presto.spectralpower(ffdot.flat)
theo_max_pow = N**2.0/4.0
frp = max(pffdot) / theo_max_pow # Fraction of recovered power
print("Fraction of recovered signal power = %f" % frp)
a = time.clock()
[maxpow, rmax, zmax, rd] = presto.maximize_rz(ft, r+norm(1)[0]/5.0,
                                              z+norm(1)[0], norm=1.0)
print("Time for rz:", time.clock()-a)
print(r, rmax, z, zmax, theo_max_pow, maxpow)
a = time.clock()
[maxpow, rmax, zmax, wmax, rd] = presto.maximize_rzw(ft, r+norm(1)[0]/5.0,
                                                     z+norm(1)[0],
                                                     w+norm(1)[0]*5.0,
                                                     norm=1.0)
print("Time for rzw:", time.clock()-a)
print(r, rmax, z, zmax, w, wmax, theo_max_pow, maxpow)
#print "Raw power should be ~%.2e" % theo_max_pow
pffdot = pffdot / theo_max_pow
pffdot.shape = (np, np)
rs = num.arange(np) * dr - np//2*dr
zs = num.arange(np) * dz - np//2*dz
rgx = num.asarray([rs[0], rs[np-1]])
rgy = num.asarray([zs[0], zs[np-1]])
freqcut = pffdot[np//2, :]
fdotcut = pffdot[:, np//2]

image='antirainbow'
device='ffdot_combined.eps/VCPS'
device='/XWIN'
labx='Fourier Frequency Offset (bins)'
laby='Fourier Frequency Derivative (bins)'
contours = num.asarray([0.1, 0.3, 0.5, 0.7, 0.9])

imfract = 0.65
margin = 0.08

ppgplot.pgopen(device)
ppgplot.pgpap(0.0, 1.0)
ppgplot.pgpage()

# Give z and w values and power change
ppgplot.pgsvp(margin+imfract, 1.0-margin/2, margin+imfract, 1.0-margin/2)
ppgplot.pgswin(0.0, 1.0, 0.0, 1.0)
ppgplot.pgtext(0.1, 0.8, "Frac Recovered" % frp)
ppgplot.pgtext(0.2, 0.65, "Power = %.3f" % frp)
ppgplot.pgtext(0.1, 0.4, "signal z = %.1f" % z)
ppgplot.pgtext(0.1, 0.25, "signal w = %.1f" % w)

# freq cut
ppgplot.pgsvp(margin, margin+imfract, margin+imfract, 1.0-margin/2)
ppgplot.pgswin(min(rs), max(rs), -0.1, 1.1)
ppgplot.pgbox("BCST", 0.0, 0, "BCNST", 0.0, 0)
ppgplot.pgline(rs, freqcut)
ppgplot.pgmtxt("L", 2.0, 0.5, 0.5, "Relative Power");

#fdot cut
ppgplot.pgsvp(margin+imfract, 1.0-margin/2, margin, margin+imfract)
ppgplot.pgswin(-0.1, 1.1, min(zs), max(zs))
ppgplot.pgbox("BCNST", 0.0, 0, "BCST", 0.0, 0)
ppgplot.pgline(fdotcut, zs)
ppgplot.pgmtxt("B", 2.4, 0.5, 0.5, "Relative Power");

# f-fdot image
ppgplot.pgsvp(margin, margin+imfract, margin, margin+imfract)
ppgplot.pgswin(min(rs), max(rs), min(zs), max(zs))
ppgplot.pgmtxt("B", 2.4, 0.5, 0.5, labx);
ppgplot.pgmtxt("L", 2.0, 0.5, 0.5, laby);
lo_col_ind, hi_col_ind = ppgplot.pgqcol()
lo_col_ind = lo_col_ind + 2
ppgplot.pgscir(lo_col_ind, hi_col_ind)
pgpalette.setpalette(image)
ppgplot.pgctab(pgpalette.l, pgpalette.r, pgpalette.g, pgpalette.b)
ppgplot.pgimag_s(pffdot, 0.0, 0.0, rgx[0], rgy[0], rgx[1], rgy[1])
ppgplot.pgsci(1)
ppgplot.pgcont_s(pffdot, len(contours), contours, rgx[0], rgy[0], rgx[1], rgy[1])
ppgplot.pgbox("BCST", 0.0, 0, "BCST", 0.0, 0)
ppgplot.pgsci(1)
ppgplot.pgbox("N", 0.0, 0, "N", 0.0, 0)

# gray axes
ppgplot.pgscr(1, 0.5, 0.5, 0.5)
ppgplot.pgsci(1)
ppgplot.pgslw(2)
ppgplot.pgline(rgx, num.asarray([0.0, 0.0]))
ppgplot.pgline(num.asarray([0.0, 0.0]), rgy)

ppgplot.pgclos()
