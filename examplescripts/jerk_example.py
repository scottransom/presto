from __future__ import print_function
import numpy as num
from presto import presto
import presto.ppgplot as ppgplot
import time
from presto.Pgplot import pgpalette

N = 2**14
r = N/4.0 # average freq over "observation"
#r = N/4.0 + 0.5 # average freq over "observation"
rint = num.floor(r)
numbetween = 8
dr = 1.0/numbetween
dz = 4.0/numbetween
dw = 20.0/numbetween
np = 256 # number of pixels across for f-fdot image
z = 0.0 # average fourier f-dot
w = 0.0 # fourier freq double deriv
#noise = 0.0
noise = 0.0*num.random.standard_normal(N)

us = num.arange(N, dtype=num.float64) / N # normalized time coordinate
r0 = r - 0.5 * z + w / 12.0 # Make symmetric for all z and w
z0 = z - 0.5 * w
phss = 2.0 * num.pi * (us * (us * (us * w/6.0 + z0/2.0) + r0))
ft = presto.rfft(num.cos(phss)+noise)

a = time.clock()
vol = presto.fdotdot_vol(ft, rint-np/2*dr, dr, np,
                         0.0-np/2*dz, dz, np,
                         0.0-np/2*dw, dw, np)
print("First jerk vol took %.3f s" % (time.clock()-a))
a = time.clock()
vol = presto.fdotdot_vol(ft, rint-np/2*dr, dr, np,
                         0.0-np/2*dz, dz, np,
                         0.0-np/2*dw, dw, np)
print("Second jerk vol took %.3f s" % (time.clock()-a))
pvol = presto.spectralpower(vol.flat)
theo_max_pow = N**2.0/4.0
frp = max(pvol) / theo_max_pow # Fraction of recovered power
print("Fraction of recovered signal power = %f" % frp)
[maxpow, rmax, zmax, rd] = presto.maximize_rz(ft, r+num.random.standard_normal(1)[0]/5.0,
                                       z+num.random.standard_normal(1)[0], norm=1.0)
print(r, rmax, z, zmax, theo_max_pow, maxpow)
# print("Raw power should be ~%.2e" % theo_max_pow)
pvol = pvol / theo_max_pow
pvol.shape = (np, np, np)
rs = num.arange(np) * dr - np/2*dr
zs = num.arange(np) * dz - np/2*dz
ws = num.arange(np) * dw - np/2*dw
rgx = num.asarray([rs[0], rs[np-1]])
rgy = num.asarray([zs[0], zs[np-1]])

# Use the following if you want frames for a movie.  See the bottom
# of this file for the other commands to generate that movie.
#device='jerk_%03d.eps/VCPS'
device='/XWIN'

image='antirainbow'
labx='Fourier Frequency Offset (bins)'
laby='Fourier Frequency Derivative (bins)'
contours = num.asarray([0.1, 0.3, 0.5, 0.7, 0.9])

imfract = 0.65
margin = 0.08

if device=="/XWIN":
    ppgplot.pgopen(device)
    ppgplot.pgpap(0.0, 1.0)
    ppgplot.pgpage()
    
for ii in range(np):
    if not device=="/XWIN":
        ppgplot.pgopen(device%ii)
        ppgplot.pgpap(0.0, 1.0)
        ppgplot.pgpage()
    freqcut = pvol[ii, np//2, :]
    fdotcut = pvol[ii, :, np//2]
    frp = pvol[ii].max() # Fraction of recovered power
    print("w = %.3f  frac pow recovered = %.3f" % (ws[ii], frp))
    # Give z and w values and power change
    ppgplot.pgsvp(margin+imfract, 1.0-margin/2, margin+imfract, 1.0-margin/2)
    ppgplot.pgswin(0.0, 1.0, 0.0, 1.0)
    ppgplot.pgtext(0.1, 0.8, "Frac Recovered")
    ppgplot.pgtext(0.2, 0.65, "Power = %.3f" % frp)
    ppgplot.pgtext(0.1, 0.4, "signal z = %.1f" % z)
    ppgplot.pgtext(0.1, 0.25, "w = %.1f" % ws[ii])

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
    ppgplot.pgimag_s(pvol[ii], 0.0, 0.0, rgx[0], rgy[0], rgx[1], rgy[1])  
    ppgplot.pgsci(1)
    ppgplot.pgcont_s(pvol[ii], len(contours), contours, rgx[0], rgy[0], rgx[1], rgy[1])  
    ppgplot.pgbox("BCST", 0.0, 0, "BCST", 0.0, 0)
    ppgplot.pgsci(1)
    ppgplot.pgbox("N", 0.0, 0, "N", 0.0, 0)

    # gray axes
    ppgplot.pgscr(1, 0.5, 0.5, 0.5)
    ppgplot.pgsci(1)
    ppgplot.pgslw(2)
    ppgplot.pgline(rgx, num.asarray([0.0, 0.0]))
    ppgplot.pgline(num.asarray([0.0, 0.0]), rgy)

    time.sleep(0.1)
    if device=="/XWIN":
        ppgplot.pgeras()
    else:
        ppgplot.pgclos()

if device=="/XWIN":
    ppgplot.pgclos()
else:
    print("""If you want to make a movie with the resulting .eps files, here are
the appropriate commands:
    
> python jerk_example.py 
> pstoimg -density 200 -antialias -crop a jerk_*eps
> ffmpeg -r 16 -f image2 -s 1000x1000 -i jerk_%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p jerk_search.mp4
""")
