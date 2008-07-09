from numpy import *
from presto import *
from ppgplot import *
from Pgplot import *

N = 2**20
ro = 2.0**18
dr = 0.03125
dz = 0.2
np = 384

data = cos(arange(N)*2*pi*ro/N)
ft = rfft(data)
ffdot = ffdot_plane(ft, ro-np/2*dr, dr, np, 0.0-np/2*dz, dz, np)
pffdot = spectralpower(ffdot.flat)
print "Max pow = %f  Min pow = %f" % (max(pffdot), min(pffdot))
pffdot = pffdot / max(pffdot)
pffdot.shape = (np,np)
r = arange(np) * dr - np/2*dr
z = arange(np) * dz - np/2*dz
print len(r), len(z)
rgx = asarray([r[0],r[np-1]])
rgy = asarray([z[0],z[np-1]])
freqcut = pffdot[np/2, :]
fdotcut = pffdot[:, np/2]

image='antirainbow'
device='ffdot_combined.eps/VCPS'
device='/XWIN'
labx='Fourier Frequency Offset (bins)'
laby='Fourier Frequency Derivative (bins)'
contours=asarray([0.1, 0.3, 0.5, 0.7, 0.9])

imfract = 0.65
margin = 0.08

pgopen(device)
pgpap(0.0, 1.0)
pgpage()

# freq cut
pgsvp(margin, margin+imfract, margin+imfract, 1.0-margin/2)
pgswin(min(r), max(r), -0.1, 1.1)
pgbox("BCST", 0.0, 0, "BCNST", 0.0, 0)
pgline(r, freqcut)
pgmtxt("L", 2.0, 0.5, 0.5, "Relative Power");

#fdot cut
pgsvp(margin+imfract, 1.0-margin/2, margin, margin+imfract)
pgswin(-0.1, 1.1, min(z), max(z))
pgbox("BCNST", 0.0, 0, "BCST", 0.0, 0)
pgline(fdotcut, z)
pgmtxt("B", 2.4, 0.5, 0.5, "Relative Power");

# f-fdot image
pgsvp(margin, margin+imfract, margin, margin+imfract)
pgswin(min(r), max(r), min(z), max(z))
pgmtxt("B", 2.4, 0.5, 0.5, labx);
pgmtxt("L", 2.0, 0.5, 0.5, laby);
lo_col_ind, hi_col_ind = pgqcol()
lo_col_ind = lo_col_ind + 2
pgscir(lo_col_ind, hi_col_ind)
pgpalette.setpalette(image)
pgctab(pgpalette.l,pgpalette.r,pgpalette.g,pgpalette.b)
pgimag_s(pffdot, 0.0, 0.0, rgx[0], rgy[0], rgx[1], rgy[1])  
pgsci(1)
pgcont_s(pffdot, len(contours), contours, rgx[0], rgy[0], rgx[1], rgy[1])  
pgbox("BCST", 0.0, 0, "BCST", 0.0, 0)
pgsci(1)
pgbox("N", 0.0, 0, "N", 0.0, 0)

# gray axes
pgscr(1, 0.5, 0.5, 0.5)
pgsci(1)
pgslw(2)
pgline(rgx, asarray([0.0, 0.0]))
pgline(asarray([0.0, 0.0]), rgy)

closeplot()


