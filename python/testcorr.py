from __future__ import print_function
from builtins import range
from numpy import *
from Pgplot import *
from presto import *

file = "testz"
displaybins = 20
numbetween = 16
numpts = displaybins * numbetween

# Read the '.mak' file
md = read_makfile(file)

# Open and read the data
fftfile = open(file+".fft","r")
filelen = chkfilelen(fftfile, 8)
data = fftfile.read()
data = fromstring(data, "F")
fftfile.close()
nph = data[0].real

# Set up some useful things
centerr = md.r + md.z / 2.0
startbin = floor(centerr - displaybins / 2.0)
ca = zeros(numpts, dtype=complex64)
cf = arange(startbin, startbin + displaybins, 1.0 / numbetween)

# Show the power spectrum without f-dot correction
kern_halfwidth = z_resp_halfwidth(0.0, HIGHACC)
numkern = 2 * numbetween * kern_halfwidth
for i in range(numpts):
    ca[i] = rz_interp(data, filelen, cf[i], 0.0, kern_halfwidth)
cpow = spectralpower(asarray(ca)) / nph
cphs = spectralphase(ca)
maxval = argmax(cpow)
plotxy(cpow, cf-1e6, labx="Fourier Frequency - 1e6", laby="Power")
print("Maximum value is at r =", startbin + maxval / float(numbetween))
print("   Power =", cpow[maxval], "  Phase =", cphs[maxval])
closeplot()

# Show the power spectrum with f-dot correction
kern_halfwidth = z_resp_halfwidth(md.z, HIGHACC)
numkern = 2 * numbetween * kern_halfwidth
for i in range(numpts):
    ca[i] = rz_interp(data, filelen, cf[i], md.z, kern_halfwidth)
cpow = spectralpower(ca) / nph
cphs = spectralphase(ca)
maxval = argmax(cpow)
plotxy(cpow, cf-1e6, labx="Fourier Frequency - 1e6", laby="Power")
print("Maximum value is at r =", startbin + maxval / float(numbetween))
print("   Power =", cpow[maxval], "  Phase =", cphs[maxval])
closeplot()

