from numpy import *
from Pgplot import *
from presto import *

kern_half_width = 10
numbetween = 10
numkern = 2 * numbetween * kern_half_width
f = arange(numkern, dtype=float64) / numbetween - kern_half_width
kern = gen_z_response(0.0, numbetween, 0.0, numkern)
pkern = spectralpower(kern)
print "Freq ", f[len(f)/2], " = ", pkern[len(f)/2]
plotxy(pkern, f, labx="Fourier Frequency Offset", \
       laby="Normalized Power", device="z_responses.eps/CPS")
#plotxy(pkern, f, labx="Fourier Frequency Offset", \
#       laby="Normalized Power")
for i in range(4):
    z = 5.0 * i
    kern = gen_z_response(0.0, numbetween, z, numkern)
    pkern = spectralpower(kern)
    plotxy(pkern, f, color=i+1)
    ppgplot.pgtext(5.0, 0.8 - i*0.1, 'z = %1.0f' % z)
closeplot()
