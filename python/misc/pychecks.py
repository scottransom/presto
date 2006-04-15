from numpy import *
from Pgplot import *
from psr_utils import *
from presto import *

# Check the image plotting software
a = distance(201)
b = exp(-a**2/1000.0)
x = arange(-100.0, 101.0, 1.0)
plot2d(b, x, x, contours=array([0.2, 0.5, 0.7]), color='black', labels='white', labelint=50, labelmin=30)
closeplot()

# Check the f-fdot plane generator
fftfile = open('../../tests/testz.fft')
m = 0
numbetween = 8
lof = 19990.0
loz = -20.0
hiz = 20.0
numz = 81.0
numfplot = 20.0
dz = (hiz-loz)/(numz-1.0)
a = corr_rz_plane_file(fftfile, numbetween, lof, loz, hiz, numz, 1024, LOWACC, m) 
a[1]
b = (a[0].real**2 + a[0].imag**2)/get_numphotons(fftfile)
x = arange(numfplot*numbetween, dtype='f') / numbetween + lof
y = arange(numz, dtype='f') * dz + loz
c = arange(6) * 10.0 + 10.0
plot2d(b[:,0:int(numfplot*numbetween)], x, y, \
       labx="Fourier Frequency", laby="Fourier F-dot", \
       contours=c, color='black', image='heat') #, \
#       device='test.ps/CPS')
nextplotpage(1)
plot2d(b[:,0:int(numfplot*numbetween)], x, y, \
       labx="Fourier Frequency", laby="Fourier F-dot", \
       contours=c, color='black', image='astro') #, \
#       device='test.ps/CPS')
closeplot()

