import umath
import Numeric
import Pgplot
import copy
from Multipack import leastsq

def smooth_phases(phases):
   for ii in xrange(len(phases)-1):
      l = phases[ii]
      h = phases[ii+1]
      if (h - l < -umath.pi):
         phases[ii+1:] = phases[ii+1:] + 2.0 * umath.pi
      if (h - l > umath.pi):
         phases[ii+1:] = phases[ii+1:] - 2.0 * umath.pi
   return phases

class dftvector:
   def __init__(self, filename=None):
      if filename:
         self.fromfile(filename)
   def fromfile(self, filename):
      """
      dftvector.fromfile(filename):
           Read a dftvector from 'filename'.
      """
      self.filename = filename
      infile = open(filename, "rb")
      dblarr = Numeric.fromstring(infile.read(), 'd')
      infile.close()
      self.n = int(dblarr[0])
      self.numvect = int(dblarr[1])
      self.dt = dblarr[2]
      self.r = dblarr[3]
      self.norm = dblarr[4]
      self.components = Numeric.zeros(self.numvect, 'D')
      self.timefract = (Numeric.arange(self.numvect) + 1.0) / self.numvect
      self.times = self.timefract * self.dt * self.n * self.numvect
      for ii in xrange(self.numvect):
         index = 2 * ii + 5
         self.components[ii] = complex(dblarr[index], dblarr[index+1])
      self.vector = self.addcomponents()
   def addcomponents(self):
      """
      dftvector.addcomponents():
           Add the dftvector.components into a vector.
      """
      return umath.add.accumulate(self.components)
   def rotate(self, angles):
      """
      dftvector.rotate(angles):
           Rotate dftvector.components by the angles in the array 
           'angles' (radians) and return the newly rotated vector.
      """
      return self.components * umath.exp(complex(0.0, 1.0) * angles)
   def maxamplitude(self):
      """
      dftvector.maxamplitude():
           Return the maximum possible amplitude for a dftvector.
      """
      return umath.add.reduce(umath.absolute(self.components))
   def amplitude(self):
      """
      dftvector.amplitude():
           Return the current amplitude for a dftvector.
      """
      return umath.absolute(self.vector[-1])
   def phases(self):
      """
      dftvector.phases():
           Return the accumulated phase (radians) of dftvector.vector.
           The phases are corrected to contain no discontinuities.
           The phase of the first vector component (0.0, 0.0) is
                _not_ returned.
      """
      phases = umath.arctan2(self.vector.imag, self.vector.real)
      return smooth_phases(phases)
   def rzw_phases(self, przw):
      """
      dftvector.rzw_phases(self, przw):
           Return an accumulated phase model that includes the
           initial Fourier phase in radians (p), the Fourier
           frequency offset (r), Fourier f-dot (z), and Fourier
           f-dotdot (w).  'przw' is a list with 4 elements.
      """
      u = self.timefract
      return (2.0 * umath.pi * u *
              (przw[1] + u * (0.5 * przw[2] + przw[3] / 6.0 * u)) + przw[0])
   def fitrzw(self, przw=None):
      """
      dftvector.fitrzw():
           Fit a model consisting of and arbitrary phase, Fourier
           frequency offset (dr), Fourier f-dot (z), and Fourier
           f-dotdot (w) to dftvector.  Return the best fitting
           parameters in a tuple of (phase, dr, z, w).  The mod
      """
      def funct(przw, self):
         x = Numeric.zeros(4, 'd')
         x[0:3] = przw
         corrvect = umath.add.accumulate(dv.rotate(-dv.rzw_phases(x)))
         phases = umath.arctan2(corrvect.imag, corrvect.real)
         #self.plotvector(corrvect)
         #Pgplot.closeplot()
         return smooth_phases(phases)
         #weights = self.timefract
         #Pgplot.plotxy(weights * (self.phases() - self.rzw_phases(x)))
         #Pgplot.closeplot()
         #return weights * (self.phases() - self.rzw_phases(x))
      if not przw:
         przw = [0.0, 0.0, 0.0]
      retval = leastsq(funct, przw, args=(self))
      print retval[1]
      x = Numeric.zeros(4, 'd')
      x[0:3] = retval[0]
      return x
   def plotvector(self, vector, amp=None, color='red', device='/XWIN'):
      """
      dftvector.plotvector(vector, amp=None, color='red', device='/XWIN'):
           Plot 'vector' in the complex plane.
      """
      x = Numeric.zeros(self.numvect+1, 'd')
      y = Numeric.zeros(self.numvect+1, 'd')
      x[1:] = vector.real
      y[1:] = vector.imag
      if not amp:
         amp = 1.1 * umath.absolute(complex(x[-1], y[-1]))
      Pgplot.plotxy(y, x, rangex=[-amp, amp], rangey=[-amp, amp],
                    labx='Real Amplitude', laby='Imaginary Amplitude',
                    color=color, aspect=1.0)
   def plot(self, amp=None, color='white', device='/XWIN'):
      """
      dftvector.plot(amp=None, color='white', device='/XWIN'):
           Plot the dftvector summation.
      """
      self.plotvector(self.vector, amp=amp, color=color, device=device)

#dv = dftvector("/raid/data/Ter5_98/Day1-2/Ter5_1-2_bary_DM236_3564593.380.dftvec")
#print dv.maxamplitude()
#print dv.amplitude()
#dv.plot()

dv = dftvector("testz_20000.000.dftvec")
print dv.maxamplitude()
print dv.amplitude()
przw = dv.fitrzw([2.0, 2.0, 6.0])
print "przw =", przw
rot = dv.rzw_phases(przw)
#rot = dv.rzw_phases([0.0, 0.345, 0.0, 0.0])
corrvect = dv.rotate(-rot)
dv.plot()
dv.plotvector(umath.add.accumulate(corrvect))
Pgplot.closeplot()


#dv = dftvector("/pulsar/data/parkes/Ter5_1_Part/Ter5_DM243_76375.000.dftvec")
#print dv.maxamplitude()
#print dv.amplitude()
#dv.plot(amp=12)
#przw = dv.fitrzw()
#rot = dv.rzw_phases(przw)
#corrvect = dv.rotate(-rot)
#dv.plotvector(umath.add.accumulate(corrvect))
#Pgplot.closeplot()
