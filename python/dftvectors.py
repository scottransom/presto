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

def rebin_components(array, factor):
   array.shape = (len(array) / factor, factor)
   return umath.add.reduce(Numeric.transpose(array))

def add_components(components):
   """
   add_components(components):
        Add (accumulate) the components into a vector.
   """
   return umath.add.accumulate(components)

def vector_amplitude(vector):
   return umath.absolute(vector[-1])

def vector_phase(vector):
   return umath.arctan2(vector[-1].imag, vector[-1].real)

def rotate_vector(vector, angles):
   return vector * umath.exp(complex(0.0, 1.0) * angles)   

def point_distances(vec1, vec2):
   return umath.sqrt((vec2.real-vec1.real)**2.0 +
                     (vec2.imag-vec1.imag)**2.0)

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
      self.vector = add_components(self.components)
   def rebin(self, factor):
      """
      dftvector.rebin(factor):
           Rebin dftvector.components and dftvector.vector by 'factor'.
           For example, if 'factor' = 2, then the resulting vector will
           have half the number of bins as the current vector.
           If 'factor' = 0, reset the vector and its components to the
           original values.
      """
      if (factor == 0):
         self.n = self.orig_n
         self.numvect = self.orig_numvect
         self.components = copy.copy(self.orig_components)
         self.vector = copy.copy(self.orig_vector)
         self.timefract = copy.copy(self.orig_timefract)
         self.times = copy.copy(self.orig_times)
      elif not (factor * (self.numvect / factor) == self.numvect):
         print "'factor' in dftvector.rebin(factor) must evenly divide"
         print "   the number of vector components (numvect)."
      else:
         try:
            self.orig_n
         except AttributeError:
            self.orig_n = self.n
            self.orig_numvect = self.numvect
            self.orig_components = copy.copy(self.components)
            self.orig_vector = copy.copy(self.vector)
            self.orig_timefract = copy.copy(self.timefract)
            self.orig_times = copy.copy(self.times)
         self.n = self.n * factor
         self.numvect = self.numvect / factor
         self.components = rebin_components(self.components, factor)
         self.vector = add_components(self.components)
         self.timefract = (Numeric.arange(self.numvect) + 1.0) / self.numvect
         self.times = self.timefract * self.dt * self.n * self.numvect
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
      return vector_amplitude(self.vector)
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
              (przw[1] + u * (0.5 * przw[2] + przw[3] / 6.0 * u))) + przw[0]
   def fitrzw(self, przw=None):
      """
      dftvector.fitrzw():
           Fit a model consisting of an arbitrary phase, Fourier
           frequency offset (dr), Fourier f-dot (z), and Fourier
           f-dotdot (w) to dftvector.  Return the best fitting
           parameters in an array of [phase, dr, z, w].
      """
      def funct(przw, self):
         x = Numeric.zeros(4, 'd')
         x[0:3] = przw
         corrvect = umath.add.accumulate(self.rotate(-self.rzw_phases(x)))
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
   def fitrzw_dist(self, aprzw=None):
      """
      dftvector.fitrzw_dist():
           Fit a model consisting of an amplitude, phase, Fourier
           frequency offset (dr), Fourier f-dot (z), and Fourier
           f-dotdot (w) to dftvector.  Return the best fitting
           parameters in an array of [amp, phase, dr, z, w].
      """
      def funct(aprzw, self):
         amplitudes = Numeric.arange(self.numvect, typecode='d') / \
                      self.numvect  * aprzw[0]
         phases = self.rzw_phases(aprzw[1:])
         corrvect = add_components(rotate_vector(amplitudes, phases))
         #self.plotvector(corrvect)
         #return weights * point_distances(self.vector, corrvect)
         #weights = self.timefract**2.0
         return point_distances(self.vector, corrvect)
      if not aprzw:
         aprzw = [0.5 * self.maxamplitude(), 0.0, 0.0, 0.0, 0.0]
      retval = leastsq(funct, aprzw, args=(self))
      print retval[1]
      return retval[0]
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
         amp = 1.1 * umath.maximum(max(umath.absolute(x)),
                                   max(umath.absolute(y)))
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

dv = dftvector("testz_1000000.000.dftvec")
print dv.maxamplitude()
print dv.amplitude()
aprzw = [dv.amplitude(), 224.657*umath.pi/180.0, 0.345, 5.789, 0.0]
print "aprzw =", aprzw
phasecorr = -dv.rzw_phases(aprzw[1:])
corrvect = umath.add.accumulate(dv.rotate(phasecorr))
dv.plotvector(corrvect, color="red")
## amplitudes = Numeric.arange(dv.numvect, typecode='d') / \
##              dv.numvect * aprzw[0]
## phases = dv.rzw_phases(aprzw[1:])
## corrvect = add_components(rotate_vector(amplitudes, phases))
## dv.plotvector(corrvect, color="yellow")

aprzw = dv.fitrzw_dist()
print "aprzw =", aprzw
phasecorr = -dv.rzw_phases(aprzw[1:])
corrvect = umath.add.accumulate(dv.rotate(phasecorr))
dv.plotvector(corrvect, color="blue")
dv.plot()
Pgplot.closeplot()

## dv = dftvector("/raid/data/NGC6544/NGC6544_DM134_bary_566514.000.dftvec")
## przw = [0.0, -6.37, 8.6, 0.0]
## print "przw =", przw
## corrvect = add_components(dv.rotate(-dv.rzw_phases(przw)))
## print "corrvect amp =", vector_amplitude(corrvect)
## print "corrvect phs =", vector_phase(corrvect)
## dv.plotvector(corrvect)
## dv.plot()
## Pgplot.closeplot()

#dv = dftvector("/pulsar/data/parkes/Ter5_1_Part/Ter5_DM243_76375.000.dftvec")
#print dv.maxamplitude()
#print dv.amplitude()
#dv.plot(amp=12)
#przw = dv.fitrzw()
#rot = dv.rzw_phases(przw)
#corrvect = dv.rotate(-rot)
#dv.plotvector(umath.add.accumulate(corrvect))
#Pgplot.closeplot()
