import umath
import Numeric
from Pgplot import plotxy, closeplot
from copy import copy
from time import sleep
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

def p_to_f(p, pd, pdd=None):
   """
   p_to_f(p, pd, pdd=None):
      Convert period, period derivative and period second
      derivative to the equivalent frequency counterparts.
      Will also convert from f to p.
   """
   f = 1.0 / p
   fd = -pd / (p * p)
   if (pdd==None):
       return [f, fd]
   else:
       if (pdd==0.0):
           fdd = 0.0
       else:
           fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
       return [f, fd, fdd]

def rzw_dist_model(przw, self):
   przw = przw * Numeric.array([1.0, 1.0, 4.0, 20.0])
   amplitudes = Numeric.ones(self.numvect, 'D')
   phases = self.rzw_phases(przw)
   corrvect = add_components(rotate_components(amplitudes, phases))
   norm = self.amplitude() / vector_amplitude(corrvect)
   weights = Numeric.arange(1.0, 0.0, -1.0/self.numvect)
   return weights * point_distances(self.vector, norm * corrvect)

def rzw_phase_model(przw, self):
   przw = przw * Numeric.array([1.0, 1.0, 4.0, 20.0])
   corrvect = add_components(self.rotate(-self.rzw_phases(przw)))
   phases = umath.arctan2(corrvect.imag, corrvect.real)
   weights = self.timefract
   return weights * smooth_phases(phases)

def vector_amplitude(vector):
   return umath.absolute(vector[-1])

def vector_phase(vector):
   return umath.arctan2(vector[-1].imag, vector[-1].real)

def rotate_components(components, angles):
   return components * umath.exp(complex(0.0, 1.0) * angles)   

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
      self.T = dblarr[5]
      self.components = Numeric.zeros(self.numvect, 'D')
      self.timefract = (Numeric.arange(self.numvect) + 1.0) / self.numvect
      self.times = self.timefract * self.dt * self.n * self.numvect
      for ii in xrange(self.numvect):
         index = 2 * ii + 6
         self.components[ii] = complex(dblarr[index], dblarr[index+1])
      self.vector = add_components(self.components)
      self.mask = Numeric.ones(self.numvect)
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
         self.components = copy(self.orig_components)
         self.vector = copy(self.orig_vector)
         self.timefract = copy(self.orig_timefract)
         self.times = copy(self.orig_times)
      elif not (factor * (self.numvect / factor) == self.numvect):
         print "'factor' in dftvector.rebin(factor) must evenly divide"
         print "   the number of vector components (numvect)."
      else:
         try:
            self.orig_n
         except AttributeError:
            self.orig_n = self.n
            self.orig_numvect = self.numvect
            self.orig_components = copy(self.components)
            self.orig_vector = copy(self.vector)
            self.orig_timefract = copy(self.timefract)
            self.orig_times = copy(self.times)
         if (factor > 1):
            self.n = self.n * factor
            self.numvect = self.numvect / factor
            self.components = rebin_components(self.components, factor)
            self.vector = add_components(self.components)
            self.timefract = (Numeric.arange(self.numvect) + 1.0) / self.numvect
            self.times = self.timefract * self.dt * self.n * self.numvect
   def mask(self, start, stop):
      """
      dftvector.mask(start, stop):
           Set the dftvector.components that are located between the
           fractional times (i.e. 0.0-1.0) 'start' and 'stop' to 0.0.
      """
      self.rebin(1)
      beg = Numeric.searchsorted(self.timefract, start)
      end = Numeric.searchsorted(self.timefract, stop)
      self.mask[beg:end] = 0
      self.components = self.components * self.mask
      self.vector = add_components(self.components)
   def unmask(self):
      """
      dftvector.unmask():
           Remove all masks from the data.
      """
      self.rebin(0)
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
      if (przw[3]==0.0):
         if (przw[2]==0.0):
            return 2.0 * umath.pi * u * przw[1] + przw[0]
         else:
            return 2.0 * umath.pi * u * (przw[1] + 0.5 * u * przw[2]) + przw[0]
      else:
         return (2.0 * umath.pi * u *
                 (przw[1] + 0.5 * u * (przw[2] + u * przw[3] / 3.0))) + przw[0]
   def rzwinfo(self, przw):
      """
      dftvector.rzwinfo(self, przw):
           Nicely print the conversions of the przw array.
      """
      modprzw = Numeric.asarray(copy(przw[1:]))
      modprzw[0] = (self.r + modprzw[0])
      [f, fd, fdd] = (modprzw / self.T**Numeric.array([1.0, 2.0, 3.0])).tolist()
      [p, pd, pdd] = p_to_f(f, fd, fdd)
      print ""
      print "      Phase (rad) = %.3f" % (przw[0])
      print "      Phase (deg) = %.3f" % (przw[0] * 180.0 / umath.pi)
      print "     Fourier freq = %.3f" % (modprzw[0])
      print "     Fourier fdot = %.5g" % (modprzw[1])
      print "  Fourier fdotdot = %.5g" % (modprzw[2])
      print "   Frequency (Hz) = %.9f" % (f)
      print "      fdot (Hz/s) = %.5g" % (fd)
      print " fdotdot (Hz/s^2) = %.5g" % (fdd)
      print "       Period (s) = %.12f" % (p)
      print "       pdot (s/s) = %.5g" % (pd)
      print "  pdotdot (s/s^2) = %.5g" % (pdd)
      print ""
   def fitrzw(self, przw=None, funct=rzw_dist_model):
      """
      dftvector.fitrzw(przw=None, funct=rzw_dist_model):):
           Fit a model consisting of an arbitrary phase, Fourier
           frequency offset (dr), Fourier f-dot (z), and Fourier
           f-dotdot (w) to dftvector.  Return the best fitting
           parameters in an array of [phase, dr, z, w].  If 'przw'
           is defined it is used as the starting guess for the
           fit.  If 'funct' is defined that is the function that
           is used during the minimization.
      """
      if przw:
         przw = Numeric.asarray(przw) / Numeric.array([1.0, 1.0, 4.0, 20.0])
      else:
         przw = Numeric.array([0.0, 0.0, 0.0, 0.0])
      retval = leastsq(funct, przw, args=(self))
      print retval[1]
      return retval[0] * Numeric.array([1.0, 1.0, 4.0, 20.0])
   def plotvector(self, vector, amp=None, color='red',
                  device='/XWIN', step=0):
      """
      dftvector.plotvector(vector, amp=None, color='red',
                           device='/XWIN', step=0)
           Plot 'vector' in the complex plane.
      """
      x = Numeric.zeros(self.numvect+1, 'd')
      y = Numeric.zeros(self.numvect+1, 'd')
      x[1:] = vector.real
      y[1:] = vector.imag
      if not amp:
         amp = 1.1 * umath.maximum(max(umath.absolute(x)),
                                   max(umath.absolute(y)))
      if step:
         for ii in xrange(self.numvect):
            sleep(0.01)
            plotxy(y[ii:ii+2], x[ii:ii+2],
                   rangex=[-amp, amp], rangey=[-amp, amp],
                   labx='Real Amplitude', laby='Imaginary Amplitude',
                   color=color, aspect=1.0)
      else:
         plotxy(y, x, rangex=[-amp, amp], rangey=[-amp, amp],
                labx='Real Amplitude', laby='Imaginary Amplitude',
                color=color, aspect=1.0)
   def plot(self, amp=None, color='white', device='/XWIN', step=0):
      """
      dftvector.plot(amp=None, color='white', device='/XWIN', step=0)
           Plot the dftvector summation.
      """
      self.plotvector(self.vector, amp=amp, color=color,
                      device=device, step=step)

x=1
if (x==0):
   dv = dftvector("/raid/data/Ter5_98/Day1-2/Ter5_1-2_bary_DM236_3564593.380.dftvec")
   phs = 4.50
   rfold = dv.r
   rmeas = 3564593.40
   z = -0.3
   przw = [phs, (rmeas - 0.5 * z) - rfold, z, 0.0]
elif (x==1):
   dv = dftvector("/raid/data/NGC6544/NGC6544_DM134_bary_566514.000.dftvec")
   phs = 2.63
   rfold = dv.r
   rmeas = 566511.93
   z = 8.6
   przw = [phs, (rmeas - 0.5 * z) - rfold, z, 0.0]
else:
   dv = dftvector("testz_1000000.000.dftvec")
   phs = 3.88
   rfold = dv.r
   rmeas = 1000003.24
   z = 5.8
   przw = [phs, (rmeas - 0.5 * z) - rfold, z, 0.0]
   dv.rebin(4)

if (1):
   print "    Max amplitude =", dv.maxamplitude()
   print "Current amplitude =", dv.amplitude()
   dv.rzwinfo(przw)
   phasecorr = dv.rzw_phases(przw)
   corrvect = add_components(dv.rotate(-phasecorr))
   dv.plotvector(corrvect, color="red")
   amplitudes = Numeric.ones(dv.numvect, 'D')
   corrvect = add_components(rotate_components(amplitudes, phasecorr))
   norm = dv.amplitude() / vector_amplitude(corrvect)
   dv.plotvector(norm * corrvect, color="blue")
   dv.plot()
   closeplot()
   
   przw = dv.fitrzw(przw)
   dv.rzwinfo(przw)
   phasecorr = dv.rzw_phases(przw)
   corrvect = add_components(dv.rotate(-phasecorr))
   dv.plotvector(corrvect, color="red")
   amplitudes = Numeric.ones(dv.numvect, 'D')
   corrvect = add_components(rotate_components(amplitudes, phasecorr))
   norm = dv.amplitude() / vector_amplitude(corrvect)
   dv.plotvector(norm * corrvect, color="blue")
   dv.plot()
   closeplot()
   
   przw = dv.fitrzw(przw, rzw_phase_model)
   dv.rzwinfo(przw)
   phasecorr = dv.rzw_phases(przw)
   corrvect = add_components(dv.rotate(-phasecorr))
   dv.plotvector(corrvect, color="red")
   amplitudes = Numeric.ones(dv.numvect, 'D')
   corrvect = add_components(rotate_components(amplitudes, phasecorr))
   norm = dv.amplitude() / vector_amplitude(corrvect)
   dv.plotvector(norm * corrvect, color="blue")
   dv.plot()
   closeplot()
else:
   dv = dftvector("/raid/data/Ter5_98/Day1-2/Ter5_1-2_bary_DM243_2600530.000.dftvec")
   dv.plot(step=0.03)
   closeplot()
   plotxy(abs(dv.components), dv.timefract)
   dv.mask(0.0, 0.33)
   dv.mask(0.41, 0.51)
   dv.mask(0.60, 1.01)
   plotxy(abs(dv.components), dv.timefract, color='red')
   closeplot()
   dv.plot(color='red')
   closeplot()
   plotxy(dv.phases(), dv.timefract)
