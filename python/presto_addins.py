

#-------------- Extra Stuff to Make Things Easier -----------------

import math, Numeric, Pgplot

def read_inffile(filename):
   """
   read_inffile(filename):
       Return an infodata 'C' structure containing the data from the
       'inf' file in 'filename'.  'filename' should not include the
       '.inf' suffix.
   """
   id = infodata()
   print "Reading information from", "\""+filename+".inf\""
   readinf(id, filename)
   return id

def read_makfile(filename):
   """
   read_makfile(filename):
       Return an makedata 'C' structure containing the data from the
       'mak' in 'filename'.  'filename' should not include the
       '.mak' suffix.
   """
   md = makedata()
   read_mak_file(filename, md)
   return md

def psrepoch(psrname, epoch):
   """
   psrepoch(psrname, epoch):
       Return a psrparams 'C' structure which includes data for
           PSR 'psrname' (a string of the B1950 or J2000 name of the
           pulsar -- without PSR, J, or B included) at epoch 'epoch'
           (in MJD format).
   """
   pp = psrparams()
   num = return_psrparams_at_epoch(pp, psrname, epoch)
   print 'Retrieved data at MJD %f for %s' % (epoch, pp.jname)
   print 'The pulsar was #%d in the database.' % num
   return pp

def collect_psrdata():
    """
    collect_psrdata():
        Return a list of all of the pulsars in the Taylor et al.
            pulsar database including their characteristics.
    """
    pdata = []
    np = num_psrs_in_database()
    print 'There are %d pulsars in the database.' % np
    for i in range(0, np):
        pdata.append(psrdata())
        get_psrdata_by_num(pdata[i], i)
    return pdata

def next2_to_n(x):
    """
    next2_to_n(x):
        Return the first value of 2^n >= x.
    """
    i = 1L
    while (i < x): i = i << 1
    return i

def rfft(data, sign=-1):
   """
   rfft(data, sign=-1):
       Return the FFT of the real-valued 'data'.
       Note:  This only returns the positive frequency half of the FFT,
              since the other half is symmetric.  The Nyquist frequency
              is stored in the complex part of frequency 0 as per
              Numerical Recipes.
       The optional value 'sign' should be positive or negative 1.
   """
   # Default to sign = -1 if the user gives a bad value
   tmp = Numeric.array(data, copy=1)
   if (sign == -1 or sign != 1):
      tmp = tofloatvector(tmp)
      realfft(tmp, len(tmp), -1)
      float_to_complex(tmp)
   else:
      complex_to_float(tmp)
      realfft(tmp, len(tmp), 1)
   return tmp

def spectralpower(fftarray):
    """
    spectralpower(fftarray):
        Return the power spectrum of a complex FFT 'fftarray'.
    """
    fftarray = Numeric.asarray(fftarray)
    if fftarray.typecode()=='F':
       return power_arr(fftarray, len(fftarray))
    elif fftarray.typecode()=='D':
       return dpower_arr(fftarray, len(fftarray))
    else:
       print 'fftarray must be complex in spectralpower()'
       return None
    
def spectralphase(fftarray):
    """
    spectralphase(fftarray):
        Return the spectral phase (deg) of a complex FFT 'fftarray'.
    """
    fftarray = Numeric.asarray(fftarray)
    if fftarray.typecode()=='F':
       return phase_arr(fftarray, len(fftarray))
    elif fftarray.typecode()=='D':
       return dphase_arr(fftarray, len(fftarray))
    else:
       print 'fftarray must be complex in spectralpower()'
       return None

def maximize_rz(data, r, z, norm = None):
   """
   maximize_rz(data, r, z, norm = None):
       Optimize the detection of a signal at location 'r', 'z' in
           the F-Fdot plane.  The routine returns a list containing
           the optimized values of the maximum normalized power, rmax,
           zmax, and an rderivs structure for the peak.
   """
   rd = rderivs()
   (maxpow, rmax, zmax) = max_rz_arr(data, len(data), r, z, rd)
   if not norm:
      maxpow = maxpow / rd.locpow
   else:
      maxpow = maxpow / norm
   return [maxpow, rmax, zmax, rd]

def search_fft(data, numcands, norm='default'):
   """
   search_fft(data, numcands):
      Search a short FFT and return a list containing the powers and
      Fourier frequencies of the 'numcands' highest candidates in 'data'.
      'norm' is the value to multiply each pow power by to get
         a normalized power spectrum (defaults to  1.0/(Freq 0) value)
   """
   if (norm=='default'): norm = 1.0/data[0].real
   hp = Numeric.zeros(numcands, 'f')
   hf = Numeric.zeros(numcands, 'f')
   search_minifft(data, len(data), norm, numcands, hp, hf) 
   cands = []
   for i in range(numcands):
      cands.append([hp[i],hf[i]])
   return cands

def ffdot_plane(data, r, dr, numr, z, dz, numz):
   """
   ffdot_plane(data, r, dr, numr, z, dz, numz):
       Generate an F-Fdot plane centered on the point 'r', 'z'.
       There will be a total of 'numr' x 'numz' points in the array.
       The F-Fdot plane will be interpolated such the points are
       separated by 'dr' in the 'r' (i.e. f) direction and 'dz'
       in the 'z' (i.e. fdot) direction.  'data' is the input FFT.
       Note:  'dr' much be the reciprocal of an integer
              (i.e. 1 / numbetween).  Also, 'r' is considered to be
              the average frequency (r = ro + z / 2).
   """
   numbetween = int(1.0 / dr)
   startbin = int(r - (numr * dr) / 2)
   startz = int(z - (numz * dz) / 2)
   maxabsz = max(abs(startz), abs(startz + numz * dz))
   fftlen = next2_to_n(numr + 2 * numbetween * \
                       z_resp_halfwidth(maxabsz, LOWACC))
   ffdp = Numeric.zeros((numz, numr), 'F')
   for i in range(numz):
      z = startz + i * dz
      (ffdraw, nextbin) = corr_rz_interp(data, len(data), numbetween, \
                                         startbin, z, fftlen, LOWACC)
      ffdp[i][0:numr] = Numeric.array(ffdraw[0:numr], copy=1)
   return ffdp

def show_ffdot_plane(data, r, z, dr = 0.125, dz = 0.5,
                     numr = 300, numz = 300, 
                     contours = None,
                     title = None, 
                     image = "astro",
                     device = "/XWIN",
                     norm = 1.0):
   """
   show_ffdot_plane(data, r, z):
       Show a color plot of the F-Fdot plane centered on the point 'r', 'z'.
   """
   ffdp = ffdot_plane(data, r, dr, numr, z, dz, numz)
   ffdpow = spectralpower(ffdp.flat)
   ffdpow.shape = (numz, numr)
   startbin = int(r - (numr * dr) / 2)
   startz = int(z - (numz * dz) / 2)
   x = Numeric.arange(numr, typecode="d") * dr + startbin
   y = Numeric.arange(numz, typecode="d") * dz + startz
   highpt = Numeric.argmax(ffdpow.flat)
   hir = highpt % numr
   hiz = highpt / numr
   print ""
   print "Fourier Freqs from ", min(x), "to", max(x), "."
   print "Fourier Fdots from ", min(y), "to", max(y), "."
   print "Maximum normalized power is ", ffdpow[hiz][hir]
   print "The max value is located at:  r =", startbin + hir * dr, \
         "  z =", startz + hiz * dz
   print ""
   Pgplot.plot2d(ffdpow, x, y, labx = "Fourier Frequency (bins)", \
                 laby = "Fourier Frequency Derivative", \
                 title = title, image = image, \
                 contours = contours, device = device)

def make_orbit(psr, npts, dt, eo=0.0):
   """
   make_orbit(psr, npts, dt, eo=0.0):
       Return a Numeric array containing the Eccentric anomaly for a
       given orbit at a string of times.
           'psr' is a psrparams instance containing info about the pulsar.
           'npts' is the number of data points to return.
           'dt' is the difference in time (s) between returned points.
           'eo' is the initial Eccentric anomaly. (default = 0.0)
   """
   e = Numeric.zeros(npts, 'd')
   t = Numeric.zeros(npts, 'd')
   dorbint(e, eo, t, 0.0, npts, dt, psr.orb)
   return e

def v_from_e(e, psr):
   """
   v_from_e(e, psr):
       Return a vector of velocities (km/s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   v = Numeric.array(e, copy=1)
   E_to_v(v, len(v), psr.orb)
   return v

def d_from_e(e, psr):
   """
   d_from_e(e, psr):
       Return a vector of time delays (s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   d = Numeric.array(e, copy=1)
   E_to_phib(d, len(d), psr.orb)
   return d

def p_from_e(e, psr):
   """
   p_from_e(e, psr):
       Return a vector of pulsar periods (s) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
   """
   p = Numeric.array(e, copy=1)
   E_to_p(p, len(p), psr.p, psr.orb)
   return p

def z_from_e(e, psr, T):
   """
   z_from_e(e, psr):
       Return a vector of Fourier F-dots (bins) from a vector of Eccentric
       anomalys.
           'e' is the vector of Eccentric anomalys.
           'psr' is a psrparams instance containing info about the pulsar.
           'T' is the total length of the observation (s).
   """
   z = Numeric.array(e, copy=1)
   E_to_z(z, len(z), psr.p, T, psr.orb)
   return z


   
