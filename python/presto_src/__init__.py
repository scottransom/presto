from __future__ import print_function
from __future__ import absolute_import
from builtins import input
from builtins import range
from .prestoswig import *
import os.path
import numpy as np
from presto import Pgplot
from presto import psr_utils


def val_with_err(value, error, length=0, digits=2, latex=0):
    """
    val_with_err(value, error, length=0, digits=2):
        Returns a string of length length (auto if 0) with 'value'
            rounded to the appropriate decimal place and the
            'error' in parenthesis as in scientific journals.
            The error has 'digits' decimal places.    
        Notes:
           'length' should be ~20 to show full double precision          
                if the base 10 exponent of the error needs to be shown.       
            If length == 0, left-justified minimum length string is returned.
            If length > 0, the string returned is right justified.       
            If length < 0, the string returned is left justified.       
            If latex=1, the string is converted into LaTeX markup.
    """
    slen = 40
    outstr = ' ' * slen
    if abs(length) > slen:
        slen = abs(length)
    if digits == 2:
        slen = nice_output_2(outstr, value, error, length)
    else:
        slen = nice_output_1(outstr, value, error, length)
    outstr = outstr[:slen].strip()  # remove null termination and any space
    if length < 0:
        outstr = outstr + (20 - len(outstr)) * ' '
    if length > 0:
        outstr = (20 - len(outstr)) * ' ' + outstr
    if latex:
        if outstr.find("x10") > 0:
            outstr = outstr.replace("x10^", r"$\times$10$^{") + "}$"
    return outstr


def read_inffile(filename, verbose=True):
    """
    read_inffile(filename, verbose=True):
        Return an infodata 'C' structure containing the data from the
           'inf' file in 'filename'.
    """
    fname = filename[:-4] if (filename[-4:] == ".inf") else filename
    id = infodata()
    if verbose:
        print("Reading information from", "\"" + fname + ".inf\"")
    readinf(id, fname)
    return id


def write_inffile(infodata, verbose=True):
    """
    wite_inffile(infodata, verbose=True):
        Write an '.inf' file based on its input structure
    """
    if verbose:
        print("Writing .inf file to '%s.inf'" % infodata.name)
    writeinf(infodata)


def psrepoch(psrname, epoch, verbose=True):
    """
    psrepoch(psrname or parname, epoch):
        Return a psrparams 'C' structure which includes data for
            PSR 'psrname' (a string of the B1950 or J2000 name of the
            pulsar -- without PSR, J, or B included) at epoch 'epoch'
            (in MJD format) from the ATNF database, or, a parfile is
            passed, read the pulsar information from it instead.
    """
    pp = psrparams()
    if os.path.isfile(psrname):
        get_psr_from_parfile("1903+0327.par", epoch, pp)
        if verbose:
            print('Retrieved data at MJD %f from "%s"' % (epoch, psrname))
    else:
        num = get_psr_at_epoch(psrname, epoch, pp)
        if verbose:
            print('Retrieved data at MJD %f for %s' % (epoch, pp.jname))
            print('The pulsar was #%d in the database.' % num)
    return pp


def read_rzwcands(filename):
    """
    read_rzwcands(filename):
        Return a list of all of the rzw search candidates from
        the file 'filename'.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = fourierprops()
    while (read_rzw_cand(infile, nextcand)):
        cands.append(nextcand)
        nextcand = fourierprops()
    infile.close()
    return cands


def read_rawbincands(filename):
    """
    read_rawbincands(filename):
        Return a list of all of the raw binary search candidates
            from the file 'filename'.
    """
    infile = open(filename, "r")
    cands = []
    nextcand = rawbincand()
    while (read_rawbin_cand(infile, nextcand)):
        cands.append(nextcand)
        nextcand = rawbincand()
    infile.close()
    return cands


def next2_to_n(x):
    """
    next2_to_n(x):
        Return the first value of 2^n >= x.
    """
    i = 1
    while (i < x): i = i << 1
    return i


def rfft(data, sign=-1):
    """
    rfft(data, sign=-1):
        Return the FFT of the real-valued, 32-bit floating point 'data'
        Note:  This only returns the positive frequency half of the FFT,
            since the other half is symmetric.  The Nyquist frequency
            is stored in the complex part of frequency 0 as per
            Numerical Recipes.
        The optional value 'sign' should be -1 (forward) or +1 (inverse).
    """
    # Default to sign = -1 if the user gives a bad value
    if (sign == -1 or sign != 1):
        tmp = np.array(data, copy=1).astype(np.float32)
        realfft(tmp, -1)
        return tmp.view(np.complex64)
    else:
        tmp = np.array(data.view(np.float32), copy=1).astype(np.float32)
        realfft(tmp, 1)
        return tmp.view(np.float32)


def spectralpower(fftarray):
    """
    spectralpower(fftarray):
        Return the power spectrum of a complex FFT 'fftarray'.
    """
    return power_arr(np.asarray(fftarray).astype(np.complex64))


def spectralphase(fftarray):
    """
    spectralphase(fftarray):
        Return the spectral phase (deg) of a complex FFT 'fftarray'.
    """
    return phase_arr(np.asarray(fftarray).astype(np.complex64))


def rzw_response(roffset, z, w, numbetween=1, numkern=None):
    """
    rzw_response(roffset, z, w, numbetween=1, numkern=None):
        Return the response of a signal offset from a Fourier bin
            by roffset bins, with a Fourier f-dot of z, and a
            Fourier f-dotdot of w.  The Fourier interpolation
            factor is the integer numbetween, and the the length
            of the resulting kernel will be auto-determined if
            numkern is None.
    """
    if numkern is None:
        numkern = w_resp_halfwidth(z, w, LOWACC)
    return gen_w_response(roffset, numbetween, numkern, z, w)


def maximize_r(data, r, norm=None):
    """
    maximize_r(data, r, norm = None):
        Optimize the detection of a signal at Fourier frequency 'r' in
            a FFT 'data'.  The routine returns a list containing
            the optimized values of the maximum normalized power, rmax,
            and an rderivs structure for the peak.
    """
    rd = rderivs()
    (rmax, maxpow) = max_r_arr(data, r, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, rd]


def maximize_rz(data, r, z, norm=None):
    """
    maximize_rz(data, r, z, norm = None):
        Optimize the detection of a signal at location 'r', 'z' in
            the F-Fdot plane.  The routine returns a list containing
            the optimized values of the maximum normalized power, rmax,
            zmax, and an rderivs structure for the peak.
    """
    rd = rderivs()
    (rmax, zmax, maxpow) = max_rz_arr(data, r, z, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, zmax, rd]


def maximize_rz_harmonics(data, r, z, numharm, norm=None):
    """
    maximize_rz_harmonics(data, r, z, numharm, norm = None):
        Optimize the detection of a signal at location 'r', 'z' in
            the F-Fdot plane, including harmonic summing of the harmonics.
            The routine returns a list containing the optimized values of 
            the maximum normalized power, rmax, zmax, and a list of 
            rderivs structures for the peak.
    """
    rds = [rderivs() for ii in range(numharm)]
    derivdata = np.zeros(7 * numharm, dtype=np.float64)
    rmax, zmax = max_rz_arr_harmonics(data, r, z, derivdata)
    maxpow = 0.0
    for ii in range(numharm):
        rds[ii].pow = derivdata[ii * 7 + 0]
        rds[ii].phs = derivdata[ii * 7 + 1]
        rds[ii].dpow = derivdata[ii * 7 + 2]
        rds[ii].dphs = derivdata[ii * 7 + 3]
        rds[ii].d2pow = derivdata[ii * 7 + 4]
        rds[ii].d2phs = derivdata[ii * 7 + 5]
        rds[ii].locpow = derivdata[ii * 7 + 6]
        maxpow += rds[ii].pow / rds[ii].locpow if norm is None else rds[ii].pow / norm
    return [maxpow, rmax, zmax, rds]


def maximize_rzw(data, r, z, w, norm=None):
    """
    maximize_rzw(data, r, z, w, norm = None):
        Optimize the detection of a signal at location 'r', 'z', 'w' in
            the F-Fdot-Fdotdot plane.  The routine returns a list containing
            the optimized values of the maximum normalized power, rmax,
            zmax, wmax, and an rderivs structure for the peak.
    """
    rd = rderivs()
    (rmax, zmax, wmax, maxpow) = max_rzw_arr(data, r, z, w, rd)
    maxpow = maxpow / rd.locpow if norm is None else maxpow / norm
    return [maxpow, rmax, zmax, wmax, rd]


def maximize_rzw_harmonics(data, r, z, w, numharm, norm=None):
    """
    maximize_rzw_harmonics(data, r, z, w, numharm, norm = None):
        Optimize the detection of a signal at location 'r', 'z', 'w' in
            the F-Fd-Fdd volume, including harmonic summing of the harmonics.
            The routine returns a list containing the optimized values of 
            the maximum normalized power, rmax, zmax, wmax, and a list of 
            rderivs structures for the peak.
    """
    rds = [rderivs() for ii in range(numharm)]
    derivdata = np.zeros(7 * numharm, dtype=np.float64)
    rmax, zmax, wmax = max_rzw_arr_harmonics(data, r, z, w, derivdata)
    maxpow = 0.0
    for ii in range(numharm):
        rds[ii].pow = derivdata[ii * 7 + 0]
        rds[ii].phs = derivdata[ii * 7 + 1]
        rds[ii].dpow = derivdata[ii * 7 + 2]
        rds[ii].dphs = derivdata[ii * 7 + 3]
        rds[ii].d2pow = derivdata[ii * 7 + 4]
        rds[ii].d2phs = derivdata[ii * 7 + 5]
        rds[ii].locpow = derivdata[ii * 7 + 6]
        maxpow += rds[ii].pow / rds[ii].locpow if norm is None else rds[ii].pow / norm
    return [maxpow, rmax, zmax, wmax, rds]


def search_fft(data, numcands, norm='default'):
    """
    search_fft(data, numcands):
        Search a short FFT and return a list containing the powers and
        Fourier frequencies of the 'numcands' highest candidates in 'data'.
        'norm' is the value to multiply each pow power by to get
             a normalized power spectrum (defaults to  1.0/(Freq 0) value)
    """
    if (norm == 'default'): norm = 1.0 / data[0].real
    hp = np.zeros(numcands, 'f')
    hf = np.zeros(numcands, 'f')
    search_minifft(data, len(data), norm, numcands, hp, hf)
    cands = []
    for i in range(numcands):
        cands.append([hp[i], hf[i]])
    return cands


def ffdot_plane(data, lor, dr, numr, loz, dz, numz):
    """
    ffdot_plane(data, lor, dr, numr, loz, dz, numz):
         Generate an F-Fdot plane with the 'lower-left' corners
         at the point 'lor', 'loz'.  The plane will have 'numr' frequency
         bins and 'numz' slices in the fdot direction, separated by 'dr'
         and 'dz' respectively.  'lor', 'numr', and 'numz' should all be
         integers.  'data' is the input FFT.
         Note:  'dr' much be the reciprocal of an integer
              (i.e. 1 / numbetween).  Also, 'r' is considered to be
              the average frequency (r = ro + z / 2).
    """
    lor = int(lor)
    numr = int(numr)
    numz = int(numz)
    numbetween = int(1.0 / dr)
    hiz = loz + (numz - 1) * dz
    maxabsz = max(abs(loz), abs(hiz))
    kern_half_width = z_resp_halfwidth(maxabsz, LOWACC)
    fftlen = next2_to_n(numr + 2 * numbetween * kern_half_width)
    ffd = corr_rz_plane(data, numbetween, lor, loz, hiz,
                        numz, fftlen, LOWACC)
    return np.array(ffd[:, 0:numr], copy=1)


def fdotdot_vol(data, lor, dr, numr, loz, dz, numz, low, dw, numw):
    """
    fdotdot_vol(data, lor, dr, numr, loz, dz, numz, low, dw, numw):
         Generate an F-Fdot-Fdotdot volume with the 'lower-left' corners
         at the point 'lor', 'loz', 'low'.  The vol will have 'numr' frequency
         bins, 'numz'/'numw' slices in the fdot/fdotdot direction, separated 
         by 'dr', 'dz', and 'dw' respectively.  'lor', 'numr', 'numz', and 
         'numw' should all be integers.  'data' is the input FFT.
         Note:  'dr' much be the reciprocal of an integer
              (i.e. 1 / numbetween).  Also, 'r' is considered to be
              the average frequency (r = r0 + w/6 + z0/2), and 'z'
              is the average fdot (z = z0 + w / 2).
    """
    lor = int(lor)
    numr, numz, numw = int(numr), int(numz), int(numw)
    numbetween = int(1.0 / dr)
    hiz = loz + (numz - 1) * dz
    maxabsz = max(abs(loz), abs(hiz))
    hiw = low + (numw - 1) * dw
    maxabsw = max(abs(low), abs(hiw))
    kern_half_width = w_resp_halfwidth(maxabsz, maxabsw, LOWACC)
    fftlen = next2_to_n(numr + 2 * numbetween * kern_half_width)
    ffd = corr_rzw_vol(data, numbetween, lor, loz, hiz,
                       numz, low, hiw, numw, fftlen, LOWACC)
    return np.array(ffd[:, :, 0:numr], copy=1)


def estimate_rz(psr, T, show=0, device='/XWIN'):
    """
    estimate_rz(psr, T, show=0, device='/XWIN'):
        Return estimates of a pulsar's average Fourier freq ('r')
        relative to its nominal Fourier freq as well as its
        Fourier f-dot ('z') in bins, of a pulsar.
           'psr' is a psrparams structure describing the pulsar.
           'T' is the length of the observation in sec.
           'show' if true, displays plots of 'r' and 'z'.
           'device' if the device to plot to if 'show' is true.
    """
    startE = keplers_eqn(psr.orb.t, psr.orb.p, psr.orb.e, 1.0E-15)
    numorbpts = int(T / psr.orb.p + 1.0) * 1024 + 1
    dt = T / (numorbpts - 1)
    E = dorbint(startE, numorbpts, dt, psr.orb)
    z = z_from_e(E, psr, T)
    r = T / p_from_e(E, psr) - T / psr.p
    if show:
        times = np.arange(numorbpts) * dt
        Pgplot.plotxy(r, times, labx='Time', \
                      laby='Fourier Frequency (r)', device=device)
        if device == '/XWIN':
            print('Press enter to continue:')
            try:
                i = raw_input()
            except NameError:
                i = input()
        Pgplot.nextplotpage()
        Pgplot.plotxy(z, times, labx='Time',
                      laby='Fourier Frequency Derivative (z)', device=device)
        Pgplot.closeplot()
    return r.mean(), z.mean()


def alias(r, rny):
    """
    alias_to_r(r, rny):
        Convert an aliased Fourier frequency into the 'true' Fourier
        frequency of a signal.  Or vise-versa -- the transformation is
        symmetric about the Nyquist Freq.
           'r' is the signal's Fourier frequency to convert.
           'rny' is the Nyquist frequency (in bins).  For an FFT
              of real data, 'rny' = number of data points FFT'd / 2.
    """
    return 2.0 * rny - r


def show_ffdot_plane(data, r, z, dr=0.125, dz=0.5,
                     numr=300, numz=300, T=None,
                     contours=None, title=None,
                     image="astro", device="/XWIN", norm=1.0):
    """
    show_ffdot_plane(data, r, z):
        Show a color plot of the F-Fdot plane centered on the point 'r', 'z'.
    """
    ffdp = ffdot_plane(data, r, dr, numr, z, dz, numz)
    ffdpow = spectralpower(ffdp.ravel())
    ffdpow.shape = (numz, numr)
    startbin = int(r - (numr * dr) / 2)
    startz = int(z - (numz * dz) / 2)
    x = np.arange(numr, dtype="d") * dr + startbin
    y = np.arange(numz, dtype="d") * dz + startz
    highpt = np.argmax(ffdpow.ravel())
    hir = highpt % numr
    hiz = highpt / numr
    print("")
    print("Fourier Freqs from ", min(x), "to", max(x), ".")
    print("Fourier Fdots from ", min(y), "to", max(y), ".")
    print("Maximum normalized power is ", ffdpow[hiz][hir])
    print("The max value is located at:  r =", startbin + hir * dr, \
          "  z =", startz + hiz * dz)
    print("")
    if not T:
        Pgplot.plot2d(ffdpow, x, y, labx="Fourier Frequency (bins)", \
                      laby="Fourier Frequency Derivative", \
                      title=title, image=image, \
                      contours=contours, device=device)
    else:
        Pgplot.plot2d(ffdpow, x / T, y / (T ** 2.0), labx="Frequency (hz)", \
                      laby="Frequency Derivative (Hz/sec)", \
                      rangex2=[x[0], x[-1]], rangey2=[y[0], y[-1]], \
                      labx2="Fourier Frequency", \
                      laby2="Fourier Frequency Derivative", \
                      title=title, image=image, \
                      contours=contours, device=device)


def v_from_e(e, psr):
    """
    v_from_e(e, psr):
        Return a vector of velocities (km/s) from a vector of Eccentric
        anomalys.
            'e' is the vector of Eccentric anomalys.
            'psr' is a psrparams instance containing info about the pulsar.
    """
    oldw = psr.orb.w
    v = np.array(e, copy=1)
    E_to_v(v, psr.orb)
    psr.orb.w = oldw
    return v


def d_from_e(e, psr):
    """
    d_from_e(e, psr):
        Return a vector of time delays (s) from a vector of Eccentric
        anomalys.
            'e' is the vector of Eccentric anomalys.
            'psr' is a psrparams instance containing info about the pulsar.
    """
    oldw = psr.orb.w
    d = np.array(e, copy=1)
    E_to_phib(d, psr.orb)
    psr.orb.w = oldw
    return d


def p_from_e(e, psr):
    """
    p_from_e(e, psr):
        Return a vector of pulsar periods (s) from a vector of Eccentric
        anomalys.
            'e' is the vector of Eccentric anomalys.
            'psr' is a psrparams instance containing info about the pulsar.
    """
    oldw = psr.orb.w
    psr.orb.w = psr.orb.w * DEGTORAD
    p = np.array(e, copy=1)
    E_to_p(p, psr.p, psr.orb)
    psr.orb.w = oldw
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
    oldw = psr.orb.w
    psr.orb.w = psr.orb.w * DEGTORAD
    z = np.array(e, copy=1)
    E_to_z(z, psr.p, T, psr.orb)
    psr.orb.w = oldw
    return z


def pcorr(data, kernel, numbetween, lo, hi):
    """
    pcorr(data, kernel, numbetween, lo, hi):
        Perform a correlation with the raw complex vectors 'data' and
        'kernel'.  The returned vector should start at frequency
        'lo' (must be an integer), and go up to but not include 'hi'
        (also an integer).
    """
    kern_half_width = len(kernel) / (2 * numbetween)
    result = np.zeros((hi - lo) * numbetween, 'F')
    corr_complex(data, len(data), RAW,
                 kernel, len(kernel), RAW,
                 result, len(result), lo,
                 numbetween, kern_half_width, CORR)
    return result


def p_to_f(p, pd, pdd):
    """
    p_to_f(p, pd, pdd):
       Convert period, period derivative and period second
       derivative to the equivalent frequency counterparts.
       Will also convert from f to p.
    """
    f = 1.0 / p
    fd = -pd / (p * p)
    if (pdd == 0.0):
        fdd = 0.0
    else:
        fdd = 2.0 * pd * pd / (p ** 3.0) - pdd / (p * p)
    return [f, fd, fdd]


def bary_to_topo(pb, pbd, pbdd, infofilenm, ephem="DE200"):
    """
    bary_to_topo(pb, pbd, pbdd, infofilenm, ephem="DE200"):
       Use least squares to calculate topocentric period
       period derivative, and period second derivative
       for the corresponding barycentric values.  The data
       for the observation must be found in the info file.
    """
    from numpy.linalg.old import linear_least_squares
    if infofilenm[-4:] == ".inf":  infofilenm = infofilenm[:-4]
    obs = read_inffile(infofilenm)
    T = obs.N * obs.dt
    dt = 10.0
    tto = obs.mjd_i + obs.mjd_f
    tts = np.arange(tto, tto + (T + dt) / SECPERDAY, dt / SECPERDAY)
    nn = len(tts)
    bts = np.zeros(nn, 'd')
    vel = np.zeros(nn, 'd')
    ra = psr_utils.coord_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
    dec = psr_utils.coord_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
    if (obs.telescope == 'Parkes'):
        tel = 'PK'
    elif (obs.telescope == 'Effelsberg'):
        tel = 'EB'
    elif (obs.telescope == 'Arecibo'):
        tel = 'AO'
    elif (obs.telescope == 'MMT'):
        tel = 'MT'
    else:
        print("Telescope not recognized.")
        return 0
    barycenter(tts, bts, vel, nn, ra, dec, tel, ephem)
    print("Topocentric start time = %17.11f" % tts[0])
    print("Barycentric start time = %17.11f" % bts[0])
    avgvel = np.add.reduce(vel) / nn
    print("Average Earth velocity = %10.5e c" % (avgvel))
    tts = np.arange(nn, dtype='d') * dt
    bts = (bts - bts[0]) * SECPERDAY
    [fb, fbd, fbdd] = p_to_f(pb, pbd, pbdd)
    b = fb * bts + fbd * bts ** 2.0 / 2.0 + fbdd * bts ** 3.0 / 6.0
    a = np.transpose(np.asarray([tts, tts ** 2.0, tts ** 3.0]))
    [ft, ftd, ftdd], residuals, rank, sv = linear_least_squares(a, b)
    [pt, ptd, ptdd] = p_to_f(ft, ftd, ftdd)
    print("    Topocentric period = %15.12f" % pt)
    print("     Topocentric p-dot = %15.9e" % ptd)
    print("  Topocentric p-dotdot = %15.9e" % ptdd)
    print("     Quick Topo period = %15.12f" % (pb * (1.0 + avgvel)))
    print("      Quick Topo p-dot = %15.9e" % (pbd * (1.0 + avgvel)))
    print("   Quick Topo p-dotdot = %15.9e" % (pbdd * (1.0 + avgvel)))
    return [pt, ptd, ptdd]


def measure_phase(profile, template, sigma, fwhm):
    """
    measure_phase(profile, template, sigma, fwhm):
       TOA measurement technique from J. H. Taylor's talk
       _Pulsar_Timing_and_Relativistic_Gravity_.  Routine
       takes two profiles, the first measured and the
       second a high S/N template and determines the phase
       offset of 'profile' from 'template'.  Both profiles
       must have the same number of points.  'sigma' denotes
       the RMS noise level of the 'profile'.  'fwhm' is the
       approximate width of the template pulse (0-1).  The phase
       returned is cyclic (i.e. from 0-1).  The routine
       returns a tuple comtaining (tau, tau_err, b, b_err, a).
       Where 'tau' is the phase, 'B' is the scaling factor,
       and 'a' is the DC offset.  The error values are
       estimates of the 1 sigma errors.
    """
    from simple_roots import newton_raphson
    N = len(profile)
    if not (N == len(template)):
        print("Lengths of 'profile' and 'template' must")
        print("  be equal in measure_phase().")
        return 0.0
    ft = rfft(profile)
    p0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    P_k = abs(ft)
    frotate(P_k, len(ft), 1)
    Theta_k = np.arctan2(-ft.imag, ft.real)
    frotate(Theta_k, len(ft), 1)
    ft = rfft(template)
    s0 = ft[0].real
    # Nyquist freq
    ft[0] = complex(ft[0].imag, 0.0)
    S_k = abs(ft)
    frotate(S_k, len(ft), 1)
    Phi_k = np.arctan2(-ft.imag, ft.real)
    frotate(Phi_k, len(ft), 1)
    # Estimate of the noise sigma (This needs to be checked)
    # Note:  Checked 10 Jul 2000.  Looks OK.
    sig = sigma * np.sqrt(N)
    k = np.arange(len(ft), dtype='d') + 1.0

    def fn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
        # Since Nyquist freq always has phase = 0.0
        k[-1] = 0.0
        return np.add.reduce(k * p * s *
                             np.sin(phi - theta + k * tau))

    def dfn(tau, k=k, p=P_k, s=S_k, theta=Theta_k, phi=Phi_k):
        # Since Nyquist freq always has phase = 0.0
        k[-1] = 0.0
        return np.add.reduce(k * k * p * s *
                             np.cos(phi - theta + k * tau))

    numphases = 200
    ddchidt = np.zeros(numphases, 'd')
    phases = np.arange(numphases, dtype='d') / \
             float(numphases - 1) * TWOPI - PI
    for i in np.arange(numphases):
        ddchidt[i] = dfn(phases[i])
    maxdphase = phases[np.argmax(ddchidt)] + \
                0.5 * TWOPI / (numphases - 1.0)
    # Solve for tau
    tau = newton_raphson(fn, dfn, maxdphase - 0.5 * fwhm * TWOPI,
                         maxdphase + 0.5 * fwhm * TWOPI)
    # Solve for b
    c = P_k * S_k * np.cos(Phi_k - Theta_k + k * tau)
    d = np.add.reduce(S_k ** 2.0)
    b = np.add.reduce(c) / d
    # tau sigma
    tau_err = sig * np.sqrt(1.0 / (2.0 * b *
                                   np.add.reduce(k ** 2.0 * c)))
    # b sigma  (Note:  This seems to be an underestimate...)
    b_err = sig * np.sqrt(1.0 / (2.0 * d))
    # Solve for a
    a = (p0 - b * s0) / float(N)
    return (tau / TWOPI, tau_err / TWOPI, b, b_err, a)


def get_baryv(ra, dec, mjd, T, obs="PK"):
    """
    get_baryv(ra, dec, mjd, T, obs="PK"):
      Determine the average barycentric velocity towards 'ra', 'dec'
      during an observation from 'obs'.  The RA and DEC are in the
      standard string format (i.e. 'hh:mm:ss.ssss' and 'dd:mm:ss.ssss').
      'T' is in sec and 'mjd' is (of course) in MJD.  The obs variable
      is the standard two character string from TEMPO:  PK, GB, AO, GM, JB, ...
    """
    tts = np.linspace(mjd, mjd + T / 86400.0, 100)
    nn = len(tts)
    bts = np.zeros(nn, dtype=np.float64)
    vel = np.zeros(nn, dtype=np.float64)
    barycenter(tts, bts, vel, ra, dec, obs, "DE421")
    return vel.mean()


def fold(indata, dt, nbins, f, fd=0.0, fdd=0.0, startphs=0.0, tlo=0.0, standard=True):
    """
    fold(indata, dt, nbins, f, fd=0.0, fdd=0.0, startphs=0.0, tlo=0.0):
      This is an interface into PRESTO's fold() code, which is what
      prepfold uses to fold data.  It will return a tuple of a
      double-precision profile of length nbins, and the ending phase
      (0-1) of the fold.
        indata is an array of floats to fold
        dt is the duration in sec of each of the indata bins
        f, fd, and fdd are the freq, freq deriv, and freq 2nd deriv to fold (Hz)
        startphs (0-1) is the phase for the beginning of the first bin
        tlo is the time (in sec) referring to the start of the first bin,
          with respect to the reference time of f, fd, and fdd (i.e. tlo=0.0).
        If standard (bool), then traditional prepfold "drizzling" will be
          used.  Otherwise, treat each input data point as a sample and put
          it fully in a single profile bin.
    """
    prof = np.zeros(nbins, dtype=np.float64)
    data = indata.astype(np.float32)
    phs = simplefold(data, dt, tlo, prof, startphs, f, fd, fdd,
                     1 if standard else 0)
    return (prof, phs)

