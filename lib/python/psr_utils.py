import Numeric, umath, FFT, Pgplot, ppgplot, bisect, sinc_interp
from Scientific.Functions.FindRoot import newtonRaphson
from Scientific.Statistics import average, standardDeviation
from Scientific.Statistics.Histogram import Histogram
from scipy.special.cephes import ndtr, ndtri, chdtrc, chdtri, i0, kolmogorov
from scipy.optimize import leastsq
from psr_constants import *

isintorlong = lambda x: type(x) == type(0) or type(x) == type(0L)

def span(Min, Max, Num):
    """
    span(Min, Max, Num):
        Create a range of 'Num' floats given inclusive 'Min' and 'Max' values.
    """
    assert isintorlong(Num)
    if isintorlong(Min) and isintorlong(Max) and \
       (Max-Min) % (Num-1) != 0:
        Max = float(Max) # force floating points
    return Min+(Max-Min)*Numeric.arange(Num)/(Num-1)

def distance(width):
    """
    distance(width):
        Return a 'width' x 'width' Numeric Python array with each
            point set to the geometric distance from the array's center.
    """
    x = Numeric.arange(-width/2.0+0.5, width/2.0+0.5, 1.0)**2
    x = Numeric.resize(x, (width,width))
    return Numeric.sqrt(x + Numeric.transpose(x))

def running_avg(arr, navg):
    """
    running_avg(arr, navg):
        Return an array of the running average of 'navg' bins from the
        input array 'arr'.
    """
    a = Numeric.asarray(arr, 'd')
    a.shape = (len(a) / navg, navg)
    return umath.add.reduce(Numeric.transpose(a)) / navg

def hist(data, bins, range=None, labx="", color=1, line=1,
         device='/XWIN'):
    """
    hist(data, bins, range=None):
    Return and plot a histogram in one variable.
      data  -- a sequence of data points
      bins  -- the number of bins into which the data is to be sorted
      range -- a tuple of two values, specifying the lower and
               the upper end of the interval spanned by the bins.
               Any data point outside this interval will be ignored.
               If no range is given, the smallest and largest
               data values are used to define the interval.
    Note:  This command also accepts some basic flags for the plot,
           like labx, color, line, and device.
     """
    rawh = Histogram(data, bins, range)
    h = Numeric.transpose(rawh[:])
    Pgplot.plotbinned(h[1], h[0], rangey=[0,int(1.1*max(h[1]))],
                      laby="Number", color=color, line=line,
                      labx=labx, device=device)
    return h

def KS_test(data, cumdist, output=0):
    """
    KS_test(data, cumdist, output=0):
        Perform a Kolmogorov-Smirnov test on data compared to the
            cumulative-distribution function cumdist.
    """
    nn = len(data)
    sdata = Numeric.sort(Numeric.asarray(data))
    D1 = umath.maximum.reduce(umath.absolute(cumdist(sdata)-
                                             Numeric.arange(nn, typecode='d')/nn))
    D2 = umath.maximum.reduce(umath.absolute(cumdist(sdata)-
                                             Numeric.arange(1,nn+1, typecode='d')/nn))
    D = max((D1, D2))
    P = kolmogorov(sqrt(nn)*D)
    if (output):
        print "Max distance between the cumulative distributions (D) = %.5g" % D
        print "Prob the data is from the specified distrbution   (P) = %.3g" % P
    return (D, P)


def MJD_to_JD(MJD):
    """
    MJD_to_JD(MJD):
       Convert Modified Julian Date (MJD) to Julian Date (JD)
    """
    return MJD+2400000.5

def JD_to_MJD(JD):
    """
    JD_to_MJD(JD):
       Convert Julian Date (JD) to Modified Julian Date (MJD)
    """
    return JD-2400000.5

def MJD_to_Julian_Epoch(MJD):
    """
    MJD_to_Julian_Epoch(MJD):
       Convert Modified Julian Date (MJD) to Julian Epoch
    """
    return 2000.0 + (MJD-51544.5)/365.25

def Julian_Epoch_to_MJD(jepoch):
    """
    Julian_Epoch_to_MJD(jepoch):
       Convert Julian Epoch to Modified Julian Date (MJD)
    """
    return 51544.5 + (jepoch-2000.0)*365.25

def MJD_to_Besselian_Epoch(MJD):
    """
    MJD_to_Besselian_Epoch(MJD):
       Convert Modified Julian Date (MJD) to Besselian Epoch
    """
    return 1900.0 + (MJD-15019.81352)/365.242198781

def Besselian_Epoch_to_MJD(bepoch):
    """
    Besselian_Epoch_to_MJD(bepoch):
       Convert Besselian Epoch to Modified Julian Date (MJD)
    """
    return 15019.81352 + (bepoch-1900.0)*365.242198781

def rad_to_dms(rad):
    """
    rad_to_dms(rad):
       Convert radians to degrees, minutes, and seconds of arc.
    """
    if (rad < 0.0): sign = -1
    else: sign = 1
    arc = RADTODEG * umath.fmod(umath.fabs(rad), PI)
    d = int(arc)
    arc = (arc - d) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (sign * d, m, s)

def dms_to_rad(deg, min, sec):
    """
    dms_to_rad(deg, min, sec):
       Convert degrees, minutes, and seconds of arc to radians.
    """
    if (deg < 0.0): sign = -1
    else: sign = 1
    return sign * ARCSECTORAD * \
           (60.0 * (60.0 * umath.fabs(deg) +
                    umath.fabs(min)) + umath.fabs(sec))

def dms_to_deg(deg, min, sec):
    """
    dms_to_deg(deg, min, sec):
       Convert degrees, minutes, and seconds of arc to degrees.
    """
    return RADTODEG * dms_to_rad(deg, min, sec)

def rad_to_hms(rad):
    """
    rad_to_hms(rad):
       Convert radians to hours, minutes, and seconds of arc.
    """
    rad = umath.fmod(rad, TWOPI)
    if (rad < 0.0): rad = rad + TWOPI
    arc = RADTOHRS * rad
    h = int(arc)
    arc = (arc - h) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (h, m, s)

def hms_to_rad(hour, min, sec):
    """
    hms_to_rad(hour, min, sec):
       Convert hours, minutes, and seconds of arc to radians
    """
    if (hour < 0.0): sign = -1
    else: sign = 1
    return sign * SECTORAD * \
           (60.0 * (60.0 * umath.fabs(hour) +
                    umath.fabs(min)) + umath.fabs(sec))

def hms_to_hrs(hour, min, sec):
    """
    hms_to_hrs(hour, min, sec):
       Convert hours, minutes, and seconds of arc to hours.
    """
    return RADTOHRS * hms_to_rad(hour, min, sec)

def coord_to_string(h_or_d, m, s):
    """
    coord_to_string(h_or_d, m, s):
       Return a formatted string of RA or DEC values as
       'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
    """
    if (s >= 10.0):
        return "%.2d:%.2d:%.4f" % (h_or_d, m, s)
    else:
        return "%.2d:%.2d:0%.4f" % (h_or_d, m, s)

def ra_to_rad(ra_string):
    """
    ra_to_rad(ar_string):
       Given a string containing RA information as
       'hh:mm:ss.ssss', return the equivalent decimal
       radians.
    """
    h, m, s = ra_string.split(":")
    return hms_to_rad(int(h), int(m), float(s))

def dec_to_rad(dec_string):
    """
    dec_to_rad(dec_string):
       Given a string containing DEC information as
       'dd:mm:ss.ssss', return the equivalent decimal
       radians.
    """
    d, m, s = dec_string.split(":")
    return dms_to_rad(int(d), int(m), float(s))

def delta_m(flux_factor):
    """
    delta_m(flux_factor):
        Return the change in magnitudes caused by a change
            in flux of flux_factor.
    """
    return -2.5*umath.log10(flux_factor)

def flux_factor(delta_m):
    """
    flux_factor(delta_m):
        Return the change in flux caused by a change
            in magnitude of delta_m magnitudes
    """
    return 10.0**(delta_m/-2.5)

def distance_modulus_to_distance(dm, absorption=0.0):
    """
    distance_modulus_to_distance(dm, absorption=0.0):
        Return the distance (kpc) given a distance modulus dm and
            an optional absorption.
    """
    return 10.0**(((dm-absorption)+5.0)/5.0)/1000.0

def distance_to_distance_modulus(d, absorption=0.0):
    """
    distance_to_distance_modulus(d, absorption=0.0):
        Return the distance modulus given a distance d and
            an optional absorption.
    """
    return 5.0*umath.log10(d*1000.0)-5.0+absorption

def true_anomaly(E, ecc):
    """
    true_anomaly(E, ecc):
        Return the True Anomaly (in radians) given the Eccentric anomaly
            (E in radians) and the eccentricity (ecc)
    """
    return 2.0*umath.arctan(umath.sqrt((1.0+ecc)/(1.0-ecc))*umath.tan(E/2.0))

def mass_funct(pb, x):
    """
    mass_funct(pb, x):
        Return the mass function of an orbit given the following:
            'pb' is the binary period in sec.
            'x' is the projected semi-major axis in lt-sec.
    """
    return 8015123.37129 * x**3.0 / (pb * pb)

def mass_funct2(mp, mc, i):
    """
    mass_funct2(mp, mc, i):
        Return the mass function of an orbit given the following:
            'mp' is the mass of the primary in solar masses.
            'mc' is the mass of the companion in solar masses.
            'i' is the orbital inclination (rad).
        Note:  An 'average' orbit has cos(i) = 0.5, or i = 60 deg
    """
    return (mc * umath.sin(i))**3 / (mc + mp)**2

def asini_c(pb, mf):
    """
    asini_c(pb, mf):
        Return the orbital projected semi-major axis (lt-sec) given:
            'pb' is the binary period in sec.
            'mf' is the mass function of the orbit.
    """
    return (mf * pb * pb / 8015123.37129)**(1.0 / 3.0)


def bins_to_accel(rdot, T, f=[1.0, 1000.0], device="/XWIN"):
    """
    bins_to_accel(rdot, T, f=[1.0, 1000.0], device="/XWIN"):
        Make a plot showing the acceleration which corresponds
        to a certain number of Fourier bins drifted 'rdot' during
        an observation of length 'T'.
    """
    fs = span(umath.log10(f[0]), umath.log10(f[1]), 1000)
    accels = rdot * SOL / (T * T * 10.0**fs)
    if (device):
        Pgplot.plotxy(umath.log10(accels), fs, logx=1, logy=1,
                      labx="Frequency (Hz)",
                      laby="Acceleration (m/s\u2\d)", device=device)
        ppgplot.pgmtxt("T", -2.0, 0.75, 0.0, "T = %.0f sec"%T)
        ppgplot.pgmtxt("T", -3.5, 0.75, 0.0, "r\B\u\.\d = %.1f bins"%rdot)
        if (device != '/XWIN'):
            Pgplot.closeplot()
    else:
        return accels

def companion_mass(pb, x, inc=60.0, mpsr=1.35):
    """
    companion_mass(pb, x, inc=60.0, mpsr=1.35):
        Return the companion mass (in solar mass units) for a binary
        system with the following characteristics:
            'pb' is the binary period in sec.
            'x' is the projected semi-major axis in lt-sec.
            'inc' is the orbital inclination in degrees.
            'mpsr' is the mass of the pulsar in solar mass units.
    """
    massfunct = mass_funct(pb, x)
    def localmf(mc, mp=mpsr, mf=massfunct, i=inc*DEGTORAD):
        return (mc*umath.sin(i))**3.0/(mp + mc)**2.0 - mf
    return newtonRaphson(localmf, 0.0, 1000.0, 0.00000001)
        
def companion_mass_limit(pb, x, mpsr=1.35):
    """
    companion_mass_limit(pb, x, mpsr=1.35):
        Return the lower limit (corresponding to i = 90 degrees) of the
        companion mass (in solar mass units) in a binary system with
        the following characteristics:
            'pb' is the binary period in sec.
            'x' is the projected semi-major axis in lt-sec.
            'mpsr' is the mass of the pulsar in solar mass units.
    """
    return companion_mass(pb, x, inc=90.0, mpsr=mpsr)
        
def OMDOT(porb, e, Mp, Mc):
    """
    OMDOT(porb, e, Mp, Mc):
        Return the predicted advance of periaston (deg/yr) given the
        orbital period (s), eccentricity, and pulsar and companion masses.
    """
    return 3.0 * (porb/TWOPI)**(-5.0/3.0) * \
           (Tsun*(Mp+Mc))**(2.0/3.0) / (1.0-e**2.0) * \
           RADTODEG * SECPERJULYR

def gamma(porb, e, Mp, Mc):
    """
    gamma(porb, e, Mp, Mc):
        Return the predicted value of relativistic gamma (ms) given the
        orbital period (s), eccentricity, and pulsar and companion masses.
    """
    return e * (porb/TWOPI)**(1.0/3.0) * Tsun**(2.0/3.0) * \
           (Mp+Mc)**(-4.0/3.0) * Mc * (Mp+2.0*Mc) * 1000.0

def PBDOT(porb, e, Mp, Mc):
    """
    PBDOT(porb, e, Mp, Mc):
        Return the predicted orbital period derivative (s/s) given the
        orbital period (s), eccentricity, and pulsar and companion masses.
    """
    return -192.0*PI/5.0 * (porb/TWOPI)**(-5.0/3.0) * \
           (1.0 + 73.0/24.0*e**2.0 + 37.0/96.0*e**4.0) * \
           (1.0-e**2.0)**(-7.0/2.0) * Tsun**(5.0/3.0) * \
           Mp * Mc * (Mp+Mc)**(-1.0/3.0)


def OMDOT_to_Mtot(OMDOT, porb, e):
    """
    OMDOT_to_Mtot(OMDOT, porb, e):
        Return the total mass (in solar units) of a system given an advance
        of periastron (OMDOT) in deg/yr.  The opbital period should be in sec.
    """
    prob /= SECPERDAY
    wd = OMDOT/SECPERJULYR*DEGTORAD # rad/s
    return (wd/3.0*(1.0-e*e)*(porb*SECPERDAY/TWOPI)**(5.0/3.0))**(3.0/2.0)/Tsun

def beam_halfwidth(obs_freq, dish_diam):
    """
    beam_halfwidth(obs_freq, dish_diam):
        Return the telescope beam halfwidth in arcmin
            'obs_freq' = the observing frqeuency in MHz
            'dish_diam' = the telescope diameter in m
    """
    return 1.2*SOL/(obs_freq*10.0**6)/dish_diam*RADTODEG*60/2

def limiting_flux_dens(Ttot, G, BW, T, P=0.01, W=0.05, polar=2, factor=15.0):
    """
    limiting_flux_dens(Ttot, G, BW, T, P=0.01, W=0.05, polar=2, factor=15.0):
        Return the approximate limiting flux density for a pulsar
        survey in mJy based of the following characteristics:
            'Ttot' = sky + system temperature (K)
            'G' = forward gain of the antenna (K/Jy)
            'BW' = observing bandwidth (MHz)
            'T' = integration time (s)
            'P' = pulsar period (s) (default = 0.01)
            'W' = duty cycle of pulsar (0-1) (default = 0.05)
            'polar' = number of polarizations (default = 2)
            'factor' = normalization factor that take into account
                limiting SNR, hardware limitations etc. (default = 15.0)
        Note:  This is a _very_ approximate calculation.  For a better
            calculation, see Cordes and Chernoff, ApJ, 482, p971, App. A.
        Observatories:
            Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy
    """
    w = W * P
    return umath.sqrt(w/((P-w)*polar*BW*T))*factor*Ttot/G

def dm_info(dm=None, dmstep=1.0, freq=1390.0, numchan=512, chanwidth=0.5):
    """
    dm_info(dm=None, dmstep=1.0, freq=1390.0, numchan=512, chanwidth=0.5):
        Return info about potential DM smearing during an observation.
    """
    BW = chanwidth * numchan
    print "      Center freq (MHz) = %.3f" % (freq)
    print "     Number of channels = %d" % (numchan)
    print "    Channel width (MHz) = %.3g" % (chanwidth)
    print "  Total bandwidth (MHz) = %.3g" % (BW)
    print "   DM offset (0.5*step) = %.3g" % (0.5 * dmstep)
    print "  Smearing over BW (ms) = %.3g" % \
          (1000.0 * dm_smear(0.5 * dmstep, BW, freq))
    if (dm):
        print " Smearing per chan (ms) = %.3g" % \
              (1000.0 * dm_smear(dm, chanwidth, freq))

def best_dm_step(maxsmear=0.1, dt=0.00080, dm=0.0, freq=1390.0, numchan=512, chanwidth=0.5):
    """
    best_dm_step(maxsmear=0.1, dt=0.00080, dm=0.0, freq=1390.0, numchan=512, chanwidth=0.5):
        Return the required DM step to keep the total smearing below 'maxsmear' (in ms).
    """
    BW = chanwidth * numchan
    tau_tot = maxsmear/1000.0
    tau_chan = dm_smear(dm, chanwidth, freq)
    tau_samp = dt
    if (tau_tot**2.0 < (tau_chan**2.0+tau_samp**2.0)):
        print "The requested total smearing is smaller than one or more of the components."
        return 0.0
    else:
        return 0.0001205*freq**3.0*2.0/BW*umath.sqrt(tau_tot**2.0-tau_chan**2.0-tau_samp**2.0)

def dm_smear(dm, BW, center_freq):
    """
    dm_smear(dm, BW, center_freq):
        Return the smearing in sec caused by a 'dm' over a bandwidth
        of 'BW' MHz centered at 'center_freq' MHz.
    """
    return dm * BW / (0.0001205 * center_freq * center_freq * center_freq)

def diagonal_DM(dt, chanBW, center_freq):
    """
    diagonal_DM(dt, chanBW, center_freq):
        Return the so-called "diagonal DM" where the smearing across
        one channel is equal to the sample time.
    """
    return (0.0001205 * center_freq * center_freq * center_freq) * dt / chanBW

def guess_DMstep(DM, dt, BW, f_ctr):
    """
    guess_DMstep(DM, dt, BW, f_ctr):
        Choose a reasonable DMstep by setting the maximum smearing across the
        'BW' to equal the sampling time 'dt'.
    """
    return dt*0.0001205*f_ctr**3.0/(0.5*BW)

def delay_from_DM(DM, freq_emitted):
    """
    Return the delay in seconds caused by dispersion, given
    a Dispersion Measure (DM) in cm-3 pc, and the emitted
    frequency (freq_emitted) of the pulsar in MHz.
    """
    if (freq_emitted > 0.0):
        return DM/(0.000241*freq_emitted*freq_emitted)
    else:
        return 0.0

def smear_plot(dm=[1.0,1000.0], dmstep=1.0, subdmstep=10.0, freq=1390.0,
               numchan=512, numsub=32, chanwidth=0.5, dt=0.000125,
               device='/xwin'):
    """
    smear_plot(dm=[0.0,1000.0], dmstep=1.0, subdmstep=10.0, freq=1390.0,
               numchan=512, numsub=32, chanwidth=0.5, dt=0.000125,
               device='/xwin'):
         Show a plot that displays the expected smearing in ms
         from various effects during a radio pulsar search.
    """
    numpts = 500
    BW = numchan * chanwidth
    subBW = numchan / numsub * chanwidth
    maxDMerror = 0.5 * dmstep
    maxsubDMerror = 0.5 * subdmstep
    ldms = span(umath.log10(dm[0]), umath.log10(dm[1]), numpts)
    dms = 10.0**ldms
    # Smearing from sample rate
    dts = Numeric.zeros(numpts) + 1000.0 * dt
    # Smearing due to the intrinsic channel width
    chan_smear = 1000.0 * dm_smear(dms, chanwidth, freq)
    # Smearing across the full BW due to max DM mismatch
    BW_smear = Numeric.zeros(numpts) + \
               1000.0 * dm_smear(maxDMerror, BW, freq)
    # Smearing in each subband due to max DM mismatch
    subband_smear = Numeric.zeros(numpts) + \
                    1000.0 * dm_smear(maxsubDMerror, subBW, freq)
    total_smear = umath.sqrt(dts**2.0 + chan_smear**2.0 +
                             subband_smear**2.0 + BW_smear**2.0)
    maxval = umath.log10(2.0 * max(total_smear))
    minval = umath.log10(0.5 * min([min(dts), min(chan_smear),
                                    min(BW_smear), min(subband_smear)]))
    Pgplot.plotxy(umath.log10(total_smear), ldms, rangey=[minval, maxval],
                  logx=1, logy=1, labx="Dispersion Measure",
                  laby="Smearing (ms)", device=device)
    ppgplot.pgsch(0.8)
    ppgplot.pgmtxt("t", 1.5, 1.0/12.0, 0.5, "\(2156)\dcenter\u = %gMHz" % freq)
    ppgplot.pgmtxt("t", 1.5, 3.0/12.0, 0.5, "N\dchan\u = %d" % numchan)
    ppgplot.pgmtxt("t", 1.5, 5.0/12.0, 0.5, "N\dsub\u = %d" % numsub)
    ppgplot.pgmtxt("t", 1.5, 7.0/12.0, 0.5, "BW\dchan\u = %gMHz" % chanwidth)
    ppgplot.pgmtxt("t", 1.5, 9.0/12.0, 0.5, "\gDDM = %g" % dmstep)
    ppgplot.pgmtxt("t", 1.5, 11.0/12.0, 0.5, "\gDDM\dsub\u = %g" % subdmstep)
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt("b", -7.5, 0.95, 1.0, "Total")
    Pgplot.plotxy(umath.log10(dts), ldms, color="green",
                  logx=1, logy=1)
    ppgplot.pgmtxt("b", -6.0, 0.95, 1.0, "Sample Rate")
    Pgplot.plotxy(umath.log10(chan_smear), ldms, color="purple",
                  logx=1, logy=1)
    ppgplot.pgmtxt("b", -4.5, 0.95, 1.0, "Channel")
    Pgplot.plotxy(umath.log10(BW_smear), ldms, color="red",
                  logx=1, logy=1)
    ppgplot.pgmtxt("b", -3.0, 0.95, 1.0, "Full BW")
    Pgplot.plotxy(umath.log10(subband_smear), ldms, color="blue",
                  logx=1, logy=1)
    ppgplot.pgmtxt("b", -1.5, 0.95, 1.0, "Subband")
    ppgplot.pgsci(1)


def search_sensitivity(Ttot, G, BW, chan, freq, T, dm, ddm, dt, Pmin=0.001,
                        Pmax=1.0, W=0.1, polar=2, factor=15.0, pts=1000):
    """
    (periods, S_min) = search_sensitivity(Ttot, G, BW, chan, freq, T, dm,
             ddm, dt, Pmin=0.001, Pmax=1.0, W=0.1, polar=2, factor=15.0, pts=1000):
        Return the approximate limiting flux density for a pulsar
        survey in mJy based of the following characteristics:
            'Ttot' = sky + system temperature (K)
            'G' = forward gain of the antenna (K/Jy)
            'BW' = observing bandwidth (MHz)
            'chan' = number of channels in the filterbank
            'freq' = central observing frequency (MHz)
            'T' = integration time (s)
            'dm' = Dispersion Measure in pc cm^-3
            'ddm' = Dispersion Measure stepsize in pc cm^-3
            'dt' = Sample time for each data point in sec
            'Pmin' = minimum pulsar period (s) (default = 0.001)
            'Pmax' = maximum pulsar period (s) (default = 1.0)
            'W' = duty cycle of pulsar (0-1) (default = 0.1)
            'polar' = number of polarizations (default = 2)
            'factor' = normalization factor that take into account
                limiting SNR, hardware limitations etc. (default = 15.0)
            'pts' = the number of points to calculate
        Note:  This is a _very_ approximate calculation.  For a better
            calculation, see Cordes and Chernoff, ApJ, 482, p971, App. A.
        Observatories:
            Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy
    """
    periods = span(Pmin, Pmax, pts)
    widths = umath.sqrt((W * periods)**2.0 +
                        dm_smear(dm, BW/chan, freq)**2.0 + \
                        dm_smear(ddm/2.0, BW, freq)**2.0 + \
                        dt**2.0) / periods
    return (periods, limiting_flux_dens(Ttot, G, BW, T, periods, widths, 
                                        polar=polar, factor=factor))

def smin_noise(Ttot, G, BW, dt):
    """
    smin_noise(Ttot, G, BW, dt):
        Return the 1 sigma Gaussian noise level (mJy) for each time
        series bin in a pulsar data simulation.  Default is for a
        sinusoidal pulse (i.e. W = P / 2) with freq << Nyquist freq.
            'Ttot' = sky + system temperature (K)
            'G' = forward gain of the antenna (K/Jy)
            'BW' = observing bandwidth (MHz)
            'dt' = time per time series bin (s)
        Observatories:
            Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy
    """
    return Ttot / (G * umath.sqrt(2 * BW * dt))

def calc_phs(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0): 
    """
    calc_phs(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
        Return the rotational phase (0-1) at MJD given a reference MJD and the
            rotational freq (f0) and optional freq derivs (f1-f5).
    """
    t = (MJD-refMJD)*SECPERDAY
    return umath.fmod(t*(f0 +t*(f1/2.0 +
                                t*(f2/6.0 +
                                   t*(f3/24.0 +
                                      t*(f4/120.0 +
                                         t*f5/720.0))))), 1.0)

def calc_freq(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
    """
    calc_freq(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
        Return the instantaneous frequency at an MJD given a reference
            MJD and the rotational freq (f0) and optional freq derivs (f1-f5).
    """
    t = (MJD-refMJD)*SECPERDAY
    return f0 + t*(f1 +
                   t*(f2/2.0 +
                      t*(f3/6.0 +
                         t*(f4/24.0 +
                            t*f5/120.0))))

def calc_t0(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
    """
    calc_t0(MJD, refMJD, f0, f1=0.0, f2=0.0, f3=0.0, f4=0.0, f5=0.0):
        Return the closest previous MJD corresponding to phase=0 of the pulse.
    """
    phs = calc_phs(MJD, refMJD, f0, f1, f2, f3, f4, f5)
    p = 1.0/calc_freq(MJD, refMJD, f0, f1, f2, f3, f4, f5)
    return MJD-phs*p/SECPERDAY


def write_princeton_toa(toa, toaerr, freq, dm, obs='@', name=' '*13):
    """
    Princeton Format

    columns     item
    1-1     Observatory (one-character code) '@' is barycenter
    2-2     must be blank
    16-24   Observing frequency (MHz)
    25-44   TOA (decimal point must be in column 30 or column 31)
    45-53   TOA uncertainty (microseconds)
    69-78   DM correction (pc cm^-3)
    """
    if dm!=0.0:
        print obs+" %13s %8.3f %19.13f %8.2f              %9.4f" % \
              (name, freq, toa, toaerr, dm)
    else:
        print obs+" %13s %8.3f %19.13f %8.2f" % \
              (name, freq, toa, toaerr)

def rotate(arr, bins):
    """
    rotate(arr, bins):
        Return an array rotated by 'bins' places to the left
    """
    bins = bins % len(arr)
    if bins==0:
        return arr
    else:
        return Numeric.concatenate((arr[bins:], arr[:bins]))

def interp_rotate(arr, bins, zoomfact=10):
    """
    interp_rotate(arr, bins, zoomfact=10):
        Return a sinc-interpolated array rotated by 'bins' places to the left.
            'bins' can be fractional and will be rounded to the closest
            whole-number of interpolated bins.  The resulting vector will
            have the same length as the oiginal.
    """
    newlen = len(arr)*zoomfact
    rotbins = int(umath.floor(bins*zoomfact+0.5)) % newlen
    newarr = sinc_interp.periodic_interp(arr, zoomfact)  
    return rotate(newarr, rotbins)[::zoomfact]

def corr(profile, template):
    """
    corr(profile, template):
        Cross-correlate (using FFTs) a 'profile' and a 'template'.
    """
    return FFT.inverse_real_fft(FFT.real_fft(template) * umath.conjugate(FFT.real_fft(profile)))

def maxphase(profile, template):
    """
    maxphase(profile, template):
        Return the phase offset required to get the 'profile' to best
            match the 'template'.
    """
    return float(Numeric.argmax(corr(profile, template))) / len(profile)

def linear_interpolate(vector, zoom=10):
    """
    linear_interpolate(vector, zoom=10):
        Linearly interpolate 'vector' by a factor of 'zoom'.
    """
    n = len(vector)
    ivect = Numeric.zeros(zoom*n, typecode='d')
    nvect = Numeric.concatenate((vector, vector[:1]))
    ivals = Numeric.arange(zoom, typecode='d')/zoom
    loy = nvect[0]
    for ii in range(n):
        hiy = nvect[ii+1]
        ivect[ii*zoom:(ii+1)*zoom] = ivals*(hiy-loy) + loy
        loy = hiy
    return ivect

def downsample(vector, factor):
    """
    downsample(vector, factor):
        Downsample (i.e. co-add consecutive numbers) a short section
            of a vector by an integer factor.
    """
    if (len(vector) % factor):
        print "Lenght of 'vector' is not divisible by 'factor'=%d!" % factor
        return 0
    newvector = Numeric.reshape(vector, (len(vector)/factor, factor))
    return umath.add.reduce(newvector, 1)

def measure_phase_corr(profile, template, zoom=10):
    """
    measure_phase_corr(profile, template, zoom=10):
        Return the phase offset required to get the 'profile' to best
            match the 'template', each of which has been interpolated
            by a factor of 'zoom'.
    """
    zoomprof = zoomtemp = zoom
    if (len(template) != len(profile)):
        if (len(template)%len(profile) == 0):
            zoomprof = zoom*len(template)/len(profile)
        else:
            print "Warning!:  The lengths of the template (%d) and profile (%d)" % \
                  (len(template), len(profile))
            print "           are not the same!"
    itemp = linear_interpolate(rotate(template, Numeric.argmax(template)), zoomtemp)
    iprof = linear_interpolate(profile, zoomprof)
    return maxphase(iprof, itemp)

def spike_profile(N, phase, fwhm):
    """
    spike_profile(N, phase, fwhm):
        Return a triangular pulse profile with 'N' bins and
        an integrated 'flux' of 1 unit.
            'N' = the number of points in the profile
            'phase' = the pulse phase (0-1)
            'fwhm' = the triangular pulses full width at half-max
    """
    phsval = Numeric.arange(N, typecode='d') / float(N)
    peakst = 0.5 - fwhm
    peakend = 0.5 + fwhm
    normalize = 1.0 / fwhm
    if (mean < 0.5):
        phsval = Numeric.where(umath.greater(phsval, mean+0.5),
                               phsval-1.0, phsval)
    else:
        phsval = Numeric.where(umath.less(phsval, mean-0.5),
                               phsval+1.0, phsval)
    return Numeric.where(umath.less_equal(phsval, 0.5),
                         Numeric.where(umath.less_equal(phsval, peakst),
                                       0.0, (phsval - peakst) *
                                       normalize * normalize),
                         Numeric.where(umath.greater(phsval, peakend),
                                       0.0, (1.0 - (phsval - 0.5) *
                                             normalize) * normalize))

def harm_to_sum(fwhm):
    """
    harm_to_sum(fwhm):
        For an MVMD profile returns the optimal number
            of harmonics to sum incoherently
    """
    fwhms = [0.0108, 0.0110, 0.0113, 0.0117, 0.0119, 0.0124, 0.0127, 0.0132,
             0.0134, 0.0140, 0.0145, 0.0151, 0.0154, 0.0160, 0.0167, 0.0173,
             0.0180, 0.0191, 0.0199, 0.0207, 0.0220, 0.0228, 0.0242, 0.0257,
             0.0273, 0.0295, 0.0313, 0.0338, 0.0366, 0.0396, 0.0437, 0.0482,
             0.0542, 0.0622, 0.0714, 0.0836, 0.1037, 0.1313, 0.1799, 0.2883]
    return len(fwhms)-bisect.bisect(fwhms, fwhm)+1

def expcos_profile(N, phase, fwhm):
    """
    expcos_profile(N, phase, fwhm):
        Return a pulse profile with 'N' bins and an integrated 'flux'
        of 1 unit based on the 'Exponentiated Sinusoid'.
            'N' = the number of points in the profile
            'phase' = the pulse phase (0-1)
            'fwhm' = pulse full width at half-max (0.0 < fwhm <= 0.5)
    """
    from simple_roots import secant
    def fwhm_func(k, fwhm=fwhm):
        if (fwhm < 0.02):
            return umath.arccos(1.0-umath.log(2.0)/k)/PI-fwhm
        else:
            return umath.arccos(umath.log(0.5*(umath.exp(k)+
                                               umath.exp(-k)))/k)/PI-fwhm
    phsval = TWOPI * Numeric.arange(N, typecode='d') / float(N)
    phi = -phase * TWOPI
    if (fwhm >= 0.5):
        return umath.cos(phsval + phi) + 1.0
    elif (fwhm < 0.02):
        # The following is from expanding of iO(x) as x->Infinity.
        k = umath.log(2.0) / (1.0 - umath.cos(PI * fwhm))
        # print "Expansion:  k = %f  FWHM = %f" % (k, fwhm_func(k, 0.0))
        phsval = umath.fmod(phsval + phi, TWOPI)
        phsval = Numeric.where(Numeric.greater(phsval, PI),
                               phsval - TWOPI, phsval)
        denom = ((1 + 1/(8*k) + 9/(128*k*k) + 75/(1024*k**3) + 
                 3675/(32768*k**4) + 59535/(262144*k**5)) / umath.sqrt(TWOPI*k))
        return Numeric.where(umath.greater(umath.fabs(phsval/TWOPI), 3.0*fwhm), 0.0, 
                             umath.exp(k*(umath.cos(phsval)-1.0))/denom)
    else:
        k = secant(fwhm_func, 1e-8, 0.5)
        norm = 1.0 / (i0(k) - umath.exp(-k))
        # print "Full Calc:  k = %f  FWHM = %f" % (k, fwhm_func(k, 0.0))
    if (k < 0.05):
        tmp = umath.cos(phsval + phi)
        tmp2 = tmp * tmp
        return norm * (k * (tmp + 1) +
                       k * k * (tmp2 - 1.0) / 2.0 +
                       k * k * k * (tmp2 * tmp + 1.0) / 6.0)
    else:
        return norm * (umath.exp(k * umath.cos(phsval + phi)) -
                       umath.exp(-k))

def gaussian_profile(N, phase, fwhm):
    """
    gaussian_profile(N, phase, fwhm):
        Return a gaussian pulse profile with 'N' bins and
        an integrated 'flux' of 1 unit.
            'N' = the number of points in the profile
            'phase' = the pulse phase (0-1)
            'fwhm' = the gaussian pulses full width at half-max
        Note:  The FWHM of a gaussian is approx 2.35482 sigma
    """
    sigma = fwhm / 2.35482
    mean = phase % 1.0
    phsval = Numeric.arange(N, typecode='d') / float(N)
    if (mean < 0.5):
        phsval = Numeric.where(umath.greater(phsval, mean+0.5),
                               phsval-1.0, phsval)
    else:
        phsval = Numeric.where(umath.less(phsval, mean-0.5),
                               phsval+1.0, phsval)
    try:
        zs = (phsval-mean)/sigma
        okzinds = Numeric.compress(umath.fabs(zs)<20.0, Numeric.arange(N))
        okzs = Numeric.take(zs, okzinds)
        retval = Numeric.zeros(N, 'd')
        Numeric.put(retval, okzinds, umath.exp(-0.5*(okzs)**2.0)/(sigma*umath.sqrt(2*PI)))
        return retval
    except OverflowError:
        print "Problem in gaussian prof:  mean = %f  sigma = %f" % \
              (mean, sigma)
        return Numeric.zeros(N, 'd')
        
def gauss_profile_params(profile, output=0):
    """
    gauss_profile_params(profile, output=0):
        Return parameters of a best-fit gaussian to a profile.
        The funtion returns a tuple containg the following values:
           ret[0] = Best-fit gaussian integrated 'flux'.
           ret[1] = Best-fit gaussian FWHM.
           ret[2] = Best-fit gaussian phase (0.0-1.0).
           ret[3] = Baseline (i.e. noise) average value.
           ret[4] = Residuals average value.
           ret[5] = Residuals standard deviation.
        If 'output' is true, the fit will be plotted and 
           the return values will be printed.
    """
    def funct(afpo, profile):
        return afpo[0] * gaussian_profile(len(profile), afpo[2], afpo[1]) \
               + afpo[3] - profile
    ret = leastsq(funct, [max(profile)-min(profile),
                          0.25, Numeric.argmax(profile)/float(len(profile)),
                          min(profile)], args=(profile))
    if (output):
        phases = Numeric.arange(0.0, 1.0,
                                1.0 / len(profile)) + 0.5 / len(profile)
        Pgplot.plotxy(profile, phases, rangex=[0.0, 1.0],
                      labx='Pulse Phase', laby='Pulse Intensity')
    bestfit = ret[0][0] * gaussian_profile(len(profile),
                                           ret[0][2], ret[0][1]) \
              + ret[0][3]
    if (output):
        Pgplot.plotxy(bestfit, phases, color='red')
        Pgplot.closeplot()
    residuals = bestfit - profile
    resid_avg = average(residuals)
    resid_std = standardDeviation(residuals)
    if (output):
        Pgplot.plotxy(residuals, phases, rangex=[0.0, 1.0],
                      rangey=[min(residuals) - 2 * resid_std,
                              max(residuals) + 2 * resid_std],
                      labx='Pulse Phase', laby='Residuals',
                      line=None, symbol=3)
        ppgplot.pgerrb(6, phases, residuals,
                       Numeric.zeros(len(residuals), 'd') + \
                       resid_std, 2)
        Pgplot.plotxy([resid_avg, resid_avg], [0.0, 1.0], line=2)
        Pgplot.closeplot()
        print ""
        print "  Best-fit gaussian integrated 'flux'  = ", ret[0][0]
        print "               Best-fit gaussian FWHM  = ", ret[0][1]
        print "    Best-fit gaussian phase (0.0-1.0)  = ", ret[0][2]
        print "        Baseline (i.e. noise) average  = ", ret[0][3]
        print "                    Residuals average  = ", resid_avg
        print "         Residuals standard deviation  = ", resid_std
        print ""
    return (ret[0][0], ret[0][1], ret[0][2], ret[0][3], resid_avg, resid_std)

def estimate_flux_density(profile, N, dt, Ttot, G, BW, prof_stdev, display=0):
    """
    estimate_flux_density(profile, N, dt, Ttot, G, BW, prof_stdev, display=0):
        Return an estimate of the flux density (mJy) of a pulsar.
            'profile' = the pulse profile you are using
            'N' = number of time series bins folded
            'dt' = time per time series bin (s)
            'Ttot' = sky + system temperature (K)
            'G' = forward gain of the antenna (K/Jy)
            'BW' = observing bandwidth (MHz)
            'prof_stdev' = profile standard deviation
            'display' = if set, the gaussian fit plots are shown
        Observatories:
            Parkes Multibeam: Tsys = 21 K, G = 0.735 K/Jy
    """
    (amp, fwhm, phase, offset, resid_avg, resid_std) = \
          gauss_profile_params(profile, display)
    T = N * dt
    norm_fact = (prof_stdev * len(profile)) / \
                 smin_noise(Ttot, G, BW, T / len(profile))
    return umath.add.reduce(profile - offset) / norm_fact

def max_spike_power(FWHM):
    """
    max_spike_power(FWHM):
        Return the (approx.) ratio of the highest power from a
        triangular spike pulse profile to the power from a
        perfect sinusoidal pulse profile.  In other words, if a
        sine gives you a power of 1, what power does a spike profile
        give you?  Both the spike and the sine are assumed to have
        an area under one full pulse of 1 unit.  Note:  A gaussian
        profile gives almost identical powers as a spike profile
        of the same width.  This expression was determined using
        a least-squares fit (Max abs error ~ 0.016).
            'FWHM' is the full width at half-max of the spike.
                (0.0 < FWHM <= 0.5)
    """
    return ((36.4165309504 * FWHM - 32.0107844537) * FWHM \
           + 0.239948319674) * FWHM + 4.00277916584

def num_spike_powers(FWHM):
    """
    num_spike_powers(FWHM):
        Return the (approx.) number of powers from a triangular spike
        pulse profile which are greater than one half the power
        perfect sinusoidal pulse profile.  Both the spike and the
        sine are assumed to have an area under one full pulse of 1 unit.
        Note:  A gaussian profile gives almost identical numbers of
        high powers as a spike profile of the same width.  This
        expression was determined using a least-squares fit.
        (Errors get large as FWHM -> 0).
            'FWHM' is the full width at half-max of the spike.
                (0.0 < FWHM <= 0.5)
    """
    return -3.95499721563e-05 / FWHM**2 + 0.562069634689 / FWHM - \
           0.683604041138

def incoherent_sum(amps):
    """
    incoherent_sum(amps):
        Given a series of complex Fourier amplitudes, return a vector
            showing the accumulated incoherently-summed powers.
    """
    return umath.add.accumulate(umath.absolute(amps)**2.0)

def coherent_sum(amps):
    """
    coherent_sum(amps):
        Given a series of complex Fourier amplitudes, return a vector
            showing the accumulated coherently-summed powers.
    """
    phss = umath.arctan2(amps.imag, amps.real)
    phs0 = phss[0]
    phscorr = phs0 - umath.fmod((Num.arange(len(amps), typecode='d')+1.0)*phs0, TWOPI)
    sumamps = umath.add.accumulate(amps*umath.exp(complex(0.0, 1.0)*phscorr))
    return umath.absolute(sumamps)**2.0

def prob_power(power):
    """
    prob_power(power):
        Return the probability for noise to exceed a normalized power
        level of 'power' in a power spectrum.
    """
    return umath.exp(-power)

def prob_sum_powers(power, nsum):
    """
    prob_sum_powers(power, nsum):
        Return the probability for noise to exceed 'power' in
        the sum of 'nsum' normalized powers from a power spectrum.
    """
    return chdtrc(2*nsum, 2.0*power)

def sigma_power(power):
    """
    sigma_power(power):
        Return the approximate equivalent Gaussian sigma for noise
        to exceed a normalized power level given as 'power'
        in a power spectrum.
    """
    if power > 36.0: return umath.sqrt(2.0 * power -
                                       umath.log(PI * power))
    else: return ndtri(1.0 - prob_power(power))

def sigma_sum_powers(power, nsum):
    """
    sigma_sum_powers(power, nsum):
        Return the approximate equivalent Gaussian sigma for noise
        to exceed a sum of 'nsum' normalized powers given by 'power'
        in a power spectrum.
    """
    return ndtri(1.0 - prob_sum_powers(power, nsum))

def power_at_sigma(sigma):
    """
    power_at_sigma(sigma):
        Return the approximate normalized power level that is
        equivalent to a detection of significance 'sigma'.
    """
    return sigma**2 / 2.0 + umath.log(umath.sqrt(PIBYTWO)
                                      * sigma)
def powersum_at_sigma(sigma, nsum):
    """
    powersum_at_sigma(sigma, nsum):
        Return the approximate sum of 'nsum' normalized powers that is
        equivalent to a detection of significance 'sigma'.
    """
    return 0.5 * chdtri(2.0 * nsum, 1.0 - ndtr(sigma))

def cand_sigma(N, power):
    """
    cand_sigma(N, power):
        Return the sigma of a candidate found in a power
        spectrum of 'N' bins after taking into account the
        number of bins searched.
    """
    return ndtri(1.0 - N * prob_power(power))

def fft_max_pulsed_frac(N, numphot, sigma=3.0):
    """
    fft_max_pulsed_frac(N, numphot, sigma=3.0):
        Return the approximate maximum pulsed fraction for a
        sinusoidal signal that _wasn't_ found in a FFT-based
        search.  'N' is the number of bins searched in the FFT.
        'numphot' is the number of photons present in the data.
        And 'sigma' is your confidence (in sigma) that you
        have in expressing this limit.
    """
    # The following is the power level required to get a
    # noise spike that would appear somewhere in N bins
    # at the 'sigma' level
    power_required = -umath.log((1.0-ndtr(sigma))/N)
    return umath.sqrt(4.0 * numphot * power_required)/N

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

def pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
    """
    pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
       Calculate the period or frequency errors and
       the pdot or fdot errors from the opposite one.
    """
    if (pdorfd==None):
        return [1.0 / porf, porferr / porf**2.0]
    else:
        forperr = porferr / porf**2.0
        fdorpderr = umath.sqrt((4.0 * pdorfd**2.0 * porferr**2.0) / porf**6.0 +
                               pdorfderr**2.0 / porf**4.0)
        [forp, fdorpd] = p_to_f(porf, pdorfd)
        return [forp, forperr, fdorpd, fdorpderr]

def psr_info(porf, pdorfd, time=None, input=None):
    """
    psr_info(porf, pdorfd, input=None):
        Print a list of standard derived pulsar parameters based
        on the period (or frequency) and its first derivative.  The
        routine will automaticallu assume you are using periods if
        'porf' <= 1.0 and frequencies otherwise.  You can override this
        by setting input='p' or 'f' appropriately.
    """
    if ((input==None and porf > 1.0) or
        (input=='f' or input=='F')):
        pdorfd = - pdorfd / (porf * porf)
        porf = 1.0 / porf
    I = 1.0e45  # Moment of Inertia in g cm^2
    Edot = 4.0 * PI * PI * I * pdorfd / porf ** 3.0
    Bo = 3.2e19 * umath.sqrt(porf * pdorfd)
    age = porf / (2.0 * pdorfd * 31557600.0) 
    [f, fd] = p_to_f(porf, pdorfd)
    print ""
    print "             Period = %f s" % porf
    print "              P-dot = %g s/s" % pdorfd
    print "          Frequency = %f Hz" % f
    print "              F-dot = %g Hz/s" % fd
    if (time):
        print "       Fourier Freq = %g bins" % (f * time)
        print "      Fourier F-dot = %g bins" % (fd * time * time)
    print "              E-dot = %g ergs/s" % Edot
    print "    Surface B Field = %g gauss" % Bo
    print " Characteristic Age = %g years" % age
    print "          Assumed I = %g g cm^2" % I
    print ""
