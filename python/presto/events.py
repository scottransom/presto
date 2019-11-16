from __future__ import print_function
import bisect
from presto.psr_constants import PI, TWOPI, PIBYTWO
from presto.simple_roots import newton_raphson
from scipy.special import iv, chdtri, ndtr, ndtri
from presto.cosine_rand import *
import numpy as np


def sine_events(pulsed_frac, Nevents, phase=0.0):
    """
    sine_events(pulsed_frac, Nevents, phase=0.0):
       Return an array of 'Nevents' of phase values [0,1)
       simulating a folded profile with a pulsed fraction
       'pulsed_frac', a phase offset 'phase', and with a
       sinusoidal pulse profile.
    """
    Nsrc = int(pulsed_frac*Nevents+0.5)
    Nbak = Nevents - Nsrc
    phases = Num.zeros(Nevents, dtype=Num.float)
    phases[:Nsrc] += cosine_rand(Nsrc) + phase
    phases[Nsrc:] += Num.random.random(Nbak)
    phases = Num.fmod(phases, 1.0)
    phases[phases<0.0] += 1.0
    return phases

def gaussian_events(pulsed_frac, Nevents, fwhm, phase=0.0):
    """
    gaussian_events(pulsed_frac, Nevents, phase=0.0):
       Return an array of 'Nevents' of phase values [0,1)
       simulating a folded profile with a pulsed fraction
       'pulsed_frac', a phase offset 'phase', and with a
       gaussian pulse profile of width 'fwhm'
    """
    sigma = fwhm / 2.35482
    Nsrc = int(pulsed_frac*Nevents+0.5)
    Nbak = Nevents - Nsrc
    phases = Num.zeros(Nevents, dtype=Num.float)
    phases[:Nsrc] += Num.random.standard_normal(Nsrc)*sigma + phase
    phases[Nsrc:] += Num.random.random(Nbak)
    phases = Num.fmod(phases, 1.0)
    phases[phases<0.0] += 1.0
    return phases

def harm_to_sum(fwhm):
    """
    harm_to_sum(fwhm):
       For an MVMD profile of width 'fwhm', returns the
       optimal number of harmonics to sum incoherently
    """
    fwhms = [0.0108, 0.0110, 0.0113, 0.0117, 0.0119, 0.0124, 0.0127, 0.0132,
             0.0134, 0.0140, 0.0145, 0.0151, 0.0154, 0.0160, 0.0167, 0.0173,
             0.0180, 0.0191, 0.0199, 0.0207, 0.0220, 0.0228, 0.0242, 0.0257,
             0.0273, 0.0295, 0.0313, 0.0338, 0.0366, 0.0396, 0.0437, 0.0482,
             0.0542, 0.0622, 0.0714, 0.0836, 0.1037, 0.1313, 0.1799, 0.2883]
    return len(fwhms)-bisect.bisect(fwhms, fwhm)+1

def DFTexact(times, f, maxnumharms=20):
    """
    DFTexact(times, f, maxnumharms=20):
       Return an array of 'maxnumharms' complex amplitudes
       corresponding to the harmonics of the 'times' (in sec)
       with a fundamental at frequency 'f' Hz.
    """
    const = -TWOPI*(Num.arange(maxnumharms, dtype=Num.float)+1.0)*f*complex(0.0, 1.0)
    return Num.add.reduce(Num.exp(Num.outerproduct(const,times)), axis=1)

def incoherent_sum(amps):
    """
    incoherent_sum(amps):
       Return the incoherent sum of an array of complex Fourier
       amplitudes.  Usually these correspond to the complex
       harmonics of a periodic signal.
    """
    return Num.add.accumulate(Num.abs(amps)**2.0)

def coherent_sum(amps):
    """
    coherent_sum(amps):
       Return the coherent sum (i.e. including phase information)
       of an array of complex Fourier amplitudes.  Usually these
       correspond to the complex harmonics of a periodic signal.
    """
    phss = Num.arctan2(amps.imag, amps.real)
    phs0 = phss[0]
    phscorr = phs0 - Num.fmod(Num.arange(1.0, len(amps)+1,
                                         dtype=Num.float)*phs0, TWOPI)
    sumamps = Num.add.accumulate(amps*Num.exp(complex(0.0, 1.0)*phscorr))
    return Num.abs(sumamps)**2.0

def Htest_exact(phases, maxnumharms=20, weights=None):
    """
    Htest_exact(phases, maxnumharms=20, weights=None):
       Return an exactly computed (i.e. unbinned) H-test statistic
       for periodicity for the events with folded phases 'phases' [0,1).
       Also return the best number of harmonics.  The H-statistic and
       harmonic number are returned as a tuple: (hstat, harmnum).
       This routine returns the Leahy normalized H-statistic, and the
       best number of harmonics summed.  If weights are set to be
       fractional photon weights, then the weighted Htest is returned
       (see Kerr 2011: http://arxiv.org/pdf/1103.2128.pdf)
    """
    N = len(phases)
    Zm2s = np.zeros(maxnumharms, dtype=np.float)
    rad_phases = 2.0*np.pi*phases
    weightfact = 1.0/(np.sum(weights**2.0) / N) if \
                 weights is not None else 1.0
    for harmnum in range(1, maxnumharms+1):
        phss = harmnum*rad_phases
        Zm2s[harmnum-1] = 2.0/N*(np.add.reduce(np.sin(phss))**2.0+
                                 np.add.reduce(np.cos(phss))**2.0)
        Zm2s[harmnum-1] *= weightfact
    hs = np.add.accumulate(Zm2s) - \
         4.0*np.arange(1.0, maxnumharms+1)+4.0
    bestharm = hs.argmax()
    return (hs[bestharm], bestharm+1)

def Hstat_prob(h):
    """
    Hstat_prob(h):
       Return the probability associated with an H-test statistic
       of value 'h'.  Uses de Jager & Busching 2010 result.
    """
    return Num.exp(-0.4 * h)

def gauss_sigma_to_prob(sigma):
    """
    gauss_sigma_to_prob(sigma):
        Returns the area under the Gaussian probability density
        function, integrated from 'sigma' to infinity.
    """
    if sigma < 5.0:
        return 1.0 - ndtr(sigma)
    else:
        # From A&S page 932, eqn 26.2.12 for Q(x)
        x = sigma
        Z = 1.0/Num.sqrt(2.0*Num.pi) * Num.exp(-0.5*x*x)
        series = Num.sum(Num.asarray([1.0, -1.0/(x*x), 3.0/(x**4.0),
                                      -15.0/(x**6.0), 105.0/(x**8.0)]))
        return Z/x*series

def prob_to_gauss_sigma(prob):
    """
    prob_to_gauss_sigma(prob):
        Returns the Gaussian sigma for which the area under the
        Gaussian probability density function (integrated from minus
        infinity to 'sigma') is equal to 'prob'.
    """
    return ndtri(prob)

def xray_time_to_detect(ctrate, pfract, dt, fpsr, bins=0, confidence=0.99,
                        detectfract=0.99):
    """
    xray_time_to_detect(ctrate, pfract, dt, fpsr, bins=0, confidence=0.99,
                        detectfract=0.99):
        Return the observation duration required (assuming no breaks
        and a sinusoidal pulse profile) to detect pulsations at
        frequency 'fpsr' while looking in a number of Fourier
        bins equal to 'bins' (Note: the default value of 0 means
        that all bins will be examined).  'dt' is the bin duration in
        sec, 'ctrate' is the total expected count rate, and 'pfract' is
        the expected pulsed fraction.  'confidence' is the confidence
        level that the signal is not caused by noise, and 'detectfract'
        is the fraction of the time that you want this observation to
        occur (i.e. if set to 0.5, 50% of observations of this duration
        would detect the specified signal at 'confidence' level).
    """
    nyquist_freq = 0.5 / dt
    factor = binning_factor(fpsr, nyquist_freq)**2.0
    A = pfract * ctrate  # Signal ct rate
    if (bins):
        P_detect = max_noise_power(bins, confidence=confidence)
        power_required = required_signal_power(P_detect, confidence=detectfract)
	# The following is from para 1, sect 3.3, of Ransom, Gaensler, and Slane, 2002
	#return  (power_required - 1.0)          4 * ctrate * dt**2.0 / (A**2.0 * factor) * 
    else:
        print("Not implemented yet...I think we need to iterate.")

# The following routines are based on the method of signal
# estimation described by Vaughan et al., 1994, ApJ, 435, p362.
# The math comes from Groth, 1975, ApJS, 29, p285.

def power_average(signal_power, n=1):
    """
    power_average(signal_power, n=1):
        Return the expectation value of the measured power given
        a signal with intrinsic power 'signal_power' and 'n'
        summed powers.  This is from equation 14 in Groth, 1975.
    """
    return signal_power + n

def power_variance(signal_power, n=1):
    """
    power_variance(signal_power, n=1):
        Return the variance of the measured power given a signal
        with intrinsic power 'signal_power' and 'n' summed
        powers.  This is from equation 14 in Groth, 1975.
    """
    return 2.0 * signal_power + n

def power_sigma(signal_power, n=1):
    """
    power_sigma(signal_power, n=1):
        Return the standard deviation of the measured power
        given a signal with intrinsic power 'signal_power' and
        'n' summed powers.  This is from equation 14 in Groth, 1975.
    """
    return Num.sqrt(power_variance(signal_power, n))

def log_fact_table(maxn):
    """
    log_fact_table(maxn):
        Return a table of the natural logarithms of the
        first 'maxn'+1 factorials.
    """
    table = Num.arange(maxn+1, dtype='d')
    table[0] = 1.0
    return Num.add.accumulate(Num.log(table))

def binning_factor(freq, nyquist_freq):
    """
    binning_factor(freq, nyquist_freq):
        Return the factor that causes high frequency Fourier
        Amplitudes to be decreased if the time series is
        made of binned events.  Square this for a power
        spectrum adjustment.  'freq' is the frequency of
        interest and 'nyquist_freq' is the Nyquist Frequency
        which can be defined as N/(2*T).
    """
    x = 0.5 * Num.asarray(freq) / nyquist_freq
    return Num.sinc(x)  # numpy sinc is defined with pi

def max_noise_power(bins, n=1, confidence=0.99):
    """
    max_noise_power(bins, n=1, confidence=0.99):
        Return the power level that gives you some
        'confidence' that spectral noise could not cause
        that level in your power spectrum.  The total number
        of independent frequencies searched is 'bins'.
        This is P_detect in Vaughan et. al, 1994, and is also
        known as P_threshold.
    """
    if (n==1):
        return -Num.log((1.0 - confidence) / bins)
    else:
        return 0.5 * chdtri(2.0 * n, (1.0 - confidence) / bins)

def prob_power_series(power, signal_power, n=1, TOL=1.0e-14):
    """
    prob_power_series(power, signal_power, n=1, TOL=1.0e-14):
        Return the integrated probability from P=0 to 'power'
        that a signal with theoretical power 'signal_power'
        will show up in a power spectrum with power 'power'.
        This method evaluates the integral using an infinite
        sum and is equation 16 in Groth, 1975.
    """
    fact = Num.exp(-(power + signal_power))
    lf = log_fact_table((power + signal_power) * 5)
    lp, lps = Num.log(power), Num.log(signal_power)
    sum = 0.0
    term = 1.0
    m = 0
    while (1):
        kmax = m + n
        term = fact * Num.add.reduce(Num.exp((Num.arange(kmax)*lp + m*lps) - \
                                             (lf[0:kmax] + lf[m])))
        sum = sum + term
        if (m > signal_power and term < TOL):  break
        m = m + 1
    return 1.0 - sum

def prob_power_integral(power, signal_power, n=1):
    """
    prob_power_integral(power, signal_power, n=1):
        Return the integrated probability from P=0 to 'power'
        that a signal with theoretical power 'signal_power'
        will show up in a power spectrum with power 'power'.
        This method evaluates the integral numerically and
        is equation 18 in Groth, 1975.
    """
    def integrand(theta, p, ps, n):
        t1 = 2 * n * theta
        t2 = Num.sin(2.0 * theta)
        A = t1 + ps * t2
        B = t1 + (ps - p) * t2
        sintheta = Num.sin(theta)
        sin2theta = sintheta**2.0
        return (Num.exp(-2.0 * ps * sin2theta) *
                (Num.sin(A - theta) - Num.exp(-2.0 * p * sin2theta) *
                 Num.sin(B - theta)) / sintheta)
    (val, err) = quad(integrand, 0.0, PIBYTWO, (power, signal_power, n))
    return val/PI

def power_probability(power, signal_power, n=1):
    """
    power_probability(power, signal_power, n=1):
        Return the probability of a signal with power
        'signal_power' actually showing up with power
        'power' in a power spectrum'  This is equation
        12 in Groth, 1975 and is the integrand of the
        prob_power_* functions (which integrate it from 0 to P)
    """
    return (power / signal_power)**(0.5 * (n - 1)) * \
           Num.exp(-(power + signal_power)) * \
           iv(n - 1.0, 2 * Num.sqrt(power * signal_power))

def required_signal_power(power, n=1, confidence=0.99):
    """
    required_signal_power(power, n=1, confidence=0.99):
        Return the required power of a signal that will cause
        at least a power 'power' in a power spectrum a fraction
        'confidence' of the time.  This is the inverse of
        equation 16 in Groth, 1975, with solves for P_signal.
        If called with 'power' = P_detect the result is
        the search sensitivity.  If called with 'power' = P_max,
        then the result is the upper limit on the signal power
        in the power spectrum.
    """
    prob = 1.0 - confidence
    def func(x, power=power, prob=prob, n=n):
        return prob_power_series(power, x, n) - prob
    def dfunc(x, power=power, n=n):
        return power_probability(power, x, n)
    P_signal = newton_raphson(func, dfunc, 0.0001, 100.0)
    return P_signal

def fft_sensitivity(N, bins=0, n=1, confidence=0.99):
    """
    fft_sensitivity(N, bins=0, n=1, confidence=0.99):
        Return a measure of the weakest signal power you can
        confidently detect in an FFT search containing 'N' data
        points (this is the number of bins in the time series -- the
        number of Frequency bins searched is usually N/2).  'bins' is
        only different from 0 if the number of independent frequencies
        searched does not equal N/2 (i.e. when an acceleration search
        is performed).  'confidence' is our fractional confidence in
        the result (i.e. 0.99 = 99% limit).  This calculation does not
        include the correction to sensitivity due to binning effects.
        These calculations are based on the Vaughan et al 1994 paper
        and compute P_sens.
    """
    if not (bins): bins = N / 2
    P_threshold = max_noise_power(bins, n, confidence)
    return required_signal_power(P_threshold, n, confidence)

def rzw_sensitivity(N, zlo=-100.0, zhi=100.0, n=1, confidence=0.99):
    """
    rzw_sensitivity(N, zlo=-100.0, zhi=100.0, n=1, confidence=0.99):
        Return a measure of the weakest signal power you can
        confidently detect in an RZW (Fourier acceleration) search
        containing 'N' data points (this is the number of bins in the
        time series) and low and high acceleration values of 'zlo'
        and 'zhi'.  'confidence' is our fractional confidence in
        the result (i.e. 0.99 = 99% limit).  This calculation does not
        include the correction to sensitivity due to binning effects.
        These calculations are based on the Vaughan et al 1994 paper
        and compute P_sens.
    """
    bins = N / 2.0 * (zhi - zlo + 1.0) / 6.95
    P_threshold = max_noise_power(bins, n, confidence)
    return required_signal_power(P_threshold, n, confidence)

def binned_fft_sensitivity(N, dt, freq, bins=0, n=1, confidence=0.99):
    """
    binned_fft_sensitivity(N, dt, freq, bins=0, n=1, confidence=0.99):
        Return a measure of the weakest signal power of frequency 'freq'
        Hz you can confidently detect in an FFT search containing 'N'
        data points (this is the number of bins in the time series --
        the number of Frequency bins searched is usually 1/2 of this
        value) each of which was binned into 'dt' sec bins.
        'bins' is only different from 0 if the number of independent
        frequencies searched does not equal N/2 (i.e. when an
        acceleration search is performed).  'confidence' is our
        fractional confidence in the result (i.e. 0.99 = 99% limit).
        This calculation includes the correction to sensitivity
        due to binning effects.  These calculations are based on
        the Vaughan et al 1994 paper and compute P_sens.
    """
    nyquist_freq = 0.5 / dt
    factor = binning_factor(freq, nyquist_freq)**2.0
    return fft_sensitivity(N, bins, n, confidence) / factor

def binned_rzw_sensitivity(N, dt, freq, zlo=-100.0, zhi=100.0,
                           n=1, confidence=0.99):
    """
    binned_rzw_sensitivity(N, dt, freq, zlo=-100.0, zhi=100.0,
                           n=1, confidence=0.99):
        Return a measure of the weakest signal power of frequency 'freq'
        Hz you can confidently detect in an RZW (Fourier acceleration)
        search containing 'N' data points (this is the number of bins in
        the time series) each of which was binned into 'dt' sec bins.
        Low and high acceleration values of 'zlo' and 'zhi' were used.
        'confidence' is our fractional confidence in the result (i.e.
        0.99 = 99% limit).  This calculation includes the correction to
        sensitivity due to binning effects.  These calculations are
        based on the Vaughan et al 1994 paper and compute P_sens.
    """
    bins = N / 2.0 * (zhi - zlo + 1.0) / 6.95
    nyquist_freq = 0.5 / dt
    factor = binning_factor(freq, nyquist_freq)**2.0
    return fft_sensitivity(N, bins, n, confidence) / factor

def pulsed_fraction_limit(Nphot, Pow):
    """
    pulsed_fraction_limit(phot, Pow):
        Return an _observational_ (i.e. not intrinsic) upper limit
        to the pulsed fraction of a signal that is in the data but
        was not detected.  By observational, I mean that some of the
        unpulsed events do not come from the source you are looking
        for pulsations in.  The data contain a total of 'Nphot'
        photons and the largest measured power (or P_sens as
        calculated using the *_sensitivity functions in this module)
        is 'Pow'.  If you want the _intrinsic_ pulsed fraction,
        you should divide the returned value by the fraction of Nphot
        that actually comes from the _source_ (i.e. the NS).
    """
    return Num.sqrt(4.0 * (Pow - 1.0) / Nphot)


if __name__=="__main__":
    from presto.psr_utils import *
    from presto.Pgplot import *
    from presto.presto import *
    from RandomArray import *

    prof = expcos_profile(128, 0.0, 0.1) + normal(0.0, 5.0, 128)
    plotxy(prof)
    closeplot()
    fprof = rfft(prof)
    fprof = fprof/Num.sqrt(fprof[0].real)
    pows = spectralpower(fprof)
    tcsum =  Num.add.accumulate(Num.sqrt(pows[1:10]))**2.0
    csum = coherent_sum(fprof[1:10])
    isum = incoherent_sum(fprof[1:10])
    print(isum)
    print(csum)
    print(tcsum)
    for ii in range(len(csum)):
        print(candidate_sigma(isum[ii], ii+1, 1), candidate_sigma(csum[ii]/(ii+1), 1, 1))
        
