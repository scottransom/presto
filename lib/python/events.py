import bisect
import umath
import Numeric as Num
from psr_constants import PI, TWOPI, PIBYTWO
from simple_roots import newton_raphson
from scipy.special.cephes import iv, chdtri, ndtr, ndtri

def gauss_sigma_to_prob(sigma):
    """
    gauss_sigma_to_prob(sigma):
        Returns the area under the Gaussian probability density
        function, integrated from minus infinity to 'sigma'.
    """
    return ndtr(sigma)

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
        return 4 * ctrate * dt**2.0 / (A**2.0 * factor) * power_required
    else:
        print "Not implemented yet...I think we need to iterate."

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
    return umath.sqrt(power_variance(signal_power, n))

def log_fact_table(maxn):
    """
    log_fact_table(maxn):
        Return a table of the natural logarithms of the
        first 'maxn'+1 factorials.
    """
    table = Num.arange(maxn+1, typecode='d')
    table[0] = 1.0
    return Num.add.accumulate(umath.log(table))

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
    x = PIBYTWO * freq / nyquist_freq
    if (x == 0.0):
        return 1.0
    else:
        return (umath.sin(x) / x)

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
        return -umath.log((1.0 - confidence) / bins)
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
    fact = umath.exp(-(power + signal_power))
    lf = log_fact_table((power + signal_power) * 5)
    lp, lps = umath.log(power), umath.log(signal_power)
    sum = 0.0
    term = 1.0
    m = 0
    while (1):
        kmax = m + n
        term = fact * Num.add.reduce(umath.exp((Num.arange(kmax)*lp + m*lps) - \
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
        t2 = umath.sin(2.0 * theta)
        A = t1 + ps * t2
        B = t1 + (ps - p) * t2
        sintheta = umath.sin(theta)
        sin2theta = sintheta**2.0
        return (umath.exp(-2.0 * ps * sin2theta) *
                (umath.sin(A - theta) - umath.exp(-2.0 * p * sin2theta) *
                 umath.sin(B - theta)) / sintheta)
    (val, err) = quad(integrand, 0.0, PIBYTWO, (power, signal_power, n))
    return val/PI

def power_probability(power, signal_power, n=1):
    """
    power_probability(power, signal_power, n=1):
        Return the probability of a signal with power
        'signal_power' actually showing up with power
        'power' in a power spectrum'  This is equation
        12 in Groth, 1975 and is the integrand of the previous
        two functions (which integrate it from 0 to P)
    """
    return (power / signal_power)**(0.5 * (n - 1)) * \
           umath.exp(-(power + signal_power)) * \
           iv(n - 1.0, 2 * umath.sqrt(power * signal_power))

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
        due to binning effects.
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
        sensitivity due to binning effects.
    """
    bins = N / 2.0 * (zhi - zlo + 1.0) / 6.95
    nyquist_freq = 0.5 / dt
    factor = binning_factor(freq, nyquist_freq)**2.0
    return fft_sensitivity(N, bins, n, confidence) / factor

def pulsed_fraction_limit(numphot, power_limit):
    """
    pulsed_fraction_limit(numphot, power_limit):
        Return an upper limit to the pulsed fraction of a signal
        that is in the data but was not detected.  The data
        contain a total of 'numphot' photons and the largest
        measured power is 'power_limit'.
    """
    return umath.sqrt(4.0 * (power_limit - 1.0) / numphot)

def exact_DFT(times, f, maxnumharms=20):
    """
    exact_DFT(times, f, maxnumharms=20):
        Return an array of maxnumharms complex amplitudes
            corresponding to the the harmonics of the times
            (in sec) with a fundamental at frequency f
    """
    const = -TWOPI*(Num.arange(maxnumharms, typecode='d')+1.0)*f*complex(0.0, 1.0)
    return umath.add.reduce(umath.exp(Num.outerproduct(const,times)), axis=1)

def exact_H_test(phases, maxnumharms=20):
    """
    exact_H_test(phases, maxnumharms=20):
        Return the value of 'h' and the corresponding harmonic
            after calculating an 'exact' H-test on the phases (0-1).
    """
    N = len(phases)
    Zm2s = Num.zeros(maxnumharms, 'd')
    rad_phases = TWOPI*phases
    for harmnum in range(1, maxnumharms+1):
        phss = harmnum*rad_phases
        Zm2s[harmnum-1] = 2.0/N * (umath.add.reduce(umath.sin(phss))**2.0 +
                                   umath.add.reduce(umath.cos(phss))**2.0)
    hs = umath.add.accumulate(Zm2s) - \
         4.0*Num.arange(1, maxnumharms+1, typecode='d') + 4.0
    bestharm = Num.argmax(hs)
    return (hs[bestharm]/2.0, bestharm+1)

def H_prob(h):
    """
    H_prob(h):
        Return the probability of getting a vaule 'h' or greater
            from the H-test.
    """
    h *= 2.0
    if (h <= 23.0):
        return 1.210597*umath.exp(-0.45901*h+0.0022900*h*h)
    elif (h < 50.0):
        return 0.9999755*umath.exp(-0.39802*h)
    else:
        return 0.0

if __name__=="__main__":
    from psr_utils import *
    from Pgplot import *
    from presto import *
    from RandomArray import *
    
    prof = expcos_profile(128, 0.0, 0.1) + normal(0.0, 5.0, 128)
    plotxy(prof)
    closeplot()
    fprof = rfft(prof)
    fprof = fprof/umath.sqrt(fprof[0].real)
    pows = spectralpower(fprof)
    tcsum =  umath.add.accumulate(umath.sqrt(pows[1:10]))**2.0
    csum = coherent_sum(fprof[1:10])
    isum = incoherent_sum(fprof[1:10])
    print isum
    print csum
    print tcsum
    for ii in range(len(csum)):
        print candidate_sigma(isum[ii], ii+1, 1), candidate_sigma(csum[ii]/(ii+1), 1, 1)
        
