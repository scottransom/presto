#!/usr/bin/env python
import Numeric
from umath import pi, log, exp, sin, sqrt
from simple_roots import newton_raphson
from cephes import iv
#from Multipack import quad

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
    return sqrt(power_variance(signal_power, n))

def log_fact_table(maxn):
    """
    log_fact_table(maxn):
        Return a table of the natural logarithms of the
        first 'maxn'+1 factorials.
    """
    table = Numeric.arange(maxn+1, typecode='d')
    table[0] = 1.0
    return Numeric.add.accumulate(log(table))

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
    x = pi / 2.0 * freq / nyquist_freq
    if (x == 0.0):
        return 1.0
    else:
        return (sin(x) / x)

def max_noise_power(bins, confidence=0.99):
    """
    max_noise_power(bins, confidence=0.99):
        Return the power level that gives you some
        'confidence' that spectral noise could not cause
        that level in your power spectrum.  The total number
        of independent frequencies searched is 'bins'.
        This is P_detect in Vaughan et. al, 1994, and is also
        know as P_threshold.
    """
    return -log((1.0 - confidence) / bins)

def prob_power_series(power, signal_power, n=1, TOL=1.0e-14):
    """
    prob_power_series(power, signal_power, n=1, TOL=1.0e-14):
        Return the integrated probability from P=0 to 'power'
        that a signal with theoretical power 'signal_power'
        will show up in a power spectrum with power 'power'.
        This method evaluates the integral using an infinite
        sum and is equation 16 in Groth, 1975.
    """
    fact = exp(-(power + signal_power))
    lf = log_fact_table((power + signal_power) * 5)
    lp, lps = log(power), log(signal_power)
    sum = 0.0
    term = 1.0
    m = 0
    while (1):
        kmax = m + n
        term = fact * Numeric.add.reduce(exp((Numeric.arange(kmax) * lp + \
                                              m * lps) - (lf[0:kmax] + lf[m])))
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
        t2 = sin(2.0 * theta)
        A = t1 + ps * t2
        B = t1 + (ps - p) * t2
        sintheta = sin(theta)
        sin2theta = sintheta**2.0
        return (exp(-2.0 * ps * sin2theta) *
                (sin(A - theta) - exp(-2.0 * p * sin2theta) *
                 sin(B - theta)) / sintheta)
    (val, err) = quad(integrand, 0.0, pi/2.0, (power, signal_power, n))
    return val/pi

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
           exp(-(power + signal_power)) * \
           iv(n - 1.0, 2 * sqrt(power * signal_power))

def required_signal_power(power, confidence=0.99):
    """
    required_signal_power(power, confidence=0.99):
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
    def func(x, power=power, prob=prob):
        return prob_power_series(power, x) - prob
    def dfunc(x, power=power):
        return power_probability(power, x)
    P_signal = newton_raphson(func, dfunc, 0.0001, 100.0)
    return P_signal

def fft_sensitivity(N, bins=0, confidence=0.99):
    """
    fft_sensitivity(N, bins=0, confidence=0.99):
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
    P_threshold = max_noise_power(bins, confidence)
    return required_signal_power(P_threshold, confidence)

def rzw_sensitivity(N, zlo=-100.0, zhi=100.0, confidence=0.99):
    """
    rzw_sensitivity(N, zlo=-100.0, zhi=100.0, confidence=0.99):
        Return a measure of the weakest signal power you can
        confidently detect in an RZW (Fourier acceleration) search
        containing 'N' data points (this is the number of bins in the
        time series) and low and high acceleration values of 'zlo'
        and 'zhi'.  'confidence' is our fractional confidence in
        the result (i.e. 0.99 = 99% limit).  This calculation does not
        include the correction to sensitivity due to binning effects.
    """
    bins = N / 2.0 * (zhi - zlo) / 6.95
    P_threshold = max_noise_power(bins, confidence)
    return required_signal_power(P_threshold, confidence)

def binned_fft_sensitivity(N, dt, freq, bins=0, confidence=0.99):
    """
    binned_fft_sensitivity(N, dt, freq, bins=0, confidence=0.99):
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
    return fft_sensitivity(N, bins, confidence) / factor

def binned_rzw_sensitivity(N, dt, freq, zlo=-100.0, zhi=100.0,
                           confidence=0.99):
    """
    binned_rzw_sensitivity(N, dt, freq, zlo=-100.0, zhi=100.0,
                           confidence=0.99):
        Return a measure of the weakest signal power of frequency 'freq'
        Hz you can confidently detect in an RZW (Fourier acceleration)
        search containing 'N' data points (this is the number of bins in
        the time series) each of which was binned into 'dt' sec bins.
        Low and high acceleration values of 'zlo' and 'zhi' were used.
        'confidence' is our fractional confidence in the result (i.e.
        0.99 = 99% limit).  This calculation includes the correction to
        sensitivity due to binning effects.
    """
    bins = N / 2.0 * (zhi - zlo) / 6.95
    nyquist_freq = 0.5 / dt
    factor = binning_factor(freq, nyquist_freq)**2.0
    return fft_sensitivity(N, bins, confidence) / factor

def pulsed_fraction_limit(numphot, power_limit):
    """
    pulsed_fraction_limit(numphot, power_limit):
        Return an upper limit to the pulsed fraction of a signal
        that is in the data but was not detected.  The data
        contain a total of 'numphot' photons and the largest
        measured power is 'power_limit'.
    """
    return sqrt(4.0 * (power_limit - 1.0) / numphot)


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print '\nUsage:  powerstats.py N [confidence=0.99] [dt] [numphot]\n'
        sys.exit(0)
    elif len(sys.argv)==2:
        N = float(sys.argv[1])
        confidence = 0.99
        binned = 0
    elif len(sys.argv)==3 or len(sys.argv)==4:
        N = float(sys.argv[1])
        confidence = float(sys.argv[2])
        binned = 0
    else:
        N = float(sys.argv[1])
        confidence = float(sys.argv[2])
        binned = 1
        dt = float(sys.argv[3])
        numphot = float(sys.argv[4])
    if (binned):
        print ""
        print "       Power Stats for Binned Data"
        print "   -----------------------------------"
        print "    Number of data points = %.0f" % N
        print "      Time per sample (s) = %g" % dt
        print "  Total number of photons = %.0f" % numphot
        print "         Confidence Level = %g%%" % (100 * confidence)
        print "                 P_detect = %.2f" % \
              max_noise_power(N/2, confidence)
        print ""
        sens = (binned_fft_sensitivity(N, dt, 1.0, 0, confidence),
                binned_fft_sensitivity(N, dt, 10.0, 0, confidence),
                binned_fft_sensitivity(N, dt, 100.0, 0, confidence),
                binned_fft_sensitivity(N, dt, 1000.0, 0, confidence))
        print "                        f=1Hz     f=10Hz    f=100Hz   f=1kHz"
        print "                        -----     ------    -------   ------"
        print "      P_sensitivity  =  %-10.2f%-10.2f%-10.2f%-10.2f" % sens
        print "       Pulsed Fract  <  %-10.3g%-10.3g%-10.3g%-10.3g" % \
              (pulsed_fraction_limit(numphot, sens[0]),
               pulsed_fraction_limit(numphot, sens[1]),
               pulsed_fraction_limit(numphot, sens[2]),
               pulsed_fraction_limit(numphot, sens[3]))
        print ""
    else:
        print ""
        print "       Power Stats for Normal Data"
        print "   -----------------------------------"
        print "    Number of data points = %.0f" % N
        print "         Confidence Level = %g%%" % (100 * confidence)
        print "                 P_detect = %.2f" % \
              max_noise_power(N/2, confidence)
        sens = fft_sensitivity(N, 0, confidence)
        print "            P_sensitivity = %.2f" % sens
        print ""
