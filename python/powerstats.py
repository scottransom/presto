#!/usr/bin/env python
import Numeric
from umath import pi, log, exp, sin, sqrt
from simple_roots import newton_raphson
from cephes import iv, chdtri, ndtr, ndtri
#from Multipack import quad

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
        return -log((1.0 - confidence) / bins)
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
    return sqrt(4.0 * (power_limit - 1.0) / numphot)

def answer_yes(question):
    yes = ['', 'Y', 'y', 'Yes', 'yes', 'YES',
           'T', 't', 'True', 'true', 'TRUE']
    return raw_input('\n'+question) in yes

def ask_float(question, default=None):
    while 1:
        ans = raw_input('\n'+question)
        if not ans:
            ans = default
        try:
            return float(ans)
        except (ValueError, TypeError):
            print "\nThat was not a valid number.  Try again...\n"

def ask_int(question, default=None):
    while 1:
        ans = raw_input('\n'+question)
        if not ans:
            ans = default
        try:
            return int(ans)
        except (ValueError, TypeError):
            print "\nThat was not a valid number.  Try again...\n"

if __name__ == '__main__':
    print "\nPower Statistics Calculation Routine"
    print "       Scott Ransom, June 2000\n"

    conf = ask_float(\
        "What confidence level would you like to use?  [0.99]  ", 0.99)
    Ntot = ask_int(\
        "How many data points were FFTd (N)?  ")
    dt = ask_float("What was the length in time (s) of each bin?  ")
    T = Ntot * dt
    P_max = ask_float(\
        "What was the maximum normalized power found?  ")
    rlo = 1
    rhi = Ntot / 2
    if answer_yes(\
        "Was this an RZW acceleration search (y/n)?  [y]  "):
        rlo = T * ask_float(\
                "What was the lowest freq searched (Hz)?   [1.0]  ", 1.0)
        rhi = T * ask_float(\
                "What was the highest freq searched (Hz)? [%.2f]  " %
                ((Ntot/2.0)/T), (Ntot/2.0)/T)
        zlo = ask_float(\
                "What was the lowest 'z' value searched?  [-100]  ", -100.0)
        zhi = ask_float(\
                "What was the highest 'z' value searched?  [100]  ", 100.0)
        Nsearch = (rhi - rlo) * (zhi - zlo + 1.0) / 6.95
    else:
        Nsearch = ask_int(\
                "How many independent bins were searched?  [N/2]  ", Ntot/2)
    if answer_yes(\
        "Was the data composed of binned counts (y/n)?  [y]  "):
        numphot = ask_int("How many counts (photons) were there?  ")
        lofreq, hifreq = rlo / T, rhi / T
        trial_freqs = (10.0**(Numeric.arange(7.0)-2.0)).tolist()
        trial_freqs = filter(lambda x:  x > lofreq and x < hifreq,
                             trial_freqs)
        print "\nThe trial frequencies (Hz) are:", trial_freqs
        if answer_yes(\
            "Would you like to add any more?  [y]  "):
            new_freq = ask_float(\
                "Enter a frequency (Hz) or '0' to stop.  ")
            while (new_freq):
                trial_freqs.append(new_freq)                
                new_freq = ask_float(\
                    "Enter a frequency (Hz) or '0' to stop.  ")
        trial_freqs.sort()
        print "\n\nCalculating...\n\n"
        print ""
        print "         Power Stats for Binned Data"
        print "     -----------------------------------"
        print "      Number of data points = %.0f" % Ntot
        print "        Time per sample (s) = %g" % dt
        print "    Total number of photons = %.0f" % numphot
        print "           Confidence Level = %g%%" % (100 * conf)
        print " Number of independent bins = %.2e" % Nsearch
        print " Threshold Power (P_detect) > %.2f" % \
              max_noise_power(Nsearch, 1, conf)
        ulim = required_signal_power(P_max, 1, conf)
        print "    Max Power Found (P_max) = %.2f" % P_max
        print " Max Signal Power (P_limit) < %.2f" % ulim
        print "  Pulsed Fraction (P_limit) < %.3g" % \
              pulsed_fraction_limit(numphot, ulim)
        print ""
        sens = []
        ulim = []
        for f in trial_freqs:
            sens.append(binned_fft_sensitivity(Ntot, dt, f, Nsearch, 1, conf))
            ulim.append(required_signal_power(P_max, 1, conf))
        print "          Freq (Hz)  = ",
        for f in trial_freqs:
            print " f=%-7g" % (f),
        print '\n                       '+'-'*len(trial_freqs)*11
        print "  Power Sensitivity  > ",
        for s in sens:
            print " %-8.2f " % (s),
        print ''
        pfract = []
        for s in sens:
            pfract.append(pulsed_fraction_limit(numphot, s))
        print "    Pulsed Fraction  < ",
        for p in pfract:
            print " %-8.3g " % (p),
        print '\n'
    else:
        print "\n\nCalculating...\n\n"
        print ""
        print "         Power Stats for Normal Data"
        print "     -----------------------------------"
        print "      Number of data points = %.0f" % Ntot
        print "           Confidence Level = %g%%" % (100 * conf)
        print " Number of independent bins = %.2e" % Nsearch
        print " Threshold Power (P_detect) > %.2f" % \
              max_noise_power(Nsearch/2, 1, conf)
        sens = fft_sensitivity(Ntot, Nsearch, 1, conf)
        print "          Power Sensitivity > %.2f" % sens
        ulim = required_signal_power(P_max, 1, conf)
        print "    Max Power Found (P_max) = %.2f" % P_max
        print " Max Signal Power (P_limit) < %.2f" % ulim
        print ""
