from umath import pi, log, exp, sin, cos
from Numeric import *
#from Multipack import quad
from simple_roots import newton_raphson
from cephes import iv

def log_fact_table(maxn):
    table = arange(maxn+1, typecode='d')
    table[0] = 1.0
    return add.accumulate(log(table))

def binning_factor(freq, nyquist_freq):
    x = pi / 2.0 * freq / nyquist_freq
    if (x == 0.0):
        return 1.0
    else:
        return (sin(x) / x)**2.0

def max_noise_power(bins, confidence=0.99):
    """
    max_noise_power(bins, confidence=0.99):
        Return the power level that gives you a certain
        'confidence' that spectral noise could not cause
        it in your power spectrum.  The total number
        of independent frequencies searched is 'bins'.
    """
    return -log((1.0 - confidence) / bins)

def prob_power_series(power, signal_power, n=1, TOL=1.0e-14):
    fact = exp(-(power + signal_power))
    lf = log_fact_table((power + signal_power) * 5)
    lp, lps = log(power), log(signal_power)
    sum = 0.0
    term = 1.0
    m = 0
    while (1):
        kmax = m + n
        term = fact * add.reduce(exp((arange(kmax) * lp + m * lps) -
                                     (lf[0:kmax] + lf[m])))
        sum = sum + term
        if (m > signal_power and term < TOL):  break
        m = m + 1
    return 1.0 - sum

def prob_power_integral(power, signal_power, n=1):
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
    # Note:  This is the integrand of the previous
    #        two functions (which integrate it from 0 to P
    return (power / signal_power)**(0.5 * (n - 1)) * \
           exp(-(power + signal_power)) * \
           iv(n - 1.0, 2 * sqrt(power * signal_power))

def fft_search_sensitivity(N, dt, freq, bins=0, confidence=0.99):
    """
    fft_search_sensitivity(N, dt, freq, bins=N/2, confidence=0.99):
        Return a measure of the weakest signal power of frequency 'freq'
        Hz you can confidently detect in an FFT search containing 'N'
        data points each of which was binned into 'dt' sec long bins.
        'bins' is only different from 0 if the number of independent
        frequencies searched does not equal N/2 (i.e. when an
        acceleration search is performed).  'confidence' is our
        fractional confidence in the result (i.e. 0.99 = 99% limit).
    """
    if not (bins): bins = N / 2
    T = N * dt
    nyquist = N / 2.0 * T
    prob = 1.0 - confidence
    P_threshold = max_noise_power(bins, confidence)
    print "P_threshold =", P_threshold
    def func(x, power=P_threshold, prob=prob):
        return prob_power_series(power, x) - prob
    def dfunc(x, power=P_threshold):
        return power_probability(power, x)
    P_signal = newton_raphson(func, dfunc, 0.0001, 100.0)
    print "P_signal = ", P_signal
    return P_signal

if __name__ == '__main__':
    import sys
    import time
    from Pgplot import plotxy, closeplot

    if len(sys.argv) < 2:
        print '\nusage:  python power_integrals.py signal_power\n'
        sys.exit(0)
    signal_power = float(sys.argv[1])
    powers = arange(signal_power * 4.0) * 0.5 + 0.1
    series_probs = []
#    integral_probs = []
#    integral_time = time.clock()
#    for power in powers:
#        integral_probs.append(prob_power_integral(power, signal_power))
#    integral_time = time.clock() - integral_time
#    print 'Integral method took', integral_time, 'seconds.'
    series_time = time.clock()
    for power in powers:
        series_probs.append(prob_power_series(power, signal_power))
    series_time = time.clock() - series_time
    print 'Series method took', series_time, 'seconds.'
#    plotxy(integral_probs, powers, labx='Power', laby='Probability',
#           color='red')
    plotxy(series_probs, powers, labx='Power', laby='Probability',
           color='blue')
    closeplot()
