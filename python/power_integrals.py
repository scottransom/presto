from umath import log, exp, sin, cos
from Numeric import *
from Multipack import quad

def log_fact_table(maxn):
    table = arange(maxn+1, typecode='d')
    table[0] = 1.0
    return add.accumulate(log(table))

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

if __name__ == '__main__':
    import sys
    import time
    from Pgplot import plotxy, closeplot

    if len(sys.argv) < 2:
        print '\nusage:  python power_integrals.py signal_power\n'
        sys.exit(0)
    signal_power = float(sys.argv[1])
    powers = arange(signal_power * 20.0) * 0.1 + 0.1
    series_probs = []
    integral_probs = []
    integral_time = time.clock()
    for power in powers:
        integral_probs.append(prob_power_integral(power, signal_power))
    integral_time = time.clock() - integral_time
    print 'Integral method took', integral_time, 'seconds.'
    series_time = time.clock()
    for power in powers:
        series_probs.append(prob_power_series(power, signal_power))
    series_time = time.clock() - series_time
    print 'Series method took', series_time, 'seconds.'
    plotxy(integral_probs, powers, labx='Power', laby='Probability',
           color='red')
    plotxy(series_probs, powers, labx='Power', laby='Probability',
           color='blue')
    closeplot()
