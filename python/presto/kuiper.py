from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import numpy as num
from presto import Pgplot
from functools import reduce

def noverk(n,k):
    # This is the combinations formula
    return float(reduce(lambda a,b: a*(n-b)/(b+1), range(k),1))

def Tt(t, z, N):
    overN = 1.0/float(N)
    y = z + t*overN
    return y**(t-3.0) * (y**3.0*N
                         - y*y*t*(3.0-2.0*overN)*overN
                         - (t*(t-1.0)*(t-2.0))*overN*overN)

def kuiper_prob(D, N):
    # From section 14.3 in Numerical Recipes
    EPS1 = 1e-6
    EPS2 = 1e-12
    en = num.sqrt(N)
    lamda = (en + 0.155 + 0.24 / en) * D
    if (lamda < 0.4): return 1.0
    probks = termbf = 0.0
    a2 = -2.0 * lamda * lamda
    for ii in range(1, 100):
        a2ii2 = a2 * ii * ii
        term = 2.0 * (-2.0 * a2ii2 - 1.0) * num.exp(a2ii2)
        probks += term
        if (num.fabs(term) <= EPS1*termbf or
            num.fabs(term) <= EPS2*probks):
            return probks
        termbf = num.fabs(term)
    return 1.0

def kuiper_prob2(D, N):
    # From Paltani 2004, eqn 3 (for large N)
    EPS1 = 1e-6
    EPS2 = 1e-12
    z = D * num.sqrt(N)
    term1bf = term2bf = 0.0
    term1 = term2 = 0.0
    for m in range(1, 1000):
        x = 4.0*m*m*z*z
        term = 2.0 * (x - 1.0) * num.exp(-0.5*x)
        term1 += term
        if (num.fabs(term) <= EPS1*term1bf or
            num.fabs(term) <= EPS2*term1):
            break
        term1bf = num.fabs(term1)
    for m in range(1, 1000):
        x = 4.0*m*m*z*z
        term = m * m * (x - 3.0) * num.exp(-0.5*x)
        term2 += term
        if (num.fabs(term) <= EPS1*term2bf or
            num.fabs(term) <= EPS2*term2):
            break
        term2bf = num.fabs(term2)
    return term1 - 8.0*z/(3.0*num.sqrt(N)) * term2

def kuiper_prob3(D, N):
    # From Paltani 2004, eqn 6 (for large D)
    # note:  this equation does not seem consistent with the other 2...
    EPS1 = 1e-6
    EPS2 = 1e-12
    prob = termbf = 0.0
    for t in range(1000):
        term = noverk(N, t) * (1.0-D-t/float(N))**(N-t-1) * Tt(t, D, N)
        prob += term
        if (num.fabs(term) <= EPS1*termbf or
            num.fabs(term) <= EPS2*prob):
            return prob
        termbf = num.fabs(term)
    return 1.0

def kuiper_uniform_test(data, output=0):
    """
    kuiper_uniform_test(data, output=0):
       Conduct a Kuiper test on the data.  The data must be values
       within [0,1) (e.g. phases from a periodicity search).  They
       will be compared to a uniform distribution.  The return value
       is the probability that the data is uniformly distributed.
    """
    sdata = num.asarray(data)
    N = sdata.size
    sdata.sort()
    f0 = num.arange(N, dtype=num.float64)/N
    fn = (num.arange(N, dtype=num.float64)+1.0)/N
    Dp = (fn - sdata).max()
    Dm = (sdata - f0).max()
    D = Dp + Dm
    P = kuiper_prob(D, N)
    if (output):
        xs = (num.arange(N+3, dtype=num.float64)/(N+2.0)).repeat(2)[1:-1]
        ys = num.concatenate((num.asarray([0.0]), sdata, num.asarray([1.0]))).repeat(2)
        Pgplot.plotxy(ys, xs, rangex=[-0.03, 1.03], rangey=[-0.03, 1.03], aspect=1.0, 
                      labx="Fraction of Data", laby="Cumulative Value", width=2)
        Pgplot.plotxy(num.asarray([0.0, 1.0]), num.asarray([0.0, 1.0]), width=1)
        Pgplot.closeplot()
        print("Max distance between the cumulative distributions (D) = %.5g" % D)
        print("Prob the data is from the specified distrbution   (P) = %.3g" % P)
    return (D, P)

if __name__=="__main__":
    if (0):
        from kstest import *
        for ii in range(4):
            N = 200
            data = num.random.random(N)
            print("-------")
            print("K-S:")
            (D, P) = KS_test(data, cum_uniform_dist, output=1)
            print("Kuiper:")
            (D, P) = kuiper_uniform_test(data, output=1)
    if (1):
        p1s = []
        p2s = []
        p3s = []
        N = 50
        Ds = num.arange(0.01, 0.6, 0.01)
        for D in Ds:
           p1s.append(kuiper_prob(D, N))
           p2s.append(kuiper_prob2(D, N))
           p3s.append(kuiper_prob3(D, N))
        Pgplot.plotxy(num.log10(num.asarray(p1s)), Ds, color='red')
        Pgplot.plotxy(num.log10(num.asarray(p2s)), Ds, color='blue')
        Pgplot.plotxy(num.log10(num.asarray(p3s)), Ds, color='green')
        Pgplot.closeplot()
