from __future__ import print_function
import numpy as Num
import cPickle, os

n = 1000

if (0):
    # Use the following to generate the xs
    from simple_roots import newton_raphson

    rval = 0.0
    rs = Num.arange(n+1, dtype=Num.float)/n
    xs = Num.zeros(n+1, dtype=Num.float)
    def func(x):
        return Num.sin(2.0*Num.pi*x)/(2.0*Num.pi) + x - rval
    def dfunc(x):
        return Num.cos(2.0*Num.pi*x) + 1
    for (ii, rval) in enumerate(rs[:-1]):
        if (ii==n/2):
            xs[ii] = 0.5
        else:
            xs[ii] = newton_raphson(func, dfunc, 0.0, 1.0)
    xs[0] = 0.0
    xs[n] = 1.0
    cPickle.dump(xs, file("cosine_rand.pickle", "w"), 1)
else:
    pfile = os.path.join(os.environ['PRESTO'],
                         "lib", "python", "cosine_rand.pickle")
    xs = cPickle.load(file(pfile))

def cosine_rand(num):
    """cosine_rand(num):  Return num phases that are randomly distributed
          as per a sinusoid with maximum at phase=0 (0 < phase < 1).
    """
    rands = n*Num.random.random(num)
    indices = rands.astype(Num.int)
    fracts = rands-indices
    lo = Num.take(xs, indices)
    hi = Num.take(xs, indices+1)
    return fracts*(hi-lo)+lo

if __name__ == '__main__':
    from psr_utils import hist
    from Pgplot import plotxy, closeplot
    import time
    
    if (0):
        numtrials = 20
        numrandnums = 1000000
        for funct in [cosine_rand1, cosine_rand2]:
            times = []
            for jj in range(numtrials):
                tt = time.clock()
                funct(numrandnums)
                times.append(time.clock()-tt)
            print("Average time = ", Num.add.reduce(Num.asarray(times))/numtrials)
    else:
        rs = Num.arange(n+1, dtype=Num.float)/n
        plotxy(xs, rs)
        closeplot()
        hist(cosine_rand(10000), 100, color='red')
        closeplot()
