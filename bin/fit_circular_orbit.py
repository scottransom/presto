#!/usr/bin/env python
import numpy as num
import psr_utils as pu
import parfile, bestprof, sys
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

period = []
time = []

def parse_eph(filenm):
    global period, time
    suffix = filenm.split(".")[-1]
    if suffix=="bestprof":
        x = bestprof.bestprof(filenm)
        f0, f1, f2 = pu.p_to_f(x.p0_bary, x.p1_bary, x.p2_bary)
        f3 = 0.0
        epoch = x.epochi_bary + x.epochf_bary
        T = x.T
    elif suffix=="par":
        x = parfile.psr_par(filenm)
        if not hasattr(x, "F1"): x.F1 = 0.0
        if not hasattr(x, "F2"): x.F2 = 0.0
        if not hasattr(x, "F3"): x.F3 = 0.0
        f0, f1, f2, f3 = x.F0, x.F1, x.F2, x.F3
        epoch = x.PEPOCH
        T = (x.FINISH - x.START) * 86400.0
    else:
        print "I don't recognize the file type for", filenm
        sys.exit()
    for minute in num.arange(int(T/10.0+0.5), dtype=num.float):
        t = epoch + minute/8640.0
        time.append(t)
        period.append(1.0 / pu.calc_freq(t, epoch, f0, f1, f2, f3))
    print "%13.7f:  %12.7f Hz  %10.3e Hz/s  %10.3e Hz/s/s  (%0.1f sec)" % \
          (epoch, f0, f1, f2, T)

def orbeqn(Ppxt, times):
    # P = Ppsr, p = Porb, x = a*sin(i)/s, t = T_o
    phi = pu.TWOPI*(times - Ppxt[3])*86400.0/Ppxt[1]
    return Ppxt[0]*(1.0+pu.TWOPI*Ppxt[2]/Ppxt[1]*num.cos(phi))

def funct(Ppxt, times, measured):
    return orbeqn(Ppxt, times) - measured

if __name__ == '__main__':
    if len(sys.argv)==1:
        print "\nusage: fit_circular_orbit.py P_psr P_orb X_orb parfiles or bestprofs"
        exit(0)
    Ppsr = float(sys.argv[1])
    Porb = float(sys.argv[2])*86400.0
    Xorb = float(sys.argv[3])
    for infile in sys.argv[4:]:
        parse_eph(infile)
    Torb = min(time)
    period = num.asarray(period, dtype=float)
    time = num.asarray(time, dtype=float)
    ret = leastsq(funct, [Ppsr, Porb, Xorb, Torb], args=(time, period))
    To = ret[0][3]

    if (ret[0][2] < 0.0):
        print "Modifying TO because of negative asini/c..."
        ret[0][3] += 0.5 * (ret[0][1]/86400.0)
        ret[0][2] = abs(ret[0][2])

    print "P_orb = %.3f hrs" % (ret[0][1]/3600.0)
    print "P0  %17.15g 1" % ret[0][0]
    print "PB  %17.15g 1" % (ret[0][1]/86400.0)
    print "A1  %17.15g 1" % ret[0][2]
    print "T0  %17.15g 1" % ret[0][3]
    print "E                 0.0"
    print "OM                0.0"

    T = max(time)-min(time)
    model_time = num.arange(min(time)-0.1*T, max(time)+0.1*T, 0.01)

    plt.figure()
    plt.plot(time-model_time[0],
             (period-ret[0][0])*1000.0, '.')
    plt.plot(model_time-model_time[0],
             (orbeqn(ret[0], model_time)-ret[0][0])*1000.0, 'r')
    plt.xlabel("Days + %.7f"%model_time[0])
    plt.ylabel("Pulsar Period - %.7f (ms)"%(ret[0][0]*1000.0))
    plt.show()

