#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import sys
import numpy as num
from presto import psr_utils as pu
from presto import psr_constants as pc
from presto import parfile
from presto import bestprof
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

period = num.asarray([])
time = num.asarray([])


def parse_eph(filenm):
    global period, time
    suffix = filenm.split(".")[-1]
    if suffix=="bestprof":
        x = bestprof.bestprof(filenm)
        fs = pu.p_to_f(x.p0_bary, x.p1_bary, x.p2_bary)
        epoch = x.epochi_bary + x.epochf_bary
        T = x.T
    elif suffix=="par":
        x = parfile.psr_par(filenm)
        # Try to see how many freq derivs we have
        fs = [x.F0]
        for ii in range(1, 20):  # hopefully 20 is an upper limit!
            attrib = "F%d"%ii
            if hasattr(x, attrib):
                fs.append(getattr(x, attrib))
            else:
                break
        epoch = x.PEPOCH
        T = (x.FINISH - x.START) * 86400.0
    else:
        print("I don't recognize the file type for", filenm)
        sys.exit()
    newts = epoch + num.arange(int(T/10.0+0.5), dtype=num.float)/8640.0
    time = num.concatenate((time, newts))
    newps = 1.0 / pu.calc_freq(newts, epoch, *fs)
    period = num.concatenate((period, newps))
    print("%13.7f (%0.1f sec): " % (epoch, T), fs)


def orbeqn(Ppxt, times):
    # P = Ppsr, p = Porb, x = a*sin(i)/s, t = T_o
    phi = pc.TWOPI*(times - Ppxt[3])*86400.0/Ppxt[1]
    return Ppxt[0]*(1.0+pc.TWOPI*Ppxt[2]/Ppxt[1]*num.cos(phi))


def funct(Ppxt, times, measured):
    return orbeqn(Ppxt, times) - measured


if __name__ == '__main__':
    if len(sys.argv)==1:
        print("\nusage: fit_circular_orbit.py P_psr P_orb X_orb parfiles or bestprofs")
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
        print("Modifying TO because of negative asini/c...")
        ret[0][3] += 0.5 * (ret[0][1]/86400.0)
        ret[0][2] = abs(ret[0][2])

    print("P_orb = %.3f hrs" % (ret[0][1]/3600.0))
    print("P0  %17.15g 1" % ret[0][0])
    print("PB  %17.15g 1" % (ret[0][1]/86400.0))
    print("A1  %17.15g 1" % ret[0][2])
    print("T0  %17.15g 1" % ret[0][3])
    print("E                 0.0")
    print("OM                0.0")

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

