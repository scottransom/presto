from __future__ import print_function
from builtins import range
from time import clock
from math import *
from Numeric import *
from presto import *
from miscutils import *
from Statistics import *
from random import expovariate
import RNG

global theo_sum_pow, b_pows, bsum_pows, newpows, noise, fftlen

# Some admin variables
parallel = 0          # True or false
showplots = 0         # True or false
debugout = 0          # True or false
outfiledir = '/home/ransom'
outfilenm = 'monte'
pmass = 1.35                                 # Pulsar mass in solar masses
cmass = {'WD': 0.3, 'NS': 1.35, 'BH': 10.0}  # Companion masses to use
ecc = {'WD': 0.0, 'NS': 0.6, 'BH': 0.6}      # Eccentricities to use
orbsperpt = {'WD': 20, 'NS': 100, 'BH': 100}   # of orbits to avg per pt
ppsr = [0.002, 0.02, 0.2, 2.0]               # Pulsar periods to test

# Simulation parameters
ctype = 'BH'             # One of 'WD', 'NS', or 'BH'
Pb = 7200.0              # Orbital period in seconds
dt = 0.0001              # The duration of each data sample (s)
searchtype = 'sideband'  # One of 'ffdot', 'sideband', 'shortffts'

scope = 'PK'  # One of the following 'PK' = Parkes Multibeam (20cm)
              #                      'EB' = Effelsberg (20cm)
scopename = {'PK':'Parkes Multibeam', 'EB':'Effelsberg'}
gain = {'PK':1.5, 'EB':0.666}  # Antenna gain (Jy/K)
tsys = {'PK':23.5, 'EB':25.0}  # System noise temp (K)
bw = {'PK':64.0, 'EB':90.0}    # Bandwidth in MHz
numpol = 2                     # Number of polarizations
dt = 0.0001                    # Length of each time series bin
P_orb = 7200.0                 # The orbital period (s) used
detect_sigma = 6.0             # What sigma constitutes a detection

# Calculated parameters
sigma_t = 1000.0 * gain[scope] * tsys[scope] / \
          sqrt(numpol * dt * bw[scope] * 1e6)
dist = RNG.ExponentialDistribution(1.0)
rng = RNG.CreateGenerator(0, dist)

##################################################
# You shouldn't need to edit anyting below here. #
##################################################

# Figure out our environment 
if showplots:
    import Pgplot
if parallel:
    import mpi
    from mpihelp import *
    myid = mpi.comm_rank()
    numprocs = mpi.comm_size()
    outfilenm = (outfiledir+'/'+outfilenm+repr(myid)+
                 '_'+searchtype+'_'+ctype+'.out')
else:
    myid = 0
    numprocs = 1
    outfilenm = (outfiledir+'/'+outfilenm+
                 '_'+searchtype+'_'+ctype+'.out')


def secant(func, oldx, x, TOL=1e-6):    # f(x)=func(x)
    """
    Summary 
       Solve for a zero of function using Secant method 

    Usage 
       real = func(real) 
       real = secant(func, real, real [, TOL=real]) 
       
    Similar to Newton's method, but the derivative is estimated by divided
    difference using only function calls.  A root is estimated by
    x = x - f(x) (x - oldx)/(f(x) - f(oldx))
    where oldx = x[i-1] and x = x[i].
    """
    oldf, f = func(oldx), func(x)
    if (abs(f) > abs(oldf)):            # swap so that f(x) is closer to 0
        oldx, x = x, oldx
        oldf, f = f, oldf
    count = 0
    while 1:
        dx = f * (x - oldx) / float(f - oldf)
        if abs(dx) < TOL * (1 + abs(x)): return x - dx
        if count > 50:
            x = average([x, oldx, x - dx])
            f = func(x)
            # print "secant(%d): x=%s, f(x)=%s" % (count, x, f)
            return x
        oldx, x = x, x - dx
        oldf, f = f, func(x)
        count = count + 1
        # print "secant(%d): x=%s, f(x)=%s" % (count, x, f)

def mini_fft_sum_pows(tryamp):
    global theo_sum_pow, b_pows, bsum_pows, newpows, noise, fftlen

    fdata = rfft(newpows * tryamp + noise)
    norm = fdata[0].real
    rpred = predict_mini_r(fftlen, psr.orb.p, T)
    [b_pows[ct], rmax, rd] = \
                 maximize_r(fdata, rpred, norm=norm)
    # print 'avg(dat) = ',average(newpows * tryamp[ct] + noise)
    # print 'avg(fft) = ',average(spectralpower(fdata)[1:]/norm)
    # print tryamp
    if debugout:
        print('Nyquist = '+repr(fftlen/2))
        print('   rpred = %10.3f  power = %10.7f' % \
              (rpred, b_pows[ct]))
    bsum_pows[ct] = b_pows[ct]
    if (TbyPb[x] > 2.0):
        for harmonic in arange(int(TbyPb[x]-1.0))+2:
            hrpred = predict_mini_r(fftlen, harmonic * \
                                    psr.orb.p, T)
            [tmppow, hrmax, rd] = \
                     maximize_r(fdata, hrpred, norm=norm)
            bsum_pows[ct] = bsum_pows[ct] + tmppow
            if debugout:
                print('  hrpred = %10.3f  power = %10.7f' % \
                      (hrpred, tmppow))
    if debugout:
        print('  r = %10.3f  meas_r = %10.3f '\
              'alias_r = %10.3f' % \
              (fftlen * psr.orb.p / T, rmax,
               alias(rmax, fftlen/2)))
        print('  p = %10.3f  meas_p = %10.3f '\
              'alias_p = %10.3f' % \
              (psr.orb.p, rmax * T / fftlen,
               alias(rmax, fftlen/2) * T / fftlen))
        print('  BigPow = %10.7f  SumPow = %10.7f' % \
              (b_pows[ct], bsum_pows[ct]))
    return bsum_pows[ct] - theo_sum_pow
    
def psrparams_from_list(pplist):
    psr = psrparams()
    psr.p = pplist[0]
    psr.orb.p = pplist[1]
    psr.orb.x = pplist[2]
    psr.orb.e = pplist[3]
    psr.orb.w = pplist[4]
    psr.orb.t = pplist[5]
    return psr

def predict_mini_r(fftlen, Pb, T):
    nyquist = fftlen / 2
    r = fftlen * Pb / T
    if (r > nyquist): rpred = alias(r, nyquist)
    else: rpred = r
    return rpred

def slice_resp(psr, T, response):
    c1 = TWOPI * psr.orb.x / (psr.orb.p * sqrt(1.0 - psr.orb.e * psr.orb.e))
    c2 = psr.orb.e * cos(psr.orb.w * DEGTORAD)
    v1 = c1 * (c2 + 1.0);
    v2 = c1 * (c2 - 1.0);
    if (v1 < v2):  lo, hi = v2, v1
    else: lo, hi = v1, v2
    lo = len(response)/2 - int(fabs(T * fabs(lo) / \
                                    (psr.p * (1.0 + fabs(lo)))))
    hi = len(response)/2 + int(fabs(T * fabs(hi) / \
                                    (psr.p * (1.0 + fabs(hi)))))
    diff = hi-lo
    newlen = (diff/100 + 1)*100
    newresp = zeros(newlen, 'f')
    newresp[0:diff] = response[lo:hi]
    return newresp

####################################################################

# Calculate the values of our X and Y axis points
TbyPb = arange(1.15, 10.15, 0.2)

# Open a file to save each orbit calculation
file = open(outfilenm,'w')

# The Simulation loops

# Loop over T / Porb
xb = asini_c(Pb, mass_funct2(pmass, cmass[ctype], pi / 3.0))
for x in range(len(TbyPb)):
    T = Pb * TbyPb[x]
    N = T / dt
    # Loop over ppsr
    for y in range(len(ppsr)):
        # Each processor calculates its own point
        z = 2 * pi * xb / ppsr[y]
        if not (y % numprocs == myid):  continue
        else:
            b_pows = zeros(orbsperpt[ctype], 'd')
            tryamp = zeros(orbsperpt[ctype], 'd')
            bsum_pows = zeros(orbsperpt[ctype], 'd')
            fftlen = 0
            # Loop over the number of tries per point
            for ct in range(orbsperpt[ctype]):
                stim = clock()
                if (ecc[ctype] == 0.0):
                    wb, tp = 0.0, ct * Pb / orbsperpt[ctype]
                else:
                    (orbf, orbi)  = modf(ct / sqrt(orbsperpt[ctype]))
                    orbi = orbi / sqrt(orbsperpt[ctype])
                    wb, tp = orbf * 180.0, Pb * orbi

                # Generate the PSR response
                psr = psrparams_from_list([ppsr[y], Pb, xb, ecc[ctype], wb, tp])
                psr_numbins = 2 * bin_resp_halfwidth(psr.p, T, psr.orb)
                psr_resp = gen_bin_response(0.0, 1, psr.p, T, psr.orb,
                                            psr_numbins)
                if debugout:
                    print('T = %9.3f  Pb = %9.3f  Ppsr = %9.7f' % \
                          (T, psr.orb.p, psr.p))

                newpows = slice_resp(psr, T, spectralpower(psr_resp))
                if showplots:
                    print("The raw response:")
                    Pgplot.plotxy(newpows)
                    Pgplot.closeplot()
                fftlen = len(newpows)
                noise = rng.sample(fftlen)
                tryamp[ct] = 500.0
                theo_sum_pow = powersum_at_sigma(detect_sigma,
                                                 int(T/psr.orb.p))
                if debugout:
                    print('theo_sum_pow = ', theo_sum_pow)
                newloop = 1
                tryamp[ct] = secant(mini_fft_sum_pows, tryamp[ct]/2,
                                    tryamp[ct], 0.01)
                # Pgplot.plotxy(spectralpower(fdata)[1:]/norm, \
                #              arange(len(fdata))*T/fftlen, \
                #              labx='Orbital Period (s))', \
                #              laby='Power')
                # Pgplot.closeplot()
                #print '  BigPow = %10.7f  SumPow = %10.7f  S(mJy) = %10.5f' % \
                #      (b_pows[ct], bsum_pows[ct]-theo_sum_pow, 2 * sigma_t * sqrt(tryamp[ct]/N))
                tim = clock() - stim
                if debugout:
                    print('Time for this point was ',tim, ' s.')
        # Note:  The output contains the average value of tryamp.  To convert this
        #        to a minimum flux density, use the formula
        #               S(mJy) = 2 * sigma_t * sqrt(tryamp / N)   
        file.write('%9.6f  %8.6f  %10d  %7d  %13.9f  %13.9f  %13.7f\n' % \
                   (TbyPb[x], ppsr[y], N, fftlen, average(b_pows),
                    average(bsum_pows), average(tryamp)))
        file.flush()
file.close()
