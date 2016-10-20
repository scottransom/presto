from __future__ import print_function
from builtins import range
from time import clock
from math import *
from Numeric import *
from presto import *
from miscutils import *
from Statistics import *

# Some admin variables
parallel = 0          # True or false
showplots = 0         # True or false
debugout = 0          # True or false
outfiledir = '/home/ransom'
outfilenm = 'monte'
pmass = 1.35                                 # Pulsar mass in solar masses
cmass = {'WD': 0.3, 'NS': 1.35, 'BH': 10.0}  # Companion masses to use
ecc = {'WD': 0.0, 'NS': 0.6, 'BH': 0.6}      # Eccentricities to use
orbsperpt = {'WD': 20, 'NS': 100, 'BH': 100} # # of orbits to avg per pt
ppsr = [0.002, 0.02, 0.2, 2.0]               # Pulsar periods to test

# Simulation parameters
numTbyPb = 100        # The number of points along the x axis
minTbyPb = 0.01       # Minimum Obs Time / Orbital Period
maxTbyPb = 10.0       # Maximum Obs Time / Orbital Period
ctype = 'BH'          # The type of binary companion: 'WD', 'NS', or 'BH'
Pb = 7200.0           # Orbital period in seconds
dt = 0.0001           # The duration of each data sample (s)
searchtype = 'ffdot'  # One of 'ffdot', 'sideband', 'shortffts'
maxTbyPb_ffdot = 11.0

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


####################################################################

# Calculate the values of our X and Y axis points
logTbyPb = span(log(minTbyPb), log(maxTbyPb), numTbyPb)
TbyPb = exp(logTbyPb)

# Open a file to save each orbit calculation
file = open(outfilenm,'w')

# The Simulation loops

# Loop over T / Porb
for x in range(numTbyPb)[51:66]:
    T = Pb * TbyPb[x]
    xb = asini_c(Pb, mass_funct2(pmass, cmass[ctype], pi / 3.0))
    eb = ecc[ctype]
    # Loop over ppsr
    for y in range(len(ppsr)):
        # Each processor calculates its own point
        if not (y % numprocs == myid):  continue
        else:
            pows = zeros(orbsperpt[ctype], 'd')
            stim = clock()
            if (TbyPb[x] <= maxTbyPb_ffdot):
                # Loop over the number of tries per point
                for ct in range(orbsperpt[ctype]):
                    if (eb == 0.0):
                        wb, tp = 0.0, ct * Pb / orbsperpt[ctype]
                    else:
                        (orbf, orbi)  = modf(ct / sqrt(orbsperpt[ctype]))
                        orbi = orbi / sqrt(orbsperpt[ctype])
                        wb, tp = orbf * 180.0, Pb * orbi
                    if debugout:
                        print('T = '+repr(T)+'  ppsr = '+repr(ppsr[y])+\
                              ' Pb = '+repr(Pb)+' xb = '+repr(xb)+' eb = '+\
                              repr(eb)+' wb = '+repr(wb)+' tp = '+repr(tp))
                    psr = psrparams_from_list([ppsr[y], Pb, xb, eb, wb, tp])
                    psr_numbins = 2 * bin_resp_halfwidth(psr.p, T, psr.orb)
                    psr_resp = gen_bin_response(0.0, 1, psr.p, T, psr.orb,
                                                psr_numbins)
                    if showplots:
                        print("The raw response:")
                        Pgplot.plotxy(spectralpower(psr_resp))
                        Pgplot.closeplot()
                    # The following places the nominative psr freq
                    # approx in bin len(data)/2
                    datalen = next2_to_n(psr_numbins * 2)
                    if datalen < 1024: datalen = 1024
                    data = zeros(datalen, 'F')
                    lo = (len(data) - len(psr_resp)) / 2
                    hi = lo + len(psr_resp)
                    data[lo:hi] = array(psr_resp, copy=1)
                    (tryr, tryz) = estimate_rz(psr, T, show=showplots)
                    tryr = tryr + len(data) / 2.0
                    numr = 200
                    numz = 200
                    dr = 0.5
                    dz = 1.0
                    if debugout:
                        #   print 'psr_numbins = '+`psr_numbins`+\
                        #   ' TbyPb[x] = '+`TbyPb[x]`+\
                        #   ' ppsr[y] = '+`ppsr[y]`+\
                        #   ' len(data) = '+`len(data)`
                        print(' tryr = %11.5f   tryz = %11.5f' % \
                              (tryr, tryz))
                    ffd = ffdot_plane(data, tryr, dr, numr,
                                      tryz, dz, numz)
                    maxarg = argmax(spectralpower(ffd.flat))
                    rind = maxarg % numr
                    zind = maxarg / numr
                    peakr = (rind * dr + int(tryr - (numr * dr) / 2.0))
                    peakz = (zind * dz + tryz - (numz * dz) / 2.0)
                    peakpow = ffd[zind][rind].real**2 + \
                              ffd[zind][rind].imag**2
                    if showplots:
                        xvals = arange(numr, typecode="d") * dr + \
                                int(tryr - (numr * dr) / 2.0)
                        yvals = arange(numz, typecode="d") * dz + \
                                tryz - (numz * dz) / 2.0
                        ffdpow = spectralpower(ffd.flat)
                        ffdpow.shape = (numz, numr)
                        Pgplot.plot2d(ffdpow, xvals, yvals, \
                                      labx="Fourier Freq (bins)", \
                                      laby="Fourier Freq Deriv", \
                                      image="astro")
                        Pgplot.closeplot()
                    if debugout:
                        print('peakr = %11.5f  peakz = %11.5f' % \
                              (peakr, peakz))
                    if (peakpow < 0.2):
                        pows[ct] = peakpow
                        rmax = peakr
                        zmax = peakz
                    else:
                        [pows[ct], rmax, zmax, rd] = \
                                   maximize_rz(data, peakr, peakz, \
                                               norm=1.0)
                    if debugout:
                        print('max_r = %11.5f  max_z = %11.5f' % \
                              (rmax, zmax))
                    if debugout:
                        print(repr(x)+'  '+repr(y)+'  '+repr(TbyPb[x])+'  ', end=' ')
                        print(repr(ppsr[y])+'  '+repr(pows[ct]))
            tim = clock() - stim
            if debugout:
                print('Time for this point was ',tim, ' s.')
            file.write('%5d  %9.6f  %8.6f  %11.9f  %11.9f  %11.9f\n' % \
                       (y * numTbyPb + x, TbyPb[x], ppsr[y], \
                        average(pows), max(pows), min(pows)))
            file.flush()
file.close()
