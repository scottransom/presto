from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import input
from builtins import range
import math, string, Numeric, presto, random, sys, pickle
from LeastSquares import leastSquaresFit
from orbitstuff import *

# Some admin variables
parallel = 0        # True or false
showplots = 0       # True or false
debugout = 0        # True or false
numloops = 200      # Number of orbits to test
numbetween = 16     # The number of bins to interpolate

# Observation parameters
dt = 0.000125       # The duration of each data sample
N = 2**28           # The number of points in the observation
T = N*dt            # The total observation time
ctype = 'NS'        # The type of binary companion: 'WD', 'NS', or 'BH'

# These are the minimum distances to measure from the true values
Dp = 0.00002        # fraction of orbital period
Dx = 0.0002         # fraction of projected orbital semi-major axis
De = 0.01           # eccentricity (absolute)
Dw = 0.5            # degrees (absolute)
Dt = 0.00002        # fraction of orbital period

if showplots:
    import Pgplot
if parallel:
    import mpi
    from mpihelp import *
    myid = mpi.comm_rank()
    numprocs = mpi.comm_size()
    if ctype=='WD':
        if numprocs!=3:
            raise SystemExit('You need 3 procs for the NS-WD simulation.')
    else:
        if numprocs!=5:
            raise SystemExit('You need 5 procs for a NS-NS or NS-BH simulation.')
else:
    myid = 0

# The mathematical model of the Fourier Peak.
def quadratic(parameters, x):
    a = parameters[0]
    b = parameters[1]
    c = parameters[2]
    return (a * x + b) * x + c

# Store important psrparams info in a list
def psrparams_to_list(psr):
    result = []
    result.append(psr.p)
    result.append(psr.orb.p)
    result.append(psr.orb.x)
    result.append(psr.orb.e)
    result.append(psr.orb.w)
    result.append(psr.orb.t)
    return result

# Get important psrparams info from a list
def psrparams_from_list(pplist):
    psr = presto.psrparams()
    psr.p = pplist[0]
    psr.orb.p = pplist[1]
    psr.orb.x = pplist[2]
    psr.orb.e = pplist[3]
    psr.orb.w = pplist[4]
    psr.orb.t = pplist[5]
    return psr

# Correlate a kernel with a data array.  Return the good parts.
def corr(data, kernel, numbetween, firsttime=0):
    m = len(kernel)/2
    lobin = m / numbetween
    numbins = len(data) - 2 * lobin
    if firsttime: samedata = 0
    else: samedata = 1
    result = Numeric.zeros(numbins * numbetween, 'F')
    presto.corr_complex(data, len(data), lobin, numbins, 
                        numbetween, kernel, m, result, samedata, 0)
    return result

# Open a file to save each orbit calculation
file = open('montebinopt_saves.txt','w')

# Decide what work we have to perform
work = ['p', 'x', 't', 'e', 'w']
if not parallel:
    if ctype=='WD':
        numjobs = 3
    else:
        numjobs = 5
else:
    numjobs = 1
    
# Begin the loop over the candidates
widths = []
for i in range(numloops):
    # Generate the fake pulsar parameters
    if myid==0:
        psr = fake_mspsr(companion = ctype)
        psrlist = psrparams_to_list(psr)
    else:
        psr = presto.psrparams()
        psrlist = None
    if parallel:
        psrlist = bcast_general(psrlist, myid, 0)
        if myid==0: psrlist = psrparams_to_list(psr)
        else: psr = psrparams_from_list(psrlist)
        if debugout:
            allproc_print(numprocs, 'Psr period =', psr.p)
    print('')
    print('Trial', i)
    if debugout:
        print('')
        print('   PSR mass              =', mpsr)
        print('   Companion mass        =', mc)
        print('   PSR period (s)        =', psr.p)
        print('   PSR frequency (hz)    =', 1.0/psr.p)
        print('   Orbit period (s)      =', psr.orb.p)
        print('   Orbit asini/c (lt-s)  =', psr.orb.x)
        print('   Orbit eccentricity    =', psr.orb.e)
        print('   Orbit angle (deg)     =', psr.orb.w)
        print('   Orbit time (s)        =', psr.orb.t)
        print('   Orbit Fourier Freq    =', T/psr.orb.p)
        print('   Orbit z               =', \
              presto.TWOPI*psr.orb.x/psr.p)
        print('')
            
    # Create the data set
    cand = presto.orbitparams()
    m = 0
    comb = presto.gen_bin_response(0.0, 1, psr.p, T, psr.orb , 
                                   presto.LOWACC, m)
    ind = len(comb)
    # The follwoing is performed automatically in gen_bin_resp() now
    # m = (ind / 2 + 10) * numbetween
    data = Numeric.zeros(3 * ind, 'F')
    data[ind:2*ind] = comb
    if showplots and not parallel:
        Pgplot.plotxy(presto.spectralpower(data), color='red',
                      title='Data', labx='Fourier Frequency',
                      laby='Relative Power')
        a = input("Press enter to continue...")
        Pgplot.nextplotpage(1)
        
    # Perform the loops over the Keplerian parameters
    for job in range(numjobs):
        if parallel:
            myjob = work[myid]
        else:
            myjob = work[job]
        if myjob=='p':
            Dd = Dp
            psrref = psr.orb.p
        if myjob=='x':
            Dd = Dx
            psrref = psr.orb.x
        if myjob=='t':
            Dd = Dt
            psrref = psr.orb.p
        if myjob=='e':
            Dd = De
            psrref = 1.0
        if myjob=='w':
            Dd = Dw
            psrref = 1.0
        firsttime = 1
        vals = []
        vals.append((0.0, 1.0))
        vals.append((0.0, 1.0))
        ddelta = Dd * psrref
        delta = ddelta
        while vals[-1][1] > 0.5 and vals[-2][1] > 0.5:
            for currentdelta in [delta, -delta]:
                # Adjust our candidate orbital period
                cand = copyorb(psr.orb, cand)
                if myjob=='p': cand.p = psr.orb.p + currentdelta
                if myjob=='x': cand.x = psr.orb.x + currentdelta
                if myjob=='t': cand.t = psr.orb.t + currentdelta
                if myjob=='e': cand.e = psr.orb.e + currentdelta
                if myjob=='w': cand.w = psr.orb.w + currentdelta
                # Generate the new correlation kernel
                kernel = presto.gen_bin_response(0.0, numbetween,
                                                 psr.p, T, cand, 
                                                 presto.LOWACC, m)
                # Perform the correlation
                result = corr(data, kernel, numbetween, firsttime)
                firsttime = 0
                # Convert to a power spectrum
                respow = presto.spectralpower(result)
                vals.append((currentdelta/psrref,
                             Numeric.maximum.reduce(respow)))
                if debugout:
                    # Print the most recent results
                    print('   %s:  Delta = %10.6f   Response = %8.5f' % \
                          (myjob, vals[-1][0], vals[-1][1]))
                if showplots and not parallel:
                    # Plot the results of the correlation
                    Pgplot.plotxy(respow, labx='Frequency',
                                  laby='Relative Power')
                    a = input("Press enter to continue...")
                    Pgplot.nextplotpage(1)
            # A very rough adaptive stepsize
            if abs(vals[-3][1] - vals[-1][1]) < 0.04:
                ddelta = ddelta * 2.0            
            delta = delta + ddelta
        # Fit a quadratic to the width values
        fit = leastSquaresFit(quadratic, (-1.0, 0.0, 1.0), vals)
        if debugout:
            print('\n   %sfit = %fx^2 + %fx + %f\n' % (myjob, fit[0][0],
                                                       fit[0][1],
                                                       fit[0][2]))
        width = 2.0*math.sqrt(-0.5/fit[0][0])
        if parallel:
            newwidths = mpi.gather_string(str(width), 0)
            if myid==0:
                for proc in range(numprocs):
                    psrlist.append(string.atof(newwidths[proc]))
        else:
            psrlist.append(width)
    widths.append(psrlist)
    if debugout:
        print('Widths are', widths[i])
    # Save our most recent orbit and width information
    pickle.dump(widths[i], file, 1)
file.close()

