import string, random, sys, cPickle
from math import *
from Numeric import *
from miscutils import *
from presto import *
from orbitstuff import *

# Some admin variables
parallel = 0          # True or false
showplots = 1         # True or false
debugout = 1          # True or false
outfilenm = '/home/ransom/montebinresp_saves'
orbsperpt = 100       # Number of orbits to test
pmass = 1.35          # Pulsar mass in solar masses
cmass = {'WD': 0.3, 'NS': 1.35, 'BH': 10.0}  # Companion masses to use
ecc = {'WD': 0.0, 'NS': 0.6, 'BH': 0.6}      # Eccentricities to use

# Simulation parameters
numTbyPb = 100        # The number of points along the x axis
minTbyPb = 0.01       # Minimum Obs Time / Orbital Period
maxTbyPb = 10.0       # Maximum Obs Time / Orbital Period
numppsr = 100         # The number of points along the y axis
minppsr = 0.0005      # Minimum pulsar period (s)
maxppsr = 5.0         # Maximum pulsar period (s)
ctype = 'WD'          # The type of binary companion: 'WD', 'NS', or 'BH'
dt = 0.0001           # The duration of each data sample (s)
searchtype = 'ffdot'  # One of 'ffdot', 'sideband', 'shortffts'
maxTbyPb_ffdot = 1.0
minTbyPb_sideband = 0.3
fftlen_shortffts = 0.05

##################################################
# You shouldn't need to edit anyting below here. #
##################################################

def psrparams_from_list(pplist):
    psr = psrparams()
    psr.p = pplist[0]
    psr.orb.p = pplist[1]
    psr.orb.x = pplist[2]
    psr.orb.e = pplist[3]
    psr.orb.w = pplist[4]
    psr.orb.t = pplist[5]
    return psr

def estimate_rz(psr, T, eo=0.0):
    """
    estimate_rz(psr, T, eo=0.0):
        Return an estimate of the average Fourier frequency ('r')
        and Fourier F-dot ('z') values in bins given the psrparams
        'psr', the length of the observation 'T' in sec,  and the
        inital Eccentric anomaly 'eo'.
    """
    import LeastSquares
    dt = 1.0
    e = make_orbit(psr, int(T+1.0), dt, eo)
    z = z_from_e(e, psr, T)
    r = T / p_from_e(e, psr)
    return (average(r), average(z))
    
outfilenm = outfilenm+'_'+searchtype+'.txt'

# Calculate the values of our X and Y axis points
logTbyPb = span(log(minTbyPb), log(maxTbyPb), numTbyPb)
logppsr = span(log(minppsr), log(maxppsr), numppsr)
TbyPb = exp(logTbyPb)
ppsr = exp(logppsr)
# Adjust the number of points in the time series so as to
# keep the ppsr in the center of the FFT freq range
# (i.e. fpsr = 1 / ppsr = Nyquist freq / 2)
N = 4 / ppsr
T = N * dt

# Figure out our environment 
if showplots:
    import Pgplot
if parallel:
    import mpi
    from mpihelp import *
    myid = mpi.comm_rank()
    numprocs = mpi.comm_size()
else:
    myid = 0
    numprocs = 1

# Open a file to save each orbit calculation
file = open(outfilenm,'w')

# The Simulation loops

# Loop over T / Porb
for x in range(numTbyPb):
    Pb = T[x] / TbyPb[x]
    xb = asini_c(Pb, mass_funct2(pmass, cmass[ctype], pi / 3.0))
    eb = ecc[ctype]
    # Loop over ppsr
    for y in range(numppsr):
        # Each processor calculates its own point
        if not (y % numprocs == myid):  continue
        else:
            avgpow = 0.0
            if searchtype == 'ffdot':
                if Pb > maxTbyPb_ffdot:
                    print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                    print `ppsr[y]`+'  '+`avgpow`
                    file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  '+
                               `ppsr[y]`+'  '+`avgpow`)
                    continue
            elif searchtype == 'sideband':
                if Pb < minTbyPb_sideband:
                    print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                    print `ppsr[y]`+'  '+`avgpow`
                    file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  '+
                               `ppsr[y]`+'  '+`avgpow`)
                    continue
            # Loop over the number of tries per point
            for ct in range(orbsperpt):
                if (eb == 0.0):
                    wb, tp = 0.0, 0.0
                else:
                    (orbf, orbi)  = modf(ct / sqrt(orbsperpt))
                    orbi = orbi / sqrt(orbsperpt)
                    wb, tp = orbf * 180.0, Pb * orbi
                psr = psrparams_from_list([ppsr[y], Pb, xb, eb, wb, tp])
                print numTbyPb, psr.p, psr.orb.p, psr.orb.x
                psr_halfwidth = bin_resp_halfwidth(psr.p, psr.orb)
                psr_numbins = 2 * numbetween * psr_halfwidth
                psr_resp = gen_bin_response(0.0, 1, psr.p, T[x], psr.orb, \
                                            psr_numbins)
                if searchtype == 'ffdot':
                    # The following places the nominative psr freq
                    # in bin # ~len(searchdata) / 2
                    search = zeros(next2_to_n(psr_numbins * 2), 'F')
                    lo = (len(search) - len(psr_resp)) / 2
                    hi = lo + len(psr_resp)
                    search[lo:hi] = psr_resp
                    eo = keplars_eqn(psr.orb.t, psr.orb.p, psr.orb.e, 1.0e-14)
                    (tryr, tryz) = estimate_fdot(psr, T[x], eo)
                    tryr = T[x] / psr.p - tryr + len(searchdata) / 2.0
                    nrz = 201
                    ffd = ffdot_plane(search, tryr, 0.5, nrz, tryz, 2.0, nrz)
                    maxarg = argmax(spectralpower(ffd.flat))
                    [maxpow, rmax, zmax, rd] = maximize_rz(search,
                                                           maxarg % nrz,
                                                           maxarg / nrz)
                    avgpow = avgpow + maxpow
                    if debugout:
                        print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                        print `ppsr[y]`+'  '+`maxpow`
                        file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  '+
                                   `ppsr[y]`+'  '+`maxpow`)                    
                elif searchtype == 'sideband':
                    psr_pows = spectralpower(psr_resp)
                    search = zeros(next2_to_n(psr_numbins * 2), 'd')
                    search[0:psr_numbins] = psr_pows
            avgpow = avgpow / orbsperpt
            print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
            print `ppsr[y]`+'  '+`avgpow`
            file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  '+
                       `ppsr[y]`+'  '+`avgpow`)                    
            # Save our most recent orbit and search information
#            cPickle.dump(widths[i], file, 1)
file.close()
