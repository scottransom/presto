from time import clock
from math import *
from Numeric import *
from miscutils import *
from Statistics import *
from presto import *
from orbitstuff import *

# Some admin variables
parallel = 0          # True or false
showplots = 1         # True or false
debugout = 1          # True or false
#outfiledir = '/tmp/scratch'
outfiledir = '/home/ransom'
outfilenm = 'montebinresp'
orbsperpt = 100       # Number of orbits to average per point
pmass = 1.35          # Pulsar mass in solar masses
cmass = {'WD': 0.3, 'NS': 1.35, 'BH': 10.0}  # Companion masses to use
ecc = {'WD': 0.0, 'NS': 0.6, 'BH': 0.6}      # Eccentricities to use

# Simulation parameters
numTbyPb = 50         # The number of points along the x axis
minTbyPb = 0.01       # Minimum Obs Time / Orbital Period
maxTbyPb = 10.0       # Maximum Obs Time / Orbital Period
numppsr = 50          # The number of points along the y axis
minppsr = 0.0005      # Minimum pulsar period (s)
maxppsr = 5.0         # Maximum pulsar period (s)
ctype = 'WD'          # The type of binary companion: 'WD', 'NS', or 'BH'
Pb = 2 * 3600.0       # Orbital period in seronds
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
    e = dorbint(eo, int(T+1.0), dt, psr.orb)
    z = z_from_e(e, psr, T)
    r = T / p_from_e(e, psr)
    return (average(r), average(z))
    
# Figure out our environment 
if showplots:
    import Pgplot
if parallel:
    import mpi
    from mpihelp import *
    myid = mpi.comm_rank()
    numprocs = mpi.comm_size()
    outfilenm = (outfiledir+'/'+outfilenm+`myid`+
                 '_'+searchtype+'_'+ctype+'.out')
else:
    myid = 0
    numprocs = 1
    outfilenm = (outfiledir+'/'+outfilenm+
                 '_'+searchtype+'_'+ctype+'.out')

# Calculate the values of our X and Y axis points
logTbyPb = span(log(minTbyPb), log(maxTbyPb), numTbyPb)
logppsr = span(log(minppsr), log(maxppsr), numppsr)
TbyPb = exp(logTbyPb)
ppsr = exp(logppsr)

# Open a file to save each orbit calculation
file = open(outfilenm,'w')

# The Simulation loops

# Loop over T / Porb
for x in range(numTbyPb):
    T = Pb * TbyPb[x]
    xb = asini_c(Pb, mass_funct2(pmass, cmass[ctype], pi / 3.0))
    eb = ecc[ctype]
    # Loop over ppsr
    for y in range(numppsr):
        # Each processor calculates its own point
        if not (y % numprocs == myid):  continue
        else:
            pows = zeros(orbsperpt, 'd')
            if ((searchtype == 'ffdot' and
                 TbyPb[x] < maxTbyPb_ffdot) or
                (searchtype == 'sideband' and
                 TbyPb[x] > minTbyPb_sideband) or
                (searchtype == 'shortffts')):
                # Loop over the number of tries per point
                for ct in range(orbsperpt):
                    stim = clock()
                    if (eb == 0.0):
                        wb, tp = ct * 180.0 / orbsperpt, 0.0
                    else:
                        (orbf, orbi)  = modf(ct / sqrt(orbsperpt))
                        orbi = orbi / sqrt(orbsperpt)
                        wb, tp = orbf * 180.0, Pb * orbi
                    if debugout:
                        print 'T = '+`T`+'  ppsr = '+`ppsr[y]`+\
                              ' Pb = '+`Pb`+' xb = '+`xb`+' eb = '+\
                              `eb`+' wb = '+`wb`+' tp = '+`tp`
                    psr = psrparams_from_list([ppsr[y], Pb, xb,
                                               eb, wb, tp])
                    psr_numbins = 2 * bin_resp_halfwidth(psr.p, T, psr.orb)
                    psr_resp = gen_bin_response(0.0, 1, psr.p, T, psr.orb,
                                                psr_numbins)
                    if showplots:
                        Pgplot.plotxy(spectralpower(psr_resp))
                        Pgplot.nextplotpage()
                    if searchtype == 'ffdot':
                        # The following places the nominative psr freq
                        # approx in bin len(data)/2
                        data = zeros(next2_to_n(psr_numbins * 2), 'F')
                        lo = (len(data) - len(psr_resp)) / 2
                        hi = lo + len(psr_resp)
                        data[lo:hi] = array(psr_resp, copy=1)
                        eo = keplars_eqn(psr.orb.t, psr.orb.p,
                                         psr.orb.e, 1.0e-14)
                        (tryr, tryz) = estimate_rz(psr, T, eo)
                        tryr = T / psr.p - tryr + len(data) / 2.0
                        width = 201
                        if debugout:
                            print 'psr_numbins = '+`psr_numbins`+\
                                  ' TbyPb[x] = '+`TbyPb[x]`+\
                                  ' ppsr[y] = '+`ppsr[y]`+' eo = '+\
                                  `eo`+' len(data) = '+`len(data)`+\
                                  ' tryr = '+`tryr`+' tryz = '+`tryz`,
                        ffd = ffdot_plane(data, tryr, 0.5, width,
                                          tryz, 2.0, width)
                        maxarg = argmax(spectralpower(ffd.flat))
                        peakr = ((maxarg % width) * 0.5 +
                                 int(tryr - (width * 0.5) / 2.0))
                        peakz = ((maxarg / width) * 2.0 +
                                 tryz - (width * 2.0) / 2.0)
                        if showplots:
                            show_ffdot_plane(data, tryr, tryz)
                            Pgplot.nextplotpage()
                        if debugout:
                            print ' peakr = '+`peakr`+' peakz = '+`peakz`
                        [pows[ct], rmax, zmax, rd] = \
                                   maximize_rz(data, peakr, peakz, \
                                               norm=1.0)
                    elif searchtype == 'sideband':
                        psr_pows = spectralpower(psr_resp)
                        data = zeros(next2_to_n(psr_numbins * 2), 'd')
                        data[0:psr_numbins] = psr_pows
                    if debugout:
                        print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                        print `ppsr[y]`+'  '+`pows[ct]`
                    tim = clock() - stim
                    if debugout:
                        print 'Time for this point was ',tim, ' s.'
            if debugout:
                print `x`+'  '+`y`+'  '+`TbyPb[x]`+'  ',
                print `ppsr[y]`+'  '+`average(pows)`+'  ',
                print `max(pows)`+'  '+`min(pows)`
            file.write(`x`+'  '+`y`+'  '+`TbyPb[x]`+'  ')
            file.write(`ppsr[y]`+'  '+`average(pows)`+'  ')
            file.write(`max(pows)`+'  '+`min(pows)`+'\n')
            file.flush()
file.close()
