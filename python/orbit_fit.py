from Numeric import *
from miscutils import *
from Scientific.Functions.Interpolation import *
from Scientific.Statistics import *
from presto import rfft, frotate, drotate, tofloatvector, measure_phase, pfd
from Pgplot import *
from Multipack import leastsq
from TableIO import readColumns
from sys import argv

fwhm = 0.04 # FWHM of expcos profile to fit for first TOA analysis
fit_template = expcos_profile(42, 0.0, fwhm)

class profdat:
    pass

def read_bestprof(proffile):
    prof = profdat()
    infile = open(proffile)
    for line in infile.readlines():
        startline = line[:22]
        if (startline == "# Epoch_topo       =  "):
            if not (line[22:25] == "N/A"):
                profdat.MJD_topo = float(line[22:])
        elif (startline == "# Epoch_bary (MJD) =  "):
            if not (line[22:25] == "N/A"):
                profdat.MJD_bary = float(line[22:])
        elif (startline == "# T_sample         =  "):
            profdat.dt = float(line[22:])
        elif (startline == "# Data Folded      =  "):
            profdat.N = int(line[22:])
        elif (startline == "# Profile Bins     =  "):
            profdat.proflen = int(line[22:])
        elif (startline == "# Profile Avg      =  "):
            profdat.profavg = float(line[22:])
        elif (startline == "# Profile StdDev   =  "):
            profdat.profstd = float(line[22:])
        elif (startline == "# P_topo (ms)      =  "):
            if not (line[22:25] == "N/A"):
                profdat.p_topo = float(line[22:40])/1000.0
                profdat.p_topo_err = float(line[43:])/1000.0
        elif (startline == "# P'_topo (s/s)    =  "):
            if not (line[22:25] == "N/A"):
                profdat.pd_topo = float(line[22:40])
                profdat.pd_topo_err = float(line[43:])
        elif (startline == "# P_bary (ms)      =  "):
            if not (line[22:25] == "N/A"):
                profdat.p_bary = float(line[22:40])/1000.0
                profdat.p_bary_err = float(line[43:])/1000.0
        elif (startline == "# P'_bary (s/s)    =  "):
            if not (line[22:25] == "N/A"):
                profdat.pd_bary = float(line[22:40])
                profdat.pd_bary_err = float(line[43:])
        elif (startline == "# P_orb (s)        =  "):
            if not (line[22:25] == "N/A"):
                profdat.porb = float(line[22:])
        elif (startline == "# asin(i)/c (s)    =  "):
            if not (line[22:25] == "N/A"):
                profdat.xorb = float(line[22:])
        elif (startline == "# eccentricity     =  "):
            if not (line[22:25] == "N/A"):
                profdat.eorb = float(line[22:])
        elif (startline == "# w (rad)          =  "):
            if not (line[22:25] == "N/A"):
                profdat.worb = float(line[22:])
        elif (startline == "# T_peri           =  "):
            if not (line[22:25] == "N/A"):
                profdat.torb = float(line[22:])
        elif (startline == "######################"):
            break
    infile.close()
    profdat.profile = readColumns(proffile, "#")[1]
    return profdat

def write_itoa(toa, toaerr, name, freq, dm, obs):
    """
    ITOA Format

    columns     item
    1-2      ignored, but must not be blank (often PSRNAME goes here)
    10-28    TOA (decimal point must be in column 15)
    29-34    TOA uncertainty (microseconds)
    35-45    Observing frequency (MHz)
    46-55    DM correction (pc cm^-3)
    58-59    Observatory (two-letter code)
    """
    print "%8s %19.13f %5.1f %10.5f %9.4f  %2s" % \
          (name, toa, toaerr, freq, dm, obs)

def corr(profile, template):
    ft = rfft(template)
    fp = rfft(profile)
    return rfft(ft*fp, 1)

def maxphase(profile, template):
    return float(argmax(corr(profile, template))) / len(profile)

def interp_prof(profile, zoom=10):
    nprofile = zeros(len(profile)+1, 'd')
    nprofile[-1] = profile[0]
    nprofile[:-1] = profile
    profphases = arange(0.0, 1.0001, 1.0/len(profile))
    iprof_func = InterpolatingFunction((profphases,), nprofile)
    iprof = zeros(zoom * len(profile), 'd')
    for i in xrange(zoom * len(profile)):
        iprof[i] = iprof_func(float(i) / (zoom * len(profile)))
    return iprof

def shift_to_phase(profile, phase, fwhm, zoom=10):
    template = gaussian_profile(fabs(arange(0.5, -0.5, -1.0/
                                            (zoom * len(profile)))), fwhm)
    iprof = interp_prof(profile, zoom)
    peak = argmax(corr(iprof, template))
    fprof = tofloatvector(iprof)
    frotate(fprof, len(fprof), peak-len(fprof)/2)
    return fprof

def old_measure_phase(profile, template, fwhm, zoom=10):
    template = gaussian_profile(len(profile)*zoom, 0.5, fwhm)
    iprof = interp_prof(profile, zoom)
    return maxphase(iprof, template)

def toaorbeqn(xpto, times):
    # x = a*sin(i)/c, p = orbital period, t = Time of Peri, o = offset
    return xpto[0] * sin(TWOPI / xpto[1] * (times - xpto[2])) + xpto[3]

def toafunct(xpto, times, delays):
    return toaorbeqn(xpto, times) - delays

def periodorbeqn(args, goodts):
    [ppsr, porb, xorb, torb] = args
    voverc = TWOPI * xorb * cos(TWOPI * (torb + goodts)/porb) / porb
    return (1.0 + voverc) * ppsr

def orbresid(args, goodts, goodps):
    return periodorbeqn(args, goodts) - goodps

# These are barycentric profs and standard deviations of the profiles

print "Reading profiles",
profs = []
for (ii, file) in indexed(argv[1:]):
    print ".",
    profs.append(read_bestprof(file))
print "\n"

ps = []
times = []
for prof in profs:
    ps.append(prof.P_bary)
    ts.append(prof.MJD_bary)

ps = asarray(ps)
ts = asarray(ts)
ts = (ts - ts[0])*86400.0

ret = leastsq(orbresid, [1.0/191.246, 5700.0, 0.04, 0.0],
              args=(ts, ps))

plotxy(ps*1000.0, ts, line=None, symbol=2,
       labx="Time (s)", laby="Pulse Period (ms)",
       device="4U_1820-303_orbit.ps/CPS")
fullts = arange(1000.0)/1000.0*max(ts)
plotxy(periodorbeqn(ret[0], fullts)*1000.0, fullts,
       line=1, symbol=None, color='red')
closeplot()

print "Best Ppsr (ms)  = %17.15g" % ret[0][0]
print "Best Porb (s)   = %17.15g" % ret[0][1]
print "Best Porb (h)   = %17.15g" % (ret[0][1]/3600.0)
print "Best xorb (s)   = %17.15g" % ret[0][2]
print "Best torb (s)   = %17.15g" % ret[0][3]
print "Best torb (MJD) = %17.15g" % (MJD-ret[0][3]/86400.0)
print "M_c (M_sun)     < %17.15g" % \
      companion_mass_limit(1.4, ret[0][1], ret[0][2])


if (0):
    for ii in xrange(numparts):
        (phases[ii], phaseerrs[ii], b, berr, a) = \
                    measure_phase(profs[ii], fit_template, stdevs[ii], fwhm)
        # phases[i] = old_measure_phase(profs[i], fit_template, fwhm)
        # phaseerrs[i] = 0.0

    plotxy(phases, times, labx='Observation Time',
           laby='Pulse Phase', line=None, symbol=3)
    closeplot()

    # fdcorr is the phase delay given by the f-dot at each time
    fdcorr = 0.5 * fd * times**2.0
    delays = (phases - fdcorr) / f
    toas = (delays + times) / 86400.0 + MJD
    toaerrs = phaseerrs / f * 1000000.0  # in microsec
    print ""

    print "Performing linear fit"
    ret = leastsq(funct, [0.01, 5700.0, 1.0, 0.0], args=(times, delays))
    x = ret[0][0]
    p = ret[0][1]
    t = -ret[0][2] / (2.0 * pi) * p
    #t = ret[0][2] * 180.0 / pi
    print "Linear fit (no weights:"
    print ret[0][0], ret[0][1], ret[0][2], ret[0][3]

    print "raw To = ", t

    t = t / 86400.0 + MJD
    resid = delays - orbeqn(ret[0], times)

    print ret[-1]
    print ""
    print "       P_pulsar = %.14f s" % (1.0 / f)
    print "        P_orbit = %.4f s" % p
    print "     a*sin(i)/c = %.7f s" % x
    print "   Time of peri = %.12f MJD" % t
    print "  mass function = %.4g solar masses" % mass_funct(p, x)
    print " companion mass > %.4f solar masses" % companion_mass_limit(1.4, p, x)
    print "  Residual mean = %.4f ms" % average(resid*1000.0)
    print " Residual stdev = %.4f ms" % standardDeviation(resid*1000.0)
    print ""
    print " Note:  Assuming pulsar mass = 1.4 solar masses"
    print ""
