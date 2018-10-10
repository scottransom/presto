#!/usr/bin/env python
# Fit an ellipse to a set of measured Periods and Accelerations to get an initial orbit estimate
# This uses the methods of Freire et al. (2001)
# This code written by Paul Ray <paul.ray@nrl.navy.mil
# Inputs are a set of .bestprof files or .par files from which the P0 and P1 (or F0 and F1) values
# and their errors are read.  It can ignore points with too large an F1 error
#
from __future__ import print_function
from numpy import *
from pylab import *
from presto.psr_utils import *
from sys import argv
from presto import parfile
from matplotlib.patches import Ellipse
from scipy.optimize import leastsq
verbose = True
cspeed = 299792458.0 # m/s


def read_bestprof(filename,f1errmax=999.0):
    infile = open(filename)
    bary = N = 0
    epoch = dt = p0 = p1 = p2 = 0.0
    for line in infile.readlines():
        if line[0]=="#":
            if line.startswith("# T_sample"):
                dt = float(line.split("=")[-1])
                continue
            if line.startswith("# Data Folded"):
                N = float(line.split("=")[-1])
                continue
            if line.startswith("# Epoch_topo"):
                try:
                    epochi = float(line.split("=")[-1].split(".")[0])
                    epochf = float("0."+line.split("=")[-1].split(".")[1])
                    epoch = epochi+epochf
                except ValueError:
                    pass
                continue
            if line.startswith("# Epoch_bary"):
                try:
                    epochi = float(line.split("=")[-1].split(".")[0])
                    epochf = float("0."+line.split("=")[-1].split(".")[1])
                    epoch = epochi+epochf
                    bary = 1
                except ValueError:
                    pass
            if ((bary and line.startswith("# P_bary")) or
                (not bary and line.startswith("# P_topo"))):
                valstr = line.split("=")[-1]
                p0 = float(valstr.split("+")[0])/1000.0
                p0err = float(valstr.split("+")[1][3:])/1000.0
                continue
            if ((bary and line.startswith("# P'_bary")) or
                (not bary and line.startswith("# P'_topo"))):
                valstr = line.split("=")[-1]
                p1 = float(valstr.split("+")[0])
                p1err = float(valstr.split("+")[1][3:])
                continue
            if ((bary and line.startswith("# P''_bary")) or
                (not bary and line.startswith("# P''_topo"))):
                p2 = float(line.split("=")[-1].split("+")[0])
                continue
        else:
                break
    f0,f0err,f1,f1err = pferrs(p0,p0err,p1,p1err)
    print("%.4f %10.9g %8.3g  %10.5g %8.3g" % (epoch,f0,f0err,f1,f1err), end=' ')
    if (f1err > f1errmax):
        print(" * Ignored *")
    else:
        print()

    
    #print "----- ",filename
    #print "PEPOCH ",epoch
    #print "F0 ", f0
    #print "F1 ", f1

    return (epoch, N*dt, f0, f0err, f1, f1err)


def read_par(pfname,f1errmax=999.0):
    pf = parfile.psr_par(pfname)
    f0 = pf.F0
    p0 = pf.P0
    try:
        f0err = pf.F0_ERR
    except:
        f0err = 2.0e-5
    if not isfinite(f0err):
        f0err = 3.0e-5
    f1 = pf.F1
    try:
        p1 = pf.P1
    except:
        p1 = 0.0
    try:
        f1err = pf.F1_ERR
    except:
        f1err = 10.0e-8
    mjd = pf.PEPOCH
    if (verbose):
#        print "%6s: %.4f F0 %10.9g +/- %8.03g   F1 %10.5g +/- %8.03g" % (pfname,mjd,f0,f0err,f1,f1err)
        print("%.4f %10.9g %8.3g  %10.5g %8.3g" % (mjd,f0,f0err,f1,f1err), end=' ')
        if (f1err > f1errmax):
            print(" * Ignored *")
        else:
            print()
#        print "          P0 = %g,    P1 = %g" % (p0,p1)

        print("----- ",pfname)
        print("PEPOCH ",mjd)
        print("F0 ", f0)
        print("F1 ", f1)

    return mjd,f0,f0err,f1,f1err


def readPeriodAccelFromPars(parfilelist,f1errmax=3.0e-6):
    mjds = []
    f0s = []
    f0errs = []
    f1s = []
    f1errs = []
    accs = []
    if (verbose):
        print("MJD        F0        F0_err     F1       F1_err")
    for fn in argv[1:]:
        if fn.endswith('.bestprof'):
            mjd,Tobs,f0,f0err,f1,f1err = read_bestprof(fn,f1errmax)
        else:
            mjd,f0,f0err,f1,f1err = read_par(fn,f1errmax)
        mjds.append(mjd)
        f0s.append(f0)
        f0errs.append(f0err)
        f1s.append(f1)
        f1errs.append(f1err)

    del f0,f1,f0err,f1err,mjd

    mjds = array(mjds)
    f0s = array(f0s)
    f0errs = array(f0errs)
    f1s = array(f1s)
    f1errs = array(f1errs)

    # Select only the values where Fdot is measured with some accuracy
    idx = where(f1errs < f1errmax)

    ps, perrs, p1s, p1errs = pferrs(f0s[idx],f0errs[idx],f1s[idx],f1errs[idx])
    selmjds = mjds[idx]

    accs = cspeed*p1s/ps
    accerrs = absolute(accs*sqrt((p1errs/p1s)**2 + (perrs/ps)**2))

    accfreqs = f0s[idx]
    accfreqerrs = f0errs[idx]

    return selmjds, ps, perrs, accs, accerrs


def parabola_funct(pars, x, y_measured, y_err):
    '''Generic parabola fitting function.
    pars is the array of parameters [p0, p1, p2].
    Fit function is y = p2*x**2 + p1*x + p0
    x, y_measured and y_err must all be same length'''
    return ((pars[2]*x**2 + pars[1]*x + pars[0]) - y_measured)/y_err

def funct(pars, ps, Asq_measured, Asq_errs):
    '''Fitting function from Eqn A1 of Freire et al. 2001.
    pars[i] is the array of 3 parameters [a_0, a_1, a_2]
    Asq_measured is the array of measures accelerations SQUARED
    ps is the array of measured pulse periods.'''
    return ((pars[2]*ps**2 + pars[1]*ps + pars[0]) - Asq_measured)/Asq_errs

    
def fitellipse(mjds, accs,accerrs,ps,P0_init,Porb_init,X_init):
    '''Fit an orbit using Eqn A1 of Freire et al. 2001, MNRAS.
    Period errors are assumed to be negligible.'''
    asq = accs**2
    asq_errs = 2*asq*(accerrs/accs)
    apar_init = array([0.0, 0.0, 0.0])
    # Init parameters based on initial orbit, using Eqs A2-A4.
    A1 = 4.0*pi**2*X_init*cspeed/(Porb_init**2)
    P1 = 2.0*pi*X_init*P0_init/Porb_init

    apar_init[2] = -A1**2/(P1**2)
    apar_init[1] = 2.0*P0_init*A1**2/(P1**2)
    apar_init[0] = A1**2 - A1**2*P0_init**2/(P1**2)
    #print "apar init = ",apar_init
    out = leastsq(funct,apar_init,args=(ps,asq,asq_errs),full_output=1)
    apar = out[0]
    covar = out[1]
    P0 = -apar[1]/(2.0*apar[2])
    Porb = (2.0*pi*cspeed)/(P0*sqrt(-apar[2]))
    X = Porb*sqrt(P0**2-apar[0]/apar[2])/(2.0*pi*P0)
    A1 = 4.0*pi**2*X*cspeed/Porb**2
    P1 = 2.0*pi*X*P0/Porb
    
    figure(7)
    errorbar(ps,asq,asq_errs,fmt='o')
    xs = linspace(ps.min(),ps.max(),100)
    plot(xs,apar_init[2]*xs**2 + apar_init[1]*xs + apar_init[0],'b--')
    plot(xs,apar[2]*xs**2 + apar[1]*xs + apar[0],'r--')
    title('Eqn A1 fit')
    ylabel('Acc Squared')
    xlabel('Period')
    grid(1)
    
    return P0, Porb, X, A1, P1

if __name__ == '__main__':

# First read the periods and accelerations from the parfiles
    parfilelist = argv[1:]
    if len(parfilelist)<1:
        print("No par files specified")
        sys.exit(1)
    mjds,ps,perrs,accs,accerrs = readPeriodAccelFromPars(parfilelist,
                                                         f1errmax=3.0e-7)
    print()

    print("MJD :",mjds)
    print("accs :",accs)
    print("accerrs :",accerrs)

# Now setup initial parameter values based on observed periods and accs
    P0_init = ps.mean()
    P1_init = (ps.max()-ps.min())/2.0
    A1_init = abs(accs).max()
    Porb_init = 2.0*pi*cspeed*P1_init/(P0_init*A1_init)
    X_init = P1_init**2*cspeed/(P0_init**2*A1_init)

    vmin = cspeed*(ps.min()/P0_init - 1)
    vmax = cspeed*(ps.max()/P0_init - 1)
    print("vmin = %.2f km/s" % (vmin/1000.0,))
    print("vmax = %.2f km/s" % (vmax/1000.0,))

    print("amin = %.4f m/s^2" % (accs.min(),))
    print("amax = %.4f m/s^2" % (accs.max(),))

    print("pmin = %.6f ms" % (1000.0*ps.min(),))
    print("pmax = %.6f ms" % (1000.0*ps.max(),))


    print("Initial Values:")
    print(" P0 = ",P0_init)
    print(" Porb = %g s (%.3f days)" % (Porb_init,Porb_init/86400.0))
    print(" X   = ",X_init)
    print(" A1  = ",A1_init)
    print(" P1  = ",P1_init)
    print()

# If enough points, do the ellipse fit
    if len(mjds)>=3:
        P0, Porb, X, A1, P1 = fitellipse(mjds,accs,accerrs,ps,P0_init,Porb_init,X_init)

        print("Fitted Values:")
        print(" P0 = ",P0)
        print(" Porb = %g s (%.3f days)" % (Porb, Porb/86400.0))
        print(" X   = ",X)
        print(" A1  = ",A1)
        print(" P1  = ",P1)
        #print "Mcomp,min = ",companion_mass_limit(Porb/86400.0,X)
    else:
        A1 = 0.0
        P0 = 0.0
        P1 = 0.0

# Plot initial and final ellipses with points
    figure(1)
    errorbar(1000.0*ps,accs,xerr=1000.0*perrs,yerr=accerrs,fmt='.')
    grid(True)
    ax = gca()
    # Current initial ellipse is commented out
    #ax.add_artist(Ellipse((P0_init*1000.0,0.0), 2*1000.0*P1_init, 2*A1_init,ls='dashed',color='blue',fc='none',lw=2))
    if (A1>0.0 and P0 > 0.0 and P1 > 0.0):
        ax.add_artist(Ellipse((P0*1000.0,0.0), 2*1000.0*P1, 2*A1,ls='dashed',color='red',fc='none',lw=2))

    title('Acceleration vs. Period')
    xlabel('Period (ms)')
    ylabel('Acceleration (m/s^2)')

    # With Porb and X determined, should look at Tasc
    # This needs to be done...
    if (A1 > 0.0) :
        Porb_days = Porb/86400.0
        phis = arctan(-accs*P1/(A1*(ps-P0)))
        Tascs = mjds - phis*Porb_days/(2.0*pi)
        # Make this do some kind of fit!
        T0 = Tascs[0]
        
        #print "phis = ",phis
        #print "Tascs = ",Tascs
        #print "Tascs in orbits = ",(Tascs-Tascs[0])/Porb_days
        #figure(3)
        #resids = (Tascs-Tascs[0])/Porb_days - floor((Tascs-Tascs[0])/Porb_days)
        #plot(mjds,resids,'s')
        #title('Tasc residuals')
        #xlabel('Time (MJD)')
        #grid(1)

        print()
        print("PAR file of fit: ")

        print("P0 %.15f" % P0)
        print("BINARY BT")
        print("PB %.8f" % Porb_days)
        print("A1 %.6f" % X)
        print("T0 %.6f" % T0)
        print("OM 0.0")
        print("E 0.0")

    show()
