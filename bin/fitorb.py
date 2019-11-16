#!/usr/bin/env python

"""
fitorb: A non-linear optimizer for solving pulsar orbits by Ryan Lynch
"""
from __future__ import print_function
from builtins import range

from numpy import *
from presto.mpfit import mpfit
from presto.psr_constants import SECPERDAY,TWOPI,DEGTORAD,SOL
from presto import psr_utils
from presto import parfile
from pylab import *
import sys
cspeed = 299792458.0  # m/s


# begin function definitions

def print_usage():
    print("""
A non-linear least-squares optimizer for solving pulsar orbits

Usage: fitorb.py [-p p] [-pb pb] [-x x] [-T0 T0] [-e e] [-w w] [-par par_file] [-nofit const_params] bestprof_files

            -p: Initial guess for pulsar spin period (s; required without -par)
           -pb: Initial guess for orbital period (days; required without -par)
            -x: Initial guess for projected semi-major axis (s; required
                without -par)
           -T0: Initial guess for time of periastron passage (MJD; required
                without -par)
            -e: Initial guess for orbital eccentricity (default = 0)
            -w: Initial guess for longitude of periastron (deg; default = 0)
          -par: A tempo .par file with values of the spin and orbital
                parameters to use as initial guesses (if given, you don't need
                to specify the above parameters)
        -nofit: A comma separated string of parameters to hold constant
               (e.g. "-nofit p,pb,x")
            -o: Root of the output file name(s) (default = "fitorb")
bestprof_files: prepfold .bestprof files containing measurements of p and p-dot
                (and possibly p-ddot)
    """)

    return None


def parse_cmd_line(args):
    """
    Parse command line argumentss
    
    Input
    -----
       args - a list of command line aruments and values
       
    Output
    ------
       user-supplied values for the given arguments
    """

    # create a list of valid command line arguments
    valid_args = ["-p", "-pb", "-x", "-T0", "-e", "-w", "-par", "-nofit","-o"]

    if len(args) == 0:
        print_usage()
        exit(0)
    
    for arg in args:
        # check to make sure all arguments are valid
        if (arg.startswith("-")) and (arg not in valid_args) and \
               not arg.strip("-").replace(".","").isdigit():
            print("ERROR: Unknown arg %s"%arg)
            print_usage()
            exit(0)
    
    # go through the given arguments and store user-supplied values
    try:
        const_params = args.pop(args.index("-nofit")+1)
        args.remove("-nofit")
    except ValueError:
        const_params = ""
        pass
        
    if "-par" in args:
        try:
            par_file_name = args.pop(args.index("-par")+1)
            args.remove("-par")
            par = parfile.psr_par(par_file_name)
            p = par.P0
            pb_days = par.PB
            x = par.A1
            T0 = par.T0
            e = par.E
            w = par.OM
        except IOError:
            print("ERROR: %s not found\n"%par_file_name)
            exit(0)
        except AttributeError:
            print("ERROR: %s does not appear to be a valid binary .par file\n" \
                  %par_file_name)
            exit(0)
    else:        
        try:
            p = float(args.pop(args.index("-p")+1))
            args.remove("-p")
        except ValueError:
            print("ERROR: You must specify a spin period\n")
            exit(0)
        
        try:
            pb_days = float(args.pop(args.index("-pb")+1))
            args.remove("-pb")
        except ValueError:
            print("ERROR: You must specify an orbital period\n")
            exit(0)
        
        try:
            x = float(args.pop(args.index("-x")+1))
            args.remove("-x")
        except ValueError:
            print("ERROR: You must specify a projected semi-major axis\n")
            exit(0)
        
        try:
            T0 = float(args.pop(args.index("-T0")+1))
            args.remove("-T0")
        except ValueError:
            print("ERROR: You must specify a time of periastron passage\n")
            exit(0)
        
        try:
            e = float(args.pop(args.index("-e")+1))
            args.remove("-e")
        except ValueError:
            print("WARNING: Orbital eccentricity not specified, assuming e = 0\n")
            e = 0.0
            const_params = const_params + ",e"
            pass
        
        try:
            w = float(args.pop(args.index("-w")+1))
            args.remove("-w")
        except ValueError:
            print("WARNING: Longitude of periastron not specified, assuming w = 0\n")
            w = 0.0
            const_params = const_params + ",w"
            pass
        
    try:
        out_file_root = args.pop(args.index("-o")+1)
        args.remove("-o")
    except ValueError:
        out_file_root = "fitorb"
        pass
        
    in_files = args
    
    return p,pb_days,x,T0,e,w,const_params,out_file_root,in_files


def read_bestprof(file_name):
    """
    Read relevant information from prepfold .bestprof files (written by
    Scott Ransom
    
    Input
    -----
        file_name - string containing the path to a .bestprof file
    
    Output
    ------
        epoch - the barycentric epoch (MJD) of the observation
        N*dt - length of observation (number of data points * sampling time)
        p0,p1,p2 - observed spin period and higher-order period derivatives
    """
    
    in_file = open(file_name)
    bary = N = 0
    epoch = dt = p0 = p1 = p2 = 0.0
    for line in in_file.readlines():
        if line[0] == "#":
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
                p0 = float(line.split("=")[-1].split("+")[0])/1000.0
                continue
            if ((bary and line.startswith("# P'_bary")) or
                (not bary and line.startswith("# P'_topo"))):
                p1 = float(line.split("=")[-1].split("+")[0])
                continue
            if ((bary and line.startswith("# P''_bary")) or
                (not bary and line.startswith("# P''_topo"))):
                p2 = float(line.split("=")[-1].split("+")[0])
                continue
        else:
            break
    return (epoch, N*dt, p0, p1, p2)

def read_par(pfname,f1errmax=999.0):
    pf = parfile.psr_par(pfname)
    # Try to see how many freq derivs we have
    fs = [pf.F0]
    for ii in range(1, 20):  # hopefully 20 is an upper limit!
        attrib = "F%d"%ii
        if hasattr(pf, attrib):
            fs.append(getattr(pf, attrib))
        else:
            break
    epoch = pf.PEPOCH
    Tobs = (pf.FINISH - pf.START) * 86400.0
    return epoch,Tobs,fs

def get_params_info(params_start, const_params):
    """
    Build a list of dictionaries with information about spin and orbital
    parameters to be passed to mpfit
    
    Input
    -----
        params_start - a list of initial guesses for parameter values
        const_params - a string containing the parameters to hold constant
            during fit

    Output
    ------
        params_info - a list of dictionaries with information on each
            parameter
    """
    
    params_info = []

    # check to see if each parameter should be helt constant
    if "p" in const_params.split(","):
        params_info.append({"parname":"p", # parameter name
                            "value":params_start[0], # initial guess
                            "limited":[True,True], # bounded above and below?
                            # upper and low limits (used if "limited" is "True"
                            "limits":[0.9*params_start[0],1.1*params_start[0]],
                            "fixed":True}) # parameter fixed?
        print("Holding spin period constant")
    else:
        params_info.append({"parname":"p",
                            "value":params_start[0],
                            "limited":[True,True],
                            "limits":[0.9*params_start[0],1.1*params_start[0]],
                            "fixed":False})

    if "pb" in const_params.split(","):
        params_info.append({"parname":"pb",
                            "value":params_start[1],
                            "limited":[True,False],
                            "limits":[0.0,0.0],
                            "fixed":True})
        print("Holding orbital period constant")
    else:
        params_info.append({"parname":"pb",
                            "value":params_start[1],
                            "limited":[True,False],
                            "limits":[0.0,0.0],
                            "fixed":False})

    if "x" in const_params.split(","):
        params_info.append({"parname":"x",
                            "value":params_start[2],
                            "limited":[True,False],
                            "limits":[0.0,0.0],
                            "fixed":True})
        print("Holding projected semi-major axis constant")
    else:
        params_info.append({"parname":"x",
                            "value":params_start[2],
                            "limited":[True,False],
                            "limits":[0.0,0.0],
                            "fixed":False})

    if "T0" in const_params.split(","): 
        params_info.append({"parname":"T0",
                            "value":params_start[3],
                            "limited":[True,True],
                            "limits":[params_start[3] - params_start[1]/SECPERDAY,
                                      params_start[3] + params_start[1]/SECPERDAY],
                            "fixed":True})
        print("Holding time of periastron passage constant")
    else:
        params_info.append({"parname":"T0",
                            "value":params_start[3],
                            "limited":[True,True],
                            "limits":[params_start[3] - params_start[1]/SECPERDAY,
                                      params_start[3] + params_start[1]/SECPERDAY],
                            "fixed":False})

    if "e" in const_params.split(","):
        params_info.append({"parname":"e",
                            "value":params_start[4],
                            "limited":[True,True],
                            "limits":[0.0,1.0],
                            "fixed":True})
        print("Holding eccentricity constant")
    else:
        params_info.append({"parname":"e",
                            "value":params_start[4],
                            "limited":[True,True],
                            "limits":[0.0,1.0],
                            "fixed":False})

    if "w" in const_params.split(","):
        params_info.append({"parname":"w",
                            "value":params_start[5],
                            "limited":[True,True],
                            "limits":[0.0,360.0],
                            "fixed":True})
        print("Holding longitude of periastron constant")
    else:
        params_info.append({"parname":"w",
                            "value":params_start[5],
                            "limited":[True,True],
                            "limits":[0.0,360.0],
                            "fixed":False})

    return params_info


def myasarray(a):
    """
    Properly format array (written by Scott Ransom)
    
    Input
    -----
        a - python array

    Output
    ------
        a - modified python array
    """
    
    if type(a) in [type(1.0),type(1),type(1),type(1j)]:
        a = asarray([a])
    if len(a) == 0:
        a = asarray([a])
    return a


def calc_omega(params, MJD):
    """
    Calculate w in at the barycentric epoch MJD (written by Scott Ransom)

    Input
    -----
        params - a list of parameter values
        MJD - barycentric epoch MJD

    Output
    ------
        w in radians
    """
    
    return params[5]*DEGTORAD


def eccentric_anomaly(params, mean_anomaly):
    """
    Calculate the eccentric anomaly using a simplte iteration to solve
    Kepler's Equations (written by Scott Ransom)

    Input
    -----
        params - a list of parameter values
        mean_anomaly - the mean anomaly

    Output
    ------
        the eccentric anomaly in radians
    """
    ma = fmod(mean_anomaly, TWOPI)
    ma = where(ma < 0.0, ma+TWOPI, ma)
    eccentricity = params[4]
    ecc_anom_old = ma
    ecc_anom = ma + eccentricity*sin(ecc_anom_old)
    # This is a simple iteration to solve Kepler's Equation
    while (maximum.reduce(fabs(ecc_anom-ecc_anom_old)) > 5e-15):
        ecc_anom_old = ecc_anom[:]
        ecc_anom = ma + eccentricity*sin(ecc_anom_old)
    return ecc_anom


def calc_anoms(params, MJD):
    """
    Calculate the mean, eccentric, and true anomalies at the barycentric
    epoch MJD (written by Scott Ransom)

    Input
    -----
        params - a list of parameter values
        MJD - the barycentric epoch MJD

    Output
    ------
        mean_anom - mean anomaly in radians
        ecc_anom - eccentric enomaly in radians
        true_anom - the true anomaly in radians
    """
    MJD = myasarray(MJD)
    difft = (MJD - params[3])*SECPERDAY
    sec_since_peri = fmod(difft, params[1])
    sec_since_peri[sec_since_peri < 0.0] += params[1]
    mean_anom = sec_since_peri/params[1]*TWOPI
    ecc_anom = eccentric_anomaly([params[0],params[1],
                                  params[2],params[3],params[4],
                                  params[5]], mean_anom)
    true_anom = psr_utils.true_anomaly(ecc_anom, params[4])
    return (mean_anom, ecc_anom, true_anom)


def radial_velocity(params, MJD):
    """
    Calculate the radial velocity of the pulsar at the given MJD
    (written by Scott Ransom)

    Input
    -----
        params - a list of parameter values
        MJD - the barycentric epoch MJD

    Output
    ------
        the radial velocity in km/s
    """
    ma, ea, ta = calc_anoms([params[0],params[1],params[2],
                             params[3],params[4],params[5]], MJD)
    ws = calc_omega([params[0],params[1],params[2],
                     params[3],params[4],params[5]], MJD)
    c1 = TWOPI*params[2]/params[1];
    c2 = cos(ws)*sqrt(1-params[4]*params[4]);
    sws = sin(ws);
    cea = cos(ea)
    return SOL/1000.0*c1*(c2*cea - sws*sin(ea)) / (1.0 - params[4]*cea)

def plot_file_panel(in_file,params):
    period = []
    time = []
    if in_file.endswith('.par'):
        (epoch,T,fs) = read_par(in_file)
        for minute in arange(int(T/60.0+0.5)):
            t = minute/1440.0
            time.append(t)
            period.append(1.0/psr_utils.calc_freq(epoch+t,epoch,*fs))
    else:
        (epoch, T, p0, p1, p2) = read_bestprof(in_file)
        for minute in arange(int(T/60.0+0.5)):
            t = minute*60.0
            time.append(minute/1440.0)
            period.append(p0 + t*(p1 + 0.5*t*p2))
    print("Plotting: file,  epoch, Tobs",in_file,epoch,T)
    period = asarray(period)
    time = asarray(time)
    plot(time,period*1000.0,'o')
    xlabel('Time (s)')
    ylabel('Period (ms)')
    title("%.3f" % (epoch,))
    model_time = arange(epoch-0.1, epoch+max(time)+0.1, 0.001)
    plot( model_time-epoch,doppler_period(params, model_time)*1000.0,'r')
    

def doppler_period(params, MJD):
    """
    Calculate the doppler modulated pulse period (written by Scott Ransom)

    Input
    -----
        params - list of parameter values
        MJD - barycentric epoch MJD

    Output
    ------
        observed pulse period in seconds
    """
    vs = radial_velocity([params[0],params[1],params[2],
                          params[3],params[4],params[5]], MJD) \
                          *1000.0 # m/s
    return params[0]*(1.0+vs/SOL)


def funct(params, fjac=None, times=None, measured=None):
    """
    Calculate the difference between the modeled and observed period

    Input
    -----
        params - list of parameter values
        fjac - function for calculating the Jacobian (if None mpfit
               will use a default method)
        times - array of MJDs when period observations were made
        measured - array of observed periods (in seconds)

    Output
    ------
        a list containing the exit status (used by mpfit) and the
        differences between the model and data
    """
    status = 0 # this will probably always be zero
    return [status,doppler_period([params[0],params[1],params[2],
                           params[3],params[4],params[5]],
                          times) - measured]


# parse the command line
p_start,pb_days_start,x_start,T0_start,e_start,w_start,const_params,out_file_root, \
    in_files = parse_cmd_line(sys.argv[1:])

pb_sec_start = pb_days_start * SECPERDAY # need orbital period in seconds

# store user-supplied initial guesses for parameter values
params_start = [p_start, pb_sec_start, x_start, T0_start, e_start, w_start]

# build the dictionary of parameter information
params_info = get_params_info(params_start, const_params)

period = []
time = []
pepochs = [] 
p0s = []
p1s = []
# get the observed periods and times from the .bestprof files
for in_file in in_files:
    if in_file.endswith('.par'):
        (epoch,T,fs) = read_par(in_file)
        if (fs[1] != 0.0):
            p0tmp,p1tmp = psr_utils.p_to_f(fs[0],fs[1])
            p0s.append(p0tmp)
            p1s.append(p1tmp)
            pepochs.append(epoch)
        for minute in arange(int(T/60.0+0.5)):
            t = epoch + minute/1440.0
            time.append(t)
            period.append(1.0/psr_utils.calc_freq(t,epoch,*fs))
    else:
        (epoch, T, p0, p1, p2) = read_bestprof(in_file)
        for minute in arange(int(T/60.0+0.5)):
            t = minute*60.0
            time.append(epoch + minute/1440.0)
            period.append(p0 + t*(p1 + 0.5*t*p2))
            if p1 != 0.0:
                p0s.append(p0)
                p1s.append(p1)
                pepochs.append(epoch)

Torb = min(time)
period = asarray(period)
time = asarray(time)
p0s = asarray(p0s)
p1s = asarray(p1s)
pepochs = asarray(pepochs)
accs = cspeed*p1s/p0s

# Plot ellipse figure
figure(1)
plot(p0s*1000.0,accs,"o")
grid(True)
title('Acceleration vs Period')
xlabel('Period (ms)')
ylabel('Acceleration (m/s^2)')

figure(2)
if len(in_files) < 5:
    nsubrows = 1
    nsubcols = len(in_files)
else:
    nsubrows = (len(in_files)-1)//5 + 1
    nsubcols = 5
ip = 1
for in_file in in_files:
    subplot(nsubrows,nsubcols,ip)
    plot_file_panel(in_file,params_start)
    ip+=1

# do the actual fitting
ret = mpfit(funct, functkw={"times":time, "measured":period}, parinfo=params_info, iterfunct=None)

print("\nmpfit exited with status %i (1--4 is OK)\n"%ret.status)

# print the parameters in a tempo .par file format
print("PEPOCH %17.15g"%epoch)
print("P0  %17.15g"%ret.params[0])
print("BINARY BT")
print("A1  %17.15g"%ret.params[2])
print("E   %17.15g"%ret.params[4])
print("T0  %17.15g"%(ret.params[3]+0.5*ret.params[0]/SECPERDAY))
print("PB  %17.15g"%(ret.params[1]/SECPERDAY))
print("OM  %17.15g"%ret.params[5])
print()

print("Generating plots...")

# make the plots...

amp = TWOPI*abs(ret.params[2])/ret.params[1]  # 2 pi x / Porb 

# plot the model for the full range of observations
if (max(time)-min(time) < 500*(ret.params[1]/SECPERDAY)): # but don't plot too many orbits...
    figure(5)
    model_time = arange(min(time)-0.1, max(time)+0.1, 0.001)
    plot(model_time-model_time[0],doppler_period(ret.params, model_time)*1000.0,'r',alpha=0.5)
    plot(time-model_time[0],period*1000.0,'o')
    xlabel("Days + %.7f"%model_time[0])
    ylabel("Pulsar Period (ms)")

# make more detailed plots around each observation
figure(3)
ip=1
for in_file in in_files:
    subplot(nsubrows,nsubcols,ip)
    plot_file_panel(in_file,ret.params)
    ip+=1

show()
