#!/usr/bin/env python
import sys
import numpy as np
from optparse import OptionParser
from presto.binary_psr import binary_psr
from presto.psr_constants import SOL
from presto.infodata import infodata
import presto.psr_utils as pu
import matplotlib.pyplot as plt

full_usage = """
usage:  binary_info.py [options] parfile(s)
  [-h, --help]   : Display this help
  [-p, --plot]   : Show plot of "v" velocity or "p" period
  [-t, --time]   : Duration of the calculation or observation (sec)
  [-s, --start]  : Barycentric MJD for the start of the calculation/observation
  [-i, --info]   : .inf file name to get barycentric MJD and obs duration

  Read a parfile of a binary pulsar and compute min/max observed barycentric
  spin periods or velocities as either a function of the orbit (default),
  or for a prescribed duration of time. Also show basic information about the
  binary. Optionally plot the velocity/periods/range.

  Copyright Scott Ransom <sransom@nrao.edu>, 2023
"""
usage = "usage: %prog [options] parfile(s)"

def main():
    parser = OptionParser(usage)
    parser.add_option("-p", "--plot", type="string", dest="plot", default=None,
                      help="Plot the results ('v' velocity or 'p' period)")
    parser.add_option("-t", "--time", type="float", dest="duration", default=0.0,
                      help="Duration of the calculation or observation (sec)")
    parser.add_option("-s", "--start", type="float", dest="T_start", default=0.0,
                      help="Barycentric MJD for the start of the calculation/observation")
    parser.add_option("-i", "--info", type="string", dest="inf", default=None,
                      help=".inf file name to get barycentric MJD and obs duration")
    (opts, args) = parser.parse_args()
    if len(args)==0:
        print(full_usage)
        sys.exit(0)

    for pfile in args:
        # Read the parfile(s)
        psr = binary_psr(pfile)
        if hasattr(psr.par, "BINARY"):
            inf = infodata(opts.inf) if opts.inf else None
            if inf:
                opts.T_start = inf.epoch
                opts.duration = inf.dt * inf.N
            if opts.T_start==0: opts.T_start = psr.T0
            if opts.duration==0: opts.duration = psr.PBsec
            T_end = opts.T_start + opts.duration / 86400.0
            norbits = (T_end - opts.T_start) / psr.par.PB
            times = np.linspace(opts.T_start, T_end, max([int(norbits * 1000), 1000]))
            rvs = psr.radial_velocity(times) # km/s
            ps = psr.par.P0*(1.0 + rvs * 1000.0 / SOL) * 1000 # ms
            print("---------------------------------")
            print(f"PSR:          {psr.par.PSR if hasattr(psr.par, 'PSR') else psr.par.PSRJ}")
            print(f"Pspin (ms):   {psr.par.P0 * 1000.0:.6f}")
            print(f"  min (ms):   {ps.min():.6f}")
            print(f"  max (ms):   {ps.max():.6f}")
            print(f" diff (ms):   {ps.max() - ps.min():.6e}")
            print(f"Fspin (Hz):   {1.0/psr.par.P0:.6f}")
            print(f"  min (Hz):   {1000.0 / ps.max():.6f}")
            print(f"  max (Hz):   {1000.0 / ps.min():.6f}")
            print(f" diff (Hz):   {1000.0 / ps.min() - 1000.0 / ps.max():.6e}")
            print(f"Porb (days):  {psr.par.PB:.6f}")
            print(f"Porb (hours): {psr.par.PB * 24:.6f}")
            print(f"Eccentricity: {psr.par.ECC:.6g}")
            print(f"Mcmin (Msun): {pu.companion_mass_limit(psr.par.PB, psr.par.A1):.6f}")
            print(f"Mcmed (Msun): {pu.companion_mass(psr.par.PB, psr.par.A1):.6f}")
            if hasattr(psr.par, "OMDOT"):
                print(f"Mtot (Msun):  {pu.OMDOT_to_Mtot(psr.par.OMDOT, psr.par.PB, psr.par.ECC):.6f}")
            if opts.duration < 0.15 * psr.PBsec:
                pvals = np.polyfit((times - times[0]) * 86400.0, ps / 1000.0, deg=2, full=False)
                # Note:  the pvals[0] * 2 is because there is a 1/2 factor for Taylor expansion
                print(f"\nFold observation with:\nprepfold -p {pvals[2]} -pd {pvals[1]} "
                      f"-pdd {pvals[0] * 2} {inf.basenm if inf else '...'}.dat")
            if opts.plot:
                if opts.plot=="v":
                    vals = rvs
                    ylab = "Radial Velocity (km/s)"
                elif opts.plot=="p":
                    vals = ps - ps.mean()
                    ylab = f"Observed Spin Period (ms) - {ps.mean()} ms)"
                else:
                    print("The only plot choices are 'v' and 'p'.")
                    sys.exit()
                if opts.duration == psr.PBsec:
                    times = np.linspace(0, 1, len(times))
                    xlab = "Orbital Phase"
                elif opts.duration < psr.PBsec:
                    times = np.linspace(0, opts.duration, len(times))
                    xlab = "Seconds"
                else:
                    xlab = "MJD"
                plt.plot(times, vals)
                plt.xlabel(xlab)
                plt.ylabel(ylab)
                plt.show()

if __name__=='__main__':
    main()
    
