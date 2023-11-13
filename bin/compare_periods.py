#!/usr/bin/env python
from presto.parfile import psr_par
import sys, getopt

def usage():
    sys.stderr.write("""
usage:  compare_periods.py [options] period(ms) parfiles
  [-h, --help]                      : Display this help
  [-t thresh, --threshold=thresh]   : Fractional error to flag a potential clash (def=1e-4)
  [-n harm, --numharm=harm]         : Highest harmonic (or fractional harmonic) to check

""")

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:n:t:",
                                   ["help", "numharm=", "threshold="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)<=2:
        usage()
        sys.exit(2)
    thresh = 0.0001
    harm = 10
    p1 = float(args[0]) * 1e-3 # now in seconds
    parfiles = args[1:]
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-t", "--threshold"):
            thresh = float(a)
        if o in ("-n", "--numharm"):
            harm = int(a)

    for parfile in parfiles:
        psr = psr_par(parfile)
        name = psr.PSR if hasattr(psr, "PSR") else psr.PSRJ
        p2 = 1.0/psr.F0 if hasattr(psr, "F0") else psr.P0
        for num in range(1, harm+1):
            for denom in range(1, harm+1):
                if num==denom and num > 1:
                    continue
                tfract = num/denom
                ofract = p1/p2
                ferr = abs(tfract-ofract)/tfract
                if ferr < thresh:
                    print(f"{name}:  {num}/{denom} = {p1*1e3:10.6f} / {p2*1e3:10.6f} (error = {ferr:.4g})")

