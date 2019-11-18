from __future__ import print_function
from builtins import zip
from builtins import range
import os

# To use this script to help you dedisperse a bunch of time series, first
# run DDplan.py with appropriate values for your data to generate a
# dedispersion plan:
#
# sransom@foops:~$ DDplan.py -d 200 -t 0.000072 -s 32 -n 96 -b 48.0 -f 820.0
#
# Minimum total smearing     : 0.102 ms
# --------------------------------------------
# Minimum channel smearing   : 3.76e-05 ms
# Minimum smearing across BW : 0.00361 ms
# Minimum sample time        : 0.072 ms
#
# Setting the new 'best' resolution to : 0.072 ms
# Best guess for optimal initial dDM is 0.199
#
#   Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract
#     0.000     38.400    0.20       1    4.80     192      24       8    0.7273
#    38.400     60.000    0.30       2    7.20      72      24       3    0.1364
#    60.000    108.000    0.50       4   12.00      96      24       4    0.09091
#   108.000    204.000    1.00       8   24.00      96      24       4    0.04545
#
#
# Now with that plan, fill in the lists below and appropriate variables
# for your data and you can then generate the subbands and time series
# using "python dedisp.py"
#


def myexecute(cmd):
    print("'%s'"%cmd)
    os.system(cmd)


# dDM steps from DDplan.py
dDMs      = [0.2, 0.3, 0.5, 1.0]
# dsubDM steps
dsubDMs   = [4.8, 7.2, 12.0, 24.0]
# downsample factors
downsamps = [1, 2, 4, 8]
# number of calls per set of subbands
subcalls  = [8, 3, 4, 4]
# The low DM for each set of DMs 
startDMs  = [0.0, 38.4, 60.0, 108.0]
# DMs/call
dmspercall = 24
# Number of subbands
nsub = 32
# The number of points in the least-downsampled time series
numout = 500000
# The basename of the output files you want to use
basename = "S0062004000"
# The name of the raw data file (or files if you use wildcards) to use
rawfiles = basename+"*.bcpm2"
# The name of the maskfile to apply (if no mask, use None)
maskfile = basename+"_rfifind.mask"

# Loop over the DDplan plans
for dDM, dsubDM, downsamp, subcall, startDM in \
        zip(dDMs, dsubDMs, downsamps, subcalls, startDMs):
    # Get our downsampling right
    subdownsamp = downsamp/2
    datdownsamp = 2
    if downsamp < 2: subdownsamp = datdownsamp = 1
    # Loop over the number of calls
    for ii in range(subcall):
        subDM = startDM + (ii+0.5)*dsubDM
        # First create the subbands
        if maskfile:
            myexecute("prepsubband -mask %s -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                      (maskfile, subDM, nsub, subdownsamp, basename, rawfiles))
        else:
            myexecute("prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                      (subDM, nsub, subdownsamp, basename, rawfiles))
        # And now create the time series
        loDM = startDM + ii*dsubDM
        subnames = basename+"_DM%.2f.sub[0-9]*"%subDM
        myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                  (loDM, dDM, dmspercall, datdownsamp, basename, subnames))
