#!/usr/bin/env python
import numpy as np
import sys

def read_weights(filename):
    """
    reat_weights(filename):
        Read the channels and weights from a .weights text file.
    """
    chans, weights = np.loadtxt(filename, unpack=True)
    return chans.astype(np.int), weights.astype(np.float)

def build_chanline(weights):
    """
    build_chanline(weights):
        Build an ignorechans line from a .weights file written by hand or rfifind.py
            The ignorechans line can be used to ignore channels with prepdata, prepfold, 
            prepsubband, rfifind, and mpiprepsubband.
    """
    N = len(weights)
    outstr = ''
    if N==0:
        return outstr
    lastbad = -2
    ii = 0
    while ii < len(weights):
        if weights[ii] == 0.0:
            if lastbad == ii-1:
                outstr += '' if outstr[-1]==':' else ':'
            else:
                outstr += str(ii)
            lastbad = ii
        else:
            if lastbad == ii-1:
                outstr += str(ii-1)+',' if outstr[-1]==':' else ','
        ii += 1
    if outstr[-1]==':':
        outstr += str(N-1)
    return outstr if outstr[-1]!=',' else outstr[:-1]

def build_pazline(chanline):
    """
    build_pazline(chanline):
        Build a command line for paz from an ignorechan line
    """
    outstr = 'paz -e zap '
    for part in chanline.split(','):
        outstr += "-Z '"+part.replace(":", " ")+"' " if ":" \
            in part else "-z "+part+" "
    return outstr

if __name__=="__main__":

    if len(sys.argv) != 2:
        print("\nusage: {} file\n".format(sys.argv[0]))
        sys.exit(1)

    # Read the channels and weights
    chans, weights = read_weights(sys.argv[1])

    # Get the chanline
    chanline = build_chanline(weights)
    print(chanline)

    # Convert it to a paz command
    pazline = build_pazline(chanline)
    print("\n"+pazline)
