#!/bin/env python
import numpy as num
import sys

if len(sys.argv)==1:
    print("usage:  combine_weights.py INPUT_rfifind.weight_FILES")
    sys.exit()

wgts = None
for name in sys.argv[1:]:
    tmpwgts = num.loadtxt(name, dtype=num.int32, comments="#",
                          usecols=(1,), unpack=True)
    if wgts is None:
        wgts = tmpwgts
        print("'%s' had %d bad channels" % \
              (name, num.equal(wgts, 0).sum()))
        chans = num.arange(len(wgts))
    else:
        gain = (tmpwgts-wgts) < 0
        print("'%s' gained %d chans:" % (name, gain.sum()), chans[gain])
        loss = (tmpwgts-wgts) > 0
        print("'%s'   lost %d chans:" % (name, loss.sum()), chans[loss])
        wgts = num.logical_and(wgts, tmpwgts)

print("There are %d channels total." % num.equal(wgts, 0).sum())
print("Writing them to 'combined.weights'")

outfile = open("combined.weights", 'w')
for ii, wgt in enumerate(wgts):
    outfile.write("%4d %d\n" % (ii, wgt))
outfile.close()

        
    
