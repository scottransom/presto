#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import sys, re

maxwappnum = 7
wappfiles = {}
for filename in sys.argv[1:]:
    for wappnum in range(1, maxwappnum + 1):
        if ((wappnum == 1 and re.search("\.wapp\.", filename)) or \
                (wappnum > 1 and re.search("\.wapp%d?\." % wappnum, filename))):
            if wappnum in wappfiles:
                wappfiles[wappnum].append(filename)
            else:
                wappfiles[wappnum] = [filename]
            break

for key in list(wappfiles.keys()):
    numfiles = len(wappfiles[key])
    wappfiles[key].sort()

for filenum in range(numfiles):
    for wappnum in range(1, maxwappnum + 1):
        if wappnum in wappfiles:
            print(wappfiles[wappnum][filenum], end=' ')
