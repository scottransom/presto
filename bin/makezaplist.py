#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import sys
from operator import attrgetter
from presto.presto import read_inffile, binary_velocity, psrepoch

class bird(object):
    def __init__(self, freq, width, bary=0):
        self.freq = freq
        self.width = width
        self.bary = bary

    def __str__(self):
        if self.bary:
            return "B %20.15g  %20.15g\n" % (self.freq, self.width)
        else:
            return "  %20.15g  %20.15g\n" % (self.freq, self.width)
            
def processbirds(filename):
    # Filename should be the rootfilename of the FFT + '.birds'
    try:
        ii = filename.index(".birds")
    except ValueError:
        print("\nThe birdie filename must end in '.birds'\n")
        sys.exit(0)
    else:
        rootname = filename[:ii]
    psrs = 0
    freqs = 0
    trains = 0
    birds = []
    print("\nProcessing the birds...")
    info = read_inffile(rootname)
    T = info.dt * info.N;
    # PSRs get 40 bins minimum zapped (overkill for most,
    # but required for the _really_ bright ones
    min_psr_width = 40.0 / T
    file = open(filename, "r")
    for line in file.readlines():
        if (len(line)<=3 or line[0]=='#'):
            continue
        elif (line[0]=='P'):
            (tmp, psrname, numharm) = line.split()
            numharm = int(numharm)
            psr = psrepoch(psrname, info.mjd_i+info.mjd_f)
            if (psr.orb.p):
                (minv, maxv) = binary_velocity(T, psr.orb)
            psrs += 1
            for harm in range(1, numharm+1):
                if (psr.orb.p):
                    midv = 0.5 * (maxv + minv)
                    midf = (1.0 + midv) * psr.f * harm
                    width = (maxv - minv) * psr.f * harm
                    if (0.1 * width < min_psr_width):
                        width = width + min_psr_width
                    else:
                        width = width * 1.1
                else:
                    midf = psr.f * harm
                    width = min_psr_width
                birds.append(bird(midf, width, bary=1))
        else:
            words = line.split()
            increase_width = 0
            bary = 0
            if (len(words) >= 3):
                freq = float(words[0])
                width = float(words[1])
                numharm = int(words[2])
                if (len(words) >= 4):
                    increase_width = int(words[3])
                    if (len(words) >= 5):
                        bary = int(words[4])
                trains += 1
                if (increase_width):
                    for harm in range(1, numharm+1):
                        birds.append(bird(freq * harm, width * harm, bary))
                else:
                    for harm in range(1, numharm+1):
                        birds.append(bird(freq * harm, width, bary))
            else:
                freqs += 1
                birds.append(bird(float(words[0]), float(words[1])))
    print("\nRead %d freqs, %d pulsars, and %d harmonic series." % \
          (freqs, psrs, trains))
    print("Total number of birdies = %d" % (len(birds))) 
    birds.sort(key=attrgetter('freq'))
    file.close()
    file = open(rootname+".zaplist", "w")
    file.write("# This file created automatically with makebirds.py\n")
    file.write("# Lines beginning with '#' are comments\n")
    file.write("# Lines beginning with 'B' are barycentric freqs (i.e. PSR freqs)\n")
    file.write("# %20s  %20s\n" % ("Freq", "Width"))
    file.write("# %s  %s\n" % ("-"*20, "-"*20))
    for birdie in birds:
        file.write(str(birdie))
    file.close()
    print("\nWrote '%s'\n" % (rootname+".zaplist"))

if __name__ == '__main__':
    if len(sys.argv)==1:
        print("\nusage:  makezaplist.py birdsfilename")
        print("       Note:  'birdsfilename' must end in '.birds'\n")
        print("               and a related infofile ('.inf') must exist.\n")
    else:
        processbirds(sys.argv[1])
