#!/usr/bin/env python

from sys import argv, exit
from string import index
from presto import *

class bird:
    def __init__(self, freq, width, bary):
        self.freq = freq
        self.width = width
        self.bary = bary
    def __cmp__(self, other):
        return cmp(self.freq, other.freq)
    def __str__(self):
        if self.bary:
            return "B %20.15g  %20.15g\n" % (self.freq, self.width)
        else:
            return "  %20.15g  %20.15g\n" % (self.freq, self.width)
            
def processbirds(filename):
    # Filename should be the rootfilename of the FFT + '.birds'
    try:
        ii = index(filename, ".birds")
    except ValueError:
        print "\nThe birdie filename must end in '.birds'\n"
        exit(0)
    else:
        rootname = filename[:ii]
    psrs = 0
    freqs = 0
    trains = 0
    birds = []
    print "\nProcessing the birds..."
    info = read_inffile(rootname)
    T = info.dt * info.N;
    # PSRs get 40 bins minimum zapped (overkill for most,
    # but required for the _really_ bright ones
    min_psr_width = 40.0 / T
    file = open(filename, "r")
    for line in file.readlines():
        if (line[0]=='#'):
            continue
        elif (line[0]=='P'):
            (tmp, psrname, numharm) = line.split()
            numharm = int(numharm)
            psr = psrepoch(psrname, info.mjd_i+info.mjd_f)
            if (psr.orb.p):
                (minv, maxv) = binary_velocity(T, psr.orb)
            psrs += 1
            for harm in xrange(1, numharm+1):
                if (psr.orb.p):
                    midv = 0.5 * (maxv + minv)
                    midf = (1.0 + midv) * psr.f * harm
                    width = 1.05 * (maxv - minv) * psr.f * harm
                    if (width < min_psr_width):
                        width = min_psr_width
                else:
                    midf = psr.f * harm
                    width = min_psr_width
                birds.append(bird(midf, width, 1))
        else:
            words = line.split()
            if (len(words)==2):
                freqs += 1
                birds.append(bird(float(words[0]), float(words[1]), 0))
            else:
                freq = float(words[0])
                width = float(words[1])
                numharm = int(words[2])
                trains += 1
                for harm in xrange(1, numharm+1):
                    birds.append(bird(freq * harm, width, 0))
    print "\nRead %d freqs, %d pulsars, and %d harmonic series." % \
          (freqs, psrs, trains)
    print "Total number of birdies = %d" % (len(birds)) 
    birds.sort()
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
    print "\nWrote '%s'\n" % (rootname+".zaplist")

if __name__ == '__main__':
    if len(argv)==1:
        print "\nusage:  makezaplist.py birdsfilename"
        print "       Note:  'birdsfilename' must end in '.birds'\n"
        print "               and a related infofile ('.inf') must exist.\n"
    else:
        processbirds(argv[1])
