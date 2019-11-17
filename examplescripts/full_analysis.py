#!/usr/bin/python
from os import system, chdir, remove
from sys import stdout, argv, exit
from glob import glob
from optparse import OptionParser
from fcntl import *

def myexecute(cmd):
    stdout.write("\n'"+cmd+"'\n")
    stdout.flush()
    system(cmd)

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-d", "--dualcpus", action="store_true", dest="dualcpus", default=0,
                      help="There are two CPUs working on the same data files.")
    parser.add_option("-f", "--fft", action="store_true", dest="fft", default=0,
                      help="Calculate (and later remove) the .fft files")
    parser.add_option("-v", "--baryvel", type="float", dest="baryv", default=0.0,
                      help="Barycentric velocity in units of c")
    parser.add_option("-o", "--outdir", type="string", dest="outdir", default=".",
                      help="Output directory to store results")
    parser.add_option("-w", "--workdir", type="string", dest="workdir", default=".",
                      help="Working directory for search")
    parser.add_option("-l", "--flo", type="float", dest="flo", default=1.0,
                      help="Low frequency (Hz) to search")
    parser.add_option("-x", "--fhi", type="float", dest="fhi", default=10000.0,
                      help="High frequency (Hz) to search")
    parser.add_option("-z", "--zmax", type="int", dest="zmax", default=170,
                      help="Maximum fourier drift (bins) to search")
    parser.add_option("-a", "--numharm", type="int", dest="numharm", default=8,
                      help="Number of harmonics to sum when searching")
    parser.add_option("-s", "--sigma", type="float", dest="sigma", default=2.0,
                      help="Cutoff sigma to consider a candidate")
    parser.add_option("-p", "--pmax", type="int", dest="pmax", default=6,
                      help="Maximum # of harmonics to sum in sideband search")
    (options, args) = parser.parse_args()
    if (options.outdir[-1]!="/"):
        options.outdir = options.outdir+"/"
    if (options.workdir!='.'):
        chdir(options.workdir)

    # Get the datafiles and determine their DMs from their names
    datanames = glob('*.dat')
    if (len(datanames)==0):
        exit(0)

    dms = []
    for dataname in datanames:
        loptr = dataname.find("_DM")+3
        hiptr = dataname.find(".dat")
        dms.append(float(dataname[loptr:hiptr]))
    dms.sort()

    # If using 2 CPUs, determine if we will analyze the first
    # half of the data files or the second half
    if (options.dualcpus):
        flag = open("PROC_firsthalf", "w+")
        try: lockf(flag.fileno(), LOCK_EX | LOCK_NB)
        except IOError:     # Can't get the lock, therefore 2nd
            flag.close()
            flag = open("PROC_lasthalf", "w+")
            firsthalf = 0
            dms = dms[len(dms)/2:]
        else:               # Got the lock, therefore 1st
            firsthalf = 1
            dms = dms[:len(dms)/2]

    # The basename of the data files
    basename = datanames[0][:loptr-3]

    # Get the bird file (the first birdie file in the directory!)
    birdname = glob("*.birds")
    if birdname:
        birdname = birdname[0]

    for dm in dms:
        filenamebase = basename+'_DM%.2f'%dm
        outnamebase = options.outdir+filenamebase
        if options.fft:
            myexecute('realfft '+filenamebase+'.dat')
        myexecute('cp '+birdname+' '+filenamebase+'.birds')
        myexecute('makezaplist.py '+filenamebase+'.birds')
        myexecute('rm '+filenamebase+'.birds')
        myexecute('zapbirds -zap -zapfile '+filenamebase+
                  '.zaplist -baryv %g '%options.baryv+filenamebase+'.fft')
        myexecute('rm '+filenamebase+'.zaplist')
        myexecute('search_bin -flo 80 -ncand 200 -harmsum 1 '+filenamebase+'.fft')
        myexecute('search_bin -flo 80 -ncand 200 -harmsum %d '%options.pmax+filenamebase+'.fft')
        myexecute('cp '+filenamebase+'_bin* '+options.outdir)
        myexecute('accelsearch -sigma %f -zmax 4 -numharm %d -flo %f -fhi %f ' % \
                  (options.sigma, options.numharm, options.flo, options.fhi)+filenamebase+'.fft')
        myexecute('cp '+filenamebase+'_ACCEL_4 '+options.outdir)
        myexecute('cp '+filenamebase+'_ACCEL_4.cand '+options.outdir)
        myexecute('cp '+filenamebase+'.inf '+options.outdir)
        myexecute('accelsearch -sigma %f -zmax %d -numharm %d -flo %f -fhi %f ' % \
                  (options.sigma, options.zmax, options.numharm, options.flo, options.fhi)+filenamebase+'.fft')
        myexecute('cp '+filenamebase+'_ACCEL_%d '%options.zmax+options.outdir)
        myexecute('cp '+filenamebase+'_ACCEL_%d.cand '%options.zmax+options.outdir)
        myexecute('single_pulse_search.py -p '+filenamebase+'.dat')
        myexecute('cp '+filenamebase+'.singlepulse '+options.outdir)
        if options.fft:
            myexecute('rm '+filenamebase+'.fft')

    if (options.dualcpus):
        if firsthalf:
            lockf(flag.fileno(), LOCK_UN)
            flag.close()
            remove("PROC_firsthalf")
        else:
            flag.close()
            remove("PROC_lasthalf")

if __name__ == "__main__":
      main()

