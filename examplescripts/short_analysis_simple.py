#!/usr/bin/python
from os import system, chdir
from sys import stdout, argv, exit
from glob import glob
from optparse import OptionParser
from presto.presto import read_inffile, writeinf, get_baryv
from presto import infodata


def myexecute(cmd):
    stdout.write("\n'"+cmd+"'\n")
    stdout.flush()
    system(cmd)


def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-n", "--number", type="int", dest="nM", default=40,
                      help="Number of points in each chunk (millions)")
    parser.add_option("-o", "--outdir", type="string", dest="outdir", default=".",
                      help="Output directory to store results")
    parser.add_option("-d", "--workdir", type="string", dest="workdir", default=".",
                      help="Working directory for search")
    parser.add_option("-l", "--flo", type="float", dest="flo", default=10.0,
                      help="Low frequency (Hz) to search")
    parser.add_option("-f", "--frac", type="float", dest="frac", default=0.5,
                      help="Fraction to overlap")
    parser.add_option("-x", "--fhi", type="float", dest="fhi", default=10000.0,
                      help="High frequency (Hz) to search")
    parser.add_option("-z", "--zmax", type="int", dest="zmax", default=160,
                      help="Maximum fourier drift (bins) to search")
    parser.add_option("-w", "--wmax", type="int", dest="wmax", default=0,
                      help="Maximum fourier drift deriv (bins) to search")
    parser.add_option("-a", "--numharm", type="int", dest="numharm", default=4,
                      help="Number of harmonics to sum when searching")
    parser.add_option("-s", "--sigma", type="float", dest="sigma", default=2.0,
                      help="Cutoff sigma to consider a candidate")
    (options, args) = parser.parse_args()
    if options.outdir[-1]!= "/":
        options.outdir = options.outdir+"/"
    if options.workdir!= '.':
        chdir(options.workdir)
    if options.nM >= 1000000:
        if options.nM % 1000000:
            print("If you specify --num nM to be > 1000000, it must be divisible by 1000000.")
            exit(1)
    else:
        options.nM *= 1000000
    short_nM = options.nM // 1000000

    # The basename of the data files
    if argv[1].endswith(".dat"):
        basename = "../"+argv[1][:-4]
    else:
        basename = "../"+argv[1]

    # Get the bird file (the first birdie file in the directory!)
    birdname = glob("../*.birds")
    if birdname:
        birdname = birdname[0]

    outnamebase = options.outdir+basename[3:]
    inf = read_inffile(basename)
    idata = infodata.infodata(basename+".inf")
    N = inf.N
    t0i = inf.mjd_i
    t0f = inf.mjd_f
    num = 0
    point = 0
    T = options.nM * inf.dt / 86400.0
    baryv = get_baryv(idata.RA, idata.DEC, idata.epoch, T, obs='GB')
    print("Baryv = ", baryv)
    inf.N = options.nM
    inf.numonoff = 0
    nM = options.nM // 1000000
    while point + options.nM < N:
        pM = point // 1000000
        outname = basename[3:]+'_%03dM'%nM+'_%02d'%num
        stdout.write('\n'+outname+'\n\n')
        inf.name = outname
        tstartf = inf.mjd_f + num * T * options.frac
        if tstartf > 1.0:
            tstartf = tstartf - 1.0
            inf.mjd_i = inf.mjd_i + 1
        inf.mjd_f = tstartf
        writeinf(inf)
        myexecute('dd if=' + basename +'.dat of=' + outname +'.dat bs=4000000 skip=' +
                  repr(pM) + ' count=' + repr(nM))
        myexecute('realfft '+outname+'.dat')
        myexecute('rm -f '+outname+'.dat')
        myexecute('cp '+birdname+' '+outname+'.birds')
        myexecute('makezaplist.py '+outname+'.birds')
        myexecute('rm -f '+outname+'.birds')
        myexecute('zapbirds -zap -zapfile '+outname+'.zaplist -baryv %g '%
                  baryv+outname+'.fft')
        myexecute('rm -f '+outname+'.zaplist')
        if options.wmax > 0:
            myexecute('accelsearch -sigma %.2f -zmax %d -wmax %d -numharm %d -flo %f -fhi %f '%
                      (options.sigma, options.zmax, options.wmax,
                       options.numharm, options.flo, options.fhi)+outname+'.fft')
            myexecute('rm '+outname+'.fft '+outname+'_JERK_%d.txtcand'%options.wmax)
            myexecute('cp '+outname+'_JERK_%d '%options.wmax + options.outdir)
            myexecute('cp '+outname+'_JERK_%d.cand '%options.wmax + options.outdir)
        else:
            myexecute('accelsearch -sigma %.2f -zmax %d -numharm %d -flo %f -fhi %f '%
                      (options.sigma, options.zmax,
                       options.numharm, options.flo, options.fhi)+outname+'.fft')
            myexecute('rm '+outname+'.fft '+outname+'_ACCEL_%d.txtcand'%options.zmax)
            myexecute('cp '+outname+'_ACCEL_%d '%options.zmax + options.outdir)
            myexecute('cp '+outname+'_ACCEL_%d.cand '%options.zmax + options.outdir)
        myexecute('cp '+outname+'.inf '+options.outdir)
        num = num + 1
        point = point + int(options.nM * options.frac)


if __name__ == "__main__":
      main()

