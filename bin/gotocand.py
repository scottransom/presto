#!/usr/bin/env python
import sys
import os
import os.path
import glob
import string
import re
from subprocess import Popen, PIPE, STDOUT
from presto.presto import fourierprops, get_rzw_cand

short_re = re.compile("_\d\d\dM_\d\d_ACCEL_")

def determine_dt(candfile):
    for line in open(candfile):
        if line.startswith(b" Width of each time series bin"):
            return float(line.split()[-1])

def short_stuff(candfile, candnum, shortinfo, nodename, datfile):
    tmp = shortinfo[0].split("_")
    ii = int(tmp[2])
    searchN = 1000000 * int(tmp[1][:-1])
    fileN = get_datfile_len(nodename, datfile)
    start = ii * float(searchN/2) / fileN
    end = start + float(searchN) / fileN
    dt = determine_dt(candfile)
    chunkT = dt * searchN
    obsT = dt * fileN
    cand = fourierprops()
    get_rzw_cand(candfile+'.cand', candnum, cand)
    # fourier props file reports average r and average z.
    # We need the starting values for this chunk.
    z0 = cand.z - 0.5 * cand.w
    r0 = cand.r - 0.5 * z0 - cand.w / 6.0
    f = r0 / chunkT
    fd = z0 / chunkT**2
    fdd = cand.w / chunkT**3
    return (" -start %.3f -end %.3f "%(start, end),
            "_%.3f-%.3f"%(start, end), f, fd, fdd)

def get_dm(filenm):
    parts = filenm.split("_")
    for part in parts:
        if part[:2]=='DM':
            return part[2:]
    return None

def get_basename(filenm):
    offset = filenm.find("_DM")
    if offset > 0:
        return filenm[:offset]
    else:
        return None

def find_node(DM):
    nodefiles = glob.glob("node*")
    for nodefile in nodefiles:
        for line in open(nodefile):
            if line[:4]=="nimr":
                if DM in line.split():
                    return line.split()[0]
    return None

def find_local_datfile(basename, DM):
    p = Popen("find .. -name \*%s\*DM%s\*dat"%(basename, DM), shell=True,
              bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    (i, o) = (p.stdin, p.stdout)
    datfile = ''
    for line in o:
        line = line.strip()
        if line.startswith(b"find:"):
            line = line.join(line.split()[1:])
        if line.endswith(b".dat"):
            datfile = line.decode("utf-8")
    print("'%s'"%datfile)
    if datfile!='':
        return datfile

def find_datfile(nodename, basename, DM):
    p = Popen("ssh %s find -L /scratch -name \*%s\*DM%s\*dat"%(nodename, basename, DM),
              shell=True, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    (i, o) = (p.stdin, p.stdout)
    datfile = ''
    for line in o:
        line = line.strip()
        if line.startswith(b"find:"):
            line = line.join(line.split()[1:])
        if line.endswith(b".dat"):
            datfile = line.decode("utf-8")
    print("'%s'"%datfile)
    if datfile!='' and datfile.startswith(b"/scratch"):
        return datfile
    return None

def get_datfile_len(nodename, datfile):
    if nodename:
        p = Popen("ssh %s ls -l %s | awk '{ print $5 };'"%(nodename, datfile),
                  shell=True, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=STDOUT,
                  close_fds=True)
        (i, o) = (p.stdin, p.stdout)
    else:
        p = Popen("ls -l %s | awk '{ print $5 };'"%(datfile),
                  shell=True, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=STDOUT,
                  close_fds=True)
        (i, o) = (p.stdin, p.stdout)
    filelen = o.readline().decode("utf-8")
    if filelen!='':
        return int(filelen)/4
    return None

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print("\nusage: gotocand.py [-local] candfile:candnum\n")
        sys.exit(0)
    
    local = 0
    if (sys.argv[1]=="-local"):
        local = 1
        sys.argv.pop(1)
    outdir = os.getcwd()
    if (len(sys.argv) > 2):
        extraargs = "".join(sys.argv[2:])
    else:
        extraargs = ""
    candfile, candnum = sys.argv[1].split(':')
        
    dm = get_dm(candfile)
    if dm is None:
        print("Error:  Could not find a DM value in '%s'!"%candfile)
        sys.exit(0)

    base = get_basename(candfile)
    if base is None:
        print("Error:  Could not find the base filename in '%s'!"%candfile)
        sys.exit(0)

    # Is the candidate from a short-chunk search?
    shortcand = short_re.findall(candfile)

    if (local):
        node = None
        datfile = find_local_datfile(base, dm)
    else:
        node = find_node(dm)
        if node is None:
            print("Error:  Could not find the node where the dat file should be!")
            sys.exit(0)

        datfile = find_datfile(node, base, dm)
        if datfile is None:
            print("Error:  Could not find .dat file on the node!")
            sys.exit(0)

    fullcandfile = os.path.join(outdir, candfile)+".cand"
    outfile = base+"_DM%s"%dm
    datfiledir, datfilenm = os.path.split(datfile)
    
    if not local:
        print("\nGoing to %s and folding candidate #%s from the file %s."%\
              (node,candnum,candfile))
    print("  Folding command:")

    if shortcand:
        shortparts, shortoutext, f, fd, fdd = short_stuff(candfile, int(candnum),
                                                          shortcand, node, datfile)
        extraargs += shortparts
        outfile += shortoutext
        foldcommand = "prepfold %s -f %.15g -fd %.15g -fdd %.15g -o %s %s"%\
                      (extraargs, f, fd, fdd, outfile, datfile)
        print(foldcommand)
        if not local:
            os.system("ssh -X %s 'cd %s ; %s'"%(node, datfiledir, foldcommand))
            os.system("scp -c blowfish %s:%s*_%.2f*.pfd* %s"% \
                      (node, os.path.join(datfiledir, outfile), f, outdir))
        else:
            os.system("%s"%(foldcommand))
    else:
        foldcommand = "prepfold %s -accelcand %s -accelfile %s -o %s %s"%\
                      (extraargs, candnum, fullcandfile, outfile, datfile)

        print("    %s"%foldcommand)
        if not local:
            os.system("ssh -X %s 'cd %s ; %s'"%(node, datfiledir, foldcommand))
            os.system("scp -c blowfish %s:%s*ACCEL_Cand_%d*.pfd* %s"% \
                      (node, os.path.join(datfiledir, outfile), int(candnum), outdir))
        else:
            os.system("%s"%(foldcommand))
