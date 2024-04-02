#!/usr/bin/env python

from builtins import str
from builtins import range
from builtins import object
import sys
import os
import numpy as np
import presto.infodata as pi
import presto.presto as pp

scopes = {"gbt": "GB", "arecibo": "AO", "vla": "VL", "parkes": "PK",
    "jodrell": "JB", "gb43m": "G1", "gb 140ft": "G1", "nrao20": "G1",
    "nancay": "NC", "effelsberg": "EF", "srt": "SR", "fast": "FA",
    "wsrt": "WT", "gmrt": "GM", "chime": "CH", "lofar": "LF",
    "lwa": "LW", "mwa": "MW", "meerkat": "MK", "ata": "AT",
    "k7": "K7", "geocenter": "0 "}


def mod_get_baryv(ra, dec, mjd, T, obs="PK", bary=True):
    """
    mod_get_baryv(ra, dec, mjd, T, obs="PK", bary=True):
      Determine the average barycentric velocity towards 'ra', 'dec'
      during an observation from 'obs'.  The RA and DEC are in the
      standard string format (i.e. 'hh:mm:ss.ssss' and 'dd:mm:ss.ssss').
      'T' is in sec and 'mjd' is (of course) in MJD.  The obs variable
      is the standard two character string from TEMPO:  PK, GB, AO, GM, JB, ...
      If bary is true, we will need to back out the topocentric times
    """
    if bary:
        tts = np.linspace(mjd, mjd + T / 86400.0, 200)
    else:
        tts = np.linspace(mjd - 500.0 /86400.0,
                          mjd + (T + 500) / 86400.0, 200)
    nn = len(tts)
    bts = np.zeros(nn, dtype=np.float64)
    vel = np.zeros(nn, dtype=np.float64)
    pp.barycenter(tts, bts, vel, ra, dec, obs, "DE421")
    if bary:
        ind0 = np.fabs(tts - mjd).argmin()
        ind1 = np.fabs(tts - (mjd + T / 86400.0)).argmin()
        # newtts = tts - (bts - tts)
        return vel[ind0:ind1].mean()
    else:
        return vel.mean()


def group_infiles(infilenms):
    """Find the common basenames of files, and sort each by numerical DM, if present

    Parameters
    ----------
    infilenms : list of strings
        These are the input filenames

    Returns
    -------
    sorted list of filenames
    """
    # Make sure that all the files are ".fft" files
    for infilenm in infilenms:
        assert(infilenm.endswith(".fft"))
    # Sort the filenames
    names = sorted(infilenms)
    basenames = []
    DMs = []
    for name in names:
        try:
            ind = name.rfind("_DM")
            if name[:ind] not in basenames:
                basenames.append(name[:ind])
            try:
                dm = float(name[ind+3:-4])
            except ValueError:
                dm = None
        except ValueError:
            if name[:-4] not in basenames:
                basenames.append(name[:-4])
                dm = None
        DMs.append(dm)
    if len(basenames)==1:
        print(f"All files have the same basename '{basenames[0]}'")
    if len(DMs)>1 and None in DMs:
        print("Not all input file names have DM values")
    # Now sort via basename first, then DM, then anything else
    outnames = []
    for basename in basenames:
        tmp = []
        nodms = []
        for ii, name in enumerate(names):
            if name.startswith(basename):
                if DMs[ii] is not None:
                    tmp.append((DMs[ii], name))
                else:
                    nodms.append(name)
        tmp = sorted(tmp) # This sorts by DM, numerically
        for fn in tmp: # These are the files with DMs
            outnames.append(fn[1])
        for fn in nodms: # These are the files without DMs
            outnames.append(fn)
    assert(len(outnames)==len(names))
    return basenames, outnames


def read_birds(birdsname):
    print(f"Reading the birds from '{birdsname}'")
    with open(birdsname, "r") as bfile:
        return bfile.readlines()


def process_birds(birdlines, T, baryv, info):
    psrs = 0
    freqs = 0
    trains = 0
    birds = []
    # PSRs get 40 bins minimum zapped (overkill for most,
    # but required for the _really_ bright ones
    min_psr_width = 40.0 / T
    for line in birdlines:
        bary = 0
        baryfact = 1.0
        line = line[:-1]
        if (len(line)<=3 or line[0]=='#'):
            continue
        elif (line[0]=='P'):
            (tmp, psrname, numharm) = line.split()
            numharm = int(numharm)
            psr = pp.psrepoch(psrname, info.epoch)
            if (psr.orb.p):
                (minv, maxv) = pp.binary_velocity(T, psr.orb)
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
                if info.bary==0:
                    midf /= (1.0 + baryv)
                birds.append((midf, width))
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
                        if info.bary:
                            baryfact = 1.0 if bary else (1.0 + baryv)
                        else:
                            baryfact = (1.0 + baryv) if bary else 1.0
                trains += 1
                if (increase_width):
                    for harm in range(1, numharm+1):
                        birds.append((freq * harm * baryfact, width * harm))
                else:
                    for harm in range(1, numharm+1):
                        birds.append((freq * harm * baryfact, width))
            else:
                freqs += 1
                birds.append((float(words[0]), float(words[1])))
    print("  Read %d freqs, %d pulsars, and %d harmonic series." % \
          (freqs, psrs, trains))
    print("  Total number of birdies = %d" % (len(birds))) 
    return sorted(birds)


def zapfile(fftfile, zaplist, info):
    """Zap the frequencies and widths in zaplist from fftfile

    Parameters
    ----------
    fftfile : file oject
        The .fft file that will be zapped (opened in "rb+" mode)
    zaplist : list of tuples
        List of (freq, width)s (in Hz) to zap
    info : infodata object
        From the .inf file describing the .fft file
    """
    # Use memory-mapping
    ft = np.memmap(fftfile, mode='r+', dtype='complex64')
    T = info.dt * info.N
    for (f, w) in zaplist:
        lor = int(np.round((f - 0.5 * w) * T))
        if lor < 1: lor = 1
        if lor > len(ft): break
        hir = int(np.round((f + 0.5 * w) * T)) + 1
        if hir > len(ft): hir = len(ft)
        # print(lor, hir, lor/T, hir/T)
        # To zap, so that median normalization works, and the Fourier
        # phases are still correct, get the median level just outside
        # the window, and use that as the target level within the
        # zap window.  Remember that median != mean in power spectra
        winlol = int(np.round((f - 2 * w) * T))
        if winlol < 1: winlol = 1
        winhir = int(np.round((f + 2 * w) * T))
        if winhir > len(ft): winhir = len(ft)
        win = np.abs(np.concatenate((ft[winlol:lor], ft[hir:winhir])))
        tgt = np.sqrt(np.median(win**2) / np.log(2)) # sqrt(window mean power)
        # the following sets each zap region aplitude to tgt
        ft[lor:hir] *= tgt / np.abs(ft[lor:hir])
    ft.flush()
    fftfile.close()

if __name__ == '__main__':
    if len(sys.argv)==1:
        print(
    """\nusage:  simple_zapbirds.py .birdsfile .fftfile(s)

  This routine does what makezaplist.py and zapbirds do, but all in 
  one command, and over multiple .fft files.  It also auto-determines
  the barycentric velocity.
  
  The format of the .birds file is a simple text file as shown below.
  Lines starting with '#' are comments and with 'P' are assumed to name a
  pulsar in the ATNF catalog.  The only columns that are required are the 
  first (which specifies a freq, in Hz) and the second, which specifies 
  the width (or, if a pulsar, the number of harmonics zapped).  All others
  are optional.

  The 'grow' flag specifies if the width for each harmonic increases in size.
  That is sometimes useful for some types of RFI or for binary pulsars.
  
  The 'bary' column tells whether the specified freq is barycentric
  or not (i.e. topocentric, like pure, local, RFI tones).
  
  Example .birds file:
#-------------------------------------
# Freq   Width #harm  grow?  bary?
#-------------------------------------
28.760   0.1   3      0      0
60.0     0.05  2      1      0
# Zaps 10 harmonics for PSR J1643-1224:
PSR J1643-1224 10
# Zap 100 Hz with a width of 0.2 Hz
100.0    0.2
""")
    else:
        birds = read_birds(sys.argv[1])
        bases, infilenms = group_infiles(sys.argv[2:])
        lastsize = 0
        lastT = 0
        lastbase = bases[0]
        baryv = 0
        for infilenm in infilenms:
            currsize = os.stat(infilenm).st_size
            with open(infilenm, "rb+") as infile:
                currbase = [x for x in bases if infilenm.startswith(x)][-1]
                if (currsize != lastsize) or (currbase != lastbase):
                    fn, ext = os.path.splitext(infilenm)
                    print(f"Reading file info from '{fn}.inf'")
                    info = pi.infodata(fn+".inf")
                    currT = info.dt * info.N
                    # Only re-compute baryv if we need to
                    if baryv==0 or (currbase != lastbase):
                        baryv = mod_get_baryv(info.RA, info.DEC, info.epoch, currT, 
                                              obs=scopes[info.telescope.lower()],
                                              bary=info.bary)
                    # Only re-compute freqs to zap if the times are also different
                    if (currT != lastT):
                        zaplist = process_birds(birds, currT, baryv, info)
                        # print(zaplist)
                # Now actually do the zapping
                print(f"Zapping '{infilenm}' ... ", end='')
                zapfile(infile, zaplist, info)
                print("done.")
                lastsize = currsize
                lastbase = currbase
                lastT = currT
