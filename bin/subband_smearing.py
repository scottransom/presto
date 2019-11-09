#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as num
import presto.psr_utils as pu

def subband_smear(DM, subDM, subBW, fctr):
    """
    subband_smear(DM, subDM, subBW, fctr):
        Return the smearing in ms caused by subbanding at DM='DM' given
        subbands of bandwidth 'subBW' (MHz) at DM='subDM'.  All values
        are computed at the frequency fctr in MHz.
    """
    return 1000.0 * pu.dm_smear(num.fabs(DM-subDM), subBW, fctr)

def chan_smear(DM, chanDM, chanBW, fctr):
    """
    chan_smear(DM, chanDM, chanBW, fctr):
        Return the smearing in ms caused by a finite channels at DM='DM'
        given channels of bandwidth 'chanBW' (MHz) at DM='chanDM'.  All
        values are computed at the frequency fctr in MHz.
    """
    return subband_smear(DM, chanDM, chanBW, fctr)

def orig_smear(DM, nchan, chanDM, BW, fctr, dt):
    """
    orig_smear(DM, nchan, chanDM, BW, fctr, dt):
        Return the total smearing in ms due to the sampling rate,
        and the smearing over each channel.
    """
    return num.sqrt((1000.0*dt)**2.0 +
                    chan_smear(DM, chanDM, BW/nchan, fctr)**2.0)

def total_smear(DM, nchan, chanDM, nsub, subDM,
                BW, fctr, dt, downsamp):
    """
    total_smear(DM, nchan, chanDM, nsub, subDM,
                BW, fctr, dt, downsamp):
        Return the total smearing in ms due to the original channel
        format and the properties of the subbands.
    """
    # the factor of two comes from integer-bin shifts when doing
    # the incoherent subbanding
    return num.sqrt(2 * (1000.0*dt*downsamp)**2.0 +
                    chan_smear(DM, chanDM, BW/nchan, fctr)**2.0 +
                    subband_smear(DM, subDM, BW/nsub, fctr)**2.0)

def usage():
    print("""
usage:  subband_smearing.py [options]
  [-l loDM, --loDM=loDM]          : Low DM
  [-h hiDM, --hiDM=HIDM]          : High DM 
  [-t dt, --dt=dt]                : Sample time (s)
  [-s subbands, --nsub=nsub]      : Number of subbands
  [-m subdm, --subDM=subDM]       : DM of each channel
  [-f fctr, --fctr=fctr]          : Center frequency in MHz
  [-b BW, --bw=bandwidth]         : Bandwidth in MHz
  [-n #chan, --nchan=#chan]       : Number of channels
  [-c chanDM, --chanDM=chanDM]    : DM in each channel (default = 0.0)
  [-d N, --downsamp=downsamp]     : Integer downsample (default = 1)

""")    

if __name__=='__main__':
    import getopt, sys

    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:h:t:s:m:f:b:n:c:d:",
                                   ["loDM=", "hiDM=", "dt=",
                                    "nsub=", "subDM="
                                    "fctr=", "bw=",
                                    "nchan=", "chanDM=", "downsamp="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    # Defaults
    chanDM = 0.0
    downsamp = 1

    for o, a in opts:
        if o in ("-l", "--loDM"):
            loDM = float(a)
        elif o in ("-h", "--hiDM"):
            hiDM = float(a)
        elif o in ("-t", "--dt"):
            dt = float(a)
        elif o in ("-s", "--nsub"):
            nsub = int(a)
        elif o in ("-m", "--subDM"):
            subDM = float(a)
        elif o in ("-f", "--fctr"):
            fctr = float(a)
        elif o in ("-b", "--bw"):
            BW = float(a)
        elif o in ("-n", "--nchan"):
            nchan = int(a)
        elif o in ("-c", "--chanDM"):
            chanDM = float(a)
        elif o in ("-d", "--downsamp"):
            downsamp = float(a)

    DMs = num.linspace(loDM, hiDM, 1000)

    samp = num.ones_like(DMs) * 1000.0 * dt
    dsamp = samp * downsamp
    chan = chan_smear(DMs, chanDM, BW/nchan, fctr)
    subband = subband_smear(DMs, subDM, BW/nsub, fctr)
    orig = orig_smear(DMs, nchan, chanDM, BW, fctr, dt)
    total = total_smear(DMs, nchan, chanDM, nsub, subDM,
                        BW, fctr, dt, downsamp)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(DMs, samp, 'g:',
                DMs, dsamp, 'g--',
                DMs, chan, 'r:',
                DMs, subband, 'r--',
                DMs, orig, 'k:',
                DMs, total, 'k')
    leg = ax.legend(('Sampling time', 'Downsampling',
                     'Channel smear', 'Subband smear',
                     'Original time res', 'Total time res'),
                    loc='upper center')
    ax.set_xlabel('Disperson Measure')
    ax.set_ylabel('Smearing (ms)')
    ax.set_xlim([DMs.min(), DMs.max()])
    ax.set_ylim([0.5*1000.0*dt, 2.0*total.max()])
    plt.show()
    
