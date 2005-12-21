#!/usr/bin/env python
import bisect, os, sys, getopt
import scipy, scipy.signal, umath, ppgplot
from Numeric import *
from presto import read_inffile
from psr_utils import coord_to_string
from TableIO import readColumns

class candidate:
    def __init__(self, DM, sigma, time, bin, downfact):
        self.DM = DM
        self.sigma = sigma
        self.time = time
        self.bin = bin
        self.downfact = downfact
    def __str__(self):
        return "%7.2f %7.2f %13.6f %10d     %3d\n"%\
               (self.DM, self.sigma, self.time, self.bin, self.downfact)
    def __cmp__(self, other):
        return cmp(self.bin, other.bin)

def cmp_sigma(self, other):
    retval = -cmp(self.sigma, other.sigma)
    return retval

def usage():
    print """
usage:  single_pulse_search.py [options] .dat files _or_ .singlepulse files
  [-h, --help]        : Display this help
  [-t, --threshold]   : Set a different threshold SNR (default=5.0)
  [-m, --maxfact]     : Set the max downsample factor (see below for default)
  [-x, --xwin]        : Don't make a postscript plot, just use an X-window

  Perform a single-pulse search (or simply re-plot the results of a
  single-pulse search) on a set of de-dispersed time series (.dat
  files).

  The search attempts to find pulses by matched-filtering the data
  with a series of different width boxcar functions.  The possible
  boxcar sizes are [1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150]
  bins.  By default the the boxcars <= 30 are used.  You can specify
  that the larger boxcars are used with the -m (or --maxfact) option.

  The matched filtering (and accounting for all the possible 'phase'
  offsets of each boxcar) is accomplished by convolving the boxcars
  with the full resolution data.  'Duplicate' candidates from this
  process are filtered, leaving only the most significant.  The time
  series are initially smoothed using a piecewise linear fit to the
  data where each piece is 2000 data points long.

  If the input files are .singlepulse files, we won't actually perform
  a search, we'll only read in the output .singlepulse files and make
  a plot using the information they contain (along with the
  corresponding .inf files).

  Copyright Scott Ransom <sransom@nrao.edu>, 2005
"""
    
def read_singlepulse_files(infiles, threshold):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    for ii, infile in enumerate(infiles):
        filenmbase = infile.rstrip(".singlepulse")
        info = read_inffile(filenmbase)
        DMstr = "%.2f"%info.dm
        DMs.append(info.dm)
        num_v_DMstr[DMstr] = 0
        if ii==0:
            info0 = info
        if os.stat(infile)[6]:
            cands = transpose(readColumns(infile, "#"))
            for cand in cands:
                if cand[1] >= threshold:
                    candlist.append(candidate(*cand))
                    num_v_DMstr[DMstr] += 1
    DMs.sort()
    return info0, DMs, candlist, num_v_DMstr

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hxm:t:",
                                   ["help", "xwin", "maxfact=", "threshold="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    if len(sys.argv)==1:
        usage()
        sys.exit(2)
    dosearch = True
    threshold = 5.0
    chunklen = 10000
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150]
    max_downfact = 30
    pgplot_device = ""
    for ii, (o, a) in enumerate(opts):
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-m", "--maxfact"):
            max_downfact = int(a)
        if o in ("-t", "--threshold"):
            threshold = float(a)
        if o in ("-x", "--xwin"):
            pgplot_device = "/XWIN"

    downfacts = [x for x in default_downfacts if x <= max_downfact]
    break_points = chunklen/5 * arange(1,5)
    overlap = max(downfacts)/2
    float_len = 4
    if args[0].endswith(".singlepulse"):
        filenmbase = args[0].rstrip(".singlepulse")
        dosearch = False
    else:
        filenmbase = args[0].rstrip(".dat")

    # Don't do a search, just read results and plot
    if not dosearch:
        info, DMs, candlist, num_v_DMstr = \
              read_singlepulse_files(args, threshold)
        orig_N, orig_dt = int(info.N), info.dt
        obstime = orig_N * orig_dt
    else:
        DMs = []
        candlist = []
        num_v_DMstr = {}

        # Loop over the input files
        for filenm in args:
            filenmbase = filenm.rstrip(".dat")
            fileblocklen = chunklen + 2*overlap
            info = read_inffile(filenmbase)
            DMstr = "%.2f"%info.dm
            DMs.append(info.dm)
            N, dt = int(info.N), info.dt
            obstime = N * dt
            if (filenm == args[0]):
                orig_N = N
                orig_dt = dt
            infile = open(filenmbase+'.dat', mode='rb')
            outfile = open(filenmbase+'.singlepulse', mode='w')
            print 'Single-pulse searching   "%s"...'%filenm
            fileptr = 0
            dataptr = 0

            # Step through the file
            chunk = zeros(fileblocklen, typecode='d')
            dm_candlist = []
            while (fileptr < N):
                if (N-fileptr < fileblocklen):
                    tmpchunk = scipy.io.fread(infile, N-fileptr, 'f')
                else:
                    tmpchunk = scipy.io.fread(infile, fileblocklen, 'f')
                fileptr += len(tmpchunk)

                # Detrend the data in a piece-wise linear fashion
                bp = compress(break_points<len(tmpchunk), break_points)
                tmpchunk = scipy.signal.detrend(tmpchunk, bp=bp)
                # Compute the standard deviation of the data
                if (len(tmpchunk) > 0.1*fileblocklen):
                    stdev = scipy.stats.std(tmpchunk)
                # Normalize so that the data has mean=0, std=1.0
                tmpchunk = tmpchunk/stdev

                # Take care of beginning and end of file overlap issues
                if fileptr==0: # Beginning of file
                    chunk[overlap:] += tmpchunk[:-overlap]
                    infile.seek(-overlap*float_len, 1)
                    fileptr -= overlap
                elif (len(tmpchunk) < fileblocklen):  # End of file
                    chunk[:len(tmpchunk)] += tmpchunk
                    chunk[len(tmpchunk):] = scipy.stats.mean(tmpchunk)
                else:
                    chunk = tmpchunk

                # This is the good part of the data (end effects removed)
                goodchunk = chunk[overlap:-overlap]

                # Search non-downsampled data first
                hibins = compress(goodchunk>threshold, arange(chunklen))
                hivals = take(goodchunk, hibins)
                hibins += dataptr
                # Add the candidates (which are sorted by bin)
                for bin, val in zip(hibins, hivals):
                    time = bin * dt
                    dm_candlist.append(candidate(info.dm, val, time, bin, 1))

                # Now do the downsampling...
                for ii, downfact in enumerate(downfacts):
                    # The normalization of this kernel keeps the
                    # post-smoothing RMS = 1
                    kernel = ones(downfact, typecode='d') / sqrt(downfact)
                    smoothed_chunk = scipy.signal.convolve(chunk, kernel, 1)
                    goodchunk = smoothed_chunk[overlap:-overlap]
                    hibins = compress(goodchunk>threshold, arange(chunklen))
                    hivals = take(goodchunk, hibins)
                    hibins += dataptr
                    hibins = hibins.tolist()
                    hivals = hivals.tolist()
                    # Now walk through the new candidates and remove those
                    # that are not the highest but are within downfact/2
                    # bins of a higher signal pulse
                    toremove = []
                    for ii in range(0, len(hibins)-1):
                        if ii in toremove:  continue
                        xbin, xsigma = hibins[ii], hivals[ii]
                        for jj in range(ii+1, len(hibins)):
                            ybin, ysigma = hibins[jj], hivals[jj]
                            if (abs(ybin-xbin) > downfact/2):
                                break
                            else:
                                if jj in toremove:
                                    continue
                                if (xsigma > ysigma):
                                    if jj not in toremove:
                                        toremove.append(jj)
                                else:
                                    if ii not in toremove:
                                        toremove.append(ii)
                    # Now zap them starting from the end
                    toremove.sort()
                    toremove.reverse()
                    for bin in toremove:
                        del(hibins[bin])
                        del(hivals[bin])

                    # Insert the new candidates into the candlist, but
                    # keep it sorted...
                    for bin, val in zip(hibins, hivals):
                        time = bin * dt
                        bisect.insort(dm_candlist,
                                      candidate(info.dm, val, time, bin, downfact))

                # Backup in the infile by the overlap
                infile.seek(-overlap*float_len, 1)
                if (fileptr < N):
                    fileptr -= overlap
                # Update the data pointer
                dataptr += chunklen

            # Now walk through the dm_candlist and remove the ones that
            # are within the downsample proximity of a higher
            # signal-to-noise pulse
            toremove = []
            for ii in range(0, len(dm_candlist)-1):
                if ii in toremove:  continue
                xx = dm_candlist[ii]
                xbin, xsigma = xx.bin, xx.sigma
                for jj in range(ii+1, len(dm_candlist)):
                    yy = dm_candlist[jj]
                    ybin, ysigma = yy.bin, yy.sigma
                    if (abs(ybin-xbin) > max(downfacts)/2):
                        break
                    else:
                        if jj in toremove:
                            continue
                        prox = max([xx.downfact/2, yy.downfact/2, 1])
                        if (abs(ybin-xbin) <= prox):
                            if (xsigma > ysigma):
                                if jj not in toremove:
                                    toremove.append(jj)
                            else:
                                if ii not in toremove:
                                    toremove.append(ii)
            # Now zap them starting from the end
            toremove.sort()
            toremove.reverse()
            for bin in toremove:
                del(dm_candlist[bin])
            print "  Found %d pulse candidates"%len(dm_candlist)

            # Write the pulses to an ASCII output file
            if len(dm_candlist):
                #dm_candlist.sort(cmp_sigma)
                outfile.write("# DM      Sigma      Time (s)     Sample    Downfact\n")
                for cand in dm_candlist:
                    outfile.write(str(cand))
            outfile.close()

            # Add these candidates to the overall candidate list
            for cand in dm_candlist:
                candlist.append(cand)
            num_v_DMstr[DMstr] = len(dm_candlist)

    # Step through the candidates to make a SNR list
    maxsnr = 0.0
    snrs = []
    for cand in candlist:
        snrs.append(cand.sigma)
        if cand.sigma > maxsnr:
            maxsnr = cand.sigma
    maxsnr = int(maxsnr) + 3

    # Generate the SNR histogram
    snrs = asarray(snrs)
    (num_v_snr, lo_snr, d_snr, num_out_of_range) = \
                scipy.stats.histogram(snrs,
                                      int(maxsnr-threshold+1),
                                      [threshold, maxsnr])
    snrs = arange(maxsnr-threshold+1, typecode='d')*d_snr + lo_snr + 0.5*d_snr
    num_v_snr = where(num_v_snr==0, 0.001, num_v_snr)

    # Generate the DM histogram
    num_v_DM = zeros(len(DMs))
    for ii, DM in enumerate(DMs):
        num_v_DM[ii] = num_v_DMstr["%.2f"%DM]
    DMs = asarray(DMs)

    # open the plot device
    if pgplot_device:
        ppgplot.pgopen(pgplot_device)
    else:
        ppgplot.pgopen(filenmbase+'_singlepulse.ps/VPS')
    ppgplot.pgpap(7.5, 1.0)  # Width in inches, aspect

    # plot the SNR histogram
    ppgplot.pgsvp(0.06, 0.31, 0.6, 0.87)
    ppgplot.pgswin(threshold, maxsnr, log10(0.5), log10(2*max(num_v_snr)))
    ppgplot.pgsch(0.8)
    ppgplot.pgbox("BCNST", 0, 0, "BCLNST", 0, 0)
    ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Signal-to-Noise")
    ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
    ppgplot.pgsch(1.0)
    ppgplot.pgbin(snrs, log10(num_v_snr), 1)

    # plot the DM histogram
    ppgplot.pgsvp(0.39, 0.64, 0.6, 0.87)
    ppgplot.pgswin(min(DMs)-0.5, max(DMs)+0.5, 0.0, 1.1*max(num_v_DM))
    ppgplot.pgsch(0.8)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
    ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "DM (pc cm\u-3\d)")
    ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
    ppgplot.pgsch(1.0)
    ppgplot.pgbin(DMs, num_v_DM, 1)
    
    # plot the SNR vs DM plot 
    ppgplot.pgsvp(0.72, 0.97, 0.6, 0.87)
    ppgplot.pgswin(min(DMs)-0.5, max(DMs)+0.5, threshold, maxsnr)
    ppgplot.pgsch(0.8)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
    ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "DM (pc cm\u-3\d)")
    ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Signal-to-Noise")
    ppgplot.pgsch(1.0)
    cand_ts = zeros(len(candlist), typecode='f')
    cand_SNRs = zeros(len(candlist), typecode='f')
    cand_DMs = zeros(len(candlist), typecode='f')
    for ii, cand in enumerate(candlist):
        cand_ts[ii], cand_SNRs[ii], cand_DMs[ii] = \
                     cand.time, cand.sigma, cand.DM
    ppgplot.pgpt(cand_DMs, cand_SNRs, 20)

    # plot the DM vs Time plot
    ppgplot.pgsvp(0.06, 0.97, 0.08, 0.52)
    ppgplot.pgswin(0.0, obstime, min(DMs)-0.5, max(DMs)+0.5)
    ppgplot.pgsch(0.8)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
    ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Time (s)")
    ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "DM (pc cm\u-3\d)")
    # Circles are symbols 20-26 in increasing order
    snr_range = 15.0-threshold
    cand_symbols = (cand_SNRs-threshold)/snr_range * 6.0 + 0.5 + 20.0
    cand_symbols = cand_symbols.astype('i')
    cand_symbols = where(cand_symbols>26, 26, cand_symbols)
    for ii in range(len(cand_symbols)):
        ppgplot.pgpt(cand_ts[ii:ii+1], cand_DMs[ii:ii+1], cand_symbols[ii:ii+1])

    # Now fill the infomation area
    ppgplot.pgsvp(0.05, 0.95, 0.87, 0.97)
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt('T', 0.5, 0.0, 0.0,
                   "Single pulse results for '%s'"%\
                   filenmbase[:filenmbase.find("_DM")])
    ppgplot.pgsch(0.8)
    # first row
    ppgplot.pgmtxt('T', -1.1, 0.02, 0.0, 'Source: %s'%\
                   info.object)
    ppgplot.pgmtxt('T', -1.1, 0.33, 0.0, 'RA (J2000):')
    ppgplot.pgmtxt('T', -1.1, 0.5, 0.0,
                   coord_to_string(info.ra_h, info.ra_m, info.ra_s))
    ppgplot.pgmtxt('T', -1.1, 0.73, 0.0, 'N samples: %.0f'%orig_N)
    # second row
    ppgplot.pgmtxt('T', -2.4, 0.02, 0.0, 'Telescope: %s'%\
                   info.telescope)
    ppgplot.pgmtxt('T', -2.4, 0.33, 0.0, 'DEC (J2000):')
    ppgplot.pgmtxt('T', -2.4, 0.5, 0.0,
                   coord_to_string(info.dec_d, info.dec_m, info.dec_s))
    ppgplot.pgmtxt('T', -2.4, 0.73, 0.0, 'Sampling time: %.2f \gms'%\
                   (orig_dt*1e6))
    # third row
    if info.instrument.find("pigot"):
        instrument = "Spigot"
    else:
        instrument = info.instrument
    ppgplot.pgmtxt('T', -3.7, 0.02, 0.0, 'Instrument: %s'%instrument)
    if (info.bary):
        ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, 'MJD\dbary\u: %.12f'%\
                       (info.mjd_i+info.mjd_f))
    else:
        ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, 'MJD\dtopo\u: %.12f'%\
                       (info.mjd_i+info.mjd_f))
    ppgplot.pgmtxt('T', -3.7, 0.73, 0.0, 'Freq\dctr\u: %.1f MHz'%\
                   ((info.num_chan/2-0.5)*info.chan_wid+info.freq))
    ppgplot.pgiden()
    ppgplot.pgend()
