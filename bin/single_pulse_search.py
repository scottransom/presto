#!/usr/bin/env python
from __future__ import print_function
from builtins import str, zip
from optparse import OptionParser
import bisect, os, sys, glob
import scipy, scipy.signal, scipy.stats
from presto.presto import rfft, next2_to_n
from presto import infodata
from presto.Pgplot import *

# This is for Python 2/3 comptibility
def mycmp(a, b):
    return ((a > b) - (a < b))

class candidate(object):
    def __init__(self, DM, sigma, time, bin, downfact):
        self.DM = DM
        self.sigma = sigma
        self.time = time
        self.bin = bin
        self.downfact = downfact
    def __str__(self):
        return "%7.2f %7.2f %13.6f %10d     %3d\n"%\
               (self.DM, self.sigma, self.time, self.bin, self.downfact)
    def __eq__(self, other):
        return (self.bin == other.bin)
    def __ne__(self, other):
        return (self.bin != other.bin)
    def __lt__(self, other):
        return (self.bin < other.bin)
    def __le__(self, other):
        return (self.bin <= other.bin)
    def __gt__(self, other):
        return (self.bin > other.bin)
    def __ge__(self, other):
        return (self.bin >= other.bin)
    def __cmp__(self, other):
	# Sort by time (i.e. bin) by default)
        return mycmp(self.bin, other.bin)

def cmp_sigma(self, other):
    #Comparison function to sort candidates by significance
    retval = -mycmp(self.sigma, other.sigma)
    return retval

def fft_convolve(fftd_data, fftd_kern, lo, hi):
    """
    fft_convolve(fftd_data, fftd_kern, lo, hi):
        Perform a convolution with the complex floating point vectors
            'fftd_data' and 'fftd_kern'.  The returned vector will start at
            at bin 'lo' (must be an integer), and go up to but not
            include bin 'hi' (also an integer).
    """
    # Note:  The initial FFTs should be done like:
    # fftd_kern = rfft(kernel, -1)
    # fftd_data = rfft(data, -1)
    prod = Num.multiply(fftd_data, fftd_kern)
    prod.real[0] = fftd_kern.real[0] * fftd_data.real[0]
    prod.imag[0] = fftd_kern.imag[0] * fftd_data.imag[0]
    return rfft(prod, 1)[lo:hi].astype(Num.float32)

def make_fftd_kerns(downfacts, fftlen):
    fftd_kerns = []
    for downfact in downfacts:
        kern = Num.zeros(fftlen, dtype=Num.float32)
        # These offsets produce kernels that give results
        # equal to scipy.signal.convolve
        if downfact % 2:  # Odd number
            kern[:downfact//2+1] += 1.0
            kern[-(downfact//2):] += 1.0
        else:             # Even number
            kern[:downfact//2+1] += 1.0
            if (downfact > 2):
                kern[-(downfact//2-1):] += 1.0
        # The following normalization preserves the
        # RMS=1 characteristic of the data
        fftd_kerns.append(rfft(kern / Num.sqrt(downfact), -1))
    return fftd_kerns

def prune_related1(hibins, hivals, downfact):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the raw 
    # candidate arrays and uses the single downfact
    # that they were selected with.
    toremove = set()
    for ii in range(0, len(hibins)-1):
        if ii in toremove:  continue
        xbin, xsigma = hibins[ii], hivals[ii]
        for jj in range(ii+1, len(hibins)):
            ybin, ysigma = hibins[jj], hivals[jj]
            if (abs(ybin-xbin) > downfact//2):
                break
            else:
                if jj in toremove:
                    continue
                if (xsigma > ysigma):
                    toremove.add(jj)
                else:
                    toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for bin in toremove:
        del(hibins[bin])
        del(hivals[bin])
    return hibins, hivals
    
def prune_related2(dm_candlist, downfacts):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the candidate 
    # instances and looks at the different downfacts of the
    # the different candidates.
    toremove = set()
    for ii in range(0, len(dm_candlist)-1):
        if ii in toremove:  continue
        xx = dm_candlist[ii]
        xbin, xsigma = xx.bin, xx.sigma
        for jj in range(ii+1, len(dm_candlist)):
            yy = dm_candlist[jj]
            ybin, ysigma = yy.bin, yy.sigma
            if (abs(ybin-xbin) > max(downfacts)//2):
                break
            else:
                if jj in toremove:
                    continue
                prox = max([xx.downfact//2, yy.downfact//2, 1])
                if (abs(ybin-xbin) <= prox):
                    if (xsigma > ysigma):
                        toremove.add(jj)
                    else:
                        toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for bin in toremove:
        del(dm_candlist[bin])
    return dm_candlist

def prune_border_cases(dm_candlist, offregions):
    # Ignore those that are located within a half-width
    # of the boundary between data and padding
    #print offregions
    toremove = set()
    for ii in range(len(dm_candlist)-1, -1, -1):
        cand = dm_candlist[ii]
        loside = cand.bin-cand.downfact/2
        hiside = cand.bin+cand.downfact/2
        if hiside < offregions[0][0]: break
        for off, on in offregions:
            if (hiside > off and loside < on):
                toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for ii in toremove:
        del(dm_candlist[ii])
    return dm_candlist

full_usage = """
usage:  single_pulse_search.py [options] .dat files _or_ .singlepulse files
  [-h, --help]        : Display this help
  [-m, --maxwidth]    : Set the max downsampling in sec (see below for default)
  [-p, --noplot]      : Look for pulses but do not generate a plot
  [-t, --threshold]   : Set a different threshold SNR (default=5.0)
  [-x, --xwin]        : Don't make a postscript plot, just use an X-window
  [-s, --start]       : Only plot events occuring after this time (s)
  [-e, --end]         : Only plot events occuring before this time (s)
  [-g, --glob]        : Use the files from these glob expressions (in quotes)
  [-f, --fast]        : Use a less-accurate but much faster method of detrending
  [-b, --nobadblocks] : Don't check for bad-blocks (may save strong pulses)
  [-d, --detrendlen]  : Chunksize for detrending (pow-of-2 in 1000s, default=1)

  Perform a single-pulse search (or simply re-plot the results of a
  single-pulse search) on a set of de-dispersed time series (.dat
  files).

  The search attempts to find pulses by matched-filtering the data with
  a series of different width boxcar functions.  The possible boxcar
  sizes are [1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]
  bins.  By default the boxcars <= 30 are used.  You can specify
  that the larger boxcars are used with the -m (or --maxwidth) option.

  The matched filtering (and accounting for all the possible 'phase'
  offsets of each boxcar) is accomplished by convolving the boxcars
  with the full resolution data.  'Duplicate' candidates from this
  process are filtered, leaving only the most significant.  The time
  series are initially smoothed (by default) using a piecewise linear
  fit to the data where each piece is 1000 data points long.

  If the input files are .singlepulse files, we won't actually perform
  a search, we'll only read in the output .singlepulse files and make
  a plot using the information they contain (along with the
  corresponding .inf files).

  Notes on usage and performance:

    -- single_pulse_search.py is tuned for finding *narrow* pulses
       (i.e. those of only a few bins width).  Because of this, you
       should always search appropriately downsampled data (as
       recommended by DDplan.py, for instance) where dispersion
       smearing is <~ 1 time series bin.

    -- the linear-piecewise detrending is very useful in long
       observations with modern instrumentation where you can see
       long timescale power fluctuations.  Strong pulses can skew the
       statistics of the 1000-bin chunks, though, and caused some
       suppression in the detection levels of bright pulses (weak
       pulses are mostly unaffected since they don't strongly change
       the statistics).  If your data have no long-timescale
       fluctuations (for instance, if you are processing old 1-bit
       analog filterbank data which is AC-coupled or if you remove
       rednoise via realfft/rednoise/(inverse-)realfft), I recommend
       using the -f/--fast flag.  And if you want to find wide
       pulses, it might be worth making the chunksize bigger (i.e.
       4000 or 8000).

    -- The bad-block detection and removal code can and does remove
       blocks that have very strong, and particularly strong and broad,
       pulses in them.  It can also quite effectively remove RFI-
       infused portions of the data.  Whether to turn it on or off
       depends on your data.  Note that if there are multiple pulses,
       only the brightest will usually be "bad-blocked" and removed.

    -- The fourier-domain matched filtering used here has no phase-
       dependent effects.  So a 15-bin pulse can be found with equal
       significance no matter which bin it starts in in the time series.

    -- The definition of "sigma" used is possibly slightly different
       from that used in other codes for S/N:
           sigma = sum(signal-bkgd_level)/RMS/sqrt(boxcar_width)
       where the bkgd_level is typically 0 after detrending and RMS=1
       after normalization.  This definition has the advantage that
       you will get (basically) the same sigma for any pulse no
       matter how much the input time series has been downsampled as
       long as the pulse is still resolved.

  Copyright Scott Ransom <sransom@nrao.edu>, 2015
"""
usage = "usage: %prog [options] .dat files _or_ .singlepulse files"
    
def read_singlepulse_files(infiles, threshold, T_start, T_end):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    for ii, infile in enumerate(infiles):
        if infile.endswith(".singlepulse"):
            filenmbase = infile[:infile.rfind(".singlepulse")]
        else:
            filenmbase = infile
        info = infodata.infodata(filenmbase+".inf")
        DMstr = "%.2f"%info.DM
        DMs.append(info.DM)
        num_v_DMstr[DMstr] = 0
        if ii==0:
            info0 = info
        if os.stat(infile)[6]:
            try:
                cands = Num.loadtxt(infile)
                if len(cands.shape)==1:
                    cands = Num.asarray([cands])
                for cand in cands:
                    if cand[2] < T_start: continue
                    if cand[2] > T_end: break
                    if cand[1] >= threshold:
                        candlist.append(candidate(*cand))
                        num_v_DMstr[DMstr] += 1
            except:  # No candidates in the file
                IndexError
    DMs.sort()
    return info0, DMs, candlist, num_v_DMstr

def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Don't make a postscript plot, just use an X-window")
    parser.add_option("-p", "--noplot", action="store_false", dest="makeplot",
                      default=True, help="Look for pulses but do not generate a plot")
    parser.add_option("-m", "--maxwidth", type="float", dest="maxwidth", default=0.0,
                      help="Set the max downsampling in sec (see below for default)")
    parser.add_option("-t", "--threshold", type="float", dest="threshold", default=5.0,
                      help="Set a different threshold SNR (default=5.0)")
    parser.add_option("-s", "--start", type="float", dest="T_start", default=0.0,
                      help="Only plot events occuring after this time (s)")
    parser.add_option("-e", "--end", type="float", dest="T_end", default=1e9,
                      help="Only plot events occuring before this time (s)")
    parser.add_option("-g", "--glob", type="string", dest="globexp", default=None,
                      help="Process the files from this glob expression")
    parser.add_option("-f", "--fast", action="store_true", dest="fast",
                      default=False, help="Use a faster method of de-trending (2x speedup)")
    parser.add_option("-b", "--nobadblocks", action="store_false", dest="badblocks",
                      default=True, help="Don't check for bad-blocks (may save strong pulses)")
    parser.add_option("-d", "--detrendlen", type="int", dest="detrendfact", default=1,
                      help="Chunksize for detrending (pow-of-2 in 1000s)")
    (opts, args) = parser.parse_args()
    if len(args)==0:
        if opts.globexp==None:
            print(full_usage)
            sys.exit(0)
        else:
            args = []
            for globexp in opts.globexp.split():
                args += glob.glob(globexp)
    useffts = True
    dosearch = True
    if opts.xwin:
        pgplot_device = "/XWIN"
    else:
        pgplot_device = ""

    fftlen = 8192     # Should be a power-of-two for best speed
    chunklen = 8000   # Must be at least max_downfact less than fftlen
    assert(opts.detrendfact in [1,2,4,8,16,32])
    detrendlen = opts.detrendfact*1000
    if (detrendlen > chunklen):
        chunklen = detrendlen
        fftlen = int(next2_to_n(chunklen))
    blocks_per_chunk = chunklen // detrendlen
    overlap = (fftlen - chunklen) // 2
    worklen = chunklen + 2*overlap  # currently it is fftlen...

    max_downfact = 30
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]

    if args[0].endswith(".singlepulse"):
        filenmbase = args[0][:args[0].rfind(".singlepulse")]
        dosearch = False
    elif args[0].endswith(".dat"):
        filenmbase = args[0][:args[0].rfind(".dat")]
    else:
        filenmbase = args[0]

    # Don't do a search, just read results and plot
    if not dosearch:
        info, DMs, candlist, num_v_DMstr = \
              read_singlepulse_files(args, opts.threshold, opts.T_start, opts.T_end)
        orig_N, orig_dt = int(info.N), info.dt
        obstime = orig_N * orig_dt
    else:
        DMs = []
        candlist = []
        num_v_DMstr = {}

        # Loop over the input files
        for filenm in args:
            if filenm.endswith(".dat"):
                filenmbase = filenm[:filenm.rfind(".dat")]
            else:
                filenmbase = filenm
            info = infodata.infodata(filenmbase+".inf")
            DMstr = "%.2f"%info.DM
            DMs.append(info.DM)
            N, dt = int(info.N), info.dt
            obstime = N * dt
            # Choose the maximum width to search based on time instead
            # of bins.  This helps prevent increased S/N when the downsampling
            # changes as the DM gets larger.
            if opts.maxwidth > 0.0:
                downfacts = [x for x in default_downfacts if x*dt <= opts.maxwidth]
            else:
                downfacts = [x for x in default_downfacts if x <= max_downfact]
            if len(downfacts) == 0:
                downfacts = [default_downfacts[0]]
            if (filenm == args[0]):
                orig_N = N
                orig_dt = dt
                if useffts:
                    fftd_kerns = make_fftd_kerns(default_downfacts, fftlen)
            if info.breaks:
                offregions = list(zip([x[1] for x in info.onoff[:-1]],
                                 [x[0] for x in info.onoff[1:]]))

                # If last break spans to end of file, don't read it in (its just padding)
                if offregions[-1][1] == N - 1:
                    N = offregions[-1][0] + 1

            outfile = open(filenmbase+'.singlepulse', mode='w')

            # Compute the file length in detrendlens
            roundN = N // detrendlen * detrendlen
            numchunks = roundN // chunklen
            # Read in the file
            print('Reading "%s"...'%filenm)
            timeseries = Num.fromfile(filenm, dtype=Num.float32, count=roundN)
            # Split the timeseries into chunks for detrending
            numblocks = roundN // detrendlen
            timeseries.shape = (numblocks, detrendlen)
            stds = Num.zeros(numblocks, dtype=Num.float64)
            # de-trend the data one chunk at a time
            print('  De-trending the data and computing statistics...')
            for ii, chunk in enumerate(timeseries):
                if opts.fast:  # use median removal instead of detrending (2x speedup)
                    tmpchunk = chunk.copy()
                    tmpchunk.sort()
                    med = tmpchunk[detrendlen//2]
                    chunk -= med
                    tmpchunk -= med
                else:
                    # The detrend calls are the most expensive in the program
                    timeseries[ii] = scipy.signal.detrend(chunk, type='linear')
                    tmpchunk = timeseries[ii].copy()
                    tmpchunk.sort()
                # The following gets rid of (hopefully) most of the 
                # outlying values (i.e. power dropouts and single pulses)
                # If you throw out 5% (2.5% at bottom and 2.5% at top)
                # of random gaussian deviates, the measured stdev is ~0.871
                # of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
                # The following is roughly .std() since we already removed the median
                stds[ii] = Num.sqrt((tmpchunk[detrendlen//40:-detrendlen//40]**2.0).sum() /
                                    (0.95*detrendlen))
            stds *= 1.148
            # sort the standard deviations and separate those with
            # very low or very high values
            sort_stds = stds.copy()
            sort_stds.sort()
            # identify the differences with the larges values (this
            # will split off the chunks with very low and very high stds
            locut = (sort_stds[1:numblocks//2+1] -
                     sort_stds[:numblocks//2]).argmax() + 1
            hicut = (sort_stds[numblocks//2+1:] -
                     sort_stds[numblocks//2:-1]).argmax() + numblocks//2 - 2
            std_stds = scipy.std(sort_stds[locut:hicut])
            median_stds = sort_stds[(locut+hicut)//2]
            print("    pseudo-median block standard deviation = %.2f" % (median_stds))
            if (opts.badblocks):
                lo_std = median_stds - 4.0 * std_stds
                hi_std = median_stds + 4.0 * std_stds
                # Determine a list of "bad" chunks.  We will not search these.
                bad_blocks = Num.nonzero((stds < lo_std) | (stds > hi_std))[0]
                print("    identified %d bad blocks out of %d (i.e. %.2f%%)" % \
                      (len(bad_blocks), len(stds),
                       100.0*float(len(bad_blocks))/float(len(stds))))
                stds[bad_blocks] = median_stds
            else:
                bad_blocks = []
            print("  Now searching...")

            # Now normalize all of the data and reshape it to 1-D
            timeseries /= stds[:,Num.newaxis]
            timeseries.shape = (roundN,)
            # And set the data in the bad blocks to zeros
            # Even though we don't search these parts, it is important
            # because of the overlaps for the convolutions
            for bad_block in bad_blocks:
                loind, hiind = bad_block*detrendlen, (bad_block+1)*detrendlen
                timeseries[loind:hiind] = 0.0
            # Convert to a set for faster lookups below
            bad_blocks = set(bad_blocks)

            # Step through the data
            dm_candlist = []
            for chunknum in range(numchunks):
                loind = chunknum*chunklen-overlap
                hiind = (chunknum+1)*chunklen+overlap
                # Take care of beginning and end of file overlap issues
                if (chunknum==0): # Beginning of file
                    chunk = Num.zeros(worklen, dtype=Num.float32)
                    chunk[overlap:] = timeseries[loind+overlap:hiind]
                elif (chunknum==numchunks-1): # end of the timeseries
                    chunk = Num.zeros(worklen, dtype=Num.float32)
                    chunk[:-overlap] = timeseries[loind:hiind-overlap]
                else:
                    chunk = timeseries[loind:hiind]

                # Make a set with the current block numbers
                lowblock = blocks_per_chunk * chunknum
                currentblocks = set(Num.arange(blocks_per_chunk) + lowblock)
                localgoodblocks = Num.asarray(list(currentblocks -
                                                   bad_blocks)) - lowblock
                # Search this chunk if it is not all bad
                if len(localgoodblocks):
                    # This is the good part of the data (end effects removed)
                    goodchunk = chunk[overlap:-overlap]

                    # need to pass blocks/chunklen, localgoodblocks
                    # dm_candlist, dt, opts.threshold to cython routine

                    # Search non-downsampled data first
                    # NOTE:  these nonzero() calls are some of the most
                    #        expensive calls in the program.  Best bet would 
                    #        probably be to simply iterate over the goodchunk
                    #        in C and append to the candlist there.
                    hibins = Num.flatnonzero(goodchunk>opts.threshold)
                    hivals = goodchunk[hibins]
                    hibins += chunknum * chunklen
                    hiblocks = hibins // detrendlen
                    # Add the candidates (which are sorted by bin)
                    for bin, val, block in zip(hibins, hivals, hiblocks):
                        if block not in bad_blocks:
                            time = bin * dt
                            dm_candlist.append(candidate(info.DM, val, time, bin, 1))

                    # Prepare our data for the convolution
                    if useffts: fftd_chunk = rfft(chunk, -1)

                    # Now do the downsampling...
                    for ii, downfact in enumerate(downfacts):
                        if useffts: 
                            # Note:  FFT convolution is faster for _all_ downfacts, even 2
                            goodchunk = fft_convolve(fftd_chunk, fftd_kerns[ii],
                                                     overlap, -overlap)
                        else:
                            # The normalization of this kernel keeps the post-smoothing RMS = 1
                            kernel = Num.ones(downfact, dtype=Num.float32) / \
                                     Num.sqrt(downfact)
                            smoothed_chunk = scipy.signal.convolve(chunk, kernel, 1)
                            goodchunk = smoothed_chunk[overlap:-overlap]
                        #hibins = Num.nonzero(goodchunk>opts.threshold)[0]
                        hibins = Num.flatnonzero(goodchunk>opts.threshold)
                        hivals = goodchunk[hibins]
                        hibins += chunknum * chunklen
                        hiblocks = hibins // detrendlen
                        hibins = hibins.tolist()
                        hivals = hivals.tolist()
                        # Now walk through the new candidates and remove those
                        # that are not the highest but are within downfact/2
                        # bins of a higher signal pulse
                        hibins, hivals = prune_related1(hibins, hivals, downfact)
                        # Insert the new candidates into the candlist, but
                        # keep it sorted...
                        for bin, val, block in zip(hibins, hivals, hiblocks):
                            if block not in bad_blocks:
                                time = bin * dt
                                bisect.insort(dm_candlist,
                                              candidate(info.DM, val, time, bin, downfact))

            # Now walk through the dm_candlist and remove the ones that
            # are within the downsample proximity of a higher
            # signal-to-noise pulse
            dm_candlist = prune_related2(dm_candlist, downfacts)
            print("  Found %d pulse candidates"%len(dm_candlist))
            
            # Get rid of those near padding regions
            if info.breaks: prune_border_cases(dm_candlist, offregions)

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

    if (opts.makeplot):

        # Step through the candidates to make a SNR list
        DMs.sort()
        snrs = []
        for cand in candlist:
            if not Num.isinf(cand.sigma):
                snrs.append(cand.sigma)
        if snrs:
            maxsnr = max(int(max(snrs)), int(opts.threshold)) + 3
        else:
            maxsnr = int(opts.threshold) + 3

        # Generate the SNR histogram
        snrs = Num.asarray(snrs)
        (num_v_snr, edges) = Num.histogram(snrs,
                                           int(maxsnr-opts.threshold+1),
                                           [opts.threshold, maxsnr])
        snrs = edges[:-1] + 0.5*(edges[1]-edges[0])
        num_v_snr = num_v_snr.astype(Num.float32)
        num_v_snr[num_v_snr==0.0] = 0.001

        # Generate the DM histogram
        num_v_DM = Num.zeros(len(DMs))
        for ii, DM in enumerate(DMs):
            num_v_DM[ii] = num_v_DMstr["%.2f"%DM]
        DMs = Num.asarray(DMs)

        # open the plot device
        short_filenmbase = filenmbase[:filenmbase.find("_DM")]
        if opts.T_end > obstime:
            opts.T_end = obstime
        if pgplot_device:
            ppgplot.pgopen(pgplot_device)
        else:
            if (opts.T_start > 0.0 or opts.T_end < obstime):
                ppgplot.pgopen(short_filenmbase+'_%.0f-%.0fs_singlepulse.ps/VPS'%
                               (opts.T_start, opts.T_end))
            else:
                ppgplot.pgopen(short_filenmbase+'_singlepulse.ps/VPS')
        ppgplot.pgpap(7.5, 1.0)  # Width in inches, aspect

        # plot the SNR histogram
        ppgplot.pgsvp(0.06, 0.31, 0.6, 0.87)
        ppgplot.pgswin(opts.threshold, maxsnr,
                       Num.log10(0.5), Num.log10(2*max(num_v_snr)))
        ppgplot.pgsch(0.8)
        ppgplot.pgbox("BCNST", 0, 0, "BCLNST", 0, 0)
        ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Signal-to-Noise")
        ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
        ppgplot.pgsch(1.0)
        ppgplot.pgbin(snrs, Num.log10(num_v_snr), 1)

        # plot the DM histogram
        ppgplot.pgsvp(0.39, 0.64, 0.6, 0.87)
        # Add [1] to num_v_DM in YMAX below so that YMIN != YMAX when max(num_v_DM)==0
        ppgplot.pgswin(min(DMs)-0.5, max(DMs)+0.5, 0.0, 1.1*max(num_v_DM+[1]))
        ppgplot.pgsch(0.8)
        ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, r"DM (pc cm\u-3\d)")
        ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
        ppgplot.pgsch(1.0)
        ppgplot.pgbin(DMs, num_v_DM, 1)

        # plot the SNR vs DM plot 
        ppgplot.pgsvp(0.72, 0.97, 0.6, 0.87)
        ppgplot.pgswin(min(DMs)-0.5, max(DMs)+0.5, opts.threshold, maxsnr)
        ppgplot.pgsch(0.8)
        ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, r"DM (pc cm\u-3\d)")
        ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "Signal-to-Noise")
        ppgplot.pgsch(1.0)
        cand_ts = Num.zeros(len(candlist), dtype=Num.float32)
        cand_SNRs = Num.zeros(len(candlist), dtype=Num.float32)
        cand_DMs = Num.zeros(len(candlist), dtype=Num.float32)
        for ii, cand in enumerate(candlist):
            cand_ts[ii], cand_SNRs[ii], cand_DMs[ii] = \
                         cand.time, cand.sigma, cand.DM
        ppgplot.pgpt(cand_DMs, cand_SNRs, 20)

        # plot the DM vs Time plot
        ppgplot.pgsvp(0.06, 0.97, 0.08, 0.52)
        ppgplot.pgswin(opts.T_start, opts.T_end, min(DMs)-0.5, max(DMs)+0.5)
        ppgplot.pgsch(0.8)
        ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Time (s)")
        ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, r"DM (pc cm\u-3\d)")
        # Circles are symbols 20-26 in increasing order
        snr_range = 12.0
        cand_symbols = (cand_SNRs-opts.threshold)/snr_range * 6.0 + 20.5
        cand_symbols = cand_symbols.astype(Num.int32)
        cand_symbols[cand_symbols>26] = 26
        for ii in [26, 25, 24, 23, 22, 21, 20]:
            inds = Num.nonzero(cand_symbols==ii)[0]
            ppgplot.pgpt(cand_ts[inds], cand_DMs[inds], ii)

        # Now fill the infomation area
        ppgplot.pgsvp(0.05, 0.95, 0.87, 0.97)
        ppgplot.pgsch(1.0)
        ppgplot.pgmtxt('T', 0.5, 0.0, 0.0,
                       "Single pulse results for '%s'"%short_filenmbase)
        ppgplot.pgsch(0.8)
        # first row
        ppgplot.pgmtxt('T', -1.1, 0.02, 0.0, 'Source: %s'%\
                       info.object)
        ppgplot.pgmtxt('T', -1.1, 0.33, 0.0, 'RA (J2000):')
        ppgplot.pgmtxt('T', -1.1, 0.5, 0.0, info.RA)
        ppgplot.pgmtxt('T', -1.1, 0.73, 0.0, 'N samples: %.0f'%orig_N)
        # second row
        ppgplot.pgmtxt('T', -2.4, 0.02, 0.0, 'Telescope: %s'%\
                       info.telescope)
        ppgplot.pgmtxt('T', -2.4, 0.33, 0.0, 'DEC (J2000):')
        ppgplot.pgmtxt('T', -2.4, 0.5, 0.0, info.DEC)
        ppgplot.pgmtxt('T', -2.4, 0.73, 0.0, 'Sampling time: %.2f \gms'%\
                       (orig_dt*1e6))
        # third row
        if info.instrument.find("pigot") >= 0:
            instrument = "Spigot"
        else:
            instrument = info.instrument
        ppgplot.pgmtxt('T', -3.7, 0.02, 0.0, 'Instrument: %s'%instrument)
        if (info.bary):
            ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, r'MJD\dbary\u: %.12f'%info.epoch)
        else:
            ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, r'MJD\dtopo\u: %.12f'%info.epoch)
        ppgplot.pgmtxt('T', -3.7, 0.73, 0.0, r'Freq\dctr\u: %.1f MHz'%\
                       ((info.numchan/2-0.5)*info.chan_width+info.lofreq))
        ppgplot.pgiden()
        ppgplot.pgend()

if __name__ == '__main__':
    if (0):
        # The following is for profiling
        import hotshot
        prof = hotshot.Profile("hotshot_edi_stats")
        prof.runcall(main)
        prof.close()
        # To see the results:
        if (0):
            from presto.hotshot import stats
            s = stats.load("hotshot_edi_stats")
            s.sort_stats("time").print_stats()
    else:
        # print(sys.version)
        main()
