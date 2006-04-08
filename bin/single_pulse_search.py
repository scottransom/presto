#!/usr/bin/env python
import bisect, os, sys, getopt, infodata
import scipy, scipy.signal, umath, ppgplot
from Numeric import *
from presto import rfft
from psr_utils import coord_to_string
from TableIO import readColumns
from optparse import OptionParser

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
	# Sort by time (i.e. bin) by default)
        return cmp(self.bin, other.bin)

def cmp_sigma(self, other):
    #Comparison function to sort candidates by significance
    retval = -cmp(self.sigma, other.sigma)
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
    prod = fftd_data * fftd_kern
    prod.real[0] = fftd_kern.real[0] * fftd_data.real[0]
    prod.imag[0] = fftd_kern.imag[0] * fftd_data.imag[0]
    return rfft(prod, 1)[lo:hi].astype('d')

def make_fftd_kerns(downfacts, fftlen):
    fftd_kerns = []
    for downfact in downfacts:
        kern = zeros(fftlen, typecode='d')
        # These offsets produce kernels that give results
        # equal to scipy.signal.convolve
        if downfact % 2:  # Odd number
            kern[:downfact/2+1] += 1.0
            kern[-(downfact/2):] += 1.0
        else:             # Even number
            kern[:downfact/2+1] += 1.0
            if (downfact > 2):
                kern[-(downfact/2-1):] += 1.0
        # The following normalization preserves the
        # RMS=1 characteristic of the data
        fftd_kerns.append(rfft(kern/sqrt(downfact), -1))
    return fftd_kerns

def prune_related1(hibins, hivals, downfact):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the raw 
    # candidate arrays and uses the single downfact
    # that they were selected with.
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
    return hibins, hivals
    
def prune_related2(dm_candlist, downfacts):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the candidate 
    # instances and looks at the different downfacts of the
    # the different candidates.
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
    return dm_candlist

def prune_border_cases(dm_candlist, offregions):
    # Ignore those that are locate in a half-width
    # of the boundary between data and padding
    #print offregions
    toremove = []
    for ii in range(len(dm_candlist))[::-1]:
        cand = dm_candlist[ii]
        loside = cand.bin-cand.downfact/2
        hiside = cand.bin+cand.downfact/2
        if hiside < offregions[0][0]: break
        for off, on in offregions:
            if (hiside > off and loside < on):
                toremove.append(ii)
    # Now zap them starting from the end
    toremove.sort()
    toremove.reverse()
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

  Perform a single-pulse search (or simply re-plot the results of a
  single-pulse search) on a set of de-dispersed time series (.dat
  files).

  The search attempts to find pulses by matched-filtering the data
  with a series of different width boxcar functions.  The possible
  boxcar sizes are [1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150]
  bins.  By default the boxcars <= 30 are used.  You can specify
  that the larger boxcars are used with the -m (or --maxwidth) option.

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
usage = "usage: %prog [options] .dat files _or_ .singlepulse files"
    
def read_singlepulse_files(infiles, threshold, T_start, T_end):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    for ii, infile in enumerate(infiles):
        filenmbase = infile.rstrip(".singlepulse")
        info = infodata.infodata(filenmbase+".inf")
        DMstr = "%.2f"%info.DM
        DMs.append(info.DM)
        num_v_DMstr[DMstr] = 0
        if ii==0:
            info0 = info
        if os.stat(infile)[6]:
            cands = transpose(readColumns(infile, "#"))
            for cand in cands:
                if cand[2] < T_start: continue
                if cand[2] > T_end: break
                if cand[1] >= threshold:
                    candlist.append(candidate(*cand))
                    num_v_DMstr[DMstr] += 1
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
    (opts, args) = parser.parse_args()
    if len(args)==0:
        print full_usage
        sys.exit(0)
    useffts = True
    dosearch = True
    max_downfact = 30
    fftlen = 8192    # Should be a power-of-two for best speed
    chunklen = 8000  # Must be at least max_downfact less than fftlen
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150]
    if opts.xwin:
        pgplot_device = "/XWIN"
    else:
        pgplot_device = ""
    break_points = chunklen/5 * arange(1,5)
    overlap = (fftlen - chunklen)/2
    fileblocklen = chunklen + 2*overlap  # currently it is fftlen...
    float_len = 4
    if args[0].endswith(".singlepulse"):
        filenmbase = args[0].rstrip(".singlepulse")
        dosearch = False
    else:
        filenmbase = args[0].rstrip(".dat")

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
            filenmbase = filenm.rstrip(".dat")
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
                    fftd_kerns = make_fftd_kerns(downfacts, fftlen)
            if info.breaks:
                offregions = zip([x[1] for x in info.onoff[:-1]],
                                 [x[0] for x in info.onoff[1:]])
            infile = open(filenmbase+'.dat', mode='rb')
            outfile = open(filenmbase+'.singlepulse', mode='w')
            print 'Single-pulse searching   "%s"...'%filenm
            fileptr = 0
            dataptr = 0
            oldstdev = 0.0
            stdev = 0.0

            # Step through the file
            dm_candlist = []
            while (fileptr < N):
                search_chunk = True
                if (N-fileptr < fileblocklen):
                    tmpchunk = scipy.io.fread(infile, N-fileptr, 'f')
                else:
                    tmpchunk = scipy.io.fread(infile, fileblocklen, 'f')
                numread = len(tmpchunk)

                # Detrend the data in a piece-wise linear fashion
                tmpchunk = tmpchunk.astype('d') # keep the precision high for detrend
                bp = compress(break_points<numread, break_points)
                tmpchunk = scipy.signal.detrend(tmpchunk, bp=bp)

                # Compute the standard deviation of the data if we are not
                # in a padding region
                padding = 0
                if info.breaks:
                    newfileptr = fileptr+numread
                    for off, on in offregions:
                        if ((off < fileptr < on) or
                            (off < newfileptr < on) or
                            (fileptr < off and newfileptr > on)):
                            padding = 1
                            break

                if (numread > 0.1*fileblocklen and not padding):
                    stdev = scipy.stats.std(tmpchunk)
                    if oldstdev==0.0:  oldstdev = stdev
                    oldstdev = 0.5*oldstdev + 0.5*stdev

                # Make sure that we are not searching a section of padding
                if padding or stdev < 1e-7 or stdev < 0.3*oldstdev:
                    search_chunk = False

                # Normalize so that the data has mean=0, std=1.0
                tmpchunk = tmpchunk/stdev

                # Take care of beginning and end of file overlap issues
                if (fileptr==0 and numread==fileblocklen): # Beginning of file
                    chunk = zeros(fileblocklen, typecode='d')
                    chunk[overlap:] = tmpchunk[:-overlap]
                    infile.seek(-overlap*float_len, 1)
                    fileptr -= overlap
                elif (numread < fileblocklen): # End of file
                    chunk = zeros(fileblocklen, typecode='d')
                    chunk[:numread] = tmpchunk
                else:
                    chunk = tmpchunk
                fileptr += numread

                if search_chunk:
                    # This is the good part of the data (end effects removed)
                    goodchunk = chunk[overlap:-overlap]

                    # Search non-downsampled data first
                    #
                    # NOTE:  these compress() calls (and the nonzeros() that result) are
                    #        currently the most expensive call in the program.  Best
                    #        bet would probably be to simply iterate over the goodchunk
                    #        in C and append to the candlist there.
                    hibins = compress(goodchunk>opts.threshold, arange(chunklen))
                    hivals = take(goodchunk, hibins)
                    hibins += dataptr
                    # Add the candidates (which are sorted by bin)
                    for bin, val in zip(hibins, hivals):
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
                            kernel = ones(downfact, typecode='d') / sqrt(downfact)
                            smoothed_chunk = scipy.signal.convolve(chunk, kernel, 1)
                            goodchunk = smoothed_chunk[overlap:-overlap]
                        hibins = compress(goodchunk>opts.threshold, arange(chunklen))
                        hivals = take(goodchunk, hibins)
                        hibins += dataptr
                        hibins = hibins.tolist()
                        hivals = hivals.tolist()
                        # Now walk through the new candidates and remove those
                        # that are not the highest but are within downfact/2
                        # bins of a higher signal pulse
                        hibins, hivals = prune_related1(hibins, hivals, downfact)
                        # Insert the new candidates into the candlist, but
                        # keep it sorted...
                        for bin, val in zip(hibins, hivals):
                            time = bin * dt
                            bisect.insort(dm_candlist,
                                          candidate(info.DM, val, time, bin, downfact))

                # Backup in the infile by the overlap
                infile.seek(-2*overlap*float_len, 1)
                if (fileptr < N): fileptr -= 2*overlap
                # Update the data pointer
                dataptr += chunklen

            # Now walk through the dm_candlist and remove the ones that
            # are within the downsample proximity of a higher
            # signal-to-noise pulse
            dm_candlist = prune_related2(dm_candlist, downfacts)
            print "  Found %d pulse candidates"%len(dm_candlist)
            
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
                                          int(maxsnr-opts.threshold+1),
                                          [opts.threshold, maxsnr])
        snrs = arange(maxsnr-opts.threshold+1, typecode='d')*d_snr + lo_snr + 0.5*d_snr
        num_v_snr = where(num_v_snr==0, 0.001, num_v_snr)

        # Generate the DM histogram
        num_v_DM = zeros(len(DMs))
        for ii, DM in enumerate(DMs):
            num_v_DM[ii] = num_v_DMstr["%.2f"%DM]
        DMs = asarray(DMs)

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
        ppgplot.pgswin(opts.threshold, maxsnr, log10(0.5), log10(2*max(num_v_snr)))
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
        ppgplot.pgswin(min(DMs)-0.5, max(DMs)+0.5, opts.threshold, maxsnr)
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
        ppgplot.pgswin(opts.T_start, opts.T_end, min(DMs)-0.5, max(DMs)+0.5)
        ppgplot.pgsch(0.8)
        ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Time (s)")
        ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "DM (pc cm\u-3\d)")
        # Circles are symbols 20-26 in increasing order
        snr_range = 12.0
        cand_symbols = (cand_SNRs-opts.threshold)/snr_range * 6.0 + 0.5 + 20.0
        cand_symbols = cand_symbols.astype('i')
        cand_symbols = where(cand_symbols>26, 26, cand_symbols)
        for ii in [26, 25, 24, 23, 22, 21, 20]:
            inds = nonzero(equal(cand_symbols, ii))
            ts = take(cand_ts, inds)
            DMs = take(cand_DMs, inds)
            ppgplot.pgpt(ts, DMs, ii)

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
            ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, 'MJD\dbary\u: %.12f'%info.epoch)
        else:
            ppgplot.pgmtxt('T', -3.7, 0.33, 0.0, 'MJD\dtopo\u: %.12f'%info.epoch)
        ppgplot.pgmtxt('T', -3.7, 0.73, 0.0, 'Freq\dctr\u: %.1f MHz'%\
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
    else:
        main()
