#!/usr/bin/env python
#
# Script to fold the profile from the dat-file and plot it
# or to plot the time series, or plot the stack of individual
# pulses or subintegrations in different modes
# Also added functionality to use the tim-file in Sigproc-style
# format, and any other raw binary (float) data files
# if corresponding info is provided with --mjd (or --chandra),
# --tsamp and --pulsar options.
# Also, possible to give ascii file with events times (in seconds).
# This is useful to make plots for the X-ray data.
#
# List of options with short help:    -h or --help
# Full list of options without help:  -l or --list
#
# Vlad, Nov 26, 2010 (c)
#
# Last major update - Dec 5, 2010
#
######################################################################

from __future__ import print_function
from builtins import str
from builtins import range
import numpy as np
import scipy
import scipy.signal
import array as ar
import os.path, stat, sys
import math
import optparse as opt
import presto.infodata as inf
import presto.polycos as poly
from presto import sigproc

# importing VariableColormap from kapteyn module of it exists
try:
        from kapteyn.mplutil import VariableColormap
        kapteyn_loaded = True
except ImportError:
        kapteyn_loaded = False

detrendlen = 1000 # length of a linear piecewise chunk of data for detrending


# get the period from the polyco-file with the fileid=pid
# using the best record in the polyco-file for the mjd=mjd
def get_period (id, mjd):
    return 1. / id.get_phs_and_freq(float(str(mjd).split(".")[0]), float("0."+str(mjd).split(".")[1]))[1]


# I took this part from the Scott's single_pulse_search.py
def detrending (data, is_fast):
    # Compute the file length in detrendlens
    roundN = np.size(data)/detrendlen * detrendlen
    data = data[:roundN] # here we redefining the arrat and loosing some samples
    # Split the data into chunks for detrending
    numblocks = roundN/detrendlen
    data.shape = (numblocks, detrendlen)
    stds = np.zeros(numblocks, dtype=np.float64)
    # de-trend the data one chunk at a time and compute statistics
    for ii, chunk in enumerate(data):
        if is_fast:  # use median removal instead of detrending (2x speedup)
            tmpchunk = chunk.copy()
            tmpchunk.sort()
            med = tmpchunk[detrendlen/2]
            chunk -= med
            tmpchunk -= med
        else:
            # The detrend calls are the most expensive in the program
            data[ii] = scipy.signal.detrend(chunk, type='linear')
            tmpchunk = data[ii].copy()
            tmpchunk.sort()
        # The following gets rid of (hopefully) most of the
        # outlying values (i.e. power dropouts and single pulses)
        # If you throw out 5% (2.5% at bottom and 2.5% at top)
        # of random gaussian deviates, the measured stdev is ~0.871
        # of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
        # The following is roughly .std() since we already removed the median
            stds[ii] = np.sqrt((tmpchunk[detrendlen/40:-detrendlen/40]**2.0).sum() / (0.95*detrendlen))

    stds *= 1.148
    # sort the standard deviations and separate those with
    # very low or very high values
    sort_stds = stds.copy()
    sort_stds.sort()
    # identify the differences with the larges values (this
    # will split off the chunks with very low and very high stds
    locut = (sort_stds[1:numblocks/2+1] - sort_stds[:numblocks/2]).argmax() + 1
    hicut = (sort_stds[numblocks/2+1:] - sort_stds[numblocks/2:-1]).argmax() + numblocks/2 - 2
    std_stds = scipy.std(sort_stds[locut:hicut])
    median_stds = sort_stds[(locut+hicut)/2]
    lo_std = median_stds - 4.0 * std_stds
    hi_std = median_stds + 4.0 * std_stds
    # Determine a list of "bad" chunks.  We will not search these.
    bad_blocks = np.nonzero((stds < lo_std) | (stds > hi_std))[0]
    stds[bad_blocks] = median_stds

    # Now normalize all of the data and reshape it to 1-D
    data /= stds[:,np.newaxis]
    data.shape = (roundN,)

    return data

# gives the short list of options without explanations
def list_options(prg):
    print("Usage: %s [options] <.dat OR .tim (use --tim option)>" % (prg))
    print()
    print("Options:")
    print("  [-h, --help]              [-n, --nbins #BINS]     [-p, --pulsar NAME]")
    print("  [--polyco FILE]           [--period PERIOD]       [-b, --block SAMPLES]")
    print("  [-s, --startphase PHASE]  [-e, --endphase PHASE]  [--start TIME]")
    print("  [-w, --window TIME]       [-f, --fast-detrend]    [--no-detrend]")
    print("  [-t, --timeseries]        [-a, --rebin FACTOR]    [-y, --profileonly]")
    print("  [-k, --stacking]          [--offset OFFSET]       [-d, --dump #PULSES | TIME]")
    print("  [--saveprof FILE]         [--saveonly]            [-i, --image FILEEXT]")
    print("  [--tim]                   [--events]              [-l, --list]")
    print("  [-2, --double]            [-m, --mjd MJD]         [--tsamp TIME]")
    print("  [--chandra]")
    print()
    print("Graphics Options:")
    print("  [--fontsize SIZE]         [--color COLOR]         [--linestyle STYLE]")
    print("  [--linewidth WIDTH]       [--marker TYPE]         [--markercolor COLOR]")
    print("  [--markerwidth WIDTH]     [--markersize SIZE]     [--facecolor COLOR]")
    print("  [--cmap COLORMAP]         [-c, --colorbar]        [--title STR]")
    print("  [--legend STR]            [--loc STR]             [-g, --grid]")
    print("  [--titlepos STR]          [--label STR]           [--labelpos STR]")
    print("  [--no-top-axis]           [--no-right-axis]")
    if kapteyn_loaded:
        print("  [--cmap-scale SCALE]")


###################################################################################################################################
#                                                      M  A  I  N 
###################################################################################################################################
if __name__=="__main__":

        #
        # Parsing the command line options
        #
    usage = "Usage: %prog [options] <.dat OR .tim (use --tim option)>"
    cmdline = opt.OptionParser(usage)
    cmdline.add_option('-n', '--nbins', dest='nbins', metavar='#BINS',
                       help="number of phase bins per pulsar period (default: number of samples)", default=-1, type='int')
    cmdline.add_option('-p', '--pulsar', dest='psrname', metavar='NAME',
                       help="pulsar name to be used for polyco instead of inf-file", type='str')
    cmdline.add_option('--polyco', dest='polycofile', metavar='FILE',
                       help="polyco file to be used for folding (default: %default)", default="polyco.dat", type='str')
    cmdline.add_option('--period', dest='period', metavar='PERIOD',
                       help="period in ms for folding. If not given then will use 'polyco.dat'", default=-1, type='float')
    cmdline.add_option('-s', '--startphase', dest='phase_start', metavar='PHASE',
                       help="start phase to plot", default=0.0, type='float')
    cmdline.add_option('-e', '--endphase', dest='phase_end', metavar='PHASE',
                       help="end phase to plot", default=1.0, type='float')
    cmdline.add_option('--start', dest='start_time', metavar='TIME',
                       help="time offset from the start in seconds (default: %default)", default=0, type='float')
    cmdline.add_option('-w', '--window', dest='window_time', metavar='TIME',
                       help="duration of the window in seconds (default: whole file)", default=-1, type='float')
    cmdline.add_option('-b', '--block', dest='blocksize', metavar='SAMPLES',
                       help="size of the block for reading the dat-file when folding profile. When using polyco file "
                            "period is updated for every block. So, for very close binaries this block size "
                            "should be probably very small. The smallest is 3000. Default: %default", default=1000000, type='int')
    cmdline.add_option('-f', '--fast-detrend', action="store_true", dest="fast_detrend",
                           help="Use a faster method of de-trending the time-series (2x speedup). "
                                "Sometimes (for strong pulsars (?), many strong pulses of which could effect calculation "
                                "of the linear trend) this method produces even better results when only median is subtracted", default=False)
    cmdline.add_option('--no-detrend', action="store_true", dest="is_no_detrend",
                       help="do not detrend the data", default=False)
    cmdline.add_option('-a', '--rebin', dest='rebin', metavar='FACTOR',
                       help="averaging time series by FACTOR (default: %default)", default=1, type='int')
    cmdline.add_option('-t', '--timeseries', action="store_true", dest="is_timeseries",
                       help="no folding. Time series will be plotted", default=False)
    cmdline.add_option('-y', '--profileonly', action="store_true", dest="is_profileonly",
                       help="only plot the average profile. No stack of pulses or subints", default=False)
    cmdline.add_option('-2', '--double', action="store_true", dest="is_period_doubled",
                       help="plot doubled-period profile (only when options -s and -e are not given)", default=False)
    cmdline.add_option('-k', '--stacking', action="store_true", dest="is_stacking",
                       help="Plot series of pulses or subints in stacking mode. Default is grey-scale mode", default=False)
    cmdline.add_option('-d', '--dump', dest='dump', metavar='#PULSES | TIME',
                       help="number of pulses or time (if . given) to dump for subintegrations", default='', type='str')
    cmdline.add_option('--offset', dest='offset', metavar='OFFSET',
                       help="Offset between individual profiles in stacking mode. Default = %default. "
                            "Offset is in the same units as for profiles' flux density. "
                            "Only positive (>=0) offsets are allowed", default=1.0, type='float')
    cmdline.add_option('--saveprof', dest='proffile', metavar='FILE',
                       help="save profile to binary file FILE", default='', type='str')
    cmdline.add_option('--saveonly', action="store_true", dest="is_saveonly",
                       help="only saves png-file and exits", default=False)
    cmdline.add_option('-i', '--image', dest='imageext', metavar='FILEEXT',
                       help="image file extension when used with --saveonly (default: %default)", default='png', type='str')
    cmdline.add_option('--tim', action="store_true", dest="is_timfile",
                       help="input file is Sigproc-style tim-file. None inf-file is necessary in this case", default=False)
    cmdline.add_option('-l', '--list', action="store_true", dest="is_printlist",
                       help="print short list of all options", default=False)
    cmdline.add_option('--events', action="store_true", dest="is_events",
                       help="input file is ascii file with event time in secs. Must use --mjd option to provide the start MJD. "
                            "Number of bins _must_ be also given with --nbins option unless --timeseries is set", default=False)
    cmdline.add_option('-m', '--mjd', dest='mjd', metavar='MJD',
                       help="start MJD of the data. By default, is read from .inf file or tim-file header. "
                            "If events file is used, it _must_ be given or --chandra is used", default='', type='str')
    cmdline.add_option('--chandra', action="store_true", dest="is_chandra",
                       help="events file is Chandra file, so start MJD is set to 50814.0 (Chandra reference MJD)", default=False)
    cmdline.add_option('--tsamp', dest='tsamp', metavar='TIME',
                       help="sampling time in sec of the data. By default, is read from .inf file or tim-file header", default='', type='str')
    group = opt.OptionGroup(cmdline, "Graphics Options")
    group.add_option('--fontsize', dest='fs', metavar='SIZE',
                       help="font size for labels (default: %default)", default=10, type='int')
    group.add_option('--color', dest='color', metavar='COLOR',
                       help="line and marker color (default: %default)", default='green', type='str')
    group.add_option('--linestyle', dest='linestyle', metavar='STYLE',
                       help="set linestyle. Default is '%default' (solid). "
                            "Other possible values: '--', '-.', ':', 'None', ' ', '' with possible combination "
                            "with 'steps', 'steps-pre', 'steps-mid', or 'steps-post'. In stacking mode only "
                            "'-', '--', '-.', ':' are possible, or their aliases: 'solid', 'dashed', 'dashdot', "
                            "and 'dotted'", default='-', type='str')
    group.add_option('--linewidth', dest='linewidth', metavar='WIDTH',
                       help="set linewidth. Default is '%default'", default='1.0', type='float')
    group.add_option('--marker', dest='marker', metavar='TYPE',
                       help="set line and marker color (default: %default). Use 'None' to not use the marker symbol", default=',', type='str')
    group.add_option('--markercolor', dest='markercolor', metavar='COLOR',
                       help="set only marker color (default: %default)", default='green', type='str')
    group.add_option('--markerwidth', dest='markerwidth', metavar='WIDTH',
                       help="set marker width (default: %default)", default='1.0', type='float')
    group.add_option('--markersize', dest='markersize', metavar='SIZE',
                       help="set marker size (default: %default)", default='6.0', type='float')
    group.add_option('--facecolor', dest='facecolor', metavar='COLOR',
                       help="marker facecolor or fill color in stacking mode (default: %default)", default='white', type='str')
    group.add_option('--cmap', dest='colormap', metavar='COLORMAP',
                      help="set colormap for plotting pulse stack. Default: %default. " 
                           "Other good colormaps: gray, gist_yarg (reversed gray), "
                           "gist_stern, hot, jet, pink, gist_heat, gist_gray, copper, ...", default='gist_earth', type='str')
    if kapteyn_loaded:
        group.add_option('--cmap-scale', dest='cmap_scaling', metavar='SCALE',
                              help="set colormap scaling. Default: %default. "
                                   "Other possible values: SQRT, LOG, EXP, SQUARE", default='LINEAR', type='str')
        group.add_option('-c', '--colorbar', action="store_true", dest="is_colorbar",
                           help="plot colorbar when plotting pulse stack", default=False)
        group.add_option('--title', dest='title', metavar='STR',
                           help="set title of the plot. Use '\\n' to have several lines", default='', type='str')
        group.add_option('--titlepos', dest='titlepos', metavar='STR',
                           help="set the position of the title in plot coordinates from 0 till 1. "
                                "Default is '%default', i.e. x=0.5, y=1.07 (top of the plot), and horizontal "
                                "alignment is 'center'", default='0.5,1.07,center', type='str')
        group.add_option('--legend', dest='legend', metavar='STR',
                           help="set label for the legend when only used with --timeseries or --profileonly options. "
                                "Use '\\n' to have several lines", default='', type='str')
        group.add_option('--loc', dest='loc', metavar='STR',
                           help="location of the legend. Default is '%default'. Other possible values: 'upper left', "
                                "'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', "
                                "'upper center', 'center'.", default='upper right', type='str')
        group.add_option('--label', dest='label', metavar='STR',
                           help="set label to put anywhere in the figure. The position is determined by --labelpos option. "
                                "Use '\\n' to have several lines", default='', type='str')
        group.add_option('--labelpos', dest='labelpos', metavar='STR',
                           help="set the position of the label in plot coordinates from 0 till 1. "
                                "Default is '%default', i.e. x=0.05, y=0.9, and horizontal "
                                "alignment is 'left'", default='0.05,0.9,left', type='str')
        group.add_option('-g', '--grid', action="store_true", dest="is_grid",
                           help="turn on the grid lines", default=False)
        group.add_option('--no-top-axis', action="store_true", dest="is_no_top_axis",
                           help="turn off the different labeling on the top axis", default=False)
        group.add_option('--no-right-axis', action="store_true", dest="is_no_right_axis",
                           help="turn off the different labeling on the right axis", default=False)
    cmdline.add_option_group(group)

    # reading cmd options
    (opts,args) = cmdline.parse_args()

    # print short list of options
    if opts.is_printlist:
        list_options(sys.argv[0].split("/")[-1])    
        sys.exit(0)

    # check if input file is given
    if len(args) != 0:
        datfile = args[0]
    else:
        #cmdline.print_help()
        cmdline.print_usage()
        sys.exit(0)

    # importing matplotlib
    if opts.is_saveonly:
        import matplotlib
        matplotlib.use("Agg")
        pngname = datfile.split(".dat")[0] + "." + opts.imageext
    else:
        import matplotlib

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.cm as cm
    import matplotlib.collections as collections
    import matplotlib.font_manager as fm

    headersize = 0  # size in bytes of the header
    # if input file is Sigproc-style tim-file
    if opts.is_timfile:
        try:
            filhdr, headersize = sigproc.read_header(datfile)
            startmjd = filhdr['tstart']
            tsamp = filhdr['tsamp']
            source = filhdr['source_name']
        except:
            print("Error: Can't open the tim-file '%s'!" % (datfile,))
            sys.exit(1)

    elif opts.is_events:
        opts.is_no_detrend = True  # we don't do detrending for events
        if opts.mjd == '' and not opts.is_chandra:
            print("Error: for events' file start MJD _must_ be given with --mjd option or --chandra option!")
            sys.exit(1)
        if opts.nbins == -1:
            print("Error: number of bins _must_ be given with --nbins option!")
            sys.exit(1)
        if opts.rebin != 1:
            print("Event data can not be re-binned")
            opts.rebin = 1
    else:
        if (opts.mjd == '' and not opts.is_chandra) or opts.tsamp == '' or opts.psrname == '':
            # reading inf-file to get corresponding info
            inffile = datfile.split(".dat")[0] + ".inf"
            try:
                id = inf.infodata(inffile)
                tsamp = id.dt        # sampling time
                startmjd = id.epoch  # start MJD
                source = id.object   # pulsar name
            except:
                print("Error: Can't read the inf-file '%s'!" % (inffile,))
                sys.exit(1)

    # overwriting MJD, tsamp, pulsarname from the command line if given
    if opts.is_chandra:
        opts.mjd = "50814.0"
        print("Chandra event file. Reference start MJD is %s" % (opts.mjd))
    if opts.mjd != '':
        startmjd = float(opts.mjd)
    if opts.tsamp != '':
        tsamp = float(opts.tsamp)
    if opts.psrname != '':
        source = opts.psrname

    # checking start, end times and adjusting them if necessary
    if opts.start_time < 0:
        print("Error: start time %.3f is < 0!" % (opts.start_time))
        sys.exit(1)

    if not opts.is_events:
        # getting the size of the file
        try:
            size=(os.stat(datfile)[stat.ST_SIZE] - headersize) / 4   # 4 bytes in float
        except:
            print("Error: Can't open the input file '%s'!" % (datfile,))
            sys.exit(1)

        # checking start, end times and adjusting them if necessary
        start_sample = int(opts.start_time / tsamp)
        if start_sample >= size:
            print("Start time %.3f is out of range!" % (opts.start_time))
            sys.exit(1)
        else:
            size = size - start_sample
    
        # adjusting start MJD
        opts.start_time = start_sample * tsamp
    
        if opts.window_time != -1:
            window = int(opts.window_time / tsamp)
            opts.window_time = window * tsamp
        else:
            window = size
            opts.window_time = window * tsamp

        if start_sample + window > size:
            print("End time %.3f is out of range. The window duration will be adjusted" % (opts.start_time + opts.window_time))
            window = size - start_sample
            opts.window_time = window * tsamp
            print("New window is %.3f s" % (opts.window_time))

    else:  # reading the file, and so, we can get the start and end events' time
        events = np.loadtxt(datfile, comments='#', usecols=(0,0), dtype=float, unpack=True)[0]
        try:
            energy = np.loadtxt(datfile, comments='#', usecols=(1,1), dtype=float, unpack=True)[0]
        except:
            energy = []
        if np.size(energy) == 0 and opts.is_timeseries:
            print("No energy column is given in the events'file '%s'. It _must_ be given to plot the timeseries!" % (datfile))
            sys.exit(1)
        duration = events[-1]-events[0]
        if opts.start_time >= duration:
            print("Start time %.3f sec is more than duration of observation of %.3f sec!" % (opts.start_time, duration))
            sys.exit(1)
        else:
            events.compress(events >= opts.start_time + events[0])
            if opts.is_timeseries:
                energy.compress(events >= opts.start_time + events[0])
        if opts.window_time == -1:
            opts.window_time = duration
        if opts.start_time + opts.window_time > duration:
            print("End time %.3f is out of range. The window duration will be adjusted" % (opts.start_time + opts.window_time))
            opts.window_time = duration - opts.start_time
            print("New window is %.3f s" % (opts.window_time))

    # checking dump settings
    if opts.dump != '':
        if "." in opts.dump:  # dump time is given
            dump_time = float(opts.dump)
            dump_pulses = 0
            if dump_time > 0:
                print("Dump time is %.3f s" % (dump_time))
            else:
                print("Dump time %.3f should be > 0!" % (dump_time))
                sys.exit(1)
        else:  # number of pulses to dump is given
            dump_pulses = int(opts.dump)
            dump_time = 0    
            if dump_pulses > 0:
                print("Number of pulses in subintegration is %d" % (dump_pulses))
            else:
                print("Number of pulses in subintegration %d should be > 0!" % (dump_pulses))
                sys.exit(1)
    else:
        dump_pulses = 0
        dump_time = 0

    # Checking the phases and correct if necessary
    if not opts.is_timeseries:
        if opts.phase_start < 0 or opts.phase_start >= 1.0:
            opts.phase_start = 0.0
            print("Start phase is out of range. Adjusted value is %.3f" % (opts.phase_start,))

        if opts.phase_end <= 0.0 or opts.phase_end > 1.0:
            opts.phase_end = 1.0
            print("End phase is out of range. Adjusted value is %.3f" % (opts.phase_end,))

        if opts.phase_end <= opts.phase_start:
            print("End phase %.3f is <= than start phase %.3f! Please adjust." % (opts.phase_end, opts.phase_start))
            sys.exit(1)

        if not opts.is_no_detrend:
            # checking the blocksize
            if opts.blocksize < 3 * detrendlen:
                print("Block size is too small: %d. Will be increased to %d." % (opts.blocksize, 3 * detrendlen))
                opts.blocksize = 3 * detrendlen
            # making blocksize to be divisible by detrendlen (in order not to loose samples between blocks)
            if int (opts.blocksize / detrendlen) * detrendlen != opts.blocksize:
                opts.blocksize = detrendlen * (1 + int (opts.blocksize / detrendlen))
                print("Adjusting block size to %d to be divisible by detrendlen=%d" % (opts.blocksize, detrendlen))

    
    # forming the array of time samples
    if not opts.is_events:
        try:
            infile = open(datfile, "rb")
        except:
            print("Error: Can't read the dat-file '%s'!" % (datfile,))
            sys.exit(1)
        dataptr = ar.array('f')  # 'f' - for float
        infile.seek(headersize + 4 * start_sample)  # position to the first byte to read; '4' - is the size of float 
    else:
        data = events / 86400.  # converting events to days
        data += startmjd        # converting to MJD
        events -= events[0]     # converting events' time relative to the start of observation
    

    # Folding the profile
    if not opts.is_timeseries:
        if opts.period == -1:   # Period is not given in the cmdline, so will use polyco file
            pid=poly.polycos(source, opts.polycofile)
            try:
                if not opts.is_events:
                            fold_period = get_period (pid, startmjd)
                else:
                            fold_period = get_period (pid, data[0])
            except:
                print("Check the name of the pulsar in polyco file '%s' and inf-file '%s'!" % (opts.polycofile, inffile))
                print("If different, try --pulsar option to set the name of pulsar the same as in polyco file.")
                sys.exit(1)
            is_update_period = True
            if fold_period <= 0:
                print("Computed fold period is bad: %f. Check your polyco and/or MJD!" % (float(fold_period)))
                sys.exit(1)
        else: # period is given
            fold_period = opts.period / 1000.
            is_update_period = False
        # if number of bins is not specified  (should always be specified for events)
        if opts.nbins == -1:
            opts.nbins = int(fold_period / tsamp)
        # if dump_time is given - checking that it is >= fold_period
        if dump_time > 0:
            if dump_time < fold_period:
                print("Dump time %.3f s is less than folding period of %f s. Adjusting to match." % (dump_time, fold_period))
                dump_time = fold_period
                print("Dump time is now %.3f s" % (dump_time)) 
            if dump_time > opts.window_time:
                print("Dump time %.3f is more than window duration of %f s. Adjusting..." % (dump_time, opts.window_time))
                # we make it a one period less than duration, because otherwise plt.imshow
                # fails to plot
                dump_time = opts.window_time - fold_period
                print("Dump time is now %.3f s" % (dump_time)) 
        else: # we need this for plotting purposes
            dump_time = dump_pulses * fold_period
            if dump_time > opts.window_time:
                print("Number of pulses per subintegration %d is more than within window duration of %f s." % (dump_pulses, opts.window_time))
                print("Adjusting to match to the closest maximum possible number.")
                dump_pulses = int((opts.window_time - fold_period) / fold_period)
                dump_time = dump_pulses * fold_period
                print("Number of pulses per subintegration is now %d" % (dump_pulses))

        bin_start = int (opts.nbins * opts.phase_start)
        bin_end = int(math.ceil((opts.nbins - 1) * opts.phase_end))
        bin_range = bin_end - bin_start + 1
        if "%f" % (opts.phase_start) != "%f" % (float(bin_start)/opts.nbins):
            opts.phase_start = float(bin_start)/opts.nbins
            print("Adjusting the start phase to %.3f to account for integer number of bins" % (opts.phase_start))
        if "%f" % (opts.phase_end) != "%f" % (float(bin_end)/(opts.nbins - 1)):
            opts.phase_end = float(bin_end)/(opts.nbins - 1)
            print("Adjusting the end phase to %.3f to account for integer number of bins" % (opts.phase_end))

        # initializing some variables and arrays
        elapsed_secs = opts.start_time
        elapsed_turns = 0.
        ibin = 0
        lbin = -1
        pulsecount = 0
        profile = np.zeros (bin_range, dtype=float)   # profile array
        counts = np.zeros (bin_range, dtype=float)    # array that keeps number of counts in each bin of the profile
            
        if not opts.is_profileonly: # initialize the 2D array to keep the stack of pulses
            if dump_pulses > 0:
                npulses_expect = int(math.ceil(int((opts.window_time) / fold_period) / float(dump_pulses)))
            elif dump_time > 0:
                npulses_expect = int(math.ceil((opts.window_time) / dump_time))
            else:
                npulses_expect = int(math.ceil((opts.window_time) / fold_period))
            # increase that value by 25% (overkill probably) in order to be safe if period is changing significantly
            # over the course of observations
            npulses_expect += int(0.25 * npulses_expect)
            if npulses_expect == 1:
                npulses_expect += 1
            pulsestack = np.zeros((npulses_expect, bin_range), dtype=float)
            if dump_pulses > 0 or dump_time > 0:
                dump = np.zeros (bin_range, dtype=float)   # current subintegration
                dumpcount = 0  # number of pulses in current subintegration
                dump_next = elapsed_secs + dump_time  # time of the next dump
                dumpcounts = np.zeros (bin_range, dtype=float)    # array that keeps number of counts in each bin of the subintegration

        # in case of events' file we always use dump_time
        if opts.is_events:
            if dump_time == 0:  # always use dump_time
                dump_time = fold_period
            dumpcount = 0  # number of pulses in current subintegration
            dump_next = elapsed_secs + dump_time  # time of the next dump
            dumpcounts = np.zeros (bin_range, dtype=float)    # array that keeps number of counts in each bin of the subintegration



        if not opts.is_events:
            # Loop over the chunks of data to read
            while 1:
                samplesleft = (headersize + (start_sample + window) * 4 - infile.tell()) / 4
                if samplesleft <= 0:  # leaving from this infinite while loop
                    break 
                if samplesleft > opts.blocksize:
                    dataptr.fromfile(infile, opts.blocksize)
                else:
                    dataptr.fromfile(infile, samplesleft)
                data = np.array(dataptr)
                dataptr = ar.array('f')  # clear the array. This is important!
                readsamples = np.size(data)  # how many samples we have read

                # detrending the data
                if not opts.is_no_detrend: 
                    # this check is necessary when reading the very last block and if its size is _very_ small
                    # then we just have to skip it when detrending
                    if readsamples < 3 * detrendlen:
                        break
                    data = detrending (data, opts.fast_detrend)
                size = np.size(data)

                # updating the period if necessary
                if is_update_period:
                    fold_period = get_period (pid, startmjd + elapsed_secs / 86400.)
                
                # main folding loop over the read samples
                turns = elapsed_turns
                secs = elapsed_secs
                for s in range(0, size):
                    # phase of each sample
                    phase=turns - math.floor(turns)
                    if phase >= opts.phase_start and phase <= opts.phase_end:
                        ibin = int (opts.nbins * phase) - bin_start
                        if ibin == opts.nbins:  ibin -= 1
                        if ibin - lbin < 0: pulsecount += 1
                        counts[ibin] += 1.
                        profile[ibin] += data[s]
                        if not opts.is_profileonly:
                            # making the subintegration
                            if dump_pulses > 0 or dump_time > 0:
                                # if number of pulses to dump is given
                                if dump_pulses > 0 and ibin - lbin < 0 and pulsecount != 0 and pulsecount%dump_pulses == 0:
                                    pulsestack[dumpcount] = np.array([dump[i]/(dumpcounts[i] == 0.0 and 1.0 or dumpcounts[i]) for i in range(0, bin_range)], dtype=float)
                                    dumpcount += 1
                                    dumpcounts[:] = 0.0
                                    dump[:] = 0.0
                                # if dump time is given
                                if dump_time > 0 and dump_pulses <= 0 and dump_next - secs <= tsamp/2.:
                                    pulsestack[dumpcount] = np.array([dump[i]/(dumpcounts[i] == 0.0 and 1.0 or dumpcounts[i]) for i in range(0, bin_range)], dtype=float)
                                    dumpcount += 1
                                    dumpcounts[:] = 0.0
                                    dump[:] = 0.0
                                    dump_next = secs + dump_time

                                # after the dumps (if happened) we still need to continue with the current sample 
                                # that belongs already to the next dump
                                dumpcounts[ibin] += 1.
                                dump[ibin] += data[s]
                            else: # no subintegrations
                                pulsestack[pulsecount][ibin] = data[s]         

                    turns += (tsamp / fold_period)
                    secs += tsamp
                    lbin = ibin
                elapsed_secs += readsamples * tsamp
                elapsed_turns += (readsamples * tsamp / fold_period)

        else:   # if events are given
            for tt in range(np.size(data)):
                # updating the period if necessary
                if is_update_period:
                    fold_period = get_period (pid, data[tt])
                turns = events[tt] / fold_period
                # phase of each sample
                phase=turns - math.floor(turns)
                if phase >= opts.phase_start and phase <= opts.phase_end:
                    ibin = int (opts.nbins * phase) - bin_start
                    if ibin == opts.nbins:  ibin -= 1
                    counts[ibin] += 1.
                    if not opts.is_profileonly:
                        # if dump time is given (always use in case of events file)
                        if ((tt < np.size(data)-1 and dump_next < events[tt+1]) or tt == np.size(data)-1):
                            pulsestack[dumpcount] = dumpcounts
                            dumpcount += 1
                            dumpcounts[:] = 0.0
                            dump_next += dump_time
                        # after the dumps (if happened) we still need to continue with the current sample 
                        # that belongs already to the next dump
                        dumpcounts[ibin] += 1.


        # normalizing the profile
        if not opts.is_events:
            profile = np.array([profile[i]/(counts[i] == 0.0 and 1.0 or counts[i]) for i in range(0, bin_range)], dtype=float)
        else:
            profile = counts
        # saving the profile to binary file
        if opts.proffile != '':
            outp = open(opts.proffile, 'wb')    
            outbins = ar.array('f')
            outbins.fromlist(profile.tolist())
            outbins.tofile(outp)
            outp.close()


    else:   
        if not opts.is_events:
            # if one wants just to plot time series, just read the whole selected block and then plot it
            # If the size of window is too large, then probably Python won't be able to read it all at once
            dataptr.fromfile(infile, window)
            data = np.array(dataptr)
            dataptr = ar.array('f')  # clear the array. Though, here it's not important or necessary
            # detrending the data
            if not opts.is_no_detrend: 
                data = detrending (data, opts.fast_detrend)
            window = np.size(data)
            opts.window_time = window * tsamp

    #
    # Plotting
    #
        fig = plt.figure()

    if opts.is_timeseries:            # plotting the time series
        if opts.rebin == 1:
            if not opts.is_events:
                flux = data
            else:   
                flux = energy
        else:
            window = int(window/opts.rebin)
            tsamp *= opts.rebin
            flux = [np.average(data[k*opts.rebin:(k+1)*opts.rebin]) for k in range(0, window)]

        if not opts.is_events:
            time = [opts.start_time + tsamp * t  for t in range(0, window)]
        else:
            time = events

            ax = fig.add_subplot(111)
            plt.xlabel("Time (s)", fontsize=opts.fs)
        if opts.is_events:
                plt.ylabel("Energy", fontsize=opts.fs)
        elif opts.is_no_detrend:
                plt.ylabel("Flux density (arb. units)", fontsize=opts.fs)
        else:
            plt.ylabel("Flux density ($\sigma$)", fontsize=opts.fs)

        ax.plot (time, flux, color="%s" % (opts.color), marker="%s" % (opts.marker),
                 markeredgecolor="%s" % (opts.markercolor), markerfacecolor="%s" % (opts.facecolor),
                 linestyle="%s" % (opts.linestyle), linewidth=opts.linewidth, markeredgewidth=opts.markerwidth,
                 markersize=opts.markersize, label="%s" % ("\n".join(opts.legend.split("\\n"))))
        # turn on grid
        if opts.is_grid:
            plt.grid(True)

        if not opts.is_no_top_axis:
            if not opts.is_events:
                axtop = plt.twiny()
                axtop.xaxis.tick_top()
                axtop.xaxis.set_label_position("top")
                axtop.set_xlim(xmin=start_sample, xmax=start_sample+window-1)
                for label in axtop.get_xticklabels(): label.set_fontsize(opts.fs)
                plt.xlabel("Samples", fontsize=opts.fs)
                plt.gca().minorticks_on()

        ax.set_xlim(xmin=time[0], xmax=time[-1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
        for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)

    elif opts.is_profileonly:   # plotting only the profile

        flux = profile        
        phase = [float(n)/opts.nbins for n in range(bin_start, bin_end + 1)]
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            flux = np.append(flux, profile)
            phase = np.append(phase, [float(opts.nbins + n)/opts.nbins for n in range(bin_start, bin_end + 1)])

            ax = fig.add_subplot(111)
            plt.xlabel("Pulse phase", fontsize=opts.fs)
        if opts.is_events:
                plt.ylabel("Counts/bin", fontsize=opts.fs)
        elif opts.is_no_detrend:
                plt.ylabel("Flux density (arb. units)", fontsize=opts.fs)
        else:
            plt.ylabel("Flux density ($\sigma$)", fontsize=opts.fs)

        ax.plot (phase, flux, color="%s" % (opts.color), marker="%s" % (opts.marker),
                 markeredgecolor="%s" % (opts.markercolor), markerfacecolor="%s" % (opts.facecolor),
                 linestyle="%s" % (opts.linestyle), linewidth=opts.linewidth, markeredgewidth=opts.markerwidth,
                 markersize=opts.markersize, label="%s" % ("\n".join(opts.legend.split("\\n"))))
        # turn on grid
        if opts.is_grid == True:
            plt.grid(True)

        if not opts.is_no_top_axis:
            axtop = plt.twiny()
            axtop.xaxis.tick_top()
            axtop.xaxis.set_label_position("top")
            if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
                axtop.set_xlim(xmin=bin_start, xmax=opts.nbins + bin_end)
            else:
                axtop.set_xlim(xmin=bin_start, xmax=bin_end)
                for label in axtop.get_xticklabels(): label.set_fontsize(opts.fs)
                plt.xlabel("Phase bin", fontsize=opts.fs)
            plt.gca().minorticks_on()

        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            ax.set_xlim(xmin=opts.phase_start, xmax=2.*opts.phase_end)
        else:
            ax.set_xlim(xmin=opts.phase_start, xmax=opts.phase_end)
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
            for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)

    elif opts.is_stacking:   # plotting the stack of pulses or subints

            ax = fig.add_subplot(111)
        if not opts.is_no_top_axis:
                plt.xlabel("Phase bin", fontsize=opts.fs)
        if not opts.is_no_right_axis:
            if opts.is_events:
                    plt.ylabel("Counts/bin", fontsize=opts.fs)
            elif opts.is_no_detrend:
                    plt.ylabel("Flux density (arb. units)", fontsize=opts.fs)
            else:
                plt.ylabel("Flux density ($\sigma$)", fontsize=opts.fs)

        dr = abs(opts.offset)
        if dump_pulses > 0 or dump_time > 0:
            ncount = dumpcount
        else:
            ncount = pulsecount
        ymin = pulsestack[0:ncount].min()
        ymax = pulsestack[0:ncount].max()
        dr = abs(opts.offset)
        dmins = []
        dmaxs = []

        t = np.arange(bin_start, bin_end + 1, dtype=float)            
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            t = np.append(t, [opts.nbins + n for n in range(bin_start, bin_end + 1)])
        pulseverts = []
        for i in np.arange(ncount-1, -1, -1):
            temp = pulsestack[i] + i * dr
            if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
                temp = np.append(temp, temp)
            dmins.append(temp.min())
            dmaxs.append(temp.max())
            pulseverts.append(np.vstack((np.hstack((t[:,np.newaxis], temp[:,np.newaxis])), [t[-1], ymin-1], [0, ymin-1])))

        ymin = np.array(dmins).min()
        ymax = np.array(dmaxs).max()
        polys = collections.PolyCollection(pulseverts, closed=True)
        polys.set_edgecolor("%s" % (opts.color))
        polys.set_facecolor("%s" % (opts.facecolor))
        polys.set_linestyle("%s" % (opts.linestyle))
        polys.set_linewidth("%s" % (opts.linewidth))
        ax.add_collection(polys)

        plt.gca().minorticks_on()

        # turn on grid
        if opts.is_grid:
            plt.grid(True)

        axbot = plt.twiny()
        axbot.xaxis.tick_bottom()
        axbot.xaxis.set_label_position("bottom")
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            axbot.set_xlim(xmin=opts.phase_start, xmax=2.*opts.phase_end)
        else:
            axbot.set_xlim(xmin=opts.phase_start, xmax=opts.phase_end)
            for label in axbot.get_xticklabels(): label.set_fontsize(opts.fs)
            plt.xlabel("Pulse phase", fontsize=opts.fs)
        plt.gca().minorticks_on()

        ayleft = plt.twinx()
        ayleft.yaxis.tick_left()
        ayleft.yaxis.set_label_position("left")
        if (dump_pulses > 0 or dump_time > 0) and not opts.is_events:
            ayleft.set_ylim(ymin=0.0, ymax=ncount*dump_time)
        else:
            ayleft.set_ylim(ymin=0.0, ymax=opts.window_time)
            for label in ayleft.get_yticklabels(): label.set_fontsize(opts.fs)
            plt.ylabel("Observing time (s)", fontsize=opts.fs)
        plt.gca().minorticks_on()

        # Determining the main (top/right) axes
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            ax.set_xlim(xmin=bin_start, xmax=opts.nbins + bin_end-0.1)
        else:
            ax.set_xlim(xmin=bin_start, xmax=bin_end-0.1)
        ax.set_ylim(ymin=ymin, ymax=ymax)
        if not opts.is_no_top_axis:
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
            plt.gca().minorticks_on()
        else:   # if don't want to show then turning everything off (I can not duplicate the limits and
            # ticks from the bottom axis because then imshow won't work properly
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.xaxis.set_minor_locator(ticker.NullLocator())

        if not opts.is_no_right_axis:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
            for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)
            plt.gca().minorticks_on()
        else:   # if don't want to show then turning everything off (I can not duplicate the limits and
            # ticks from the left axis because then imshow won't work properly
            ax.yaxis.set_major_formatter(ticker.NullFormatter())
            ax.yaxis.set_major_locator(ticker.NullLocator())
            ax.yaxis.set_minor_locator(ticker.NullLocator())


    else:   # plotting the stack of pulses or subints in grey-scale mode

        ax = fig.add_subplot(111)
        if not opts.is_no_top_axis:
            plt.xlabel("Phase bin", fontsize=opts.fs)

        if not kapteyn_loaded:
            colormap = cm.get_cmap(opts.colormap)
        else:
            colormap=VariableColormap(cm.get_cmap(opts.colormap))
            colormap.set_scale(opts.cmap_scaling)

        if dump_pulses > 0 or dump_time > 0:
            ncount = dumpcount
            if not opts.is_no_right_axis: plt.ylabel("Sub-integration", fontsize=opts.fs)
        else:
            ncount = pulsecount
            if not opts.is_no_right_axis: plt.ylabel("Pulse number", fontsize=opts.fs)

        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            pulsestack = np.concatenate((pulsestack[:ncount], pulsestack[:ncount]), axis=1)
            cax = ax.imshow(pulsestack[:ncount], interpolation=None, aspect='auto', extent=(bin_start, opts.nbins + bin_end, 0, ncount), origin='lower', cmap=colormap)
        else:
            cax = ax.imshow(pulsestack[:ncount], interpolation=None, aspect='auto', extent=(bin_start, bin_end, 0, ncount), origin='lower', cmap=colormap)
        if opts.is_colorbar:
            cbar = fig.colorbar(cax, orientation='horizontal', spacing='proportional')
            if opts.is_events:
                    cbar.ax.set_xlabel("Counts/bin", fontsize=opts.fs)
            elif opts.is_no_detrend:
                    cbar.ax.set_xlabel("Flux density (arb. units)", fontsize=opts.fs)
            else:
                    cbar.ax.set_xlabel("Flux density ($\sigma$)", fontsize=opts.fs)
            for label in cbar.ax.get_xticklabels(): label.set_fontsize(opts.fs)

        plt.gca().minorticks_on()

        # turn on grid
        if opts.is_grid == True:
            plt.grid(True)

        axbot = plt.twiny()
        axbot.xaxis.tick_bottom()
        axbot.xaxis.set_label_position("bottom")
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            axbot.set_xlim(xmin=opts.phase_start, xmax=2.*opts.phase_end)
        else:
            axbot.set_xlim(xmin=opts.phase_start, xmax=opts.phase_end)
            for label in axbot.get_xticklabels(): label.set_fontsize(opts.fs)
            plt.xlabel("Pulse phase", fontsize=opts.fs)
        plt.gca().minorticks_on()

        ayleft = plt.twinx()
        ayleft.yaxis.tick_left()
        ayleft.yaxis.set_label_position("left")
        if (dump_pulses > 0 or dump_time > 0) and not opts.is_events:
            ayleft.set_ylim(ymin=0.0, ymax=ncount*dump_time)
        else:
            ayleft.set_ylim(ymin=0.0, ymax=opts.window_time)
            for label in ayleft.get_yticklabels(): label.set_fontsize(opts.fs)
            plt.ylabel("Observing time (s)", fontsize=opts.fs)
        plt.gca().minorticks_on()

        # Determining the main (top/right) axes
        if opts.is_period_doubled and bin_start == 0 and bin_end == opts.nbins - 1:
            ax.set_xlim(xmin=bin_start, xmax=opts.nbins + bin_end)
        else:
            ax.set_xlim(xmin=bin_start, xmax=bin_end)
        ax.set_ylim(ymin=0, ymax=ncount)
        if not opts.is_no_top_axis:
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
                for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
            plt.gca().minorticks_on()
        else:   # if don't want to show then turning everything off (I can not duplicate the limits and
                        # ticks from the bottom axis because then imshow won't work properly
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.xaxis.set_minor_locator(ticker.NullLocator())

        if not opts.is_no_right_axis:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
                for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)
            plt.gca().minorticks_on()
        else:   # if don't want to show then turning everything off (I can not duplicate the limits and
                        # ticks from the left axis because then imshow won't work properly
            ax.yaxis.set_major_formatter(ticker.NullFormatter())
            ax.yaxis.set_major_locator(ticker.NullLocator())
            ax.yaxis.set_minor_locator(ticker.NullLocator())


    # Making the title
    if opts.title != '':
        ax.set_title("\n".join(opts.title.split("\\n")), fontsize=opts.fs, x=opts.titlepos.split(",")[0],
                     y=opts.titlepos.split(",")[1], ha='%s' % (opts.titlepos.split(",")[-1]))

    # Putting the label
    if opts.label != '':
        ax.annotate("\n".join(opts.label.split("\\n")), fontsize=opts.fs, xycoords='axes fraction', xy=(0,0),
                    xytext=(opts.labelpos.split(",")[0], opts.labelpos.split(",")[1]), ha='%s' % (opts.labelpos.split(",")[-1]))

    # turn on the Legend (only when is used in profile-only mode or timeseries mode)
    if opts.legend != '' and (opts.is_timeseries or opts.is_profileonly):
        prop=fm.FontProperties(size=opts.fs)
        plt.legend(prop=prop, loc="%s" % (opts.loc))

    # end of plotting
           if opts.is_saveonly: 
        plt.savefig(pngname)
    else:
        plt.show()

    # closing input-file
    if not opts.is_events:
        infile.close()
