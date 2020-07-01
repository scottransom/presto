from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
import numpy as np
from scipy.signal import medfilt
from presto import infodata
from presto.Pgplot import *

"""
class rfifind:
    Read in stats and mask information from an rfifind run.
    Use that information to make a list of bad channels to
    zap and/or weight and offset properly to (hopefully) improve
    your overall signal-to-noise ratio.

    A typical usage would be something like:
        In [1]: import rfifind

        In [2]: a = rfifind.rfifind('SPIGOT_040424_1614-22_rfifind.mask')

        In [3]: a.set_zap_chans(power=1000)
        Recommending to zap 75 channels:
          55 channels from Fourier power levels > 1000.0
          6 channels from avgs levels being off by > 5.0 sigma
          2 channels from stds levels being off by > 2.0 sigma
          20 channels for being within 0.01 of the band edge
          0 channels were specified by the user
          0 channels for having no variation
          The mask recommends 116 additional bad channels
          adding them to the zap_chans list.

              << Plot is shown here >>

           Type <RETURN> for next page:

        In [4]: a.write_zap_chans()

        In [5]: a.set_weights_and_offsets()

        In [6]: a.write_weights_and_offsets()
        
"""

class rfifind(object):
    def __init__(self, filename):
        self.basename = filename[:filename.find("_rfifind.")+8]
        self.idata = infodata.infodata(self.basename+".inf")
        self.read_stats()
        self.read_mask()
        self.get_bandpass()
        if len(self.goodints):
            self.get_median_bandpass()
            self.determine_padvals()

    def read_stats(self):
        x = open(self.basename+".stats")
        self.nchan, self.nint, self.ptsperint, self.lobin, self.numbetween = \
                    np.fromfile(x, dtype=np.int32, count=5)
        count = self.nchan * self.nint
        self.pow_stats = np.fromfile(x, dtype=np.float32, count=count)
        self.avg_stats = np.fromfile(x, dtype=np.float32, count=count)
        self.std_stats = np.fromfile(x, dtype=np.float32, count=count)
        self.pow_stats.shape = (self.nint, self.nchan)
        self.avg_stats.shape = (self.nint, self.nchan)
        self.std_stats.shape = (self.nint, self.nchan)
        x.close()

    def read_mask(self):
        x = open(self.basename+".mask")
        self.time_sig, self.freq_sig, self.MJD, self.dtint, self.lofreq, self.df = \
                       np.fromfile(x, dtype=np.float64, count=6)
        self.freqs = np.arange(self.nchan)*self.df + self.lofreq
        self.times = np.arange(self.nint)*self.dtint
        self.MJDs = self.times/86400.0 + self.MJD
        nchan, nint, ptsperint = np.fromfile(x, dtype=np.int32, count=3)
        nzap = np.fromfile(x, dtype=np.int32, count=1)[0]
        if nzap:
            self.mask_zap_chans = np.fromfile(x, dtype=np.int32, count=nzap)
        else:
            self.mask_zap_chans = np.asarray([])
        self.mask_zap_chans = set(self.mask_zap_chans)
        if len(self.mask_zap_chans)==self.nchan:
            print("WARNING!:  All channels recommended for masking!")
        nzap = np.fromfile(x, dtype=np.int32, count=1)[0]
        if nzap:
            self.mask_zap_ints = np.fromfile(x, dtype=np.int32, count=nzap)
        else:
            self.mask_zap_ints = np.asarray([])
        if len(self.mask_zap_ints)==self.nint:
            print("WARNING!:  All intervals recommended for masking!")
        nzap_per_int = np.fromfile(x, dtype=np.int32, count=nint)
        self.mask_zap_chans_per_int = []
        for nzap in nzap_per_int:
            if nzap:
                if nzap == nchan:
                    tozap = np.arange(0, nchan, dtype=np.int32)
                else:
                    tozap = np.fromfile(x, dtype=np.int32, count=nzap)
            else:
                tozap = np.asarray([])
            self.mask_zap_chans_per_int.append(tozap)
        x.close()

    def get_bandpass(self, plot=False):
        """
        get_bandpass():
            This routine returns a 'good' bandpass based on an average
                of the average bandpasses, with the exception of the
                intervals that were recommended for zapping in the mask.
        """
        ints = np.arange(self.nint)
        badints = self.mask_zap_ints
        goodints = set(ints) - set(badints)
        goodints = np.asarray(list(goodints))
        self.goodints = goodints
        if not len(goodints):
            print("WARNING!:  Cannot get bandpass because all intervals zapped.")
            return 0.0
        self.bandpass_avg = self.avg_stats[goodints,:].mean(0)
        self.bandpass_std = self.std_stats[goodints,:].mean(0)
        self.bandpass_pow = self.pow_stats[goodints,:].mean(0)
        if plot:
            plotxy(self.bandpass_avg, self.freqs, labx="Frequency (MHz)")
            plotxy(self.bandpass_avg+self.bandpass_std, self.freqs, color="red")
            plotxy(self.bandpass_avg-self.bandpass_std, self.freqs, color="red")
            closeplot()
        return self.bandpass_avg

    def get_median_bandpass(self, medlen=21, plot=False):
        self.median_bandpass_avg = medfilt(self.bandpass_avg, medlen)
        self.median_bandpass_std = medfilt(self.bandpass_std, medlen)
        if plot:
            plotxy(self.median_bandpass_avg, self.freqs,
                   labx="Frequency (MHz)")
            plotxy(self.median_bandpass_avg+self.median_bandpass_std,
                   self.freqs, color="red")
            plotxy(self.median_bandpass_avg-self.median_bandpass_std,
                   self.freqs, color="red")
            closeplot()
        return self.median_bandpass_avg

    def determine_padvals(self, frac_to_keep=0.8):
        """
        determine_padvals():
            This routines determines padding values to use for each
            channel.
        """
        # NOTE: Casting to 64/32-bit floats are done to mimick 'mask.c'.
        num = int(np.round(self.nint*frac_to_keep))
        start = (self.nint - num) // 2
        self.padvals = np.zeros(self.nchan, dtype='float32')
        for ichan in range(self.nchan):
            isort = np.argsort(self.avg_stats[:,ichan])
            self.padvals[ichan] = np.mean(self.avg_stats.astype('float64')[isort,ichan][start:start+num])

    def get_pow_zap_chans(self, pow_threshold=100.0):
        return np.where(self.bandpass_pow > pow_threshold)[0]

    def get_avg_zap_chans(self, sigma=5.0):
        resid = np.fabs(self.bandpass_avg - self.median_bandpass_avg)
        has_var = self.median_bandpass_std != 0.0
        no_var = self.median_bandpass_std == 0.0
        resid[has_var] /= self.median_bandpass_std[has_var]
        resid[no_var] = np.inf
        return np.where(resid > sigma)[0]

    def get_std_zap_chans(self, sigma=5.0):
        resid = np.fabs(self.bandpass_std - self.median_bandpass_std)
        has_var = self.median_bandpass_std != 0.0
        no_var = self.median_bandpass_std == 0.0
        resid[has_var] /= self.median_bandpass_std[has_var]
        resid[no_var] = np.inf
        return np.where(resid > sigma)[0]

    def get_no_signal_chans(self):
        med_std = np.median(self.median_bandpass_std)
        return np.where(self.bandpass_std < 0.15 * med_std)[0]

    def get_edge_chans(self, edges=0.01):
        nedge = int(self.nchan * edges)
        return np.concatenate((np.arange(nedge),
                               np.arange(nedge) + self.nchan-nedge))

    def set_zap_chans(self, power=100.0, asigma=5.0, ssigma=2.0, plot=True,
                      edges=0.01, usemask=True, chans=[]):
        self.user_zap_chans = set(chans)
        self.pow_zap_chans = set(self.get_pow_zap_chans(power))
        self.avg_zap_chans = set(self.get_avg_zap_chans(asigma))
        self.std_zap_chans = set(self.get_std_zap_chans(ssigma))
        self.no_signal_chans = set(self.get_no_signal_chans())
        self.edge_chans = set(self.get_edge_chans(edges))
        self.zap_chans = self.pow_zap_chans | \
                         self.avg_zap_chans | \
                         self.std_zap_chans | \
                         self.user_zap_chans | \
                         self.edge_chans | \
                         self.no_signal_chans
        print("Recommending to zap %d channels:" % len(self.zap_chans))
        print("  %d channels from Fourier power levels > %.1f" % \
              (len(self.pow_zap_chans), power))
        print("  %d channels from avgs levels being off by > %.1f sigma" % \
              (len(self.avg_zap_chans), asigma))
        print("  %d channels from stds levels being off by > %.1f sigma" % \
              (len(self.std_zap_chans), ssigma))
        print("  %d channels for being within %.3g of the band edge" % \
              (len(self.edge_chans), edges))
        print("  %d channels were specified by the user" % len(self.user_zap_chans))
        print("  %d channels for having no variation" % len(self.no_signal_chans))
        if usemask:
            onlymask = self.mask_zap_chans - self.zap_chans
            print("  The mask recommends %d additional bad channels" % len(onlymask))
            if len(onlymask):
                print("  adding them to the zap_chans list.")
            self.zap_chans = self.zap_chans | self.mask_zap_chans
        self.zap_chans = np.asarray(sorted(list(self.zap_chans)))
        if plot: self.plot_zapped_bandpass()

    def plot_zapped_chans(self, device="/xwin"):
        plotxy(self.bandpass_avg, self.freqs,
               labx="Frequency (MHz)", device=device)
        plotxy(self.median_bandpass_avg+self.median_bandpass_std,
               self.freqs, color='blue')
        plotxy(self.median_bandpass_avg-self.median_bandpass_std,
               self.freqs, color='blue')
        plotxy(self.bandpass_avg[self.zap_chans],
               self.freqs[self.zap_chans],
               line=None, symbol=16, color="red")
        closeplot()

    def plot_zapped_bandpass(self, device="/xwin"):
        not_zapped = set(np.arange(self.nchan)) - set(self.zap_chans)
        not_zapped = np.asarray(list(not_zapped))
        if len(not_zapped):
            yhi = (self.median_bandpass_avg+1.5*self.median_bandpass_std).max()
            ylo = (self.median_bandpass_avg-1.5*self.median_bandpass_std).min()
            plotxy(self.median_bandpass_avg, self.freqs, rangey=[ylo, yhi],
                   labx="Frequency (MHz)", color='light gray', device=device)
            plotxy(self.median_bandpass_avg+self.median_bandpass_std,
                   self.freqs, color='blue')
            plotxy(self.median_bandpass_avg-self.median_bandpass_std,
                   self.freqs, color='blue')
            plotxy(self.bandpass_avg[not_zapped], self.freqs[not_zapped], color='white')
            plotxy(self.median_bandpass_avg[self.zap_chans], self.freqs[self.zap_chans],
                   line=None, symbol=16, color='red')
            closeplot()
        else:
            print("WARNING!:  All channels recommended for masking!")

    def write_zap_chans(self, filename=None):
        if filename is None:
            filename = self.basename+".zapchans"
        outfile = open(filename, "w")
        for chan in np.unique(sorted(self.zap_chans)):
            outfile.write("%d\n" % chan)
        outfile.close()

    def set_weights_and_offsets(self):
        # The channel with the highest, non-zapped std has a weight of 1.0
        # That channels std is called std_norm.  All the others are set to
        # std_norm / std[i], where i is the channel number
        not_zapped = set(np.arange(self.nchan)) - set(self.zap_chans)
        not_zapped = np.asarray(list(not_zapped))
        if len(not_zapped):
            std_norm = self.bandpass_std[not_zapped].max()
            has_var = self.bandpass_std != 0.0
            # weights for channels without variance will automatically be 0
            self.weights = np.zeros_like(self.bandpass_std)
            self.weights[has_var] = std_norm / self.bandpass_std[has_var]
            self.weights[self.zap_chans] = 0.0
            self.offsets = self.bandpass_avg
        else:
            print("WARNING!:  All channels recommended for masking!")

    def write_weights_and_offsets(self, filename=None, invertband=False):
        if filename is None:
            filename = self.basename+".weights"
        outfile = open(filename, "w")
        outfile.write("# Chan    Weight    Offset\n")
        weights = self.weights[::-1] if invertband else self.weights
        offsets = self.offsets[::-1] if invertband else self.offsets
        for c, w, o in zip(np.arange(self.nchan), self.weights, self.offsets):
            outfile.write("%5d     %7.5f   %7.5f\n" % (c, w, o))
        outfile.close()

    def write_bandpass(self, filename=None, invertband=False):
        if filename is None:
            filename = self.basename+".bandpass"
        outfile = open(filename, "w")
        avg = self.bandpass_avg[::-1] if invertband else self.bandpass_avg
        std = self.bandpass_std[::-1] if invertband else self.bandpass_std
        freqs = self.freqs[::-1] if invertband else self.freqs
        outfile.write("# Chan   Freq(MHz)     Mean       StDev\n")
        for ii, (freq, mean, stdev) in enumerate(zip(freqs, avg, std)):
            outfile.write("%6d  %9.3f  %9.3f  %9.3f\n" % (ii, freq, mean, stdev))
        outfile.close()

    def write_median_bandpass(self, filename=None, invertband=False):
        if filename is None:
            filename = self.basename+".median_bandpass"
        outfile = open(filename, "w")
        avg = self.median_bandpass_avg[::-1] if invertband else self.median_bandpass_avg
        std = self.median_bandpass_std[::-1] if invertband else self.median_bandpass_std
        freqs = self.freqs[::-1] if invertband else self.freqs
        outfile.write("# Chan   Freq(MHz)     Mean       StDev\n")
        for ii, (freq, mean, stdev) in enumerate(zip(freqs, avg, std)):
            outfile.write("%6d  %9.3f  %9.3f  %9.3f\n" % (ii, freq, mean, stdev))
        outfile.close()

    def write_weights(self, threshold=0.05, filename=None, invertband=False):
        # This weights file works for psrfits_subband
        if filename is None:
            filename = self.basename+".weights"
        outfile = open(filename, "w")
        outfile.write("# Chan    Weight\n")
        weights = self.weights[::-1] if invertband else self.weights
        for c, w in zip(np.arange(self.nchan), weights):
            if w > threshold:
                outfile.write("%5d     1\n" % (c))
            else:
                outfile.write("%5d     0\n" % (c))
        outfile.close()

if __name__=="__main__":
    import sys
    a = rfifind(sys.argv[1])
    if len(a.goodints)==0:
        exit(0)
    sys.stderr.write("WARNING!:  If raw data have channels in decreasing freq\n")
    sys.stderr.write("           order, the channel ordering as given will be\n")
    sys.stderr.write("           inverted!  Use 'invertband=True' in \n")
    sys.stderr.write("           write_weights() in that case!\n")
    if (a.idata.telescope=='GBT' and a.idata.lofreq < 1000.0):
        sys.stderr.write("Data is from GBT Prime Focus, auto-flipping the weights/offsets...\n\n")
        invert = True
    else:
        invert = False
    a.set_zap_chans(power=200.0,
                    edges=0.01,
                    asigma=2.0,
                    ssigma=2.0,
                    usemask=True,
                    plot=True,
                    chans=[])
    a.write_zap_chans()
    a.set_weights_and_offsets()
    a.write_weights(invertband=invert)
    a.write_bandpass(invertband=invert)
    #a.write_weights_and_offsets(invertband=invert)
