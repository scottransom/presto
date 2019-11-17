#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip, str, range, object
from operator import attrgetter
import sys, re, os, copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path
import glob
from presto import infodata
from presto.presto import candidate_sigma

# Note: the following are global variables that can
#       (and should) be set in whatever module
#       imports this file.  The following values are
#       "OK" defaults for some searches....

# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
r_err = 1.1
# Longest period candidates to consider (s)
long_period = 15.0
# Shortest period candidates to consider (s)
short_period = 0.0005
# Ignore candidates with a sigma (from incoherent power summation) less than this
sigma_threshold = 6.0
# Ignore candidates with a coherent power less than this
c_pow_threshold = 100.0
# Ignore any candidates where at least one harmonic does not exceed this power
harm_pow_cutoff = 8.0

# If the birds file works well, the following shouldn't
# be needed at all...
#                (ms, err)
known_birds_p = []
#                (Hz, err)
known_birds_f = []
#---------------------------------------------------

fund_re = re.compile("^\d")
harms_re = re.compile("^[ ]\d")
DM_re = re.compile("DM(\d+\.\d{2})")

# Add some functions to maintain support for the old
# sifting API
def remove_duplicate_candidates(candlist, *args, **kwargs):
    """Remove duplicate candidates. The candlist is modified
        **in-place**.

        Note: This function is defined to maintain support
            for old code. It simply calls the 
            'remove_duplicate_candidates' method of candlist.

        Inputs:
            ** All arguments are passed onto the
            'remove_duplicate_candidates' method of candlist.

        Output:
            candlist: The modified candidate list.
    """
    candlist.remove_duplicate_candidates(*args, **kwargs)
    return candlist


def remove_DM_problems(candlist, *args, **kwargs):
    """Remove candidates with DM problems. The candlist is modified
        **in-place**.

        Note: This function is defined to maintain support
            for old code. It simply calls the 
            'remove_DM_problems' method of candlist.

        Inputs:
            ** All arguments are passed onto the
            'remove_DM_problems' method of candlist.

        Output:
            candlist: The modified candidate list.
    """
    candlist.remove_DM_problems(*args, **kwargs)
    return candlist


def remove_harmonics(candlist, *args, **kwargs):
    """Remove harmonics. The candlist is modified
        **in-place**.

        Note: This function is defined to maintain support
            for old code. It simply calls the 
            'remove_harmonics' method of candlist.

        Inputs:
            ** All arguments are passed onto the
            'remove_harmonics' method of candlist.

        Output:
            candlist: The modified candidate list.
    """
    candlist.remove_harmonics(*args, **kwargs)
    return candlist


def write_candlist(candlist, *args, **kwargs):
    candlist.to_file(*args, **kwargs)


def sigma_to_size(sigmas):
    """Given a numpy array of sigma values, return an array
        of same size with sizes of markers to plot.

        Inputs:
            sigmas: Numpy array of sigma values.

        Output:
            sizes: Numpy array of marker sizes.
    """
    # return 8+sigmas**1.7
    return np.clip(20**(sigmas/6), 5, 400)


def print_sift_globals():
    print("r_err =", r_err)
    print("short_period =", short_period)     
    print("long_period =", long_period)     
    print("sigma_threshold =", sigma_threshold) 
    print("c_pow_threshold =", c_pow_threshold) 
    print("harm_pow_cutoff =", harm_pow_cutoff) 
    print("known_birds_p =", known_birds_p)   
    print("known_birds_f =", known_birds_f)


def parse_power(pow):
    power = float(pow.split("(")[0])
    if ("^" in pow):  # add exponent...
        try:
            expon = float(pow.split("^")[1])
        except ValueError:
            expon = 5 # power gets chopped off if this large
        power *= 10.0**(expon)
    return power


class Candidate(object):
    def __init__(self, candnum, sigma, numharm, ipow, cpow, bin, z,
                 DMstr, filename, T):
        self.path, self.filename = os.path.split(filename)
        self.candnum = candnum
        self.sigma = sigma
        self.numharm = numharm
        self.ipow_det = ipow
        self.cpow = cpow
        self.r = bin
        self.f = bin/T
        self.z = z
        self.T = T
        self.p = 1.0/self.f
        self.DMstr = DMstr
        self.DM = float(DMstr)
        self.harm_pows = None
        self.harm_amps = None
        self.snr = 0.0
        self.hits = []
        self.note = ""

    def add_as_hit(self, other):
        self.hits.extend(other.hits)

    def __str__(self):
        cand = self.filename + ':' + repr(self.candnum)
        return "%-65s   %7.2f  %6.2f  %6.2f  %s   %7.1f  %7.1f  " \
                "%12.6f  %10.2f  %8.2f " % \
                        (cand, self.DM, self.snr, self.sigma, \
                        "%2d".center(7)%self.numharm,
                        self.ipow_det, self.cpow, self.p*1000, self.r, self.z)

    def harms_to_snr(self):
        # Remove the average power level
        harmamps = np.asarray(self.harm_pows) - 1.0
        # Set the S/N to 0.0 for harmonics with "negative" amplitudes
        harmamps[harmamps < 0.0] = 0.0
        self.snr = np.sum(np.sqrt(harmamps))


class Candlist(object):
    def __init__(self, cands=None, trackbad=False, trackdupes=False):
        if cands is None:
            self.cands = []
        else:
            self.cands = cands
        self.trackbad = trackbad # Should we keep track of bad candidates
        self.trackdupes = trackdupes # Should we keep track of duplicates
        # Set default badlists
        self.badlists = {'knownbirds': [], \
                         'longperiod': [], \
                         'shortperiod': [], \
                         'threshold': [], \
                         'harmpowcutoff': [], \
                         'rogueharmpow': [], \
                         'harmonic': [], \
                         'dmproblem': []}
        self.duplicates = []

    def __iter__(self):
        return iter(self.cands)

    def __getitem__(self, key):
        return self.cands[key]

    def __delitem__(self, key):
        del(self.cands[key])

    def sort(self, *args, **kwargs):
        self.cands.sort(*args, **kwargs)

    def plot_summary(self, usefreqs=True):
        """Produce a plot summarizing the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Get all candidates and sort by sigma
        allcands = self.get_all_cands()
        sigmas = np.array([c.sigma for c in allcands])
        isort = sigmas.argsort()
        sigmas = sigmas[isort]

        if usefreqs:
            xdata = np.array([c.f for c in allcands])[isort]
            xlabel = "Freq (Hz)"
            xscale = "log"
        else:
            xdata = np.array([c.p for c in allcands])[isort]
            xlabel = "Period (s)"
            xscale = "loglin"

        dms = np.array([c.DM for c in allcands])[isort]
        numharms = np.array([c.numharm for c in allcands])[isort]

        # Plot the all candidates 
        scatt = plt.scatter(xdata, dms, s=sigma_to_size(sigmas), \
                                c=np.log2(numharms), \
                                marker='o', alpha=0.7, zorder=-1) 
        plt.set_cmap("Spectral") 
  
        # Add colorbar 
        fmtr = matplotlib.ticker.FuncFormatter(lambda x, pos: "%d" % 2**x)
        cax = plt.axes((0.18, 0.06, 0.67, 0.035))
        cb = plt.colorbar(scatt, cax=cax, ticks=(0,1,2,3,4), format=fmtr, \
                            orientation="horizontal")
        cb.set_label("Num harmonics summed") 
        
        plt.axes(ax) # Set scatter plot's axes as current
        plt.xscale(xscale)
        plt.xlabel(xlabel)
        mindm = np.min(dms)
        maxdm = np.max(dms)
        dmrange = np.ptp(dms)

        # Use log-scale y-axis if max DM > 2000
        yscale = "log" if maxdm > 2000.0 else "linear"
        plt.yscale(yscale)

        if yscale is "log":
            plt.ylim(1.0, maxdm+0.1*dmrange)
        else:
            plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)

        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        if not usefreqs:
            plt.gca().xaxis.set_ticks(np.concatenate((\
                                        np.logspace(-4,0,4, endpoint=False), \
                                        np.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(np.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(xdata)/5.0), \
                        min(long_period+0.5, max(xdata)+0.5))
        else:
            plt.xlim(min(xdata)/5.0, max(xdata)*2.0)
            
        ax.format_coord = lambda x,y: "x=%g, y=%g" % (x,y)
        return fig

    def plot_rejects(self, usefreqs=True):
        """Produce a plot showing why candidates were rejected by
            the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Plot bad candidates
        candlists = [self.badlists['knownbirds'], self.badlists['longperiod'], \
                     self.badlists['shortperiod'], self.badlists['threshold'], \
                     self.badlists['harmpowcutoff'], self.badlists['rogueharmpow'], \
                     self.badlists['harmonic'], self.badlists['dmproblem'], \
                     self.cands, self.duplicates]
        labels = ['Known birdires', 'Long period', 'Short period', \
                    'Threshold', 'Harm power cutoff', 'Rogue harm power', \
                    'Harmonic cand', 'DM problem', 'Good cands', 'Hits']
        colours = ['#FF0000', '#800000', '#008000', '#00FF00', \
                    '#00FFFF', '#0000FF', '#FF00FF', '#800080', 'r', 'k']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'x', 's']
        zorders = [-2, -2, -2, -2, -2, -2, -2, -2, 0, 0]
        sizes = [50, 50, 50, 50, 50, 50, 50, 50, 100, 10]
        fixedsizes = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
        lws = [1,1,1,1,1,1,1,1,2,1]
        handles = []
        for cands, colour, marker, zorder, size, fixedsize, lw in \
                zip(candlists, colours, markers, zorders, sizes, fixedsizes, lws):
            if len(cands):
                sigmas = np.array([c.sigma for c in cands])
                isort = sigmas.argsort()
                sigmas = sigmas[isort]
                if usefreqs:
                    xdata = np.array([c.f for c in cands])[isort]
                    xlabel = "Freq (Hz)"
                    xscale = "log"
                else:
                    xdata = np.array([c.p for c in cands])[isort]
                    xlabel = "Period (s)"
                    xscale = "loglin"
                dms = np.array([c.DM for c in cands])[isort]
                
                # Plot the candidates
                if fixedsize:
                    plt.scatter(xdata, dms, s=size, lw=lw, \
                                c=colour, marker=marker, alpha=0.7, zorder=zorder)
                else:
                    plt.scatter(xdata, dms, s=sigma_to_size(sigmas), lw=lw, \
                                c=colour, marker=marker, alpha=0.7, zorder=zorder)
            handles.append(plt.scatter([0], [0], s=size, c=colour, \
                                    marker=marker, alpha=0.7))

        fig.legend(handles, labels, 'lower center', \
                        prop={'size':'x-small'}, ncol=4)

        plt.xscale(xscale) 
        plt.xlabel(xlabel)

        alldms = np.array([c.DM for c in self.get_all_cands()])
        mindm = np.min(alldms)
        maxdm = np.max(alldms)
        dmrange = np.ptp(alldms)

        # Use log-scale y-axis if max DM > 2000
        yscale = "log" if maxdm > 2000.0 else "linear"
        plt.yscale(yscale)

        if yscale is "log":
            plt.ylim(1.0, maxdm+0.1*dmrange)
        else:
            plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)

        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        if not usefreqs:
            all_xdata = np.array([c.p for c in self.get_all_cands()])
            plt.gca().xaxis.set_ticks(np.concatenate((\
                                        np.logspace(-4,0,4, endpoint=False), \
                                        np.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(np.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(all_xdata)/5.0), \
                        min(long_period+0.5, max(all_xdata)+0.5))
        else:
            all_xdata = np.array([c.f for c in self.get_all_cands()])
            plt.xlim(min(all_xdata)/5.0, max(all_xdata)*2.0)

        return fig

    def plot_goodcands(self, usefreqs=True):
        """Produce a plot highlighting good candidates as selected by
            the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Plot candidates
        labels = []
        candlists = []
        for key in self.badlists:
            labels.append(key.title())
            candlists.append(self.badlists[key])
        candlists.append(self.cands)
        labels.append('Good cands')
        colours = ['#FF0000', '#800000', '#008000', '#00FF00', \
                    '#00FFFF', '#0000FF', '#FF00FF', '#800080', 'r']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
        zorders = [-2, -2, -2, -2, -2, -2, -2, -2, 0]
        sizes = [10, 10, 10, 10, 10, 10, 10, 10, 50]
        fixedsizes = [1, 1, 1, 1, 1, 1, 1, 1, 0]
        lws = [1,1,1,1,1,1,1,1,1,1]
        ecs = ['none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'k']
        alphas = [1,1,1,1,1,1,1,1,0.7]
        handles = []
        for cands, colour, marker, zorder, size, fixedsize, lw, alpha, ec in \
                zip(candlists, colours, markers, zorders, sizes, fixedsizes, lws, alphas, ecs):
            sigmas = []
            dms = []
            xdata = []
            for c in cands:
                sigmas.extend([h.sigma for h in c.hits])
                dms.extend([h[0] for h in c.hits])
                if usefreqs:
                    xval = c.f
                else:
                    xval = c.p
                xdata.extend([xval]*len(c.hits))
            sigmas = np.array(sigmas)
            dms = np.array(dms)
            xdata = np.array(xdata)

            isort = sigmas.argsort()
            sigmas = sigmas[isort]
            dms = dms[isort]
            xdata = xdata[isort]
            if usefreqs:
                xlabel = "Freq (Hz)"
                xscale = "log"
            else:
                xlabel = "Period (s)"
                xscale = "loglin"
            
            # Plot the candidates
            if fixedsize:
                plt.scatter(xdata, dms, s=size, lw=lw, edgecolors=ec, \
                            c=colour, marker=marker, alpha=alpha, zorder=zorder)
            else:
                plt.scatter(xdata, dms, s=sigma_to_size(sigmas), lw=lw, edgecolors=ec, \
                            c=colour, marker=marker, alpha=alpha, zorder=zorder)
            handles.append(plt.scatter([], [], s=size, c=colour, \
                                    marker=marker, alpha=0.7))

        fig.legend(handles, labels, 'lower center', \
                        prop={'size':'x-small'}, ncol=4)

        plt.xscale(xscale) 
        plt.xlabel(xlabel) 
        mindm = np.min(dms)
        maxdm = np.max(dms)
        dmrange = np.ptp(dms)
        plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)
        plt.ylabel(r"DM (pc cm$^{-3}$)")
        if not usefreqs:
            plt.gca().xaxis.set_ticks(np.concatenate((\
                                        np.logspace(-4,0,4, endpoint=False), \
                                        np.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(np.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(xdata)/5.0), \
                        min(long_period+0.5, max(xdata)+0.5))
        return fig

    def mark_as_bad(self, icand, badlistname):
        cand = self.cands.pop(icand)
        if self.trackbad:
            badlist = self.badlists.setdefault(badlistname, [])
            badlist.append(cand)

    def mark_as_duplicate(self, icand):
        cand = self.cands.pop(icand)
        if self.trackdupes:
            self.duplicates.append(cand)

    def get_all_cands(self):
        cands = self.get_all_goodcands()
        return self.get_all_goodcands() + self.get_all_badcands()

    def get_all_goodcands(self):
        return self.cands + self.duplicates

    def get_all_badcands(self):
        cands = []
        for key in list(self.badlists.keys()):
            cands += self.badlists[key]
        return cands

    def reject_longperiod(self, long_period=None):
        """Find and remove very long period candidates.
            Long period candidates are moved to the 'longperiod' bad-list.

            Inputs:
                long_period: The longest allowed period for a 'good' cand.
                    (Default: Globally defined limit, "long_period")

            Outputs:
                None
        """
        if long_period is None:
            long_period = globals()['long_period']
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            if (cand.p > long_period):
                cand.note = "Period is too long (%g ms > %g ms)" % \
                            (cand.p*1000, long_period*1000)
                self.mark_as_bad(ii, 'longperiod')
    
    def reject_shortperiod(self, short_period=None):
        """Find and remove very short period candidates.
            Short period candidates are moved to the 'shortperiod' bad-list.

            Inputs:
                short_period: The shortest allowed period for a 'good' cand.
                    (Default: Globally defined limit, "short_period")

            Outputs:
                None
        """
        if short_period is None:
            short_period = globals()['short_period']
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            if (cand.p < short_period):
                cand.note = "Period is too short (%g ms < %g ms)" % \
                            (cand.p*1000, short_period*1000)
                self.mark_as_bad(ii, 'shortperiod')

    def reject_knownbirds(self, known_birds_f=[], known_birds_p=[]):
        """Find and remove candidates conincident with known birds.

            Inputs:
                known_birds_f: A list of tuples containing bad frequencies
                    and widths. The tuples should contain
                        (<bad freq (Hz)>, <one-sided width (Hz)>)
                    (Default: Globally defined "known_birds_f")
                known_birds_p: A list of tuples containing bad peridocities
                    and widths. The tuples should contain
                        (<bad freq (ms)>, <one-sided width (ms)>)
                    (Default: Globally defined "known_birds_p")

            Outputs:
                None
        """
        if known_birds_f is None:
            known_birds_f = globals()['known_birds_f']
        if known_birds_p is None:
            known_birds_p = globals()['known_birds_p']
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            known_bird = 0
            for bird, err in known_birds_f:
                if (np.fabs(cand.f-bird) < err):
                    known_bird = 1
                    cand.note = "Freq (%.2f Hz) is within %g Hz " \
                                    "of a known birdie centred at %.2f Hz" % \
                                    (cand.f, err, bird)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue
            for bird, err in known_birds_p:
                if (np.fabs(cand.p*1000.0-bird) < err):
                    known_bird = 1
                    cand.note = "Period (%.2f ms) is within %g ms " \
                                    "of a known birdie centred at %.2f ms" % \
                                    (cand.f*1000, err, bird)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue

    def reject_threshold(self, sigma_threshold=None, \
                    c_pow_threshold=None):
        """Find and remove candidates that don't pass our threshold.
            The conditions for rejection are different for candidates
            with a single harmonic, and multiple harmonics.

            Inputs:
                sigma_threshold: The threshold for sigma.
                    (Default: Globally defined threshold, "sigma_threshold")
                c_pow_threshold: The threshold for coherent power.
                    (Default: Globally defined threshold, "c_pow_threshold")

            Outputs:
                None
        """
        if sigma_threshold is None:
            sigma_threshold = globals()['sigma_threshold']
        if c_pow_threshold is None:
            c_pow_threshold = globals()['c_pow_threshold']
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            
            if cand.numharm == 1:
                # Single harmonic case
                if (cand.sigma < sigma_threshold) and \
                                    (cand.cpow < c_pow_threshold):
                    cand.note = "Only 1 harmonic and both sigma " \
                                "(%g < %g) and coherent power (%g < %g) are " \
                                "too low." % (cand.sigma, sigma_threshold, \
                                                cand.cpow, c_pow_threshold)
                    self.mark_as_bad(ii, 'threshold')
            else:
                # Multiple harmonic case
                if cand.sigma < sigma_threshold:
                    cand.note = "%d harmonics and sigma " \
                                "(%g < %g) is too low." % \
                                (cand.numharm, cand.sigma, sigma_threshold)
                    self.mark_as_bad(ii, 'threshold')
        

    def reject_harmpowcutoff(self, harm_pow_cutoff=None):
        """Find and remove the candidates where the harmonic with the
            highest power is not more than harm_pow_cutoff.

            Inputs:
                harm_pow_cutoff: Minimum power for a good harmonic.
                    (Default: Globally defined "harm_pow_cutoff")

            Outputs:
                None
        """
        if harm_pow_cutoff is None:
            harm_pow_cutoff = globals()['harm_pow_cutoff']
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            maxharm = np.argmax(cand.harm_pows)
            maxpow = cand.harm_pows[maxharm]
            if maxpow < harm_pow_cutoff:
                cand.note = "All harmonics have power < %g" % harm_pow_cutoff
                self.mark_as_bad(ii, 'harmpowcutoff')

    def reject_rogueharmpow(self):
        """Find and remove candidates which are dominated by a single
            high-power but high-numbered harmonic.

            Inputs:
                None
                
            Ouputs:
                None
        """
        for ii in reversed(list(range(len(self.cands)))):
            cand = self.cands[ii]
            maxharm = np.argmax(cand.harm_pows)
            maxpow = cand.harm_pows[maxharm]
            
            # Sort the harmonics by power
            sortedpows = np.sort(cand.harm_pows)

            if (cand.numharm >= 8 and maxharm > 4 and \
                                        maxpow > 2*sortedpows[-2]):
                # Max-power harmonic is at least 2x more powerful 
                # than the next highest-power harmonic, and is the
                # 4+th harmonic our of 8+ harmonics
                cand.note = "High-numbered harmonic (%d) has too " \
                            "much power" % maxharm
                self.mark_as_bad(ii, 'rogueharmpow')
            elif (cand.numharm >= 4 and maxharm > 2 and \
                                        maxpow > 3*sortedpows[-2]):
                # Max-power harmonic is at least 3x more powerful 
                # than the next highest-power harmonic, and is the
                # 2+th harmonic our of 4+ harmonics
                cand.note = "High-numbered harmonic (%d) has too " \
                            "much power" % maxharm
                self.mark_as_bad(ii, 'rogueharmpow')
        
    def default_rejection(self):
        """Run all rejection methonds with default arguments.

            Inputs:
                None

            Outputs:
                None
        """
        self.reject_longperiod()
        self.reject_shortperiod()
        self.reject_knownbirds()
        self.reject_threshold()
        self.reject_harmpowcutoff()
        self.reject_rogueharmpow()

    def remove_duplicate_candidates(self, verbosity=1):
        """Remove lower-significance 'duplicate' (i.e. same period)
            candidates from a list of candidates.  For the highest
            significance candidate, include a list of the DMs (and SNRs)
            of all the other detections.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        if verbosity >= 1:
            print("  Sorting the %d candidates by frequency..." % \
                        self.get_numcands())
        self.cands.sort(key=attrgetter('r'))
        if verbosity >= 1:
            print("  Searching for dupes...")
        ii = 0
        # Find any match
        while ii < self.get_numcands():
            jj = ii + 1
            if jj < self.get_numcands() and \
                        np.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                # Find others that match
                jj += 1
                while jj < self.get_numcands() and \
                        np.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                    jj += 1
                matches = self.cands[ii:jj]
                matches.sort(key=attrgetter('sigma'), reverse=True)
                bestindex = self.cands.index(matches[0])
                #sigmas = [c.sigma for c in matches]
                #bestindex = np.argmax(sigmas)+ii
                # flag the duplicates
                bestcand = self.cands[bestindex]
                # Add other matching cands as hit of highest-sigma cand
                for matchind in reversed(list(range(ii, jj))):
                    if matchind == bestindex:
                        # The current candidate is the highest-sigma cand
                        # Don't remove it
                        continue
                    match = self.cands[matchind]
                    bestcand.add_as_hit(match)
                    match.note = "This candidate is a duplicate of %s:%d" % \
                                (bestcand.filename, bestcand.candnum)
                    self.mark_as_duplicate(matchind)
                    if verbosity >= 2:
                        print("Removing %s:%d (index: %d)" % \
                                (match.filename, match.candnum, matchind))
                        print("    %s" % match.note)
                # If the best candidate isn't at the same freq
                # as ii, then it's possible even more hits should
                # be added. So we don't increment the index
                # (note that the best cand has moved into position ii).
            else:
                ii += 1 # No candidates to be added as hits, move on
        if verbosity >= 1:
            print("Found %d candidates.\n" % self.get_numcands())
        self.cands.sort(key=attrgetter('sigma'), reverse=True)

    def remove_harmonics(self, verbosity=1):
        """Remove the candidates that are lower significance harmonics
            of other candidates from the candlist.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Note:  should probably put the harmonics into the fundamental as hits (use sets)
        numremoved = 0
        self.cands.sort(key=attrgetter('sigma'), reverse=True)
        f_err = r_err/self.cands[0].T
        if verbosity >= 1:
            print("\nSearching for duplicate harmonics...")
        ii = 0
        while 1:
            fundcand = self.cands[ii]
            jj = len(self.cands) - 1
            zapj = 0
            while 1:
                harmcand = self.cands[jj]
                if zapj:  print("Hey!")
                for factor in np.arange(1.0, 17.0):
                    if np.fabs(fundcand.f - harmcand.f*factor) < f_err*factor:
                        zapj = 1
                        harmstr = "1/%dth" % factor
                    elif np.fabs(fundcand.f - harmcand.f/factor) < f_err/factor:
                        zapj = 1
                        if factor==2.0:
                            harmstr = "%dnd" % factor
                        else:
                            harmstr = "%dth" % factor
                    if zapj:
                        if verbosity >= 2:
                            print("Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%s) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        harmstr, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f))
                        break
                # Check a few other common ratios
                for numer,denom in zip([3.0, 5.0, 2.0, 4.0, 5.0, \
                                        3.0, 5.0, 2.0, 3.0, 4.0],
                                       [2.0, 2.0, 3.0, 3.0, 3.0, \
                                        4.0, 4.0, 5.0, 5.0, 5.0]):
                    factor = numer/denom
                    if np.fabs(fundcand.f-harmcand.f*factor) < f_err*factor:
                        if verbosity >= 2:
                            print("Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%d/%dth) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        denom, \
                                        numer, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f))
                        harmstr = "%d/%dth" % (denom, numer)
                        zapj = 1
                        break
                if zapj:
                    harmcand.note = "This candidate (P=%.4f s, DM=%.2f) is " \
                                    "a harmonic (%s) of %s:%d " \
                                    "(P=%.4f s, DM=%.2f)." % \
                                (harmcand.p, harmcand.DM, harmstr, \
                                    fundcand.filename, fundcand.candnum, \
                                    fundcand.p, fundcand.DM)
                    numremoved += 1
                    self.mark_as_bad(jj, 'harmonic')
                    if verbosity >= 2:
                        print("Removing %s:%d (index: %d)" % \
                                (harmcand.filename, harmcand.candnum, jj))
                        print("    %s" % harmcand.note)
                    zapj = 0
                jj -= 1
                if jj == ii:
                    break
            ii += 1
            if ii >= len(self.cands) - 1:
                break
        if verbosity >= 1:
            print("Removed a total of %d harmonics.\n" % numremoved)

    def remove_DM_problems(self, numdms, dmlist, low_DM_cutoff, verbosity=1):
        """Remove the candidates where any of the following are true:
            1) The number of hits is < numdms
            2) The highest S/N candidate occurs below a DM of low_DM_cutoff
            3) The minimum difference in DM indices between the hits is > 1

            Inputs:
                numdms: The minimum number of hits for a good candidate.
                dmlist: List of DMs.
                low_DM_cutoff: The lowest DM possible for a good candidate.
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Create a dictionary where the key is the dmstr 
        # and the values are the index
        dmdict = {}
        dms = np.unique([float(dm) for dm in dmlist])
        dmstrs = ['%.2f'%dm for dm in dms]
        dmdict = dict(list(zip(dmstrs, list(range(len(dms))))))
        numremoved = 0
        num_toofew = 0
        num_toolow = 0
        num_gaps = 0
        self.cands.sort(key=attrgetter('sigma'), reverse=True)
        for ii in reversed(list(range(len(self.cands)))):
            currcand = self.cands[ii]
            # Remove all the candidates without enough DM hits
            if len(currcand.hits) < numdms:
                numremoved += 1
                num_toofew += 1
                currcand.note = "Candidate has only %d DM hits. This is less " \
                                "than minimum for 'good' cands (%d hits)" % \
                                (len(currcand.hits), numdms)
                self.mark_as_bad(ii, 'dmproblem')
                if verbosity >= 2:
                    print("Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii))
                    print("    %s" % currcand.note)
                continue

            # Remove all the candidates where the max sigma DM is 
            # less than the cutoff DM
            # Recall - A hit is a 3-tuple: (DM, SNR, sigma)
            imax = np.argmax(np.array([hit[2] for hit in currcand.hits]))
            hitdm, hitsnr, hitsigma = currcand.hits[imax]
            if float(hitdm) <= low_DM_cutoff:
                numremoved += 1
                num_toolow += 1
                currcand.note = "Hit with max sigma (%g) has dm (%.2f) " \
                                "<= low DM cutoff (%.2f) " % \
                                    (hitsigma, hitdm, low_DM_cutoff)
                self.mark_as_bad(ii, 'dmproblem')
                if verbosity >= 2:
                    print("Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii))
                    print("    %s" % currcand.note)
                continue

            # Remove all the candidates where there are no hits at consecutive DMs
            if len(currcand.hits) > 1:
                currcand.hits.sort(key=lambda cand: float(cand[0]))
                dm_indices = np.asarray([dmdict["%.2f"%currcand.hits[jj][0]]
                                          for jj in range(len(currcand.hits))])
                min_dmind_diff = min(dm_indices[1:] - dm_indices[:-1])
                if min_dmind_diff > 1:
                    numremoved += 1
                    num_gaps += 1
                    currcand.note = "DM list of hits has gaps (i.e. " \
                                    "consecutive DMs don't have hits)."
                    self.mark_as_bad(ii, 'dmproblem')
                    if verbosity >= 2:
                        print("Removing %s:%d (index: %d)" % \
                                (currcand.filename, currcand.candnum, ii))
                        print("    %s" % currcand.note)
                    continue

        if verbosity >= 1:
            print("Removed %d candidates with DM problems.\n" % numremoved)
        if verbosity >= 2:
            print("  # with too few hits:", num_toofew)
            print("  # with peak SNR too low:", num_toolow)
            print("  # with gaps in DM hits:", num_gaps)

    def print_cand_summary(self, summaryfilenm=None):
        """Write a summary of all candidates to file (or stdout).

            Input:
                summaryfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if summaryfilenm is None:
            summaryfile = sys.stdout
        elif summaryfilenm in [sys.stdout, sys.stderr]:
            summaryfile = summaryfilenm
        else:
            summaryfile = open(summaryfilenm, "w")
        summaryfile.write("   Candlist contains %d 'good' candidates\n" % \
                            len(self.cands))
        summaryfile.write("      # Known RFI rejects:           %d\n" % \
              len(self.badlists['knownbirds']))
        summaryfile.write("      # Short period rejects:        %d\n" % \
              len(self.badlists['shortperiod']))
        summaryfile.write("      # Long period rejects:         %d\n" % \
              len(self.badlists['longperiod']))
        summaryfile.write("      # Missed threshold:            %d\n" % \
              len(self.badlists['threshold']))
        summaryfile.write("      # No good harmonics:           %d\n" % \
              len(self.badlists['harmpowcutoff']))
        summaryfile.write("      # One bad harmonic:            %d\n" % \
              len(self.badlists['rogueharmpow']))
        summaryfile.write("      # Duplicate candidates:        %d\n" % \
              len(self.duplicates))
        summaryfile.write("      # Harmonic candidates:         %d\n" % \
              len(self.badlists['harmonic']))
        summaryfile.write("      # Candidates with DM problems: %d\n" % \
              len(self.badlists['dmproblem']))
        if summaryfilenm not in [None, sys.stdout, sys.stderr]:
            summaryfile.close()
  
    def write_cand_report(self, reportfilenm=None):
        """Write a report of all bad candidates to file (or stdout).

            Input:
                reportfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if reportfilenm is None:
            reportfile = sys.stdout
        else:
            reportfile = open(reportfilenm, "w")
        reportfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                       "SNR".center(8) + "sigma".center(8) + "numharm".center(9) +
                       "ipow".center(9) + "cpow".center(9) +  "P(ms)".center(14) +
                       "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
        badcands = self.get_all_badcands()
        for badcand in badcands:
            reportfile.write("%s (%d)\n" % (str(badcand), len(badcand.hits)))
            reportfile.write("    Note: %s\n\n" % badcand.note)
        if reportfilenm is not None:
            reportfile.close()
        

    def __add__(self, other):
        copy_of_self = copy.deepcopy(self)
        copy_of_self.extend(other)
        return copy_of_self

    def get_numcands(self):
        """Get the number of good candidates (i.e. len(self.cands)).

            Inputs:
                None

            Outputs:
                None
        """
        return len(self)

    def __len__(self):
        # return the number of good candidates
        return len(self.cands)

    def extend(self, other):
        """Extend Candlist with another. This combines
            the candidates, as well as the lists of bad cands.
        
            Inputs:
                other: A second Candlist object to extend from.

            Outputs:
                None - the original Candlist object is extended in place.
        """
        self.cands.extend(other.cands)
        self.duplicates.extend(other.duplicates)
        for key in other.badlists:
            bad = self.badlists.setdefault(key, [])
            bad.extend(other.badlists[key])

    def to_file(self, candfilenm=None):
        """Write Candlist to file (or stdout).
            
            Input:
                candfilenm: Name of file to write to. If None,
                    write to stdout. (Default: write to stdout).

            Outputs:
                None
        """
        if candfilenm is None:
            candfile = sys.stdout
        else:
            candfile = open(candfilenm, "w")
        candfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                       "SNR".center(8) + "sigma".center(8) + "numharm".center(9) +
                       "ipow".center(9) + "cpow".center(9) +  "P(ms)".center(14) +
                       "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
        for goodcand in self.cands:
            candfile.write("%s (%d)\n" % (str(goodcand), len(goodcand.hits)))
            if (len(goodcand.hits) > 1):
                goodcand.hits.sort(key=lambda cand: float(cand[0]))
                for hit in goodcand.hits:
                    numstars = int(hit[2]/3.0)
                    candfile.write("  DM=%6.2f SNR=%5.2f Sigma=%5.2f   "%hit + \
                                    numstars*'*' + '\n')
        if candfilenm is not None:
            candfile.close()


def candlist_from_candfile(filename, trackbad=False, trackdupes=False):
    candfile = open(filename, 'r')
    # First identify the length of the observation searched
    for line in candfile:
        if line.startswith(" Number of bins in the time series"):
            numsamp = int(line.split()[-1])
        if line.startswith(" Width of each time series bin (sec)"):
            dt = float(line.split()[-1])
    tobs = numsamp * dt

    # Go back to the start of the file to read the candidates
    candfile.seek(0)

    cands = []
    candnums = []
    current_goodcandnum = 0
    last_candnum = 0
    last_goodcandnum = 0

    for line in candfile:
        # Identify the candidates in the top of the file
        if fund_re.match(line):
            split_line = line.split()
            candnum   = int(split_line[0])
            if len(cands) and len(split_line[0])==4 and cands[-1].candnum >= 9999:
                candnum = last_candnum + 1
            sigma     = float(split_line[1])
            i_pow_det = float(split_line[2])
            c_pow     = float(split_line[3])
            numharm   = int(split_line[4])
            bin       = float(split_line[7].split("(")[0])
            z         = float(split_line[9].split("(")[0])
            f = bin / tobs    # Spin freq in hz
            p = 1.0 / f       # Spin period in sec

            # Add it to the candidates list
            DMstr = DM_re.search(filename).groups()[0]
            cands.append(Candidate(candnum, sigma, numharm,
                                          i_pow_det, c_pow, bin, z, 
                                          DMstr, filename, tobs))
            candnums.append(candnum)
            last_candnum = candnum
            continue

        # Parse the harmonic powers
        elif harms_re.match(line):
            split_line = line.split()
            candnum = int(split_line[0])
            # For rare cases where >10K cands, candnum is wrong in ACCEL file
            if len(split_line[0])==4 and last_goodcandnum >= 9999:
                candnum = last_goodcandnum + 1
            if candnum in candnums:
                cand = cands[candnums.index(candnum)]
                cand.harm_pows = np.zeros(cand.numharm, dtype=np.float64)
                cand.harm_amps = np.zeros(cand.numharm, dtype=np.complex64) 
                power = parse_power(split_line[3])
                phase = float(split_line[9].split("(")[0])
                cand.harm_pows[0] = power
                cand.harm_amps[0] = np.sqrt(power) * np.exp(phase*1.0j)
                if (cand.numharm > 1):
                    current_goodcandnum = candnum
                    current_harmnum = 1
                else:
                    last_goodcandnum = candnum
                    current_goodcandnum = 0
                    # Compute the S/N
                    cand.harms_to_snr()
                    # These are the "optimized" power...
                    opt_ipow = cand.harm_pows[0]
                    # and sigma (calculated assuming _1_ trial!)
                    opt_sigma = candidate_sigma(opt_ipow, 1, 1)
                    cand.sigma = opt_sigma
                    # Now that S/N and sigma are available
                    # List candidate as a hit of itself
                    cand.hits = [(cand.DM, cand.snr, cand.sigma)]
                    cand.ipow_det = opt_ipow
            continue

        # Parse the higher (than the first) harmonic powers
        if current_goodcandnum:
            cand = cands[candnums.index(current_goodcandnum)]
            power = parse_power(line.split()[2])
            phase = float(line.split()[8].split("(")[0])
            cand.harm_pows[current_harmnum] = power
            cand.harm_amps[current_harmnum] = np.sqrt(power) * np.exp(phase*1.0j)
            current_harmnum += 1
            # Calculate other stats after all harmonics have been read in
            if (current_harmnum==cand.numharm):
                # Compute the S/N
                cand.harms_to_snr()
                # Compute sigma and incoherent power
                opt_ipow = sum(cand.harm_pows)
                opt_sigma = candidate_sigma(opt_ipow, cand.numharm, 1)
                cand.sigma = opt_sigma
                # Now that S/N and sigma are available
                # List candidate as a hit of itself
                cand.hits = [(cand.DM, cand.snr, cand.sigma)]
                cand.ipow_det = opt_ipow
                last_goodcandnum = candnum
                current_goodcandnum = 0
    candfile.close()
    return Candlist(cands, trackbad=trackbad, trackdupes=trackdupes)


def read_candidates(filenms, prelim_reject=True, track=False):
    """Read in accelsearch candidates from the test ACCEL files.
        Return a Candlist object of Candidate instances.

        Inputs:
            filenms: A list of files to read candidates from.
            prelim_reject: If True, perform preliminary rejection of
                candidates. (Default: True)
            track: If True, keep track of bad/duplicate candidates.
                (Default: False)

    """
    candlist = Candlist(trackbad=track, trackdupes=track)
    numfiles = len(filenms)
    if filenms:
        print("\nReading candidates from %d files...." % len(filenms))
        for ii, filenm in enumerate(filenms):
            curr_candlist = candlist_from_candfile(filenm, trackbad=track, trackdupes=track)
            if prelim_reject:
                curr_candlist.default_rejection()
            candlist.extend(curr_candlist)
            sys.stdout.write(" Read %d of %d files (%d cands)\r" % (ii+1, numfiles, len(candlist)))
            sys.stdout.flush()
        print("\nDone")
    else:
        print("Error:  There are no candidate files to read!")
    return candlist


def sift_directory(dir, outbasenm):
    """Sift candidates in given directory.
        
        Inputs:
            dir: The directory containing candidates to sift.
            outbasenm: The base name of the output files.

        Outputs:
            None
    """

  
    # TODO: Remove hard-coded values in this function
    #       replace with command line options.

    global sigma_threshold
    sigma_threshold = 5.0

    print_sift_globals()
    # Get list of DMs from *.inf files
    inffns = glob.glob(os.path.join(dir, '*.inf'))
    dmstrs = ['%.2f'%infodata.infodata(inffn).DM for inffn in inffns]

    # Read candidates found in low acceleration searching
    lo_accel_fns = glob.glob(os.path.join(dir, '*ACCEL_0'))
    lo_accel_cands = read_candidates(lo_accel_fns, track=True)
    print("Read %d candidates from %d files" % \
                (len(lo_accel_cands), len(lo_accel_fns)))
    print("%d candidates passed default rejection" % len(lo_accel_cands))
    if len(lo_accel_cands):
        lo_accel_cands.remove_duplicate_candidates()
    if len(lo_accel_cands):
        lo_accel_cands.remove_DM_problems(2, dmstrs, 2.0)
    lo_accel_cands.print_cand_summary()

    # Read candidates found in high acceleration searching
    hi_accel_fns = glob.glob(os.path.join(dir, '*ACCEL_50'))
    hi_accel_cands = read_candidates(hi_accel_fns, track=True)
    print("Read %d candidates from %d files" % \
                (len(hi_accel_cands), len(hi_accel_fns)))
    print("%d candidates passed default rejection" % len(hi_accel_cands))
    if len(hi_accel_cands):
        hi_accel_cands.remove_duplicate_candidates()
    if len(hi_accel_cands):
        hi_accel_cands.remove_DM_problems(2, dmstrs, 2.0)
    hi_accel_cands.print_cand_summary()

    all_accel_cands = lo_accel_cands + hi_accel_cands
    if len(all_accel_cands):
        all_accel_cands.remove_harmonics()
        # Note:  the candidates will be sorted in _sigma_ order, not _SNR_!
        all_accel_cands.cands.sort(key=attrgetter('sigma'), reverse=True)
        print("Found %d good candidates" % len(all_accel_cands))
        all_accel_cands.to_file(outbasenm+".accelcands")
    all_accel_cands.write_cand_report(outbasenm+".accelcands.report")
    all_accel_cands.print_cand_summary()
    all_accel_cands.plot_rejects(usefreqs=True)
    plt.savefig(outbasenm+".accelcands.rejects.png")
    #all_accel_cands.plot_goodcands(usefreqs=False)
    #plt.savefig(outbasenm+".accelcands.goodcands.png")
    all_accel_cands.plot_summary(usefreqs=True)
    plt.savefig(outbasenm+".accelcands.summary.png")
    plt.show()


def LogLinScaleFactory(b):
    class LogLinScale(matplotlib.scale.ScaleBase):
        name = 'loglin'

        def __init__(self, axis, **kwargs):
            matplotlib.scale.ScaleBase.__init__(self)
            self.thresh = kwargs.pop("thresh", 1e-5)
            if self.thresh <= 0.0:
                raise ValueError("thresh must be larger than 0")

        def get_transform(self):
            return self.LogLinTransform(self.thresh)

        def set_default_locators_and_formatters(self, axis):
            pass

        def limit_range_for_scale(self, vmin, vmax, minpos):
            return max(vmin, self.thresh), vmax

        class LogLinTransform(matplotlib.transforms.Transform):
            input_dims = 1
            output_dims = 1
            is_separable = True
            brk = b

            def __init__(self, thresh):
                matplotlib.transforms.Transform.__init__(self)
                self.thresh = thresh

            def transform(self, a):
                aa = np.ma.masked_where(a<self.thresh, a)
                if aa.mask.any():
                    aa[a<self.brk] = np.ma.log10(a[a<self.brk]) - \
                                        np.log10(self.brk)+self.brk
                else:
                    aa[a<self.brk] = np.log10(a[a<self.brk]) - \
                                        np.log10(self.brk)+self.brk
                return aa

            def inverted(self):
                return LogLinScale.InvertedLogLinTransform(self.thresh)

        class InvertedLogLinTransform(matplotlib.transforms.Transform):
            input_dims = 1
            output_dims = 1
            is_separable = True
            brk = b

            def __init__(self, thresh):
                matplotlib.transforms.Transform.__init__(self)
                self.thresh = thresh

            def transform(self, a):
                aa = a.copy()
                aa[a<self.brk] = np.ma.power(10, a[a<self.brk]-self.brk + \
                                                np.log10(self.brk))
                return aa

            def inverted(self):
                return LogLinScale.LogLinTransform(self.thresh)
    return LogLinScale

matplotlib.scale.register_scale(LogLinScaleFactory(1))

def main():
    # Sift candidates in PWD
    outbasenm = sys.argv[1]
    sift_directory(os.getcwd(), outbasenm)


if __name__ == '__main__':
    main()
