#!/usr/bin/env python
import sys
import re
import os
import copy

import numpy as Num
import matplotlib.pyplot as plt

from presto import candidate_sigma

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
    copy_of_candlist = copy.deepcopy(candlist)
    copy_of_candlist.remove_duplicate_candidates(*args, **kwargs)
    return copy_of_candlist


def remove_DM_problems(candlist, *args, **kwargs):
    copy_of_candlist = copy.deepcopy(candlist)
    copy_of_candlist.remove_DM_problems(*args, **kwargs)
    return copy_of_candlist


def remove_harmonics(candlist, *args, **kwargs):
    copy_of_candlist = copy.deepcopy(candlist)
    copy_of_candlist.remove_harmonics(*args, **kwargs)
    return copy_of_candlist


def write_candlist(candlist, *args, **kwargs):
    candlist.to_file(*args, **kwargs)


def print_sift_globals():
    print "r_err =", r_err
    print "short_period =", short_period     
    print "long_period =", long_period     
    print "sigma_threshold =", sigma_threshold 
    print "c_pow_threshold =", c_pow_threshold 
    print "harm_pow_cutoff =", harm_pow_cutoff 
    print "known_birds_p =", known_birds_p   
    print "known_birds_f =", known_birds_f


def parse_power(pow):
    power = float(pow.split("(")[0])
    if ("^" in pow):  # add exponent...
        try:
            expon = float(pow.split("^")[1])
        except ValueError:
            expon = 5 # power gets chopped off if this large
        power *= 10.0**(expon)
    return power


def cmp_sigma(self, other):
    retval = -cmp(self.sigma, other.sigma)
    if retval==0:
        return -cmp(self.ipow_det, other.ipow_det)
    else:
        return retval


def cmp_snr(self, other):
    retval = -cmp(self.snr, other.snr)
    if retval==0:
        return -cmp(self.ipow_det, other.ipow_det)
    else:
        return retval


def cmp_dms(self, other):
    return cmp(float(self[0]), float(other[0]))


def cmp_freq(self, other):
    return cmp(self.r, other.r)


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
        cand = self.filename + ':' + `self.candnum`
        return "%-65s   %7.2f  %6.2f  %6.2f  %s   %7.1f  %7.1f  " \
                "%12.6f  %10.2f  %8.2f " % \
                        (cand, self.DM, self.snr, self.sigma, \
                        "%2d".center(7)%self.numharm,
                        self.ipow_det, self.cpow, self.p*1000, self.r, self.z)

    def harms_to_snr(self):
        # Remove the average power level
        harmamps = Num.asarray(self.harm_pows) - 1.0
        # Set the S/N to 0.0 for harmonics with "negative" amplitudes
        harmamps[harmamps < 0.0] = 0.0
        self.snr = Num.sum(Num.sqrt(harmamps))


class Candlist(object):
    def __init__(self, cands=None):
        if cands is None:
            self.cands = []
        else:
            self.cands = cands
        self.badcands_knownbirds = []
        self.badcands_longperiod = []
        self.badcands_shortperiod = []
        self.badcands_threshold = []
        self.badcands_harmpowcutoff = []
        self.badcands_rogueharmpow = []
        self.duplicate_cands = []
        self.harmonic_cands = []
        self.dmproblem_cands = []

    def sort(self, *args, **kwargs):
        self.cands.sort(*args, **kwargs)

    def plot_summary(self):
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.08, 0.87, 0.80)) 
        plt.set_cmap("Spectral") 
        
        # Get all candidates and sort by sigma
        allcands = self.get_all_cands()
        sigmas = Num.array([c.sigma for c in allcands])
        isort = sigmas.argsort()
        sigmas = sigmas[isort]
        freqs = Num.array([c.f for c in allcands])[isort]
        dms = Num.array([c.DM for c in allcands])[isort]
        numharms = Num.array([c.numharm for c in allcands])[isort]

        # Plot the all candidates 
        plt.scatter(freqs, dms, s=8+sigmas**1.7, c=Num.log2(numharms), \
                                marker='o', alpha=0.7, zorder=-1) 
  
        # Add colorbar 
        fmtr = matplotlib.ticker.FuncFormatter(lambda x, pos: "%d" % 2**x) 
        cb = plt.colorbar(ticks=(0,1,2,3,4), format=fmtr) 
        cb.set_label("Num harmonics summed") 
         
        plt.xscale('log', base=10.0) 
        plt.xlim(0.5, 10000) 
        plt.ylim(-10, 1200) 
        plt.xlabel("Freq (Hz)") 
        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        return fig

    def plot_goodbad(self):
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Plot bad candidates
        candlists = [self.badcands_knownbirds, self.badcands_longperiod, \
                     self.badcands_shortperiod, self.badcands_threshold, \
                     self.badcands_harmpowcutoff, self.badcands_rogueharmpow, \
                     self.harmonic_cands, self.dmproblem_cands, self.cands, \
                     self.duplicate_cands]
        labels = ['Known birdires', 'Long period', 'Short period', \
                    'Threshold', 'Harm power cutoff', 'Rogue harm power', \
                    'Harmonic cand', 'DM problem', 'Good cands', 'Hits']
        colours = ['#FF0000', '#800000', '#008000', '#00FF00', \
                    '#00FFFF', '#0000FF', '#FF00FF', '#800080', 'r', 'k']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'x', ',']
        zorders = [-2, -2, -2, -2, -2, -2, -2, -2, 0, 0]
        sizes = [50, 50, 50, 50, 50, 50, 50, 50, 100, 10]
        fixedsizes = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
        lws = [1,1,1,1,1,1,1,1,2,1]
        handles = []
        for cands, colour, marker, zorder, size, fixedsize, lw in \
                zip(candlists, colours, markers, zorders, sizes, fixedsizes, lws):
            sigmas = Num.array([c.sigma for c in cands])
            isort = sigmas.argsort()
            sigmas = sigmas[isort]
            freqs = Num.array([c.f for c in cands])[isort]
            dms = Num.array([c.DM for c in cands])[isort]
            
            # Plot the candidates
            if fixedsize:
                plt.scatter(freqs, dms, s=size, lw=lw, \
                            c=colour, marker=marker, alpha=0.7, zorder=zorder)
            else:
                plt.scatter(freqs, dms, s=8+sigmas**1.7, lw=lw, \
                            c=colour, marker=marker, alpha=0.7, zorder=zorder)
            handles.append(plt.scatter([], [], s=size, c=colour, \
                                    marker=marker, alpha=0.7))

        fig.legend(handles, labels, 'lower center', \
                        prop={'size':'x-small'}, ncol=4)

        plt.xscale('log', base=10.0) 
        plt.xlim(0.5, 10000) 
        plt.ylim(-10, 1200) 
        plt.xlabel("Freq (Hz)") 
        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        return fig

    def get_all_cands(self):
        return self.cands + self.badcands_knownbirds + self.badcands_longperiod + \
                self.badcands_shortperiod + self.badcands_threshold + \
                self.badcands_harmpowcutoff + self.badcands_rogueharmpow + \
                self.duplicate_cands + self.harmonic_cands + self.dmproblem_cands

    def get_all_goodcands(self):
        return self.cands + self.duplicate_cands

    def get_all_badcands(self):
        return self.badcands_knownbirds + self.badcands_longperiod + \
                self.badcands_shortperiod + self.badcands_threshold + \
                self.badcands_harmpowcutoff + self.badcands_rogueharmpow + \
                self.harmonic_cands + self.dmproblem_cands

    def reject_longperiod(self, long_period=None):
        """Find and remove very long period candidates.
            Long period candidates are moved to the 'badcands_longperiod' list.

            Inputs:
                long_period: The longest allowed period for a 'good' cand.
                    (Default: Globally defined limit, "long_period")

            Outputs:
                None
        """
        if long_period is None:
            long_period = globals()['long_period']
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            if (cand.p > long_period):
                cand.note = "Period is too long (%g ms > %g ms)" % \
                            (cand.p*1000, long_period*1000)
                self.badcands_longperiod.append(self.cands.pop(ii))
    
    def reject_shortperiod(self, short_period=None):
        """Find and remove very short period candidates.
            Short period candidates are moved to the 'badcands_shortperiod' list.

            Inputs:
                short_period: The shortest allowed period for a 'good' cand.
                    (Default: Globally defined limit, "short_period")

            Outputs:
                None
        """
        if short_period is None:
            short_period = globals()['short_period']
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            if (cand.p < short_period):
                cand.note = "Period is too short (%g ms > %g ms)" % \
                            (cand.p*1000, short_period*1000)
                self.badcands_shortperiod.append(self.cands.pop(ii))

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
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            known_bird = 0
            for bird, err in known_birds_f:
                if (Num.fabs(cand.f-bird) < err):
                    known_bird = 1
                    cand.note = "Freq (%.2f Hz) is within %g Hz " \
                                    "of a known birdie centred at %.2f Hz" % \
                                    (cand.f, err, bird)
                    break
            if known_bird:
                self.badcands_knownbirds.append(self.cands.pop(ii))
                continue
            for bird, err in known_birds_p:
                if (Num.fabs(cand.p*1000.0-bird) < err):
                    known_bird = 1
                    cand.note = "Period (%.2f ms) is within %g ms " \
                                    "of a known birdie centred at %.2f ms" % \
                                    (cand.f*1000, err, bird)
                    break
            if known_bird:
                self.badcands_knownbirds.append(self.cands.pop(ii))
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
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            
            if cand.numharm == 1:
                # Single harmonic case
                if (cand.sigma < sigma_threshold) and \
                                    (cand.cpow < c_pow_threshold):
                    cand.note = "Only 1 harmonic and both sigma " \
                                "(%g < %g) and coherent power (%g < %g) are " \
                                "too low." % (cand.sigma, sigma_threshold, \
                                                cand.cpow, c_pow_threshold)
                    self.badcands_threshold.append(self.cands.pop(ii))
            else:
                # Multiple harmonic case
                if cand.sigma < sigma_threshold:
                    cand.note = "%d harmonics and sigma " \
                                "(%g < %g) is too low." % \
                                (cand.numharm, cand.sigma, sigma_threshold)
                    self.badcands_threshold.append(self.cands.pop(ii))
        

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
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            maxharm = Num.argmax(cand.harm_pows)
            maxpow = cand.harm_pows[maxharm]
            if maxpow < harm_pow_cutoff:
                cand.note = "All harmonics have power < %g" % harm_pow_cutoff
                self.badcands_harmpowcutoff.append(self.cands.pop(ii))

    def reject_rogueharmpow(self):
        """Find and remove candidates which are dominated by a single
            high-power but high-numbered harmonic.

            Inputs:
                None
                
            Ouputs:
                None
        """
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            maxharm = Num.argmax(cand.harm_pows)
            maxpow = cand.harm_pows[maxharm]
            
            # Sort the harmonics by power
            sortedpows = Num.sort(cand.harm_pows)

            if (cand.numharm >= 8 and maxharm > 4 and \
                                        maxpow > 2*sortedpows[-2]):
                # Max-power harmonic is at least 2x more powerful 
                # than the next highest-power harmonic, and is the
                # 4+th harmonic our of 8+ harmonics
                cand.note = "High-numbered harmonic (%d) has too " \
                            "much power" % maxharm
                self.badcands_rogueharmpow.append(self.cands.pop(ii))
            elif (cand.numharm >= 4 and maxharm > 2 and \
                                        maxpow > 3*sortedpows[-2]):
                # Max-power harmonic is at least 3x more powerful 
                # than the next highest-power harmonic, and is the
                # 2+th harmonic our of 4+ harmonics
                cand.note = "High-numbered harmonic (%d) has too " \
                            "much power" % maxharm
                self.badcands_rogueharmpow.append(self.cands.pop(ii))
        
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
            print "  Sorting the %d candidates by frequency..." % \
                        self.get_numcands()
        self.cands.sort(cmp_freq)
        if verbosity >= 1:
            print "  Searching for dupes..."
        ii = 0
        # Find any match
        while ii < self.get_numcands():
            jj = ii + 1
            if jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                # Find others that match
                jj += 1
                while jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                    jj += 1
                matches = self.cands[ii:jj]
                matches.sort(cmp_sigma)
                bestindex = self.cands.index(matches[0])
                #sigmas = [c.sigma for c in matches]
                #bestindex = Num.argmax(sigmas)+ii
                # flag the duplicates
                bestcand = self.cands[bestindex]
                # Add other matching cands as hit of highest-sigma cand
                for matchind in reversed(range(ii, jj)):
                    if matchind == bestindex:
                        # The current candidate is the highest-sigma cand
                        # Don't remove it
                        continue
                    match = self.cands[matchind]
                    bestcand.add_as_hit(match)
                    match.note = "This candidate is a duplicate of %s:%d" % \
                                (bestcand.filename, bestcand.candnum)
                    self.duplicate_cands.append(self.cands.pop(matchind))
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (match.filename, match.candnum, matchind)
                        print "    %s" % match.note
                # If the best candidate isn't at the same freq
                # as ii, then it's possible even more hits should
                # be added. So we don't increment the index
                # (note that the best cand has moved into position ii).
            else:
                ii += 1 # No candidates to be added as hits, move on
        if verbosity >= 1:
            print "Found %d candidates.\n" % self.get_numcands()
        self.cands.sort(cmp_sigma)

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
        self.cands.sort(cmp_sigma)
        f_err = r_err/self.cands[0].T
        if verbosity >= 1:
            print "\nSearching for duplicate harmonics..."
        ii = 0
        while 1:
            fundcand = self.cands[ii]
            jj = len(self.cands) - 1
            zapj = 0
            while 1:
                harmcand = self.cands[jj]
                if zapj:  print "Hey!"
                for factor in Num.arange(1.0, 17.0):
                    if (Num.fabs(fundcand.f - harmcand.f*factor) < f_err or
                        Num.fabs(fundcand.f - harmcand.f/factor) < f_err):
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        zapj = 1
                        break
                # Check a few other common ratios
                for factor in [3.0/2.0, 5.0/2.0,
                               2.0/3.0, 4.0/3.0, 5.0/3.0,
                               3.0/4.0, 5.0/4.0,
                               2.0/5.0, 3.0/5.0, 4.0/5.0]:
                    if Num.fabs(fundcand.f-harmcand.f*factor) < f_err:
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        zapj = 1
                        break
                if zapj:
                    harmcand.note = "This candidate (P=%.4f s, DM=%.2f) is " \
                                    "a harmonic of %s:%d (P=%.4f s, DM=%.2f)." % \
                                (harmcand.p, harmcand.DM, fundcand.filename, \
                                    fundcand.candnum, fundcand.p, fundcand.DM)
                    numremoved += 1
                    self.harmonic_cands.append(self.cands.pop(jj))
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (harmcand.filename, harmcand.candnum, jj)
                        print "    %s" % harmcand.note
                    zapj = 0
                jj -= 1
                if jj == ii:
                    break
            ii += 1
            if ii >= len(self.cands) - 1:
                break
        if verbosity >= 1:
            print "Removed a total of %d harmonics.\n" % numremoved

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
        dms = Num.unique([float(dm) for dm in dmlist])
        dmstrs = ['%.2f'%dm for dm in dms]
        dmdict = dict(zip(dmstrs, range(len(dms))))
        numremoved = 0
        num_toofew = 0
        num_toolow = 0
        num_gaps = 0
        self.cands.sort(cmp_sigma)
        for ii in reversed(range(len(self.cands))):
            currcand = self.cands[ii]
            # Remove all the candidates without enough DM hits
            if len(currcand.hits) < numdms:
                numremoved += 1
                num_toofew += 1
                currcand.note = "Candidate has only %d DM hits. This is less " \
                                "than minimum for 'good' cands (%d hits)" % \
                                (len(currcand.hits), numdms)
                self.dmproblem_cands.append(self.cands.pop(ii))
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            # Remove all the candidates where the max SNR DM is 
            # less than the cutoff DM
            imax = Num.argmax(Num.array([hit[1] for hit in currcand.hits]))
            hitdm, hitsnr = currcand.hits[imax]
            if float(hitdm) <= low_DM_cutoff:
                numremoved += 1
                num_toolow += 1
                currcand.note = "Hit with max SNR (%g) has dm (%.2f) " \
                                "<= low DM cutoff (%.2f) " % \
                                    (hitsnr, hitdm, low_DM_cutoff)
                self.dmproblem_cands.append(self.cands.pop(ii))
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            # Remove all the candidates where there are no hits at consecutive DMs
            if len(currcand.hits) > 1:
                currcand.hits.sort(cmp_dms)
                dm_indices = Num.asarray([dmdict["%.2f"%currcand.hits[jj][0]]
                                          for jj in range(len(currcand.hits))])
                min_dmind_diff = min(dm_indices[1:] - dm_indices[:-1])
                if min_dmind_diff > 1:
                    numremoved += 1
                    num_gaps += 1
                    currcand.note = "DM list of hits has gaps (i.e. " \
                                    "consecutive DMs don't have hits)."
                    self.dmproblem_cands.append(self.cands.pop(ii))
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (currcand.filename, currcand.candnum, ii)
                        print "    %s" % currcand.note
                    continue

        if verbosity >= 1:
            print "Removed %d candidates with DM problems.\n" % numremoved
        if verbosity >= 2:
            print "  # with too few hits:", num_toofew
            print "  # with peak SNR too low:", num_toolow
            print "  # with gaps in DM hits:", num_gaps

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
        else:
            summaryfile = open(summaryfilenm, "w")
        summaryfile.write("   Candlist contains %d 'good' candidates\n" % \
                            len(self.cands))
        summaryfile.write("      # Known RFI rejects:           %d\n" % \
              len(self.badcands_knownbirds))
        summaryfile.write("      # Short period rejects:        %d\n" % \
              len(self.badcands_shortperiod))
        summaryfile.write("      # Long period rejects:         %d\n" % \
              len(self.badcands_longperiod))
        summaryfile.write("      # Missed threshold:            %d\n" % \
              len(self.badcands_threshold))
        summaryfile.write("      # No good harmonics:           %d\n" % \
              len(self.badcands_harmpowcutoff))
        summaryfile.write("      # One bad harmonic:            %d\n" % \
              len(self.badcands_rogueharmpow))
        summaryfile.write("      # Duplicate candidates:        %d\n" % \
              len(self.duplicate_cands))
        summaryfile.write("      # Harmonic candidates:         %d\n" % \
              len(self.harmonic_cands))
        summaryfile.write("      # Candidates with DM problems: %d\n" % \
              len(self.dmproblem_cands))
        if summaryfilenm is not None:
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
        self.badcands_knownbirds.extend(other.badcands_knownbirds)
        self.badcands_longperiod.extend(other.badcands_longperiod)
        self.badcands_shortperiod.extend(other.badcands_shortperiod)
        self.badcands_threshold.extend(other.badcands_threshold)
        self.badcands_harmpowcutoff.extend(other.badcands_harmpowcutoff)
        self.badcands_rogueharmpow.extend(other.badcands_rogueharmpow)
        self.duplicate_cands.extend(other.duplicate_cands)
        self.harmonic_cands.extend(other.harmonic_cands)
        self.dmproblem_cands.extend(other.dmproblem_cands)

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
                goodcand.hits.sort(cmp_dms)
                for hit in goodcand.hits:
                    numstars = int(hit[1]/3.0)
                    candfile.write("  DM=%6.2f SNR=%5.2f   "%hit + numstars*'*' + '\n')
        if candfilenm is not None:
            candfile.close()


def candlist_from_candfile(filename):
    # First identify the length of the observation searched
    for line in open(filename, 'r'):
        if line.startswith(" Number of bins in the time series"):
            numsamp = int(line.split()[-1])
        if line.startswith(" Width of each time series bin (sec)"):
            dt = float(line.split()[-1])
    tobs = numsamp * dt

    cands = []
    candnums = []
    current_goodcandnum = 0
    for line in open(filename, 'r'):
        # Identify the candidates in the top of the file
        if fund_re.match(line):
            split_line = line.split()
            candnum   = int(split_line[0])
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
            continue

        # Parse the harmonic powers
        elif harms_re.match(line):
            split_line = line.split()
            candnum = int(split_line[0])

            if candnum in candnums:
                cand = cands[candnums.index(candnum)]
                cand.harm_pows = Num.zeros(cand.numharm, dtype=Num.float64)
                cand.harm_amps = Num.zeros(cand.numharm, dtype=Num.complex64) 
                power = parse_power(split_line[3])
                phase = float(split_line[9].split("(")[0])
                cand.harm_pows[0] = power
                cand.harm_amps[0] = Num.sqrt(power) * Num.exp(phase*1.0j)
                if (cand.numharm > 1):
                    current_goodcandnum = candnum
                    current_harmnum = 1
                else:
                    current_goodcandnum = 0
                    # Compute the S/N
                    cand.harms_to_snr()
                    # Now that S/N is available
                    # List candidate as a hit of itself
                    cand.hits = [(cand.DM, cand.snr)]
                    # These are the "optimized" power...
                    opt_ipow = cand.harm_pows[0]
                    # and sigma (calculated assuming _1_ trial!)
                    opt_sigma = candidate_sigma(opt_ipow, 1, 1)
                    cand.sigma = opt_sigma
                    cand.ipow_det = opt_ipow
            continue

        # Parse the higher (than the first) harmonic powers
        if current_goodcandnum:
            cand = cands[candnums.index(current_goodcandnum)]
            power = parse_power(line.split()[2])
            phase = float(line.split()[8].split("(")[0])
            cand.harm_pows[current_harmnum] = power
            cand.harm_amps[current_harmnum] = Num.sqrt(power) * Num.exp(phase*1.0j)
            current_harmnum += 1
            # Calculate other stats after all harmonics have been read in
            if (current_harmnum==cand.numharm):
                # Compute the S/N
                cand.harms_to_snr()
                # Now that S/N is available
                # List candidate as a hit of itself
                cand.hits = [(cand.DM, cand.snr)]
                # Compute sigma and incoherent power
                opt_ipow = sum(cand.harm_pows)
                opt_sigma = candidate_sigma(opt_ipow, cand.numharm, 1)
                cand.sigma = opt_sigma
                cand.ipow_det = opt_ipow
                current_goodcandnum = 0
    return Candlist(cands)


def read_candidates(filenms, prelim_reject=True):
    """Read in accelsearch candidates from the test ACCEL files.
        Return a Candlist object of Candidate instances.

        Inputs:
            filenms: A list of files to read candidates from.
            prelim_reject: If True, perform preliminary rejection of
                candidates. (Default: True)

    """
    if not len(filenms):
        print "Error:  There are no candidate files to read!"
        return None
    print "\nReading candidates from %d files...." % len(filenms)
    candlist = Candlist()
    for filenm in filenms:
        curr_candlist = candlist_from_candfile(filenm)
        candlist.extend(curr_candlist)
    if prelim_reject:
        candlist.default_rejection()
    return candlist


def sift_directory(dir, outbasenm):
    """Sift candidates in given directory.
        
        Inputs:
            dir: The directory containing candidates to sift.
            outbasenm: The base name of the output files.

        Outputs:
            None
    """
    import os.path
    import glob
    import infodata
  
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
    lo_accel_cands = read_candidates(lo_accel_fns)
    print "Read %d candidates from %d files" % \
                (len(lo_accel_cands), len(lo_accel_fns))
    print "%d candidates passed default rejection" % len(lo_accel_cands)
    if len(lo_accel_cands):
        lo_accel_cands.remove_duplicate_candidates()
    if len(lo_accel_cands):
        lo_accel_cands.remove_DM_problems(2, dmstrs, 2.0)
    lo_accel_cands.print_cand_summary()

    # Read candidates found in high acceleration searching
    hi_accel_fns = glob.glob(os.path.join(dir, '*ACCEL_50'))
    hi_accel_cands = read_candidates(hi_accel_fns)
    print "Read %d candidates from %d files" % \
                (len(hi_accel_cands), len(hi_accel_fns))
    print "%d candidates passed default rejection" % len(hi_accel_cands)
    if len(hi_accel_cands):
        hi_accel_cands.remove_duplicate_candidates()
    if len(hi_accel_cands):
        hi_accel_cands.remove_DM_problems(2, dmstrs, 2.0)
    hi_accel_cands.print_cand_summary()

    all_accel_cands = lo_accel_cands + hi_accel_cands
    if len(all_accel_cands):
        all_accel_cands.remove_harmonics()
        # Note:  the candidates will be sorted in _sigma_ order, not _SNR_!
        all_accel_cands.cands.sort(cmp_sigma)
        print "Found %d good candidates" % len(all_accel_cands)
        all_accel_cands.to_file(outbasenm+".accelcands")
    all_accel_cands.write_cand_report(outbasenm+".accelcands.report")
    all_accel_cands.print_cand_summary()
    all_accel_cands.plot_goodbad()
    plt.savefig(outbasenm+".accelcands.rejects.png")
    all_accel_cands.plot_summary()
    plt.savefig(outbasenm+".accelcands.summary.png")


def main():
    # Sift candidates in PWD
    outbasenm = sys.argv[1]
    sift_directory(os.getcwd(), outbasenm)


if __name__ == '__main__':
    main()
