#!/usr/bin/env python

"""
rrattrap.py 

Single pulse sifting code: performs grouping and ranking of single pulses identified by PRESTO's single_pulse_search.py.
Usage:
    python rrattrap.py [OPTIONS] --inffile <.inf file> <.singlepulse files>

Output:
    - groups.txt : a file listing all single pulse groups and their ranking.
    - several colourized DM vs. time single-pulse plots, for different DM ranges, with colours corresponding to group ratings.

Chen Karako May 7, 2014
Updated by Chitrang Patel June 10, 2016.
"""
from __future__ import print_function
#from past.builtins import cmp
from time import strftime
from presto import infodata
import matplotlib.pyplot as plt
import numpy as np
from presto.Pgplot import *
import optparse
from presto.singlepulse import spio

FRACTIONAL_SIGMA = 0.9 # change to 0.8?
ALL_RANKS_ORDERED = [1,2,0,3,4,5,6]
DEBUG = True # if True, will be verbose


def dmthreshold(dm, use_dmplan, min_group=45):
    if use_dmplan:
        import singlepulse.rrattrap_config as rrattrap_config
        dmt, min_group = rrattrap_config.use_dmplan(dm)
    else:
        dmt = 1
        min_group = min_group 
    return dmt, min_group


class SinglePulseGroup(object): # Greg's modification
    """Define single pulse group
    """
    __slots__ = ['min_dm', 'max_dm', 'max_sigma', 'center_time', 
                 'min_time', 'max_time', 'duration', 
                 'singlepulses', 'numpulses', 'rank'] # Greg's modification
    
    def __init__(self, dm, sigma, time, sample, downfact):
        """SinglePulseGroup constructor.
            Takes as input one single pulse (creates a group of one)
            inputs DM,Sigma,Time,Sample,Downfact.
        """
        self.min_dm = dm
        self.max_dm = dm
        self.max_sigma = sigma
        self.center_time = time
        if sample == 0:
            dt = 0 # this will ignore events with sample=0. 
                   # better would be to use the inf files 
                   # to determine the dt for these events
        else:
            dt = time/sample
        self.min_time = time-downfact/2.0*dt
        self.max_time = time+downfact/2.0*dt
        self.duration = self.max_time - self.min_time
        self.singlepulses = [(dm,sigma,time,sample,downfact)]
        self.numpulses = 1
        self.rank = 0


    # Python 3 no longer honors the __cmp__ special method, so another way of comparing is required
    # see: https://portingguide.readthedocs.io/en/latest/comparisons.html
    # def __cmp__(self, other):
    #     return cmp(ALL_RANKS_ORDERED.index(self.rank),
    #                ALL_RANKS_ORDERED.index(other.rank))

    def timeisclose(self,other,use_dmplan,time_thresh=0.5):
        """Checks whether the overlap in time of self and other is within
            time_thresh. Takes as input other, a SinglePulseGroup object,
            as well as the optional input time_thresh (in s).
        """
        if self.duration < other.duration:
            narrow = self
            wide = other
        else:
            narrow = other
            wide = self
        time_thresh = dmthreshold(self.min_dm, use_dmplan)[0]*time_thresh
        dt = max(time_thresh, narrow.duration/2.0) # always group groups within time_thresh (or duration/2, if longer) of one another
        timeisclose = (wide.max_time >= (narrow.center_time - dt)) and\
                        (wide.min_time <= (narrow.center_time + dt))
        
        return timeisclose

    def dmisclose(self,other,use_dmplan,dm_thresh=0.1):
        """Checks whether the DM of self and other is within dm_thresh of one
            another. Takes as input other, a SinglePulseGroup object, as well as the optional input dm_thresh (in pc cm-3).
        """
        dm_thresh = dmthreshold(self.min_dm, use_dmplan)[0]*dm_thresh
        dmisclose = (other.max_dm >= (self.min_dm-dm_thresh)) and\
                    (other.min_dm <= (self.max_dm+dm_thresh))

        return dmisclose

    def combine(self,other):
        """combines self and other SinglePulseGroup objects.
            takes as input other, a SinglePulseGroup object.
            combines in place; nothing returned.
        """
        self.min_dm = min(self.min_dm, other.min_dm)
        self.max_dm = max(self.max_dm, other.max_dm)
        self.min_time = min(self.min_time, other.min_time)
        self.max_time = max(self.max_time, other.max_time)
        self.max_sigma = max(self.max_sigma, other.max_sigma)
        self.duration = self.max_time - self.min_time
        self.center_time = (self.min_time + self.max_time)/2.0
        self.numpulses = self.numpulses + other.numpulses
        self.singlepulses.extend(other.singlepulses)
    
    def __str__(self):
        s = ["Group of %d single pulses: " % len(self.singlepulses), \
             "\tMin DM (cm-3 pc): %f" % self.min_dm, \
             "\tMax DM (cm-3 pc): %f" % self.max_dm, \
             "\tCenter time (s):  %f" % self.center_time, \
             "\tDuration (s):     %f" % self.duration, \
             "\tMax sigma:        %f" % self.max_sigma, \
             "\tRank:             %f" % self.rank]
        return '\n'.join(s)

def create_groups(sps, inffile, min_nearby=1, time_thresh=0.5, \
                    dm_thresh=0.1, ignore_obs_end=0, use_dmplan=False):
    """Given a recarray of singlepulses return a list of
        SinglePulseGroup objects.

        Inputs:
            sps: A recarray of single pulse info.
            min_nearby: Minimum number of nearby single pulse events
                to bother creating a group.
            time_thresh: Time-range within which another event must be found
            dm_thresh: DM-range within which another event must be found
            ignore_obs_end: if non-zero, the time (in seconds) to ignore from 
                the end of the observation. Useful for beams on which zero-DMing 
                was applied and caused artifacts in the sp output at the end of 
                the obs.

            *** NOTE: time_thresh and dm_thresh are used together

        Outputs:
            groups: A list of SinglePulseGroup objects.
    """

    Tobs = get_obs_info(inffile)['T'] # duration of observation
    if not (0 <= ignore_obs_end < Tobs):
        print("Invalid ignore_obs_end value. Value must be: \
            0 <= ignore_obs_end < Tobs. Setting ignore_obs_end to 0.")
        ignore_obs_end = 0
    Tignore = Tobs - ignore_obs_end # sps with t>=Tignore will be ignored

    numsps = len(sps)
    # Sort by time
    sps.sort(order='time')
    groups = []
    for ii in range(numsps):
        # Time and DM of current SP event
        ctime = sps[ii]['time']
        if ignore_obs_end and (ctime > Tignore):
            continue
        cdm = sps[ii]['dm']
        ngood = 0 # number of good neighbours
        time_thresh = dmthreshold(cdm, use_dmplan)[0]*time_thresh
        dm_thresh = dmthreshold(cdm, use_dmplan)[0]*dm_thresh
        
        jj = ii+1
        while (ngood < min_nearby) and (jj < numsps) and \
                    (sps[jj]['time'] < (ctime+time_thresh)):
            if abs(sps[jj]['dm'] - cdm) < dm_thresh:
                ngood += 1
            jj += 1
        # Look backward as well
        jj = ii-1
        while (ngood < min_nearby) and (jj >= 0) and \
                    (sps[jj]['time'] > (ctime-time_thresh)):
            if abs(sps[jj]['dm'] - cdm) < dm_thresh:
                ngood += 1
            jj -= 1
        if ngood >= min_nearby:
            # At least min_nearby nearby SP events
            grp = SinglePulseGroup(*sps[ii])
            groups.append(grp)
    return groups


def grouping_sp_dmt(groups, use_dmplan=False, time_thresh=0.5, dm_thresh=0.1):
    """Groups SinglePulse objects based on proximity in time, DM. 
        Outputs list of Single Pulse Groups.
    """
    didcombine = True
    while didcombine:
        didcombine = False
        groups.sort(key=lambda group: group.min_time) #Sort by time
        for i, grp1 in enumerate(groups):
            j=i+1
            while (j<len(groups) and groups[i].center_time+0.2 > groups[j].center_time): #Only look at groups that are close in time
               if grp1.dmisclose(groups[j], use_dmplan, dm_thresh):
                    if grp1.timeisclose(groups[j], use_dmplan, time_thresh):
                        grp1.combine(groups.pop(j))
                        didcombine = True
               j=j+1


def grouping_rfi(groups, use_dmplan=False, time_thresh=0.5, dm_thresh=0.1):
    """
    Groups together close groups of RFI, and considers as RFI other groups
    that are close to RFI.
    """
    didcombine = True
    while didcombine:
        didcombine = False
        # If a group is very close to a group of rfi, set it as rfi  
        for i in reversed(list(range(len(groups)))):
            grp1 = groups[i]
            for j in reversed(list(range(len(groups)))):
                if j <= i:
                    continue
                grp2 = groups[j]
                if (grp1.rank != 2) and (grp2.rank != 2):
                    continue
                if grp1.dmisclose(grp2,use_dmplan,10) and grp1.timeisclose(grp2, use_dmplan, time_thresh): # use bigger time thresh?
                    grp1.combine(groups.pop(j))
                    # FIXME: Should we set as RFI without checking
                    #        sigma behaviour (ie re-check rank) for group?
                    grp1.rank = 2 # Set as rfi
                    didcombine = True


def grouping_sp_t(groups, use_dmplan=False, time_thresh=0.5, dm_thresh=0.1):
    """Groups SinglePulse objects based on proximity in time, assuming 
        the DM difference is no more than DMDIFF=10.

        Inputs:
            groups: A list of SinglePulseGroup objects.

        Outputs:
            groups: A list of SinglePulseGroup objects.
    """
    DMDIFF = 10 #max allowed DM difference between groups that will be grouped in time
    didcombine = True
    while didcombine:
        didcombine = False
        for i, grp1 in enumerate(groups):
            for j in range(len(groups)-1,i,-1):
                if grp1.timeisclose(groups[j], use_dmplan, time_thresh) and \
                    grp1.dmisclose(groups[j],use_dmplan,DMDIFF): # We check if two events
                                                      # have similar time and 
                                                      # a DM difference < DMDIFF
                    grp1.combine(groups.pop(j)) # Note group rank is not 
                                                # updated when combine groups,
                                                # need to re-run ranking after.
                    didcombine = True
    return groups


def flag_noise(groups, use_dmplan=False, min_group=45):
    """Flag groups as noise based on group size.
        If the number of sp events in a group is < min_group,
        this group is marked as noise.

        Inputs:
            groups: A list of SinglePulseGroup objects.
            min_group: The minimum group size that a group must have
                        in order not to be considered as noise. The
                        default min_group is MIN_GROUP.

        Outputs:
            None
    """
    for grp in groups:
        min_group = dmthreshold(grp.min_dm, use_dmplan, min_group)[1]
        if grp.numpulses < min_group:
            grp.rank = 1
    return groups


def flag_rfi(groups, close_dm = 2.0):
    """Flag groups as RFI based on sigma behavior.
        Takes as input list of Single Pulse Groups.
        The ranks of the groups are updated in-place.

        Inputs:
            groups: A list of SinglePulseGroup objects.

        Outputs:
            None
    """
    for grp in groups:
        if (grp.rank != 2) and (grp.min_dm <= close_dm): # if grp has not 
                                                         # yet been marked RFI
            for sp in grp.singlepulses:
                if (sp[0] <= close_dm) and \
                    (sp[1] >= (FRACTIONAL_SIGMA*grp.max_sigma)):
                    # if any sp in the group has low dm, and its sigma is >= frac sigma*grp.max_sigma, call that grp rfi
                    grp.rank = 2 
                    break


def rank_groups(groups, use_dmplan=False, min_group=45, min_sigma=8.0):
    """Rank groups based on their sigma vs. DM behaviour. 
        Takes as input list of Single Pulse Groups.
        The ranks of the groups are updated in-place.

        Inputs:
            groups: A list of SinglePulseGroup objects.

        Outputs:
            None
    """
#   divide groups into 5 parts (based on number events) to examine sigma behaviour
    for grp in groups:
        min_group = dmthreshold(grp.min_dm, use_dmplan, min_group)[1]
        if len(grp.singlepulses) < min_group:
            grp.rank = 1
        elif grp.rank != 2: # don't overwrite ranks of rfi groups
            numsps = len(grp.singlepulses)
            # sort list by increasing DM
            idmsort = np.argsort([sp[0] for sp in grp.singlepulses])
            
            sigmas = np.ma.zeros(int(np.ceil(numsps/5.0)*5))
            sigmas[-numsps:] = np.asarray([sp[1] for sp in grp.singlepulses])[idmsort]
            # Mask sigma=0. These are elements added to pad size of array
            # to have multiple of 5 elements
            # (there should never be actual SPs with sigma=0)
            sigmas = np.ma.masked_equal(sigmas, 0.0)
         
            sigmas.shape = (5, int(np.ceil(numsps/5.0)))
           
            maxsigmas = sigmas.max(axis=1)
            avgsigmas = sigmas.mean(axis=1)
            # standard deviation in signal to noise values in the group. 
            stdsigmas = sigmas.std(axis=1)
            # The largest std deviation
            maxstd = np.ma.max(stdsigmas)
            # The smallest std deviation
            minstd = np.ma.min(stdsigmas)
            # The largest maxsigma
            maxsigma = np.ma.max(maxsigmas)
            # The smallest maxsigma
            minsigma = np.ma.min(maxsigmas)
            # The largest avgsigma
            maxavgsigma = np.ma.max(avgsigmas)
            # The smallest avgsigma
            minavgsigma = np.ma.min(avgsigmas)
                       
            #if maxavgsigma<1.05*minavgsigma:
            if all(stdsigma < 0.1 for stdsigma in stdsigmas): 
                # Sigmas pretty much constant. Group is RFI
                grp.rank = 2
            if maxsigmas[2] > maxsigmas[1]:
                if maxsigmas[2] > maxsigmas[3]:
                    # nearest neighbour subgroups both have smaller sigma
                    grp.rank = 3
                    if (maxsigmas[3] > maxsigmas[4]) and (maxsigmas[1] > maxsigmas[0]): 
                        #next-nearest subgps have sigma < nearest neighbours
                        grp.rank = 4
                        if maxsigmas[2] > min_sigma:  
                            # We want the largest maxsigma to be at least 
                            # 1.15 times bigger than the smallest
                            grp.rank = 5
                            if (avgsigmas[2] > avgsigmas[0]) and \
                                (avgsigmas[2] > avgsigmas[4]) and \
                                maxsigma>1.15*minsigma:
                                    grp.rank = 6 
                else: #ie. maxsigmas[2] <= maxsigmas[3], allowing for asymmetry:
                    if maxsigmas[1] > maxsigmas[0]:
                        grp.rank = 3
                        if maxsigmas[3] > maxsigmas[4]:
                            grp.rank = 4
                            if maxsigmas[3] > min_sigma:
                                grp.rank = 5
                                if (avgsigmas[3] > avgsigmas[0]) and \
                                    (avgsigmas[3] > avgsigmas[4]) and \
                                    maxsigma>1.15*minsigma:
                                   grp.rank = 6 
            else: #ie. maxsigma2 >= maxsigma3, allowing for asymmetry:
                if (maxsigmas[1] > maxsigmas[0]) and (maxsigmas[2] > maxsigmas[3]):
                    grp.rank = 3
                    if maxsigmas[3] > maxsigmas[4]:
                        grp.rank = 4
                        if maxsigmas[1] > min_sigma:
                            grp.rank = 5
                            if (avgsigmas[1] >= avgsigmas[0]) and \
                                (avgsigmas[1] > avgsigmas[4]) and \
                                maxsigma>1.15*minsigma:
                                grp.rank = 6
            if any(stdsigma < 0.1 for stdsigma in stdsigmas) and (grp.max_sigma < 5.5): # if max sigma of the group is less than 5.5 and the sigma distribution is mostly flat, then it is not likely to be astrophysical.
                grp.rank = 0
            if grp.rank == 0:
                pass 

def check_dmspan(groups, dt, lofreq, hifreq):
    """Read in groups and check whether each group's DM span exceeds the threshold.
    """
    for grp in groups:
        for sp in grp.singlepulses:
            if sp[1] == grp.max_sigma:
                downsamp = (sp[2]/dt)/sp[3]
                width_ms = 1000.0*sp[4]*dt*downsamp
                break
        if (grp.max_dm-grp.min_dm > 5*spio.theoritical_dmspan(grp.max_sigma, 5.0, width_ms, lofreq, hifreq)): 
            # checks if the DM span is more than 5 times theoritical dm value.
            if not ((grp.rank == 5) or (grp.rank == 6)): #if group is not good or excellent
                grp.rank = 2                             # then its most likely RFI.

def get_obs_info(inffile):
    """Read in an .inf file to extract observation information.
        Return observation RA, Dec, duration, and source name.
    """
    inf = infodata.infodata(inffile)
    T = inf.dt * inf.N # total observation time (s)
    RA = inf.RA
    dec = inf.DEC
    src = inf.object
    MJD = inf.epoch
    telescope = inf.telescope
    freq = (inf.numchan/2-0.5)*inf.chan_width+inf.lofreq # center freq
    return {'T': T, 'RA': RA, 'dec': dec, 'src': src, 'MJD': MJD, 'telescope': telescope, 'freq': freq}


def plot_sp_rated_all(groups, ranks, inffile, ylow=0, yhigh=100, xlow=0, xhigh=120):
    """Take in dict of Single Pulse Group lists and 
        plot the DM vs. t for all, with the plotted 
        colour corresponding to group rank. 
        The DM range to plot can also be specified.
    """
    rank_to_color = {2:'darkgrey', 0:'k', 3:'c', 4:'royalblue', 5:'b', 6:'m'}

    # Prepare data to plot
    dm = [] 
    time = []
    size = []
    colors = []
    for grp in groups:
        if grp.rank not in ranks:
            continue
        if grp.min_dm < yhigh and grp.max_dm > ylow:
            for sp in grp.singlepulses:
                dm.append(sp[0])
                size.append(np.clip(3*sp[1]-14,0,50))
                time.append(sp[2])
                colors.append(rank_to_color[grp.rank])

    # Plot
    plt.axes()
    if len(time): # check if there are points to plot
        plt.scatter(time, dm, c=colors, marker='o', s=size, edgecolor='none')
    plt.xlabel('Time (s)')
    plt.ylabel('DM (pc cm$^{-3}$)')
    # if inf file exists, will override xlow and xhigh 
    # specified when function is called
    if get_obs_info(inffile) is not None: # if inf files exist, can get obs info
        obsinfo = get_obs_info(inffile)
        plt.title('Single Pulse Results for %s\nRA: %s Dec: %s' % 
                  (obsinfo['src'], obsinfo['RA'], obsinfo['dec']))
        xhigh = obsinfo['T'] # set xhigh to observation duration
    plt.xlim((xlow, xhigh)) 
    plt.ylim((ylow, yhigh))

    print_debug("Saving figure...")
    plt.savefig('grouped_sps_DMs%s-%s.png' % (ylow, yhigh), dpi=300)


def plot_sp_rated_pgplot(groups, ranks, inffile, ylow=0, yhigh=100, xlow=0, xhigh=120):
    """Plot groups according to their ranks. Uses pgplot rather 
        than matplotlib for faster, more memory-efficient plotting.

        Inputs:
            groups: A list of SinglePulseGroup objects.
            ylow (optional): lower y limit to plot. Default: ylow=0.
            yhigh (optional): higher y limit to plot. Default: yhigh=100.
            xlow (optional): lower x limit to plot. Default: xlow=0.
            xhigh (optional): higher x limit to plot. Default: use inf file
                to find observation duration, or if inf file unavailable,
                use xhigh=120s.

        Outputs:
            None; saves a colorized sp plot.
    """
    if get_obs_info(inffile) is not None: # if inf files exist, can get obs info
        obsinfo = get_obs_info(inffile)
        #plt.title('Single Pulse Results for %s\nRA: %s Dec: %s' % 
                  #(obsinfo['src'], obsinfo['RA'], obsinfo['dec']))
        xhigh = obsinfo['T'] # set xhigh to observation duration
    ppgplot.pgopen('grouped_pgplot_DMs%s-%s.ps/vcps' % (ylow, yhigh))
    # copying from single_pulse_search.py plotting
    ppgplot.pgpap(9, 0.7) # Width in inches, aspect
#    ppgplot.pgsvp(0.06, 0.97, 0.08, 0.52) # single_pulse_search.py
    ppgplot.pgsvp(0.87, 0.97, 0.80, 0.97)
    ppgplot.pgsch(0.9)
    ppgplot.pgsci(15)
    ppgplot.pgpt(np.array([0.01]), np.array([0.12]), 9)
    ppgplot.pgptxt(0.1, 0.1, 0.0, 0.0, '2: RFI')
    ppgplot.pgsci(1)
    ppgplot.pgpt(np.array([0.01]), np.array([0.27]), 9)
    ppgplot.pgptxt(0.1, 0.25, 0.0, 0.0, '0: Noise')
    ppgplot.pgsci(5)
    ppgplot.pgpt(np.array([0.01]), np.array([0.42]), 9)
    ppgplot.pgptxt(0.1, 0.4, 0.0, 0.0, '3: Ok')
    ppgplot.pgsci(11)
    ppgplot.pgpt(np.array([0.01]), np.array([0.57]), 9)
    ppgplot.pgptxt(0.1, 0.55, 0.0, 0.0, '4: Good')
    ppgplot.pgsci(4)
    ppgplot.pgpt(np.array([0.01]), np.array([0.72]), 9)
    ppgplot.pgptxt(0.1, 0.7, 0.0, 0.0, '5:Very Good')
    ppgplot.pgsci(6)
    ppgplot.pgpt(np.array([0.01]), np.array([0.87]), 9)
    ppgplot.pgptxt(0.1, 0.85, 0.0, 0.0, '6: Excellent')
    ppgplot.pgsci(1)

    ppgplot.pgsvp(0.06, 0.97, 0.08, 0.80)
    #ppgplot.pgenv(xlow, xhigh, ylow, yhigh, 0, 1) #not sure if need 0,1
    ppgplot.pgswin(xlow, xhigh, ylow, yhigh)
    ppgplot.pgsch(0.8)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0) # redundant with pgenv
    ppgplot.pgmtxt('B', 2.5, 0.5, 0.5, "Time (s)")
    ppgplot.pgmtxt('L', 1.8, 0.5, 0.5, "DM (pc cm\\u-3\\d)")
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt('T', 2.5, 0.3, 0.0,
                   "Single Pulse Results for %s" % obsinfo['src'])
    ppgplot.pgmtxt('T', 1.5, 0.3, 0.0, 'RA: %s  Dec: %s' \
                    % (obsinfo['RA'], obsinfo['dec']))
    ppgplot.pgmtxt('T', 0.5, 0.2, 0.0, 'Telescope: %s  MJD: %s Freq: %.1f MHz'\
                    % (obsinfo['telescope'], obsinfo['MJD'], obsinfo['freq']))
    ppgplot.pgsch(0.8)
    rank_to_color = {2:15, # grey
                     0:1, # black
                     3:5, # cyan
                     4:11, # dim blue
                     5:4, # dark blue
                     6:6} # magenta
    
    # Plotting scheme taken from single_pulse_search.py
    # Circles are symbols 20-26 in increasing order
    snr_range = 12.0
    spthresh = 5.0 # 5 for gbncc (5.5 plotting), 6 for palfa...

    for grp in groups:
        cand_symbols = []
        dm = []
        time = []
        if grp.rank not in ranks:
            continue
        if grp.min_dm < yhigh and grp.max_dm > ylow:
            ppgplot.pgsci(rank_to_color[grp.rank])
            for sp in grp.singlepulses:
                dm.append(sp[0])
                time.append(sp[2])
                cand_SNR = sp[1]
                # DEBUG: UNCOMMENT first line, then remove next 2 lines
                if np.isfinite(cand_SNR):
                    cand_symbol = int((cand_SNR-spthresh)/snr_range * 6.0 + 20.5)
                else:
                    cand_symbol = 26
                cand_symbols.append(min(cand_symbol, 26)) # biggest circle is 26
        cand_symbols = np.array(cand_symbols)
        dm = np.array(dm)
        time = np.array(time)
        for ii in [26, 25, 24, 23, 22, 21, 20]:
            inds = np.nonzero(cand_symbols==ii)[0]
            ppgplot.pgpt(time[inds], dm[inds], ii)
    ppgplot.pgclos()


def print_debug(msg):
    if DEBUG:
        print(msg)
#        print h.heap()


def pop_by_rank(groups, rank):
    """Remove groups with specified rank from a list of groups.
        Removes the groups in place; nothing returned.

        Inputs:
            groups: A list of SinglePulseGroup objects.
            rank: The rank of groups to be removed.

        Outputs:
            None
    """
    for j in reversed(list(range(len(groups)))):
        if groups[j].rank == rank:
           del groups[j] 
    

def rank_occur(groups):
    """Return a dict of the number of groups of each rank in the groups list.

        Inputs:
            groups: A list of SinglePulseGroup objects.

        Outputs:
            rank_occur: A dict of ranks and the number of their occurrences
                        in the groups list.
    """
    rank_occur = {}
    for grp in groups:
        nn = rank_occur.setdefault(grp.rank, 0)
        rank_occur[grp.rank] = nn+1

    return rank_occur


#@profile
def main():
    parser = optparse.OptionParser(prog="rrattrap.py", \
                         version="Chen Karako, updated by Chitrang Patel(June 23, 2015)",\
                         usage="%prog --inffile <.inf file> [options] *.singlepulse",\
                         description="Group single pulse events and rank them based on the sigma behavior. \
                                       Plot DM vs time with different colours for different ranks as follows:\
                                       \t\tRank 1 (Other)      : Grey\
                                       \t\tRank 2 (RFI)        : Red\
                                       \t\tRank 3 (ok)         : Cyan\
                                       \t\tRank 4 (good)       : dim blue\
                                       \t\tRank 5 (very good)  : dark blue\
                                       \t\tRank 6 (excellent)  : Magenta")

    parser.add_option('--CLOSE-DM', dest='close_dm', type='float', \
                        help="DM to below which the signalis considered RFI(Default: 2", \
                        default=2.0)
    parser.add_option('--use-configfile', dest='use_configfile', action='store_true', \
                        help="If this flag is set - import the config file for selecting grouping" \
                        "parameters.(Default: do not use a config file.)", default=False)
    parser.add_option('--use-DMplan', dest='use_DMplan', action='store_true', \
                        help="If this flag is set - Use the ddplan for selecting grouping" \
                        "parameters. Make sure that you have a corresponding config file containing" \
                        "the DDplan.  (Default: do not use ddplan)", default=False)
    parser.add_option('--min-group', dest='min_group', type='int', \
                        help="minimum number of events in a group to no be considered noise." \
                             "(Default: 45)", \
                         default=45)
    parser.add_option('--dm-thresh', dest='dm_thresh', type='float', \
                        help="DM threshold to use for nearest neighbour. Suggest a value greater " \
                              " than the DM step size(Default: 0.5 pc/cm^3 - will not work if DM " \
                              "step size is greater than 0.5)", default=0.5)
    parser.add_option('--time-thresh', dest='time_thresh', type='float', \
                        help="Time threshold to use for nearest neighbour. Suggest a value that " \
                             " is a few times the max pulse width(Default: 0.1 s)", default=0.1)
    parser.add_option('--min-sigma', dest='min_sigma', type='float', \
                        help="minimum signal-to-noise above which the group is highly considered" \
                        "to be astrophysical. (Default: 8.0)", \
                         default=8.0)
    parser.add_option('--no-plot', dest='plot', action='store_false', \
                        help="Do not plot the groups in the DM time plot." \
                                "(Default: Make a plot)", default=True)
    parser.add_option('--plottype', dest='plottype', type = 'string',\
                       help="Make a plot using : 'matplotlib' or 'pgplot'."\
                       , default='pgplot')
    parser.add_option('--min-rank-to-plot', dest='min_ranktoplot', type = 'int',\
                       help="Only groups with rank upto this will plotted.(default: plot \
                       all except rank 1)", default=0)
    parser.add_option('--min-rank-to-write', dest='min_ranktowrite', type = 'int',\
                       help="Only info of the groups with rank upto this will written." \
                       "(default: write all except rank 1)", default=0)
    parser.add_option('--inffile', dest='inffile', type = 'string',\
                       help="A .inf file. I suggest a .rfifind.inf file."\
                       , default=None)
    parser.add_option('-o', dest='outbasenm', type = 'string',\
                       help="outfile base name. .groups.txt will be added to the given name."\
                       , default='')
    options, args = parser.parse_args()

    if not hasattr(options, 'inffile'):
        raise ValueError("You must supply a .inf file. I suggest .rfifind.inf")
    
    if not options.inffile.endswith(".inf"):
        raise ValueError("Cannot recognize file type from extension. "
                         " Only '.inf' types are supported.")
    
    if options.use_DMplan or options.use_configfile:
        import singlepulse.rrattrap_config as rrattrap_config

    RANKS = np.asarray([2,0,3,4,5,6])
    
    if options.use_configfile:
        CLOSE_DM = rrattrap_config.CLOSE_DM
        MIN_GROUP = rrattrap_config.MIN_GROUP
        TIME_THRESH = rrattrap_config.TIME_THRESH
        DM_THRESH = rrattrap_config.DM_THRESH
        MIN_SIGMA = rrattrap_config.MIN_SIGMA
        PLOT = rrattrap_config.PLOT
        PLOTTYPE = rrattrap_config.PLOTTYPE
        RANKS_TO_WRITE = rrattrap_config.RANKS_TO_WRITE
        RANKS_TO_PLOT = rrattrap_config.RANKS_TO_PLOT
    else:
        CLOSE_DM = options.close_dm
        MIN_GROUP = options.min_group
        TIME_THRESH = options.time_thresh
        DM_THRESH = options.dm_thresh
        MIN_SIGMA = options.min_sigma
        PLOT = options.plot
        PLOTTYPE = options.plottype
        RANKS_TO_WRITE = list(RANKS[RANKS>options.min_ranktowrite]) 
        RANKS_TO_PLOT = list(RANKS[RANKS>options.min_ranktoplot])
    
    inffile = options.inffile
    inf = infodata.infodata(inffile)    
    print_debug("Beginning read_sp_files... "+strftime("%Y-%m-%d %H:%M:%S"))
    
    groups = spio.read_sp_files(args[1:])[0]
    print_debug("Finished read_sp_files, beginning create_groups... " +
                strftime("%Y-%m-%d %H:%M:%S"))
    print_debug("Number of single pulse events: %d " % len(groups))
    
    groups = create_groups(groups, inffile, min_nearby=1, ignore_obs_end=10, time_thresh=TIME_THRESH, dm_thresh=DM_THRESH, use_dmplan=options.use_DMplan) # ignore the last 10 seconds of the obs, for palfa
    print_debug("Number of groups: %d " % len(groups))
    print_debug("Finished create_groups, beginning grouping_sp_dmt... " +
                    strftime("%Y-%m-%d %H:%M:%S"))
    
    grouping_sp_dmt(groups, use_dmplan=options.use_DMplan, time_thresh=TIME_THRESH, dm_thresh=DM_THRESH)
    print_debug("Number of groups (after initial grouping): %d " % len(groups))
    print_debug("Finished grouping_sp_dmt, beginning flag_noise... " + 
                strftime("%Y-%m-%d %H:%M:%S"))
    flag_noise(groups, use_dmplan=options.use_DMplan, min_group=MIN_GROUP) # do an initial coarse noise flagging and removal
    pop_by_rank(groups, 1)
    print_debug("Number of groups (after removed noise gps w <10 sps): %d " % len(groups))
    print_debug("Beginning grouping_sp_t... " +
                strftime("%Y-%m-%d %H:%M:%S"))
    
    # Regroup good groups based on proximity in time only (compensate for missing middles):
    groups = grouping_sp_t(groups, use_dmplan=options.use_DMplan, time_thresh=TIME_THRESH, dm_thresh=DM_THRESH)
    print_debug("Finished grouping_sp_t. " + strftime("%Y-%m-%d %H:%M:%S"))
    
    # Flag RFI groups, noise
    flag_rfi(groups, close_dm=CLOSE_DM)
    # Rank groups and identify noise (<45/40/35/30 sp events) groups
    
    print_debug("Ranking groups...")
    rank_groups(groups, use_dmplan=options.use_DMplan, min_group=MIN_GROUP, min_sigma=MIN_SIGMA)
    # Remove noise groups
    print_debug("Before removing noise, len(groups): %s" % len(groups))
    pop_by_rank(groups, 1)
    print_debug("After removing noise, len(groups): %s" % len(groups))
    
    # Group rfi with very close groups
    print_debug("len(groups) before grouping_rfi: %s" % len(groups))
    print_debug("Beginning grouping_rfi... " + strftime("%Y-%m-%d %H:%M:%S"))
    grouping_rfi(groups, use_dmplan=options.use_DMplan, time_thresh=TIME_THRESH, dm_thresh=DM_THRESH)
    print_debug("Finished grouping_rfi. " + 
                strftime("%Y-%m-%d %H:%M:%S"))
    
    # Rank groups
    print_debug("Finished rank_groups, beginning DM span check... " + 
                strftime("%Y-%m-%d %H:%M:%S"))
    # Remove groups that are likely RFI, based on their large span in DM
    print_debug("Beginning DM span check...")
    check_dmspan(groups, inf.dt, inf.lofreq, inf.lofreq+inf.BW)
    print_debug("Finished DM span check, beginning writing to outfile... " + 
                strftime("%Y-%m-%d %H:%M:%S"))

    outfile = open(options.outbasenm+'groups.txt', 'w')
    summaryfile = open(options.outbasenm+'spsummary.txt', 'w')
    
    rank_dict = rank_occur(groups)
    for rank in sorted(ALL_RANKS_ORDERED):
        if rank != 1:
            outfile.write("Number of rank %d groups: %d \n" % 
                      (rank, rank_dict.get(rank, 0)))
            summaryfile.write("Number of rank %d groups: %d \n" % 
                      (rank, rank_dict.get(rank, 0)))
    outfile.write("\n")
    summaryfile.close()

    # Reverse sort lists so good groups are written at the top of the file
    groups.sort(key=lambda x: ALL_RANKS_ORDERED.index(x.rank), reverse=True)

    # write list of events in each group
    for grp in groups:
        if grp.rank in RANKS_TO_WRITE:
            outfile.write(str(grp) + '\n') #print group summary
            outfile.write('\n')
            outfile.write("# DM      Sigma     Time (s)    Sample    Downfact \n")
            for sp in grp.singlepulses:
                outfile.write("%7.2f %7.2f %13.6f %10d   %3d \n" % sp)
            outfile.write('\n')
    outfile.close()

    print_debug("Finished writing to outfile, now plotting... " + 
                strftime("%Y-%m-%d %H:%M:%S"))
    
    if PLOT:
        ranks = RANKS_TO_PLOT 
        # Sort groups so better-ranked groups are plotted on top of worse groups
        groups.sort(key=lambda x: ALL_RANKS_ORDERED.index(x.rank))
        # create several DM vs t plots, splitting up DM in overlapping intervals 
        # DMs 0-30, 20-110, 100-300, 300-1000 
        if PLOTTYPE.lower() == 'pgplot':
            # Use PGPLOT to plot
            plot_sp_rated_pgplot(groups, ranks, inffile, 0, 30)
            print_debug("Finished PGplotting DMs0-30 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_pgplot(groups, ranks, inffile, 20, 110)
            print_debug("Finished PGplotting DMs20-110 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_pgplot(groups, ranks, inffile, 100, 310)
            print_debug("Finished PGplotting DMs100-310 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_pgplot(groups, ranks, inffile, 300, 1000)
            print_debug("Finished PGplotting DMs100-310 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_pgplot(groups, ranks, inffile, 1000, 10000)
            print_debug("Finished PGplotting DMs100-310 "+strftime("%Y-%m-%d %H:%M:%S"))
        elif PLOTTYPE.lower() == 'matplotlib':
            # Use matplotlib to plot
            plot_sp_rated_all(groups, ranks, inffile, 0, 30)
            print_debug("Finished plotting DMs0-30 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_all(groups, ranks, inffile, 20, 110)
            print_debug("Finished plotting DMs20-110 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_all(groups, ranks, inffile, 100, 310)
            print_debug("Finished plotting DMs100-310 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_all(groups, ranks, inffile, 300, 1000)
            print_debug("Finished plotting DMs300-1000 "+strftime("%Y-%m-%d %H:%M:%S"))
            plot_sp_rated_all(groups, ranks, inffile, 1000, 10000)
            print_debug("Finished plotting DMs1000-10000 "+strftime("%Y-%m-%d %H:%M:%S"))
        else:
            print("Plot type must be one of 'matplotlib' or 'pgplot'. Not plotting.")


if __name__ == '__main__':
    main()
