import numpy as Num
from presto import candidate_sigma
import sifting, sys, re

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

class candidate:
    def __init__(self, candnum, sigma, numharm, ipow, cpow, bin, z,
                 DMstr, filename, T):
        self.filename = filename
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
    def __str__(self):
        cand = self.filename + ':' + `self.candnum`
        return "%-65s   %7.2f  %6.2f  %6.2f  %s   %7.1f  %7.1f  %12.6f  %10.2f  %8.2f "%\
               (cand, self.DM, self.snr, self.sigma, "%2d".center(7)%self.numharm,
                self.ipow_det, self.cpow, self.p*1000.0, self.r, self.z)
    def harms_to_snr(self):
        # Remove the average power level
        harmamps = Num.asarray(self.harm_pows) - 1.0
        # Set the S/N to 0.0 for harmonics with "negative" amplitudes
        harmamps[harmamps < 0.0] = 0.0
        self.snr = Num.sum(Num.sqrt(harmamps))

class file_candidates:
    def __init__(self, filename):
        self.filename = filename
        self.cands = {}
        self.badcands_knownbirds = []
        self.badcands_longperiod = []
        self.badcands_shortperiod = []
        self.badcands_threshold = []
        self.badcands_harmpowcutoff = []
        self.badcands_rogueharmpow = []
        current_goodcandnum = 0

        # First identify the length of the observation searched
        for line in open(filename, 'r'):
            if line.startswith(" Number of bins in the time series"):
                self.N = int(line.split()[-1])
            if line.startswith(" Width of each time series bin (sec)"):
                self.dt = float(line.split()[-1])
        self.T = self.N * self.dt

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
                f = bin / self.T  # Spin freq in hz
                p = 1.0 / f       # Spin period in sec

                # Reject very long period candidates
                if (p > long_period):
                    self.badcands_longperiod.append(candnum)
                    continue

                # Reject very short period candidates
                if (p < short_period):
                    self.badcands_shortperiod.append(candnum)
                    continue

                # Check to see if the candidate is in the known birds list
                known_bird = 0
                for bird, err in known_birds_f:
                    if (Num.fabs(f-bird) < err):
                        known_bird = 1
                        break
                if known_bird:
                    self.badcands_knownbirds.append(candnum)
                    continue
                for bird, err in known_birds_p:
                    if (Num.fabs(p*1000.0-bird) < err):
                        known_bird = 1
                        break
                if known_bird:
                    self.badcands_knownbirds.append(candnum)
                    continue

                # Add it to the candidates list
                DMstr = DM_re.search(filename).groups()[0]
                self.cands[candnum]=candidate(candnum, sigma, numharm,
                                              i_pow_det, c_pow, bin, z, 
                                              DMstr, filename, self.T)
                continue

            # Parse the harmonic powers
            elif harms_re.match(line):

                split_line = line.split()
                candnum = int(split_line[0])

                # Only read the harmonics for the candidates that weren't
                # rejected in the initial pass
                if self.cands.has_key(candnum):
                    self.cands[candnum].harm_pows = Num.zeros(self.cands[candnum].numharm, dtype=Num.float64)
                    self.cands[candnum].harm_amps = Num.zeros(self.cands[candnum].numharm, dtype=Num.complex64) 
                    power = parse_power(split_line[3])
                    phase = float(split_line[9].split("(")[0])
                    self.cands[candnum].harm_pows[0] = power
                    self.cands[candnum].harm_amps[0] = Num.sqrt(power) * Num.exp(phase*1.0j)
                    if (self.cands[candnum].numharm > 1):
                        current_goodcandnum = candnum
                        current_harmnum = 1
                    else:
                        current_goodcandnum = 0
                        # Compute the S/N
                        self.cands[candnum].harms_to_snr()
                        # These are the "optimized" power...
                        opt_ipow = self.cands[candnum].harm_pows[0]
                        # and sigma (calculated assuming _1_ trial!)
                        opt_sigma = candidate_sigma(opt_ipow, 1, 1)
                        self.cands[candnum].sigma = opt_sigma
                        self.cands[candnum].ipow_det = opt_ipow
                        # Remove the single harmonic candidates that
                        # don't pass our threshold
                        if (opt_sigma < sigma_threshold and
                            self.cands[candnum].cpow < c_pow_threshold):
                            self.badcands_threshold.append(candnum)
                            del(self.cands[candnum])
                continue

            # Parse the higher (than the first) harmonic powers
            if current_goodcandnum:
                cand = self.cands[current_goodcandnum]
                power = parse_power(line.split()[2])
                phase = float(line.split()[8].split("(")[0])
                cand.harm_pows[current_harmnum] = power
                cand.harm_amps[current_harmnum] = Num.sqrt(power) * Num.exp(phase*1.0j)
                
                current_harmnum += 1

                # Do the more advanced sifting after all the harmonics
                # have been read in
                if (current_harmnum==cand.numharm):

                    # Compute the S/N
                    cand.harms_to_snr()

                    # Insure that the sum of the optimized powers is > threshold
                    opt_ipow = sum(cand.harm_pows)
                    # Try to correct for the fact that by optimizing each
                    # harmonic power, we get a power that is slightly higher
                    # than it should be.  Simulations suggest that the average
                    # increase in power is ~2 per hamonic.  For single harmonics,
                    # though, the optimized power should be correct.  So the
                    # correction should be approx -2*(cand.numharm-1)
                    opt_sigma = candidate_sigma(opt_ipow-2.0*(cand.numharm-1),
                                                cand.numharm, 1)
                    self.cands[current_goodcandnum].sigma = opt_sigma
                    self.cands[current_goodcandnum].ipow_det = opt_ipow
                    if (opt_sigma < sigma_threshold):
                        self.badcands_threshold.append(current_goodcandnum)
                        del(self.cands[current_goodcandnum])
                        current_goodcandnum = 0
                        continue

                    # Remove the candidates where the harmonic with the
                    # highest power is not more than harm_pow_cutoff
                    maxharm = Num.argmax(cand.harm_pows)
                    maxpow = cand.harm_pows[maxharm]
                    if maxpow < harm_pow_cutoff:
                        self.badcands_harmpowcutoff.append(current_goodcandnum)
                        del(self.cands[current_goodcandnum])
                        current_goodcandnum = 0
                        continue

                    # Sort the harmonics by power
                    sortedpows = Num.sort(cand.harm_pows)

                    # Remove cands which are dominated by a single high-power
                    # but high-numbered harmonic
                    if (cand.numharm >= 8 and
                        maxharm > 4 and
                        maxpow > 2*sortedpows[-2]):
                        self.badcands_rogueharmpow.append(current_goodcandnum)
                        del(self.cands[current_goodcandnum])
                        current_goodcandnum = 0
                        continue
                    elif (cand.numharm >= 4 and
                          maxharm > 2 and
                          maxpow > 3*sortedpows[-2]):
                        self.badcands_rogueharmpow.append(current_goodcandnum)
                        del(self.cands[current_goodcandnum])
                        current_goodcandnum = 0
                        continue
                    current_goodcandnum = 0
                continue

    def print_cand_summary(self):
        print "   %s contains %d 'good' candidates" % \
              (self.filename, len(self.cands))
        print "      Known RFI rejects    (%d):" % \
              (len(self.badcands_knownbirds)), self.badcands_knownbirds
        print "      Short period rejects (%d):" % \
              (len(self.badcands_shortperiod)), self.badcands_shortperiod
        print "      Long period rejects  (%d):" % \
              (len(self.badcands_longperiod)), self.badcands_longperiod
        print "      Missed threhold      (%d):" % \
              (len(self.badcands_threshold)), self.badcands_threshold
        print "      No good harmonics    (%d):" % \
              (len(self.badcands_harmpowcutoff)), self.badcands_harmpowcutoff
        print "      One bad harmonic     (%d):" % \
              (len(self.badcands_rogueharmpow)), self.badcands_rogueharmpow

def read_candidates(filenms):
    """
    read_candidates(filenms):
        Read in accelsearch candidates from the test ACCEL files.
        Return a list of file_candidates instances.
    """
    if not len(filenms):
        print "Error:  There are no candidate files to read!"
        return None
    print "\nReading candidates...."
    candlist = []
    for filenm in filenms:
        filecands = file_candidates(filenm)
        #filecands.print_cand_summary()
        candlist += filecands.cands.values()
    return candlist

def remove_duplicate_candidates(candlist):
    """
    remove_duplicate_candidates(candlist):
        Remove lower-significance 'duplicate' (i.e. same period)
        candidates from a list of candidates.  For the highest
        significance candidate, include a list of the DMs (and SNRs)
        of all the other detections.
    """
    n = len(candlist)
    print "  Sorting the %d candidates by frequency..." % n
    candlist.sort(cmp_freq)
    print "  Searching for dupes..."
    ii = 0
    # Find any match
    while ii < n:
        jj = ii + 1
        if jj < n and Num.fabs(candlist[ii].r-candlist[jj].r) < r_err:
            # Find the rest that match
            jj += 1
            while jj < n and Num.fabs(candlist[ii].r-candlist[jj].r) < r_err:
                jj += 1
            matches = candlist[ii:jj]
            matches.sort(cmp_sigma)
            # flag the duplicates
            bestindex = candlist.index(matches[0])
            candlist[bestindex].hits = []
            for match in matches:
                candlist[bestindex].hits.append((match.DM, match.snr))
        else:
            candlist[ii].hits = [(candlist[ii].DM, candlist[ii].snr)]
        ii = jj
    # Remove the duplicate candidates (reverse order is necessary for indexing!)
    print "  Removing duplicates..."
    for ii in range(len(candlist)-1, -1, -1):
        if not hasattr(candlist[ii], 'hits'):
            del(candlist[ii])
    print "Found %d candidates.  Sorting them by significance...\n" % len(candlist)
    candlist.sort(cmp_sigma)
    return candlist

def remove_harmonics(candlist):
    """
    remove_harmonics(candlist):
        Remove the candidates that are lower significance harmonics
        of other candidates from the candlist.  Return a new candlist.
    """
    # Note:  should probably put the harmonics into the fundamental as hits (use sets)
    numcands = 0
    candlist.sort(cmp_sigma)
    f_err = r_err/candlist[0].T
    print "\nSearching for duplicate harmonics..."
    ii = 0
    while 1:
        jj = len(candlist) - 1
        zapj = 0
        while 1:
            if zapj:  print "Hey!"
            for factor in Num.arange(1.0, 17.0):
                if (Num.fabs(candlist[ii].f - candlist[jj].f*factor) < f_err or
                    Num.fabs(candlist[ii].f - candlist[jj].f/factor) < f_err):
                    #print "  Removing ", candlist[jj].filename, candlist[jj].candnum, \
                    #      candlist[jj].f, "because of", candlist[ii].filename, \
                    #      candlist[ii].candnum, candlist[ii].f, \
                    #      candlist[ii].f/candlist[jj].f, candlist[jj].f/candlist[ii].f
                    zapj = 1
                    break
            # Check a few other common ratios
            for factor in [3.0/2.0, 5.0/2.0,
                           2.0/3.0, 4.0/3.0, 5.0/3.0,
                           3.0/4.0, 5.0/4.0,
                           2.0/5.0, 3.0/5.0, 4.0/5.0]:
                if Num.fabs(candlist[ii].f-candlist[jj].f*factor) < f_err:
                    #print "  Removing ", candlist[jj].filename, candlist[jj].candnum, \
                    #      candlist[jj].f, "because of", candlist[ii].filename, \
                    #      candlist[ii].candnum, candlist[ii].f, \
                    #      candlist[ii].f/candlist[jj].f, candlist[jj].f/candlist[ii].f
                    zapj = 1
                    break
            if zapj:
                numcands += 1
                del(candlist[jj])
                zapj = 0
            jj -= 1
            if jj == ii:
                break
        ii += 1
        if ii >= len(candlist) - 1:
            break
    print "  Removed a total of %d harmonics.\n" % numcands
    return candlist

def remove_DM_problems(candlist, numdms, dmlist, low_DM_cutoff):
    """
    remove_DM_problems(candlist, numdms, dmlist, low_DM_cutoff):
        Remove the candidates where any of the following are true:
            1) The number of hits is < numdms
            2) The highest S/N candidate occurs below a DM of low_DM_cutoff
            3) The minimum difference in DM indices between the hits is > 1
    """
    # Create a dictionary where the key is the dmstr and the values are the index
    dmdict = {}
    for ii in range(len(dmlist)):
        dmdict[dmlist[ii]] = ii
    numcands = 0
    candlist.sort(cmp_sigma)
    for ii in range(len(candlist)-1, -1, -1):
        # Remove all the candidates without enough DM hits
        if len(candlist[ii].hits) < numdms:
            numcands += 1
            del(candlist[ii])
            continue
        # Remove all the candidates where the max SNR DM is less than the cutoff DM
        if float(candlist[ii].hits[0][0]) <= low_DM_cutoff:
            numcands += 1
            del(candlist[ii])
            continue
        # Remove all the candidates where there are no hits at consecutive DMs
        if len(candlist[ii].hits) > 1:
            candlist[ii].hits.sort(cmp_dms)
            dm_indices = Num.asarray([dmdict["%.2f"%candlist[ii].hits[jj][0]]
                                      for jj in range(len(candlist[ii].hits))])
            min_dmind_diff = min(dm_indices[1:] - dm_indices[:-1])
            if min_dmind_diff > 1:
                numcands += 1
                del(candlist[ii])
                continue

    print "Removed a total of %d candidates with DM problems.\n" % numcands
    return candlist

def write_candlist(candlist, candfilenm=None):
    if candfilenm is None:
        candfile = sys.stdout
    else:
        candfile = open(candfilenm, "w")
    candfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                   "SNR".center(8) + "sigma".center(8) + "numharm".center(9) +
                   "ipow".center(9) + "cpow".center(9) +  "P(ms)".center(14) +
                   "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
    for goodcand in candlist:
        candfile.write("%s (%d)\n" % (str(goodcand), len(goodcand.hits)))
        if (len(goodcand.hits) > 1):
            goodcand.hits.sort(cmp_dms)
            for hit in goodcand.hits:
                numstars = int(hit[1]/3.0)
                candfile.write("  DM=%6.2f SNR=%5.2f   "%hit + numstars*'*' + '\n')
    if candfilenm is not None:
        candfile.close()

