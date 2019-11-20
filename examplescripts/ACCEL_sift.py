from __future__ import absolute_import
from builtins import map
import re
import glob
import presto.sifting as sifting
from operator import itemgetter, attrgetter

# Note:  You will almost certainly want to adjust
#        the following variables for your particular search

# glob for ACCEL files
globaccel = "*ACCEL_*0"
# glob for .inf files
globinf = "*DM*.inf"
# In how many DMs must a candidate be detected to be considered "good"
min_num_DMs = 2
# Lowest DM to consider as a "real" pulsar
low_DM_cutoff = 2.0
# Ignore candidates with a sigma (from incoherent power summation) less than this
sifting.sigma_threshold = 4.0
# Ignore candidates with a coherent power less than this
sifting.c_pow_threshold = 100.0

# If the birds file works well, the following shouldn't
# be needed at all...  If they are, add tuples with the bad
# values and their errors.
#                (ms, err)
sifting.known_birds_p = []
#                (Hz, err)
sifting.known_birds_f = []

# The following are all defined in the sifting module.
# But if we want to override them, uncomment and do it here.
# You shouldn't need to adjust them for most searches, though.

# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
sifting.r_err = 1.1
# Shortest period candidates to consider (s)
sifting.short_period = 0.0005
# Longest period candidates to consider (s)
sifting.long_period = 15.0
# Ignore any candidates where at least one harmonic does exceed this power
sifting.harm_pow_cutoff = 8.0

#--------------------------------------------------------------

# Try to read the .inf files first, as _if_ they are present, all of
# them should be there.  (if no candidates are found by accelsearch
# we get no ACCEL files...
inffiles = glob.glob(globinf)
candfiles = glob.glob(globaccel)
# Check to see if this is from a short search
if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
    dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
else:
    dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
dms = list(map(float, dmstrs))
dms.sort()
dmstrs = ["%.2f"%x for x in dms]

# Read in all the candidates
cands = sifting.read_candidates(candfiles)

# Remove candidates that are duplicated in other ACCEL files
if len(cands):
    cands = sifting.remove_duplicate_candidates(cands)

# Remove candidates with DM problems
if len(cands):
    cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

# Remove candidates that are harmonically related to each other
# Note:  this includes only a small set of harmonics
if len(cands):
    cands = sifting.remove_harmonics(cands)

# Write candidates to STDOUT
if len(cands):
    cands.sort(key=attrgetter('sigma'), reverse=True)
    sifting.write_candlist(cands)
