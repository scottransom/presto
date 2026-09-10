import glob
from operator import attrgetter

import presto.sifting as sifting

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

# --------------------------------------------------------------


def dms_from_filenames(filenames):
    """Return the sorted, unique DMs encoded in a list of PRESTO filenames.

    The DMs are pulled out with ``sifting.DM_re``, which finds a positive
    float with 2 or 3 digits after the decimal point following "_DM"
    *anywhere* in the name.  That handles both the usual "..._DM12.34.inf"
    naming and the "..._DM12.34_500M_2_ACCEL_0" naming of a short search,
    without caring where in the filename the DM happens to sit.  Names with
    no recognizable DM are simply ignored.

    Parameters
    ----------
    filenames : list of str
        The filenames to examine.

    Returns
    -------
    list of float
        The unique DMs found, in increasing order.
    """
    matches = (sifting.DM_re.search(filename) for filename in filenames)
    return sorted({float(match.group(1)) for match in matches if match is not None})


# Try to read the .inf files first, as _if_ they are present, all of
# them should be there.  (if no candidates are found by accelsearch
# we get no ACCEL files...
inffiles = glob.glob(globinf)
candfiles = glob.glob(globaccel)

# Get the list of DMs that were searched.  Fall back on the ACCEL
# filenames if the .inf files aren't around (or aren't named usefully).
dms = dms_from_filenames(inffiles) or dms_from_filenames(candfiles)
if not dms:
    raise SystemExit(
        f"Error:  could not determine any DMs from the '{globinf}' "
        f"or '{globaccel}' files!"
    )
dmstrs = [f"{dm:.2f}" for dm in dms]

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
    cands.sort(key=attrgetter("sigma"), reverse=True)
    sifting.write_candlist(cands)
