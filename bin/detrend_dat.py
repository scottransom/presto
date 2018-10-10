#!/usr/bin/env python
import os, sys
import numpy as np
from astropy.stats import sigma_clip
import scipy.signal

if len(sys.argv) != 2:
    print("\nusage: {} file\n".format(sys.argv[0]))
    sys.exit(1)

os.rename(sys.argv[1], sys.argv[1]+".bak")
data = np.fromfile(sys.argv[1]+".bak", dtype=np.float32)
N = len(data)

nblocks = 10000
data.shape = (nblocks, N/nblocks)
block_stds = np.std(data, axis=1)

good_stds = sigma_clip(block_stds, sigma=4.0)
stds_inds = np.arange(nblocks)[~good_stds.mask]

# zero-out the bad blocks
data[good_stds.mask,:] *= 0.0

print("Found %d bad blocks out of %d" % (good_stds.mask.sum(), nblocks))

# Now detrend the good blocks
for ii in stds_inds:
    data[ii] = scipy.signal.detrend(data[ii], type='linear')

data.ravel().tofile(sys.argv[1])

