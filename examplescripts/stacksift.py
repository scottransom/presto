import sys
import numpy as np
import matplotlib.pyplot as plt


class candidate:  # just an empty container for candidate parameters
    pass


def read_stackout(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    # Following might need to be modified based on what you use for filenames!
    dm = float(filename[:-4].split("_")[1][2:])
    candspart = False
    data = []
    for line in lines:
        if line.startswith("#"):
            if candspart is False:
                candspart = True
            continue
        if line.strip() and candspart:  # Skip empty lines
            x = candidate()
            x.sigma, x.f, x.p, x.bin, x.power, x.nharm = list(map(float, line.split()))
            x.dm = dm
            if x.f < 1000:  # Filter out candidates with frequency > 1000 Hz
                data.append(x)
    return dm, data


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python stacksift.py <stackout_files...>")
        sys.exit(1)

    allcands = []
    cands = {}
    for filename in sys.argv[1:]:
        dm, data = read_stackout(filename)
        cands[dm] = data
        allcands += data

    dms = np.asarray(list(cands.keys()))
    numfiles = len(dms)
    allcands.sort(key=lambda x: x.f)  # sort candidates by frequency
    print(f"Total candidates: {len(allcands)} from {numfiles} files")

    # Find candidates that are close to each other in frequency (within 1 bin) and group them
    bins = np.asarray([x.bin for x in allcands])
    group_barriers = np.arange(1, len(allcands))[np.diff(bins) > 1.0]
    groups = np.split(np.arange(len(allcands)), group_barriers)  # indices of groups

    # Now look at the candidates in each group and see if the highest sigma
    # is in the middle of the DM range of the group.  If so, make a plot
    # of sigma vs DM for that group.
    for group in groups:
        grp_dms = np.asarray([allcands[ii].dm for ii in group])
        grp_freqs = np.asarray([allcands[ii].f for ii in group])
        grp_sigmas = np.asarray([allcands[ii].sigma for ii in group])
        bins_in_group = np.asarray([allcands[ii].bin for ii in group])
        assert np.all(
            np.diff(bins_in_group) <= 1.0
        )  # sanity check that these are indeed close in frequency
        bestdm = grp_dms[np.argmax(grp_sigmas)]
        if bestdm > dms.min() and bestdm < dms.max() and len(group) <= numfiles:
            dmorder = np.argsort(grp_dms)
            minsig, maxsig = grp_sigmas.min(), grp_sigmas.max()
            print(
                f"Freq: {grp_freqs.mean():.12f} Hz  {1000.0 / grp_freqs.mean():.12f} ms:"
            )
            for dm, sigma in zip(grp_dms[dmorder], grp_sigmas[dmorder]):
                stars = "*" * max(1, int(sigma - minsig))
                print(f"  {dm:.2f}: {stars}")
