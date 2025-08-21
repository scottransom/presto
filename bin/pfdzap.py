# %%
import sys
import copy
import numpy as np
from presto import presto
import presto.Pgplot as pP
import presto.prepfold as ppfd

ppg = pP.ppg  # shorthand


def plot_normprof(psr):
    pP.reset_colors()
    ppg.pgeras()
    pP.ppgplot_dev_prep_ = False
    normprof = psr.sumprof - min(psr.sumprof)
    normprof /= max(normprof)
    ppg.pgsvp(0.09, 0.92, 0.75, 0.99)
    ppg.pgswin(0.0, 1.0, 0.0, 1.0)
    ppg.pgline(np.linspace(0.0, 1.0, len(normprof)), normprof)
    ppg.pgsvp(0.09, 0.92, 0.06, 0.7)


def calc_simple_SNR(psr):
    prof = psr.sumprof - np.median(psr.sumprof)
    prof /= np.sqrt(psr.varprof)
    return prof.sum()


psr = ppfd.pfd(sys.argv[1])
psr.dedisperse()
profsbak = copy.copy(psr.profs)
pP.ppgplot_dev_open_ = True
ppg.pgopen("/xw")
ppg.pgpap(5.0, 1.6)
ppg.pgpage()  # clears panel
ppg.pgask(0)

killsubs = set([])
killints = set([])
mode = "freq"

print(f"SNR ~ {calc_simple_SNR(psr):.3f}  reduced-chi^2 = {psr.calc_redchi2():.3f}")
while True:
    plot_normprof(psr)
    if mode == "freq":
        psr.plot_subbands()
    else:
        psr.plot_intervals()
    x, y, key = ppg.pgband(7, 0)
    if key in [b"q", b"Q"]:
        pP.closeplot()
        break
    if key in [b"m", b"M"]:
        if mode == "freq":
            print("Switching to subint (time) mode.")
            mode = "time"
        else:
            print("Switching to subband (freq) mode.")
            mode = "freq"
        continue
    if key in [b"f", b"F"]:
        print("Switching to subband (freq) mode.")
        mode = "freq"
        continue
    if key in [b"t", b"T"]:
        print("Switching to subint (time) mode.")
        mode = "time"
        continue
    if key in [b"a", b"A"]:
        ii = int(y)
        if mode == "freq":
            if 0 <= ii < psr.nsub:
                if ii not in killsubs:
                    killsubs.add(ii)
                else:  # Remove the selection from the kill list
                    killsubs.remove(ii)
                    psr.profs[:, ii, :] = profsbak[:, ii, :]
                print(f"killsubs = {sorted(killsubs)}")
                psr.kill_subbands(killsubs)

        else:
            if 0 <= ii < psr.npart:
                if ii not in killints:
                    killints.add(ii)
                else:  # Remove the selection from the kill list
                    killints.remove(ii)
                    psr.profs[ii, :, :] = profsbak[ii, :, :]
                print(f"killints = {sorted(killints)}")
                psr.kill_intervals(killints)
        print(
            f"SNR ~ {calc_simple_SNR(psr):.3f}  reduced-chi^2 = {psr.calc_redchi2():.3f}"
        )
        continue
    if key in [b"x", b"X"]:
        # for this we will need to save the original profs and add them back
        pass

sys.exit()
