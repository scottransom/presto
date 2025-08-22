#!/usr/bin/env python
import sys
import copy
import numpy as np
import presto.Pgplot as pP
import presto.ppgplot as ppg
import presto.prepfold as ppfd

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

def usage():
    print(
        """
usage:  pfdzap.py pfd_file

  Interactively zap a .pfd file, and generate appropriate lines to use 
  for show_pfd, sum_profiles.py, and get_TOAs.py.

  Interactive keys:
    h or H               :  Print this help.
    a or A or left-mouse :  Kill interval or subband (or re-add it if killed)
    q or Q               :  Quit
    m or M               :  Switch frequency and time modes
    f or F               :  Switch to frequency (i.e. subband) mode
    t or T               :  Switch to time (i.e. subint) mode
    r or R               :  Zap a range rather than a single slice (via follow-on click)
    s or S               :  Show the current lists of subbands and intervals to kill
    c or C               :  Clear the kill lists and reset the profiles
    w or W               :  Write out get_TOAs.py or sum_profiles.py lines

"""
    )


if __name__ == "__main__":
    usage()
    if len(sys.argv)==1:
        sys.exit()
    psr = ppfd.pfd(sys.argv[1])
    print(f"Successfully read {sys.argv[1]}.\n")
    psr.dedisperse()
    profsbak = copy.copy(psr.profs)
    pP.ppgplot_dev_open_ = True
    ppg.pgopen("/xw")
    ppg.pgpap(6.0, 1.6)
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
        if key in [b"h", b"H"]:
            usage()
            continue
        if key in [b"q", b"Q"]:
            pP.closeplot()
            break
        if key in [b"c", b"C"]:
            psr.profs = copy.copy(profsbak)
            killsubs = set([])
            killints = set([])
            psr.kill_subbands(killsubs)
            psr.kill_intervals(killints)
            continue
        if key in [b"r", b"R"]:
            start = int(y)
            if mode == "freq":
                if 0 <= start < psr.nsub:
                    x, y, key = ppg.pgband(3, 0, x, y)
                    end = int(y)
                    if 0 <= end < psr.nsub:
                        if start < end:
                            killsubs = killsubs.union([i for i in range(start, end+1)])
                        else:
                            killsubs = killsubs.union([i for i in range(end, start+1)])
                        psr.kill_subbands(killsubs)
            else:
                if 0 <= start < psr.npart:
                    x, y, key = ppg.pgband(3, 0, x, y)
                    end = int(y)
                    if 0 <= end < psr.npart:
                        if start < end:
                            killints = killints.union([i for i in range(start, end+1)])
                        else:
                            killints = killints.union([i for i in range(end, start+1)])
                        psr.kill_intervals(killints)
            continue
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
                    psr.kill_subbands(killsubs)
            else:
                if 0 <= ii < psr.npart:
                    if ii not in killints:
                        killints.add(ii)
                    else:  # Remove the selection from the kill list
                        killints.remove(ii)
                        psr.profs[ii, :, :] = profsbak[ii, :, :]
                    psr.kill_intervals(killints)
            print(
                f"SNR ~ {calc_simple_SNR(psr):.3f}  reduced-chi^2 = {psr.calc_redchi2():.3f}"
            )
            continue
        if key in [b"s", b"S"]:
            print(f"killsubs = {sorted(killsubs)}")
            print(f"killints = {sorted(killints)}")
            continue
        if key in [b"w", b"W"]:
            subs = ",".join([str(i) for i in sorted(killsubs)]) if len(killsubs) else "X"
            ints = ",".join([str(i) for i in sorted(killints)]) if len(killints) else "X"
            print("#-----------------------")
            print("# For show_pfd:")
            print("#-----------------------")
            print(f"show_pfd -killsubs {subs} -killparts {ints} {sys.argv[1]}")
            print("#-----------------------")
            print("# For sum_profiles.py:")
            print("#-----------------------")
            print(f"{sys.argv[1]} {subs} {ints}")
            print("#-----------------------")
            print("# For get_TOAs.py:")
            print("#-----------------------")
            print(f"get_TOAs.py -n 1 -g TEMPLATE -2 -p -k {subs} -i {ints} {sys.argv[1]}")
            print("#-----------------------")
            continue

    sys.exit()
