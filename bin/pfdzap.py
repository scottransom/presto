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


def print_stats(psr):
    print(
        f"SNR ~ {calc_simple_SNR(psr):6.2f}  reduced-chi^2 = {psr.calc_redchi2():6.2f}"
    )


def usage():
    print(
        """
usage:  pfdzap.py pfd_file

  Interactively zap a .pfd file, and generate appropriate lines to use 
  for show_pfd, sum_profiles.py, and get_TOAs.py.

  Interactive keys:
    h or H                :  Print this help.
    a or A or left-mouse  :  Kill interval or subband (or re-add it if killed)
    q or Q                :  Quit
    m or M                :  Switch between frequency (subband) and time (subint) modes
    r or R or right-mouse :  Zap a range (determined via follow-on left-mouse click)
    s or S                :  Show the current lists of subbands and intervals to kill
    c or C                :  Clear the kill lists and reset the profiles
    w or W                :  Write out get_TOAs.py or sum_profiles.py lines

"""
    )


if __name__ == "__main__":
    usage()
    if len(sys.argv) == 1:
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

    while True:
        print_stats(psr)
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
            psr.kill_subbands(killsubs)  # These re-compute
            psr.kill_intervals(killints)  # the profile stats
            continue
        if key in [b"r", b"R", b"x", b"X"]:
            if 0 <= y < (psr.nsub if mode == "freq" else psr.npart):
                x, y2, key = ppg.pgband(3, 0, x, y)
                if 0 <= y2 < (psr.nsub if mode == "freq" else psr.npart):
                    start, end = int(min(y, y2)), int(max(y, y2))
                    if mode == "freq":
                        killsubs = killsubs.union([i for i in range(start, end + 1)])
                        psr.kill_subbands(killsubs)
                    else:
                        killints = killints.union([i for i in range(start, end + 1)])
                        psr.kill_intervals(killints)
            continue
        if key in [b"m", b"M"]:
            mode = "freq" if mode == "time" else "time"
            print(
                f"Switching to sub{'int' if mode == 'time' else 'band'} ({mode}) mode."
            )
            continue
        if key in [b"a", b"A"]:
            ii = int(y)
            if 0 <= ii < (psr.nsub if mode == "freq" else psr.npart):
                kill = killsubs if mode == "freq" else killints
                if ii not in kill:
                    kill.add(ii)
                else:  # Remove the selection from the kill list
                    kill.remove(ii)
                    if mode == "freq":
                        psr.profs[:, ii, :] = profsbak[:, ii, :]
                    else:
                        psr.profs[ii, :, :] = profsbak[ii, :, :]
                psr.kill_subbands(killsubs)  # We have to re-kill both since we
                psr.kill_intervals(killints)  # potentially added some data back
            continue
        if key in [b"s", b"S"]:
            print(f"killsubs = {sorted(killsubs)}")
            print(f"killints = {sorted(killints)}")
            continue
        if key in [b"w", b"W"]:
            subs = (
                ",".join([str(i) for i in sorted(killsubs)]) if len(killsubs) else "X"
            )
            ints = (
                ",".join([str(i) for i in sorted(killints)]) if len(killints) else "X"
            )
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
            print(
                f"get_TOAs.py -n 1 -g TEMPLATE -2 -p -k {subs} -i {ints} {sys.argv[1]}"
            )
            print("#-----------------------")
            continue

    sys.exit()
