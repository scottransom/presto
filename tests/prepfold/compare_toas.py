#!/usr/bin/env python3
"""Check that two .pfd files of the same data give the same TOAs.

Used by prepfold_tests.sh to compare a raw-data fold made with bary2topo()'s
cubic topocentric approximation against one made with tempo2 polycos
(prepfold -usepoly).  Both are timed with get_TOAs.py (with any extra
options given via --getoas-args, e.g. "-u" for folds that prepfold searched),
using the same template, and the TOAs must agree to within --tol microseconds.

Exits 0 if they agree, 1 if not, and 2 if get_TOAs.py fails.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys


def get_toas(
    pfdfile: str, template: str, numtoas: int, extra: list[str]
) -> list[tuple[int, float]]:
    """Return the TOA MJDs get_TOAs.py finds for a .pfd file.

    A double holds an MJD to only ~1 us, so each TOA comes back as its
    integer and fractional days.

    Parameters
    ----------
    pfdfile : str
        The .pfd file to time.
    template : str
        The template .bestprof file.
    numtoas : int
        The number of TOAs to make.
    extra : list of str
        Any other get_TOAs.py options.

    Returns
    -------
    list of (int, float)
        The integer and fractional parts of the TOA MJDs.
    """
    cmd = ["get_TOAs.py", "-2", "-n", str(numtoas), "-t", template, *extra, pfdfile]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode:
        sys.stderr.write(f"'{' '.join(cmd)}' failed:\n{out.stderr}")
        sys.exit(2)
    # Tempo2 format TOA lines:  name freq MJD err site ...
    toas = []
    for line in out.stdout.splitlines():
        if line.startswith(pfdfile):
            mjdi, mjdf = line.split()[2].split(".")
            toas.append((int(mjdi), float("0." + mjdf)))
    return toas


def main() -> int:
    """Compare the TOAs from two .pfd files and report the result."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pfd1")
    parser.add_argument("pfd2")
    parser.add_argument("-t", "--template", required=True, help="template .bestprof")
    parser.add_argument("-n", "--numtoas", type=int, default=4)
    parser.add_argument("--tol", type=float, default=1.0, help="tolerance (us)")
    parser.add_argument("--getoas-args", default="", help="extra get_TOAs.py options")
    args = parser.parse_args()
    extra = shlex.split(args.getoas_args)

    toas1 = get_toas(args.pfd1, args.template, args.numtoas, extra)
    toas2 = get_toas(args.pfd2, args.template, args.numtoas, extra)
    if len(toas1) != args.numtoas or len(toas2) != args.numtoas:
        print(f"  expected {args.numtoas} TOAs, got {len(toas1)} and {len(toas2)}")
        return 1
    diffs = [
        ((t2[0] - t1[0]) + (t2[1] - t1[1])) * 86400e6 for t1, t2 in zip(toas1, toas2)
    ]
    worst = max(abs(d) for d in diffs)
    print("  TOA differences (us): " + " ".join(f"{d:.3f}" for d in diffs))
    return 0 if worst < args.tol else 1


if __name__ == "__main__":
    sys.exit(main())
