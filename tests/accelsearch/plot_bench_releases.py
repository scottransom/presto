#!/usr/bin/env python
"""Plot accelsearch speed across PRESTO releases from bench_releases.sh's results.tsv.

Usage:  plot_bench_releases.py RESULTS_TSV [OUTFILE]   (default OUTFILE: beside it,
        accel_speed_releases.png; a .pdf or .svg name works too)

The top panel is the speedup over the first version listed (the fastest of the runs,
at the same thread count); the bottom panel is the fastest wall time.  Each case is a
cluster of narrow bars, one per thread count.  Variant cases are hatched clusters next
to their version's: "<version>-median" (the new version with block-median
normalization) separates the code speedups from the normalization change, and
"<version>-native" (built with -march=native) shows what users who built for their CPU
got.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def read_best(path: Path) -> tuple[list[str], dict[tuple[str, int], float]]:
    """Read results.tsv; return the cases in file order and the best wall time per (case, ncpus).

    Parameters
    ----------
    path : Path
        The results.tsv written by bench_releases.sh.

    Returns
    -------
    cases : list of str
        Case labels in the order they first appear.
    best : dict
        The minimum wall time (s) for each (case, ncpus).
    """
    cases: list[str] = []
    best: dict[tuple[str, int], float] = defaultdict(lambda: float("inf"))
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            case = row["case"]
            if case not in cases:
                cases.append(case)
            key = (case, int(row["ncpus"]))
            best[key] = min(best[key], float(row["wall_s"]))
    return cases, dict(best)


# Variant cases ("<version>-<kind>"): bar hatch and legend text
VARIANT_HATCHES = {
    "median": ("///", "block-median normalization"),
    "native": ("...", "built with -march=native"),
}

BAR_W = 0.22  # width of one bar
SUB_GAP = 0.10  # gap between a version's cases
VER_GAP = 0.45  # gap between versions


def version_key(label: str) -> tuple[int, ...]:
    """Sort key for 'vX.Y.Z' labels (any '-suffix' ignored)."""
    return tuple(int(x) for x in label.lstrip("v").split("-")[0].split("."))


def main() -> None:
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    results = Path(sys.argv[1])
    out = (
        Path(sys.argv[2])
        if len(sys.argv) == 3
        else results.with_name("accel_speed_releases.png")
    )
    cases, best = read_best(results)
    variants = [c for c in cases if "-" in c]
    versions = sorted((c for c in cases if "-" not in c), key=version_key)
    ncpus = sorted({n for _, n in best})
    base = versions[0]
    colors = {n: f"C{i}" for i, n in enumerate(ncpus)}

    # Lay out one cluster of len(ncpus) bars per case, the version's own build first
    # and its variants beside it.
    cluster_w = BAR_W * len(ncpus)
    left: dict[str, float] = {}  # left edge of each case's cluster
    centers: list[float] = []  # tick position for each version
    x = 0.0
    for v in versions:
        group = [v] + [c for c in variants if c.split("-", 1)[0] == v]
        x0 = x
        for c in group:
            left[c] = x
            x += cluster_w + SUB_GAP
        x -= SUB_GAP
        centers.append(0.5 * (x0 + x))
        x += VER_GAP

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for c, x0 in left.items():
        kind = c.split("-", 1)[1] if "-" in c else None
        hatch = VARIANT_HATCHES.get(kind, ("xx", ""))[0] if kind else None
        for i, n in enumerate(ncpus):
            if (c, n) not in best:
                continue
            t = best[(c, n)]
            speed = best[(base, n)] / t if (base, n) in best else float("nan")
            xc = x0 + (i + 0.5) * BAR_W
            kw = dict(width=BAR_W, edgecolor="white", linewidth=0.5)
            if kind:
                kw.update(color="none", edgecolor=colors[n], hatch=hatch, linewidth=1.0)
            else:
                kw.update(color=colors[n])
            ax1.bar(xc, speed, **kw)
            ax2.bar(xc, t, **kw)

    handles = [
        Patch(color=colors[n], label=f"{n} thread{'s' if n > 1 else ''}") for n in ncpus
    ]
    for kind in sorted({c.split("-", 1)[1] for c in variants}):
        hatch, what = VARIANT_HATCHES.get(kind, ("xx", kind))
        names = ", ".join(c.split("-")[0] for c in variants if c.endswith("-" + kind))
        handles.append(
            Patch(
                facecolor="none", edgecolor="0.4", hatch=hatch, label=f"{names}: {what}"
            )
        )
    ax1.axhline(1.0, color="0.5", lw=0.8, zorder=0)
    ax1.set_ylabel(f"speedup over {base}")
    ax1.legend(handles=handles, frameon=False, loc="upper left")
    ax2.set_ylabel("wall time (s)")
    ax2.set_yscale("log")
    ax2.set_xticks(centers, versions, rotation=30)
    ax2.set_xlabel("PRESTO release")
    for ax in (ax1, ax2):
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
    ax1.set_title("accelsearch speed across releases")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
