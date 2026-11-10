"""End-to-end search of real data: GBT observations of the globular cluster Terzan 5.

The data are 859 s of 2-bit, 64-channel GBT PSRFITS at 2 GHz (`Ter5_080912_short2bits.fits`,
173 MB, also used by the prepfold suites in `tests/prepfold/`). They are not in the repo:
`tests/prepfold/get_test_data.sh` fetches them from the repo's `test-data-v1` GitHub release,
and these tests skip if the file is absent. `pixi run test-realdata` fetches the data and runs
this file, and CI runs that task weekly (`.github/workflows/realdata.yml`).

The pipeline is rfifind -> prepsubband (8 DMs, barycentered, masked) -> realfft -> rednoise ->
accelsearch -> sifting -> prepfold, plus single_pulse_search.py and coherentsearch (at DM
238, where it must label Ter5A, Ter5C and Ter5N and reject the combs that Ter5A makes at
slow trial frequencies). The truth comes from the
timing solutions in `tests/prepfold/`: `binary_psr` predicts each pulsar's Doppler-shifted
frequency and frequency derivative during the observation. Two pulsars are checked:

- Ter5A (11.56 ms, DM 242.3, 0.076-day eclipsing binary): very bright, with z ~ -5 bins of
  acceleration and detectable single pulses.
- Ter5N (8.66 ms, DM 238.3, 0.385-day binary): weak, sigma ~ 8 in the ACCEL file.

As in test_pipeline_synthetic.py, thresholds sit below the measured values (noted at each
assert), so that FFTW plan differences across platforms cannot flip a test.
"""

from __future__ import annotations

import re
import subprocess
from operator import attrgetter
from pathlib import Path

import numpy as np
import pytest

from presto import binary_psr, infodata, rfifind, sifting

PFDIR = Path(__file__).resolve().parent / "prepfold"
RAW = PFDIR / "Ter5_080912_short2bits.fits"
DMS = [234.0 + 2.0 * i for i in range(8)]  # 234 ... 248
ZMAX = 50  # Ter5A's 8th harmonic drifts by 8 x 5 = 40 bins
DT = 2 * 81.92e-6  # prepsubband -downsamp 2
NUMOUT = 5280000  # what prepsubband writes for this file at that sampling
T = NUMOUT * DT

pytestmark = pytest.mark.skipif(
    not RAW.exists(),
    reason=f"{RAW.name} not present: run tests/prepfold/get_test_data.sh",
)


def run(cmd: list[str], cwd: Path) -> str:
    """Run a PRESTO tool, failing the test with its output if it exits non-zero."""
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    assert res.returncode == 0, (
        f"{' '.join(cmd)} exited {res.returncode}\n--- stdout ---\n{res.stdout[-4000:]}"
        f"\n--- stderr ---\n{res.stderr[-4000:]}"
    )
    return res.stdout


def base(dm: float) -> str:
    return f"ter5_DM{dm:.2f}"


class Ephem:
    """Doppler-shifted spin of a binary pulsar during the observation, from its par file.

    Times are seconds from the start of the (barycentered) time series.  The spin-down
    over the observation is negligible here (< 1e-6 Hz) and ignored.
    """

    def __init__(self, name: str, epoch: float):
        self.psr = binary_psr.binary_psr(str(PFDIR / f"{name}.par"))
        self.epoch = epoch

    def freq(self, t: float | np.ndarray) -> np.ndarray:
        return 1.0 / self.psr.doppler_period(self.epoch + np.asarray(t) / 86400.0)

    def z(self) -> float:
        """Total frequency drift over the observation, in Fourier bins (accelsearch's z)."""
        return float((self.freq(T) - self.freq(0.0))[0] * T)

    def phase(self, t: np.ndarray) -> np.ndarray:
        """Pulse phase (in turns, arbitrary zero) at times t."""
        grid = np.linspace(0.0, T, 100001)
        f = self.freq(grid)
        ph = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(grid))])
        return np.interp(t, grid, ph)


@pytest.fixture(scope="module")
def obs(tmp_path_factory) -> dict:
    """Run the whole pipeline once in a temporary directory; return what the tests need."""
    d = tmp_path_factory.mktemp("ter5")
    (d / RAW.name).symlink_to(RAW)
    rfi_out = run(
        ["rfifind", "-noscales", "-nooffsets", "-time", "2.0", "-o", "ter5", RAW.name],
        d,
    )
    run(
        [
            "prepsubband",
            "-noscales",
            "-nooffsets",
            "-nsub",
            "32",
            "-lodm",
            f"{DMS[0]}",
            "-dmstep",
            "2",
            "-numdms",
            f"{len(DMS)}",
            "-downsamp",
            "2",
            "-mask",
            "ter5_rfifind.mask",
            "-o",
            "ter5",
            RAW.name,
        ],
        d,
    )
    for dm in DMS:
        run(["realfft", f"{base(dm)}.dat"], d)
        run(["rednoise", f"{base(dm)}.fft"], d)
        run(
            ["accelsearch", "-zmax", f"{ZMAX}", "-numharm", "8", f"{base(dm)}_red.fft"],
            d,
        )
    run(["single_pulse_search.py", "-b", *[f"{base(dm)}.dat" for dm in DMS]], d)
    epoch = infodata.infodata(str(d / f"{base(DMS[0])}.inf")).epoch
    return {
        "dir": d,
        "rfifind_stdout": rfi_out,
        "epoch": epoch,
        "A": Ephem("Ter5A", epoch),
        "N": Ephem("Ter5N", epoch),
    }


def accel_cands(obs: dict, dm: float) -> list[sifting.Candidate]:
    fn = obs["dir"] / f"{base(dm)}_red_ACCEL_{ZMAX}"
    return list(sifting.candlist_from_candfile(str(fn)).cands) if fn.exists() else []


def find(cands, eph: Ephem, tol_bins: float = 1.0):
    """The strongest candidate within tol_bins of the predicted mid-observation frequency."""
    fmid = eph.freq(T / 2)[0]
    near = [c for c in cands if abs(c.f - fmid) * T < tol_bins]
    return max(near, key=attrgetter("sigma")) if near else None


def fold_fmid(path: Path) -> tuple[float, dict[str, float]]:
    """Mid-fold barycentric spin frequency from a .bestprof, plus the parsed values.

    prepfold quotes P and P' at the start of the fold.  Over 13% of Ter5A's orbit a
    P, P', P'' fit cannot match the curved Doppler track at the quoted precision, so the
    comparisons use the mean (mid-fold) frequency instead.
    """
    vals: dict[str, float] = {}
    for line in path.read_text().splitlines():
        key, _, rest = line[2:].partition("=")
        key = key.strip()
        if key in (
            "P_bary (ms)",
            "P'_bary (s/s)",
            "T_sample",
            "Data Folded",
            "Best DM",
        ):
            vals[key] = float(rest.split()[0])
        elif line.startswith("# Prob(Noise)"):
            vals["sigma"] = float(line.split("~")[1].split()[0])
    tfold = vals["Data Folded"] * vals["T_sample"]
    pmid = 1e-3 * vals["P_bary (ms)"] + vals["P'_bary (s/s)"] * tfold / 2
    return 1.0 / pmid, vals


def test_rfifind(obs: dict) -> None:
    """rfifind writes a readable mask for the whole file, masking some but not much."""
    r = rfifind.rfifind(str(obs["dir"] / "ter5_rfifind.stats"))
    r.read_mask()
    assert r.nchan == 64
    assert r.nint * r.dtint == pytest.approx(859.0, abs=3.0)
    line = next(
        ln for ln in obs["rfifind_stdout"].splitlines() if "bad   intervals" in ln
    )
    badpct = float(line.split("(")[1].split("%")[0])
    assert 1.0 < badpct < 10.0  # measured 4.36%


def test_prepsubband_inf(obs: dict) -> None:
    """The DM trials are barycentered, downsampled, and carry their DMs."""
    for dm in DMS:
        inf = infodata.infodata(str(obs["dir"] / f"{base(dm)}.inf"))
        assert inf.N == NUMOUT
        assert inf.dt == pytest.approx(DT)
        assert inf.DM == pytest.approx(dm)
        assert inf.bary == 1
        assert inf.telescope == "GBT"
    assert obs["epoch"] == pytest.approx(54721.8543, abs=1e-3)


def test_accelsearch_ter5A(obs: dict) -> None:
    """Ter5A is found at its predicted frequency and acceleration, strongest near DM 242."""
    eph = obs["A"]
    sigmas = {}
    for dm in DMS:
        cand = find(accel_cands(obs, dm), eph)
        assert cand is not None, f"Ter5A not found at DM {dm}"
        sigmas[dm] = cand.sigma
        # measured: 0.08-0.15 bins off in frequency, |dz| <= 0.15
        assert abs(cand.f - eph.freq(T / 2)[0]) * T < 0.5
        assert abs(cand.z - eph.z()) < 1.5
    # measured sigma ~ 89 at DM 242, falling to ~44 at DM 234 and ~53 at DM 248
    assert max(sigmas, key=sigmas.get) in (240.0, 242.0, 244.0), sigmas
    assert sigmas[242.0] > 50.0
    assert sigmas[242.0] > sigmas[234.0] + 15.0 and sigmas[242.0] > sigmas[248.0] + 15.0


def test_accelsearch_ter5N(obs: dict) -> None:
    """The much weaker Ter5N is found at DM 238 at its predicted frequency and z."""
    eph = obs["N"]
    cand = find(accel_cands(obs, 238.0), eph)
    assert cand is not None, "Ter5N not found at DM 238"
    assert abs(cand.f - eph.freq(T / 2)[0]) * T < 0.5  # measured 0.03 bins
    assert abs(cand.z - eph.z()) < 1.5  # measured 0.05
    assert (
        cand.sigma > 6.0
    )  # measured 11.1 (single-trial; see test_pipeline_synthetic.py)


@pytest.fixture(scope="module")
def sifted(obs: dict) -> sifting.Candlist:
    files = [str(obs["dir"] / f"{base(dm)}_red_ACCEL_{ZMAX}") for dm in DMS]
    cands = sifting.read_candidates(files)
    cands = sifting.remove_duplicate_candidates(cands)
    cands = sifting.remove_DM_problems(cands, 2, [f"{dm:.2f}" for dm in DMS], 2.0)
    cands = sifting.remove_harmonics(cands)
    cands.sort(key=attrgetter("sigma"), reverse=True)
    return cands


def test_sifting(obs: dict, sifted: sifting.Candlist) -> None:
    """After sifting, Ter5A is the top candidate, at DM 242, and Ter5N survives."""
    top = sifted[0]
    assert abs(top.f - obs["A"].freq(T / 2)[0]) * T < 0.5
    assert top.DM == pytest.approx(242.0)
    assert find(sifted, obs["N"]) is not None


def test_prepfold_raw_ter5A(obs: dict) -> None:
    """Folding the raw data from Ter5A's accelsearch candidate recovers its DM and period."""
    d, eph = obs["dir"], obs["A"]
    cand = find(accel_cands(obs, 242.0), eph)
    run(
        [
            "prepfold",
            "-noxwin",
            "-noscales",
            "-nooffsets",
            "-mask",
            "ter5_rfifind.mask",
            "-nsub",
            "16",
            "-npfact",
            "1",
            "-ndmfact",
            "1",
            "-n",
            "64",
            "-accelcand",
            f"{cand.candnum}",
            "-accelfile",
            f"{base(242.0)}_red_ACCEL_{ZMAX}.cand",
            "-dm",
            "242",
            "-o",
            "foldA",
            RAW.name,
        ],
        d,
    )
    fmid, bp = fold_fmid(d / f"foldA_ACCEL_Cand_{cand.candnum}.pfd.bestprof")
    assert abs(bp["Best DM"] - 242.34) < 0.5  # measured 242.225 (timing DM 242.342)
    assert abs(fmid - eph.freq(T / 2)[0]) * T < 0.25  # measured 0.10 bins
    assert bp["sigma"] > 100.0  # measured 244


def test_prepfold_dat_ter5N(obs: dict) -> None:
    """Folding the DM-238 time series from Ter5N's candidate recovers its period."""
    d, eph = obs["dir"], obs["N"]
    cand = find(accel_cands(obs, 238.0), eph)
    run(
        [
            "prepfold",
            "-noxwin",
            "-accelcand",
            f"{cand.candnum}",
            "-accelfile",
            f"{base(238.0)}_red_ACCEL_{ZMAX}.cand",
            "-o",
            "foldN",
            f"{base(238.0)}.dat",
        ],
        d,
    )
    fmid, bp = fold_fmid(d / f"foldN_ACCEL_Cand_{cand.candnum}.pfd.bestprof")
    assert abs(fmid - eph.freq(T / 2)[0]) * T < 0.25  # measured 0.005 bins
    assert bp["sigma"] > 6.0  # measured 9.5


def rayleigh_z(obs: dict, dm: float, minsigma: float = 6.0) -> float:
    """Rayleigh Z^2_1 of single-pulse event times folded with Ter5A's ephemeris.

    For events unrelated to Ter5A, Z is exponentially distributed with mean 1.
    """
    ev = np.loadtxt(obs["dir"] / f"{base(dm)}.singlepulse", ndmin=2)
    t = ev[ev[:, 1] > minsigma, 2]
    ph = obs["A"].phase(t)
    return len(t) * abs(np.mean(np.exp(2j * np.pi * ph))) ** 2


def test_single_pulse_search_ter5A(obs: dict) -> None:
    """single_pulse_search.py finds Ter5A's bright single pulses at its DM.

    Above 6 sigma, the DM-242 events are phase-locked to Ter5A's ephemeris, and the events
    at DMs well away from it are not.  (The many events above ~8 sigma appear at every DM
    and are RFI, so they dilute the signal at higher thresholds.)  The data have periodic
    dropouts (~6% of 1000-sample blocks are constant, from 273 s on), so this also checks
    that single_pulse_search.py skips dead blocks: before that fix, macOS found ~8x more
    events (wide-boxcar noise from rounding residuals), and Z at DM 242 fell to 5-7.
    """
    for dm in DMS:
        assert (obs["dir"] / f"{base(dm)}.singlepulse").exists()
    assert (obs["dir"] / "ter5_singlepulse.ps").exists()
    zs = {dm: rayleigh_z(obs, dm) for dm in DMS}
    counts = {
        dm: len(np.loadtxt(obs["dir"] / f"{base(dm)}.singlepulse", ndmin=2))
        for dm in DMS
    }
    info = f"Z by DM: {zs}; events by DM: {counts}"
    assert zs[242.0] > 15.0, info  # measured 30.6
    for dm in (234.0, 238.0, 248.0):
        assert zs[dm] < 8.0, info  # measured 0.9, 1.0, 0.9 (2.2-4.2 at DM 240, 244)


CS_DM = 238.0  # between Ter5N (238.3) and Ter5C (237.1); Ter5A (242.3) is bright enough


@pytest.fixture(scope="module")
def cohout(obs: dict) -> dict:
    """coherentsearch on the DM-238 spectrum: its candidates and its stdout summary."""
    d = obs["dir"]
    res = subprocess.run(
        [
            "coherentsearch",
            "-noprogress",
            "-ncpus",
            "2",
            "-threshold",
            "7",
            f"{base(CS_DM)}_red.fft",
        ],
        cwd=d,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr[-4000:]
    cands = []
    for ln in (d / f"{base(CS_DM)}_red.cohout").read_text().splitlines():
        if ln.startswith("#"):
            continue
        w = ln.split(maxsplit=6)
        cands.append(
            {"snr": float(w[1]), "f": float(w[2]), "note": w[6] if len(w) > 6 else ""}
        )
    return {"cands": cands, "stdout": res.stdout}


def test_coherentsearch_known_pulsars(obs: dict, cohout: dict) -> None:
    """coherentsearch finds Ter5A, Ter5C and Ter5N at DM 238 and labels each from the
    pulsar database; A and N also sit at their ephemeris frequencies.

    coherentsearch searches no acceleration, so a candidate lies within the |z|/2 bins
    its pulsar drifts either side of the mid-observation frequency.
    """
    cands, info = cohout["cands"], f"candidates: {cohout['cands']}"

    def labelled(name: str) -> dict:
        hits = [c for c in cands if c["note"].startswith(f"PSR {name}")]
        assert hits, f"no candidate labelled PSR {name}; {info}"
        return max(hits, key=lambda c: c["snr"])

    a = labelled("B1744-24A")
    assert a["snr"] > 30.0, info  # measured 40.6
    n = labelled("J1748-2446N")
    assert n["snr"] > 7.5, info  # measured 9.6
    c = labelled("J1748-2446C")
    assert c["snr"] > 7.5, info  # measured 9.1
    for cand, eph in ((a, obs["A"]), (n, obs["N"])):
        tol = abs(eph.z()) / 2 + 1.0
        assert abs(cand["f"] - eph.freq(T / 2)[0]) * T < tol, info
        # measured 0.90 bins for A (tol 3.6), 0.01 for N (tol 1.3)


def test_coherentsearch_rejects_periodic(obs: dict, cohout: dict) -> None:
    """Ter5A's fundamental and 2nd harmonic land on high harmonics of many slow trials,
    which fold into combs of identical peaks (S/N ~10 at 1.6-5.1 Hz); coherentsearch
    drops them and says how many."""
    m = re.search(r"(\d+) removed as periodic", cohout["stdout"])
    assert m, cohout["stdout"]
    assert int(m.group(1)) >= 4, cohout["stdout"]  # measured 8
    # No survivor carries Ter5A (or its 2nd harmonic) on its own harmonic 7..60.
    fa = obs["A"].freq(T / 2)[0]
    for c in cohout["cands"]:
        for n in range(7, 61):
            for fh in (fa, 2 * fa):
                assert abs(n * c["f"] - fh) * T > 4.0, (n, c, cohout["cands"])
