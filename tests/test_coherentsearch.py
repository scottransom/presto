"""Pin the coherentsearch port's reference numbers (``tests/coherentsearch/``).

The fixtures there were written by CoherentSearch.jl itself (``gen_fixtures.jl``) on
``synth_cs.fft``.  These tests recompute them with ``cs_oracle.py``, a plain NumPy
statement of the same algorithm, so the specification the C port is written against
is checked against the Julia code it ports, stage by stage:

interpolated harmonic stacks -> which harmonics are used -> inverse-FFT profiles ->
width banks, noise scales and shape factors -> per-trial metrics (every decimation)
-> chunk boundaries of a full search -> the metric of every raw candidate it reports.

The C code is then checked against the same fixtures through ``cs_dump`` (a test
helper built by meson, not installed), which writes one chunk's stages in the fixture
layout.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent / "coherentsearch"
FIX = HERE / "fixtures"
sys.path.insert(0, str(HERE))

import cs_oracle as cs  # noqa: E402

REPO = Path(__file__).resolve().parent.parent

PARAMS = cs.Params()
PROV = tomllib.loads((FIX / "provenance.toml").read_text())
CASES = PROV["cases"]


@pytest.fixture(scope="module")
def fft() -> cs.FFTData:
    return cs.read_fft(HERE / PROV["input"])


def load_case(name: str) -> dict:
    return tomllib.loads((FIX / name / "meta.toml").read_text())


def load(name: str, fn: str, dtype: str, ncols: int | None = None) -> np.ndarray:
    a = np.fromfile(FIX / name / fn, dtype=dtype)
    return a if ncols is None else a.reshape(-1, ncols)


@pytest.fixture(scope="module")
def oracle(fft):
    """The oracle's run of every fixture chunk, computed once."""
    out = {}
    for name in CASES:
        meta = load_case(name)
        out[name] = cs.chunk_metrics(fft, PARAMS, meta["r_lo"], meta["t0"], meta["n"])
    return out


def test_provenance():
    """The fixtures were made with the parameters the oracle assumes."""
    assert PROV["nharms"] == PARAMS.nharms
    assert PROV["m"] == PARAMS.m
    assert PROV["hidr"] == PARAMS.hidr
    assert PROV["maxdecim"] == PARAMS.maxdecim
    assert PROV["sigma_cap"] == cs.SIGMA_CAP


@pytest.mark.parametrize("name", CASES)
def test_harmonic_stacks(name, oracle):
    """Filled harmonics and the interpolated amplitudes (Eqn. 30, m = 16)."""
    meta = load_case(name)
    ftprofs, filled, _ = oracle[name]
    assert list(np.flatnonzero(filled)) == meta["filled"]
    ref = load(name, "ftprofs.c64", "<c8", PARAMS.nharms + 1)
    assert ref.shape == ftprofs.shape
    assert np.all(ref[:, ~filled] == 0)
    # The fixture is the Float64-weight Julia path (~1e-10 from the exact kernel)
    # rounded to complex64, so this is complex64 rounding plus nothing (measured
    # 5.4e-8 worst over the four chunks).
    err = np.max(np.abs(ftprofs - ref)) / np.max(np.abs(ref))
    assert err < 2e-7, err


@pytest.mark.parametrize("name", CASES)
def test_profiles(name, oracle):
    """The unnormalised c2r inverse FFT of the k = 1 stacks (stored as float32)."""
    _, _, res = oracle[name]
    ref = load(name, "profs_k1.f32", "<f4", 2 * PARAMS.nharms)
    err = np.max(np.abs(res[1].profs - ref)) / np.max(np.abs(ref))
    assert err < 2e-7, err  # measured 4.5e-8: float32 storage rounding


@pytest.mark.parametrize("name", CASES)
def test_metric_stage(name, oracle):
    """Width banks, noise scales, shape factors and every per-trial metric."""
    meta = load_case(name)
    _, _, res = oracle[name]
    for k in PARAMS.decimations:
        km = meta[f"k{k}"]
        r = res[k]
        assert r.Hk == km["Hk"]
        assert r.widths == km["widths"]
        assert r.nvalid == km["nvalid"]
        assert r.sigma_analytic == pytest.approx(km["sigma_analytic"], rel=1e-14)
        # measured 2.5e-15 (same samples, same median, different FFT rounding)
        assert r.sigma_measured == pytest.approx(km["sigma_measured"], rel=1e-12)
        np.testing.assert_allclose(r.shape, km["shape"], rtol=1e-12)
        for kind in ("analytic", "measured"):
            ref = load(name, f"metric_k{k}_{kind}.f64", "<f8")
            got = getattr(r, f"metric_{kind}")
            assert len(ref) == r.nvalid
            # measured 3.0e-14 worst: the oracle and Julia's Float64 path agree to
            # FFT/summation-order rounding.
            np.testing.assert_allclose(
                got, ref, rtol=0, atol=1e-11, err_msg=f"k={k} {kind}"
            )


@pytest.mark.parametrize("name", CASES)
def test_prod_precision(name, oracle):
    """How far Julia's own single-precision production path is from the exact one.

    Float32 interpolation weights and a Float32 profile stage: the scale of
    difference a single-precision C port should also show.
    """
    worst = 0.0
    for k in PARAMS.decimations:
        exact = load(name, f"metric_k{k}_analytic.f64", "<f8")
        prod = load(name, f"metric_k{k}_analytic_prod.f64", "<f8")
        worst = max(worst, float(np.max(np.abs(prod - exact))))
    assert worst < 2e-5, worst  # measured 2.4e-6 (metric values up to ~18)


@pytest.fixture(scope="module")
def search_meta() -> dict:
    return tomllib.loads((FIX / "search" / "meta.toml").read_text())


def test_chunk_starts(fft, search_meta):
    """Chunk boundaries of the full default search: the regular grid plus a split
    wherever a harmonic starts or stops fitting.  Exact integers."""
    m = search_meta
    got = cs.chunk_starts(fft, PARAMS, m["r_lo"], m["total"], m["blocksize"])
    assert got == m["chunk_starts"]


def test_raw_candidates(fft, search_meta):
    """Recompute the metric of every raw (threshold 6) candidate of the full search.

    Each candidate's trial is located on the global grid, its chunk found from the
    chunk boundaries (which decide the filled harmonics, hence the analytic sigma),
    and its metric recomputed exactly.  The fixture values are Julia's production
    (single-precision) numbers.
    """
    m = search_meta
    r_lo, lodr = m["r_lo"], m["lodr"]
    starts = np.array(m["chunk_starts"] + [m["total"]])
    cands = np.loadtxt(FIX / "search" / "raw_cands_thresh6.txt", ndmin=2)
    assert len(cands) == m["n_raw_thresh6"]
    cache = {}
    worst = 0.0
    for r, _freq, metric, nharm in cands:
        k = PARAMS.nharms // int(nharm)
        assert PARAMS.nharms // k == nharm
        t = int(round((r / k - r_lo) / lodr))
        assert abs(r - k * (r_lo + t * lodr)) < 1e-9 * r
        c = int(np.searchsorted(starts, t, side="right")) - 1
        t0, t1 = int(starts[c]), int(starts[c + 1])
        if (c, k) not in cache:
            filled = cs.chunk_filled(fft, PARAMS, r_lo, t0, t1 - t0)
            Hk = PARAMS.nharms // k
            widths = cs.ladder_boxcar_widths(2 * Hk, k, PARAMS)
            cache[(c, k)] = (
                filled,
                widths,
                cs.analytic_sigma(filled, k, Hk),
                cs.boxcar_shape(filled, k, Hk, widths),
            )
        filled, widths, sigma, shape = cache[(c, k)]
        stack = cs.trial_stack(fft, PARAMS, r_lo, t, filled)
        prof = cs.profiles(stack[None, :], k, PARAMS.nharms)[0]
        got = cs.boxcar_metric(prof, widths, sigma, shape)
        worst = max(worst, abs(got - metric))
    print(f"raw candidates: {len(cands)}, worst |oracle - julia prod| = {worst:.2e}")
    # measured 2.5e-6, i.e. Julia's single-precision rounding
    assert worst < 2e-5, worst


# --- the C port, through cs_dump -----------------------------------------------------


def _cs_dump() -> Path:
    """The ``cs_dump`` helper: ``$CS_DUMP``, else the pixi build's."""
    if "CS_DUMP" in os.environ:
        return Path(os.environ["CS_DUMP"])
    build = REPO / "build-pixi"
    exe = build / "src" / "cs_dump"
    if not build.is_dir():
        pytest.skip("no build-pixi/ (run `pixi run build`) and $CS_DUMP is not set")
    # A pixi build that lacks it is a build problem, not a reason to skip.
    assert exe.is_file(), f"{exe} missing: rebuild with `pixi run build`"
    return exe


@pytest.fixture(scope="module")
def cdump(tmp_path_factory) -> Path:
    """Run ``cs_dump`` on every fixture chunk; return the directory of outputs."""
    exe = _cs_dump()
    out = tmp_path_factory.mktemp("cs_dump")
    for name in CASES:
        meta = load_case(name)
        lofreq = meta["r_lo"] / cs.read_fft(HERE / PROV["input"]).T
        subprocess.run(
            [
                str(exe),
                str(HERE / PROV["input"]),
                str(out / name),
                str(meta["t0"]),
                str(meta["n"]),
                repr(lofreq),
                str(PROV["maxdecim"]),
            ],
            check=True,
        )
    return out


def cload(
    cdir: Path, name: str, fn: str, dtype: str, ncols: int | None = None
) -> np.ndarray:
    a = np.fromfile(cdir / name / fn, dtype=dtype)
    return a if ncols is None else a.reshape(-1, ncols)


@pytest.mark.parametrize("name", CASES)
def test_c_harmonic_stacks(name, cdump):
    """C interpolation (float32 weights and sums) against the exact Julia stacks."""
    meta = load_case(name)
    cmeta = tomllib.loads((cdump / name / "meta.toml").read_text())
    assert cmeta["filled"] == meta["filled"]
    ref = load(name, "ftprofs.c64", "<c8", PARAMS.nharms + 1)
    got = cload(cdump, name, "ftprofs.c64", "<c8", PARAMS.nharms + 1)
    assert got.shape == ref.shape
    err = np.max(np.abs(got - ref)) / np.max(np.abs(ref))
    assert err < 1e-6, err  # measured 2.8e-7 worst


@pytest.mark.parametrize("name", CASES)
def test_c_profiles(name, cdump):
    ref = load(name, "profs_k1.f32", "<f4", 2 * PARAMS.nharms)
    got = cload(cdump, name, "profs_k1.f32", "<f4", 2 * PARAMS.nharms)
    err = np.max(np.abs(got - ref)) / np.max(np.abs(ref))
    assert err < 1e-6, err  # measured 1.8e-7 worst


@pytest.mark.parametrize("name", CASES)
def test_c_metric_stage(name, cdump):
    """C widths, noise scales, shapes and per-trial metrics against Julia's."""
    meta = load_case(name)
    cmeta = tomllib.loads((cdump / name / "meta.toml").read_text())
    for k in PARAMS.decimations:
        km, ck = meta[f"k{k}"], cmeta[f"k{k}"]
        for key in ("Hk", "nbins", "widths", "nvalid"):
            assert ck[key] == km[key], (k, key)
        assert ck["sigma_analytic"] == pytest.approx(km["sigma_analytic"], rel=1e-14)
        np.testing.assert_allclose(ck["shape"], km["shape"], rtol=1e-12)
        # float32 MAD against Julia's float32 one: measured 2e-7 worst
        assert ck["sigma_measured"] == pytest.approx(
            km["sigma_measured_prod"], rel=1e-6
        )
        for kind in ("analytic", "measured"):
            ref = load(name, f"metric_k{k}_{kind}.f64", "<f8")
            got = cload(cdump, name, f"metric_k{k}_{kind}.f64", "<f8")
            assert len(got) == len(ref)
            # Against the exact metrics: measured 2.9e-6 (analytic) and 3.6e-6
            # (measured) worst, the same scale as Julia's own float32 path (2.4e-6).
            np.testing.assert_allclose(
                got, ref, rtol=0, atol=2e-5, err_msg=f"k={k} {kind}"
            )


def test_c_full_search(tmp_path, search_meta):
    """The whole default search: chunk boundaries exactly, and the same raw
    candidates (threshold 6, nothing collapsed) on the same trials."""
    exe = _cs_dump()
    out = tmp_path / "search.toml"
    subprocess.run(
        [str(exe), "-search", str(HERE / PROV["input"]), str(out), "6.0"], check=True
    )
    c = tomllib.loads(out.read_text())
    assert c["total"] == search_meta["total"]
    assert c["chunk_starts"] == search_meta["chunk_starts"]
    jc = np.loadtxt(FIX / "search" / "raw_cands_thresh6.txt", ndmin=2)
    r_lo, lodr = search_meta["r_lo"], search_meta["lodr"]

    def by_trial(rows) -> dict:
        """Key each candidate by (global trial, harmonics folded).  Not by r itself:
        where the compiler fuses r_lo + t*lodr into an FMA (aarch64, Apple
        Silicon) r moves by an ulp from Julia's unfused value."""
        out = {}
        for r, _f, m, h in rows:
            k = PARAMS.nharms // int(h)
            t = int(round((r / k - r_lo) / lodr))
            assert abs(r - k * (r_lo + t * lodr)) <= 1e-12 * r
            out[(t, int(h))] = m
        return out

    jd, cd = by_trial(jc), by_trial(c["cands"])
    # A trial within float32 rounding of the threshold may land on either side.
    for key in set(jd) ^ set(cd):
        m = jd.get(key, cd.get(key))
        assert abs(m - 6.0) < 1e-4, (key, m)
    common = set(jd) & set(cd)
    assert len(common) >= len(jd) - 2, (len(common), len(jd))
    worst = max(abs(cd[k] - jd[k]) for k in common)
    assert worst < 2e-5, worst  # measured 3.8e-6 (linux-64)


def test_c_cohout(tmp_path):
    """coherentsearch itself, with every default, against the Julia CLI's output."""
    exe = _cs_dump().parent / "coherentsearch"
    assert exe.is_file(), exe
    for ext in (".fft", ".inf"):
        (tmp_path / f"synth_cs{ext}").write_bytes(
            (HERE / f"synth_cs{ext}").read_bytes()
        )
    subprocess.run([str(exe), "-noprogress", "synth_cs.fft"], cwd=tmp_path, check=True)

    def rows(path: Path) -> list[list[str]]:
        return [
            ln.split() for ln in path.read_text().splitlines() if not ln.startswith("#")
        ]

    got = rows(tmp_path / "synth_cs.cohout")
    ref = rows(FIX / "search" / "default.cohout")
    assert len(got) == len(ref) == 2
    for g, r in zip(got, ref):
        assert g[0] == r[0]  # rank
        assert float(g[1]) == pytest.approx(
            float(r[1]), abs=0.0101
        )  # S/N, printed %.2f
        assert g[2:] == r[2:]  # frequency, period, harmonics, duty cycle


def test_c_threads_deterministic(tmp_path):
    """The raw candidates (order and values) do not depend on -ncpus."""
    exe = _cs_dump()
    outs = []
    for nt in (1, 3):
        out = tmp_path / f"search{nt}.toml"
        subprocess.run(
            [str(exe), "-search", str(HERE / PROV["input"]), str(out), "6.0", str(nt)],
            check=True,
        )
        outs.append(out.read_bytes())
    assert outs[0] == outs[1]


def test_c_plots(tmp_path):
    """-plotpng writes the .ps plus (via Ghostscript) one PNG per 15 candidates, and
    -replot redraws them from the .cohout without searching again."""
    exe = _cs_dump().parent / "coherentsearch"
    for ext in (".fft", ".inf"):
        (tmp_path / f"synth_cs{ext}").write_bytes(
            (HERE / f"synth_cs{ext}").read_bytes()
        )
    args = [str(exe), "-noprogress", "-ncpus", "2", "-threshold", "6", "-noharmremove"]
    subprocess.run([*args, "-plotpng", "synth_cs.fft"], cwd=tmp_path, check=True)
    ncand = sum(
        1
        for ln in (tmp_path / "synth_cs.cohout").read_text().splitlines()
        if not ln.startswith("#")
    )
    assert (tmp_path / "synth_cs_cohout.ps").read_bytes()[:4] == b"%!PS"
    if shutil.which("gs"):
        pngs = sorted(tmp_path.glob("synth_cs_cohout_*.png"))
        assert len(pngs) == -(-ncand // 15) > 1
        assert all(p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" for p in pngs)
    (tmp_path / "synth_cs_cohout.ps").unlink()
    before = (tmp_path / "synth_cs.cohout").read_bytes()
    subprocess.run([*args, "-replot", "synth_cs.fft"], cwd=tmp_path, check=True)
    assert (tmp_path / "synth_cs.cohout").read_bytes() == before  # read, not rewritten
    ps = tmp_path / "synth_cs_cohout.ps"
    assert ps.read_bytes()[:4] == b"%!PS"


# The synthetic data relabelled as a GBT pointing at PSR B1804-08 (J1807-0847,
# P = 163.727 ms, isolated), whose fundamental lies 0.25 bins from the synthetic
# 6.1153 Hz pulsar.  The small beam keeps every other catalogue pulsar (e.g.
# J1809-0743, ~1 deg away) outside comp_psr_to_cand's 10-beam search radius.
KNOWN_INF = """\
 Data file name without suffix          =  known
 Telescope used                         =  GBT
 Instrument used                        =  None
 Object being observed                  =  B1804-08 field
 J2000 Right Ascension (hh:mm:ss.ssss)  =  18:07:38.0319
 J2000 Declination     (dd:mm:ss.ssss)  =  -08:47:43.2979
 Data observed by                       =  Unknown
 Epoch of observation (MJD)             =  60000.000000000000000
 Barycentered?           (1 yes, 0 no)  =  1
 Number of bins in the time series      =  65536
 Width of each time series bin (sec)    =  0.0005
 Any breaks in the data? (1 yes, 0 no)  =  0
 Type of observation (EM band)          =  Radio
 Beam diameter (arcsec)                 =  100
 Dispersion measure (cm-3 pc)           =  112.4
 Central freq of low channel (MHz)      =  1400
 Total bandwidth (MHz)                  =  100
 Number of channels                     =  1
 Channel bandwidth (MHz)                =  100
 Data analyzed by                       =  test
 Any additional notes:
    Synthetic data relabelled as a B1804-08 observation.
"""


def test_c_known_pulsar(tmp_path):
    """Candidates matching a catalogue pulsar near the pointing are labelled, and
    the label survives -replot (which recomputes it)."""
    exe = _cs_dump().parent / "coherentsearch"
    (tmp_path / "known.fft").write_bytes((HERE / "synth_cs.fft").read_bytes())
    (tmp_path / "known.inf").write_text(KNOWN_INF)
    subprocess.run([str(exe), "-noprogress", "known.fft"], cwd=tmp_path, check=True)
    lines = (tmp_path / "known.cohout").read_text().splitlines()
    assert lines[0].endswith("Known PSR")
    assert lines[1].split()[2] == "6.115323893229"
    assert lines[1].endswith("PSR B1804-08")
    assert len(lines[2].split()) == 6  # the 88.84 Hz pulsar matches nothing
    # PGPLOT draws text as strokes, so the name is not greppable in the .ps;
    # check that -replot found the match again instead.
    res = subprocess.run(
        [str(exe), "-replot", "known.fft"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "1 of them match known pulsars" in res.stderr
    assert (tmp_path / "known_cohout.ps").is_file()
