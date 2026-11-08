"""End-to-end tests of the PRESTO search pipeline on a synthetic pulsar.

A dispersed 23.4567 Hz pulsar (10% duty cycle, matched-filter S/N 12, DM 50) is written to a
64-channel SIGPROC filterbank file (`synthpulsar.py`), then run through prepdata/prepsubband,
realfft, rednoise, accelsearch and prepfold. Each test checks what a tool finds against the
known injected values, with thresholds set below the measured results (noted at each assert)
rather than at them, so that FFTW plan differences across platforms cannot flip a test.

A second group writes short files in every supported sample format and checks that prepdata
reproduces exactly the channel sums of the samples that were written.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

from presto import filterbank, infodata, sifting, sigproc
from synthpulsar import FakePulsar, write_filterbank

PSR = FakePulsar()
NUMOUT = PSR.nspec  # 2**17
TRUE_R = PSR.f0 * PSR.T
TRUE_P_MS = 1e3 / PSR.f0


def run(cmd: list[str], cwd: Path) -> str:
    """Run a PRESTO tool, failing the test with its output if it exits non-zero."""
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    assert res.returncode == 0, (
        f"{' '.join(cmd)} exited {res.returncode}\n--- stdout ---\n{res.stdout[-4000:]}"
        f"\n--- stderr ---\n{res.stderr[-4000:]}"
    )
    return res.stdout


def search_dat(base: str, cwd: Path) -> Path:
    """realfft + rednoise + accelsearch (-zmax 0) on base.dat; return the ACCEL_0 path."""
    run(["realfft", f"{base}.dat"], cwd)
    run(["rednoise", f"{base}.fft"], cwd)
    run(["accelsearch", "-zmax", "0", "-numharm", "16", f"{base}_red.fft"], cwd)
    return cwd / f"{base}_red_ACCEL_0"


def cand_near_truth(accelfile: Path, tol_r: float = 0.5):
    """Return the candidate at the injected fundamental (within tol_r bins), or None.

    Note that ``sifting`` replaces the ACCEL file's multi-trial sigma with a single-trial
    sigma recomputed from the optimized power, so ``.sigma`` here is higher than the sigma
    column in the ACCEL file (e.g. 8.06 against 7.19 for the default observation).
    """
    if not accelfile.exists():  # accelsearch writes no file when nothing passes its cut
        return None
    cands = sifting.candlist_from_candfile(str(accelfile)).cands
    near = [c for c in cands if abs(c.r - TRUE_R) < tol_r]
    return max(near, key=lambda c: c.sigma) if near else None


def read_bestprof(path: Path) -> dict[str, float]:
    """Parse the values used here from a prepfold .bestprof file."""
    vals: dict[str, float] = {}
    for line in path.read_text().splitlines():
        if line.startswith("# P_topo (ms)"):
            fields = line.split("=")[1].split()
            vals["P_ms"], vals["P_err_ms"] = float(fields[0]), float(fields[2])
        elif line.startswith("# Best DM"):
            vals["DM"] = float(line.split("=")[1])
        elif line.startswith("# Prob(Noise)"):
            vals["sigma"] = float(line.split("~")[1].split()[0])
    return vals


@pytest.fixture(scope="module")
def obs(tmp_path_factory) -> Path:
    """Directory holding fake.fil, its DM-50 prepdata series, and that series' search."""
    d = tmp_path_factory.mktemp("synthpsr")
    write_filterbank(PSR, str(d / "fake.fil"))
    run(["prepdata", "-nobary", "-dm", f"{PSR.dm}", "-numout", f"{NUMOUT}",
         "-o", "fake_DM50.00", "fake.fil"], d)
    search_dat("fake_DM50.00", d)
    return d


def test_prepdata_inf(obs: Path) -> None:
    """The .inf written by prepdata carries the filterbank's observing parameters."""
    inf = infodata.infodata(str(obs / "fake_DM50.00.inf"))
    assert inf.N == NUMOUT
    assert inf.dt == pytest.approx(PSR.tsamp)
    assert inf.numchan == PSR.nchan
    assert inf.DM == pytest.approx(PSR.dm)
    assert inf.BW == pytest.approx(abs(PSR.foff) * PSR.nchan)
    assert inf.lofreq == pytest.approx(PSR.fch1 + PSR.foff * (PSR.nchan - 1))
    assert inf.epoch == pytest.approx(PSR.tstart)


def test_accelsearch_finds_pulsar(obs: Path) -> None:
    """accelsearch detects the fundamental at the injected frequency."""
    cand = cand_near_truth(obs / "fake_DM50.00_red_ACCEL_0")
    assert cand is not None, "pulsar not found within 0.5 bins of the injected frequency"
    assert abs(cand.r - TRUE_R) < 0.1  # measured 0.008
    assert cand.sigma > 5.5  # measured 8.06 (single-trial; see cand_near_truth)
    assert cand.numharm >= 4  # measured 8


def test_prepsubband_dm_trials(obs: Path) -> None:
    """Across prepsubband DM trials, the pulsar is strongest at the injected DM."""
    run(["prepsubband", "-nobary", "-nsub", "32", "-lodm", "0", "-dmstep", "25",
         "-numdms", "4", "-numout", f"{NUMOUT}", "-o", "sub", "fake.fil"], obs)
    sigmas = {}
    for dm in (0, 25, 50, 75):
        cand = cand_near_truth(search_dat(f"sub_DM{dm:.2f}", obs))
        sigmas[dm] = cand.sigma if cand is not None else 0.0
    # measured (single-trial sigma): DM 0 and 25 not found, 50 -> 7.97, 75 -> 4.82
    assert sigmas[50] > 5.5
    assert all(sigmas[50] > sigmas[dm] + 2.0 for dm in (0, 25, 75)), sigmas


@pytest.mark.parametrize("raw", [False, True], ids=["dat", "filterbank"])
def test_prepfold(obs: Path, raw: bool) -> None:
    """prepfold, seeded from the accelsearch candidate, recovers the period (and DM)."""
    accel = ["-accelcand", "1", "-accelfile", "fake_DM50.00_red_ACCEL_0.cand"]
    # No P-dot was injected; searching it lets P (quoted at the start epoch) trade off
    # against P-dot, which makes the period comparison much looser.
    opts = ["-noxwin", "-nopdsearch"]
    if raw:
        cmd = ["prepfold", "-topo", *opts, *accel, "-dm", f"{PSR.dm}", "-o", "rawfold",
               "fake.fil"]
    else:
        cmd = ["prepfold", *opts, *accel, "-o", "datfold", "fake_DM50.00.dat"]
    run(cmd, obs)
    bp = read_bestprof(obs / f"{cmd[-2]}_ACCEL_Cand_1.pfd.bestprof")
    # measured |dP| = 0.2 (dat) and 0.5 (filterbank) of the quoted error
    assert abs(bp["P_ms"] - TRUE_P_MS) < 3.0 * bp["P_err_ms"]
    assert bp["sigma"] > 5.5  # measured 7.6 and 7.7
    if raw:
        assert abs(bp["DM"] - PSR.dm) < 2.0  # measured 50.58; DM resolution here is ~1


# Short files in every format prepdata reads, checked sample by sample. nspec is a multiple
# of the SIGPROC reader's 2400-spectrum block, so prepdata pads nothing at the end.
FORMATS = [(8, False, np.uint8), (8, True, np.int8), (16, False, np.uint16),
           (32, False, np.float32)]


@pytest.mark.parametrize("nbits,signed,dtype", FORMATS,
                         ids=["8bit", "8bit-signed", "16bit", "32bit-float"])
def test_filterbank_formats_roundtrip(tmp_path: Path, nbits: int, signed: bool, dtype) -> None:
    """prepdata -dm 0 output equals the channel sum of exactly the samples written.

    The signed case is the aarch64 regression test: plain char is unsigned there, which
    once made signed 8-bit data read back wrapped to 0..255 (fixed in v6.1.1).
    """
    psr = FakePulsar(nspec=4800, snr=50.0)
    fn = tmp_path / "fmt.fil"
    write_filterbank(psr, str(fn), nbits=nbits, signed=signed)
    run(["prepdata", "-nobary", "-noclip", "-dm", "0", "-o", "fmt", "fmt.fil"], tmp_path)
    _, hdrlen = sigproc.read_header(str(fn))
    raw = np.fromfile(fn, dtype=dtype, offset=hdrlen).reshape(-1, psr.nchan)
    expected = raw.astype(np.float64).sum(axis=1)
    dat = np.fromfile(tmp_path / "fmt.dat", dtype=np.float32)
    assert len(dat) == psr.nspec
    if signed:
        assert raw.min() < 0  # the test only means something if negatives are present
    np.testing.assert_allclose(dat, expected, rtol=0, atol=1e-3 * np.abs(expected).max())



@pytest.mark.parametrize("nbits,signed,dtype", FORMATS,
                         ids=["8bit", "8bit-signed", "16bit", "32bit-float"])
def test_python_filterbank_reader(tmp_path: Path, nbits: int, signed: bool, dtype) -> None:
    """presto.filterbank reads back exactly the samples written, signed 8-bit included."""
    psr = FakePulsar(nspec=2400, snr=50.0)
    fn = str(tmp_path / "py.fil")
    write_filterbank(psr, fn, nbits=nbits, signed=signed)
    _, hdrlen = sigproc.read_header(fn)
    raw = np.fromfile(fn, dtype=dtype, offset=hdrlen).reshape(-1, psr.nchan)
    fb = filterbank.FilterbankFile(fn)
    assert fb.dtype == np.dtype(dtype).name
    assert fb.signed == signed
    np.testing.assert_array_equal(fb.get_spectra(0, psr.nspec).data, raw.T)

def test_filterbank_long_header_strings(tmp_path: Path) -> None:
    """Header strings longer than the reader's buffers are truncated, not fatal.

    Before the fix, any string over 80 characters (e.g. a long rawdatafile path, as
    macOS temporary directories produce) left the header reader out of sync, and it
    died with "unknown parameter: ERROR".
    """
    psr = FakePulsar(nspec=2400)
    name = "J" + "x" * 99  # 100 characters
    write_filterbank(psr, str(tmp_path / "long.fil"), rawdatafile="/" + "d" * 299,
                     source_name=name)
    run(["prepdata", "-nobary", "-dm", "0", "-o", "long", "long.fil"], tmp_path)
    inf = infodata.infodata(str(tmp_path / "long.inf"))
    assert inf.N == psr.nspec
    assert inf.dt == pytest.approx(psr.tsamp)
    assert name.startswith(inf.object) and len(inf.object) >= 70



@pytest.mark.parametrize("nbits", [1, 2, 4])
def test_lowbit_filterbank_channels(tmp_path: Path, nbits: int) -> None:
    """1/2/4-bit files written by presto.filterbank read back channel by channel.

    Independent random values in every channel are checked against both the Python
    reader and PRESTO's C reader (via prepsubband -sub, one subband per channel at DM 0),
    so any error in the in-byte channel order shows up. The packing (first channel in the
    lowest bits) was also checked against the real 4-bit GBT_Lband_PSR.fil: the Python and
    C readers agree on every sample, and rewriting it reproduces its data bytes exactly.
    """
    nchan, nspec = 64, 24000
    psr = FakePulsar(nchan=nchan, nspec=nspec)
    rng = np.random.default_rng(nbits)
    vals = rng.integers(0, 2**nbits, size=(nspec, nchan))
    header = {"telescope_id": 0, "machine_id": 0, "data_type": 1, "source_name": "LOWBIT",
              "src_raj": 0.0, "src_dej": 0.0, "tstart": psr.tstart, "tsamp": psr.tsamp,
              "fch1": psr.fch1, "foff": psr.foff, "nchans": nchan, "nifs": 1}
    fn = str(tmp_path / "low.fil")
    filterbank.create_filterbank_file(fn, header, spectra=vals, nbits=nbits)
    assert (tmp_path / "low.fil").stat().st_size > nspec * nchan * nbits // 8
    fb = filterbank.FilterbankFile(fn)
    np.testing.assert_array_equal(fb.get_spectra(0, nspec).data, vals.T)
    run(["prepsubband", "-nobary", "-noclip", "-sub", "-subdm", "0", "-nsub", f"{nchan}",
         "-o", "low", "low.fil"], tmp_path)
    subs = sorted(tmp_path.glob("low_DM0.00.sub[0-9]*"))
    assert len(subs) == nchan
    S = np.array([np.fromfile(f, dtype=np.int16) for f in subs])
    n = S.shape[1]  # prepsubband drops the final partial block
    assert n > nspec // 2
    # sub00 is the lowest frequency, i.e. the last file channel (foff < 0)
    np.testing.assert_array_equal(S, vals.T[::-1, :n])


@pytest.mark.parametrize("nbits", [4, 2, 1])
def test_lowbit_pipeline(tmp_path: Path, nbits: int) -> None:
    """The pulsar is still found in 4-, 2- and 1-bit versions of the observation."""
    write_filterbank(PSR, str(tmp_path / "low.fil"), nbits=nbits)
    run(["prepdata", "-nobary", "-dm", f"{PSR.dm}", "-numout", f"{NUMOUT}",
         "-o", "low_DM50.00", "low.fil"], tmp_path)
    cand = cand_near_truth(search_dat("low_DM50.00", tmp_path))
    assert cand is not None
    # measured: |dr| 0.13 / 0.13 / 0.008 and single-trial sigma 7.23 / 6.68 / 6.54 at
    # 4 / 2 / 1 bits (8-bit: 8.06), the expected quantization losses
    assert abs(cand.r - TRUE_R) < 0.25
    assert cand.sigma > 5.5
