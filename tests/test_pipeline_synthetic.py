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

from presto import filterbank, infodata, prepfold, psr_utils, sifting, sigproc
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


def test_accelsearch_nonorm(obs: Path) -> None:
    """accelsearch -nonorm finds the pulsar in the rednoised spectrum, and refuses a raw one.

    -nonorm skips the per-block median normalization (and allows shorter FFTs), so it
    relies on the input already being normalized: on the raw realfft output, whose
    powers are far from 1, it must exit with an error rather than pass every trial.
    """
    d = obs / "nonorm"  # its own directory: the ACCEL files would overwrite obs's
    d.mkdir(exist_ok=True)
    for ext in ("_red.fft", "_red.inf", ".fft", ".inf"):
        link = d / f"fake_DM50.00{ext}"
        if not link.exists():
            link.symlink_to(obs / f"fake_DM50.00{ext}")
    out = run(["accelsearch", "-zmax", "0", "-numharm", "16", "-nonorm",
               "fake_DM50.00_red.fft"], d)
    assert "no block normalization" in out
    cand = cand_near_truth(d / "fake_DM50.00_red_ACCEL_0")
    assert cand is not None, "pulsar not found with -nonorm"
    assert abs(cand.r - TRUE_R) < 0.1  # measured 0.008
    assert cand.sigma > 5.5  # measured 8.06, as without -nonorm
    assert cand.numharm >= 4  # measured 8

    res = subprocess.run(["accelsearch", "-zmax", "0", "-numharm", "16", "-nonorm",
                          "fake_DM50.00.fft"], cwd=d, capture_output=True, text=True)
    assert res.returncode != 0
    assert "more than 2x off" in res.stdout + res.stderr


def test_accelsearch_auto_nonorm(obs: Path) -> None:
    """A *_red.fft that passes the normalization check is searched as with -nonorm.

    -median forces the block median back, and a *_red.fft whose powers are not
    normalized falls back to the block median (with a note) instead of exiting.
    """
    d = obs / "autonorm"
    d.mkdir(exist_ok=True)
    links = {"fake_DM50.00_red": "fake_DM50.00_red", "raw_red": "fake_DM50.00"}
    for name, target in links.items():
        for ext in (".fft", ".inf"):
            link = d / f"{name}{ext}"
            if not link.exists():
                link.symlink_to(obs / f"{target}{ext}")
    search = ["accelsearch", "-zmax", "0", "-numharm", "16"]

    out = run([*search, "fake_DM50.00_red.fft"], d)
    assert "as with -nonorm" in out
    cand = cand_near_truth(d / "fake_DM50.00_red_ACCEL_0")
    assert cand is not None, "pulsar not found with automatic -nonorm"
    assert abs(cand.r - TRUE_R) < 0.1  # measured 0.008
    assert cand.sigma > 5.5  # measured 8.06

    out = run([*search, "-median", "fake_DM50.00_red.fft"], d)
    assert "using median-blocks" in out and "-nonorm" not in out

    out = run([*search, "raw_red.fft"], d)  # raw realfft powers under a _red name
    assert "not treated as" in out and "as with -nonorm" not in out


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


@pytest.mark.parametrize("nspec", [1000, 2400])
def test_prepfold_short_dat(tmp_path: Path, nspec: int) -> None:
    """prepfold refuses a .dat shorter than one 1024-point block, and cuts -npart to fit.

    It once folded zero points in both cases, and then printed NaNs and overflowed a
    buffer formatting them.
    """
    write_filterbank(FakePulsar(nspec=nspec), str(tmp_path / "short.fil"))
    run(["prepdata", "-nobary", "-dm", "0", "-o", "short", "short.fil"], tmp_path)
    res = subprocess.run(["prepfold", "-nosearch", "-noxwin", "-p", "0.05", "-n", "16",
                          "-npart", "4", "-o", "short", "short.dat"],
                         cwd=tmp_path, capture_output=True, text=True)
    if nspec < 1024:
        assert res.returncode == 1
        assert "prepfold needs at least 1024" in res.stdout
    else:
        assert res.returncode == 0, res.stdout[-4000:] + res.stderr[-4000:]
        assert "Overriding -npart to be 2" in res.stdout
        assert "Total number of data points  =  2048" in res.stdout
        assert "Amount Complete = 100%" in res.stdout
        assert "Amount Complete = 101%" not in res.stdout
        bestprof = read_bestprof(next(tmp_path.glob("short*.pfd.bestprof")))
        assert bestprof["P_ms"] == pytest.approx(50.0)


# SIGPROC tstart and the DATE-OBS-style string readfile must print for it.  The old
# formatter's "%.6g" seconds gave 60 seconds, "010" for 9.9999999, and exponents.
DATES = {
    "midnight": (60000.0, "2023-02-25T00:00:00.0000"),
    "normal": (60000.0 + (12 * 3600 + 34 * 60 + 56.789) / 86400, "2023-02-25T12:34:56.7890"),
    "tiny-sec": (60000.0 + 2e-5 / 86400, "2023-02-25T00:00:00.0000"),
    "carry-to-10s": (60000.0 + 9.9999999 / 86400, "2023-02-25T00:00:10.0000"),
    "carry-to-min": (60000.0 + (11 * 60 + 59.99996) / 86400, "2023-02-25T00:12:00.0000"),
    "carry-to-year": (60310.0 - 4e-5 / 86400, "2024-01-01T00:00:00.0000"),
}


@pytest.mark.parametrize("tstart,datestr", DATES.values(), ids=DATES.keys())
def test_filterbank_date_string(tmp_path: Path, tstart: float, datestr: str) -> None:
    """readfile prints a SIGPROC tstart as a properly rounded date and time."""
    write_filterbank(FakePulsar(nspec=2400, tstart=tstart), str(tmp_path / "date.fil"))
    out = run(["readfile", "date.fil"], tmp_path)
    assert f"Obs Date String = {datestr}\n" in out


# SIGPROC (src_raj, src_dej) and the RA/Dec strings they must give.  The cases between 0
# and -1 deg are the classic lost-sign trap: the degrees field is zero, so the sign has to
# live on the minutes or seconds (fixed for SIGPROC in v6.3.1).  The last two check that
# rounding to 4 decimals carries into the minutes, rather than printing 60.0000 seconds.
POSITIONS = {
    "missing": (None, "00:00:00.0000", "00:00:00.0000", False),
    "normal": ((123456.789, -54321.5), "12:34:56.7890", "-05:43:21.5000", False),
    "north": ((10203.04, 101214.966), "01:02:03.0400", "10:12:14.9660", False),
    "dec-0to-1-arcmin": ((10203.04, -1214.966), "01:02:03.0400", "-00:12:14.9660", False),
    "dec-0to-1-arcsec": ((10203.04, -30.5), "01:02:03.0400", "-00:00:30.5000", False),
    "dec-0to-1-sub-arcsec": ((10203.04, -0.25), "01:02:03.0400", "-00:00:00.2500", False),
    "dec+0to1-arcsec": ((10203.04, 30.5), "01:02:03.0400", "00:00:30.5000", False),
    "carry-to-deg": ((10203.04, -5959.99999), "01:02:03.0400", "-01:00:00.0000", False),
    "carry-to-min": ((105959.99996, -1259.99996), "11:00:00.0000", "-00:13:00.0000", False),
    "bad-ra": ((-4321.0, -4321.0), "00:00:00.0000", "-00:43:21.0000", True),
    "huge": ((1e300, 1e300), "00:00:00.0000", "00:00:00.0000", True),
}


@pytest.mark.parametrize("position,ra,dec,warns", POSITIONS.values(), ids=POSITIONS.keys())
def test_filterbank_position(tmp_path: Path, position, ra: str, dec: str,
                             warns: bool) -> None:
    """The SIGPROC position reaches the .inf and the .pfd, with the right signs.

    A missing or bad position gives 00:00:00.  sigproc's ``fake`` writes no
    src_raj/src_dej. The reader once left them as stack garbage, and a large value
    overflowed the 40-byte RA/Dec strings (PR #234).  prepfold on the .dat then reads
    the .inf back with the C parser and writes the position into the .pfd again.
    """
    psr = FakePulsar(nspec=2400)
    write_filterbank(psr, str(tmp_path / "pos.fil"), position=position)
    res = subprocess.run(["prepdata", "-nobary", "-dm", "0", "-o", "pos", "pos.fil"],
                         cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout[-4000:] + res.stderr[-4000:]
    assert ("Warning:  bad src_" in res.stderr) == warns
    inf = infodata.infodata(str(tmp_path / "pos.inf"))
    assert inf.RA == ra
    assert inf.DEC == dec
    if position is not None and not warns:
        # The strings agree with the python SIGPROC decoder (and so carry the sign), to
        # their printed precision of 1e-4 s (RA, 7.3e-9 rad) or 1e-4" (Dec, 4.8e-10 rad)
        assert psr_utils.ra_to_rad(ra) == pytest.approx(sigproc.ra2radians(position[0]),
                                                        abs=4e-9)
        assert psr_utils.dec_to_rad(dec) == pytest.approx(sigproc.dec2radians(position[1]),
                                                          abs=3e-10)
    run(["prepfold", "-nosearch", "-noxwin", "-p", "0.05", "-n", "16", "-npart", "4",
         "-o", "pos", "pos.dat"], tmp_path)
    pfd = prepfold.pfd(str(next(tmp_path.glob("pos*.pfd"))))
    assert (pfd.rastr.decode(), pfd.decstr.decode()) == (ra, dec)



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


# -flatten: a band with a 10x gain slope, four dead channels and three noisy ones.
FLAT_DEAD = [0, 1, 2, 3]  # file channels (the highest frequencies, as foff < 0)
FLAT_NOISY = {10: 10.0, 11: 10.0, 40: 4.0}  # file channel: extra noise (no pulsar in it)


@pytest.fixture(scope="module")
def flatobs(tmp_path_factory) -> Path:
    """Directory holding bandpass.fil (see FLAT_DEAD/FLAT_NOISY) and its rfifind mask.

    The channel means are set so that sigma/mean equals the radiometer value
    1/sqrt(|foff| * tsamp), as for real total-power data, so -flatten uses its
    sigma/mean test as well.
    """
    d = tmp_path_factory.mktemp("flatten")
    gain = np.geomspace(3.0, 0.3, PSR.nchan)
    gain[FLAT_DEAD] = 0.0
    extra = np.zeros(PSR.nchan)
    for ch, val in FLAT_NOISY.items():
        extra[ch] = val
    offset = gain * np.sqrt(abs(PSR.foff) * 1e6 * PSR.tsamp)
    write_filterbank(PSR, str(d / "bandpass.fil"), nbits=32, scale=1.0, offset=offset,
                     chan_gain=gain, chan_extra_noise=extra)
    run(["rfifind", "-time", "2.0", "-o", "bandpass", "bandpass.fil"], d)
    return d


def read_flatten(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the (freq, weight, flag) columns of a .flatten file."""
    rows = [ln.split() for ln in path.read_text().splitlines() if not ln.startswith("#")]
    freq = np.array([float(r[1]) for r in rows])
    wt = np.array([float(r[4]) for r in rows])
    flag = np.array([r[5] for r in rows])
    return freq, wt, flag


@pytest.mark.parametrize("use_mask", [False, True], ids=["prepass", "rfifind-stats"])
def test_flatten(flatobs: Path, use_mask: bool) -> None:
    """-flatten zero-weights the dead and noisy channels, keeps the rest, and helps S/N."""
    base = f"flat{'mask' if use_mask else 'nomask'}_DM50.00"
    maskopt = ["-mask", "bandpass_rfifind.mask"] if use_mask else []
    out = run(["prepdata", "-nobary", "-dm", f"{PSR.dm}", "-numout", f"{NUMOUT}",
               "-flatten", *maskopt, "-o", base, "bandpass.fil"], flatobs)
    assert ("rfifind .stats" in (flatobs / f"{base}.flatten").read_text()) == use_mask
    assert "using the sigma/mean test" in out
    freq, wt, flag = read_flatten(flatobs / f"{base}.flatten")
    filefreqs = PSR.freqs
    zeroed = {int(np.argmin(np.abs(filefreqs - f))) for f in freq[wt == 0.0]}
    assert set(FLAT_DEAD) | set(FLAT_NOISY) == zeroed, flag
    for ch in FLAT_DEAD:
        assert flag[np.argmin(np.abs(freq - filefreqs[ch]))] == "dead"
    # Each good channel has unit variance after flattening, so the dedispersed series
    # has a mean near 0 and a standard deviation near sqrt(number of good channels).
    dat = np.fromfile(flatobs / f"{base}.dat", dtype=np.float32)
    ngood = int(np.sum(wt > 0.0))
    assert abs(dat.mean()) < 0.1 * np.sqrt(ngood)
    assert dat.std() == pytest.approx(np.sqrt(ngood), rel=0.05)
    cand = cand_near_truth(search_dat(base, flatobs))
    assert cand is not None
    # measured: |dr| 0.13 / 0.13, single-trial sigma 6.73 / 6.72 (prepass / stats)
    assert abs(cand.r - TRUE_R) < 0.25
    assert cand.sigma > 5.5


def test_flatten_beats_plain(flatobs: Path) -> None:
    """Without -flatten the noisy channels swamp the pulsar.

    With equal S/N in every channel, 1/sigma weights should also beat equal weights by
    ~17% on the clean channels alone, but one noise realization can't show that reliably.
    """
    run(["prepdata", "-nobary", "-dm", f"{PSR.dm}", "-numout", f"{NUMOUT}",
         "-o", "plain_DM50.00", "bandpass.fil"], flatobs)
    plain = cand_near_truth(search_dat("plain_DM50.00", flatobs))
    run(["prepdata", "-nobary", "-dm", f"{PSR.dm}", "-numout", f"{NUMOUT}",
         "-flatten", "-o", "flat_DM50.00", "bandpass.fil"], flatobs)
    flat = cand_near_truth(search_dat("flat_DM50.00", flatobs))
    assert flat is not None
    # measured: flat sigma 6.73; plain: not detected
    assert plain is None or flat.sigma > plain.sigma + 1.0


def test_flatten_prepsubband(flatobs: Path) -> None:
    """prepsubband -flatten finds the pulsar at the right DM trial."""
    run(["prepsubband", "-nobary", "-flatten", "-lodm", f"{PSR.dm - 2}", "-dmstep", "1",
         "-numdms", "5", "-nsub", "16", "-numout", f"{NUMOUT}", "-o", "flatsub",
         "bandpass.fil"], flatobs)
    assert (flatobs / "flatsub.flatten").exists()
    cand = cand_near_truth(search_dat(f"flatsub_DM{PSR.dm:.2f}", flatobs))
    assert cand is not None
    # measured: |dr| 0.008, single-trial sigma 7.53
    assert cand.sigma > 5.5


def test_flatten_prepfold(flatobs: Path) -> None:
    """prepfold -flatten on the raw file beats a plain fold, which the noisy channels swamp."""
    sig = {}
    for name, opt in (("pfplain", []), ("pfflat", ["-flatten"])):
        run(["prepfold", "-topo", "-noxwin", "-nosearch", "-p", f"{TRUE_P_MS / 1e3:.12f}",
             "-dm", f"{PSR.dm}", "-nsub", "16", *opt, "-o", name, "bandpass.fil"], flatobs)
        bp = read_bestprof(next(flatobs.glob(f"{name}_*.pfd.bestprof")))
        sig[name] = bp["sigma"]
    assert (flatobs / "pfflat.flatten").exists()
    # measured: sigma 7.1 with -flatten, 0.8 without
    assert sig["pfflat"] > 5.5
    assert sig["pfflat"] > sig["pfplain"] + 1.0
