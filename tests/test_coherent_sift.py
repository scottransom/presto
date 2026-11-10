"""Tests for bin/coherent_sift.py on synthetic .cohout files with known answers.

Two observations, each searched at DMs 10-16, hold:
- an unknown pulsar at 12.3456789 Hz, peaking at DM 13 (seen at DMs 12-14);
- its second harmonic, weaker, at the same DMs;
- a known pulsar at 33.3 Hz, labelled by coherentsearch, in obs A only;
- broadband RFI at 60 Hz, at every DM, strongest at the lowest;
- a 0.5 Hz signal found only by the deep pass (``*_deep.cohout``), in obs A.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location(
    "coherent_sift", HERE.parent / "bin" / "coherent_sift.py"
)
cs = importlib.util.module_from_spec(SPEC)
sys.modules["coherent_sift"] = cs  # dataclasses look the module up by name
SPEC.loader.exec_module(cs)

DMS = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
HEADER = "#Num    'S/N'      Frequency (Hz)        Period (ms)    #Harm  Ducy(%)"
PSR_F, KNOWN_F, RFI_F, DEEP_F = 12.3456789, 33.3, 60.0, 0.5


def row(snr: float, f: float, nharm: int = 30, known: str = "") -> tuple:
    """One candidate: (S/N, frequency, harmonics, label)."""
    return (snr, f, nharm, known)


def write_cohout(path: Path, rows: list[tuple]) -> None:
    """Write rows (strongest first) in coherentsearch's format."""
    lines = [HEADER + ("   Known PSR" if any(r[3] for r in rows) else "")]
    for i, (snr, f, nharm, known) in enumerate(sorted(rows, reverse=True), 1):
        line = "%-4d  %8.2f  %18.12f  %18.12f   %3d   %6s" % (
            i,
            snr,
            f,
            1000.0 / f,
            nharm,
            "5.00",
        )
        lines.append(line + ("   " + known if known else ""))
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def cands(tmp_path: Path) -> Path:
    """The synthetic candidate files."""
    psr = {12.0: 8.0, 13.0: 12.0, 14.0: 9.0}
    for obs in ("A", "B"):
        for dm in DMS:
            rows = [row(20.0 - dm, RFI_F * (1 + 1e-7 * dm), 10)]
            if dm in psr:
                rows.append(row(psr[dm], PSR_F))
                rows.append(row(psr[dm] - 2.0, 2 * PSR_F, 20))
            if obs == "A" and dm in (11.0, 12.0, 13.0):
                rows.append(row(10.0 + dm - 11.0, KNOWN_F, 20, "PSR B1234+56"))
            write_cohout(tmp_path / f"obs{obs}_DM{dm:.2f}_red.cohout", rows)
        if obs == "A":
            for dm, snr in ((13.0, 9.0), (14.0, 8.0)):
                write_cohout(
                    tmp_path / f"obs{obs}_DM{dm:.2f}_red_deep.cohout",
                    [row(snr, DEEP_F, 120)],
                )
    return tmp_path


def sift(*args: str) -> dict[float, list[str]]:
    """Run the sift; return {rounded freq: [class, rank, flags]}."""
    out = {}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        cs.main(list(args))
    for line in buf.getvalue().splitlines():
        if line.startswith("#"):
            continue
        p = line.split(None, 12)
        out[round(float(p[2]), 3)] = [p[1], int(p[0]), p[12] if len(p) > 12 else ""]
    return out


def test_classes(cands: Path) -> None:
    """Each kind of signal gets the right class, and the pulsar ranks first."""
    res = sift(str(cands))
    assert res[round(PSR_F, 3)][:2] == ["candidate", 1]
    assert res[round(2 * PSR_F, 3)][0] == "harmonic"
    assert "harm(x2 of 12.345679Hz)" in res[round(2 * PSR_F, 3)][2]
    assert res[round(KNOWN_F, 3)][0] == "known"
    assert "PSR B1234+56" in res[round(KNOWN_F, 3)][2]
    assert res[round(RFI_F, 3)][0] == "RFI-like"
    assert "full-DM-span" in res[round(RFI_F, 3)][2]
    assert "deep" in res[DEEP_F][2]


def test_options(cands: Path) -> None:
    """--no-deep drops the deep-pass files, and --no-known the labelled pulsar."""
    assert DEEP_F not in sift("--no-deep", str(cands))
    assert round(KNOWN_F, 3) not in sift("--no-known", str(cands))


def test_single_observation(cands: Path) -> None:
    """With one observation, a DM-extended interior peak is still a candidate,
    and nothing is flagged single-obs."""
    res = sift(*map(str, sorted(cands.glob("obsA_*.cohout"))))
    assert res[round(PSR_F, 3)][:2] == ["candidate", 1]
    assert not any("single-obs" in v[2] for v in res.values())


def test_html(cands: Path, tmp_path: Path) -> None:
    """The HTML report is written and names the signals."""
    page = tmp_path / "sift.html"
    sift("--html", str(page), str(cands))
    text = page.read_text()
    assert "12.345679" in text and "PSR B1234+56" in text
