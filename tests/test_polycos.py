"""Tests for presto.polycos, the reader for TEMPO1-format polycos."""

import shutil
from pathlib import Path

import pytest

from presto import polycos

PARFILE = Path(__file__).parent / "prepfold" / "M28A.par"

# Two entries of "tempo2 -tempo1 -polyco" output for M28A at MeerKAT, as
# tempo2 writes them when given the site name "meerkat".  It prints the name
# with "%5s", so the 7-character name runs straight into F0.
FUSED = """\
1824-2452A 27-Jun-23  110000.00   60122.45833333330           119.879594 -0.009 -7.499
 178848880090.968488  327.405606051598meerkat   60   12  1284.000                
  4.05578682047085231e-09  1.19079031868709623e-02  4.05828381848839662e-05
 -6.89427088469414018e-09 -8.15716733169101330e-11 -7.85223115090086463e-14
 -2.01471707475872991e-15  1.99288439019170116e-17  1.39093388208051126e-18
  6.64885694481268117e-20 -2.83802734541787669e-22 -3.86935606154776053e-23
1824-2452A 27-Jun-23  120000.00   60122.49999999990           119.879594 -0.011 -7.448
 178850058752.006399  327.405606051598meerkat   60   12  1284.000                
 -6.43730616801905115e-08  1.66321767484815045e-02  3.75727119576921067e-05
 -2.63643173165391361e-08 -7.77119901490081160e-11  9.55638284115987202e-14
 -5.30047148030313050e-16  3.14739299200799202e-17  7.37938941454610344e-20
 -1.18246027690716400e-19  1.60747216064121929e-22  5.46191884551204881e-23
"""


def test_site_name_fused_to_f0(tmp_path) -> None:
    """A long site name stuck to F0 parses the same as a short, separate one."""
    fused = tmp_path / "fused.dat"
    fused.write_text(FUSED)
    spaced = tmp_path / "spaced.dat"
    spaced.write_text(FUSED.replace("598meerkat", "598    m"))
    pf = polycos.polycos("1824-2452A", str(fused))
    ps = polycos.polycos("1824-2452A", str(spaced))
    assert len(pf.polycos) == len(ps.polycos) == 2
    for a, b in zip(pf.polycos, ps.polycos):
        assert a.obs == "meerkat" and b.obs == "m"
        assert a.F0 == b.F0 == 327.405606051598
        assert (a.dataspan, a.numcoeff, a.obsfreq) == (60, 12, 1284.0)
        assert (b.dataspan, b.numcoeff, b.obsfreq) == (60, 12, 1284.0)
        assert list(a.coeffs) == list(b.coeffs)
    assert pf.get_phs_and_freq(60122, 0.47) == ps.get_phs_and_freq(60122, 0.47)


def test_plain_f0_not_split(tmp_path) -> None:
    """An F0 followed by whitespace is left alone (no digits taken as the site)."""
    fn = tmp_path / "plain.dat"
    fn.write_text(FUSED.replace("598meerkat", "598  gbt"))
    p = polycos.polycos("1824-2452A", str(fn)).polycos[0]
    assert (p.F0, p.obs, p.dataspan) == (327.405606051598, "gbt", 60)


def test_tables_consistent() -> None:
    """Every TEMPO id leads to a tempo2 site name and a track length."""
    for code, name in polycos.id_to_telescope.items():
        assert name in polycos.telescope_to_tempo2_name, code
        assert name in polycos.telescope_to_maxha, code
    t2 = polycos.telescope_to_tempo2_name
    for name, code in polycos.telescope_to_id.items():  # "a" has two names
        assert t2[polycos.id_to_telescope[code]] == t2[name]
    # A space only separates the site from F0 if the name is <= 4 characters
    assert polycos.telescope_to_tempo2_name["MeerKAT"] == "m"


NEEDS_TEMPO2 = pytest.mark.skipif(shutil.which("tempo2") is None, reason="needs tempo2")


def check_create_polycos(code: str, tmp_path, monkeypatch) -> None:
    """create_polycos() runs tempo2 for TEMPO site `code`; its polycos read back."""
    monkeypatch.chdir(tmp_path)
    p = polycos.create_polycos(str(PARFILE), code, 1284.0, 60123, 60124)
    t2name = polycos.telescope_to_tempo2_name[polycos.id_to_telescope[code]]
    assert len(p.polycos) > 0
    assert all(pc.obs == t2name for pc in p.polycos)
    assert all((pc.dataspan, pc.numcoeff) == (60, 12) for pc in p.polycos)
    # The topocentric spin frequency is within the Earth's Doppler shift of F0.
    # (Sites with short tracks only get polycos while M28A is up, so use the
    # first one's epoch.)
    f = p.get_freq(int(p.polycos[0].TMIDi), p.polycos[0].TMIDf)
    assert abs(f / 327.4056060515984 - 1.0) < 1.1e-4


@NEEDS_TEMPO2
def test_create_polycos_meerkat(tmp_path, monkeypatch) -> None:
    """MeerKAT polycos from tempo2 read back correctly.

    tempo2 takes ~2 s to start, so every other site is checked in the weekly
    real-data run instead (tests/test_polycos_sites.py).
    """
    check_create_polycos("m", tmp_path, monkeypatch)
