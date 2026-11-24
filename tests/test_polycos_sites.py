"""Polycos from tempo2 for every site PRESTO knows.

tempo2 takes ~2 s to start, so this runs ~30 of them in about a minute: too slow
for every push, so it is part of the weekly real-data run (`pixi run
test-polycos-sites`).  tests/test_polycos.py checks MeerKAT on every push.
"""

import subprocess

import pytest
from test_polycos import NEEDS_TEMPO2, PARFILE, check_create_polycos

from presto import polycos

pytestmark = NEEDS_TEMPO2


@pytest.mark.parametrize("code", sorted(polycos.id_to_telescope))
def test_create_polycos(code, tmp_path, monkeypatch) -> None:
    """create_polycos() works for every TEMPO site id."""
    check_create_polycos(code, tmp_path, monkeypatch)


# Sites with a tempo2 name but no TEMPO id, so create_polycos() can't reach them
NO_ID = sorted(
    set(polycos.telescope_to_tempo2_name.values())
    - {polycos.telescope_to_tempo2_name[n] for n in polycos.id_to_telescope.values()}
)


@pytest.mark.parametrize("t2name", NO_ID)
def test_tempo2_name(t2name, tmp_path) -> None:
    """tempo2 knows the site, and its polycos read back."""
    cmd = ["tempo2", "-tempo1", "-f", str(PARFILE), "-polyco"]
    cmd.append(f"60123 60124 60 12 12 {t2name} 1284.0")
    subprocess.run(cmd, cwd=tmp_path, capture_output=True, check=True)
    p = polycos.polycos("1824-2452A", str(tmp_path / "polyco_new.dat"))
    assert len(p.polycos) > 0
    assert all(pc.obs == t2name for pc in p.polycos)
