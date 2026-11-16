"""Check accelsearch's harmonic-sum gather kernels against the scalar loop.

``accel_gather_test`` (built from ``tests/accelsearch/gather_test.c``, not
installed) runs every version of ``add_gathered_row()`` the CPU supports (scalar,
AVX2, AVX-512, NEON and the one the runtime dispatch picks) on the r index
tables accelsearch builds for every harmonic up to ``-numharm 32``, and checks
that each gives the scalar sums bit for bit.
"""

import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


def _gather_test() -> Path:
    """The ``accel_gather_test`` helper: ``$ACCEL_GATHER_TEST``, else the pixi build's."""
    if "ACCEL_GATHER_TEST" in os.environ:
        return Path(os.environ["ACCEL_GATHER_TEST"])
    build = REPO / "build-pixi"
    exe = build / "src" / "accel_gather_test"
    if not build.is_dir():
        pytest.skip(
            "no build-pixi/ (run `pixi run build`) and $ACCEL_GATHER_TEST is not set"
        )
    # A pixi build that lacks it is a build problem, not a reason to skip.
    assert exe.is_file(), f"{exe} missing: rebuild with `pixi run build`"
    return exe


def test_gather_kernels_match_scalar() -> None:
    """Every gather kernel gives the scalar sums exactly."""
    res = subprocess.run(
        [str(_gather_test())], capture_output=True, text=True, check=False
    )
    print(res.stdout)
    assert res.returncode == 0, res.stdout + res.stderr
    # The scalar reference and the dispatched kernel always run.
    assert "scalar " in res.stdout
    assert "dispatched " in res.stdout
