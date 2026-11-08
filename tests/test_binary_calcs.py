"""Checks of PRESTO's binary-orbit routines (keplers_eqn, dorbint, E_to_v, E_to_phib).

These began as plotting demos reproducing Figure 1 of Hulse & Taylor (1975) and Figure 1 of
Champion et al. (2008); they now assert what those plots showed.
"""

from pathlib import Path

import numpy as np

from presto import binary_psr, presto

HERE = Path(__file__).resolve().parent
N = 1000  # points per calculation, spanning two orbits
MA = np.arange(float(N)) * 2.0 / N  # orbital phase


def _orbit_E(orb) -> np.ndarray:
    """Return the eccentric anomaly over two orbits, starting at orb.t = 0."""
    orb.t = 0.0
    Eo = presto.keplers_eqn(orb.t, orb.p, orb.e, 1e-15)
    return presto.dorbint(Eo, N, 2.0 * orb.p / N, orb)


def test_velocity_curve_matches_closed_form() -> None:
    """B1913+16 radial velocity spans K(-1 + e cos w) .. K(1 + e cos w)."""
    psr = presto.psrepoch("B1913+16", 42320.0)
    # OMDOT is not handled correctly in the database's binary entry, so set w by hand
    psr.orb.w = 179.0
    vs = _orbit_E(psr.orb)
    presto.E_to_v(vs, psr.orb)
    o = psr.orb
    K = 2.0 * np.pi * o.x * 299792.458 / (o.p * np.sqrt(1.0 - o.e**2))  # km/s
    ecw = o.e * np.cos(np.radians(o.w))
    # 1e-3: N samples only approximate the sharp periastron extremum
    np.testing.assert_allclose(vs.min(), K * (-1.0 + ecw), rtol=1e-3)
    np.testing.assert_allclose(vs.max(), K * (1.0 + ecw), rtol=1e-3)


def test_C_delays_match_binary_psr() -> None:
    """E_to_phib (C) agrees with binary_psr.position (Python) for J1903+0327."""
    parfile = str(HERE / "1903+0327.par")
    bpsr = binary_psr.binary_psr(parfile)
    xs, _ = bpsr.position(bpsr.T0 + MA * bpsr.par.PB)
    psr = presto.psrepoch(parfile, bpsr.T0)
    delays = _orbit_E(psr.orb)
    presto.E_to_phib(delays, psr.orb)
    # Measured agreement is 2.5e-4 s out of ~128 s
    assert np.abs(delays + xs).max() < 1e-3
