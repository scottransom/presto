"""Tests for presto.sifting."""

import warnings

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def test_loglin_scale(tmp_path) -> None:
    """The 'loglin' axis scale registers and draws without warnings, and inverts exactly.

    It is log10(x) + 1 below the break at 1 and linear above it. It raised TypeError on
    matplotlib 3.10 and 3.11 before being fixed, and warned on import with 3.11.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        from presto import sifting  # noqa: F401  (registers the scale)

        fig, ax = plt.subplots()
        ax.plot([0.01, 0.1, 1.0, 10.0, 100.0], [1, 2, 3, 4, 5])
        ax.set_xscale("loglin")
        fig.savefig(tmp_path / "loglin.png")
    t = ax.xaxis.get_transform()
    x = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 10.0])
    fwd = t.transform(x[:, None])[:, 0]
    np.testing.assert_allclose(fwd, [-1.0, 0.0, np.log10(0.5) + 1.0, 1.0, 2.0, 10.0])
    np.testing.assert_allclose(t.inverted().transform(fwd[:, None])[:, 0], x)
    plt.close(fig)
