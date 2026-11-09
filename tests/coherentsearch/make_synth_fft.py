"""Write the small synthetic PRESTO ``.fft`` that the coherentsearch fixtures use.

The output (``synth_cs.fft`` plus ``synth_cs.inf``) is committed beside this script, so
the Julia fixture generator and the C port read the very same bytes; this script only
documents how they were made.  Re-running it with a different NumPy may not reproduce
them bit for bit (NumPy does not promise that across versions), and that is fine:
regenerate the fixtures afterwards with ``gen_fixtures.jl``.

The time series is white Gaussian noise (unit variance) plus two pulsars made of
wrapped-Gaussian pulses:

* a strong, narrow one at Fourier frequency r = 200.37 bins, all of whose 60
  harmonics lie below Nyquist (``PSR_SLOW``);
* a weaker one at r = 2911.13 bins, so that its harmonics above the 11th lie past
  Nyquist (``PSR_FAST``), which exercises the partly-filled harmonic stacks.

Its FFT is normalised the way PRESTO's ``rednoise`` leaves a spectrum (Fourier powers
of mean 1) and packed the PRESTO way: element 0 holds the DC term in its real part and
the Nyquist term in its imaginary part.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

N = 2**16  # time-series points; the .fft holds N/2 complex64 amplitudes (256 KiB)
DT = 5e-4  # s, so T = 32.768 s and Nyquist is 1000 Hz
SEED = 20261109

# (Fourier frequency in bins, FWHM duty cycle, peak amplitude per sample in noise sigma)
PSR_SLOW = (200.37, 0.03, 0.5)
PSR_FAST = (2911.13, 0.08, 0.2)

# PRESTO's .inf layout for artificial data (as makedata writes it), which readinf()
# reads without the telescope, position and band fields of a real observation.
INF_TEMPLATE = """\
 Data file name without suffix          =  {basenm}
 Telescope used                         =  None (Artificial Data Set)
 Number of bins in the time series      =  {N}
 Width of each time series bin (sec)    =  {dt}
 Any breaks in the data? (1 yes, 0 no)  =  0
 Data analyzed by                       =  make_synth_fft.py
 Any additional notes:
    White noise plus two pulsars (r = {r1} and {r2} bins), FFT normalised to unit
    mean power.  Seed {seed}.
"""


def pulse_train(r: float, duty: float, amp: float) -> np.ndarray:
    """Wrapped-Gaussian pulses at Fourier frequency ``r`` (cycles per N samples).

    Parameters
    ----------
    r : float
        Pulse frequency in Fourier bins (cycles over the whole time series).
    duty : float
        Pulse FWHM as a fraction of the period.
    amp : float
        Peak amplitude per sample, in units of the noise sigma.

    Returns
    -------
    numpy.ndarray
        The ``N``-sample signal, float64.
    """
    phase = (r * np.arange(N) / N) % 1.0
    dphi = np.minimum(phase, 1.0 - phase)  # distance from the pulse peak at phase 0
    sigma = duty / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amp * np.exp(-0.5 * (dphi / sigma) ** 2)


def make_amplitudes() -> np.ndarray:
    """Return the packed, normalised ``N/2`` complex64 amplitudes.

    Returns
    -------
    numpy.ndarray
        complex64 array of length ``N // 2``, PRESTO ``.fft`` layout.
    """
    rng = np.random.Generator(np.random.PCG64(SEED))
    x = rng.standard_normal(N)
    x += pulse_train(*PSR_SLOW)
    x += pulse_train(*PSR_FAST)
    X = np.fft.rfft(x) / np.sqrt(N)  # unit-variance noise -> Fourier powers of mean 1
    amps = X[: N // 2].copy()
    amps[0] = complex(X[0].real, X[N // 2].real)
    return amps.astype(np.complex64)


def main() -> None:
    """Write ``synth_cs.fft`` and ``synth_cs.inf`` into the output directory."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "outdir",
        nargs="?",
        default=str(Path(__file__).resolve().parent),
        help="where to write synth_cs.{fft,inf} (default: beside this script)",
    )
    args = parser.parse_args()
    out = Path(args.outdir)
    basenm = "synth_cs"
    make_amplitudes().astype("<c8").tofile(out / f"{basenm}.fft")
    (out / f"{basenm}.inf").write_text(
        INF_TEMPLATE.format(
            basenm=basenm, N=N, dt=DT, r1=PSR_SLOW[0], r2=PSR_FAST[0], seed=SEED
        )
    )


if __name__ == "__main__":
    main()
