"""Generate synthetic SIGPROC filterbank files containing a dispersed pulsar.

Used by the end-to-end pipeline tests. Everything is deterministic for a given seed, and the
injected signal is specified by its true spin frequency, DM, duty cycle and S/N, so tests can
check what the tools find against the known answer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from presto import filterbank

KDM = 4.148808e3  # dispersion constant, s MHz^2 / (pc cm^-3)


@dataclass(frozen=True)
class FakePulsar:
    """Parameters of a synthetic observation and the pulsar injected into it.

    Attributes
    ----------
    f0 : float
        Topocentric spin frequency (Hz), constant over the observation.
    dm : float
        Dispersion measure (pc cm^-3).
    duty : float
        Pulse FWHM as a fraction of the period. The pulse is a wrapped Gaussian.
    snr : float
        Matched-filter S/N of the fully dedispersed, folded pulse, i.e.
        sqrt(sum_t S(t)^2) / sigma for the channel-summed signal S(t) and its noise
        sigma. Intra-channel smearing and quantization lower what a code can recover
        slightly.
    nchan, nspec : int
        Number of channels and of spectra (time samples).
    tsamp : float
        Sampling time (s).
    fch1, foff : float
        Centre frequency of the first channel and the channel spacing (MHz). A
        negative ``foff`` gives the usual descending SIGPROC order.
    tstart : float
        Start MJD.
    seed : int
        Seed for the noise.
    """

    f0: float = 23.4567
    dm: float = 50.0
    duty: float = 0.10
    snr: float = 12.0
    nchan: int = 64
    nspec: int = 2**17
    tsamp: float = 5e-4
    fch1: float = 1498.4375
    foff: float = -3.125
    tstart: float = 60000.0
    seed: int = 12345

    @property
    def T(self) -> float:
        """Observation length (s)."""
        return self.nspec * self.tsamp

    @property
    def freqs(self) -> np.ndarray:
        """Channel centre frequencies (MHz), in file order."""
        return self.fch1 + self.foff * np.arange(self.nchan)


def pulse_signal(psr: FakePulsar) -> np.ndarray:
    """Return the dispersed pulse train, shape ``(nspec, nchan)``, scaled to the requested S/N.

    Units are those of the per-channel noise (unit variance). Delays are relative to
    infinite frequency.

    Parameters
    ----------
    psr : FakePulsar
        The observation and pulsar parameters.

    Returns
    -------
    numpy.ndarray
        The float32 signal to add to unit-variance noise.
    """
    t = np.arange(psr.nspec) * psr.tsamp
    sig_phase = psr.duty / np.sqrt(8.0 * np.log(2.0))  # Gaussian sigma in phase
    delays = KDM * psr.dm / psr.freqs**2
    sig = np.empty((psr.nspec, psr.nchan), dtype=np.float32)
    for ii, delay in enumerate(delays):
        dphi = np.mod(psr.f0 * (t - delay) + 0.5, 1.0) - 0.5
        sig[:, ii] = np.exp(-0.5 * (dphi / sig_phase) ** 2)
    # Dedispersed and summed over channels, the signal is nchan * a * g(t) with noise
    # sigma sqrt(nchan), so snr = a * sqrt(nchan * sum_t g(t)^2).
    a = psr.snr / np.sqrt(psr.nchan * np.sum(sig[:, 0].astype(np.float64) ** 2))
    sig *= a
    return sig


def write_filterbank(
    psr: FakePulsar,
    path: str,
    nbits: int = 8,
    signed: bool = False,
    offset: float | None = None,
    scale: float = 12.0,
) -> None:
    """Write noise plus the pulsar to a SIGPROC filterbank file.

    Parameters
    ----------
    psr : FakePulsar
        The observation and pulsar parameters.
    path : str
        Output file name.
    nbits : int, optional
        Bits per sample: 8 or 16 (integers) or 32 (float). Default 8.
    signed : bool, optional
        Write signed 8-bit samples (sets the SIGPROC ``signed`` header key). Default False.
    offset : float, optional
        Mean level in output units. Defaults to the middle of the integer range (0 for
        signed or float data).
    scale : float, optional
        Noise standard deviation in output units. Default 12.
    """
    rng = np.random.default_rng(psr.seed)
    data = rng.standard_normal((psr.nspec, psr.nchan), dtype=np.float32)
    data += pulse_signal(psr)
    if offset is None:
        offset = 0.0 if (signed or nbits == 32) else float(2 ** (nbits - 1))
    data = data * scale + offset
    if nbits == 32:
        out = data.astype(np.float32)
    else:
        lo, hi = (-(2 ** (nbits - 1)), 2 ** (nbits - 1) - 1) if signed else (0, 2**nbits - 1)
        out = np.clip(np.rint(data), lo, hi)
        if signed:
            out = out.astype(np.int8 if nbits == 8 else np.int16)
    header = {
        "telescope_id": 0,  # "Fake"
        "machine_id": 0,  # "FAKE"
        "data_type": 1,
        "rawdatafile": path,
        "source_name": "FAKEPSR",
        "src_raj": 0.0,
        "src_dej": 0.0,
        "tstart": psr.tstart,
        "tsamp": psr.tsamp,
        "fch1": psr.fch1,
        "foff": psr.foff,
        "nchans": psr.nchan,
        "nifs": 1,
    }
    if signed:
        header["signed"] = 1
    # For signed data the writer's cast to an unsigned dtype keeps the same bytes, and the
    # header's "signed" key tells readers how to interpret them.
    filterbank.create_filterbank_file(path, header, spectra=out, nbits=nbits)

