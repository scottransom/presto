#!/usr/bin/env python
import sys
import os
import argparse
import numpy as np
import presto.presto as pp
import presto.psr_utils as pu
import presto.infodata as pi
import presto.harmonic_sum as ph
from pathlib import Path
from typing import Union


class fftfile:
    "A PRESTO FFT file (i.e. with suffix '.fft') and associated metadata"

    def __init__(self, ff: Union[str, os.PathLike]) -> None:
        """Initialize a PRESTO fftfile class instance.

        Parameters
        ----------
        ff : file or str or Path
            The PRESTO .fft file to open
        """
        self.ff: os.PathLike = ff if isinstance(ff, os.PathLike) else Path(ff)
        self.amps = np.memmap(self.ff, dtype=np.complex64)
        self.inf = pi.infodata(f"{str(self.ff)[:-4]}.inf")
        self.N: int = self.inf.N
        self.T: float = self.N * self.inf.dt
        self.dereddened = True if "_red.fft" in str(self.ff) else False
        self.detrended = True if self.dereddened else False
        self.DC, self.Nyquist = self.amps[0].real, self.amps[0].imag
        self.df: float = 1.0 / self.T
        self._powers = None
        self._phases = None

    @property
    def powers(self) -> np.ndarray:
        """The powers for the FFT amplitudes."""
        if self._powers is None:
            self._powers = pp.spectralpower(self.amps)
            # self._powers[0] = self.DC**2 # this throws away the Nyquist power
        return self._powers

    @property
    def phases(self) -> np.ndarray:
        """The phases (in degrees) for the FFT amplitudes."""
        if self._phases is None:
            self._phases = pp.spectralphase(self.amps)
            # self._phases[0] = 0.0 # this throws away the Nyquist phase (which is also 0)
        return self._phases

    @property
    def freqs(self) -> np.ndarray:
        """The frequencies (in Hz) for the FFT amplitudes."""
        self._freqs = np.linspace(0.0, self.N // 2 * self.df, self.N // 2)
        return self._freqs

    def rednoise_normalize(self, b0=50, bmax=10000, get_medians=False):
        """Normalise the power spectrum while removing rednoise. Based on PRESTO's method
        of rednoise removal, in which a logarithmically increasing window is used at low
        frequencies to compute the local median. The median value to divide for each
        frequency bin is identified by a linear fit between adjacent local medians.

        Parameters
        ----------
        b0: int
            The size of the first window to normalize the power spectrum. Default = 50

        bmax: int
            The maximum size of the largest window to normalize the power spectrum. Default = 10000

        get_medians: bool
            Whether to get the medians out of the rednoise normalization process. Default = False

        Returns
        -------
        normalized_power_spectrum: np.ndarray
            An ndarray of the normalized power spectrum with rednoise removed
        """
        power_spectrum = self.powers
        ps_len = self.N // 2
        scale = []
        for n in range(0, ps_len):
            # create the log range for normalization
            new_window = np.exp(1 + n / 3) * b0 / np.exp(1)
            if new_window > bmax:
                pass
            else:
                window = int(new_window)
            scale.append(window)
            if np.sum(scale) > ps_len:
                scale[-1] = ps_len - np.sum(scale[:-1])
                break
        # check if sum of scale is equal to ps_len
        if np.sum(scale) < ps_len:
            scale[-1] += ps_len - np.sum(scale)
        start = 0
        old_mid_bin = 0
        old_median = 1
        normalized_power_spectrum = np.zeros(shape=np.shape(power_spectrum))
        medians = []

        for bins in scale:
            mid_bin = int(start + bins / 2)
            new_median = np.nanmedian(power_spectrum[start : start + bins])
            medians.append(new_median)
            i = 0
            while np.isnan(new_median):
                i += 1
                new_median = np.nanmedian(
                    power_spectrum[start + (i * bins) : start + ((i + 1) * bins)]
                )
                if not np.isnan(new_median):
                    if start == 0:
                        new_median = new_median * (2**i)
                    else:
                        computed_median = new_median + (
                            (old_median - new_median) * (2**i / (2 ** (i + 1) - 1))
                        )
                        if computed_median - new_median < 0:
                            new_median = old_median
                        else:
                            new_median = computed_median
            if start == 0:
                normalized_power_spectrum[start : start + mid_bin] = power_spectrum[
                    start : start + mid_bin
                ] / (new_median / np.log(2))
            elif start + bins >= ps_len:
                median_slope = np.linspace(
                    old_median, new_median, num=ps_len - old_mid_bin
                )
                normalized_power_spectrum[old_mid_bin:] = power_spectrum[
                    old_mid_bin:
                ] / (median_slope / np.log(2))
            else:
                # compute slope of the power spectra
                median_slope = np.linspace(
                    old_median, new_median, num=mid_bin - old_mid_bin
                )
                normalized_power_spectrum[old_mid_bin:mid_bin] = power_spectrum[
                    old_mid_bin:mid_bin
                ] / (median_slope / np.log(2))

            start += bins
            old_mid_bin = mid_bin
            old_median = new_median

        if get_medians:
            return normalized_power_spectrum, medians, scale
        else:
            return normalized_power_spectrum


class stack:
    "A PRESTO stack of FFT power spectra, normalized as chi^2/2"

    def __init__(self, ffts: list) -> None:
        """Initialize a PRESTO FFT power spectra stack class instance.

        Parameters
        ----------
        ffts : list of file or str or Path
            The PRESTO .fft files to open and add to the stack
        """
        self.ffts = [ff if isinstance(ff, os.PathLike) else Path(ff) for ff in ffts]
        self.nstacked = 0
        self.nharms = 1
        self.stack = None
        self.hstack = None
        self.filenms = []
        for fft in self.ffts:
            self.add_to_stack(fft)

    def add_to_stack(self, ff: Union[str, os.PathLike], blocklen: int = 1000):
        """Add a PRESTO FFT to a power spectrum stack class instance.

        Parameters
        ----------
        ff : file or str or Path
            The PRESTO .fft file to open and add to the stack
        blocklen : int
            The block length for linear detrending of the powers (if detrending is needed)
        """
        ft = fftfile(ff)
        if self.stack is None:
            self.N = ft.N
            self.freqs = ft.freqs
            self.df = ft.df
            self.filenms.append(str(ft.ff))
            self.stack = np.zeros(ft.N // 2)
        else:
            assert self.N == ft.N
        self.stack += ft.powers if ft.detrended else ft.rednoise_normalize()
        self.nstacked += 1

    def get_stats(self, lobin=10):
        """Return the mean, median, and standard deviation of the stack

        Parameters
        ----------
        lobin : int
            The lowest frequency bin to include in the stats (useful for avoiding rednoise)
        """
        if self.stack is not None:
            return (self.stack[lobin:].mean(),
                    np.median(self.stack[lobin:]),
                    self.stack[lobin:].std())
        else:
            return (None, None, None)

    def show_chi2_stack_stats(self, lobin=10):
        """Show the converted stats of the stack for the proper chi^2 distribution

        Parameters
        ----------
        lobin : int
            The lowest frequency bin to include in the stats (useful for avoiding rednoise)
        """
        if self.stack is not None:
            mean, med, std = self.get_stats(lobin=lobin)
            mean *= 2 # since we normalize a single power spectrum to 1
            var = (2 * std)**2 
            print(f"Mean     = {mean:.3f} (expect {2 * self.nstacked})")
            print(f"Variance = {var:.3f} (expect {4 * self.nstacked})")

    def show_chi2_hstack_stats(self, lobin=10):
        """Show the converted stats of the harmonic stack for the proper chi^2 distribution

        Parameters
        ----------
        lobin : int
            The lowest frequency bin to include in the stats (useful for avoiding rednoise)
        """
        if self.hstack is not None:
            mean, std = (self.hstack[lobin:].mean(), self.hstack[lobin:].std())
            mean *= 2 # since we normalize a single power spectrum to 1
            var = (2 * std)**2 
            print(f"Mean     = {mean:.3f} (expect {2 * self.nstacked * self.nharms})")
            print(f"Variance = {var:.3f} (expect {4 * self.nstacked * self.nharms})")

    def sum_next_harmonics(self):
        """Generate a stack with summed harmonics based on the next power of two"""
        pstack = self.stack if self.nharms==1 else self.hstack
        self.hstack = ph.harmonic_sum(self.nharms*2, self.stack, partial=pstack, partialN=self.nharms)
        self.nharms *= 2
