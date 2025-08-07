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
        """Normalize the power spectrum while removing rednoise. Based on PRESTO's method
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
            scale.append(window)  # type: ignore
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
        self.filenms = []
        for fft in self.ffts:
            self.add_to_stack(fft)

    def add_to_stack(
        self, ff: Union[str, os.PathLike], detrended: bool = False, blocklen: int = 1000
    ):
        """Add a PRESTO FFT to a power spectrum stack class instance.

        Parameters
        ----------
        ff : file or str or Path
            The PRESTO .fft file to open and add to the stack
        detrended : bool
            Whether the FFT has already been detrended or not (default is False)
        blocklen : int
            The block length for linear detrending of the powers (if detrending is needed)
        """
        ft = fftfile(ff)
        ft.detrended = True if detrended else ft.detrended
        fname = str(ft.ff)
        print(f"Adding '{fname}' to stack.")
        if not hasattr(self, "stack"):
            self.N = ft.N
            self.freqs = ft.freqs
            self.T = ft.T
            self.df = ft.df
            self.filenms.append(fname)
            self.stack = np.zeros(ft.N // 2)
        else:
            assert np.isclose(
                self.df, ft.df
            ), f"{fname} has freq spacing {ft.df:.12g} rather than {self.df:.12g}"
        if self.N == ft.N:
            self.stack += ft.powers if ft.detrended else ft.rednoise_normalize()
        elif self.N < ft.N:  # truncate
            print(f"{fname} is shorter than current stack. Truncating.")
            self.stack += (
                ft.powers[: self.N]
                if ft.detrended
                else ft.rednoise_normalize()[: self.N]
            )
        else:  # pad
            print(f"{fname} is longer than current stack. Padding.")
            self.stack[: ft.N] += ft.powers if ft.detrended else ft.rednoise_normalize()
            self.stack[ft.N :] += 1.0
        self.nstacked += 1

    def expected_stack_stats(self, hstack=False) -> tuple[float, float, float]:
        """The mean, median, and stdev expected in the stack or hstack from chi^2 stats

        Parameters
        ----------
        hstack : boolean
            Use hstack (if True) or stack (if False)
        """
        DOF = 2 * self.nstacked * self.nharms if hstack else 2 * self.nstacked
        # The divides by two below are because we normalize powers to mean=std=1 for
        # a single power spectrum, which is off by a factor of 2 for pure chi^2
        return (DOF / 2, DOF / 2 * (1 - 2 / (9 * DOF)) ** 3, np.sqrt(DOF / 2))

    def get_stats(self, lobin=100, hstack=False):
        """Return the mean, median, and standard deviation of the stack or hstack

        Parameters
        ----------
        lobin : int
            The lowest frequency bin to include in the stats (useful for avoiding rednoise)
        hstack : boolean
            Use hstack (if True) or stack (if False)
        """
        tmpstack: np.ndarray = self.hstack if hstack else self.stack  # type: ignore
        return (
            tmpstack[lobin:].mean(),
            np.median(tmpstack[lobin:]),
            tmpstack[lobin:].std(),
        )

    def show_stats(self, lobin=100, hstack=False):
        """Show the stats of the stack or hstack as compared to what is expected

        Parameters
        ----------
        lobin : int
            The lowest frequency bin to include in the stats (useful for avoiding rednoise)
        """
        mean, med, std = self.get_stats(lobin=lobin, hstack=hstack)
        xmean, xmed, xstd = self.expected_stack_stats(hstack=hstack)
        extra = f" and nharms={self.nharms}:" if hstack else ":"
        print(
            f"For {"harmonic " if hstack else ""}stack with nstacked={self.nstacked}{extra}"
        )
        print(f"  Mean   = {mean:7.3f} (expect {xmean:7.3f})")
        print(f"  Median = {med:7.3f} (expect {xmed:7.3f})")
        print(f"  StdDev = {std:7.3f} (expect {xstd:7.3f})")

    def sum_next_harmonics(self):
        """Generate a stack with summed harmonics based on the next power of two"""
        pstack = self.stack if self.nharms == 1 else self.hstack
        self.hstack = ph.harmonic_sum(
            self.nharms * 2, self.stack, partial=pstack, partialN=self.nharms
        )
        self.nharms *= 2

    def power_threshold_from_sigma(self, sigma: float = 5.0) -> float:
        """Return the threshold power required for the harmonic stack to exceed sigma"""
        return pu.powersum_at_sigma(sigma, self.nstacked * self.nharms)

    def search_hstack(self, pthresh: float, lobin: int = 50) -> np.ndarray:
        """Return a sorted array of the hstack indices exceeding pthresh (most->least significant)

        Parameters
        ----------
        pthresh : float
            The summed power threshold required to return a candidate
                (usually set by using the power_threshold_from_sigma() method)
        lobin : int
            The lowest bin to search for significant signals
        """
        tmpstack: np.ndarray = self.hstack if self.nharms > 1 else self.stack  # type: ignore
        inds = np.arange(self.N // 2)
        # following is a boolean array for the bins above pthresh and with index above lobin
        good = np.logical_and(inds > lobin, tmpstack > pthresh)
        # Select the powers that made the cut
        goodpows = tmpstack[good]
        # And figure out what their indices are
        goodinds = inds[good]
        # return the indices with the highest powers first
        return goodinds[np.argsort(goodpows)[::-1]]

    def calc_freqs(self, bins: np.ndarray) -> np.ndarray:
        """Return the harmonic-adjusted (i.e. fundamental) frequencies (Hz) for the bins

        Parameters
        ----------
        bins : np.ndarray
            An array of indices in the hstack to compute fundamental frequencies for
        """
        return self.freqs[bins] / self.nharms


def mycmp(a, b):
    return (a > b) - (a < b)


class candidate:
    "An individual stack candidate."

    def __init__(self, bin, freq, nharm, nstack, power):
        self.bin = bin
        self.freq = freq
        self.nharm = nharm
        self.power = power
        self.sigma = pp.candidate_sigma(power, nharm * nstack, 1)

    def __str__(self):
        return f" {self.sigma:7.2f} {self.freq:13.6f} {self.bin:13.3f} {self.power:8.2f} {self.nharm:5d}"

    def __eq__(self, other):
        return self.sigma == other.sigma

    def __ne__(self, other):
        return self.sigma != other.sigma

    def __lt__(self, other):
        return self.sigma < other.sigma

    def __le__(self, other):
        return self.sigma <= other.sigma

    def __gt__(self, other):
        return self.sigma > other.sigma

    def __ge__(self, other):
        return self.sigma >= other.sigma

    def __cmp__(self, other):
        # Sort by sigma by default)
        return mycmp(self.sigma, other.sigma)


class stackcands:
    "Candidates from a stack search of power spectra."

    def __init__(self, stack: stack, indices: np.ndarray, powers: np.ndarray) -> None:
        """Initialize a stackcands class which will contain a list of stack cands.

        Parameters
        ----------
        stack : stack
            The stack class instance for the search.
        indices : np.ndarray
            The Fourier indices which were selected.
        powers : np.ndarray
            The corresponding stack powers for the candidate indices.
        """
        self.cands = []
        self.add_candidates(stack, indices, powers)

    def add_candidates(
        self, stack: stack, indices: np.ndarray, powers: np.ndarray
    ) -> None:
        """Add more candidates to the current stack candidates list.

        Parameters
        ----------
        stack : stack
            The stack class instance for the search.
        indices : np.ndarray
            The Fourier indices which were selected.
        powers : np.ndarray
            The corresponding stack powers for the candidate indices.
        """
        for ind, pow in zip(indices, powers):
            self.cands.append(
                candidate(ind, ind / stack.T / stack.nharms, stack.nharms, stack.nstacked, pow)
            )

    def output_candidates(self, outfile=None, maxncands=100) -> None:
        """Write text-formatted stack candidates to stdout or to a file.

        Parameters
        ----------
        outfile : string, optional
            Output file name, by default None
        maxncands : int, optional
            Max number of candidates to output, by default 100
        """
        self.cands = sorted(self.cands, reverse=True)
        if outfile is None:
            out = sys.stdout
        else:
            out = open(outfile, "w")
        out.write(
            f"# {"Sigma":^7} {"Freq (Hz)":^13} {"Fourier Bin":^13} {"Power":^8} {"#Harm":^5}\n"
        )
        out.write(f"#{'-'*(7+13+13+8+5+5)}\n")
        for ii, cand in enumerate(self.cands):
            if ii > maxncands:
                break
            out.write(str(cand)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read multiple PRESTO-style *.fft files and conduct a stack search for periodicities.",
        epilog="""In general, the FFT files should probably be barycentered and have known
RFI zapped (i.e. barycentering happens by default if you use `prepdata` or
`prepsubband` and zapping can be done using, for instance, simple_zapbirds.py`).
The sigma threshold is single-trial and based on equivalent gaussian sigma.
If no output candidate file name is given, the results will be written to stdout.
""",
    )
    parser.add_argument("fftfiles", nargs="*", help="PRESTO .fft files to be stacked.")
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=8,
        help="Sigma cutoff for picking candidates (default=8)",
    )
    parser.add_argument(
        "-o", "--outputfilenm", type=str, help="Output filename to record candidates"
    )
    parser.add_argument(
        "-c",
        "--ncands",
        type=int,
        default=100,
        help="Maximum number of candidates to return for each harmonic fold (default=100)",
    )
    parser.add_argument(
        "-l",
        "--lobin",
        type=int,
        default=100,
        help="Lowest frequency bin to search or to use for statistics (to avoid rednoise, default=100)",
    )
    parser.add_argument(
        "-f",
        "--lofreq",
        type=float,
        default=0.1,
        help="Lowest frequency (in Hz) to search or to use for statistics (to avoid rednoise, default=0.1)",
    )
    parser.add_argument(
        "-n",
        "--nharms",
        type=int,
        default=16,
        help="Maximum number of harmonics to sum. A power-of-two. (default=16)",
    )
    args = parser.parse_args()
    if not args.fftfiles:
        parser.print_help()
        sys.exit(1)

    print("PRESTO Stack Search by Scott Ransom (copyright 2025)\n")

    # Create the stack
    ss = stack(args.fftfiles)
    lobin = int(max(args.lobin, args.lofreq * ss.T))
    print(
        f"\nLowest bin to use for searches or stats is {lobin} ({lobin/ss.T:.4f} Hz)\n"
    )
    ss.show_stats(lobin=lobin, hstack=False)

    # Search the stack for single-harmonic candidates
    print("\nSearching the stack with no harmonics summed:")
    pthresh = ss.power_threshold_from_sigma(args.threshold)
    inds = ss.search_hstack(pthresh=pthresh, lobin=lobin)
    cands = stackcands(ss, inds, ss.stack[inds])

    while ss.nharms < pp.next2_to_n(args.nharms):
        ss.sum_next_harmonics()
        print(f"\nSearching the stack with {ss.nharms:2d} harmonics summed:", end=" ")
        pthresh = ss.power_threshold_from_sigma(args.threshold)
        inds = ss.search_hstack(pthresh=pthresh, lobin=lobin)
        print(f"({len(inds)} cands)")
        cands.add_candidates(ss, inds, ss.hstack[inds])  # type: ignore

    # Now output the candidates
    if args.outputfilenm is None:
        print("")
    cands.output_candidates(outfile=args.outputfilenm, maxncands=args.ncands)
