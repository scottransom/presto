#!/usr/bin/env python
import sys
import os
import getopt
import numpy as np
import presto.infodata as pi
import presto.presto as pp
import presto.psr_utils as pu
import matplotlib.pyplot as plt

def get_fourier_prof(fft, rr, zz=0.0, Nbins=None):
    """Generate a pulse profile from the Fourier amplitudes in a .fft file

    Parameters
    ----------
    fft : Numpy array of type complex64 which contains FFT amplitudes
    rr : Fractional Fourier bin number for the fundamental frequency
    zz : Fractional Fourier f-dot for fundamental (defaults to 0.0)
    Nbins : The number of bins to use in the pulse profile. By default None,
            which means that you will get the number corresponding to the
            number of harmonics that fits within the Nyquist frequency, and
            then rounded up to the nearest power-of-two.

    Returns
    -------
    A numpy float array containing the pulse profile.
    """
    if Nbins is None:
        Nbins = min(256, pp.next2_to_n(len(fft) / rr))
    proffft = np.zeros(Nbins // 2 + 1, dtype=np.complex64)
    for ii in range(1, Nbins // 2 + 1):
        nr, nz = ii * rr, ii * zz
        if nr > len(fft):
            break
        rl, im = pp.rz_interp(fft, nr, nz, pp.z_resp_halfwidth(nz, pp.LOWACC))
        proffft[ii] = complex(rl, im)
    return np.fft.irfft(proffft)


def estimate_profile_variance(fft, rr, Nbins=None, Ntrials=10):
    """Estimate the variance of a pulse profile based on the FFT harmonics

    Parameters
    ----------
    fft : Numpy array of type complex64 which contains FFT amplitudes
    rr : Fractional Fourier bin number for the fundamental frequency
    Nbins : The number of bins to use in the pulse profile. By default None,
            which means that you will get the number corresponding to the
            number of harmonics that fits within the Nyquist frequency, and
            then rounded up to the nearest power-of-two.
    Ntrials : int, optional. The number of random nearby frequencies to use
            to statistically compute the variance. By default 20.

    Returns
    -------
    A float estimate of the off-pulse profile variance.
    """
    minoff, maxoff = 10, 30  # in bins
    drs = np.random.uniform(minoff, maxoff, Ntrials)
    drs[::2] *= -1  # make half of the offsets negative
    return np.median([get_fourier_prof(fft, rr + dr, 0.0, Nbins).var() for dr in drs])


def optimize_freq(fft, rr, var, Nbins=None, dr=0.01, maxoff=1.0):
    """Optimize the frequency of a periodic signal in an FFT

    Parameters
    ----------
    fft : Numpy array of type complex64 which contains FFT amplitudes
    rr : Fractional Fourier bin number for the fundamental frequency
    var : The float variance of the profile
    Nbins : The number of bins to use in the pulse profile. By default None,
            which means that you will get the number corresponding to the
            number of harmonics that fits within the Nyquist frequency, and
            then rounded up to the nearest power-of-two.
    dr : float, optional (by default 0.01). The Fourier frequency step-size
            over which to search.
    maxoff : float, optional (by default 1.0). The maximum Fourier frequency
            deviation from rr over which to search.

    Returns
    -------
    A tuple of an array of reduced chi^2 values of the profiles, and a float
    of the highest significance Fourier frequency of the pulsations.
    """
    rs = np.linspace(rr - maxoff, rr + maxoff, int(2 * maxoff / dr) + 1)
    chis = np.zeros(len(rs))
    for ii, r in enumerate(rs):
        prof = get_fourier_prof(fft, r, 0.0, Nbins)
        chis[ii] = pp.compute_chi2(prof, np.median(prof), var) / (len(prof) - 1.0)
    return chis, rs[chis.argmax()]


def profile_for_plot(prof):
    "Return a centered and doubled pulse profile"
    newprof = prof - np.median(prof)
    return np.tile(pu.rotate(newprof, newprof.argmax()), 2)


def usage():
    print(
        """
usage:  fourier_fold.py [options which must include -f or -c] fft_file
  [-h, --help]                          : Display this help
  [-o, --optimize]                      : Optimize the candidates over local nearby frequencies
  [-f freq, --freq=freq]                : Frequency (Hz) to use for profile
  [-c accelfile, --accelfile=accelfile] : A PRESTO accel- or jerk-search candidate file to use.
                                          By default all candidates will be folded.
  [-n accelcand, --accelcand=accelcand  : Only generate this candidate number from candfile (1 offset)

  This program uses the complex amplitudes in a PRESTO .fft file to generate
  pulse profiles without having to do time-domain folding. It either reads in 
  all of the accelsearch candidates from a ".cand" file, or folds a single
  candidate from that file (via -n option) or one specified on the command line
  (via -f option). PNG plot or plots showing the resulting profiles are generated
  and saved. If the optimizing option is requested, then the significance of the 
  candidate is optimized before the plots are generated (along with additional
  plots of the reduced-chi^2 significance of the pulsations.
  
"""
    )


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hof:c:n:",
            ["help", "optimize", "freq=", "accelfile=", "accelcand="],
        )
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    if len(sys.argv) == 1:
        usage()
        sys.exit(2)

    if len(sys.argv) < 3:
        usage()
        sys.exit(2)

    optimize = False
    freq = None
    zz = 0.0
    accelfile = None
    accelcand = None
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--optimize"):
            optimize = True
        if o in ("-f", "--freq"):
            freq = float(a)
            filenm = f"{args[0][:-4]}_{freq:.4f}Hz.png"
        if o in ("-c", "--accelfile"):
            accelfile = a
        if o in ("-n", "--accelcand"):
            accelcand = int(a)

    #fft = np.fromfile(args[0], dtype=np.complex64)
    fft = np.memmap(args[0], dtype=np.complex64, mode='r')
    idata = pi.infodata(args[0][:-4] + ".inf")
    idata.T = idata.N * idata.dt

    if accelfile:
        ncands = os.stat(accelfile).st_size // 88

    if accelcand is not None and accelfile is None:
        print("Error: If you specify accelcand, you must specify an accelfile.")
        usage()
        sys.exit(1)

    if accelcand is not None and accelfile is not None and accelcand > ncands:
        print("Error: accelcand number is more than the candidates in the accelfile.")
        usage()
        sys.exit(1)

    if accelcand is not None and accelfile is not None:
        fprops = pp.fourierprops()
        pp.get_rzw_cand(accelfile, accelcand, fprops)
        freq = fprops.r / idata.T
        zz = fprops.z
        filenm = f"{accelfile}.{accelcand}.png"
        ncands = 1

    # Step over the candidates
    for ii in range(ncands):
        if accelfile is not None and accelcand is None:
            fprops = pp.fourierprops()
            pp.get_rzw_cand(accelfile, ii + 1, fprops)
            freq = fprops.r / idata.T
            zz = fprops.z
            filenm = f"{accelfile}.{ii+1}.png"

        rr = freq * idata.T

        if optimize:
            var = estimate_profile_variance(fft, rr, Nbins=None, Ntrials=10)
            maxoff = 0.5
            dr = 1.0 / min(len(fft)/rr, 100)
            chis, bestr = optimize_freq(fft, rr, var, Nbins=None, dr=dr, maxoff=maxoff)
            prof = get_fourier_prof(fft, bestr, zz, Nbins=None)
            # Now plot it
            fig, (ax1, ax2) = plt.subplots(
                1, 2, sharex=False, sharey=False, layout="constrained", figsize=(8, 4)
            )
            ax1.plot(np.linspace(0.0, 2.0, 2*len(prof)), profile_for_plot(prof))
            ax1.set_xlabel("Pulse phase (bins)")
            ax1.set_ylabel("Relative intensity")
            ax2.plot(np.linspace(-maxoff, maxoff, int(2 * maxoff / dr) + 1), chis)
            ax2.vlines(0.0, chis.min(), chis.max() * 1.05, colors='grey', alpha=0.4)
            ax2.set_xlabel("Fourier frequency offset (bins)")
            ax2.set_ylabel(r"Reduced-$\chi^2$")
            fig.suptitle(f"Best freq = {bestr / idata.T:.12f} Hz  ({bestr:.2f} bins)")
            fig.savefig(filenm, dpi=400)
            plt.close()
        else:
            prof = get_fourier_prof(fft, rr, zz, Nbins=None)
            # Now plot it
            plt.plot(np.linspace(0.0, 2.0, 2*len(prof)), profile_for_plot(prof))
            plt.xlabel("Pulse phase (bins)")
            plt.ylabel("Relative intensity")
            plt.title(f"Freq = {rr / idata.T:.12f} Hz  ({rr:.2f} bins)")
            plt.savefig(filenm, dpi=400)
            plt.close()
