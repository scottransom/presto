#!/usr/bin/env python
import sys
import glob
import argparse
from presto import rfifind


def find_raw_files(basename):
    """Best-effort glob for the raw PSRFITS or SIGPROC filterbank file(s)
    that produced this rfifind basename, assuming the typical
    `rfifind -o <name> <name>*fits` (or `*.fil`) convention (rfifind
    appends '_rfifind' to whatever -o was given)."""
    stem = basename
    if stem.endswith("_rfifind"):
        stem = stem[: -len("_rfifind")]
    candidates = sorted(glob.glob(stem + "*.fits") + glob.glob(stem + "*.fil"))
    return [c for c in candidates if "_subs" not in c]


def raw_band_is_inverted(rawfile):
    """Check a raw PSRFITS or SIGPROC filterbank file directly for
    decreasing-frequency channel order (CHAN_BW < 0 for PSRFITS, or
    foff < 0 for filterbank -- both equivalent to the channel frequencies
    running high-to-low). This is ground truth straight from the header,
    unlike guessing from the backend name -- e.g. GUPPI only writes
    decreasing-frequency data when the receiver in use is GBT Prime Focus;
    on other receivers it writes the band correctly."""
    if rawfile.lower().endswith(".fil"):
        from presto.filterbank import read_header

        header, _ = read_header(rawfile)
        return header["foff"] < 0

    from astropy.io import fits

    with fits.open(rawfile, memmap=True) as f:
        chan_bw = f["SUBINT"].header.get("CHAN_BW")
        if chan_bw is not None:
            return chan_bw < 0
        freqs = f["SUBINT"].data[0]["DAT_FREQ"]
        return bool(freqs[0] > freqs[-1])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "statsfile", type=str, help="rfifind .stats file to compute the channel mask"
    )
    parser.add_argument(
        "--fpower",
        type=float,
        default=200,
        help="Fourier power threshold above which entire channel will be ignored, default = 200",
    )
    parser.add_argument(
        "--band_edge_frac",
        type=float,
        default=0.01,
        help="Fraction of band edge channels to ignore, default = 0.01",
    )
    parser.add_argument(
        "--avgsigma",
        type=float,
        default=2.0,
        help="Channel mean threshold above which entire channel will be ignored, default = 2",
    )
    parser.add_argument(
        "--stdsigma",
        type=float,
        default=2.0,
        help="Channel std threshold above which entire channel will be ignored, default = 2",
    )
    parser.add_argument(
        "--invert",
        default=False,
        dest="invert",
        action="store_true",
        help="Flag to invert weights in case raw data has decreasing freq channels",
    )
    parser.add_argument(
        "--noxwin",
        default=False,
        dest="noxwin",
        action="store_true",
        help="Don't show the zapped-bandpass plot on screen (/xwin)",
    )
    parser.add_argument(
        "--rawfile",
        type=str,
        default=None,
        help="Raw PSRFITS file (or glob pattern) to check directly for true "
        "band inversion, instead of guessing from the backend name. "
        "Defaults to auto-locating a matching raw file next to the rfifind "
        "output (basename with '_rfifind' stripped).",
    )
    args = parser.parse_args()

    a = rfifind.rfifind(args.statsfile)
    sys.stderr.write("\nWARNING!:  If raw data have channels in decreasing freq\n")
    sys.stderr.write("           order, the channel ordering as given will be\n")
    sys.stderr.write("           inverted!  Use 'invertband=True' in \n")
    sys.stderr.write("           write_weights() in that case!\n")

    raw_files = (
        sorted(glob.glob(args.rawfile)) if args.rawfile else find_raw_files(a.basename)
    )
    raw_inverted = None
    if raw_files:
        try:
            raw_inverted = raw_band_is_inverted(raw_files[0])
            sys.stderr.write(
                "Checked raw file '%s' directly: band is %s.\n\n"
                % (raw_files[0], "INVERTED" if raw_inverted else "NOT inverted")
            )
        except Exception as e:
            sys.stderr.write(
                "Could not inspect raw file '%s' (%s).\n"
                "Falling back to the backend-name/Prime-Focus heuristic.\n\n"
                % (raw_files[0], e)
            )

    decreasing_freq_backend = a.idata.instrument.strip().upper() in ("VEGAS", "GUPPI")

    if raw_inverted is not None:
        # Ground truth from the raw file itself -- authoritative, so this
        # takes priority over the backend-name/Prime-Focus guesses below.
        invert = args.invert
        write_dual_output = raw_inverted and not invert
        if raw_inverted and not invert:
            sys.stderr.write(
                "Raw data is band-inverted. Standard (non-inverted) output files\n"
                "will be written as usual, plus a second, band-inverted set of\n"
                "the bandpass/weights files for use with psrfits_subband.\n\n"
            )
    elif a.idata.telescope == "GBT" and a.idata.lofreq < 1000.0:
        sys.stderr.write(
            "No raw file found to check directly; falling back to the GBT\n"
            "Prime Focus heuristic, auto-flipping the weights/offsets...\n\n"
        )
        invert = True
        write_dual_output = False
    elif decreasing_freq_backend:
        sys.stderr.write(
            "No raw file found to check directly; falling back to guessing\n"
            "from the backend name (%s), which writes raw channels in\n"
            "decreasing frequency order for some receivers (e.g. GBT Prime\n"
            "Focus). Standard (non-inverted) output files will be written as\n"
            "usual, plus a second, band-inverted set for psrfits_subband.\n\n"
            % a.idata.instrument.strip()
        )
        invert = args.invert
        write_dual_output = decreasing_freq_backend and not invert
    else:
        invert = args.invert
        write_dual_output = False

    # Write the bandpass before we zap things
    a.write_bandpass(invertband=invert)
    # Now do the zapping and set the weights
    a.set_zap_chans(
        power=args.fpower,
        edges=args.band_edge_frac,
        asigma=args.avgsigma,
        ssigma=args.stdsigma,
        usemask=True,
        plot=not args.noxwin,
        chans=[],
    )
    a.write_zap_chans()
    a.set_weights_and_offsets()
    a.write_weights(invertband=invert)
    # a.write_weights_and_offsets(invertband=invert)

    if write_dual_output:
        # Also write a band-inverted set for psrfits_subband, which (unlike
        # rfifind's own downstream PRESTO tools) applies the .weights file to
        # the raw channels with no frequency-order check of its own.
        invert_bandpass_file = a.basename + "_inverted_for_psrfits.bandpass"
        invert_weights_file = a.basename + "_inverted_for_psrfits.weights"
        a.write_bandpass(filename=invert_bandpass_file, invertband=True)
        a.write_weights(filename=invert_weights_file, invertband=True)
        sys.stderr.write(
            "Wrote band-inverted bandpass/weights for psrfits_subband:\n"
            "  %s\n  %s\n"
            "Pass these to psrfits_subband's -bandpass/-weights options.\n\n"
            % (invert_bandpass_file, invert_weights_file)
        )
