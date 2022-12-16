#!/usr/bin/env python
import sys
import argparse
from presto import rfifind

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("statsfile", type=str, help="rfifind .stats file to compute the channel mask")
    parser.add_argument("--fpower", type=float,
                        default=200,
                        help="Fourier power threshold above which entire channel will be ignored, default = 200")
    parser.add_argument("--band_edge_frac", type=float,
                        default=0.01, help="Fraction of band edge channels to ignore, default = 0.01")
    parser.add_argument("--avgsigma", type=float,
                        default=2.,
                        help="Channel mean threshold above which entire channel will be ignored, default = 2")
    parser.add_argument("--stdsigma", type=float,
                        default=2.,
                        help="Channel std threshold above which entire channel will be ignored, default = 2")
    parser.add_argument("--invert", default=False, dest="invert", action="store_true",
                        help="Flag to invert weights in case raw data has decreasing freq channels")
    args = parser.parse_args()

    a = rfifind.rfifind(args.statsfile)
    sys.stderr.write("\nWARNING!:  If raw data have channels in decreasing freq\n")
    sys.stderr.write("           order, the channel ordering as given will be\n")
    sys.stderr.write("           inverted!  Use 'invertband=True' in \n")
    sys.stderr.write("           write_weights() in that case!\n")
    if (a.idata.telescope == 'GBT' and a.idata.lofreq < 1000.0):
        sys.stderr.write("Data is from GBT Prime Focus, auto-flipping the weights/offsets...\n\n")
        invert = True
    else:
        invert = args.invert
    # Write the bandpass before we zap things   
    a.write_bandpass(invertband=invert)
    # Now do the zapping and set the weights
    a.set_zap_chans(power=args.fpower,
                    edges=args.band_edge_frac,
                    asigma=args.avgsigma,
                    ssigma=args.stdsigma,
                    usemask=True,
                    plot=True,
                    chans=[])
    a.write_zap_chans()
    a.set_weights_and_offsets()
    a.write_weights(invertband=invert)
    #a.write_weights_and_offsets(invertband=invert)
