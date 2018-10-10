#!/usr/bin/env python
import sys
from presto import rfifind

if len(sys.argv) != 2:
    print("\nusage: {} file\n".format(sys.argv[0]))
    sys.exit(1)



if __name__=="__main__":
    a = rfifind.rfifind(sys.argv[1])
    sys.stderr.write("\nWARNING!:  If raw data have channels in decreasing freq\n")
    sys.stderr.write("           order, the channel ordering as given will be\n")
    sys.stderr.write("           inverted!  Use 'invertband=True' in \n")
    sys.stderr.write("           write_weights() in that case!\n")
    if (a.idata.telescope=='GBT' and a.idata.lofreq < 1000.0):
        sys.stderr.write("Data is from GBT Prime Focus, auto-flipping the weights/offsets...\n\n")
        invert = True
    else:
        invert = False
    # Write the bandpass before we zap things
    a.write_bandpass(invertband=invert)
    # Now do the zapping and set the weights
    a.set_zap_chans(power=200.0,
                    edges=0.01,
                    asigma=2.0,
                    ssigma=2.0,
                    usemask=True,
                    plot=True,
                    chans=[])
    a.write_zap_chans()
    a.set_weights_and_offsets()
    a.write_weights(invertband=invert)
    #a.write_weights_and_offsets(invertband=invert)
