#include "backend_common.h"
#include "vectors.h"

void remove_zerodm(float *fdata, struct spectra_info *s)
// Remove the channel-weighted zero-DM from the raw data
// This implements a bandpass-weighted version of zero-DMing
// from Eatough, Keane & Lyne 2009
{
    int ii, jj;
    float zerodm, padvalsum = 0.0, *chanwts, *fptr;
    static int firsttime = 1;
    static float *bandpass = NULL;

    // Make a static copy of the bandpass to use as DC offsets for each
    // channel. This is necessary when masking as the newpadvals are
    // used to mask. So our post-zero-DM-subtracted band had better look
    // kinda like newpadvals. We want it static so that we don't get a
    // different DC offset for each row of the FITS file (which will put
    // a nice periodicity in the data).
    if (firsttime) {
        float bpsum = 0.0;
        bandpass = gen_fvect(s->num_channels);
        for (ii = 0; ii < s->num_channels; ii++) {
            bandpass[ii] = s->padvals[ii];
            bpsum += bandpass[ii];
        }
        if (bpsum == 0.0) {     // i.e. no padding is set
            float favg;
            printf("\nUsing first block channel averages for zeroDM bandpass.\n");
            printf("Would be better to use statistics from an rfifind mask...\n\n");
            for (ii = 0; ii < s->num_channels; ii++) {
                favg = 0.0;
                for (jj = 0; jj < s->spectra_per_subint; jj++)
                    favg += fdata[ii + jj * s->num_channels];
                bandpass[ii] = favg / s->spectra_per_subint;
            }
        }
    }
    // Determine what weights to use for each channel for
    // the zero-DM subtraction.  Use the padding if available, 
    // otherwise, use the initially-set bandpass
    chanwts = gen_fvect(s->num_channels);
    //if (s->clip_sigma > 0.0 && !firsttime) fptr = s->padvals;
    if (s->clip_sigma > 0.0)
        fptr = s->padvals;
    else
        fptr = bandpass;
    // Determine the sum of the padding/bandpass for weighting of the channels
    for (ii = 0; ii < s->num_channels; ii++)
        padvalsum += *(fptr + ii);
    // Set the channel weights for this block
    for (ii = 0; ii < s->num_channels; ii++)
        chanwts[ii] = *(fptr + ii) / padvalsum;

    // Now loop over all the spectra to do the zero-DMing
    for (ii = 0; ii < s->spectra_per_subint; ii++) {
        // Determine the DM=0 total power for this spectra
        zerodm = 0.0;
        fptr = fdata + ii * s->num_channels;
        for (jj = 0; jj < s->num_channels; jj++)
            zerodm += *fptr++;
        // Subtract off the correct amount of power from each point
        fptr = fdata + ii * s->num_channels;
        for (jj = 0; jj < s->num_channels; jj++) {
            // Put the constant bandpass back in since we are subtracting 
            // comparable sized numbers and don't (usually) want power < 0
            fptr[jj] -= (chanwts[jj] * zerodm - bandpass[jj]);
        }
    }
    vect_free(chanwts);
    firsttime = 0;
}
