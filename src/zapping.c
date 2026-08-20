#include "presto.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Number of bins used to measure the local power level on each side */
/* of a zapped region.                                               */
#define MEDIANBINS  200
/* Number of bins skipped on either side of the zapped region before */
/* those measurements start.  The bins just outside a birdie are     */
/* contaminated by its spectral leakage, so they are not usable.     */
#define GUARDBINS   (MEDIANBINS / 2)
/* The replacement noise is drawn from a fixed seed so that zapping   */
/* the same FFT with the same zaplist always gives the same result.   */
#define ZAPSEED     1

float calc_median_powers(fcomplex * amplitudes, int numamps)
{
    int ii;
    float *powers, powargr, powargi, med;

    /* Calculate the median power */
    powers = gen_fvect(numamps);
    for (ii = 0; ii < numamps; ii++)
        powers[ii] = POWER(amplitudes[ii].r, amplitudes[ii].i);
    med = median(powers, numamps);
    vect_free(powers);
    return med;
}

fcomplex *get_rawbins(FILE * fftfile, double bin,
                      int numtoget, float *med, long *lobin)
{
    fcomplex *result;
    *lobin = (int) bin - numtoget / 2;
    result = read_fcomplex_file(fftfile, *lobin, numtoget);
    *med = calc_median_powers(result, numtoget);
    return result;
}

static float local_median(FILE * fftfile, fcomplex * fft, long startbin)
/* Return the median power of the MEDIANBINS bins starting at 'startbin' */
{
    float med;

    if (fftfile) {              /* If we are reading a file */
        fcomplex *data;
        long lodatabin;

        data = get_rawbins(fftfile, startbin + MEDIANBINS / 2,
                           MEDIANBINS, &med, &lodatabin);
        vect_free(data);
    } else {                    /* If we are working from memory */
        med = calc_median_powers(fft + startbin, MEDIANBINS);
    }
    return med;
}

void zapbirds(double lobin, double hibin, long numbins, FILE * fftfile,
              fcomplex * fft, int constamp)
/* Replace the Fourier amplitudes between 'lobin' and 'hibin' so that a  */
/* birdie is removed from a power spectrum.  Either 'fftfile' (a file    */
/* opened for update) or 'fft' (an in-core array) must be given, and     */
/* 'numbins' is the total number of complex bins that it holds.          */
/*                                                                      */
/* The local power level is measured on both sides of the zapped region  */
/* (skipping GUARDBINS on each side, since the bins nearest a birdie are */
/* contaminated by it) and the zapped bins are then replaced by noise    */
/* drawn from that level:  chi^2 with 2 degrees of freedom powers, with  */
/* uniformly random phases.  That leaves the region statistically        */
/* indistinguishable from the surrounding spectrum, which matters when   */
/* the region is wide.  If 'constamp' is true, the behavior of PRESTO    */
/* 6.0.1 and earlier is used instead:  every zapped bin is set to a      */
/* constant amplitude equal to the square root of the local mean power,  */
/* with zero phase.                                                     */
{
    static gsl_rng *rng = NULL;
    long ii, ilobin, ihibin, binstozap, lostart, histart;
    int uselo, usehi, numsides = 0;
    float medpow = 0.0, meanpow, avgamp;
    fcomplex *data;

    ilobin = (long) floor(lobin);
    ihibin = (long) ceil(hibin);
    binstozap = ihibin - ilobin;
    if (binstozap < 1)
        return;

    /* Where the local power measurements will be made */

    lostart = ilobin - GUARDBINS - MEDIANBINS;
    histart = ihibin + GUARDBINS;
    uselo = (lostart >= 1);
    usehi = (histart + MEDIANBINS <= numbins);
    if (!uselo && !usehi) {
        printf("\nWarning:  cannot measure the power level around bins "
               "%ld-%ld.  Not zapping them.\n", ilobin, ihibin);
        return;
    }

    /* Measure the local power level outside of the zapped region */

    if (uselo) {
        medpow += local_median(fftfile, fft, lostart);
        numsides++;
    }
    if (usehi) {
        medpow += local_median(fftfile, fft, histart);
        numsides++;
    }
    medpow /= numsides;
    /* The median of chi^2 with 2 DOF powers is ln(2) times their mean */
    meanpow = medpow / -log(0.5);
    avgamp = sqrt(meanpow);

    /* Get a buffer holding the bins to zap */

    if (fftfile) {              /* If we are reading a file */
        data = gen_cvect(binstozap);
    } else {                    /* If we are working from memory */
        data = fft + ilobin;
    }

    /* Replace the amplitudes */

    if (constamp) {             /* The old behavior:  a constant amplitude */
        for (ii = 0; ii < binstozap; ii++) {
            data[ii].r = avgamp;
            data[ii].i = 0.0;
        }
    } else {                    /* Noise with the local statistics */
        double amp, phase;

        if (rng == NULL) {
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, ZAPSEED);
        }
        for (ii = 0; ii < binstozap; ii++) {
            /* Powers that are chi^2 distributed with 2 DOF (i.e. */
            /* exponential) with the local mean power             */
            amp = sqrt(gsl_ran_exponential(rng, meanpow));
            phase = TWOPI * gsl_rng_uniform(rng);
            data[ii].r = amp * cos(phase);
            data[ii].i = amp * sin(phase);
        }
    }

    /* Write the modified data */

    if (fftfile) {              /* If we are reading a file */
        chkfileseek(fftfile, ilobin, sizeof(fcomplex), SEEK_SET);
        chkfwrite(data, sizeof(fcomplex), binstozap, fftfile);
        vect_free(data);
    }
}
