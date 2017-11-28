#include "presto.h"

#define MEDIANBINS  200

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

void zapbirds(double lobin, double hibin, FILE * fftfile, fcomplex * fft)
{
    long ii, ilobin, ihibin, binstozap, lodatabin;
    float median_lo, median_hi, avgamp;
    /* double phase, radargr, radargi, radtmp; */
    fcomplex *data;

    ilobin = (long) floor(lobin);
    ihibin = (long) ceil(hibin);
    binstozap = ihibin - ilobin;
    if (lobin - 1.5 * MEDIANBINS > 1) {
        if (fftfile) {          /* If we are reading a file */
            data = get_rawbins(fftfile, lobin - MEDIANBINS,
                               MEDIANBINS, &median_lo, &lodatabin);
            vect_free(data);
            data = get_rawbins(fftfile, hibin, MEDIANBINS, &median_hi, &lodatabin);
            vect_free(data);
        } else {                /* If we are working from memory */
            lodatabin = lobin - 3 * MEDIANBINS / 2;
            median_lo = calc_median_powers(fft + lodatabin, MEDIANBINS);
            lodatabin = hibin - MEDIANBINS / 2;
            median_hi = calc_median_powers(fft + lodatabin, MEDIANBINS);
        }
        avgamp = sqrt(0.5 * (median_lo + median_hi) / -log(0.5));
        data = gen_cvect(binstozap);
        /* Read the data to zap */
        if (fftfile) {          /* If we are reading a file */
            chkfileseek(fftfile, ilobin, sizeof(fcomplex), SEEK_SET);
            chkfread(data, sizeof(fcomplex), binstozap, fftfile);
        } else {                /* If we are working from memory */
            data = fft + ilobin;
        }
        for (ii = 0; ii < binstozap; ii++) {
            /* Change the amplitudes but not the phases */
            /*
               phase = RADIAN_PHASE(data[ii].r, data[ii].i);
               data[ii].r = avgamp * cos(phase);
               data[ii].i = avgamp * sin(phase);
             */
            /* Just set the amplitudes to the avgvalue */
            data[ii].r = avgamp;
            data[ii].i = 0.0;
        }
        /* Write the modified data */
        if (fftfile) {          /* If we are reading a file */
            chkfileseek(fftfile, ilobin, sizeof(fcomplex), SEEK_SET);
            chkfwrite(data, sizeof(fcomplex), binstozap, fftfile);
            vect_free(data);
            //printf("Set bins %9d through %9d to amplitude of %.3g\n",
            //       ilobin, ihibin, avgamp);
        }
    }
}
