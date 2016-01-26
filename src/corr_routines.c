#include "presto.h"

/* define DEBUGPRINT */

int corr_complex(fcomplex * data, int numdata, presto_datainf datainf,
                 fcomplex * kern, int numkern, presto_datainf kerninf,
                 fcomplex * result, int numresult, int lobin,
                 int numbetween, int kern_half_width, presto_optype optype)
  /* This routine is a general correlation or convolution routine    */
  /* for complex data.  It can perform convolutions or correlations  */
  /* on raw complex data, data that is prepared for a convolution/   */
  /* correlation but not FFTd, or already FFTd data.  The kernel     */
  /* that it uses can also be raw, prepped, or FFTd.  If you call    */
  /* the routine multiple times with either the same kernel or data  */
  /* array, it uses a saved version of the array from the previous   */
  /* call to cut down on many processing steps. The return value     */
  /* tells how many usable (i.e.  non-contaminated) points were      */
  /* returned in the result array (the first value will be that of   */
  /* 'lobin').  This routine will _not_ perform in-place             */
  /* correlations or convolutions (i.e. it ignores those choices     */
  /* for 'optype').                                                  */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points in 'data'.          */
  /*   'datainf' is one of the following that describes the data:    */
  /*              RAW = Normal un-altered complex data.              */
  /*              PREPPED = Data has been padded and spread based    */
  /*                        on 'kern_half_width' and 'numbetween'    */
  /*                        and is ready to be FFTd.                 */
  /*              FFT = Data has already been prepared and FFTd.     */
  /*              SAME = Data is the same as the previous call.      */
  /*                        The routine uses its saved data.         */
  /*   'kern' is the correlation kernel.                             */
  /*   'numkern' is the number of complex points in 'kern'.          */
  /*   'kerninf' is one of the same choices as 'datainf' above.      */
  /*   'result' is the resulting complex array (must already exist). */
  /*   'numresult' is the number of complex points in 'result'.      */
  /*   'lobin' is the lowest fourier bin to convolve/correlate.      */
  /*   'numbetween' is the number of bins to spread the data points. */
  /*   'kern_half_width' is half the width (bins) of the raw kernel. */
  /*   'optype' is either CORR or CONV (correlation or convolution). */
  /* Notes:                                                          */
  /*   If either 'datainf' or 'kerninf' are of type PREPPED or FFT,  */
  /*   then the length of the FFTs used in the correlation/          */
  /*   convolution calculations will be of length 'numdata' or       */
  /*   'numkern'.  If both 'datainf' and 'kerninf' are of type       */
  /*   PREPPED or FFT then 'numdata' and 'numkern' must have the     */
  /*   same value.  In order for SAME values of 'datainf' and        */
  /*   'kerninf' to help out, the routine must be called with the    */
  /*   same values for 'kern_half_width' and 'numbetween' as well.   */
{
    static fcomplex *kernarray, *dataarray;
    static int firsttime = 1, oldnumbetween, oldkern_half_width;
    static int oldfftlen, oldnumdata, oldlobin;
    fcomplex *tmpdata = NULL, zeros = { 0.0, 0.0 };
    int ii, fftlen = 1, numbins, beginbin, endbin, padlo, padhi;

    /* Check entry arguments for validity */

    if (numdata <= 0) {
        printf("\n  'numdata' = %d (out of bounds) in corr_complex().\n\n", numdata);
        exit(1);
    }
    if (numkern <= 0 || numkern > MAXREALFFT / 4) {
        printf("\n  'numkern' = %d (out of bounds) in corr_complex().\n\n", numkern);
        exit(1);
    }
    if (numresult <= 0 || numresult > MAXREALFFT / 4) {
        printf("\n  'numresult' = %d (out of bounds) in corr_complex().\n\n",
               numresult);
        exit(1);
    }
    if (numbetween < 0 || numbetween > 20000) {
        printf("\n  'numbetween' = %d (out of bounds) in corr_complex().\n\n",
               numbetween);
        exit(1);
    }
    if (numresult % numbetween != 0) {
        printf("\n  'numresult' = %d must be a multiple of ", numresult);
        printf("'numbetween' in corr_complex().\n\n");
        exit(1);
    }
    if (lobin < 0 || lobin >= numdata) {
        printf("\n  'lobin' = %d (out of bounds) in corr_complex().\n\n", lobin);
        exit(1);
    }
    if ((datainf == SAME || kerninf == SAME) && firsttime) {
        printf("\n  Can't call corr_complex() with 'datainf' or 'kerninf'\n");
        printf("  being SAME if this is the 1st time calling the routine.\n\n");
    }
    if (datainf == SAME && kerninf == SAME && lobin == oldlobin) {
        printf("\n  Doesn't make sense to call corr_complex() with SAME for\n");
        printf("  both 'datainf' and 'kerninf' if 'lobin' hasn't changed.\n\n");
        exit(1);
    }
    if (datainf == PREPPED || datainf == FFT || kerninf == PREPPED || kerninf == FFT) {
        if (datainf == PREPPED || datainf == FFT)
            fftlen = numdata;
        if (kerninf == PREPPED || kerninf == FFT)
            fftlen = numkern;
        if ((datainf == PREPPED || datainf == FFT) &&
            (kerninf == PREPPED || kerninf == FFT)) {
            if (numdata != numkern) {
                printf("\n  'numdata' must equal 'numkern' if the data array and\n");
                printf
                    ("  the kernel are either PREPPED or FFT in corr_complex().\n\n");
                exit(1);
            }
        }
    } else if (datainf == SAME || kerninf == SAME) {
        fftlen = oldfftlen;
    } else {
#ifdef DEBUGPRINT
        printf("Yikes!!\n");
#endif
        fftlen = next2_to_n(numresult + 2 * kern_half_width * numbetween);
    }
    if (fftlen > MAXREALFFT / 4) {
        printf("\n  'fftlen' = %d (out of bounds) in corr_complex().\n\n", fftlen);
        exit(1);
    }
    if (kerninf == RAW || kerninf == PREPPED) {
        if (2 * kern_half_width * numbetween > numkern) {
            printf("\n  'numkern = %d (out of bounds in corr_complex().\n", numkern);
            printf("  If 'kerninf' == RAW or PREPPED, 'numkern' must be >=\n");
            printf("  to 2 * 'kern_half_width' * 'numbetween'.\n\n");
            exit(1);
        }
    } else if (kerninf == SAME) {
        if (lobin != oldlobin || kern_half_width != oldkern_half_width) {
            printf("\n  When 'kerninf' = SAME, 'lobin' and 'kern_half_width'\n");
            printf("  must match their previous values in corr_complex().\n\n");
            exit(1);
        }
    }
#ifdef DEBUGPRINT
    printf("numdata = %d  fftlen = %d\n", numdata, fftlen);
#endif

    /* Prep the data array */

    if (firsttime || !(datainf == SAME &&
                       lobin == oldlobin &&
                       fftlen == oldfftlen &&
                       numdata == oldnumdata && numbetween == oldnumbetween)) {

        numbins = fftlen / numbetween;
        beginbin = lobin - kern_half_width;
        endbin = beginbin + numbins;
        tmpdata = data + beginbin;

        if (datainf == RAW) {

#ifdef DEBUGPRINT
            printf("Prepping, spreading, and FFTing the data...\n");
#endif
            /* Zero-pad if necessary (when looking at the beginning or */
            /* the end of the data array)                              */

            if ((beginbin < 0) || (endbin > numdata)) {
                tmpdata = gen_cvect(numbins);
                for (ii = 0; ii < numbins; ii++)
                    tmpdata[ii] = zeros;
                padlo = (beginbin < 0) ? abs(beginbin) : 0;
                padhi = ((endbin - numdata) > 0) ? endbin - numdata : 0;
                memcpy(tmpdata + padlo, data + beginbin + padlo,
                       sizeof(fcomplex) * (numbins - padhi - padlo));
            }

            /* Spread the data */

            if (!firsttime)
                vect_free(dataarray);
            dataarray = gen_cvect(fftlen);
            spread_no_pad(tmpdata, numbins, dataarray, fftlen, numbetween);

            if ((beginbin < 0) || (endbin > numdata))
                vect_free(tmpdata);

            /* FFT the Data */

            COMPLEXFFT(dataarray, fftlen, -1);

        } else if (datainf == PREPPED) {

#ifdef DEBUGPRINT
            printf("FFTing and copying the data...\n");
#endif

            if (!firsttime)
                vect_free(dataarray);
            dataarray = gen_cvect(fftlen);
            memcpy(dataarray, data, sizeof(fcomplex) * fftlen);

            /* FFT the Data */

            COMPLEXFFT(dataarray, fftlen, -1);

        } else if (datainf == FFT) {

#ifdef DEBUGPRINT
            printf("Just copying the data...\n");
#endif

            if (!firsttime)
                vect_free(dataarray);
            dataarray = gen_cvect(fftlen);
            memcpy(dataarray, data, sizeof(fcomplex) * fftlen);

        }
    }

    /* Prep the kernel array */

    if (firsttime || !(kerninf == SAME && fftlen == oldfftlen)) {

        if (!firsttime)
            vect_free(kernarray);
        kernarray = gen_cvect(fftlen);

        if (kerninf == RAW) {

#ifdef DEBUGPRINT
            printf("Placing and FFTing the kernel...\n");
#endif

            place_complex_kernel(kern, numkern, kernarray, fftlen);
            COMPLEXFFT(kernarray, fftlen, -1);

        } else if (kerninf == PREPPED) {

#ifdef DEBUGPRINT
            printf("FFTing and copying the kernel...\n");
#endif

            memcpy(kernarray, kern, sizeof(fcomplex) * fftlen);
            COMPLEXFFT(kernarray, fftlen, -1);

        } else if (kerninf == FFT) {

#ifdef DEBUGPRINT
            printf("Just copying the kernel...\n");
#endif

            memcpy(kernarray, kern, sizeof(fcomplex) * fftlen);
        }
    }

    /* Perform the correlations */

    tmpdata = complex_corr_conv(dataarray, kernarray, fftlen, NOFFTS, optype);

    /* Chop off the contaminated ends and/or the extra data */

    chop_complex_ends(tmpdata, fftlen, result, numresult,
                      kern_half_width * numbetween);
    vect_free(tmpdata);

    /* Set variables for next time... */

    if (firsttime)
        firsttime = 0;
    oldlobin = lobin;
    oldfftlen = fftlen;
    oldnumdata = numdata;
    oldnumbetween = numbetween;
    oldkern_half_width = kern_half_width;
    if (numresult < fftlen - 2 * numbetween * kern_half_width)
        return numresult;
    else
        return fftlen - 2 * numbetween * kern_half_width;
}


void stretch_fft(fcomplex * data, int numdata, fcomplex * result, int numresult)
  /* This routine stretches and/or interpolates an FFT of length    */
  /* numdata.  It zeros 'result' where end-effects have occurred.   */
  /* This routine is usually used to co-add stretched FFTs to       */
  /* increase the signal-to-noise ratios of a detection.            */
  /* Arguments:                                                     */
  /*   'data' is a pointer to a complex FFT.                        */
  /*   'numdata' is the number of complex points in 'data'.         */
  /*   'result' is a pointer to the complex stretched array.        */
  /*   'numresult' is the number of complex points in 'result'.     */
  /* Notes:                                                         */
  /*   The ratio of 'numresult' to 'numdata' determines the amount  */
  /*   of stretching that will take place.  For example, if         */
  /*   'numresult' is twice 'numdata', then the data will be        */
  /*   stretched by a factor of two (i.e. interbinned).             */
{

    int numkern, numbetween, kern_half_width;
    fcomplex *kernel;

    numbetween = numresult / numdata;
    kern_half_width = r_resp_halfwidth(LOWACC);
    numkern = 2 * numbetween * kern_half_width;
    kernel = gen_r_response(0.0, numbetween, numkern);

    /* Perform the correlation */

    corr_complex(data, numdata, RAW,
                 kernel, numkern, RAW,
                 result, numresult, 0, numbetween, kern_half_width, CORR);
    vect_free(kernel);
}


float *corr_loc_pow(float *powers, int numpowers)
  /* This routine determines the local power levels for every         */
  /* frequency in an FFT containing 'numpowers' complex frequencies.  */
  /* It sets the areas where end-effects are a problem to the         */
  /* local power level of the closest bin without end effect          */
  /* problems.  It returns a vector with the local power levels.      */
  /* Arguments:                                                       */
  /*   'powers' is a pointer to a float vector of powers.             */
  /*   'numpowers' is the number of complex points in 'powers'.       */
{
    int ii, kern_half_width;
    float *local_powers, normal;
    static float *kernel;
    static int firsttime = 1;
    static int old_numpowers;

    normal = (numpowers / 2) / NUMLOCPOWAVG;
    kern_half_width = NUMLOCPOWAVG / 2 + DELTAAVGBINS;

    /* Make sure our powers array is long enough. */

    if (numpowers < 2 * kern_half_width) {
        printf("\n  numpowers = %d (out of bounds) in corr_loc_pow().\n\n",
               numpowers);
        exit(1);
    }

    /* Set up the kernel if we need to */

    if (firsttime || numpowers != old_numpowers) {
        if (!firsttime)
            vect_free(kernel);
        kernel = gen_fvect(numpowers);
        for (ii = 0; ii < numpowers; ii++)
            kernel[ii] = 0.0;
        for (ii = DELTAAVGBINS + 1; ii < (kern_half_width + 1); ii++) {
            kernel[ii] = normal;
            kernel[numpowers - ii] = normal;
        }
        realfft(kernel, numpowers, -1);
        old_numpowers = numpowers;
        firsttime = 0;
    }

    /* Do the correlation */

    local_powers = real_corr_conv(powers, kernel, numpowers, FFTD, CORR);

    /* Correct the nasty end-effects by setting these values equal to */
    /* the last good values.                                          */

    for (ii = 0; ii < kern_half_width; ii++) {
        local_powers[ii] = local_powers[kern_half_width];
        local_powers[numpowers - 1 - ii] = local_powers[numpowers - 1 -
                                                        kern_half_width];
    }
    return local_powers;
}
