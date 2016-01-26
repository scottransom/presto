#include "presto.h"

fcomplex *complex_corr_conv(fcomplex * data, fcomplex * kernel,
                            int numdata, presto_ffts ffts, presto_optype type)
  /* Perform and return a complex correlation or convolution.       */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */
{
    int ii;
    float normal, tmpd, tmpk;
    fcomplex *tmpdat;

    /* Get the normalization factor */

    normal = 1.0 / numdata;

    /* Set up all our parameters */

    if (ffts > 3) {
        printf("\nIllegal 'ffts' option (%d) in complex_corr_conv().\n", ffts);
        printf("Exiting.\n\n");
        exit(1);
    }
    if (type > 3) {
        printf("\nIllegal 'type' option (%d) in complex_corr_conv().\n", type);
        printf("Exiting.\n\n");
        exit(1);
    }
    if (type == INPLACE_CONV || type == INPLACE_CORR) {
        tmpdat = data;
    } else {
        tmpdat = gen_cvect(numdata);
        memcpy(tmpdat, data, sizeof(fcomplex) * numdata);
    }
    if (ffts == FFTDK || ffts == FFTD) {
        COMPLEXFFT(tmpdat, numdata, -1);
    }
    if (ffts == FFTDK || ffts == FFTK) {
        COMPLEXFFT(kernel, numdata, -1);
    }

    /* Do the complex multiplications */

    if (type == CORR || type == INPLACE_CORR) {
        for (ii = 0; ii < numdata; ii++) {
            tmpd = tmpdat[ii].r;
            tmpk = kernel[ii].r;
            tmpdat[ii].r = (tmpd * tmpk + tmpdat[ii].i * kernel[ii].i)
                * normal;
            tmpdat[ii].i = (tmpdat[ii].i * tmpk - kernel[ii].i * tmpd)
                * normal;
        }
    } else {
        for (ii = 0; ii < numdata; ii++) {
            tmpd = tmpdat[ii].r;
            tmpk = kernel[ii].r;
            tmpdat[ii].r = (tmpd * tmpk - tmpdat[ii].i * kernel[ii].i)
                * normal;
            tmpdat[ii].i = (tmpdat[ii].i * tmpk + kernel[ii].i * tmpd)
                * normal;
        }
    }

    /* Perform the inverse FFT on the result and return */

    COMPLEXFFT(tmpdat, numdata, 1);
    return tmpdat;
}


float *real_corr_conv(float *data, float *kernel, int numdata,
                      presto_ffts ffts, presto_optype type)
  /* Perform and return a real-valued correlation or convolution.   */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */
{
    int ii;
    float normal, tmpd, tmpk, *tmpdat;
    fcomplex *fcdata, *fckern;

    /* Get the normalization factor */

    normal = 2.0 / numdata;

    /* Set up all our parameters */

    if (ffts > 3) {
        printf("\nIllegal 'ffts' option (%d) in real_corr_conv().\n", ffts);
        printf("Exiting\n\n");
        exit(1);
    }
    if (type > 3) {
        printf("\nIllegal 'type' option (%d) in real_corr_conv().\n", ffts);
        printf("Exiting\n\n");
        exit(1);
    }
    if (type == INPLACE_CONV || type == INPLACE_CORR) {
        tmpdat = data;
    } else {
        tmpdat = gen_fvect(numdata);
        memcpy(tmpdat, data, sizeof(float) * numdata);
    }
    if (ffts == FFTDK || ffts == FFTD) {
        realfft(tmpdat, numdata, -1);
    }
    if (ffts == FFTDK || ffts == FFTK) {
        realfft(kernel, numdata, -1);
    }
    // Act like our packed-complex floats are really fcomplex values

    fcdata = (fcomplex *) tmpdat;
    fckern = (fcomplex *) kernel;

    // Do the complex multiplications

    if (type == CORR || type == INPLACE_CORR) {
        for (ii = 1; ii < numdata / 2; ii++) {
            tmpd = fcdata[ii].r;
            tmpk = fckern[ii].r;
            fcdata[ii].r = (tmpd * tmpk + fcdata[ii].i * fckern[ii].i)
                * normal;
            fcdata[ii].i = (fcdata[ii].i * tmpk - fckern[ii].i * tmpd)
                * normal;
        }
    } else {
        for (ii = 1; ii < numdata / 2; ii++) {
            tmpd = fcdata[ii].r;
            tmpk = fckern[ii].r;
            fcdata[ii].r = (tmpd * tmpk - fcdata[ii].i * fckern[ii].i)
                * normal;
            fcdata[ii].i = (fcdata[ii].i * tmpk + fckern[ii].i * tmpd)
                * normal;
        }
    }

    // Now handle bin zero

    fcdata[0].r *= fckern[0].r * normal;
    fcdata[0].i *= fckern[0].i * normal;

    // Perform the inverse FFT on the result and return

    realfft(tmpdat, numdata, 1);
    return tmpdat;
}
