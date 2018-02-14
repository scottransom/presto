#include "presto.h"

int next_good_fftlen(int N)
/* Return one of the shortest, yet best performing, FFT lengths larger
 * than N.  This assumes FFTW. */
{
    int fftlens[17] = {128, 192, 256, 384, 512, 768, 1024, 1280, 2048, 4096,
                       5120, 7680, 10240, 12288, 15360, 16384, 25600};
    int ii = 0;
    if (N <= fftlens[0])
        return fftlens[0];
    if (N > fftlens[16])
        return next2_to_n(N);
    while (N > fftlens[ii]) ii++;
    return fftlens[ii];
}

int fftlen_from_kernwidth(int kernwidth)
/* return the length of the optimal FFT to use for correlations with
 * some kernel width kernwidth.  This assumes FFTW. */
{
    // The following nummbers were determined using FFTW 3.3.7 on an
    // AVX-enabled processor.  Metric used was max throughput of good
    // correlated data.
    if (kernwidth < 6) return 128;
    else if (kernwidth < 52) return 256;
    else if (kernwidth < 67) return 512;
    else if (kernwidth < 378) return 1024;
    else if (kernwidth < 664) return 2048;
    else if (kernwidth < 1672) return 4096;
    else if (kernwidth < 3015) return 10240;
    else if (kernwidth < 3554) return 15360;
    else if (kernwidth < 6000) return 25600;
    else return next2_to_n(kernwidth*5);
}

void spread_with_pad(fcomplex * data, int numdata,
                     fcomplex * result, int numresult, int numbetween, int numpad)
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array and padding it.                   */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */
  /*   'numpad' is the number of bins to use as zero padding     */
{
    int ii, jj, numtoplace;
    fcomplex zeros = { 0.0, 0.0 };

    for (ii = 0; ii < numresult; ii++)
        result[ii] = zeros;
    numtoplace = (numresult - numpad) / numbetween;
    if (numtoplace > numdata)
        numtoplace = numdata;
    for (ii = 0, jj = 0; ii < numtoplace; ii++, jj += numbetween)
        result[jj] = data[ii];
}


void spread_no_pad(fcomplex * data, int numdata,
                   fcomplex * result, int numresult, int numbetween)
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array.                                  */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */
{
    spread_with_pad(data, numdata, result, numresult, numbetween, 0);
}


void paddata(fcomplex * data, int numdata, int numpad)
  /* Pad the last 'numpad' bins of 'data' with zeros.         */
  /* Arguments:                                               */
  /*   'data' is the FFT array to be padded                   */
  /*   'numdata' is the number of complex points in 'data'    */
  /*   'numpad' is the number of bins to use as zero padding  */
{
    int ii;
    fcomplex zeros = { 0.0, 0.0 };

    for (ii = numdata - numpad; ii < numdata; ii++)
        data[ii] = zeros;
}


void place_complex_kernel(fcomplex * kernel, int numkernel,
                          fcomplex * result, int numresult)
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a complex response function.  Bin zero  */
  /*      response is in bin numkernel/2.                  */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */
{
    int ii, halfwidth;
    fcomplex zeros = { 0.0, 0.0 };

    halfwidth = numkernel / 2;
    for (ii = 0; ii < numresult; ii++)
        result[ii] = zeros;
    memcpy(result, kernel + halfwidth, sizeof(fcomplex) * halfwidth);
    memcpy(result + numresult - halfwidth, kernel, sizeof(fcomplex) * halfwidth);
}


void place_real_kernel(float *kernel, int numkernel, float *result, int numresult)
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a real-valued response function.  Bin   */
  /*      zero response is in bin numkernel/2.             */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */
{
    int ii, halfwidth;

    halfwidth = numkernel / 2;
    for (ii = 0; ii < numresult; ii++)
        result[ii] = 0.0;
    memcpy(result, kernel + halfwidth, sizeof(float) * halfwidth);
    memcpy(result + numresult - halfwidth, kernel, sizeof(float) * halfwidth);
}


void chop_complex_ends(fcomplex * data, int numdata,
                       fcomplex * result, int numresult, int chopbins)
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */
{
    int ii, numtocopy;
    fcomplex zeros = { 0.0, 0.0 };

    if (numdata < 2 * chopbins) {
        printf("\n  'chopbins' is too large in chop_complex_ends()\n\n");
        exit(1);
    }
    numtocopy = numdata - 2 * chopbins;
    if (numresult < numtocopy)
        numtocopy = numresult;
    for (ii = 0; ii < numresult; ii++)
        result[ii] = zeros;
    memcpy(result, data + chopbins, sizeof(fcomplex) * numtocopy);
}


void chop_real_ends(float *data, int numdata,
                    float *result, int numresult, int chopbins)
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */
{
    if (numdata < 2 * chopbins) {
        printf("\n  'chopbins' is too large in chop_complex_ends()\n\n");
        exit(1);
    }
    if (numresult < numdata - 2 * chopbins) {
        printf("\n  'numresult' is too small in chop_complex_ends()\n\n");
        exit(1);
    }
    memcpy(result, data + chopbins, sizeof(float) * (numdata - 2 * chopbins));
}
