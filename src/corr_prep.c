#include "presto.h"

int next_good_fftlen(int N)
/* Return the FFT length >= N that is cheapest for accelsearch's      */
/* correlations (complex multiply + inverse FFT).  This assumes FFTW. */
{
    // Up to 8192, these are the lengths that are cheaper than every
    // longer one on both a Sapphire Rapids and a Skylake-SP Xeon (their
    // system FFTWs, 3.3.8 and 3.3.5, single thread), timed with
    // accel_fftlen_bench (tests/accelsearch/fftlen_bench.c), plus 5376
    // (cheapest on Skylake-SP, neutral on Sapphire Rapids).  On each
    // machine, this list is within ~1% of that machine's own best lengths
    // on average, and within 0.6% with conda-forge's FFTW 3.3.10.
    // Powers of 2 are much faster per point than other lengths, so the
    // others only pay off just above a power of 2.  The lengths above
    // 8192 are from the older, FFTW 3.3.7-era list.
    static const int fftlens[] = {
        128, 160, 192, 240, 256, 320, 384, 400, 512, 640, 768, 1024,
        1280, 1344, 1536, 2048, 2304, 2400, 2560, 2880, 3072, 3200, 4096,
        4800, 5120, 5376, 6144, 6400, 6912, 7680, 8192,
        10240, 12288, 15360, 16384, 25600
    };
    const int numlens = sizeof(fftlens) / sizeof(fftlens[0]);
    int ii = 0;
    if (N <= fftlens[0])
        return fftlens[0];
    if (N > fftlens[numlens - 1])
        return next2_to_n(N);
    while (N > fftlens[ii])
        ii++;
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
