#include "presto.h"

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
