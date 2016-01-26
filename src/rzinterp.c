#include "presto.h"

fcomplex **corr_rz_plane(fcomplex * data, int numdata, int numbetween,
                         int startbin, double zlo, double zhi,
                         int numz, int fftlen,
                         presto_interp_acc accuracy, int *nextbin)
  /* This routine uses the correlation method to do Fourier          */
  /* complex interpolations of the f-fdot plane.                     */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'zlo' is the lowest fdot to use (z=f-dot*T^2)                 */
  /*   'zhi' is the highest fdot to use (z=f-dot*T^2)                */
  /*   'numz' is the number of z values to use to make the plane     */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */
{
    double maxz, dz, dtmp;
    int ii, kern_half_width, numkern, numgoodbins, numgoodbins2;
    static int firsttime = 1, old_numbetween;
    static int old_kern_half_width, old_fftlen, old_numz;
    static double old_zlo, old_zhi;
    static fcomplex **kernel = NULL, **result, *response;
    static presto_interp_acc old_accuracy;
    presto_datainf datainf;

    /* Make sure our z's are in the correct order */

    if (zlo > zhi) {
        dtmp = zlo;
        zlo = zhi;
        zhi = dtmp;
    }

    /* Determine some required parameters */

    if (numz == 1)
        dz = 0.0;
    else
        dz = (zhi - zlo) / (numz - 1);
    maxz = (fabs(zlo) < fabs(zhi)) ? zhi : zlo;
    kern_half_width = z_resp_halfwidth(maxz, accuracy);
    numkern = 2 * numbetween * kern_half_width;

    /* If this is our first time or the parameters have changed, */
    /* set up and FFT the correlation kernel array               */

    if (firsttime
        || numbetween != old_numbetween
        || kern_half_width != old_kern_half_width
        || fftlen != old_fftlen
        || fabs(zlo - old_zlo) > 1.0E-5
        || fabs(zhi - old_zhi) > 1.0E-5
        || numz != old_numz || accuracy != old_accuracy) {

        /* Cleanup old data if we need to */

        if (kernel) {
            vect_free(kernel[0]);
            vect_free(kernel);
        }
        kernel = gen_cmatrix(numz, fftlen);

        /* Generate the new responses and the kernel matrix */

        for (ii = 0; ii < numz; ii++) {
            response = gen_z_response(0.0, numbetween, zlo + ii * dz, numkern);
            place_complex_kernel(response, numkern, kernel[ii], fftlen);
            vect_free(response);
            COMPLEXFFT(kernel[ii], fftlen, -1);
        }

        /* Set the old_ variables equal to the new variables  */

        old_numbetween = numbetween;
        old_kern_half_width = kern_half_width;
        old_zlo = zlo;
        old_zhi = zhi;
        old_numz = numz;
        old_fftlen = fftlen;
        old_accuracy = accuracy;
        firsttime = 0;
    }

    /* Allocate the result matrix */

    numgoodbins = fftlen - 2 * kern_half_width * numbetween;
    if (numgoodbins % numbetween) {
        while (numgoodbins % numbetween)
            numgoodbins--;
    }
    *nextbin = startbin + numgoodbins / numbetween;
    result = gen_cmatrix(numz, numgoodbins);

    /* Perform the correlations. */

    datainf = RAW;
    for (ii = 0; ii < numz; ii++) {
        numgoodbins2 = corr_complex(data, numdata, datainf,
                                    kernel[ii], fftlen, FFT,
                                    result[ii], numgoodbins, startbin,
                                    numbetween, kern_half_width, CORR);
        if (numgoodbins != numgoodbins2) {
            printf("\n  Calculated 'numgoodbins' = %d, does not equal\n",
                   numgoodbins);
            printf("\n  returned 'numgoodbins' = %d.  Exiting.\n\n", numgoodbins2);
        }
        datainf = SAME;
    }
    return result;
}


fcomplex *corr_rz_interp(fcomplex * data, int numdata, int numbetween,
                         int startbin, double z, int fftlen,
                         presto_interp_acc accuracy, int *nextbin)
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation of a slice of the f-fdot plane.           */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'z' is the fdot to use (z=f-dot*T^2).                         */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */
{
    fcomplex **result, *tempreturn;

    result = corr_rz_plane(data, numdata, numbetween, startbin, z, z, 1,
                           fftlen, accuracy, nextbin);
    tempreturn = result[0];
    vect_free(result);
    return tempreturn;
}


void rz_interp(fcomplex * data, int numdata, double r, double z,
               int kern_half_width, fcomplex * ans)
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation at a single point in the f-fdot plane.    */
  /* It does the correlations manually. (i.e. no FFTs)               */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'r' is the Fourier frequency in data that we want to          */
  /*      interpolate.  This can (and should) be fractional.         */
  /*   'z' is the fdot to use (z=f-dot*T^2 (T is integration time)). */
  /*   'kern_half_width' is the half-width of the kernel in bins.    */
  /*   'ans' is the complex answer.                                  */
{
    float *dataptr, *respptr;
    int ii, numkern, nsum, intfreq, lodata, hidata, loresp, hiresp;
    double fracfreq, dintfreq, tmpd, tmpr;
    fcomplex *response;

    /* Check 'r' and return 0.0 + 0.0i if out of bounds.        */
    /* Should this return an error and exit instead?            */

    if (r > numdata - 1.0 || r < 0.0) {
        ans->r = 0.0;
        ans->i = 0.0;
        return;
    }

    /* Split 'r' into integer and fractional parts */

    fracfreq = modf(r, &dintfreq);
    intfreq = (int) dintfreq;

    /* Return immediately if 'r' is close to an           */
    /* integer frequency and z is very close to zero      */

    if (fabs(z) < 1E-4) {
        if (fracfreq < 1E-5) {
            ans->r = data[intfreq].r;
            ans->i = data[intfreq].i;
            return;
        }
        if ((1.0 - fracfreq) < 1E-5) {
            ans->r = data[intfreq + 1].r;
            ans->i = data[intfreq + 1].i;
            return;
        }
    }

    /* Generate the response function */

    numkern = 2 * kern_half_width;
    response = gen_z_response(fracfreq, 1, z, numkern);

    /* Determine the summation boundaries */

    lodata = intfreq - kern_half_width;
    if (lodata < 0) {
        loresp = abs(lodata);
        lodata = 0;
    } else {
        loresp = 0;
    }
    hidata = intfreq + kern_half_width - 1;
    if (hidata > numdata - 1) {
        hiresp = numkern - hidata + numdata - 1;
    } else {
        hiresp = numkern;
    }
    nsum = hiresp - loresp;

    /* Set up our pointers */

    dataptr = (float *) (data + lodata);
    respptr = (float *) (response + loresp);

    /* Do the summation */

    ans->r = 0.0;
    ans->i = 0.0;

    for (ii = 0; ii < nsum; ii++) {
        tmpd = *(dataptr++);
        tmpr = *(respptr++);
        ans->r += tmpd * tmpr + (*dataptr) * (*respptr);
        ans->i += (*dataptr) * tmpr - (*respptr) * tmpd;
        dataptr++;
        respptr++;
    }

    vect_free(response);
    return;
}
