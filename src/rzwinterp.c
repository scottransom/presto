#include "presto.h"

fcomplex ***corr_rzw_vol(fcomplex * data, int numdata, int numbetween,
                         int startbin, double zlo, double zhi, int numz,
                         double wlo, double whi, int numw, int fftlen,
                         presto_interp_acc accuracy, int *nextbin)
  /* This routine uses the correlation method to do Fourier          */
  /* complex interpolations of the f-fdot-fdotdot volume.            */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'zlo' is the lowest fdot to use (z=f-dot*T^2)                 */
  /*   'zhi' is the highest fdot to use (z=f-dot*T^2)                */
  /*   'numz' is the number of z values to use to make the volume    */
  /*   'wlo' is the lowest fdotdot to use (w=f-dotdot*T^3)           */
  /*   'whi' is the highest fdotdot to use (w=f-dotdot*T^3)          */
  /*   'numw' is the number of w values to use to make the volume    */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */
{
    double maxz, dz, maxw, dw, dtmp;
    int ii, jj, kern_half_width, numkern, numgoodbins, numgoodbins2;
    static int firsttime = 1, old_numbetween;
    static int old_kern_half_width, old_fftlen, old_numz, old_numw;
    static double old_zlo, old_zhi, old_wlo, old_whi;
    static fcomplex ***kernel = NULL, ***result, *response;
    static presto_interp_acc old_accuracy;
    presto_datainf datainf;

    /* Make sure our z's and w's are in the correct order */

    if (zlo > zhi) {
        dtmp = zlo;
        zlo = zhi;
        zhi = dtmp;
    }
    if (wlo > whi) {
        dtmp = wlo;
        wlo = whi;
        whi = dtmp;
    }

    /* Determine some required parameters */

    if (numz == 1)
        dz = 0.0;
    else
        dz = (zhi - zlo) / (numz - 1);
    maxz = (fabs(zlo) < fabs(zhi)) ? zhi : zlo;
    if (numw == 1)
        dw = 0.0;
    else
        dw = (whi - wlo) / (numw - 1);
    maxw = (fabs(wlo) < fabs(whi)) ? whi : wlo;
    kern_half_width = w_resp_halfwidth(maxz, maxw, accuracy);
    numkern = 2 * numbetween * kern_half_width;

    /* If this is our first time or the parameters have changed, */
    /* set up and FFT the correlation kernel array               */

    if (firsttime
        || numbetween != old_numbetween
        || kern_half_width != old_kern_half_width
        || fftlen != old_fftlen
        || fabs(zlo - old_zlo) > 1.0E-5
        || fabs(zhi - old_zhi) > 1.0E-5
        || numz != old_numz
        || fabs(wlo - old_wlo) > 1.0E-5
        || fabs(whi - old_whi) > 1.0E-5
        || numw != old_numw
        || accuracy != old_accuracy) {

        /* Cleanup old data if we need to */

        if (kernel) {
            vect_free(kernel[0][0]);
            vect_free(kernel[0]);
            vect_free(kernel);
        }
        kernel = gen_c3Darr(numw, numz, fftlen);

        /* Generate the new responses and the kernel matrix */

        for (ii = 0; ii < numw; ii++) {
            for (jj = 0; jj < numz; jj++) {
                response = gen_w_response(0.0, numbetween,
                                          zlo + jj * dz,
                                          wlo + ii * dw, numkern);
                place_complex_kernel(response, numkern, kernel[ii][jj], fftlen);
                vect_free(response);
                COMPLEXFFT(kernel[ii][jj], fftlen, -1);
            }
        }

        /* Set the old_ variables equal to the new variables  */

        old_numbetween = numbetween;
        old_kern_half_width = kern_half_width;
        old_zlo = zlo;
        old_zhi = zhi;
        old_numz = numz;
        old_wlo = wlo;
        old_whi = whi;
        old_numw = numw;
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
    result = gen_c3Darr(numw, numz, numgoodbins);

    /* Perform the correlations. */

    datainf = RAW;
    for (ii = 0; ii < numw; ii++) {
        for (jj = 0; jj < numz; jj++) {
            numgoodbins2 = corr_complex(data, numdata, datainf,
                                        kernel[ii][jj], fftlen, FFT,
                                        result[ii][jj], numgoodbins, startbin,
                                        numbetween, kern_half_width, CORR);
            if (numgoodbins != numgoodbins2) {
                printf("\n  Calculated 'numgoodbins' = %d, does not equal\n",
                       numgoodbins);
                printf("\n  returned 'numgoodbins' = %d.  Exiting.\n\n", numgoodbins2);
            }
            datainf = SAME;
        }
    }
    return result;
}


void rzw_interp(fcomplex * data, int numdata, double r, double z,
                double w, int kern_half_width, fcomplex * ans)
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation at a single point in the f-fdot plane.    */
  /* It does the correlations manually. (i.e. no FFTs)               */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'r' is the Fourier frequency in data that we want to          */
  /*      interpolate.  This can (and should) be fractional.         */
  /*   'z' is the fdot to use (z=f-dot*T^2 (T is integration time)). */
  /*   'w' is the fdotdot to use (z=f-dotdot*T^3).                   */
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

    /* Return immediately if 'w' is close to zero  */

    if (fabs(w) < 1E-4) {
        rz_interp(data, numdata, r, z, kern_half_width, ans);
        return;
    }

    /* Generate the response function */

    numkern = 2 * kern_half_width;
    response = gen_w_response(fracfreq, 1, z, w, numkern);

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
