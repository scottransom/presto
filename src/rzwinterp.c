#include "presto.h"

/*
fcomplex *corr_rzw_interp(fcomplex *data, int numdata, int numbetween, \
			  int startbin, double z, double w, int fftlen, \
			  presto_interp_acc accuracy, int *nextbin)
*/
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation of part of the f-fdot-fdotdot volume.     */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'z' is the fdot to use (z=f-dot*T^2).                         */
  /*   'w' is the fdotdot to use (z=f-dotdot*T^3).                   */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

  /* Obviously not implemented yet! */


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
