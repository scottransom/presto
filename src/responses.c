#include "presto.h"

#define MIN_NUMDATA 131072
#define MIN_NUMORBPTS 2049      /* This should be a power-of-two + 1 */
#define NUM_PTS_WDAT 131072

/* Function declarations */
int fresnl(double xxa, double *ssa, double *cca);


int r_resp_halfwidth(presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a standard Fourier signal.                  */
  /*  Arguments:                                                       */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
   if (accuracy == HIGHACC) {
      return ((NUMFINTBINS * 3) + (NUMLOCPOWAVG >> 1) + DELTAAVGBINS);
   } else {
      return NUMFINTBINS;
   }
}


int z_resp_halfwidth(double z, presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
  /*  a constant frequency derivative)                                 */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
   int m;

   z = fabs(z);

   if (accuracy == HIGHACC) {
      m = (long) (z * (0.002057 * z + 0.0377) + NUMFINTBINS * 3);
      m += ((NUMLOCPOWAVG >> 1) + DELTAAVGBINS);

      /* Prevent the equation from blowing up in large z cases */

      if (z > 100 && m > 1.2 * z)
         m = 1.2 * z;

   } else {
      m = (long) (z * (0.00089 * z + 0.3131) + NUMFINTBINS);
      m = (m < NUMFINTBINS) ? NUMFINTBINS : m;

      /* Prevent the equation from blowing up in large z cases */

      if (z > 100 && m > 0.6 * z)
         m = 0.6 * z;
   }
   return m;
}


int w_resp_halfwidth(double z, double w, presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with an f-dot that (i.e    */
  /*  varies linearly in time -- a constant f-dotdot)                  */
  /*  Arguments:                                                       */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
   /* This needs to be made more robust... */

   if (fabs(w) < 1.0e-7)
      return z_resp_halfwidth(z, accuracy);
   else
      return (int) fabs(z) + r_resp_halfwidth(accuracy);
}


void binary_velocity(double T, orbitparams * orbit, double *minv, double *maxv)
  /*  Return the minimum and maximum orbital velocities of a pulsar    */
  /*  during an observation as a fraction of the speed of light.       */
  /*  Arguments:                                                       */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
{
   *minv = 1.0;
   *maxv = -1.0;
   if (T >= orbit->p) {

      double c1, c2;

      c1 = TWOPI * orbit->x / (orbit->p * sqrt(1.0 - orbit->e * orbit->e));
      c2 = orbit->e * cos(orbit->w * DEGTORAD);
      *maxv = c1 * (c2 + 1.0);
      *minv = c1 * (c2 - 1.0);

   } else {

      double dtb, startE, *E;
      int ii, numpoints = 1025;
      orbitparams orb;

      dtb = T / (double) (numpoints - 1);
      orb.p = orbit->p;
      orb.x = orbit->x;
      orb.e = orbit->e;
      orb.w = orbit->w * DEGTORAD;
      orb.t = orbit->t;
      startE = keplars_eqn(orb.t, orb.p, orb.e, 1.0E-15);
      E = dorbint(startE, numpoints, dtb, &orb);
      E_to_v(E, numpoints, &orb);
      for (ii = 0; ii < numpoints; ii++) {
         E[ii] *= 1000.0 / SOL;
         if (E[ii] < *minv)
            *minv = E[ii];
         if (E[ii] > *maxv)
            *maxv = E[ii];
      }
      vect_free(E);
   }
   if (*maxv < *minv) {
      printf("Something is wrong in binary_velocity()\n");
   }
}

int bin_resp_halfwidth(double ppsr, double T, orbitparams * orbit)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a pulsar in a binary orbit.                 */
  /*  Arguments:                                                       */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2 * 'numbetween' to get the   */
  /*    length of the array required to hold such a kernel.            */
{
   double maxv, minvel, maxvel, maxdevbins;
   int retval;

   binary_velocity(T, orbit, &minvel, &maxvel);
   maxv = (fabs(minvel) > fabs(maxvel)) ? minvel : maxvel;
   maxdevbins = fabs(T * maxv / (ppsr * (1.0 + maxv)));
   retval = (int) floor(1.1 * maxdevbins + 0.5);
   return (retval < NUMFINTBINS) ? NUMFINTBINS : retval;
}


fcomplex *gen_r_response(double roffset, int numbetween, int numkern)
  /*  Generate a complex response function for Fourier interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{
   int ii;
   double tmp, sinc, s, c, alpha, beta, delta, startr, r;
   fcomplex *response;

   /* Check that the arguments are OK */

   if (roffset < 0.0 || roffset >= 1.0) {
      printf("\n  roffset = %f (out of bounds) in gen_r_response().\n\n", roffset);
      exit(-1);
   }
   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_r_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_r_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_r_response().\n\n");
      exit(-1);
   }

   /* Prep the recursion */

   response = gen_cvect(numkern);
   startr = PI * (numkern / (double) (2 * numbetween) + roffset);
   delta = -PI / numbetween;
   tmp = sin(0.5 * delta);
   alpha = -2.0 * tmp * tmp;
   beta = sin(delta);
   c = cos(startr);
   s = sin(startr);

   /* Generate the points */

   for (ii = 0, r = startr; ii < numkern; ii++, r += delta) {
      if (r == 0.0)
         sinc = 1.0;
      else
         sinc = s / r;
      response[ii].r = c * sinc;
      response[ii].i = s * sinc;
      c = alpha * (tmp = c) - beta * s + c;
      s = alpha * s + beta * tmp + s;
   }

   /* Correct for divide by zero when the roffset is close to zero */

   if (roffset < 1E-3) {
      response[numkern / 2].r = 1 - 6.579736267392905746 * (tmp = roffset * roffset);
      response[numkern / 2].i = roffset * (PI - 10.335425560099940058 * tmp);
   }
   return response;
}


fcomplex *gen_z_response(double roffset, int numbetween, double z, int numkern)
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{
   int ii, signz, numkernby2;
   double absz, zd, tmp, r, xx, yy, zz, startr, startroffset;
   double fressy, frescy, fressz, frescz, tmprl, tmpim;
   double s, c, pibyz, cons, delta;
   fcomplex *response;

   /* Check that the arguments are OK */

   if (roffset < 0.0 || roffset >= 1.0) {
      printf("\n  roffset = %f (out of bounds) in gen_z_response().\n\n", roffset);
      exit(-1);
   }
   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_z_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_z_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_z_response().\n\n");
      exit(-1);
   }

   /* If z~=0 use the normal Fourier interpolation kernel */

   absz = fabs(z);
   if (absz < 1E-4) {
      response = gen_r_response(roffset, numbetween, numkern);
      return response;
   }

   response = gen_cvect(numkern);

   /* Begin the calculations */

   startr = roffset - (0.5 * z);
   startroffset = (startr < 0) ? 1.0 + modf(startr, &tmprl) : modf(startr, &tmprl);
   signz = (z < 0.0) ? -1 : 1;
   zd = signz * SQRT2 / sqrt(absz);
   cons = zd / 2.0;
   pibyz = PI / z;
   startr += numkern / (double) (2 * numbetween);
   delta = -1.0 / numbetween;

   for (ii = 0, r = startr; ii < numkern; ii++, r += delta) {
      yy = r * zd;
      zz = yy + z * zd;
      xx = pibyz * r * r;
      c = cos(xx);
      s = sin(xx);
      fresnl(yy, &fressy, &frescy);
      fresnl(zz, &fressz, &frescz);
      tmprl = signz * (frescz - frescy);
      tmpim = fressy - fressz;
      response[ii].r = ((tmp = tmprl) * c - tmpim * s) * cons;
      response[ii].i = -(tmp * s + tmpim * c) * cons;
   }

   /* Correct for divide by zero when the roffset and z is close to zero */

   if (startroffset < 1E-3 && absz < 1E-3) {
      zz = z * z;
      xx = startroffset * startroffset;
      numkernby2 = numkern / 2;
      response[numkernby2].r = 1.0 - 0.16449340668482264365 * zz;
      response[numkernby2].i = -0.5235987755982988731 * z;
      response[numkernby2].r += startroffset * 1.6449340668482264365 * z;
      response[numkernby2].i += startroffset * (PI - 0.5167712780049970029 * zz);
      response[numkernby2].r += xx * (-6.579736267392905746
                                      + 0.9277056288952613070 * zz);
      response[numkernby2].i += xx * (3.1006276680299820175 * z);
   }
   return response;
}


fcomplex *gen_w_response(double roffset, int numbetween, double z,
                         double w, int numkern)
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{

   int fftlen, ii, beginbin, numintkern;
   float *data;
   double amp, f, fd, fdd, dt, t, phase;
   static int old_numbetween = 0, old_numkern = 0, old_fftlen = 0, firsttime = 1;
   static fcomplex *kernelarray = NULL;
   fcomplex *response, *tmpresponse, *rresp, *dataarray;

   /* Check that the arguments are OK */

   if (roffset < 0.0 || roffset >= 1.0) {
      printf("\n  roffset = %f (out of bounds) in gen_w_response().\n\n", roffset);
      exit(-1);
   }
   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_w_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_w_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_w_response().\n\n");
      exit(-1);
   }

   /* If w~=0 use the normal F-dot Fourier interpolation kernel */

   if (fabs(w) < 1E-4) {
      response = gen_z_response(roffset, numbetween, z, numkern);
      return response;
   }

   /* Otherwise initialize some data */

   dt = 1.0 / (double) NUM_PTS_WDAT;
   amp = 2.0 * dt;
   f = (double) (NUM_PTS_WDAT / 4);
   fd = (z - 0.5 * w) / 2.0;
   fdd = w / 6.0;

   /* Generate the data set */

   data = gen_fvect(NUM_PTS_WDAT);
   for (ii = 0; ii < NUM_PTS_WDAT; ii++) {
      t = ii * dt;
      phase = TWOPI * (t * (t * (t * fdd + fd) + f));
      data[ii] = amp * cos(phase);
   }

   /* FFT the data */

   realfft(data, NUM_PTS_WDAT, -1);

   /* The following block saves us from having to re-compute */
   /* the Fourier interpolation kernels if 'numkern' is the  */
   /* same length as on prior calls.                         */

   fftlen = next2_to_n(numkern);
   beginbin = NUM_PTS_WDAT / 4 - numkern / numbetween;
   if (firsttime ||
       old_numkern != numkern ||
       old_numbetween != numbetween || old_fftlen != fftlen) {

      /* Generate an interpolation kernel for the data */

      numintkern = 2 * numbetween * r_resp_halfwidth(HIGHACC);
      rresp = gen_r_response(0.0, numbetween, numintkern);

      /* Free the old kernelarray if one exists */

      if (!firsttime)
         vect_free(kernelarray);

      /* Generate the interpolating kernel array */

      kernelarray = gen_cvect(fftlen);
      place_complex_kernel(rresp, numintkern, kernelarray, fftlen);
      vect_free(rresp);

      /* FFT the kernel array */

      COMPLEXFFT(kernelarray, fftlen, -1);

      /* Set our new static variables */

      old_numkern = numkern;
      old_numbetween = numbetween;
      old_fftlen = fftlen;
      firsttime = 0;
   }

   /* Generate the data array */

   dataarray = gen_cvect(fftlen);
   spread_no_pad(((fcomplex *) data) + beginbin, fftlen / numbetween,
                 dataarray, fftlen, numbetween);
   vect_free(data);

   /* Generate the final response */

   response = gen_cvect(numkern);
   tmpresponse = complex_corr_conv(dataarray, kernelarray, fftlen, FFTD, CORR);

   /* Chop off the contaminated ends and/or the extra data */

   memcpy(response, tmpresponse + numkern / 2, sizeof(fcomplex) * numkern);
   vect_free(tmpresponse);
   vect_free(dataarray);

   return response;
}


fcomplex *gen_bin_response(double roffset, int numbetween, double ppsr,
                           double T, orbitparams * orbit, int numkern)
  /*  Generate the Fourier response function for a sinusoidal PSR      */
  /*  signal from a binary orbit.                                      */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{

   int fftlen, ii, beginbin, numintkern, index, datar;
   int numdata = MIN_NUMDATA;
   int numorbpts = MIN_NUMORBPTS;       /* This should be a power-of-two + 1 */
   float *data;
   double *phi = NULL, startE;
   double amp, f, dt, dtb, t, tp, fpart, ipart, dtemp;
   static int old_numbetween = 0, old_numkern = 0, old_fftlen = 0, firsttime = 1;
   static fcomplex *kernelarray = NULL;
   fcomplex *response, *tmpresponse, *rresp, *dataarray;
   orbitparams orb;

   /* Check that the arguments are OK */

   if (roffset < 0.0 || roffset >= 1.0) {
      printf("\n  roffset = %f (out of bounds) in gen_bin_response().\n\n", roffset);
      exit(-1);
   }
   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_bin_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_bin_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_bin_response().\n\n");
      exit(-1);
   }

   /* Initialize some data */

   datar = numdata / 4;
   if (numkern > datar) {
      numdata = next2_to_n(numkern * 4);
      datar = numdata / 4;
   }
   dt = 1.0 / numdata;
   orb.p = orbit->p / T;
   if (orb.p < 1.0)
      numorbpts = next2_to_n(MIN_NUMORBPTS / orb.p) + 1;
   dtb = 1.0 / (numorbpts - 1);
   amp = 2.0 * dt;
   f = TWOPI * datar;
   orb.x = orbit->x / (ppsr * datar);
   orb.e = orbit->e;
   orb.w = orbit->w * DEGTORAD;
   orb.t = orbit->t / T;

   /* Generate the orbit */

   startE = keplars_eqn(orb.t, orb.p, orb.e, 1.0E-15);
   phi = dorbint(startE, numorbpts, dtb, &orb);
   E_to_phib(phi, numorbpts, &orb);

   /* Generate the data set */

   data = gen_fvect(numdata);
   for (ii = 0; ii < numdata; ii++) {
      t = ii * dt;
      /* The following 4 lines simply linearly interpolate   */
      /* the orbital solution and add the delay to the time. */
      fpart = modf(t / dtb, &ipart);
      index = (int) (ipart + DBLCORRECT);
      dtemp = phi[index];
      tp = t - ((fpart * (phi[index + 1] - dtemp)) + dtemp);
      data[ii] = amp * cos(f * tp);
   }
   vect_free(phi);

   /* FFT the data */

   realfft(data, numdata, -1);
   beginbin = numdata / 4 - numkern / (2 * numbetween);

   /* The following block saves us from having to re-compute */
   /* the Fourier interpolation kernels if 'numkern' is the  */
   /* same length as on prior calls.                         */

   if (numbetween > 1) {
      numintkern = 2 * numbetween * r_resp_halfwidth(HIGHACC);
      fftlen = next2_to_n(numkern + numintkern);
      if (fftlen > numdata) {
         printf("WARNING:  fftlen > numdata in gen_bin_response().\n");
      }
      if (firsttime ||
          old_numkern != numkern ||
          old_numbetween != numbetween || old_fftlen != fftlen) {

         /* Generate an interpolation kernel for the data */

         rresp = gen_r_response(0.0, numbetween, numintkern);

         /* Free the old kernelarray if one exists */

         if (!firsttime)
            vect_free(kernelarray);

         /* Generate the interpolating kernel array */

         kernelarray = gen_cvect(fftlen);
         place_complex_kernel(rresp, numintkern, kernelarray, fftlen);
         vect_free(rresp);

         /* FFT the kernel array */

         COMPLEXFFT(kernelarray, fftlen, -1);

         /* Set our new static variables */

         old_numkern = numkern;
         old_numbetween = numbetween;
         old_fftlen = fftlen;
         firsttime = 0;
      }

      /* Generate the data array */

      dataarray = gen_cvect(fftlen);
      if (fftlen / numbetween >= numdata - beginbin) {
         printf("WARNING:  fftlen too large in gen_bin_response().\n");
      }
      spread_no_pad(((fcomplex *) data) + beginbin, numkern / numbetween,
                    dataarray, fftlen, numbetween);
      vect_free(data);

      /* Generate the final response */

      response = gen_cvect(numkern);
      tmpresponse = complex_corr_conv(dataarray, kernelarray, fftlen, FFTD, CORR);

      /* Chop off the extra data */

      memcpy(response, tmpresponse, sizeof(fcomplex) * numkern);
      vect_free(tmpresponse);
      vect_free(dataarray);

   } else {

      /* Generate the final response */

      response = gen_cvect(numkern);
      memcpy(response, ((fcomplex *) data) + beginbin, sizeof(fcomplex) * numkern);
      vect_free(data);
   }
   return response;
}
