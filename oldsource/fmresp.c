#include "presto.h"

void fmresp(double q, double z, double rb, double w, int tol, \
	    double *ar, double *ai)
/* This routine generates the predicted phase modulation response of a */
/* signal given the frequency offset, q; the modulation amplitude, z;  */
/* the modulation frequency and phase, rb and w, respectively; and     */
/* a tolerance showing the order of magnitude limit between the max    */
/* value and the minimum value to keep summing for a solution.         */
/* Amplitude is returned as a complex value in 'ar', 'ai'.             */
{
  int m, N, i;
  double *bes, tmp, vtheta, vdelta, wtheta, wdelta;
  double vr, vi, va, vb, wr, wi, wa, wb;

  /* call Bessel function routine */

  bes = bessjtable(-z, &N, tol);
  m = -(N - 1) / 2;
  *ar = *ai = 0.0;

  /* Initialize trig recursion for modulation phase (w) term */

  vdelta = w;
  vtheta = m * w;
  tmp = sin(0.5 * vdelta);
  va = -2.0 * tmp * tmp;
  vb = sin(vdelta);
  vr = cos(vtheta);
  vi = sin(vtheta);

  /* Initialize trig recursion for modulation frequency term */

  wdelta = PI * rb;
  wtheta = PI * q + wdelta * m;
  tmp = sin(0.5 * wdelta);
  wa = -2.0 * tmp * tmp;
  wb = sin(wdelta);
  wr = cos(wtheta);
  wi = sin(wtheta);

  for (i = 0; i < N; i++, wtheta += wdelta) {

    /* Avoid divide by zero */

    if (fabs(wtheta) < 0.005) {
      tmp = (1.0 - wtheta * wtheta * 0.16666666666666) * bes[i];
    } else {
      tmp = wi / wtheta * bes[i];
    }

    /* Perform the summation */

    *ar += (wr * vr - wi * vi) * tmp;
    *ai += (wi * vr + wr * vi) * tmp;

    /* Increment our trig recursions */

    tmp = vr;
    vr += tmp * va - vi * vb;
    vi += vi * va + tmp * vb;
    tmp = wr;
    wr += tmp * wa - wi * wb;
    wi += wi * wa + tmp * wb;
  }

  /* Clean-up Bessel function array */

  free(bes);
}



void fmresp_dq(double qlo, double dq, int nq, double z, double rb, \
	       double w, int tol, float resp[])
/* This routine functions just like fmresp above, but it returns an   */
/* array of 'nq' complex values in the float vector 'resp'.  'resp'   */
/* must be defined before calling this routine.                       */
{
  int m, N, N2, i, j, ptrr, ptri, ptrr2, ptri2;
  double *bes, *tmparr, tmp, q;
  double vtheta, vdelta, wtheta, wdelta;
  double va, vb, vr, vi, wa, wb, wr, wi;

  /* Call Bessel function generator */

  bes = bessjtable(-z, &N, tol);
  m = -(N - 1) / 2;
  q = qlo;

  /* Initialize trig recursion for modulation phase term */

  vdelta = w;
  vtheta = m * w;
  tmp = sin(0.5 * vdelta);
  va = -2.0 * tmp * tmp;
  vb = sin(vdelta);
  vr = cos(vtheta);
  vi = sin(vtheta);

  /* Generate array of values that are independent of q */

  N2 = N << 1;
  tmparr = gen_dvect(N2);

  for (i = 0 ; i < N; i++) {
    ptrr = 2 * i;
    ptri = ptrr + 1;
    tmparr[ptrr] = vr * bes[i];
    tmparr[ptri] = vi * bes[i];
    tmp = vr;
    vr += tmp * va - vi * vb;
    vi += vi * va + tmp * vb;
  }

  /* Prep some trig recursion stuff */

  wdelta = PI * rb;
  tmp = sin(0.5 * wdelta);
  wa = -2.0 * tmp * tmp;
  wb = sin(wdelta);

  /* Loop over the q's... */

  for (i = 0 ; i < nq ; i++, q += dq) {
    ptrr = 2 * i;
    ptri = ptrr + 1;
    resp[ptrr] = 0.0;
    resp[ptri] = 0.0;

    /* Initialize trig recursion for modulation frequency term */
    
    wtheta = PI * q + wdelta * m;
    wr = cos(wtheta);
    wi = sin(wtheta);

    /* Compute the summation */

    for (j = 0 ; j < N ; j++, wtheta += wdelta) {
      ptrr2 = 2 * j;
      ptri2 = ptrr2 + 1;

      /* Avoid divide by zero */

      if (fabs(wtheta) < 0.005) {
	tmp = 1.0 - wtheta * wtheta * 0.16666666666666;
      } else {
	tmp = wi / wtheta;
      }
      resp[ptrr] += tmp * (tmparr[ptrr2] * wr - tmparr[ptri2] * wi); 
      resp[ptri] += tmp * (tmparr[ptri2] * wr + tmparr[ptrr2] * wi); 

      /* Increment our trig recursion */
      
      tmp = wr;
      wr += tmp * wa - wi * wb;
      wi += wi * wa + tmp * wb;
    }
  }

  /* Clean-up */

  free(bes);
  free(tmparr);

}


fcomplex *gen_bin_response_math(double roffset, int numbetween, \
				double ppsr, double T, \
				orbitparams * orbit, int numkern, int tol)
  /*  Generate the Fourier response function for a sinusoidal PSR      */
  /*  signal from a binary orbit.                                      */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplarian orbital parameters of the binary system.          */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
  /*    'tol' is the order-of-magnitude difference between the min and */
  /*       max Bessel functions to return during the summation.  See   */
  /*       bessjtable() for more information.                          */
{
  long m2;
  float *response;

  /* Check that roffset and numbetween are in bounds */

  if (roffset < 0.0 || roffset >= 1.0) {
    printf("   roffset out of bounds in gen_bin_response.\n");
    exit(1);
  }
  if (numbetween < 1 || numbetween >= 20000) {
    printf("   numbetween out of bounds in gen_bin_response.\n");
    exit(1);
  }

  /* If we automatically generate the width m */

  if (*m == 0) {
    *m = get_bin_resp_m(numbetween, ppsr, T, orbit);
  }
  m2 = 2 * (*m);
  response = gen_cvect(m2);

  /* Generate the response */

  fmresp_dq(((*m) / (double) numbetween + roffset), \
	    -1.0 / (double) numbetween, \
	    m2, \
	    TWOPI * orbit->x / ppsr, \
	    T / orbit->p, \
	    orbit->w * DEGTORAD, \
	    tol, \
	    response);                   

  return response;
}




