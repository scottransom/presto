#include "presto.h"

static float percolate_pows_and_freqs(float *highpows, float *highfreqs, \
				      int numcands);

void search_minifft(fcomplex *minifft, int numminifft, \
		    float norm, int numcands, float *highpows, \
		    float *highfreqs)
  /* This routine searches a short FFT (usually produced using the   */
  /* MiniFFT binary search method) and returns two vectors which     */
  /* contain the highest powers found and their Fourier frequencies. */
  /* The routine uses interbinning to help find the highest peaks.   */
  /* Arguments:                                                      */
  /*   'minifft' is the FFT to search (complex valued)               */
  /*   'numminifft' is the number of complex points in 'minifft'     */
  /*   'norm' is the value to multiply each pow power by to get      */
  /*      a normalized power spectrum.                               */
  /*   'numcands' is the length of the returned vectors.             */
  /*   'highpows' a vector containing the 'numcands' highest powers. */
  /*   'highfreqs' a vector containing the 'numcands' frequencies    */
  /*      where 'highpows' were found.                               */
  /* Notes:  The returned vectors must have already been allocated.  */
  /*   The returned vectors will be sorted by decreasing power.      */
{
  int ii, numspread, kern_half_width, numkern = 0, numbetween = 2;
  float minpow = 0.0, pwr, powargr, powargi;
  static int firsttime = 1, old_numspread = 0;
  static fcomplex *kernel;
  fcomplex *fftcopy, *spread, *kern;

  /* Copy and then normalize the minifft */

  numspread = numminifft * numbetween;
  fftcopy = gen_cvect(numminifft);
  memcpy(fftcopy, minifft, numminifft * sizeof(fcomplex));
  fftcopy[0].r = 1.0;
  fftcopy[0].i = 1.0;
  for (ii = 0; ii < numminifft; ii++){
    fftcopy[ii].r *= norm;
    fftcopy[ii].i *= norm;
  }

  /* Prep the kernel if needed */

  if (firsttime || old_numspread != numspread){
    kern_half_width = r_resp_halfwidth(LOWACC);
    numkern = 2 * numbetween * kern_half_width;
    kern = gen_r_response(0.0, numbetween, numkern);
    kernel = gen_cvect(numspread);
    place_complex_kernel(kern, numkern, kernel, numspread);
    COMPLEXFFT(kernel, numspread, -1);
    free(kern);
    firsttime = 0;
    old_numspread = numspread;
  }

  /* Interpolate the minifft */
  
  spread = gen_cvect(numspread);
  spread_no_pad(fftcopy, numminifft, spread, numspread, numbetween);
  spread = complex_corr_conv(spread, kernel, numspread, \
			     FFTD, INPLACE_CORR);
  free(fftcopy);

  /* Search the interpolated minifft */

  for (ii = 0; ii < numcands; ii++){
    highpows[ii] = 0.0;
    highfreqs[ii] = 0.0;
  }
  for (ii = 1; ii < numspread; ii++) {
    pwr = POWER(spread[ii].r, spread[ii].i);
    if (pwr > minpow) {
      highpows[numcands-1] = pwr;
      highfreqs[numcands-1] = 0.5 * (float) ii; 
      minpow = percolate_pows_and_freqs(highpows, highfreqs, numcands);
    }
  }
}


static float percolate_pows_and_freqs(float *highpows, float *highfreqs, \
				      int numcands)
  /*  Pushes a power and its corresponding frequency as far up their  */
  /*  respective sorted lists as they shoud go to keep the power list */
  /*  sorted. Returns the new low power in 'highpows'                 */
{
  int ii;
  float tempzz;

  for (ii = numcands - 2; ii >= 0; ii--) {
    if (highpows[ii] < highpows[ii + 1]) {
      SWAP(highpows[ii], highpows[ii + 1]);
      SWAP(highfreqs[ii], highfreqs[ii + 1]);
    } else {
      break;
    }
  }
  return highpows[numcands - 1];
}
