#include "presto.h"

/* Number of bins on each side of a freq to use for interpolation */
#define INTERPBINS 5

int padfftlen(int minifftlen, int numbetween, int *padlen);
static float percolate_pows_and_freqs(float *highpows, float *highfreqs, \
				      int numcands);

void search_minifft(fcomplex *minifft, int numminifft, float norm, \
		    int harmsum, int numcands, float *highpows, \
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
  /*   'harmsum' the number of harmonics to sum during the search.   */
  /*   'numcands' is the length of the returned vectors.             */
  /*   'highpows' a vector containing the 'numcands' highest powers. */
  /*   'highfreqs' a vector containing the 'numcands' frequencies    */
  /*      where 'highpows' were found.                               */
  /* Notes:  The returned vectors must have already been allocated.  */
  /*   The returned vectors will be sorted by decreasing power.      */
{
  int ii, numspread, kern_half_width, numkern = 0;
  int numbetween = 2, lowaccbins;
  float minpow = 0.0, pwr, powargr, powargi, sqrtnorm, nyquist;
  static int firsttime = 1, old_numminifft = 0, old_khw = 0;
  static fcomplex *kernel;
  fcomplex *fftcopy, *spread, *kern;

  /* Copy and then normalize the minifft */
  
  fftcopy = gen_cvect(numminifft);
  memcpy(fftcopy, minifft, numminifft * sizeof(fcomplex));
  sqrtnorm = sqrt(norm);
  nyquist = fftcopy[0].i * sqrtnorm;
  fftcopy[0].r = 1.0;
  fftcopy[0].i = 1.0;
  for (ii = 0; ii < numminifft; ii++){
    fftcopy[ii].r *= sqrtnorm;
    fftcopy[ii].i *= sqrtnorm;
  }

  /* Prep the kernel if needed */

  if (firsttime || old_numminifft != numminifft){
    if (!firsttime) free(kernel);
    numspread = padfftlen(numminifft, numbetween, &kern_half_width);
    numkern = 2 * numbetween * kern_half_width;
    kern = gen_r_response(0.0, numbetween, numkern);
    kernel = gen_cvect(numspread);
    place_complex_kernel(kern, numkern, kernel, numspread);
    COMPLEXFFT(kernel, numspread, -1);
    free(kern);
    firsttime = 0;
    old_numminifft = numminifft;
  }

  /* Interpolate the minifft */
  
  spread = gen_cvect(numspread);
  spread_no_pad(fftcopy, numminifft, spread, numspread, numbetween);
  spread = complex_corr_conv(spread, kernel, numspread, \
			     FFTD, INPLACE_CORR);
  free(fftcopy);

  /* Prep the arrays for harmonic summing */

  if (harmsum > 1){
    fullpows = gen_fvect(2 * numminifft + 1);
    fullpows[0] = 1.0;
    fullpows[2 * numminifft] = nyquist;
    for (ii = 0, jj = 2 * numminifft - 1; 
	 ii < numminifft / 2; 
	 ii++, jj--){
      
    }

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


int padfftlen(int minifftlen, int numbetween, int *padlen)
/* Choose a good (easily factorable) FFT length and an */
/* appropriate padding length (for low accuracy work). */
/* We assume that minifftlen is a power-of-2...        */
{
  int lowaccbins, newlen;

  /* First choose an appropriate number of full pad bins */

  *padlen = minifftlen / 8;
  lowaccbins = r_resp_halfwidth(LOWACC);
  if (*padlen > lowaccbins) *padlen = lowaccbins;

  /* Now choose the FFT length (This requires an FFT that */
  /* can perform non-power-of-two FFTs -- USE FFTW!!!     */

  newlen = (minifftlen + *padlen) * numbetween;

  if (newlen <= 144) return newlen;
  else if (newlen <= 288) return 288;
  else if (newlen <= 540) return 540;
  else if (newlen <= 1080) return 1080;
  else if (newlen <= 2100) return 2100;
  else if (newlen <= 4200) return 4200;
  else if (newlen <= 8232) return 8232;
  else if (newlen <= 16464) return 16464;
  else if (newlen <= 32805) return 32805;
  else if (newlen <= 65610) return 65610;
  else if (newlen <= 131220) return 131220;
  else if (newlen <= 262440) return 262440;
  else if (newlen <= 525000) return 525000;
  else if (newlen <= 1050000) return 1050000;
  /* The following might get optimized out and give garbage... */
  else return int(int((newlen + 1000)/1000) * 1000.0);
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
