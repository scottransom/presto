#include "presto.h"

/* Number of bins on each side of a freq to use for interpolation */
#define INTERPBINS 5

/* Routines defined at the bottom */

static int padfftlen(int minifftlen, int numbetween, int *padlen);
float percolate_rawbincands(rawbincand *cands, int numcands);
void print_rawbincand(rawbincand cand);
void search_minifft(fcomplex *minifft, int numminifft, \
		    rawbincand *cands, int numcands, int numharmsum, \
		    int numbetween, double numfullfft, double timefullfft, \
		    double lorfullfft, presto_interptype interptype, \
		    presto_checkaliased checkaliased)
  /* This routine searches a short FFT (usually produced using the   */
  /* MiniFFT binary search method) and returns a candidte vector     */
  /* containing information about the best binary candidates found.  */
  /* The routine uses either interbinning or interpolation as well   */
  /* as harmonic summing during the search.                          */
  /* Arguments:                                                      */
  /*   'minifft' is the FFT to search (complex valued)               */
  /*   'numminifft' is the number of complex points in 'minifft'     */
  /*   'cands' is a pre-allocated vector of rawbincand type in which */
  /*      the sorted (in decreasing sigma) candidates are returned   */
  /*   'numcands' is the length of the 'cands' vector                */
  /*   'numharmsum' the number of harmonics to sum during the search */
  /*   'numbetween' the points to interpolate per bin                */
  /*   'numfullfft' the number of points in the original long FFT    */
  /*   'timefullfft' the duration of the original time series (s)    */
  /*   'lorfullfft' the 1st bin of the long FFT that was miniFFT'd   */
  /*   'interptype' is either INTERBIN or INTERPOLATE.               */
  /*      INTERBIN = (interbinning) is fast but less sensitive.      */
  /*      INTERPOLATE = (Fourier interpolation) is slower but more   */
  /*        sensitive.                                               */
  /*   'checkaliased' is either CHECK_ALIASED or NO_CHECK_ALIASED.   */
  /*      NO_CHECK_ALIASED = harmonic summing does not include       */
  /*        aliased freqs making it faster but less sensitive.       */
  /*      CHECK_ALIASED = harmonic summing includes aliased freqs    */
  /*        making it slower but more sensitive.                     */
{
  int ii, jj, fftlen, offset, numtosearch;
  int numspread = 0, kern_half_width, numkern = 0;
  float powargr, powargi, *fullpows = NULL, *sumpows;
  double twobypi, minpow = 0.0, minsig, dr;
  static int firsttime = 1, old_numminifft = 0;
  static fcomplex *kernel;
  fcomplex *spread, *kern;

  /* Override the value of numbetween if interbinning */

  if (interptype == INTERBIN)
    numbetween = 2;

  /* Prep some other values we will need */

  dr = 1.0 / (double) numbetween;
  twobypi = 1.0 / PIBYTWO;
  fftlen = numminifft * numbetween;
  numspread = padfftlen(numminifft, numbetween, &kern_half_width);
  for (ii = 0; ii < numcands; ii++){
    cands[ii].mini_sigma = 0.0;
    cands[ii].mini_power = 0.0;
  }

  /* Prep the interpolation kernel if needed */

  if (interptype == INTERPOLATE){
    if (firsttime || (old_numminifft != numminifft)){
      if (!firsttime) free(kernel);
      numkern = 2 * numbetween * kern_half_width;
      kern = gen_r_response(0.0, numbetween, numkern);
      kernel = gen_cvect(numspread);
      place_complex_kernel(kern, numkern, kernel, numspread);
      COMPLEXFFT(kernel, numspread, -1);
      free(kern);
      firsttime = 0;
      old_numminifft = numminifft;
    }
  }
  
  /* Spread and interpolate the minifft */
  
  spread = gen_cvect(numspread);
  spread_with_pad(minifft, numminifft, spread, numspread, numbetween, 0);
  /* Nyquist is in spread[0].i, but it is usually */
  /* _big_ so we won't use it.                    */
  spread[0].r = spread[fftlen].r = 1.0;
  spread[0].i = spread[fftlen].i = 0.0;
  if (interptype == INTERPOLATE){  /* INTERPOLATE */
    spread = complex_corr_conv(spread, kernel, numspread, \
			       FFTD, INPLACE_CORR);
  } else {                         /* INTERBIN */
    for (ii = 1; ii < fftlen; ii += 2){
      spread[ii].r = twobypi * (spread[ii-1].r - spread[ii+1].r);
      spread[ii].i = twobypi * (spread[ii-1].i - spread[ii+1].i);
    }
  }

  numtosearch = (checkaliased == CHECK_ALIASED) ? 2 * fftlen : fftlen;
  fullpows = gen_fvect(numtosearch);
  fullpows[0] = 1.0;
  fullpows[fftlen] = 1.0;  /* used to be nyquist^2 */

  /* The following wraps the data around the Nyquist freq such that */
  /* we consider aliased frequencies as well (If CHECK_ALIASED).    */
  
  if (checkaliased == CHECK_ALIASED)
    for (ii = 1, jj = numtosearch - 1; ii < fftlen; ii++, jj--)
      fullpows[ii] = fullpows[jj] = POWER(spread[ii].r, spread[ii].i);
  else
    for (ii = 1; ii < numtosearch; ii++)
      fullpows[ii] = POWER(spread[ii].r, spread[ii].i);
  free(spread);

  /* Search the raw powers */

  for (ii = 1; ii < numtosearch; ii++) {
    if (fullpows[ii] > minpow) {
      cands[numcands-1].mini_r = dr * (double) ii; 
      cands[numcands-1].mini_power = fullpows[ii];
      cands[numcands-1].mini_numsum = 1.0;
      cands[numcands-1].mini_sigma = 
	sigma_from_sumpows(fullpows[ii], 1);
      minsig = percolate_rawbincands(cands, numcands);
      minpow = cands[numcands-1].mini_power;
    }
  }

  /* If needed, sum and search the harmonics */
  
  if (numharmsum > 1){
    sumpows = gen_fvect(numtosearch);
    memcpy(sumpows, fullpows, sizeof(float) * numtosearch);
    for (ii = 2; ii <= numharmsum; ii++){
      offset = ii / 2;
      minpow = sumpows_from_sigma(cands[numcands-1].mini_sigma, ii);
      for (jj = 0; jj < numtosearch; jj++){
	sumpows[jj] += fullpows[(jj + offset) / ii];
	if (sumpows[jj] > minpow) {
	  cands[numcands-1].mini_r = dr * (double) jj; 
	  cands[numcands-1].mini_power = sumpows[jj];
	  cands[numcands-1].mini_numsum = (double) ii;
	  cands[numcands-1].mini_sigma = 
	    sigma_from_sumpows(sumpows[jj], ii);
	  minsig = percolate_rawbincands(cands, numcands);
	  minpow = 
	    sumpows_from_sigma(cands[numcands-1].mini_sigma, ii);
	}
      }
    }
    free(sumpows);
  }
  free(fullpows);

  /* Add the rest of the rawbincand data to the candidate array */

  for (ii = 0; ii < numcands; ii++){
    cands[ii].full_N = numfullfft;
    cands[ii].full_T = timefullfft;
    cands[ii].full_lo_r = lorfullfft;
    cands[ii].mini_N = fftlen;
    cands[ii].psr_p = timefullfft / (lorfullfft + numminifft);
    cands[ii].orb_p = timefullfft * cands[ii].mini_r / fftlen;
  }
}


static int padfftlen(int minifftlen, int numbetween, int *padlen)
/* Choose a good (easily factorable) FFT length and an */
/* appropriate padding length (for low accuracy work). */
/* We assume that minifftlen is a power-of-2...        */
{
  int lowaccbins, newlen;

  /* First choose an appropriate number of full pad bins */

  *padlen = minifftlen / 8;
  lowaccbins = r_resp_halfwidth(LOWACC) * (numbetween / 2);
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
  else return (int)((int)((newlen + 1000)/1000) * 1000.0);
}

void print_rawbincand(rawbincand cand){
  printf("  Sigma       =  %-7.3f\n", cand.mini_sigma);
  printf("  Orbit p     =  %-8.2f\n", cand.orb_p);
  if (cand.psr_p < 0.001)
    printf("  Pulsar p    =  %-12.5e\n", cand.psr_p);
  else
    printf("  Pulsar p    =  %-12.9f\n", cand.psr_p);
  printf("  rlo (full)  =  %-10.0f\n", cand.full_lo_r);
  printf("  N (mini)    =  %-6.0f\n", cand.mini_N);
  printf("  r (detect)  =  %-9.3f\n", cand.mini_r);
  printf("  Power       =  %-8.3f\n", cand.mini_power);
  printf("  Numsum      =  %-2.0f\n", cand.mini_numsum);
  printf("  N (full)    =  %-10.0f\n", cand.full_N);
  printf("  T (full)    =  %-13.6f\n\n", cand.full_T);
}

float percolate_rawbincands(rawbincand *cands, int numcands)
  /*  Pushes a rawbincand candidate as far up the array of   */
  /*  candidates as it shoud go to keep the array sorted in  */
  /*  indecreasing significance.  Returns the new lowest     */
  /*  sigma in the array.                                    */
{
  int ii;
  rawbincand tempzz;

  for (ii = numcands - 2; ii >= 0; ii--) {
    if (cands[ii].mini_sigma < cands[ii + 1].mini_sigma) {
      SWAP(cands[ii], cands[ii + 1]);
    } else {
      break;
    }
  }
  return cands[numcands - 1].mini_sigma;
}


