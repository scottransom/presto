#include "presto.h"
#include "multibeam.h"

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 3

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 100.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3.0

/* Bins to ignore at the beginning and end of the big FFT */
#define BINSTOIGNORE 0

/* Factor to overlap the miniFFTs (power-of-two only) */
/*    1 = no overlap                    */
/*    2 = overlap half of each miniFFT  */
/*    4 = overlap 3/4 of each miniFFT   */
#define OVERLAPFACT 4

/* Blocks of maxfft (+ 0.5) to work with at a time */
#define WORKBLOCK 3

int PMsurv_phasemod_search(char *header, int N, fcomplex *bigfft, 
			   double dm, int minfft, int maxfft)
{
  int ii, jj, kk, worklen, fftlen, binsleft;
  int bigfft_pos=0, havecand=0, powers_offset;
  float *powers, *minifft, *powers_pos;
  
  /* Check our input values */

  maxfft = next2_to_n(maxfft);
  minfft = next2_to_n(minfft);
  if (N < minfft){
    printf("\nLength of input array in  PMsurv_phasemod_search()\n");
    printf("is too short or less than 0:  N = %d\n\n", N);
    exit(-1);
  }

  /* Allocate the arrays that will store the powers from */
  /* the bigFFT as well as the miniFFTs.                 */

  worklen = (2 * WORKBLOCK + 1) * (maxfft / 2);
  powers = gen_fvect(worklen);
  minifft = gen_fvect(maxfft);

  /* Loop through the bigFFT */

  while (bigfft_pos < N) {

    /* How close are we to the end of the bigFFT? */
    
    binsleft = N - bigfft_pos;

    /* Adjust our search parameters if close to end of zone to search */

    if (binsleft < worklen){
      worklen = 3 * (maxfft / 2);
      while (binsleft < worklen){
	worklen /= 2;
	maxfft /= 2;
      }
      if (worklen < minfft)
	break;
    }
    fftlen = maxfft;

    /* Get the powers from the bigFFT */

    for (ii = 0, jj = bigfft_pos; ii < worklen; ii++, jj++)
      powers[ii] = POWER(bigfft[jj].r, bigfft[jj].i);
    if (bigfft_pos == 0) powers[0] = 1.0;

    /* Chop the powers that are way above the median.  */
    /* This is a crude way of removing strong coherent */
    /* pulsations or RFI from the power spectrum.      */

    prune_powers(powers, worklen, 1);

    /* Loop through the different small FFT sizes */

    while (fftlen >= minfft) {
      powers_pos = powers;
      powers_offset = 0;

      /* Perform miniffts at each section of the powers array */

      while ((worklen - powers_offset) >
	     (int) ((1.0 - cmd->overlap) * maxfft + DBLCORRECT)){

	/* Copy the proper amount and portion of powers into minifft */

	memcpy(minifft, powers_pos, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Normalize and search the miniFFT */

	norm = sqrt((double) fftlen * (double) numsumpow) / minifft[0];
	for (ii = 0; ii < fftlen; ii++) minifft[ii] *= norm;
	search_minifft((fcomplex *)minifft, fftlen / 2, tmplist, \
		       MININCANDS, cmd->harmsum, cmd->numbetween, idata.N, \
		       T, (double) (powers_offset + bigfft_pos), \
		       INTERPOLATE, NO_CHECK_ALIASED);
		       
	/* Check if the new cands should go into the master cand list */

	for (ii = 0; ii < MININCANDS; ii++){
	  if (tmplist[ii].mini_sigma > minsig) {

	    /* Insure the candidate is semi-realistic */

	    if (tmplist[ii].orb_p > MINORBP && 
		tmplist[ii].mini_r > MINMINIBIN &&
		tmplist[ii].mini_r < tmplist[ii].mini_N - MINMINIBIN &&
		fabs(tmplist[ii].mini_r - 0.5 * tmplist[ii].mini_N) > 2.0){

	      /* Check to see if another candidate with these properties */
	      /* is already in the list.                                 */
	      
	      if (not_already_there_rawbin(tmplist[ii], list, ncand2)) {
		list[ncand2 - 1] = tmplist[ii];
		minsig = percolate_rawbincands(list, ncand2);
	      }
	    } else {
	      continue;
	    }
	  } else {
	    break;
	  }
	  /* Mini-fft search for loop */
	}

	totnumsearched += fftlen;
	powers_pos += (int)(cmd->overlap * fftlen);
	powers_offset = powers_pos - powers;

	/* Position of mini-fft in data set while loop */
      }

      fftlen >>= 1;

      /* Size of mini-fft while loop */
    }

    bigfft_pos += (worklen - (int)((1.0 - cmd->overlap) * maxfft));
    loopct++;

    /* File position while loop */
  }
  free(powers);
}
