#include "presto.h"

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 1

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 100.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3.0

/* Bins to ignore at the beginning and end of the big FFT */
#define BINSTOIGNORE 0

int PMsurv_phasemod(char *header, int N, fcomplex *fft, double dm){
  int ii, jj, kk, fftpos;
  
  /* Loop through fftfile */

  while ((fftpos + cmd->lobin) < cmd->rhi) {

    /* Calculate percentage complete */

    newper = (int) (loopct / numchunks * 100.0);

    if (newper > oldper) {
      newper = (newper > 99) ? 100 : newper;
      printf("\r   Amount complete = %3d%%", newper);
      oldper = newper;
      fflush(stdout);
    }

    /* Adjust our search parameters if close to end of zone to search */

    binsleft = cmd->rhi - (fftpos + cmd->lobin);
    if (binsleft < cmd->minfft) 
      break;
    if (binsleft < numtoread) {  /* Change numtoread */
      numtoread = cmd->maxfft;
      while (binsleft < numtoread){
	cmd->maxfft /= 2;
	numtoread = cmd->maxfft;
      }
    }
    fftlen = cmd->maxfft;

    /* Read from fftfile */

    if (cmd->stack == 0){
      data = read_fcomplex_file(fftfile, fftpos, numtoread);
      for (ii = 0; ii < numtoread; ii++)
	powers[ii] = POWER(data[ii].r, data[ii].i);
      numsumpow = 1;
    } else {
      powers = read_float_file(fftfile, fftpos, numtoread);
      numsumpow = cmd->stack;
    }
    if (fftpos == 0) powers[0] = 1.0;
      
    /* Chop the powers that are way above the median level */

    prune_powers(powers, numtoread, numsumpow);

    /* Loop through the different small FFT sizes */

    while (fftlen >= cmd->minfft) {

      halffftlen = fftlen / 2;
      powers_pos = powers;
      powers_offset = 0;

      /* Perform miniffts at each section of the powers array */

      while ((numtoread - powers_offset) >
	     (int) ((1.0 - cmd->overlap) * cmd->maxfft + DBLCORRECT)){

	/* Copy the proper amount and portion of powers into minifft */

	memcpy(minifft, powers_pos, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Normalize and search the miniFFT */

	norm = sqrt((double) fftlen * (double) numsumpow) / minifft[0];
	for (ii = 0; ii < fftlen; ii++) minifft[ii] *= norm;
	search_minifft((fcomplex *)minifft, halffftlen, tmplist, \
		       MININCANDS, cmd->harmsum, cmd->numbetween, idata.N, \
		       T, (double) (powers_offset + fftpos + cmd->lobin), \
		       cmd->interbinP ? INTERBIN : INTERPOLATE, \
		       cmd->noaliasP ? NO_CHECK_ALIASED : CHECK_ALIASED);
		       
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

    if (cmd->stack == 0) free(data);
    else free(powers);
    fftpos += (numtoread - (int)((1.0 - cmd->overlap) * cmd->maxfft));
    loopct++;

    /* File position while loop */
  }
}
