#include "presto.h"
#include "multibeam.h"

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 3

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 100.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3.0

/* Bins to ignore at the beginning and end of the bigFFT */
#define BINSTOIGNORE 0

/* Number of candidates in the stack per bigFFT */
#define NUMCANDS 10

/* Factor to overlap the miniFFTs (power-of-two only) */
/*    1 = no overlap of successive miniFFTs           */
/*    2 = overlap 1/2 of each miniFFT                 */
/*    4 = overlap 3/4 of each miniFFT                 */
#define OVERLAPFACT 4

/* Blocks of length maxfft to work with at a time */
#define WORKBLOCK 4

/* Number of harmonics to sum in each miniFFT   */
/* Usually this is greater than 1, but for the  */
/* PMsurv, 1 is probably fine since we will be  */
/* pushing the limit on orbital periods (which  */
/* determines the number of harmonics present   */
/* along with the eccentricity -- which we      */
/* expect would be 0 in a 20 minute binary).    */
#define NUMHARMSUM 1

/* Level of Fourier interpolation to use in the search      */
/*     1 = no interpolation (just the standard FFT bins)    */
/*     2 = interbinning (1 interpolated point between bins) */
/*     4 = 3 interpolated points between bins               */
#define NUMBETWEEN 2

/* If the following is defined, the program will print */
/* debugging info to stdout.                           */
/*#define DEBUGOUT*/

/* Output threshold */
/* If any candidates have a sigma greater than the following  */
/* then print the candidate and observation info to the       */
/* permanent candfile.                                        */
#define THRESH 10.0

/* Candidate file directory */
#define OUTDIR "/home/ransom"

/* Function definitions */
int not_already_there_rawbin(rawbincand newcand, 
 			    rawbincand *list, int nlist);
int comp_rawbin_to_cand(rawbincand *cand, infodata * idata,
 		       char *output, int full);
void compare_rawbin_cands(rawbincand *list, int nlist, 
 			 char *notes);
void file_rawbin_candidates(rawbincand *cand, char *notes,
				   int numcands, char name[]);
float percolate_rawbincands(rawbincand *cands, int numcands);

/******************************************************************/

int PMsurv_phasemod_search(char *header, int N, fcomplex *bigfft, 
			   float dm, int minfft, int maxfft)
{
  int ii, jj, worklen, fftlen, binsleft, overlaplen;
  int bigfft_pos=0, powers_offset, wrkblk=WORKBLOCK;
  float *powers, *minifft, *powers_pos, powargr, powargi;
  double T, norm, minsig=0.0;
  rawbincand list[NUMCANDS], tmplist[MININCANDS];
  multibeam_tapehdr *hdr=NULL;
  infodata idata;
  
  /* Prep our candidate lists */

  for (ii = 0; ii < NUMCANDS; ii++)
    list[ii].mini_sigma = 0.0;
  for (ii = 0; ii < MININCANDS; ii++)
    tmplist[ii].mini_sigma = 0.0;

  /* Copy the header information into *hdr */
    
  memcpy(hdr, header, HDRLEN);

  /* Convert the Header into usable info... */

  multibeam_hdr_to_inf(hdr, &idata);
  T = N * idata.dt;
#ifdef DEBUGOUT
      print_multibeam_hdr(hdr);
#endif

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

  worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
  powers = gen_fvect(worklen);
  minifft = gen_fvect(maxfft);

  /* Loop through the bigFFT */

  while (bigfft_pos < N) {

    /* How close are we to the end of the bigFFT? */
    
    binsleft = N - bigfft_pos;

    /* Adjust our search parameters if close to end of zone to search */

    if (binsleft < worklen){
      wrkblk = 1;
      worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
      while (binsleft < worklen){
	maxfft /= 2;
	worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
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
#ifdef DEBUGOUT
      printf("  powers_offset = %d\n", powers_offset);
#endif
      overlaplen = fftlen / OVERLAPFACT;

      /* Perform miniffts at each section of the powers array */

      while (powers_offset < wrkblk * maxfft){

	/* Copy the proper amount and portion of powers into minifft */

	memcpy(minifft, powers_pos, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Normalize and search the miniFFT */

	norm = sqrt((double) fftlen) / minifft[0];
	for (ii = 0; ii < fftlen; ii++) minifft[ii] *= norm;
	search_minifft((fcomplex *)minifft, fftlen / 2, tmplist, \
		       MININCANDS, NUMHARMSUM, NUMBETWEEN, (double) N, \
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
	      
	      if (not_already_there_rawbin(tmplist[ii], list, NUMCANDS)) {
		list[NUMCANDS - 1] = tmplist[ii];
		minsig = percolate_rawbincands(list, NUMCANDS);
	      }
	    } else {
	      continue;
	    }
	  } else {
	    break;
	  }
	  /* Mini-fft search for loop */
	}

	powers_pos += overlaplen;
	powers_offset += overlaplen;

	/* Position of mini-fft in data set while loop */
      }

      /* Switch to the next smaller miniFFT size */

      fftlen /= 2;
#ifdef DEBUGOUT
      printf("\nfftlen = %d\n\n", fftlen);
#endif

      /* Size of mini-fft while loop */
    }

    bigfft_pos += wrkblk * maxfft;
#ifdef DEBUGOUT
      printf("\nbigfft_pos = %d\n\n", bigfft_pos);
#endif

    /* BigFFT position while loop */
  }

  /* Free up our data arrays */

  free(powers);
  free(minifft);

  /* Check to see if our best candidate is worth saving */

#ifdef DEBUGOUT
  {
#else
  if (list[0].mini_sigma > THRESH){
#endif
    int newnumcands;
    char *notes, outfilenm[100];
    FILE *outfile;

    sprintf(idata.name, "%s/%s", OUTDIR, hdr->pname);
    strncpy(idata.observer, "PMSurv", 90);
    strncpy(idata.analyzer, 
	    "Scott Ransom <ransom@cfa.harvard.edu>", 90);
    idata.dm = dm;
    writeinf(&idata);

    /* Count how many candidates we actually have */
    
    ii = 0;
    while (ii < NUMCANDS && list[ii].mini_sigma != 0)
      ii++;
    newnumcands = (ii > NUMCANDS) ? NUMCANDS : ii;
    
    /* Set our candidate notes to all spaces */
    
    notes = malloc(sizeof(char) * newnumcands * 18 + 1);
    for (ii = 0; ii < newnumcands; ii++)
      strncpy(notes + ii * 18, "                     ", 18);
    
    /* Compare the candidates with each other */
    
    compare_rawbin_cands(list, newnumcands, notes);
    
    /* Send the candidates to the text file */
    
    file_rawbin_candidates(list, notes, newnumcands, idata.name);
    
    /* Write the binary candidate file */
    
    outfile = chkfopen(outfilenm, "wb");
    chkfwrite(list, sizeof(rawbincand), newnumcands, outfile);
    fclose(outfile);
    free(notes);
    return newnumcands;
  }
  return 0;
}

#undef MININCANDS
#undef MINORBP
#undef MINMINIBIN
#undef BINSTOIGNORE
#undef NUMCANDS
#undef OVERLAPFACT
#undef WORKBLOCK
#undef NUMHARMSUM
#undef NUMBETWEEN
#undef THRESH
#ifdef DEBUGOUT
#undef DEBUGOUT
#endif
