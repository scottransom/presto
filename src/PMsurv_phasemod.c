#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tape_hdr.h"

typedef struct FCOMPLEX {
  float r; /* real */
  float i; /* imag */
} fcomplex;

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 100.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3.0

/* Minimum pulsar frequency (Hz) to search */
#define MINFREQ 1.0

/* Maximum pulsar frequency (Hz) to search */
#define MAXFREQ 900.0

/* Factor to overlap the miniFFTs (power-of-two only) */
/*    1 = no overlap of successive miniFFTs           */
/*    2 = overlap 1/2 of each miniFFT                 */
/*    4 = overlap 3/4 of each miniFFT                 */
#define OVERLAPFACT 2

/* Blocks of length maxfft to work with at a time */
#define WORKBLOCK 4

/* If the following is defined, the program will print        */
/* debugging info to stdout or the results to the output dir. */
/*#define DEBUGOUT*/
#define CANDOUT

/* Candidate file directory */
#define OUTDIR "/home/ransom"
/*#define OUTDIR "/tmp_mnt/usr/users/sransom/results"*/

/* The number of miniffts to use  */
#define NUMMINIFFTS 7

/* The lengths of the miniffts to use */
static int minifftlen[NUMMINIFFTS] = \
  {64, 128, 256, 512, 1024, 2048, 4096};

/* Output thresholds */
/* If any candidates have a power level greater than the   */
/* following (which corresponds to a sigma of ~7.5 when    */
/* you take into account the number of bins searched in    */
/* each minifft), then print the candidate and observation */
/* info to the output candfile.                            */
static int threshold[NUMMINIFFTS] = \
  {34.54, 35.23, 35.93, 36.62, 37.31, 38.01, 38.70};

/* The binary candidate structure to write if we find candidates */

typedef struct RAWBINCAND {
  double full_N;       /* Number of points in original time series  */
  double full_T;       /* Length (s) of original time series        */
  double full_lo_r;    /* Lowest Fourier bin that was miniFFTd      */
  double mini_N;       /* Number of points in short (mini) FFT      */
  double mini_r;       /* Candidate Fourier bin in miniFFT          */
  double mini_power;   /* Candidate normalized power in miniFFT     */
  double mini_numsum;  /* Number of powers summed to get candidate  */
  double mini_sigma;   /* Equivalent candidate sigma (for sum pow)  */
  double psr_p;        /* Approx PSR period (miniFFT center bin)    */
  double orb_p;        /* Approx orbital period (s)                 */
} rawbincand;

/* A short structure containing useful info */

typedef struct INFO {
  double dt;
  double dm;
  int file_cntr;
  int ibeam;
  char tape_lbl[7];
  char pname[17];
  char outfilebase[50];
} info;  

void info_from_header(struct tphdr *header, double dt, double dm, 
		      info *idata)
{
  idata->dt = dt;
  idata->dm = dm;
  strncpy(idata->tape_lbl, header->tape_lbl, 6);
  idata->tape_lbl[6] = '\0';
  strncpy(idata->pname, header->pname, 16);
  idata->pname[16] = '\0';
  idata->ibeam = strtol(hdr->ibeam, NULL, 10);
  idata->file_cntr = strtol(hdr->file_cntr, NULL, 10);
  sprintf(idata->outfilebase, "%s_%d_Beam%2d_DM%.1f_%s", 
	  idata->tape_lbl, idata->file_cntr, idata->ibeam,
	  idata->idata->dm, idata->pname);
}


/******************************************************************/

int PMsurv_phasemod_search(struct tphdr *header, int N, 
			   fcomplex *bigfft, double dt, double dm)
/*
 * This routine searches an FFT (assumed to be single precision        
 * complex values in increasing order by Fourier Freq) from the        
 * Parkes Multibeam Pulsar Survey for binary pulsars in very short     
 * orbits.  The method used looks for periodic phase-modulation        
 * induced sidelobes around a pulsar spin freq in the power spectrum   
 * (see http://cfa160.harvard.edu/~ransom/readable_poster.ps.gz for    
 * more details).  The inputs are:                                     
 *   *header - the 640 byte raw header for the observation             
 *   N - the number of frequencies in the input FFT 
 *       (should be 2**22 for a full length pointing)  
 *   *bigfft - the single precision full-length FFT to be searched     
 *   dt - the sample time for the FFT (should be 0.00025 for
 *       a non-decimated pointing)
 *   dm - the DM used to prepare the data in the bigfft                
 * The routine as now stands takes approximately 10-15 sec to run
 * on a full-length FFT (i.e. a 2**22 point pointing).  
 * We should be sensitive to pulsars in very-low to very-high mass     
 * binaries with orbital periods <~ 20min and flux densities of ~1 mJy  
 * or a little less (if the duty cycle of the pulsar is short enough). 
 * Return value is 0 if the code didn't find any candidates. If the 
 * return value is positive, then a significant candidate (or more) 
 * was found in the FFT (this could be the trigger to save the bigFFT 
 * to a file for a more in depth analysis...)
 */
{
  int ii, jj, worklen, fftlen, binsleft, overlaplen;
  int bigfft_offset, powers_offset, wrkblk=WORKBLOCK;
  float *powers, *minifft, *powers_pos, powargr, powargi;
  double T, norm, dr=0.5;
  rawbincand cand;
  PKMB_tapehdr *hdr;
  info idata;
  
  /* Convert the Header into usable info... */

  info_from_header(header, dt, dm, &idata);
  T = N * dt;

  /* Allocate the arrays that will store the powers from */
  /* the bigFFT as well as the miniFFTs.                 */

  worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
  powers = (float *)malloc(sizeof(float) * worklen);
  minifft = (float *)malloc(sizeof(float) * maxfft);

  /* Loop through the bigFFT */

  while (bigfft_offset < N) {

    /* How close are we to the end of the bigFFT? */
    
    binsleft = N - bigfft_offset;

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

    for (ii = 0, jj = bigfft_offset; ii < worklen; ii++, jj++)
      powers[ii] = POWER(bigfft[jj].r, bigfft[jj].i);
    if (bigfft_offset == 0) powers[0] = 1.0;

    /* Chop the powers that are way above the median.  */
    /* This is a crude way of removing strong coherent */
    /* pulsations or RFI from the power spectrum.      */

    prune_powers(powers, worklen, 1);

    /* Loop through the different small FFT sizes */

    while (fftlen >= minfft) {
      powers_pos = powers;
      powers_offset = 0;
      overlaplen = fftlen / OVERLAPFACT;

      /* Perform miniffts at each section of the powers array */

      while (powers_offset < wrkblk * maxfft){

#ifdef DEBUGOUT
	printf("  %d:  offset = %d\n", fftlen, 
	       bigfft_offset + powers_offset);
#endif
	/* Copy the proper amount and portion of powers into minifft */

	memcpy(minifft, powers_pos, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Normalize and search the miniFFT */

	norm = sqrt((double) fftlen) / minifft[0];
	for (ii = 0; ii < fftlen; ii++) minifft[ii] *= norm;
	search_minifft((fcomplex *)minifft, fftlen / 2, 
		       min_orb_p, max_orb_p, tmplist,
		       MININCANDS, NUMHARMSUM, 2, idata.N, T,
		       (double) (powers_offset + bigfft_offset),
		       INTERBIN, NO_CHECK_ALIASED);
		       
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
	    } else
	      continue;
	  } else
	    break;

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

    bigfft_offset += wrkblk * maxfft;
#ifdef DEBUGOUT
      printf("\nbigfft_offset = %d\n\n", bigfft_offset);
#endif

    /* BigFFT position while loop */
  }

  /* Free up our data arrays */

  free(powers);
  free(minifft);

  /* Check to see if our best candidate is worth saving */

#ifdef CANDOUT
  {
#else
  if (list[0].mini_sigma > THRESH){
#endif
    int newnumcands;
    char *notes, outfilenm[100];
    FILE *outfile;

    sprintf(idata.name, "%s/%s", OUTDIR, idata.object);
    strtofilename(idata.name);
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
    
    sprintf(outfilenm, "%s_bin.cand", idata.name);
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
#undef THRESH
#ifdef CANDOUT
#undef CANDOUT
#endif
#ifdef DEBUGOUT
#undef DEBUGOUT
#endif
