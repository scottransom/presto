#include "presto.h"
#include "search_bin_cmd.h"

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 6

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 300.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3.0

/* Bins to ignore at the beginning and end of the big FFT */
#define BINSTOIGNORE 0

/* Function definitions */
static int not_already_there_rawbin(rawbincand newcand, 
				    rawbincand *list, int nlist);
static int comp_rawbin_to_cand(rawbincand *cand, infodata * idata,
			       char *output, int full);
static void compare_rawbin_cands(rawbincand *list, int nlist, 
				 char *notes);
static void file_rawbin_candidates(rawbincand *cand, char *notes,
				   int numcands, char name[]);
float percolate_rawbincands(rawbincand *cands, int numcands);

static char num[41][5] =
{"0th", "1st", "2nd", "3rd", "4th", "5th", "6th", \
 "7th", "8th", "9th", "10th", "11th", "12th", \
 "13th", "14th", "15th", "16th", "17th", "18th", \
 "19th", "20th", "21st", "22nd", "23rd", "24th", \
 "25th", "26th", "27th", "28th", "29th", "30th", \
 "31st", "32nd", "33rd", "34th", "35th", "36th", \
 "37th", "38th", "39th", "40th"};

/* Main routine */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *fftfile, *candfile;
  float powargr, powargi, *powers = NULL, *minifft;
  float norm, numchunks, *powers_pos;
  int nbins, newncand, nfftsizes, fftlen, halffftlen, binsleft;
  int numtoread, filepos = 0, loopct = 0, powers_offset, ncand2;
  int ii, ct, newper = 0, oldper = 0, numsumpow = 1;
  double T, totnumsearched = 0.0, minsig = 0.0;
  char filenm[200], candnm[200], binnm[200], *notes;
  fcomplex *data = NULL;
  rawbincand tmplist[MININCANDS], *list;
  infodata idata;
  struct tms runtimes;
  double ttim, utim, stim, tott;
  Cmdline *cmd;

  /* Prep the timer */

  tott = times(&runtimes) / (double) CLK_TCK;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
    Program = argv[0];
    printf("\n");
    usage();
    exit(1);
  }

  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("          Binary Pulsar Search Routine\n");
  printf("              by Scott M. Ransom\n");
  printf("                 23 Sept, 1999\n\n");

  /* Initialize the filenames: */

  sprintf(filenm, "%s.fft", cmd->argv[0]);
  sprintf(candnm, "%s_bin.cand", cmd->argv[0]);
  sprintf(binnm, "%s_bin", cmd->argv[0]);

  /* Read the info file */

  readinf(&idata, cmd->argv[0]);
  T = idata.N * idata.dt;
  if (idata.object) {
    printf("Analyzing %s data from '%s'.\n\n", idata.object, filenm);
  } else {
    printf("Analyzing data from '%s'.\n\n", filenm);
  }

  /* open the FFT file and get its length */

  fftfile = chkfopen(filenm, "rb");
  nbins = chkfilelen(fftfile, sizeof(fcomplex));

  /* Check that cmd->maxfft is an acceptable power of 2 */

  ct = 4;
  ii = 1;
  while (ct < MAXREALFFT || ii) {
    if (ct == cmd->maxfft)
      ii = 0;
    ct <<= 1;
  }
  if (ii) {
    printf("\n'maxfft' is out of range or not a power-of-2.\n\n");
    exit(1);
  }

  /* Check that cmd->minfft is an acceptable power of 2 */

  ct = 4;
  ii = 1;
  while (ct < MAXREALFFT || ii) {
    if (ct == cmd->minfft)
      ii = 0;
    ct <<= 1;
  }
  if (ii) {
    printf("\n'minfft' is out of range or not a power-of-2.\n\n");
    exit(1);
  }

  /* Low and high Fourier freqs to check */

  if (!cmd->rloP || (cmd->rloP && cmd->rlo < cmd->lobin)){
    if (cmd->lobin == 0)
      cmd->rlo = BINSTOIGNORE;
    else
      cmd->rlo = cmd->lobin;
  }
  if (cmd->rhiP){
    if (cmd->rhi < cmd->rlo) cmd->rhi = cmd->rlo + cmd->maxfft;
    if (cmd->rhi > cmd->lobin + nbins) 
      cmd->rhi = cmd->lobin + nbins - BINSTOIGNORE;
  } else {
    cmd->rhi = cmd->lobin + nbins - BINSTOIGNORE;
  }

  /* Determine how many different mini-fft sizes we will use */

  nfftsizes = 1;
  ii = cmd->maxfft;
  while (ii > cmd->minfft) {
    ii >>= 1;
    nfftsizes++;
  }

  /* Allocate some memory and prep some variables.             */
  /* For numtoread, the 6 just lets us read extra data at once */

  numtoread = 6 * cmd->maxfft;
  if (cmd->stack == 0)
    powers = gen_fvect(numtoread);
  minifft = gen_fvect(cmd->maxfft);
  ncand2 = 2 * cmd->ncand;
  list = (rawbincand *)malloc(sizeof(rawbincand) * ncand2);
  for (ii = 0; ii < ncand2; ii++)
    list[ii].mini_sigma = 0.0;
  for (ii = 0; ii < MININCANDS; ii++)
    tmplist[ii].mini_sigma = 0.0;
  filepos = cmd->rlo - cmd->lobin;
  numchunks = (float) (cmd->rhi - cmd->rlo) / numtoread;
  printf("Searching...\n");
  printf("   Amount complete = %3d%%", 0);
  fflush(stdout);

  /* Loop through fftfile */

  while ((filepos + cmd->lobin) < cmd->rhi) {

    /* Calculate percentage complete */

    newper = (int) (loopct / numchunks * 100.0);

    if (newper > oldper) {
      newper = (newper > 99) ? 100 : newper;
      printf("\r   Amount complete = %3d%%", newper);
      oldper = newper;
      fflush(stdout);
    }

    /* Adjust our search parameters if close to end of zone to search */

    binsleft = cmd->rhi - (filepos + cmd->lobin);
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
      data = read_fcomplex_file(fftfile, filepos, numtoread);
      for (ii = 0; ii < numtoread; ii++)
	powers[ii] = POWER(data[ii].r, data[ii].i);
      numsumpow = 1;
    } else {
      powers = read_float_file(fftfile, filepos, numtoread);
      numsumpow = cmd->stack;
    }
    if (filepos == 0) powers[0] = 1.0;
      
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
		       T, (double) (powers_offset + filepos + cmd->lobin), \
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
    filepos += (numtoread - (int)((1.0 - cmd->overlap) * cmd->maxfft));
    loopct++;

    /* File position while loop */
  }

  /* Print the final percentage update */

  printf("\r   Amount complete = %3d%%\n\n", 100);

  /* Print the number of frequencies searched */

  printf("Searched %.0f pts (including interbins).\n\n", totnumsearched);

  printf("Timing summary:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n", \
	 ttim, utim, stim);
  printf("  Total time: %.3f sec\n\n", tott);

  printf("Writing result files and cleaning up.\n");

  /* Count how many candidates we actually have */

  ii = 0;
  while (ii < ncand2 && list[ii].mini_sigma != 0)
    ii++;
  newncand = (ii > cmd->ncand) ? cmd->ncand : ii;

  /* Set our candidate notes to all spaces */

  notes = malloc(sizeof(char) * newncand * 18 + 1);
  for (ii = 0; ii < newncand; ii++)
    strncpy(notes + ii * 18, "                     ", 18);
  
  /* Check the database for possible known PSR detections */
  
  if (idata.ra_h && idata.dec_d) {
    for (ii = 0; ii < newncand; ii++) {
      comp_rawbin_to_cand(&list[ii], &idata, notes + ii * 18, 0);
    }
  }

  /* Compare the candidates with each other */

  compare_rawbin_cands(list, newncand, notes);

  /* Send the candidates to the text file */

  file_rawbin_candidates(list, notes, newncand, cmd->argv[0]);

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(list, sizeof(rawbincand), (unsigned long) newncand, candfile);
  fclose(candfile);

  /* Free our arrays and close our files */

  if (cmd->stack == 0)
    free(powers);
  free(list);
  free(minifft);
  free(notes);
  fclose(fftfile);
  printf("Done.\n\n");
  return (0);
}


static int not_already_there_rawbin(rawbincand newcand, 
				    rawbincand *list, int nlist)
{
  int ii;

  /* Loop through the candidates already in the list */

  for (ii = 0; ii < nlist; ii++) {
    if (list[ii].mini_sigma == 0.0)
      break;

    /* Do not add the candidate to the list if it is a lower power */
    /* version of an already listed candidate.                     */

    if (list[ii].mini_N == newcand.mini_N) {
      if (fabs(list[ii].mini_r - newcand.mini_r) < 0.6) {
	if (list[ii].mini_sigma > newcand.mini_sigma) {
	  return 0;
	}
      }
    }
  }
  return 1;
}


static void compare_rawbin_cands(rawbincand *list, int nlist, 
				 char *notes)
{
  double perr;
  int ii, jj, kk, ll;
  char tmp[30];

  /* Loop through the candidates (reference cands) */

  for (ii = 0; ii < nlist; ii++) {

    /* Loop through the candidates (referenced cands) */

    for (jj = 0; jj < nlist; jj++) {
      if (ii == jj)
	continue;
      perr = 0.5 * list[jj].full_T / list[jj].mini_N;

      /* Loop through the possible PSR period harmonics */
      
      for (kk = 1; kk < 41; kk++) {
	
	/* Check if the PSR Fourier freqs are close enough */

	if (fabs(list[ii].full_lo_r - list[jj].full_lo_r / kk) < 
	    list[ii].mini_N) {

	  /* Loop through the possible binary period harmonics */

	  for (ll = 1; ll < 10; ll++) {

	    /* Check if the binary Fourier freqs are close enough */

	    if (fabs(list[ii].orb_p - list[jj].orb_p / ll) < perr) {

	      /* Check if the note has already been written */

	      sprintf(tmp, "%.18s", notes + jj * 18);
	      if (!strcmp("                  ", tmp)) {

		/* Write the note */

		if (ll == 1 && kk == 1)
		  sprintf(notes + jj * 18, "Same as #%d?", ii + 1);
		else 
		  sprintf(notes + jj * 18, "MH=%d H=%d of #%d", ll, kk, ii + 1);

		break;
	      }
	    }
	  }
	}
      }
    }
  }
}


static int comp_rawbin_to_cand(rawbincand *cand, infodata * idata,
			       char *output, int full)
  /* Compares a binary PSR candidate defined by its props found in    */
  /*   *cand, and *idata with all of the pulsars in the pulsar        */
  /*   database.  It returns a string (verbose if full==1) describing */
  /*   the results of the search in *output.                          */
{
  int i, j, k;
  static int np = 0;
  static psrdatabase pdata;
  double T, theop, ra, dec, beam2, difft = 0.0, epoch;
  double bmod, pmod, orbperr, psrperr;
  char tmp1[80], tmp2[80], tmp3[80], tmp4[20];

  /* If calling for the first time, read the database. */

  if (!np) {
    np = read_database(&pdata);
  }
  /* Convert the beam width to radians */

  beam2 = 2.0 * ARCSEC2RAD * idata->fov;
  
  /* Convert RA and DEC to radians  (Use J2000) */

  ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
  dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);

  /* Calculate the time related variables  */

  T = idata->N * idata->dt;
  epoch = (double) idata->mjd_i + idata->mjd_f;
  
  /* Calculate the approximate error in our value of orbital period */
  
  orbperr = 0.5 * cand->full_T / cand->mini_N;
  
  /* Calculate the approximate error in our value of spin period */

  if (cand[k].full_lo_r == 0.0)
    psrperr = cand[k].psr_p;
  else 
    psrperr = fabs(cand->full_T / 
		   (cand->full_lo_r + 0.5 * cand->mini_N) -
		   cand->full_T / cand->full_lo_r);

  /* Run through RAs in database looking for things close  */
  /* If find one, check the DEC as well (the angle between */
  /* the sources < 2*beam diam).  If find one, check its   */
  /* period.  If this matches within 2*perr, return the    */
  /* number of the pulsar.  If no matches, return 0.       */

  for (i = 0; i < np; i++) {
    
    /* See if we're close in RA */
    
    if (fabs(pdata.ra2000[i] - ra) < 5.0 * beam2) {
      
      /* See if we're close in RA and DEC */
      
      if (sphere_ang_diff(pdata.ra2000[i], pdata.dec2000[i], \
			  ra, dec) < 5.0 * beam2) {

	/* Check that the psr in the database is in a binary   */

	if (pdata.ntype[i] & 8) {

	  /* Predict the period of the pulsar at the observation MJD */

	  difft = SECPERDAY * (epoch - pdata.epoch[i]);
	  theop = pdata.p[i] + pdata.pdot[i] * difft;
	  sprintf(tmp4, "%.8s", pdata.bname + i * 8);

	  /* Check the predicted period and its harmonics against the */
	  /* measured period.  Use both pulsar and binary periods.    */

	  for (j = 1, pmod = 1.0; j < 41; j++, pmod = 1.0 / (double) j) {
	    if (fabs(theop * pmod - cand->psr_p) < psrperr) {
	      for (k = 1, bmod = 1.0; k < 10; \
		   k++, bmod = 1.0 / (double) k) {
		if (fabs(pdata.pb[i] * bmod - cand->orb_p / SECPERDAY) < orbperr) {
		  if (!strcmp("        ", tmp4)) {
		    if (j > 1) {
		      if (full) {
			sprintf(tmp1, "Possibly the %s phasemod harmonic ", num[k]);
			sprintf(tmp2, "of the %s harmonic of PSR ", num[j]);
			sprintf(tmp3, "J%.12s (p = %11.7f s, pbin = %9.4f d).\n", \
				pdata.jname + i * 12, theop, pdata.pb[i]);
			sprintf(output, "%s%s%s", tmp1, tmp2, tmp3);
		      } else {
			sprintf(output, "%s H J%.12s", num[k], pdata.jname + i * 12);
		      }
		    } else {
		      if (full) {
			sprintf(tmp1, "Possibly the %s phasemod harmonic ", num[k]);
			sprintf(tmp2, "of PSR J%.12s (p = %11.7f s, pbin = %9.4f d).\n", \
				pdata.jname + i * 12, theop, pdata.pb[i]);
			sprintf(output, "%s%s", tmp1, tmp2);
		      } else {
			sprintf(output, "PSR J%.12s", pdata.jname + i * 12);
		      }
		    }
		  } else {
		    if (j > 1) {
		      if (full) {
			sprintf(tmp1, "Possibly the %s modulation harmonic ", num[k]);
			sprintf(tmp2, "of the %s harmonic of PSR ", num[j]);
			sprintf(tmp3, "B%s (p = %11.7f s, pbin = %9.4f d).\n", \
				tmp4, theop, pdata.pb[i]);
			sprintf(output, "%s%s%s", tmp1, tmp2, tmp3);
		      } else {
			sprintf(output, "%s H B%s", num[k], tmp4);
		      }
		    } else {
		      if (full) {
			sprintf(tmp1, "Possibly the %s phasemod harmonic ", num[k]);
			sprintf(tmp2, "of PSR B%s (p = %11.7f s, pbin = %9.4f d).\n", \
				tmp4, theop, pdata.pb[i]);
			sprintf(output, "%s%s", tmp1, tmp2);
		      } else {
			sprintf(output, "PSR B%s", tmp4);
		      }
		    }
		  }
		}
		return i + 1;
	      }
	    }
	  }
	}
      }
    }
  }

  /* Didn't find a match */

  if (full) {
    sprintf(output, "I don't recognize this candidate in the pulsar database.\n");
  } else {
    sprintf(output, "                  ");
  }
  return 0;
}


static void file_rawbin_candidates(rawbincand *cand, char *notes,
				   int numcands, char name[])
/* Outputs a .ps file describing all the binary candidates from a    */
/*   binary search. */
{
  FILE *fname;
  int i, j, k = 0;
  int nlines = 87, pages, extralines, linestoprint;
  char filenm[100], infonm[100], command[200];
  double orbperr, psrperr;
  
  sprintf(filenm, "%s_bin", name);
  sprintf(infonm, "%s.inf", name);
  fname = chkfopen(filenm, "w");

  if (numcands <= 0) {
    printf(" Must have at least 1 candidate in ");
    printf("file_bin_candidates().\n\n");
    exit(1);
  }
  pages = numcands / nlines + 1;
  extralines = numcands % nlines;

  for (i = 1; i <= pages; i++) {

    /*                       1         2         3         4         5         6         7         8         9         0         1    */  
    /*              123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234*/
    fprintf(fname, "#               P_orbit +/- Error   P_pulsar +/- Error   FullFFT   MiniFFT   MiniFFT  Num   Sum                   \n");
    fprintf(fname, "# Cand  Sigma         (sec)                (sec)         Low Bin   Length      Bin    Sum  Power  Notes           \n");
    fprintf(fname, "#------------------------------------------------------------------------------------------------------------------\n");
    
    if (i == pages) {
      linestoprint = extralines;
    } else {
      linestoprint = nlines;
    }
    
    for (j = 0; j < linestoprint; j++, k++) {
      
      /* Calculate the approximate error in our value of orbital period */
      orbperr = 0.5 * cand[k].full_T / cand[k].mini_N;
      
      /* Calculate the approximate error in our value of spin period */

      if (cand[k].full_lo_r == 0.0)
	psrperr = cand[k].psr_p;
      else 
	psrperr = fabs(cand[k].full_T / (cand[k].full_lo_r + 
					 0.5 * cand[k].mini_N) -
		       cand[k].full_T / cand[k].full_lo_r);
      
      /*  Now output it... */

      fprintf(fname, " %4d %7.3f  ", k + 1, cand[k].mini_sigma);
      fprintf(fname, " %8.2f", cand[k].orb_p);
      fprintf(fname, " %-7.2g ", orbperr);
      if (cand[k].psr_p < 0.001)
	fprintf(fname, " %12.5e", cand[k].psr_p);
      else
	fprintf(fname, " %12.9f", cand[k].psr_p);
      fprintf(fname, " %-7.2g ", psrperr);
      fprintf(fname, " %9.0f  ", cand[k].full_lo_r);
      fprintf(fname, " %6.0f ", cand[k].mini_N);
      fprintf(fname, " %8.1f ", cand[k].mini_r);
      fprintf(fname, " %2.0f ", cand[k].mini_numsum);
      fprintf(fname, "%7.2f ", cand[k].mini_power);
      fprintf(fname, " %.18s\n", notes + k * 18);
      fflush(fname);
    }
  }
  fprintf(fname, "\n Notes:  MH = Modulation harmonic.  ");
  fprintf(fname, "H = Pulsar harmonic.  # indicates the candidate number.\n\n");
  fclose(fname);
  sprintf(command, "cat %s >> %s", infonm, filenm);
  system(command);
  sprintf(command, \
	  "$PRESTO/bin/a2x -c1 -n90 -title -date -num %s > %s.ps", \
	  filenm, filenm);
  system(command);
}

