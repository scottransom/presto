#include "presto.h"
#include "search_bin_cmd.h"

/* The number of candidates we would like to examine more closely for    */
/*    every maxfft number of freqs we search.  Used to determine cutoff  */
#define NUMEXPECT 10

/* The number of harmonics to sum in the search */
#define HARMSUM 3

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 10

/* Function definitions */

int not_already_there_bin(binaryprops * binprops, \
			  binaryprops * list, int nlist);
void compare_bin_cands(binaryprops * list, \
		       int nlist, char *output);
int comp_bin_pow(const void *ca, const void *cb);
   /*  Used as compare function for qsort() */
int comp_bin_nfftbins(const void *ca, const void *cb);
   /*  Used as compare function for qsort() */
int remove_dupes_bin(binaryprops * list, int nlist);
   /*  Removes list values that are within 1 Fourier bin of the PSR freq */
   /*  from a higher power candidate. Returns # removed.                 */
int remove_other_bin(binaryprops * list, int nlist);
   /*  Removes list values whose binary parameters are unrealistic.   */
   /*  For example, orbital periods under 200 sec.                    */
   /*  Returns # removed.                                             */


/* Main routine */

int main(int argc, char *argv[])
{
  FILE *fftfile, *candfile;
  float powargr, powargi, *powr, *minifft;
  float highpows[MININCANDS], highfreqs[MININCANDS];
  float norm, numchunks, *fptr1;
  int nbins, newncand, N, nfftsizes, fftlen, halffftlen, binsleft;
  int numtoread, filepos = 0, loopct = 0;
  int ii, ct, newper = 0, oldper = 0, numsumpow = 1;
  double maxpow = 0.0, maxfreq = 0.0, dt, T;
  double totnumsearched = 0.0, minpow = 0.0;
  char filenm[200], candnm[200], binnm[200], *notes;
  fcomplex *data = NULL;
  rderivs derivs;
  fourierprops props;
  binaryprops binprops, *list;
  infodata idata;
  struct tms runtimes;
  double ttim, utim, stim, tott;
  Cmdline *cmd;

  /* Prep out timer */

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
  N = idata.N;
  dt = idata.dt;
  T = N * dt;
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

  if (cmd->rlo < cmd->lobin) cmd->rlo = cmd->lobin;
  if (cmd->rhiP){
    if (cmd->rhi < cmd->rlo) cmd->rhi = cmd->rlo + cmd->maxfft;
    if (cmd->rhi > cmd->lobin + nbins) cmd->rhi = cmd->lobin + nbins;
  } else {
    cmd->rhi = cmd->lobin + nbins;
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

  numtoread = 6 * cmd->maxfft + cmd->maxfft / 2;
  if (cmd->sum == 0){
    data = gen_cvect(numtoread);
  }
  powr = gen_fvect(numtoread);
  minifft = gen_fvect(cmd->maxfft);
  list = (binaryprops *)malloc(sizeof(binaryprops) * cmd->ncand);
  for (ii = 0; ii < cmd->ncand; ii++)
    list[ii].pow = 0.0;

  /* Determine the cutoff level to examine candidates more closely */

  if (cmd->plevP == 0){
    cmd->plev = -log(NUMEXPECT / (double) cmd->maxfft);
  }
  filepos = cmd->rlo - cmd->lobin;
  numchunks = (float) (cmd->rhi - cmd->rlo) / cmd->maxfft;
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

    binsleft = cmd->rhi - cmd->lobin - filepos;
    if (binsleft < cmd->minfft) 
      break;
    if (binsleft < numtoread)     /* 1st try decreasing the numtoread  */
      numtoread = cmd->maxfft + cmd->maxfft / 2;
    if (binsleft <= numtoread) {  /* Next try changing our maxfft size */
      for (cmd->maxfft = cmd->minfft; \
	   cmd->maxfft <= binsleft; \
	   cmd->maxfft *= 2);
      cmd->maxfft /= 2;
      numtoread = cmd->maxfft + cmd->maxfft / 2;
    }
    fftlen = cmd->maxfft;

    /* Read from fftfile */

    if (cmd->sum == 0){
      data = read_fcomplex_file(fftfile, filepos, numtoread);
      for (ii = 0; ii < numtoread; ii++)
	powr[ii] = POWER(data[ii].r, data[ii].i);
    } else {
      powr = read_float_file(fftfile, filepos, numtoread);
    }
    if (filepos == 0) powr[0] = 1.0;
    filepos += numtoread - cmd->maxfft / 2;
      
    /* Chop the powers that are way above the average level */

    prune_powers(powr, numtoread, numsumpow);

    /* Loop through the different small FFT sizes */

    while (fftlen >= cmd->minfft) {

      halffftlen = fftlen / 2;
      fptr1 = powr;

      /* Perform miniffts at each section of the powr array */

      while ((fptr1 - powr) < cmd->maxfft) {

	/* Copy the proper amount and portion of powr into minifft */

	memcpy(minifft, fptr1, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Calculate the normalization constant and search the miniFFT */

	norm = sqrt((double) fftlen * (double) numsumpow) / minifft[0];

	search_minifft((fcomplex *)minifft, halffftlen, norm, \
		       HARMSUM, MININCANDS, highpows, highfreqs);

	/* Check if the measured powers are greater than cmd->plev */

	for (ii = 0; ii < MININCANDS; ii++){
	  
	  if (highpows[ii] > cmd->plev) {
	    maxpow = max_r_arr((fcomplex *)minifft, halffftlen, \
			       highfreqs[ii], &maxfreq, &derivs);

	    /* Check if the measured power is greater than the lowest */
	    /* power we already have in the list.                     */

	    if (maxpow > minpow) {
	      calc_props(derivs, maxfreq, 0.0, 0.0, &props);
	      calc_binprops(&props, T, filepos + cmd->lobin + (fptr1 - powr), \
			    fftlen, &binprops);

	      /* Insure the candidate is semi-realistic */

	      if (binprops.pbin > 300.0) {

		/* Check to see if another candidate with these properties */
		/* is already in the list.                                 */

		if (not_already_there_bin(&binprops, list, cmd->ncand)) {
		  list[cmd->ncand - 1] = binprops;
		  minpow = percolate_bin(list, cmd->ncand);
		}
	      }
	    }
	  }
	  /* Mini-fft search for loop */
	}

	totnumsearched += fftlen;
	fptr1 += halffftlen;

	/* Position of mini-fft in data set while loop */
      }

      fftlen >>= 1;

      /* Size of mini-fft while loop */
    }

    loopct++;

    /* File position while loop */
  }

  /* Print the final percentage update */

  printf("\r   Amount complete = %3d%%\n\n", 100);

  /* Print the number of frequencies searched */

  printf("Searched %12.0f pts (including interbins).\n\n", totnumsearched);

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
  while (ii < cmd->ncand && list[ii].pow != 0)
    ii++;
  newncand = ii;

  /* Sort the results in decreasing power levels */

  qsort(list, (unsigned long) newncand, sizeof(binaryprops), \
	comp_bin_pow);

  /* Remove duplicate candidates */

  newncand -= remove_dupes_bin(list, newncand);

  /* Set our candidate notes to all spaces */

  notes = malloc(sizeof(char) * newncand * 18);
  for (ii = 0; ii < newncand; ii++) {
    sprintf(notes + ii * 18, "                  ");
  }

  /* Check the database for possible known PSR detections */

  if (idata.ra_h && idata.dec_d) {
    for (ii = 0; ii < newncand; ii++) {
      comp_bin_to_cand(&list[ii], &idata, notes + ii * 18, 0);
    }
  }
  /* Compare the candidates with each other */

  compare_bin_cands(list, newncand, notes);

  /* Send the candidates to the text file */

  file_bin_candidates(list, notes, newncand, cmd->argv[0]);

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(list, sizeof(binaryprops), (unsigned long) newncand, candfile);
  fclose(candfile);

  /* Free our arrays and close our files */

  if (cmd->sum == 0) free(data);
  free(list);
  free(powr);
  free(minifft);
  free(notes);
  fclose(fftfile);
  printf("Done.\n\n");
  return (0);
}


int not_already_there_bin(binaryprops * binprops, binaryprops * list, int nlist)
{
  int ii;

  /* Loop through the candidates already in the list */

  for (ii = 0; ii < nlist; ii++) {
    if (list[ii].pow == 0.0)
      break;

    /* Do not add the candidate to the list if it is a lower power */
    /* version of an already listed candidate.                     */

    if (binprops->lowbin == list[ii].lowbin) {
      if (binprops->nfftbins == list[ii].nfftbins) {
	if (fabs(binprops->rdetect - list[ii].rdetect) < 0.6) {
	  if (binprops->pow < list[ii].pow) {
	    return 0;
	  }
	}
      }
    }
  }
  return 1;
}


void compare_bin_cands(binaryprops * list, int nlist, char *notes)
{
  int ii, jj, kk, ll;
  char tmp[30];

  /* Loop through the candidates (reference cands) */

  for (ii = 0; ii < nlist; ii++) {

    /* Loop through the candidates (referenced cands) */

    for (jj = 0; jj < nlist; jj++) {
      if (ii == jj)
	continue;

      /* Loop through the possible PSR period harmonics */

      for (kk = 1; kk < 41; kk++) {

	/* Check if the PSR Fourier freqs are close enough */

	if (fabs(list[ii].rpsr - list[jj].rpsr / kk) < list[ii].nfftbins) {

	  /* Loop through the possible binary period harmonics */

	  for (ll = 1; ll < 10; ll++) {

	    /* Check if the binary Fourier freqs are close enough */

	    if (fabs(list[ii].rbin - list[jj].rbin * ll) < \
		list[jj].rbinerr * 20.0) {

	      /* Check if the note has already been written */

	      sprintf(tmp, "%.18s", notes + jj * 18);
	      if (!strcmp("                  ", tmp)) {

		/* Write the note */

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
