#include "presto.h"

/* The number of neighboring bins on each side of a frequency to use for */
/*    Fourier interpolation during the search phase                      */
#define SEARCHM   10

/* The number of candidates we would like to examine more closely for    */
/*    every maxfft number of freqs we search.  Used to determine cutoff  */
#define NUMEXPECT 10

/* Function definitions */

void usage(void);
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
  float powargr, powargi;
  float *powr, *minifft;
  float norm, numchunks, *fptr1;
  fcomplex *kernels[30], *data, *resp, *spread;
  unsigned long nbins;
  int ii, jj, ct, kern_half_width = SEARCHM, numkern;
  int numbetween = 2, ncand, newncand, N, nfftsizes;
  int rlo = 0, rhi, lofreq = 0;
  int minfft = 32, maxfft, fftlen;
  int numtoread, filepos = 0, ltmp, loopct = 0;
  int newper = 0, oldper = 0, numsumpow = 1;
  double maxpow = 0.0, maxfreq = 0.0, dt, T, cutoff;
  double totnumsearched = 0.0, minpow = 0.0;
  char filenm[100], candnm[100], binnm[100], *notes;
  rderivs derivs;
  fourierprops props;
  binaryprops binprops, *list;
  infodata idata;
  struct tms runtimes;
  double ttim, utim, stim, tott;

  tott = times(&runtimes) / (double) CLK_TCK;

  if (argc < 4 || argc > 9) {
    usage();
    exit(1);
  }
  printf("\n\n");
  printf("          Binary Pulsar Search Routine\n");
  printf("              by Scott M. Ransom\n");
  printf("                 23 June, 1999\n\n");

  /* Initialize the filenames: */

  sprintf(filenm, "%s.fft", argv[1]);
  sprintf(candnm, "%s_bin.cand", argv[1]);
  sprintf(binnm, "%s_bin", argv[1]);

  /* Read the info file */

  readinf(&idata, argv[1]);
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
  rhi = nbins;

  /* Read the other command line parameters */

  /* The number of candidates to record */

  ncand = atoi(argv[2]);
  if ((ncand < 1) || (ncand > 5000)) {
    printf("\n'ncand' is out of range.\n\n");
    exit(1);
  }
  /* The largest FFT to use in the search */

  maxfft = atoi(argv[3]);
  if ((maxfft <= 2) || ((unsigned long) maxfft > nbins - 1)) {
    printf("\n'maxfft' is out of range.\n\n");
    exit(1);
  }

  /* Check that maxfft is an acceptable power of 2 */

  ct = 4;
  ii = 0;
  while (ct < MAXREALFFT || !ii) {
    if (ct == maxfft)
      ii = 1;
    ct <<= 1;
  }
  if (!ii) {
    printf("\n'maxfft' is not a power of 2.\n\n");
    exit(1);
  }

  /* The smallest FFT to use in the search */

  if (argc >= 5) {
    minfft = atoi(argv[4]);

    /* Check that minfft is in range */

    if ((minfft <= 2) || ((unsigned long) minfft > nbins - 1)) {
      printf("\n'minfft' is out of range.\n\n");
      exit(1);
    }
    /* Check that minfft is an acceptable power of 2 */

    ct = 4;
    ii = 0;
    while (ct < MAXREALFFT || !ii) {
      if (ct == maxfft)
	ii = 1;
      ct <<= 1;
    }
    if (!ii) {
      printf("\n'maxfft' is out of range.\n\n");
      exit(1);
    }
  }
  /* The lowest freq present in the FFT file */

  if (argc >= 6) {
    lofreq = atoi(argv[5]);
    if ((lofreq < 0) || ((unsigned long) lofreq > nbins - 1)) {
      printf("\n'lofreq' is out of range.\n\n");
      exit(1);
    }
  }
  /* rlo and rhi */

  if (argc >= 7) {
    rlo = atoi(argv[6]);
    if ((rlo < lofreq) || ((unsigned long) rlo > nbins - 1)) {
      printf("\nLow frequency to search (rlo) is out of range.\n\n");
      exit(1);
    }
  }
  if (argc >= 8) {
    rhi = atoi(argv[7]);
    if ((rhi < rlo) || ((unsigned long) rhi > nbins)) {
      printf("\nHigh frequency to search (rhi) is out of range.\n\n");
      exit(1);
    }
  }
  /* Is the original FFT a sum of other FFTs with the amplitudes added */
  /* in quadrature?  (i.e. an incoherent sum)                          */

  if (argc >= 9) {
    numsumpow = atoi(argv[8]);
    if (numsumpow < 1) {
      printf("\nNumber of summed powers must be at least one.\n\n");
      exit(1);
    }
  }
  /* Determine how many different mini-fft sizes we will use */

  nfftsizes = 1;
  ii = maxfft;
  while (ii > minfft) {
    ii >>= 1;
    nfftsizes++;
  }

  /* Allocate and initialize our Fourier interpolation kernel arrays */

  numkern = 2 * numbetween * kern_half_width;
  resp = gen_r_response(0.0, numbetween, numkern);

  for (ii = 0; ii < nfftsizes; ii++) {
    jj = maxfft / (1 << ii);
    kernels[ii] = gen_cvect(jj);
    place_complex_kernel(resp, numkern, kernels[ii], jj);
    COMPLEXFFT(kernels[ii], jj, -1);
  }

  /* Allocate some memory */

  numtoread = maxfft + maxfft / 2;
  data = gen_cvect(numtoread);
  powr = gen_fvect(numtoread);
  minifft = gen_fvect(maxfft);
  spread = gen_cvect(maxfft);
  list = (binaryprops *)malloc(sizeof(binaryprops) * ncand);
  for (ii = 0; ii < ncand; ii++)
    list[ii].pow = 0.0;

  /* Determine the cutoff level to examine candidates more closely */

  cutoff = -log(NUMEXPECT / (double) maxfft);
  filepos = rlo - lofreq;
  numchunks = (float) (rhi - rlo) / maxfft;
  printf("Searching...\n");
  printf("   Amount complete = %3d%%", 0);
  fflush(stdout);

  /* Loop through fftfile */

  while ((filepos + lofreq) < rhi) {

    /* Calculate percentage complete */

    newper = (int) (loopct / numchunks * 100.0);

    if (newper > oldper) {
      newper = (newper > 99) ? 100 : newper;
      printf("\r   Amount complete = %3d%%", newper);
      oldper = newper;
      fflush(stdout);
    }
    /* Count which mini-FFT size we are on */

    ct = 0;

    /* Adjust our search parameters if close to end of zone to search */

    if ((ltmp = rhi - filepos - lofreq) <= numtoread) {
      maxfft = minfft;
      ct = nfftsizes - 1;
      while (2 * maxfft <= ltmp) {
	maxfft <<= 1;
	ct--;
      }
      numtoread = maxfft + maxfft / 2;
      if (ltmp >= numtoread) {
      } else {
	numtoread -= maxfft / 2;
      }
    }
    fftlen = maxfft;

    /* Read from fftfile */

    data = read_fcomplex_file(fftfile, filepos, numtoread);

    for (ii = 0; ii < numtoread; ii++)
      powr[ii] = POWER(data[ii].r, data[ii].i);
    if (filepos == 0) powr[0] = 1.0;

    /* Chop the powers that are way above the average level */

    prune_powers(powr, numtoread, numsumpow);

    /* Loop through the different small FFT sizes */

    while (fftlen >= minfft) {

      fptr1 = powr;

      /* Perform miniffts at each section of the powr array */

      while ((fptr1 - powr) < maxfft) {

	/* Copy the proper amount and portion of powr into minifft */

	memcpy(minifft, fptr1, fftlen * sizeof(float));

	/* Perform the minifft */

	realfft(minifft, fftlen, -1);

	/* Calculate the normalization constant */

	norm = sqrt((double) fftlen * (double) numsumpow) / minifft[0];

	/* Now normalize the miniFFT */

	minifft[0] = 1.0;
	minifft[1] = 1.0;
	for (ii = 2; ii < fftlen; ii++)
	  minifft[ii] *= norm;

	/* Interpolate the minifft */

	spread_no_pad((fcomplex *)minifft, fftlen / 2, spread, 
		      fftlen, numbetween);
	COMPLEXFFT(spread, fftlen, -1);
	spread = complex_corr_conv(spread, kernels[ct], \
				   fftlen, NOFFTS, INPLACE_CORR);

	/* Search the interpolated minifft */

	for (ii = 0; ii < fftlen; ii++) {

	  /* Check if the measured power is greater than cutoff */

	  if (POWER(spread[ii].r, spread[ii].i) > cutoff) {
	    maxpow = max_r_arr((fcomplex *)minifft, fftlen / 2, \
			       ii / (double) numbetween, \
			       &maxfreq, &derivs);
	    /* Check if the measured power is greater than the lowest */
	    /* power we already have in the list.                     */

	    if (maxpow > minpow) {
	      calc_props(derivs, maxfreq, 0.0, 0.0, &props);
	      calc_binprops(&props, T, filepos + lofreq + (fptr1 - powr), \
			    fftlen, &binprops);

	      /* Insure the candidate is semi-realistic */

	      if (binprops.pbin > 300.0) {

		/* Check to see if another candidate with these properties */
		/* is already in the list.                                 */

		if (not_already_there_bin(&binprops, list, ncand)) {
		  list[ncand - 1] = binprops;
		  minpow = percolate_bin(list, ncand);
		}
	      }
	    }
	  }
	  /* Mini-fft search for loop */
	}

	totnumsearched += fftlen;
	fptr1 += fftlen / 2;
	/* Position of mini-fft in data set while loop */
      }

      fftlen >>= 1;
      ct++;
      /* Size of mini-fft while loop */
    }

    filepos += maxfft;
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
  while (ii < ncand && list[ii].pow != 0)
    ii++;
  newncand = ii;

  /* Sort the results:  first by power, second by fftlen */

  qsort(list, (unsigned long) newncand, sizeof(binaryprops), \
	comp_bin_pow);
  qsort(list, (unsigned long) newncand, sizeof(binaryprops), \
	comp_bin_nfftbins);

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

  file_bin_candidates(list, notes, newncand, argv[1]);

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(list, sizeof(binaryprops), (unsigned long) newncand, candfile);
  fclose(candfile);

  /* Free our arrays and close our files */

  for (ii = 0; ii < nfftsizes; ii++)
    free(kernels[ii]);
  free(list);
  free(data);
  free(powr);
  free(minifft);
  free(notes);
  free(resp);
  free(spread);
  if (idata.onoff) free(idata.onoff);
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


void usage(void)
{
  printf("\nUsage:  search_bin file ncand maxfft [minfft] ");
  printf("[lofreq] [rlo] [rhi] [sum]\n\n");
  printf("  Mandatory arguments:\n");
  printf("       'file' = a string containing the FFT file's name.\n");
  printf("                You must have an '.inf' file of the same\n");
  printf("                name as well.  Do not add a suffix.\n");
  printf("                Candidates will be returned in a file called\n");
  printf("                'file_bin.cand'.  A Postscript format\n");
  printf("                candidate list will be in 'file_bin.ps'.\n");
  printf("      'ncand' = (int) The routine will return 'ncand' candidates.\n");
  printf("                Must be less than or equal to 5000.\n");
  printf("     'maxfft' = (int) This is the maximum length of the short\n");
  printf("                FFTs that will be used to search 'file'.\n");
  printf("  Optional arguments:\n");
  printf("     'minfft' = (int) This is the minimum length of the short\n");
  printf("                FFTs that will be used to search 'file'.\n");
  printf("                Default size is 32 points.\n");
  printf("     'lofreq' = (int) Lowest Fourier bin in FFT file.  This is\n");
  printf("                useful if you chop a long FFT for space reasons.\n");
  printf("                If 'lofreq' is present and not equal to '0', \n");
  printf("                we will assume the 0th frequency bin = 1 for\n");
  printf("                normalization purposes.\n");
  printf("        'rlo' = (int) lowest Fourier bin to search.\n");
  printf("        'rhi' = (int) highest Fourier bin to search.\n\n");
  printf("        'sum' = Number of summed FFTs (added in quadrature)\n");
  printf("                that make up the original FFT to search.\n\n");
  printf("  'search_bin' will search the Fourier transformed time \n");
  printf("  series found in 'file' for pulsars which are Doppler shifted\n");
  printf("  due to an orbital companion to such an extent that their\n");
  printf("  spread-out Fourier power can be treated as the frequency\n");
  printf("  modulated response of a signal at the intrinsic spin \n");
  printf("  frequency of the pulsar.  The frequency of modulation is\n");
  printf("  simply the orbital frequency.\n");
  printf("  The routine outputs formatted statistics for the 'ncand'\n");
  printf("  best candidates found in the search region.  If the \n");
  printf("  optional arguments are ommitted, the routine will search\n");
  printf("  the whole FFT file and assume the first bin is freq=0.\n");
  printf("                                        23 June 1999\n\n");
}
