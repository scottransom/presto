#include "presto.h"

#define SHORTESTFFT 32768	/* The shortest data set we will look at */
#define LOSKIP      200		/* The min # of low freq bins to skip    */

/* To do:  - Make an MPI version.                                        */
/*           Should allow the saving of specific output files that       */
/*           can be combined by some other routine later.                */
/*         - Make a more intelligent FFT size selection routine          */
/*           It should be based on the number of z's to search, the      */
/*           amount of the correlation we have to throw away, and        */
/*           the speeds of the FFTs.  Compare points searched/sec.       */

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

/* Define a couple functions specific to this program */

void compare_rzw_cands(fourierprops * list, int nlist, char *notes);
int not_already_there_rzw(position * newpos, position * list, int nlist);
int compare_fourierprops(const void *ca, const void *cb);
/*  Used as compare function for qsort() */
int remove_dupes(position * list, int nlist);
/*  Removes list values that are 1 unit of search away from a higher */
/*  power candidate (dr = 0.5 or dz = 2.0).  Returns # removed.      */
int remove_dupes2(fourierprops * list, int nlist);
/*  Removes list values that are within measurement error away from  */
/*  a higher power candidate.  Returns # removed.                    */
int remove_other(fourierprops * list, int nlist, long rlo, \
		 long rhi, double locpow);
/*  Removes list values whose frequencies fall outside rlo and rhi   */
/*  and candidates whose local power levels are below locpow.        */ 
/*  Returns # removed.                                               */


int main(int argc, char *argv[])
{
  FILE *fftfile, *candfile, *poscandfile;
  double dt, nph, T, N, bigz, hir, hiz, dz = 2.0;
  double powavg, powsdev, powvar, powskew, powkurt;
  float powargr, powargi, locpow = 1.0, *powlist;
  float powdiff, hipowchop, lowpowlim;
  float chkpow = 0.0, hipow = 0.0, minpow = 0.0, numr = 0.0;
  fcomplex *response, **kernels, *corrdata, *filedata; 
  unsigned long totnumsearched = 0, nreal;
  int numbetween = 2, startbin, lofreq = 0.0, lobinskip, hibinskip;
  int nbins, nr = 1, nz, corrsize = 0, mincorrsize, worknumbins = 0;
  int rlo, rhi, zlo = 0, zhi = 0, numkern, kern_half_width;
  int ii, ct, zct, nextbin, filedatalen;
  int ncand, newncand, oldper = 0, newper = 0;
  char filenm[200], candnm[200], poscandnm[200], *notes;
  char rzwnm[200];
  presto_datainf datainf;
  position *list, newpos;
  rderivs *derivs;
  fourierprops *props;
  infodata idata;
  struct tms runtimes;
  double ttim, utim, stim, tott;

  tott = times(&runtimes) / (double) CLK_TCK;
  if ((argc < 5) || (argc > 8)) {
    printf("\nUsage:  search_rzw filename ncand zlo zhi [lofreq] [rlo] [rhi]\n\n");
    printf("  Mandatory arguments:\n");
    printf("   'filename' = a string containing the FFT file's name.\n");
    printf("                You must have an '.inf' file of the same\n");
    printf("                name as well.  Do not add a suffix.\n");
    printf("                Candidates will be returned in a file called\n");
    printf("                'filename_rzw.cand'.  A Postscript format\n");
    printf("                candidate list will be in 'filename_rzw.ps'.\n");
    printf("      'ncand' = (int) The routine will return 'ncand' candidates.\n");
    printf("                Must be less than or equal to 5000.\n");
    printf("        'zlo' = (int) lowest Fourier freq deriv to search.\n");
    printf("        'zhi' = (int) highest Fourier freq deriv to search.\n\n");
    printf("  Optional arguments:\n");
    printf("     'lofreq' = (int) Lowest Fourier bin in FFT file.  This is\n");
    printf("                useful if you chop a long FFT for space reasons.\n");
    printf("                If 'lofreq' is present and not equal to '0', \n");
    printf("                we will assume the 0th frequency bin = 1 for\n");
    printf("                normalization purposes.\n");
    printf("        'rlo' = (int) lowest Fourier bin to search.\n");
    printf("        'rhi' = (int) highest Fourier bin to search.\n\n");
    printf("  'search_rzw' will search a region of the f-fdot plane for \n");
    printf("  pulsations in a file containing a long, single precision FFT\n");
    printf("  using the Correlation method (i.e. Ransom and \n");
    printf("  Eikenberry, 1997, unpublished as of yet).\n");
    printf("  The search uses a spacing of 0.5 frequency bins in\n");
    printf("  the fourier frequency (r) direction and 2 'bins' in\n");
    printf("  the fdot (z) direction.\n");
    printf("  The routine outputs formatted statistics for the 'ncand'\n");
    printf("  best candidates found in the search region.  If the\n");
    printf("  optional arguments are ommitted, the routine will search\n");
    printf("  the whole FFT file and assume the first bin is freq=0.\n\n");
    printf("  The routine currently cannot search the 'w' dimension.\n");
    printf("  This may be fixed shortly.\n");
    printf("                                        7 Nov 1997\n\n");
    exit(0);
  }
  printf("\n\n");
  printf("       Pulsar Acceleration Search Routine\n");
  printf("              by Scott M. Ransom\n");
  printf("                  7 Nov, 1997\n\n");

  /* Initialize the input filename: */

  sprintf(filenm, "%s.fft", argv[1]);

  /* Read the info file */

  readinf(&idata, argv[1]);
  if (idata.object) {
    printf("Analyzing %s data from '%s'.\n\n", idata.object, filenm);
  } else {
    printf("Analyzing data from '%s'.\n\n", filenm);
  }

  /* open the FFT file and get its length */

  fftfile = chkfopen(filenm, "rb");
  nph = get_numphotons(fftfile);
  nreal = chkfilelen(fftfile, sizeof(float));

  /* # of fourier frequencies */

  nbins = nreal >> 1;
  if (nreal < SHORTESTFFT) {
    printf("\nFFT is too short to use this routine.\n\n");
    exit(1);
  }
  rewind(fftfile);

  /* zlo and zhi */

  zlo = atoi(argv[3]);
  zhi = atoi(argv[4]);

  /* insure we have a spacing of dz = 2 */

  if ((zhi - zlo) & 1)
    zhi++;
  if ((zlo < -2000000) || (zhi > 2000000)) {
    printf("\nFrequency derivatives to search are out of range.\n\n");
    exit(1);
  }
  bigz = DMAX(fabs((double) zlo), fabs((double) zhi));
  kern_half_width = z_resp_halfwidth(bigz, LOWACC);
  nz = (int) ((zhi - zlo) / dz) + 1;

  /* Initialize the input filename: */

  sprintf(candnm, "%s_rzw_z:%d_%d.cand", argv[1], zlo, zhi);
  sprintf(poscandnm, "%s_rzw_z:%d_%d.pos", argv[1], zlo, zhi);
  sprintf(rzwnm, "%s_rzw_z:%d_%d", argv[1], zlo, zhi);

  /* Skip the lowest LOSKIP bins do to low frequency errors  */
  /* unless you specifically give rlo on the command line.   */
  /* Do the same at the high freqs.                          */

  lobinskip = bigz / 4;
  hibinskip = lobinskip;
  if (lobinskip < LOSKIP)
    lobinskip = LOSKIP;
  rlo = lobinskip;
  rhi = nbins - 1 - hibinskip;

  /* Check other arguments */

  /* The number of candidates to save */

  ncand = atoi(argv[2]);
  ncand = (int) (ncand * 1.5);
  if ((ncand < 1) || (ncand > 15000)) {
    printf("\n'ncand' must be >= 1 and <= 10000.\n\n");
    exit(1);
  }
  /* Determine the correlation sizes we will use: */

  /* The following mincorrsize ensures that we never waste more */
  /* than 25% of any correlation's points due to end effects.   */
  /* (Note: We throw away 2*m points from a correlation)        */

  mincorrsize = 6 * kern_half_width * numbetween;

  /* We need to keep nz + 2 arrays of corrsize complex numbers in */
  /* memory at all times. (nz kernels, 1 data set, 1 result)      */
  /* Therefore worknumbins is the largest corrsize we can have.   */

  worknumbins = (int) ((MAXREALFFT >> 1) / (nz + 2));

  /* Determine corrsize */

  corrsize = next2_to_n(worknumbins);

  /* Insure smaller than worknumbins then *0.5 to be sure... */

  corrsize >>= 2;

  if (mincorrsize > corrsize) {
    printf("\nYou are asking for too much memory.  Specify\n");
    printf("fewer z values to search, or lower z values.\n\n");
    exit(1);
  }
  /* Generate the correlation kernels */

  printf("Generating fdot kernels for the correlations...\n");

  kernels = gen_cmatrix(nz, corrsize);
  numkern = 2 * numbetween * kern_half_width;
  for (ii = 0; ii < nz; ii++) {
    response = gen_z_response(0.0, numbetween, zlo + ii * dz, numkern);
    place_complex_kernel(response, numkern, kernels[ii], corrsize);
    free(response);
    COMPLEXFFT(kernels[ii], corrsize, -1);
  }

  printf("Done generating kernels.\n\n");

  /* The lowest freq present in the FFT file */

  if (argc >= 6) {
    lofreq = atoi(argv[5]);
    if ((lofreq < 0) || (lofreq > nbins - 1)) {
      printf("\n'lofreq' is out of range.\n\n");
      exit(1);
    }
    if (lofreq != 0)
      nph = 1.0;
  }
  /* rlo and rhi */

  if (argc >= 7) {
    rlo = atoi(argv[6]);
    if ((rlo < lofreq) || (rlo > nbins - 1)) {
      printf("\nLow frequency to search (rlo) is out of range.\n\n");
      exit(1);
    }
  }
  if (argc == 8) {
    rhi = atoi(argv[7]);
    if ((rhi < rlo) || (rhi > nbins - 1)) {
      printf("\nHigh frequency to search (rhi) is out of range.\n\n");
      printf("\nFrequencies to search are out of range.\n\n");
      exit(1);
    }
  }
  /* Allocate some memory */

  list = malloc(sizeof(position) * ncand);
  derivs = malloc(sizeof(rderivs) * ncand);
  props = malloc(sizeof(fourierprops) * ncand);
  corrdata = gen_cvect(corrsize);

  dt = idata.dt;
  N = idata.N;
  T = N * dt;
  numr = (rhi - rlo + 1) * nz * 2;
  filedatalen = corrsize / numbetween;

  /* We will automatically get rid of any candidates that have local */
  /* powers much lower than what we would expect to find in the      */
  /* search from purely statistical reasons alone                    */
  /* Note:  6.95 is the full width in z of a signal response         */

  lowpowlim = -0.8 * log(ncand / (numr * 0.5 * dz / 6.95));

  /* Initialize the candidate list */

  for (ii = 0; ii < ncand; ii++) {
    list[ii].pow = 0.0;
    list[ii].p1 = 0.0;
    list[ii].p2 = 0.0;
    list[ii].p3 = 0.0;
  }

  /* Start the main search loop */

  nextbin = rlo;

  do {

    startbin = (unsigned long) nextbin;

    /* Get the data from the file */

    filedata = read_fcomplex_file(fftfile, startbin - kern_half_width, \
				  filedatalen);

    /*  Do the f-fdot plane correlations: */

    for (zct = 0; zct < nz; zct++) {

      /* Calculate percentage complete */

      newper = (int) (totnumsearched / (numr * 0.01)) + 1;

      if (newper > oldper) {
	newper = (newper > 99) ? 100 : newper;
	printf("\rAmount of search complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }

      if (zct == 0) datainf = RAW;
      else datainf = SAME;

      /* Perform the correlation */

      nr = corr_complex(filedata, filedatalen, datainf, \
			kernels[zct], corrsize, FFT, \
			corrdata, corrsize, kern_half_width, \
			numbetween, kern_half_width, CORR);
      nextbin = startbin + nr / numbetween;

      if (zct == 0) {

	/* Get approximate local power statistics */

	powlist = gen_fvect(nr);
	worknumbins = (nextbin > rhi) ? \
	  (rhi - startbin) * numbetween : nr;
	for (ii = 0; ii < worknumbins; ii++) 
	  powlist[ii] = POWER(corrdata[ii].r, corrdata[ii].i);
	stats(powlist, worknumbins, &powavg, &powvar, &powskew, &powkurt);
	powsdev = sqrt(powvar);

	/* Throw powers away that are greater than 4.0 sdev above mean  */
	/* Also throw powers away that are very small.  This could show */
	/* that we have come into a zero-padded area.                   */

	ct = worknumbins;
	locpow = powavg * worknumbins;
	hipowchop = 4.0 * powsdev + powavg;
	for (ii = 0; ii < worknumbins; ii++) {
	  powdiff = powlist[ii] - powavg;
	  if (powdiff > hipowchop) {
	    locpow -= powlist[ii];
	    ct--;
	  }
	}
	locpow /= (double) ct;
	free(powlist);
      }

      /* This loop is the heart of the search */

      for (ii = 0; ii < worknumbins; ii++) {
	chkpow = POWER(corrdata[ii].r, corrdata[ii].i) / locpow;

	/* Check if the measured power is greater than cutoff */

	if (chkpow > minpow) {
	  newpos.pow = chkpow;
	  newpos.p1 = startbin + 0.5 * ii;
	  if (newpos.p1 > rhi)
	    break;
	  newpos.p2 = zlo + zct * dz;
	  newpos.p3 = 0.0;

	  /* Check to see if another candidate with these properties */
	  /* is already in the list.                                 */

	  if (not_already_there_rzw(&newpos, list, ncand)) {
	    list[ncand - 1] = newpos;
	    minpow = percolate(list, ncand);
	  }
	}
      }
      totnumsearched += worknumbins;
    }
    free(filedata);
  } while (nextbin <= rhi);

  /* Free the memory used by the correlation kernels */

  free(kernels[0]);
  free(kernels);

  printf("\rAmount of search complete = %3d%%", 100);
  fflush(stdout);
  printf("\nDone searching.  ");
  printf("Now optimizing each candidate and sorting.\n\n");

  /* Do rough duplicate removal (probably not necessary) */

  newncand = ncand - remove_dupes(list, ncand);

  /* Save the list of 'rough' candidates to a file */

  poscandfile = chkfopen(poscandnm, "w");
  chkfwrite(list, sizeof(position), newncand, poscandfile);
  fclose(poscandfile);

  /* Now maximize each good candidate */

  newper = 0;
  oldper = 0;

  for (ii = 0; ii < newncand; ii++) {

    /* Calculate percentage complete */

    newper = (int) (ii / (float) (newncand) * 100.0) + 1;
    if (newper > oldper) {
      printf("\rAmount of optimization complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
    }
    hipow = max_rz_file(fftfile, list[ii].p1, list[ii].p2, \
			&hir, &hiz, &derivs[ii]);
    calc_props(derivs[ii], hir + lofreq, hiz, 0.0, &props[ii]);
  }
  printf("\rAmount of optimization complete = %3d%%\n\n", 100);

  qsort(props, (unsigned long) newncand, sizeof(fourierprops), \
	compare_fourierprops);

  /* Do fine scale duplicate removal and other cleaning */

  newncand -= remove_dupes2(props, newncand);
  newncand -= remove_other(props, newncand, rlo, rhi, lowpowlim);

  /* Set our candidate notes to all spaces */

  notes = malloc(sizeof(char) * newncand * 20);
  for (ii = 0; ii < newncand; ii++) {
    sprintf(notes + ii * 20, "                    ");
  }

  /* Compare the candidates with the pulsar database */

  if (idata.ra_h && idata.dec_d) {
    for (ii = 0; ii < newncand; ii++) {
      comp_psr_to_cand(&props[ii], &idata, notes + ii * 20, 0);
    }
  }

  /* Compare the candidates with themselves */

  compare_rzw_cands(props, newncand, notes);

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(props, sizeof(fourierprops), (unsigned long) newncand, \
	    candfile);
  fclose(candfile);

  /* Send the candidates to the text file */

  file_reg_candidates(props, notes, newncand, dt, \
		      (long) (N + DBLCORRECT), nph, argv[1], rzwnm);

  /* Finish up */

  printf("Done.\n\n");
  printf("Searched %ld pts (approximately %ld were independent).\n\n", \
	 totnumsearched, (long) (totnumsearched * 0.5 * dz / 6.95));

  printf("Timing summary:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n", \
	 ttim, utim, stim);
  printf("  Total time: %.3f sec\n\n", tott);

  printf("Candidates in binary format are stored in '%s'.\n", candnm);
  printf("A candidate Postscript table is stored in '%s.ps'.\n\n", rzwnm);

/*     readinf(&idata, infonm); */
/*     realpsr = comp_psr_to_cand(&props, &idata, compare); */
/*     printf("%s\n",compare); */

  fclose(fftfile);
  free(list);
  free(derivs);
  free(props);
  free(notes);
  if (idata.onoff) free(idata.onoff);
  return (0);
}


int not_already_there_rzw(position * newpos, position * list, int nlist)
{
  int ii;

  /* Loop through the candidates already in the list */

  for (ii = 0; ii < nlist; ii++) {
    if (list[ii].pow == 0.0)
      break;

    /* Do not add the candidate to the list if it is a lower power */
    /* version of an already listed candidate.                     */

    if (fabs(newpos->p1 - list[ii].p1) < 0.6) {
      if (fabs(newpos->p2 - list[ii].p2) < 2.1) {
	if (fabs(newpos->p3 - list[ii].p3) < 5.0) {
	  if (newpos->pow < list[ii].pow) {
	    return 0;
	  }
	}
      }
    }
  }
  return 1;
}


void compare_rzw_cands(fourierprops * list, int nlist, char *notes)
{
  int ii, jj, kk;
  char tmp[30];

  /* Loop through the candidates (reference cands) */

  for (ii = 0; ii < nlist; ii++) {

    /* Loop through the candidates (referenced cands) */

    for (jj = 0; jj < nlist; jj++) {
      if (ii == jj)
	continue;

      /* Look for standard sidelobes */

      if (fabs(list[ii].r - list[jj].r) < 15.0 && \
	  fabs(list[ii].z - list[jj].z) > 1.0 && \
	  list[ii].pow > list[jj].pow) {

	/* Check if the note has already been written */

	sprintf(tmp, "%.20s", notes + jj * 20);
	if (!strcmp("                    ", tmp)) {

	  /* Write the note */

	  sprintf(notes + jj * 20, "SL? of Cand %d", ii + 1);
	}
	continue;
      }
      /* Loop through the possible PSR period harmonics */

      for (kk = 1; kk < 61; kk++) {

	/* Check if the PSR Fourier freqs and z's are close enough */

	if ((fabs(list[ii].r - list[jj].r / kk) < list[jj].rerr * 3) && \
	    (fabs(list[ii].z - list[jj].z / kk) < list[jj].zerr * 2)) {

	  /* Check if the note has already been written */

	  sprintf(tmp, "%.20s", notes + jj * 20);
	  if (!strcmp("                    ", tmp)) {

	    /* Write the note */

	    sprintf(notes + jj * 20, "H %d of Cand %d", kk, ii + 1);

	    break;
	  }
	}
      }
    }
  }
}
#undef SHORTESTFFT
#undef LOSKIP
