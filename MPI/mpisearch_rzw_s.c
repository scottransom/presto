#include "fftapps.h"
#include "mpi.h"
#ifdef LOGGING
#include "mpe.h"
#endif

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

int main(int argc, char *argv[])
{
  FILE *fftfile, *candfile;
  double dt, nph, T, N, bigz, hir, hiz, dz = 2.0;
  float locpow = 1.0, powavg, powadev, powsdev, powvar, powskew, powcurt;
  float **kernels, *data = &locpow, *powlist, *response, *ptr;
  float powdiff, hipowchop, lowpowlim;
  float chkpow = 0.0, hipow = 0.0, minpow = 0.0, numr = 0.0;
  long numbetween = 2, startbin, lofreq = 0.0, lobinskip, hibinskip;
  long nreal, nbins, nr = 1, nz, corrsize = 0, mincorrsize, tmp = 0;
  long rlo, rhi, zlo = 0, zhi = 0;
  long i, ct, zct, ndupes, m, nextbin;
  unsigned long totnumsearched = 0;
  int ncand, newncand, oldper = 0, newper = 0;
  char filenm[100], candnm[100], *notes;
  char rzwnm[100];
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

  /* Initialize the filenames: */

  sprintf(filenm, "%s.fft", argv[1]);
  sprintf(candnm, "%s_rzw.cand", argv[1]);
  sprintf(rzwnm, "%s_rzw", argv[1]);

  /* Read the info file */

  idata = readinf(argv[1]);
  if (idata.object) {
    printf("Analyzing %s data from '%s'.\n\n", idata.object, filenm);
  } else {
    printf("Analyzing data from '%s'.\n\n", filenm);
  }

  /* open the FFT file and get its length */

  fftfile = chkfopen(filenm, "rb");
  nph = get_numphotons(fftfile);
  chkfseek(fftfile, 0, SEEK_END);

  /* # of real data points    */

  nreal = ftell(fftfile) / sizeof(float);

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
  bigz = DMAX(fabs(zlo), fabs(zhi));
  m = get_z_resp_m(numbetween, bigz, LOWACC);
  nz = (int) ((zhi - zlo) / dz) + 1;

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

  mincorrsize = 6 * m;

  /* We need to keep nz + 2 arrays of corrsize complex numbers in */
  /* memory at all times. (nz kernels, 1 data set, 1 result)      */
  /* Therefore tmp is the largest corrsize we can have.           */

  tmp = (int) ((MAXREALFFT >> 1) / (nz + 2));

  /* Determine corrsize */

  corrsize = 1;
  while (corrsize < tmp) {
    corrsize <<= 1;
  }

  /* Insure smaller than tmp then *0.5 to be sure... */

  corrsize >>= 2;

  if (mincorrsize > corrsize) {
    printf("\nYou are asking for too much memory.  Specify\n");
    printf("fewer z values to search, or lower z values.\n\n");
    exit(1);
  }
  /* Generate the correlation kernels */

  printf("Generating fdot kernels for the correlations...\n");

  kernels = gen_cmatrix(nz, corrsize);
  for (i = 0; i < nz; i++) {
    response = gen_z_response(0.0, numbetween, zlo + i * dz, \
			      LOWACC, &m);
    place_complex_kernel(response, m, corrsize, kernels[i]);
    free(response);
    COMPLEXFFT(kernels[i], corrsize, -1);
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
  notes = malloc(sizeof(char) * ncand * 20);

  dt = idata.bindt;
  N = idata.numbins;
  T = N * dt;
  numr = (rhi - rlo + 1) * nz * 2;

  /* We will automatically get rid of any candidates that have local */
  /* powers much lower than what we would expect to find in the      */
  /* search from purely statistical reasons alone                    */
  /* Note:  6.95 is the full width in z of a signal response         */

  lowpowlim = -0.8 * log(ncand / (numr * 0.5 * dz / 6.95));

  /* Initialize the candidate list */

  for (i = 0; i < ncand; i++) {
    list[i].pow = 0.0;
    list[i].p1 = 0.0;
    list[i].p2 = 0.0;
    list[i].p3 = 0.0;
  }

  /* Start the main search loop */

  nextbin = rlo;

  do {

    startbin = (unsigned long) nextbin;

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
      data = corr_complex_raw_file(fftfile, numbetween, startbin, \
				   kernels[zct], corrsize, m, &nextbin);
      nr = (int) (nextbin - startbin) * numbetween;

      if (zct == 0) {

	/* Get approximate local power statistics */

	powlist = gen_fvect(nr);
	ptr = data;
	tmp = (nextbin > rhi) ? (rhi - startbin) * numbetween : nr;
	for (i = 0; i < tmp; i++) {
	  powlist[i] = power(*(ptr++), *(ptr++));
	}
	moment(powlist - 1, tmp, &powavg, &powadev, &powsdev, \
	       &powvar, &powskew, &powcurt);

	/* Throw powers away that are greater than 4.0 sdev above mean  */
	/* Also throw powers away that are very small.  This could show */
	/* that we have come into a zero-padded area.                   */

	ct = tmp;
	locpow = powavg * tmp;
	hipowchop = 4.0 * powsdev + powavg;
	for (i = 0; i < tmp; i++) {
	  powdiff = powlist[i] - powavg;
	  if (powdiff > hipowchop) {
	    locpow -= powlist[i];
	    ct--;
	  }
	}
	locpow /= (double) ct;

	free(powlist);
      }
      /* This loop is the heart of the search */

      ptr = data;
      for (i = 0; i < tmp; i++) {
	chkpow = power(*(ptr++), *(ptr++)) / locpow;
	if (chkpow > minpow) {
	  newpos.pow = chkpow;
	  newpos.p1 = startbin + 0.5 * i;
	  if (newpos.p1 > rhi)
	    break;
	  newpos.p2 = zlo + zct * dz;
	  newpos.p3 = 0.0;
	  list[ncand - 1] = newpos;
	  minpow = percolate(list, ncand);
	}
      }

      totnumsearched += tmp;

      /* Free the memory containing the correlation result */

      free(data);

    }

  } while (nextbin <= rhi);

  /* Free the memory used by the correltation kernels */

  free(kernels[0]);
  free(kernels);

  printf("\rAmount of search complete = %3d%%", 100);
  fflush(stdout);
  printf("\nDone searching.  Now optimizing each candidate and sorting.\n\n");

  /* Do rough duplicate removal */

  ndupes = remove_dupes(list, ncand);
  newncand = ncand - ndupes;

  /* Now maximize each good candidate */

  newper = 0;
  oldper = 0;

  for (i = 0; i < newncand; i++) {

    /* Calculate percentage complete */

    newper = (int) (i / (float) (newncand) * 100.0) + 1;
    if (newper > oldper) {
      printf("\rAmount of optimization complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
    }
    hipow = max_rz_file(fftfile, list[i].p1, list[i].p2, \
			&hir, &hiz, &derivs[i]);
    props[i] = calc_props(derivs[i], hir + lofreq, hiz, 0.0);
  }
  printf("\rAmount of optimization complete = %3d%%\n\n", 100);

  qsort(props, newncand, sizeof(fourierprops), compare_fourierprops);

  /* Do fine scale duplicate removal and other cleaning */

  ndupes = remove_dupes2(props, newncand);
  newncand -= ndupes;
  ndupes = remove_other(props, newncand, rlo, rhi, lowpowlim);
  newncand -= ndupes;

  /* Need to check the database, and remove harmonics */

  for (i = 0; i < newncand; i++) {
    comp_psr_to_cand(&props[i], &idata, notes + i * 20, 0);
  }

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(props, sizeof(fourierprops), newncand, candfile);
  fclose(candfile);

  /* Send the candidates to the text file */

  file_reg_candidates(props, notes, newncand, dt, N, nph, argv[1]);

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
  printf("A candidate Postscript table is stored in '%s_rzw.ps'.\n\n", \
	 argv[1]);

/*     idata = readinf(infonm); */
/*     realpsr = comp_psr_to_cand(&props, &idata, compare); */
/*     printf("%s\n",compare); */

  fclose(fftfile);
  free(list);
  free(derivs);
  free(props);
  free(notes);
  return (0);
}

#undef SHORTESTFFT
#undef LOSKIP
