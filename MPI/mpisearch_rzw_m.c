#include "fftapps.h"
#include "mpi.h"
#ifdef LOGGING
#include "mpe.h"
#endif

#define SHORTESTFFT 32768	/* The shortest data set we will look at */
#define LOSKIP      200		/* The min # of low freq bins to skip    */

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

typedef struct FFTINFO {
    double N;       /* Number of points in original time series      */
    double dt;      /* Sample duration per point of time series      */
    long nbins;     /* Number of points in the FFT file              */
    float nph;      /* Frequency 0 bin amplitude of the FFT          */
} fftinfo;

/* Some global datatypes */

extern MPI_Datatype fftinfotype;
extern MPI_Datatype positiontype;
extern MPI_Datatype fourierpropstype;

/* Function definitions */

extern void make_position_struct(void);
extern void make_fftinfo_struct(void);
extern void make_fourierpropstype_struct(void);

void master(int numprocs, int argc, char *argv[])
{
  int myid=0;
  FILE *fftfile, *candfile;
  double dt, nph, T, N, bigz, hir, hiz, dz = 2.0;
  float locpow = 1.0, powavg, powadev, powsdev, powvar, powskew, powcurt;
  float *data = &locpow, *powlist, *ptr;
  float powdiff, hipowchop, lowpowlim;
  float chkpow = 0.0, hipow = 0.0, minpow = 0.0, numr = 0.0;
  long numbetween = 2, startbin, lofreq = 0.0, lobinskip, hibinskip;
  long numbins, beginbin, endbin, numtmplist;
  long nreal, nbins, nr = 1, nz, corrsize = 0, mincorrsize, tmp = 0;
  long rlo, rhi, zlo = 0, zhi = 0;
  long i, ct, zct, ndupes, m, nextbin;
  unsigned long totnumsearched = 0;
  int ncand, newncand, oldper = 0, newper = 0;
  char filenm[100], candnm[100], rzwnm[100], *notes;
  position *list, **tmplist, newpos;
  fourierprops *props;
  infodata idata;
  fftinfo fdata;
  MPI_Request req[];
  MPI_Status status;


#ifdef LOGGING
  MPE_Describe_state(1,   2, "Broadcast",     "red:vlines3");
  MPE_Describe_state(3,   4, "Compute",       "blue:gray3");
  MPE_Describe_state(5,   6, "Correlate",     "green:light_gray");
  MPE_Describe_state(7,   8, "Search",        "yellow:gray");
  MPE_Describe_state(9,  10, "Collect Cands", "coral:gray");
  MPE_Describe_state(11, 12, "Send Cands",    "magenta:gray");
  MPE_Describe_state(13, 14, "Send Opt",      "forestgreen:gray");
  MPE_Describe_state(15, 16, "Optimizing",    "pink:gray");
#endif  
  
  printf("\n\n");
  printf("       Pulsar Acceleration Search Routine\n");
  printf("    MPI based.  Fully parallel.  Wicked fast.\n"); 
  printf("              by Scott M. Ransom\n");
  printf("                  1 Dec, 1997\n\n");

  /* Define a few structures we will need */

  make_fftinfo_struct();
  make_position_struct();
  make_fourierpropstype_struct();

  /* Initialize the fft filename: */

  sprintf(filenm, "%s.fft", argv[1]);

  /* Read the info file */

  idata = readinf(argv[1]);
  if (idata.object) {
    printf("Analyzing %s data from '%s' using %d processors.\n\n", \
	   idata.object, filenm, numprocs);
  } else {
    printf("Analyzing data from '%s' using %d processors.\n\n", \
	   filenm, numprocs);
  }
  dt = idata.bindt;
  N = idata.numbins;
  T = N * dt;

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
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  rewind(fftfile);

  /* Define and broadcast the essential FFT info */

  fdata.N = N;
  fdata.dt = dt;
  fdata.nbins = nbins;
  fdata.nph = nph;
  MPI_Bcast(&fdata, 1, fftinfotype, 0, MPI_COMM_WORLD);

  /* Determine the correlation sizes we will use: */

  /* The following ensures that we never waste more          */
  /* than 5% of any correlation's points due to end effects. */
  /* (Note: We throw away 2*m points from a correlation)     */
  /* We need to keep nz + 3 arrays of corrsize complex numbers in */
  /* memory at all times. (nz kernels, 1 data set, 1 result, and  */
  /* 1 array of the constants needed to conduct an FFT)           */

  if (argc > 3) {
    zlo = atoi(argv[3]);
  }
  corrsize = SHORTESTFFT >> 2;
  do {
    corrsize <<= 1;
    nz = Ramsize * 0.8 / (2 * sizeof(float) * corrsize) - 3;
    nz *= numprocs;
    if (argc == 3) {
      zlo = - nz;
      if (zlo & 1) zlo++;
    }
    zhi = zlo + 2 * (nz - 1);
    bigz = DMAX(fabs(zlo), fabs(zhi));
    m = get_z_resp_m(numbetween, bigz, LOWACC);
  } while (corrsize < 40 * m); 

  /* Initialize the output filenames: */

  sprintf(candnm, "%s_rzw_z%ld-%ld.cand", argv[1], zlo, zhi);
  sprintf(rzwnm, "%s_rzw_z%ld-%ld", argv[1], zlo, zhi);

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
  if ((ncand < 1) || (ncand > 150000)) {
    printf("\n'ncand' must be >= 1 and <= 100000.\n\n");
    exit(1);
  }

  /* The lowest freq present in the FFT file */

  if (argc >= 5) {
    lofreq = atoi(argv[4]);
    if ((lofreq < 0) || (lofreq > nbins - 1)) {
      printf("\n'lofreq' is out of range.\n\n");
      exit(1);
    }
    if (lofreq != 0)
      nph = 1.0;
  }

  /* rlo and rhi */

  if (argc >= 6) {
    rlo = atoi(argv[5]);
    if ((rlo < lofreq) || (rlo > nbins - 1)) {
      printf("\nLow frequency to search (rlo) is out of range.\n\n");
      exit(1);
    }
  }
  if (argc == 7) {
    rhi = atoi(argv[6]);
    if ((rhi < rlo) || (rhi > nbins - 1)) {
      printf("\nHigh frequency to search (rhi) is out of range.\n\n");
      printf("\nFrequencies to search are out of range.\n\n");
      exit(1);
    }
  }

  numr = (rhi - rlo + 1) * nz * numbetween;

  /* We will automatically get rid of any candidates that have local */
  /* powers much lower than what we would expect to find in the      */
  /* search from purely statistical reasons alone                    */
  /* Note:  6.95 is the full width in z of a signal response         */

  lowpowlim = -0.8 * log(ncand / ((numr / numbetween) * (dz / 6.95)));

  /* Allocate some memory */

  list = (position *)malloc(sizeof(position) * ncand);
  props = (fourierprops *)malloc(sizeof(fourierprops) * ncand);
  notes = (char *)malloc(sizeof(char) * ncand * 20);
  req = (MPI_Request *)malloc(sizeof(MPI_Request) * numprocs);
  numbins = corrsize / numbetween;
  data = gen_cvect(numbins);
  powlist = gen_fvect(numbins);

  /* Allocate the receive buffers for the candidate positions */

  numtmplist = 10 * ncand * (corrsize - 2 * m) / \
    ((rhi - rlo + 1) * numbetween * numprocs);
  tmplist = (position **)malloc(sizeof(position *) * (numprocs - 1));
  tmplist[0] = (position *)malloc(sizeof(position) * (numprocs - 1) * \
				  numtmplist);
  for (i = 1; i < numprocs - 1; i++)
    m[i] = m[i - 1] + numtmplist;

  /* Initialize the candidate list */

  for (i = 0; i < ncand; i++) {
    list[i].pow = 0.0;
    list[i].p1 = 0.0;
    list[i].p2 = 0.0;
    list[i].p3 = 0.0;
  }

  printf("Done prepping the search.\n\n");

  /* Start the main search loop */

  nextbin = rlo;

  do {

    startbin = (unsigned long) nextbin;

    /* Calculate percentage complete */
    
    newper = (int) (totnumsearched / (numr * 0.01)) + 1;
    
    if (newper > oldper) {
      newper = (newper > 99) ? 100 : newper;
      printf("\rAmount of search complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
    }

    /* Find which FFT bins we need  */
    
    beginbin = startbin - m / numbetween;
    endbin = beginbin + numbins - 1;
    
    /* Read the file to get the data */
    
    readfftarray(fftfile, beginbin, endbin, data);
    nextbin = startbin + (corrsize - 2 * m) / numbetween;
    
    /* Get approximate local power statistics */
    
    ptr = data;
    tmp = (nextbin > rhi) ? (rhi - startbin) * numbetween : numbins;
    for (i = 0; i < tmp; i++) {
      powlist[i] = power(*(ptr++), *(ptr++));
    }
    moment(powlist - 1, tmp, &powavg, &powadev, &powsdev, \
	   &powvar, &powskew, &powcurt);
    
    /* Throw powers away that are greater than 6.0 sdev above mean  */
    /* Also throw powers away that are very small.  This could show */
    /* that we have come into a zero-padded area.                   */
    
    ct = tmp;
    locpow = powavg * tmp;
    hipowchop = 6.0 * powsdev + powavg;
    for (i = 0; i < tmp; i++) {
      powdiff = powlist[i] - powavg;
      if (powdiff > hipowchop) {
	locpow -= powlist[i];
	ct--;
      }
    }
    locpow /= (double) ct;

    /* Broadcast the local power and the dataset */

    MPI_Bcast(&locpow, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data, 2*numbins, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /* Prepare to receive the lists of candidates */

    /* Initialize the temporary candidate list */
    
    for (i = 0 ; i < numtmplist ; i++) {
      tmplist[i].pow = 0.0;
      tmplist[i].p1 = 0.0;
      tmplist[i].p2 = 0.0;
      tmplist[i].p3 = 0.0;
    }
    
    /* Post the non-blocking receives for the candidates */

    for (i = 1 ; i < numprocs ; i++){
      req[i] = MPI_REQUEST_NULL;
      MPI_Irecv(&tmplist[i], numtmplist, positiontype, i, \
		MPI_ANY_TAG, MPI_COMM_WORLD, &req[i]);
    }

    /* Complete the receives */

    do{
      MPI_Waitany(numprocs, req, &j, &status);

      if (j != MPI_UNDEFINED){

	/* Process the candiates */

      }

    }while(j != MPI_UNDEFINED){

    totnumsearched += tmp * nz;
    
  } while (nextbin <= rhi);

  /* Free some memory */

  free(data);
  free(powlist);

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

  printf("Candidates in binary format are stored in '%s'.\n", candnm);
  printf("A candidate Postscript table is stored in '%s_rzw.ps'.\n\n", \
	 argv[1]);

/*     idata = readinf(infonm); */
/*     realpsr = comp_psr_to_cand(&props, &idata, compare); */
/*     printf("%s\n",compare); */

  fclose(fftfile);
  free(list);
  free(props);
  free(notes);
}

#undef SHORTESTFFT
#undef LOSKIP
