#include "presto.h"
#include "bincand_cmd.h"
#include "orbint.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[]){
/* bincand: tries to optimize a binary PSR candidate based  */
/* on the mini-FFT method in a large FFT.  Uses an orbit    */
/* generator to determine trial modulation specta and then  */
/* correlates these with the initial FFT.                   */
  FILE *fftfile;
  char fftnm[200];
  float powr, bestpowr=0.0, binr, bini, nextbinr, nextbini;
  double N, T, rlo, rhi, rorb=0.0, norm, bestpsrp=0.0;
  double phiorb=0.0, ro=0.0, plo, phi;
  double epoch=0.0, freq=0.0, ppsr=0.0;
  double powargr, powargi, interbinfact;
  int pnum, ii, resplen, resphw;
  int pct, xct, tct, datalen, lodata;
  double dp=0.0, dx=0.0, dt=0.0, de=0.016, dw=0.8;
  double lop=0.0, lox=0.0, lot=0.0, loe, low;
  int np=0, nx=0, nt=0, ne, nw, numgood;
  unsigned long filelen;
  double numsearched=0.0, numtosearch=0.0, tott;
  struct tms runtimes;
  fcomplex *data=NULL, *resp=NULL, *corrdata=NULL;
  orbitparams trialorb, orb, bestorb;
  infodata idata;
  makedata mdata;
  psrdatabase pdata;
  psrparams psr;
  Cmdline *cmd;

  /* Call usage() if we have no command line arguments */

  if (argc <= 2) {
    usage();
    exit(1);
  }
  interbinfact = 1.0 / (PIBYTWO * PIBYTWO);

  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("     Binary Candidate Optimization Routine\n");
  printf("            by Scott M. Ransom\n");
  printf("               19 June, 2000\n\n");

  /* Initialize the filenames and some data: */

  sprintf(fftnm, "%s.fft", cmd->argv[0]);
  trialorb.p = 0.0; trialorb.x = 0.0; trialorb.e = 0.0;
  trialorb.w = 0.0; trialorb.t = 0.0; trialorb.wd = 0.0;
  trialorb.pd = 0.0;

  /* Open the fft file */

  fftfile = chkfopen(fftnm, "rb");

  /* The number of points in datafile */

  filelen = chkfilelen(fftfile, sizeof(float));

  /* Read the info file */

  printf("Reading observation information from \n\t'%s.inf'.\n\n", 
	 fftnm);
  readinf(&idata, cmd->argv[0]);

  /* The MJD of the beginning of the observation */

  if (idata.mjd_i && idata.mjd_f) {
    epoch = (double) idata.mjd_i + idata.mjd_f;
  }
  N = idata.N;
  T = idata.N * idata.dt;

  /* Calculate some of the parameters */

  printf("Total integration time  =  %g sec\n", T);
  printf("Original data file had  =  %.0f points\n", N);
  printf("Observation epoch (MJD) =  %.12f\n\n", epoch);

  /* Read the pulsar database if needed */

  if (cmd->psrnameP) {
    int numpsrs;

    numpsrs = read_database(&pdata);
    pnum = get_psr_at_epoch(cmd->psrname, epoch, &pdata, &psr);
    if (!pnum) {
      printf("The pulsar '%s' is not in the database.  Exiting.\n\n",
	     cmd->psrname);
      exit(1);
    }
    if (psr.ntype & 8){
      trialorb = psr.orb;
    } else {
      printf("The pulsar '%s' is not in a binary.  Exiting.\n\n",
	     cmd->psrname);
      exit(1);
    }
    ppsr = psr.p;
    if (cmd->ploP && cmd->phiP){
      plo = cmd->plo;
      phi = cmd->phi;
      rlo = T/phi;
      rhi = T/plo;
    } else if (cmd->rloP && cmd->rhiP){
      rlo = cmd->rlo;
      rhi = cmd->rhi;
      plo = T/rhi;
      phi = T/rlo;
    } else {
      rlo = psr.f*T - 1;
      rhi = psr.f*T + 1;
      plo = T/rhi;
      phi = T/rlo;
    }

    phiorb = TWOPI * trialorb.x / ppsr;
    rorb = T / trialorb.p;
    ro = T / ppsr;
    if (strcmp("PSR B        ", psr.bname)==0)
      printf("Using parameters from pulsar %s\n\n", psr.jname);
    else
      printf("Using parameters from pulsar %s\n\n", psr.bname);
    printf("    Pulsar period = %g sec\n", ppsr);
    printf("  Pulsar FFT freq = %.3f cycles\n", ro);
    printf("   Orbital period = %.3f sec\n", trialorb.p);
    printf(" Orbital FFT freq = %g cycles\n", rorb);
    printf("       a*sin(i)/c = %.3f lt-sec\n", trialorb.x);
    printf("     eccentricity = %.3f\n", trialorb.e);
    printf("  angle of periap = %.3f deg\n", trialorb.w);
    printf("   time of periap = %.3f sec\n", trialorb.t);
    printf(" Modulation 'Phi' = %.3f\n\n", phiorb);

    dp = trialorb.p * \
      exp(0.9792168 * log(trialorb.p/phiorb) - 10.9658871);
    dx = trialorb.x * \
      exp(0.9572412 * log(1.0/phiorb) + 0.7110553);
    dt = trialorb.p * \
      exp(0.9420009 * log(1.0/phiorb) - 1.1676730);
    np = nx = nt = nw = ne = 6;
    lop = trialorb.p - np / 2 * dp;
    lox = trialorb.x - nx / 2 * dx;
    lot = trialorb.t - nt / 2 * dt;
    loe = trialorb.e - ne / 2 * de;
    if (loe < 0.0) loe = 0.0;
    low = fmod(trialorb.w - nw / 2 * dw, 360.0) * 360.0;

    /* If the user specifies all of the binaries parameters */
    
  } else if (cmd->makefileP) {

    read_mak_file(cmd->argv[0], &mdata);
    ppsr = mdata.p;
    trialorb = mdata.orb;
    freq = mdata.f;
    if (cmd->ploP && cmd->phiP){
      plo = cmd->plo;
      phi = cmd->phi;
      rlo = T/phi;
      rhi = T/plo;
    } else if (cmd->rloP && cmd->rhiP){
      rlo = cmd->rlo;
      rhi = cmd->rhi;
      plo = T/rhi;
      phi = T/rlo;
    } else {
      rlo = mdata.r - 1;
      rhi = mdata.r + 1;
      plo = T/rhi;
      phi = T/rlo;
    }

    phiorb = TWOPI * trialorb.x / ppsr;
    rorb = T / trialorb.p;
    ro = T / ppsr;
    printf("Using parameters from\n\t'%s.mak':\n\n", cmd->argv[0]);
    printf("    Pulsar period = %g sec\n", ppsr);
    printf("  Pulsar FFT freq = %.3f cycles\n", ro);
    printf("   Orbital period = %.3f sec\n", trialorb.p);
    printf(" Orbital FFT freq = %g cycles\n", rorb);
    printf("       a*sin(i)/c = %.3f lt-sec\n", trialorb.x);
    printf("     eccentricity = %.3f\n", trialorb.e);
    printf("  angle of periap = %.3f deg\n", trialorb.w);
    printf("   time of periap = %.3f sec\n", trialorb.t);
    printf(" Modulation 'Phi' = %.3f\n\n", phiorb);

    dp = trialorb.p * \
      exp(0.9792168 * log(trialorb.p/phiorb) - 10.9658871);
    dx = trialorb.x * \
      exp(0.9572412 * log(1.0/phiorb) + 0.7110553);
    dt = trialorb.p * \
      exp(0.9420009 * log(1.0/phiorb) - 1.1676730);
    np = nx = nt = nw = ne = 6;
    lop = trialorb.p - np / 2 * dp;
    lox = trialorb.x - nx / 2 * dx;
    lot = trialorb.t - nt / 2 * dt;
    loe = trialorb.e - ne / 2 * de;
    if (loe < 0.0) loe = 0.0;
    low = fmod(trialorb.w - nw / 2 * dw, 360.0) * 360.0;

    /* Determine the candidate parameters to examine from a _bin.cand file */

  } else if (cmd->candnumP) {
    double hir, hipow;
    float *powers;
    int miniN;
    fourierprops fp;
    rderivs rd;
    rawbincand bincand;
    binaryprops binprops;

    if (!cmd->candfileP) {
      printf("\nYou must enter a name for the bin candidate ");
      printf("file (-candfile filename)\n");
      printf("Exiting.\n\n");
      exit(1);
    }

    printf("Optimizing candidate %d from \n\t'%s':\n", cmd->candnum, 
	   cmd->candfile);

    /* Get the candidate */

    get_rawbin_cand(cmd->candfile, cmd->candnum, &bincand);
    miniN = bincand.mini_N;

    if (cmd->ploP && cmd->phiP){
      plo = cmd->plo;
      phi = cmd->phi;
      rlo = T/phi;
      rhi = T/plo;
    } else if (cmd->rloP && cmd->rhiP){
      rlo = cmd->rlo;
      rhi = cmd->rhi;
      plo = T/rhi;
      phi = T/rlo;
    } else {
      rlo = bincand.full_lo_r - 2.0 * miniN;
      rhi = bincand.full_lo_r + 3.0 * miniN;
      plo = T/rhi;
      phi = T/rlo;
    }

    /* Read in the part of the FFT used in the search so that */
    /* we can calculate a binaryprops structure.              */

    data = read_fcomplex_file(fftfile, bincand.full_lo_r, miniN);
    powers = gen_fvect(miniN);
    for (ii = 0; ii < miniN; ii++)
      powers[ii] = POWER(data[ii].r, data[ii].i);
    free(data);
    realfft(powers, miniN, -1);
    norm = sqrt(miniN) / powers[0];
    for (ii = 0; ii < miniN; ii++) 
      powers[ii] *= norm;
    powers[0] = 1.0;
    hipow = max_r_arr((fcomplex *)powers, miniN/2, 
		      bincand.mini_r, &hir, &rd);
    free(powers);
    calc_props(rd, hir, 0.0, 0.0, &fp);
    calc_binprops(&fp, bincand.full_T, bincand.full_lo_r, miniN, 
		  &binprops);
    ppsr = binprops.ppsr;
    trialorb.p = binprops.pbin;
    trialorb.e = 0.0;
    trialorb.x = binprops.asinic;
    trialorb.w = 0.0;
    trialorb.t = 0.0;
    trialorb.wd = 0.0;
    print_bin_candidate(&binprops, 2);

    dp = trialorb.p * \
      exp(0.9792168 * log(trialorb.p/phiorb) - 10.9658871);
    dx = trialorb.x * \
      exp(0.9572412 * log(1.0/phiorb) + 0.7110553);
    dt = trialorb.p * \
      exp(0.9420009 * log(1.0/phiorb) - 1.1676730);
    np = 8.0 * binprops.pbinerr / dp + 1;
    nx = 8.0 * binprops.asinicerr / dx + 1;
    nt = trialorb.p / dt + 1;
    ne = 1;  /* This only works for circular orbs */
    nw = 1;  /* This only works for circular orbs */
    lop = trialorb.p - np / 2 * dp;
    lox = trialorb.x - nx / 2 * dx;
    lot = 0.0;
    loe = 0.0;  /* This only works for circular orbs */
    low = 0.0;  /* This only works for circular orbs */
  }

  /* Determine the pulsar parameters to check if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */

  if (!cmd->candnumP && !cmd->makefileP && !cmd->psrnameP) {
    
    if (!cmd->usrP || \
	(cmd->usrP && \
	 (!cmd->pbP || !cmd->asinicP || !cmd->eP || !cmd->ToP || !cmd->wP))){
      printf("\nIf you do not specify a pulsar, a binary candidate, \n");
      printf("of a makefile, you _must_ specify '-usr' and all of the\n");
      printf("candidate properties on the command line.  Exiting.\n\n");
      exit(1);
    } else {
      trialorb.p = cmd->pb;
      trialorb.e = cmd->e;
      trialorb.x = cmd->asinic;
      trialorb.w = cmd->w;
      trialorb.t = cmd->To;
      trialorb.wd = cmd->wdot;
    }
    if (cmd->ploP && cmd->phiP){
      plo = cmd->plo;
      phi = cmd->phi;
      rlo = T/phi;
      rhi = T/plo;
    } else if (cmd->rloP && cmd->rhiP){
      rlo = cmd->rlo;
      rhi = cmd->rhi;
      plo = T/rhi;
      phi = T/rlo;
    } else {
      printf("\nYou must specify limits on the pulsar period to search.\n");
      printf("You can do this using the '-rlo' and '-rhi' or  '-plo' and\n");
      printf("'-phi' flags.  Exiting.\n\n");
      exit(1);
    }
    ppsr = 0.5 * (plo + phi);
    phiorb = TWOPI * trialorb.x / ppsr;
    rorb = T / trialorb.p;
    ro = T / ppsr;
    printf("Using the following parameters:\n\n");
    printf("    Pulsar period = %g sec\n", ppsr);
    printf("  Pulsar FFT freq = %.3f cycles\n", ro);
    printf("   Orbital period = %.3f sec\n", trialorb.p);
    printf(" Orbital FFT freq = %g cycles\n", rorb);
    printf("       a*sin(i)/c = %.3f lt-sec\n", trialorb.x);
    printf("     eccentricity = %.3f\n", trialorb.e);
    printf("  angle of periap = %.3f deg\n", trialorb.w);
    printf("   time of periap = %.3f sec\n", trialorb.t);
    printf(" Modulation 'Phi' = %.3f\n\n", phiorb);
  }

  /* Extract the necessary data from the fftfile */

  datalen = next2_to_n(2.0 * rorb * phiorb);
  lodata = (int) (ro - datalen / 2);
  printf("Extracting frequencies from\n\t'%s.fft':\n", 
	 cmd->argv[0]);
  printf("\tNumber of frequencies = %d\n", datalen);
  printf("\tLow Fourier frequency = %d\n", lodata);
  data = read_fcomplex_file(fftfile, lodata, datalen);
  fclose(fftfile);

  /* Determine and normalize the Fourier amplitudes */
  /* using the local power level, then FFT.         */

  norm = 0.0;
  for (ii = 0; ii < datalen; ii++)
    norm += POWER(data[ii].r, data[ii].i);
  norm = 1.0 / sqrt(norm / datalen);
  for (ii = 0; ii < datalen; ii++){
    data[ii].r *= norm;
    data[ii].i *= norm;
  }
  COMPLEXFFT(data, datalen, -1);

  /* Begin the loop over the possible orbital params */
  /* So far only circular orbits are implemented...  */

  numtosearch = np * nx * nt;
  corrdata = gen_cvect(datalen);
  printf("\nWill search:\n");
  printf("\tOrbital periods = %d\n", np);
  printf("\tSemi-major axes = %d\n", nx);
  printf("\tPeriapsis times = %d\n", nt);
  printf("  For a total of %.0f candidates.\n\n", numtosearch);

  for (pct = 0; pct < np; pct++){
    orb.p = lop + pct * dp;
    for (xct = 0; xct < nx; xct++){
      orb.x = lox + xct * dx;
      for (tct = 0; tct < nt; tct++){
	orb.t = lot + tct * dt;
	printf("p = %f, x = %f, t = %f\n", orb.p, orb.x, orb.t);
	tott = times(&runtimes) / (double) CLK_TCK;

	/* Generate the response */

	resphw = bin_resp_halfwidth(ppsr, T, &orb);
	resplen = 2 * resphw;
	resp = gen_bin_response(0.0, 1, ppsr, T, &orb, resplen);

	/* Perform the correlation */

	numgood = corr_complex(data, datalen, FFT, \
			       resp, resplen, RAW, \
			       corrdata, datalen, resphw, \
			       1, resphw, CORR);

	/* Search the correlation */

	binr = corrdata[resphw].r;
	bini = corrdata[resphw].i;
	for (ii = resphw; ii < resphw + numgood; ii++){
	  /* Check the integer bin */
	  powr = POWER(binr, bini);
	  if (powr > bestpowr){
	    bestpowr = powr;
	    bestpsrp = lodata + ii;
	    bestorb = orb;
	  }
	  /* Check the interbin */
	  nextbinr = corrdata[ii+1].r;
	  nextbini = corrdata[ii+1].i;
	  binr -= nextbinr;
	  bini -= nextbini;
	  powr = interbinfact * POWER(binr, bini);
	  if (powr > bestpowr){
	    bestpowr = powr;
	    bestpsrp = lodata + ii + 0.5;
	    bestorb = orb;
	  }
	  binr = nextbinr;
	  binr = nextbini;
	}
	free(resp);
	tott = times(&runtimes) / (double) CLK_TCK - tott;
	numsearched = numsearched + 1.0;
	printf("\tt = %fs per cand.  Estimate %.3f hours for full search.\n", 
	       tott, (numtosearch - numsearched) * tott / 3600.0);
	printf("\tBestpow = %f  BestPsrP = %f  BestOrbP = %f\n", 
	       bestpowr, bestpsrp, bestorb.p);
      }
    }
  }
  free(data);
  exit(0);
}
