#include "presto.h"
#include "bincand_cmd.h"
#include "orbint.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char **argv)
/* bincand: tries to optimize a binary PSR candidate based  */
/* on the mini-FFT method in a large FFT.  Uses an orbit    */
/* generator to determine trial modulation specta and then  */
/* correlates these with the initial FFT.                  */
{
  FILE *fftfile, *candfile;
  char fftnm[200], candfilenm[200];
  char tmp1[100], tmp2[100], psrjname[30], psrbname[30], pname[30];
  float *resp=NULL, *data=NULL, *freqs=NULL, norm, nph;
  double N, T, dt, difft = 0.0, rlo, rhi, p_psr, pdot_psr;
  double rint, z, roffset, qlo, dq, ro, plo, phi;
  double epoch = 0.0, freq = 0.0;
  int nq, tol=5, pnum;
  long i, numbetween = 3, m = 0, numout = 200;
  unsigned long filelen;
  long pct, ect, tct, wct, xct;
  struct tms runtimes;
  double ttim, stim, utim, tott;
  orbitparams trialorb, bestorbs[50];
  binaryprops bincand;
  infodata idata;
  makedata mdata;
  psrdatabase pdata;
  Cmdline *cmd;

  /* Call usage() if we have no command line arguments */

  if (argc <= 2) {
    Program = argv[0];
    usage();
    exit(1);
  }
  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("     Binary Candidate Optimization Routine\n");
  printf("            by Scott M. Ransom\n");
  printf("               29 June, 1998\n\n");

  /* Initialize the filenames and some data: */

  sprintf(fftnm, "%s.fft", cmd->argv[0]);
  sprintf(candfilenm, "%s_bin.cand", cmd->argv[0]);
  trialorb.p = 0.0; trialorb.x = 0.0; trialorb.e = 0.0;
  trialorb.w = 0.0; trialorb.t = 0.0; trialorb.wd = 0.0;
  trialorb.pd = 0.0;

  /* Open the fft file */

  fftfile = chkfopen(fftnm, "rb");

  /* The number of points in datafile */

  filelen = chkfilelen(FILE *file, sizeof(float));

  /* Read the info file */

  readinf(&idata, cmd->argv[0]);
  if (idata.object) {
    printf("Optimizing a %s candidate from '%s'.\n", idata.object, datanm);
  } else {
    printf("Optimizing a candidate from '%s'.\n", datanm);
  }

  /* The MJD of the beginning of the observation */

  if (idata.mjd_i && idata.mjd_f) {
    epoch = (double) idata.mjd_i + idata.mjd_f;
  }
  dt = idata.dt;
  N = idata.N;
  T = N * dt;

  /* Read the pulsar database if needed */

  if (cmd->psrnameP) {

    np = read_database(&pdata);
    pnum = get_psr_at_epoch(cmd->psrname, epoch, &pdata, &psr);
    if (!pnum) {
      printf("The pulsar is not in the database.  Exiting.\n\n");	
      exit(1);							
    }
    if (psr.ntype & 8){
      trialorb = psr.orb;
    } else {
      printf("The pulsar is not in a binary.  Exiting.\n\n");	
      exit(1);							
    } 
    p_psr = psr.p;
    pdot_psr = psr.pd;
    strcpy(pname, psr.jname);
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

    /* If the user specifies all of the binaries parameters */
    
  } else if (cmd->makefileP) {

    read_mak_file(cmd->argv[0], &mdata);
    N = mdata.N;
    dt = mdata.dt;
    p_psr = mdata.p;
    pdot_psr = mdata.pd;
    trialorb = mdata.orb;
    freq = mdata.f;
    dfdt = mdata.fd;
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

    /* Determine the candidate parameters to examine from a _bin.cand file */

  } else if (cmd->candnumP) {

    /*  For possible use later.......
    if (cmd->binfileP) {
      printf("\nYou must enter a name for the bin candidate ");
      printf("file (-binfile filename)\n");
      printf("Exiting.\n\n");
      exit(1);
    }
    */

    /* Open the fft file */

    get_bin_cand(candfilenm, cmd->candnum, &bincand);

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
      rlo = bincand.rpsr - 10.0 * bincand.nfftbins;
      rhi = bincand.rpsr + 10.0 * bincand.nfftbins;
      plo = T/rhi;
      phi = T/rlo;
    }

  }

  /* Determine the pulsar parameters to fold if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */

  if (!cmd->candnumP && !cmd->makefileP && !cmd->psrnameP) {
    
    if (!cmd->usrP || \
	(cmd->usrP && \
	 (!cmd->pbP || !cmd->xP || !cmd->eP || !cmd->ToP || !cmd->wP)){
      printf("\nIf you do not specify a pulsar, a binary candidate, \n");
      printf("of a makefile, you _must_ specify '-usr' and all of the\n");
      printf("candidate properties on the command line.  Exiting.\n\n");
      exit(1);
    } else {
      trialorb.p = cmd->pb;
      trialorb.e = cmd->e;
      trialorb.x = cmd->x;
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
      rlo = bincand.rpsr - 10.0 * bincand.nfftbins;
      rhi = bincand.rpsr + 10.0 * bincand.nfftbins;
      plo = T/rhi;
      phi = T/rlo;
    }
    }
    
 /* Calculate some of the parameters */

  z = TWOPI * orb.x / ppsr;
  ro = T / ppsr;
  roffset = modf(ro, &rint);

  printf("Integration time  = %f sec\n",T);
  printf("\nPulsar period     = %f sec\n", ppsr);
  printf("Pulsar FFT freq   = %f cycles\n", ro);
  printf("Pulsar amplitude  = %f \n", apsr);
  printf("Response roffset  = %f cycles\n", roffset);
  printf("\nBinary period     = %f sec\n", orbit.p);
  printf("Binary FFT freq   = %f cycles\n", T / orbit.p);
  printf("Binary a*sin(i)/c = %f lt-sec\n", orbit.x);
  printf("Binary e          = %f\n", orbit.e);
  printf("Binary w          = %f deg\n", orbit.w);
  printf("Binary To         = %f sec\n", orbit.t);
  printf("\nModulation z      = %f\n", z);


  /* Generate the theoretical response */

  printf("\n\nGenerating the theoretical response...\n");

  tott = times(&runtimes) / (double) CLK_TCK;

  roffset = 0.0;

  resp = gen_bin_response(roffset, numbetween, ppsr, T, \
			  &orbit, tol, &m);
/*   resp = gen_r_response(roffset, numbetween, LOWACC, &m);  */

  qlo = m / (double) numbetween + roffset;
  dq = -1.0 / (double) numbetween;
  nq = 2 * m;
  printf("m = %ld\n",m);
  
  printf("Timing summary for theoretical response creation:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	 ttim, utim, stim);
  printf("Total time elapsed:  %.3f sec.\n\n", tott);
