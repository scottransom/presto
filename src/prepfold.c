#include "presto.h"
#include "plot2d.h"
#include "prepfold_cmd.h"
#include "multibeam.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* This causes the barycentric motion to be calculated once per second */

#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Some function definitions */

int (*readrec_ptr)(FILE * infile, float *data, int numpts,
		   double *dispdelays, int numsubbands, int numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, int numpts,
		double *dispdelays, int numsubbands, int numchan);
void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);
int dgels_(char *trans, int *mm, int *nn, int *nrhs, 
	   double *aa, int *lda, double *bb, int *ldb, 
	   double *work, int *lwork, int *info);

double calc_phase(double t, double f, double fd, double fdd)
/* Calculate the instantaneout phase at time t given f, fd and fdd */
{
  return t * (f + t * (0.5 * fd  + fdd / 6.0 * t));
}

int bary2topo(double *topotimes, double *barytimes, int numtimes, 
	      double fb, double fbd, double fbdd, 
	      double *ft, double *ftd, double *ftdd)
/* Convert a set of barycentric pulsar spin parameters (fb, fbd, fbdd) */
/* into topocentric spin parameters (ft, ftd, ftdd) by performing      */
/* a linear least-squares fit (using LAPACK routine DGELS).  The       */
/* routine equates the pulse phase using topcentric parameters and     */
/* times to the pulse phase using barycentric parameters and times.    */
{
  double *work, *aa, *bb, dtmp;
  int ii, mm=3, nn, nrhs=1, lwork, info, index;
  char trans='T';

  if (numtimes < 4){
    printf("\n'numtimes' < 4 in bary2topo():  Cannot solve.\n\n");
    exit(0);
  }
  nn = numtimes; 
  lwork = mm + nn * 9;
  aa = gen_dvect(mm * nn);
  bb = gen_dvect(nn);
  work = gen_dvect(lwork);
  for (ii = 0; ii < nn; ii++){
    index = ii * 3;
    dtmp = (topotimes[ii] - topotimes[0]) * SECPERDAY;
    aa[index] = dtmp;
    aa[index+1] = dtmp * dtmp;
    aa[index+2] = dtmp * dtmp * dtmp;
    dtmp = (barytimes[ii] - barytimes[0]) * SECPERDAY;
    bb[ii] = dtmp * (fb + dtmp * (0.5 * fbd + fbdd * dtmp / 6.0));
  }
  dgels_(&trans, &mm, &nn, &nrhs, aa, &mm, bb, &nn, work, &lwork, &info);
  *ft = bb[0];
  *ftd = bb[1];
  *ftdd = bb[2];
  free(aa);
  free(bb);
  free(work);
  return info;
}


void quick_plot(double *data, int numdata)
{
  int ii;
  double *phases;

  phases = gen_dvect(numdata);
  for (ii=0; ii<numdata; ii++)
    phases[ii] = ii / (double) numdata;
  cpgstart_x("landscape");
  dxyline(numdata, phases, data, "Profile Phase", "Intensity", 1);
  cpgend();
  free(phases);
}

/* The main program */

int main(int argc, char *argv[])
{
  FILE *infile=NULL, *filemarker;
  float *data=NULL;
  double p=0.0, pd=0.0, pdd=0.0, f=0.0, fd=0.0, fdd=0.0;
  double bestdm=0.0, bestp=0.0, bestpd=0.0;
  double difft, tt, recdt=0.0, *dispdts=NULL;
  double *periods, *pdots, *parttimes, *dms=NULL, *bestprof;
  double orb_baryepoch=0.0, topoepoch=0.0, baryepoch=0.0, barydispdt;
  double dtmp, *Ep=NULL, *tp=NULL, startE=0.0, orbdt=1.0;
  double N=0.0, dt=0.0, T=0.0, endtime=0.0, dtdays, avgvoverc=0.0;
  double *voverc=NULL, *obsf=NULL, foldf=0.0, foldfd=0.0, foldfdd=0.0;
  double *profs=NULL, *barytimes=NULL, *topotimes=NULL, proftime;
  char obs[3], ephem[10], *outfilenm, *rootfilenm;
  char pname[30], rastring[50], decstring[50], *cptr;
  int numchan=1, binary=0, np, pnum, numdelays=0;
  int info, slen, ptsperrec=1, flags=1;
  long ii, jj, kk, numbarypts=0, worklen=0, numread=0, reads_per_part=0;
  long totnumfolded=0, lorec=0, hirec=0;
  long numbinpoints=0, proflen, currentrec=0;
  unsigned long numrec=0, arrayoffset=0;
  multibeam_tapehdr hdr;
  fourierprops rzwcand;
  orbitparams orb;
  psrparams psr;
  psrdatabase pdata;
  infodata idata, rzwidata;
  Cmdline *cmd;
  foldstats *stats, beststats;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
    Program = argv[0];
    usage();
    exit(0);
  }
  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("        Pulsar Raw-Data Folding Search Routine\n");
  printf(" Used for DM and/or period determination of PSR candidates.\n");
  printf("                 by Scott M. Ransom\n");
  printf("                    23 Dec, 1999\n\n");

  /* Determine the root input file name and the input info file name */

  cptr = strrchr(cmd->argv[0], '.');
  if (cptr==NULL){
    printf("\nThe input filename (%s) must have a suffix!\n\n", 
	   cmd->argv[0]);
    exit(1);
  }
  slen = cptr - cmd->argv[0];
  rootfilenm = (char *)malloc(slen+1);
  rootfilenm[slen] = '\0';
  strncpy(rootfilenm, cmd->argv[0], slen);
  if (!cmd->pkmbP && !cmd->ebppP){
    printf("Reading input data from '%s'\n", cmd->argv[0]);
    printf("Reading information from '%s.inf'\n\n", rootfilenm);

    /* Read the info file if available */

    readinf(&idata, rootfilenm);
  } else if (cmd->pkmbP){
    printf("Reading Parkes Multibeam data from '%s'\n\n", 
	   cmd->argv[0]);
  } else if (cmd->pkmbP){
    printf("Reading New Effelsberg PSR data from '%s'\n\n", 
	   cmd->argv[0]);
  }

  /* Open the raw data file */

  infile = chkfopen(cmd->argv[0], "rb");

  /* Determine the output file name */

  slen = strlen(rootfilenm)+5;
  outfilenm = (char *)malloc(slen);
  outfilenm[slen-1] = '\0';
  sprintf(outfilenm, "%s.msv", rootfilenm);

  /* What ephemeris will we use?  (Default is DE200) */
  
  if (cmd->de405P)
    strcpy(ephem, "DE405");
  else
    strcpy(ephem, "DE200");
    
  /* Set-up values if we are using the Parkes Multibeam System */

  if (cmd->pkmbP){

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    rewind(infile);
    multibeam_hdr_to_inf(&hdr, &idata);
    idata.dm = cmd->dm;
    if (idata.object) {						
      printf("Folding a %s candidate from '%s'.\n", \
	     idata.object, cmd->argv[0]);
    } else {
      printf("Folding a candidate from '%s'.\n", cmd->argv[0]);
    }

    /* Some information about the size of the records */

    numrec = chkfilelen(infile, RECLEN);
    ptsperrec = DATLEN * 8 / idata.num_chan;

    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* Define the RA and DEC of the observation */
  
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Define some variables */

    dt = idata.dt;
    dtdays = idata.dt / SECPERDAY;
    recdt = dt * ptsperrec;

    /* Determine the number of records to use from the command line */

    lorec = (long) (cmd->startT * numrec + DBLCORRECT);
    hirec = (long) (cmd->endT * numrec + DBLCORRECT);
    numrec = hirec - lorec;

    /* The number of reads from the file we need for */
    /* each sub-integration.                         */ 

    reads_per_part = numrec / cmd->npart;

    /* Correct numrec so that each part will contain */
    /* the same number of records.                   */

    numrec = reads_per_part * cmd->npart;
    T = numrec * recdt;
    N = numrec * ptsperrec;
    endtime = T + 2 * TDT;

    /* Topocentric and barycentric times of folding epoch data */

    if (idata.mjd_i && idata.mjd_f) {
      topoepoch = (double) idata.mjd_i + idata.mjd_f + 
	lorec * recdt / SECPERDAY;
      barycenter(&topoepoch, &baryepoch, &dtmp, 1, rastring,
		 decstring, obs, ephem);

      /* Correct the barycentric time for the dispersion delay.     */
      /* This converts the barycentric time to infinite frequency.  */

      barydispdt = delay_from_dm(cmd->dm, idata.freq + 
				 (idata.num_chan - 1) * idata.chan_wid);
      baryepoch -= (barydispdt / SECPERDAY);
    }

    /* The data collection routine to use */

    readrec_ptr = read_multibeam_subbands;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    worklen = ptsperrec;

    /* How many sub-bands will we de-disperse to     */
    /* This should be made more robust by selecting  */
    /* a sub-bandwidth that keeps the smearing below */
    /* an acceptable level.                         */
    
    if (!cmd->nsubP)
      cmd->nsub = numchan / 8;
    else if (cmd->nsub > numchan)
      cmd->nsub = numchan;
  }

  /* Using the Effelsberg-Berkeley Pulsar Processor routines   */
  /*   NOTE:  This code is not yet implemented.                */

  if (cmd->ebppP) {

    /* OBS code for TEMPO */

    strcpy(obs, "EF");
  }

  /* Raw floating point data (already de-dispersed if radio data) */

  if (!cmd->ebppP && !cmd->pkmbP){

    /* Read the first header file and generate an infofile from it */

    if (idata.object) {
      printf("Folding a %s candidate from '%s'.\n", \
	     idata.object, cmd->argv[0]);
    } else {
      printf("Folding a candidate from '%s'.\n", cmd->argv[0]);
    }

    /* Some information about the size of the records */

    cmd->nsub = 1;
    numchan = 1;
    worklen = 1024;
    dt = idata.dt;
    dtdays = idata.dt / SECPERDAY;
    N = chkfilelen(infile, sizeof(float));

    /* Determine the number of records to use from the command line */

    lorec = (long) (cmd->startT * N + DBLCORRECT);
    hirec = (long) (cmd->endT * N + DBLCORRECT);
    N = hirec - lorec;
    numrec = N / worklen;

    /* The number of reads from the file we need for */
    /* each sub-integration.                         */ 

    reads_per_part = numrec / cmd->npart;

    /* Correct numrec so that each part will contain */
    /* the same number of records.                   */

    numrec = reads_per_part * cmd->npart;
    N = numrec * worklen;
    T = N * dt;
    endtime = T + 2 * TDT;

    if (cmd->nobaryP){
      if (idata.mjd_i && idata.mjd_f)
	baryepoch = (double) idata.mjd_i + 
	  idata.mjd_f + lorec * dt / SECPERDAY;
    } else {
      
      /* OBS code for TEMPO */
      
      strcpy(obs, "PK");
      
      /* Define the RA and DEC of the observation */
      
      ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
      ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);
      
      /* Topocentric and barycentric times of folding epoch data */
      
      if (idata.mjd_i && idata.mjd_f) {
	topoepoch = (double) idata.mjd_i + 
	  idata.mjd_f + lorec * dt / SECPERDAY;
	barycenter(&topoepoch, &baryepoch, &dtmp, 1, rastring,
		   decstring, obs, ephem);

      /* Correct the barycentric time for the dispersion delay.     */
      /* This converts the barycentric time to infinite frequency.  */

	barydispdt = delay_from_dm(idata.dm, idata.freq + 
				   (idata.num_chan - 1) * idata.chan_wid);
	baryepoch -= (barydispdt / SECPERDAY);
      }
    }

    /* The data collection routine to use */

    readrec_ptr = read_floats;
  }

  /* Read the pulsar database if needed */

  if (cmd->psrnameP) {
    np = read_database(&pdata);
    pnum = get_psr_at_epoch(cmd->psrname, baryepoch, &pdata, &psr);
    if (!pnum) {
      printf("The pulsar is not in the database.  Exiting.\n\n");
      exit(1);
    }
    if (psr.ntype & 8){  /* Checks if the pulsar is in a binary */
      binary = 1;
      orb = psr.orb;
      orb_baryepoch = psr.orb.t / SECPERDAY;
    }
    p = psr.p;
    pd = psr.pd;
    f = psr.f;
    fd = psr.fd;
    fdd = psr.fdd;
    strcpy(pname, psr.jname);
    
    /* If the user specifies all of the binaries parameters */	
    
  } else if (cmd->binaryP) {				
							
    /* Assume that the psr characteristics were measured at the time */
    /* of periastron passage (cmd->To)                               */
    
    difft = SECPERDAY * (baryepoch - cmd->To);				
    orb.p = cmd->pb;							
    orb.x = cmd->asinic;						
    orb.e = cmd->e;							
    orb.t = fmod(difft, orb.p);					
    if (orb.t < 0.0)							
      orb.t += orb.p;							
    orb.w = (cmd->w + difft * cmd->wdot / SECPERJULYR);		
    binary = 1;
    
  } else if (cmd->rzwcandP) {						

    if (!cmd->rzwfileP) {						
      printf("\nYou must enter a name for the rzw candidate ");	
      printf("file (-rzwfile filename)\n");				
      printf("Exiting.\n\n");						
      exit(1);							
    } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))){
      ii = (long) (cptr - cmd->rzwfile);
      cptr = (char *)malloc(ii + 1);
      cptr[ii] = '\0';
      strncpy(cptr, cmd->rzwfile, ii);
      fprintf(stderr, "\nAttempting to read '%s.inf'.  ", cptr);
      readinf(&rzwidata, cptr);
      free(cptr);
      fprintf(stderr, "Successful.\n");
      get_rzw_cand(cmd->rzwfile, cmd->rzwcand, &rzwcand);	
      f = (rzwcand.r - 0.5 * rzwcand.z) / 
	(rzwidata.dt * rzwidata.N);
      fd = rzwcand.z / ((rzwidata.dt * rzwidata.N) * 
			  (rzwidata.dt * rzwidata.N));

      /* Now correct for the fact that we may not be starting */
      /* to fold at the same start time as the rzw search.    */

      if (cmd->pkmbP)
	f += lorec * recdt * fd;
      else if (!cmd->pkmbP && !cmd->ebppP)
	f += lorec * dt * fd;
      p = 1.0 / f;
      pd = -fd / (f * f);
    } else {
      printf("\nCould not read the rzwfile.\nExiting.\n\n");
      exit(1);
    }

  }
  
  /* Determine the pulsar parameters to fold if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */
  
  if (!cmd->rzwcandP && !cmd->psrnameP) {

    if (cmd->pP) {
      p = cmd->p;
      f = 1.0 / p;
    }
    if (cmd->fP) {
      f = cmd->f;
      p = 1.0 / f;
    }
    if (cmd->pd != 0.0) {
      pd = cmd->pd;
      fd = -pd / (p * p);
    }
    if (cmd->fd != 0.0) {
      fd = cmd->fd;
      pd = -fd / (f * f);
    }
    if (cmd->pdd != 0.0) {
      pdd = cmd->pdd;
      fdd = 2 * pd * pd / (p * p * p) -
	pdd / (p * p);
    }
    if (cmd->fdd != 0.0) {
      fdd = cmd->fdd;
      pdd = 2 * fd * fd / (f * f * f) - 
	fdd / (f * f);
    }
  }
  
  /* Determine the length of the profile */				
  
  if (cmd->proflenP)
    proflen = cmd->proflen;
  else
    proflen = (long) (p / dt + 0.5);

  /* Determine the phase delays caused by the orbit if needed */

  if (binary) {
    
    /* Save the orbital solution every half second               */
    /* The times in *tp are now calculated as barycentric times. */
    /* Later, we will change them to topocentric times after     */
    /* applying corrections to Ep using TEMPO.                   */
    
    startE = keplars_eqn(orb.t, orb.p, orb.e, 1.0E-15);
    if (endtime > 2048) orbdt = 0.5;
    else orbdt = endtime / 4096.0;
    numbinpoints = (long) floor(endtime/orbdt + 0.5) + 1;
    Ep = dorbint(startE, numbinpoints, orbdt, &orb);
    tp = gen_dvect(numbinpoints);
    for (ii = 0; ii < numbinpoints; ii++) tp[ii] = ii * orbdt;

    /* Convert Eccentric anomaly to time delays */			
    
    orb.w *= DEGTORAD;
    E_to_phib(Ep, numbinpoints, &orb);
    numdelays = numbinpoints;
  }

  /* Output some informational data on the screen and to the */
  /* output file.                                            */
  
  fprintf(stdout, "\n");
  filemarker = stdout;
  for (ii = 0 ; ii < 1 ; ii++){
    if (cmd->psrnameP)
      fprintf(filemarker, 
	      "Pulsar                       =  %s\n", pname);
    if (topoepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (topo) epoch  (MJD)  =  %-17.11f\n", topoepoch);
    if (baryepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (bary) epoch  (MJD)  =  %-17.11f\n", baryepoch);
    fprintf(filemarker, 
	    "Data pt duration (dt)   (s)  =  %-.12f\n", dt);
    fprintf(filemarker, 
	    "Total number of data points  =  %-.0f\n", N);
    fprintf(filemarker, 
	    "Number of profile bins       =  %ld\n", proflen);
    fprintf(filemarker, 
	    "Folding period          (s)  =  %-.15f\n", p);
    if (pd != 0.0)
      fprintf(filemarker, 
	      "Folding p-dot         (s/s)  =  %-.10e\n", pd);
    if (pdd != 0.0)
      fprintf(filemarker, 
	      "Folding p-dotdot    (s/s^2)  =  %-.10e\n", pdd);
    fprintf(filemarker, 
	    "Folding frequency      (hz)  =  %-.12f\n", f);
    if (fd != 0.0)
      fprintf(filemarker, 
	      "Folding f-dot        (hz/s)  =  %-.8e\n", fd);
    if (fdd != 0.0)
      fprintf(filemarker, 
	      "Folding f-dotdot   (hz/s^2)  =  %-.8e\n", fdd);
    if (binary) {							
      fprintf(filemarker, 
	      "Orbital period          (s)  =  %-.10f\n", orb.p);
      fprintf(filemarker, 
	      "a*sin(i)/c (x)     (lt-sec)  =  %-.10f\n", orb.x);
      fprintf(filemarker, 
	      "Eccentricity                 =  %-.10f\n", orb.e);
      fprintf(filemarker, 
	      "Longitude of peri (w) (deg)  =  %-.10f\n", orb.w/DEGTORAD); 
      if (cmd->psrnameP){
	fprintf(filemarker, 
		"Epoch of periapsis    (MJD)  =  %-17.11f\n",
		baryepoch + orb_baryepoch);
      }
    }
  }

  printf("\nStarting work on '%s'...\n\n", cmd->argv[0]);
    
  /* Allocate and initialize some arrays and other information */
  
  data = gen_fvect(cmd->nsub * worklen);
  profs = gen_dvect(cmd->nsub * cmd->npart * proflen);
  stats = (foldstats *)malloc(sizeof(foldstats) * cmd->nsub * cmd->npart);
  for (ii = 0; ii < cmd->nsub * cmd->npart; ii++){
    for (jj = 0 ; jj < proflen; jj++)
      profs[ii * proflen + jj] = 0.0;
    stats[ii].numdata = 0.0;
    stats[ii].data_avg = 0.0;
    stats[ii].data_var = 0.0;
  }
  if (numdelays == 0) flags = 0;
    
  /* Move to the correct starting record */
  
  currentrec = lorec;
  if (cmd->pkmbP)
    skip_to_multibeam_rec(infile, lorec);
  else
    chkfileseek(infile, lorec, sizeof(float), SEEK_SET);

  /* Correct our fold parameters if we are barycentering */

  if (cmd->nobaryP) {
    foldf = f;  foldfd = fd;  foldfdd = fdd;
  } else {

    /* The number of topo to bary points to generate with TEMPO */
    
    numbarypts = endtime / TDT + 1;
    barytimes = gen_dvect(numbarypts);
    topotimes = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);

    /* topocentric times in days from data start */

    for (ii = 0 ; ii < numbarypts ; ii++)
      topotimes[ii] = topoepoch + (double) ii * TDT / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(topotimes, barytimes, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    printf("   Insure you check the file %s.tempo_out for\n", \
	   cmd->argv[0]);
    printf("   errors from TEMPO when complete.\n\n");

    /* Determine the avg v/c of the Earth's motion during the obs */

    avgvoverc = 0.0;
    for (ii = 0 ; ii < numbarypts - 1 ; ii++)
      avgvoverc += voverc[ii];
    avgvoverc /= (numbarypts - 1.0);
    free(voverc);
    printf("The average topocentric velocity is %.3g (units of c).\n", 
	   avgvoverc);

    /* Convert the barycentric folding parameters into topocentric */

    if ((info=bary2topo(topotimes, barytimes, numbarypts, 
		   f, fd, fdd, &foldf, &foldfd, &foldfdd))<0)
      printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
    printf("Topocentric folding frequency    (hz)  =  %-.12f\n", foldf);
    if (foldfd != 0.0)
      printf("Topocentric folding f-dot      (hz/s)  =  %-.8e\n", foldfd);
    if (foldfdd != 0.0)
      printf("Topocentric folding f-dotdot (hz/s^2)  =  %-.8e\n", foldfdd);
    printf("\n");

    /* Modify the binary delay times so they refer to */
    /* topocentric reference times.                   */
      
    if (binary){
      for (ii = 0; ii < numbinpoints; ii++){
	arrayoffset++;  /* Beware nasty NR zero-offset kludges! */
	dtmp = baryepoch + tp[ii] / SECPERDAY;
	hunt(barytimes, numbarypts, dtmp, &arrayoffset);
	arrayoffset--;
	tp[ii] = LININTERP(dtmp, barytimes[arrayoffset], 
			   barytimes[arrayoffset+1], 
			   topotimes[arrayoffset], 
			   topotimes[arrayoffset+1]);    
      }
      numdelays = numbinpoints;
      dtmp = tp[0];
      for (ii = 0 ; ii < numdelays ; ii++)
	tp[ii] = (tp[ii] - dtmp) * SECPERDAY;
    }
    free(barytimes);
    free(topotimes);
  }

  /* If this is 'raw' radio data, determine the dispersion delays */

  if (!strcmp(idata.band, "Radio")) {
    
    /* The observation frequencies */
    
    obsf = gen_dvect(numchan);
    obsf[0] = idata.freq;
    for (ii = 0; ii < numchan; ii++)
      obsf[ii] = obsf[0] + ii * idata.chan_wid;
    if (!cmd->nobaryP){
      for (ii = 0; ii < numchan; ii++)
	obsf[ii] = doppler(obsf[ii], avgvoverc);
    } 
    dispdts = subband_search_delays(numchan, cmd->nsub, cmd->dm,
				    idata.freq, idata.chan_wid, avgvoverc); 

    /* Convert the delays in seconds to delays in bins */

    for (ii = 0; ii < numchan; ii++)
      dispdts[ii] /= dt;
  }
  
  /* 
   *   Perform the actual folding of the data
   */

  proftime = worklen * dt;
  parttimes = gen_dvect(cmd->npart);
  printf("Folded %ld points of %.0f", totnumfolded, N);
  for (ii = 0; ii < cmd->npart; ii++){        /* sub-integrations in time  */
    parttimes[ii] = ii * reads_per_part * proftime;
    for (jj = 0; jj < reads_per_part; jj++){  /* reads per sub-integration */
      numread = readrec_ptr(infile, data, worklen, dispdts, 
			    cmd->nsub, numchan);
      /* tt is (topocentric) time from first point in sec */
      tt = parttimes[ii] + jj * proftime;
      for (kk = 0; kk < cmd->nsub; kk++)      /* frequency sub-bands */
	fold(data + kk * worklen, numread, dt, tt, 
	     profs + (ii * cmd->nsub + kk) * proflen, proflen, 
	     cmd->phs, foldf, foldfd, foldfdd, flags, Ep, tp, 
	     numdelays, NULL, &(stats[ii * cmd->nsub + kk]));
      totnumfolded += numread;
    }
    printf("\rFolded %ld points of %.0f", totnumfolded, N);
/* for (kk = 0; kk < cmd->nsub; kk++) */
/* quick_plot(profs + (ii * cmd->nsub + kk) * proflen, proflen); */
  }
  fclose(infile);

  /*
   *   Write the raw (unsummed) profiles
   */

  printf("\n\nWriting %s.\n", outfilenm);
  infile = chkfopen(outfilenm,"wb");
  chkfwrite(&topoepoch, sizeof(double), 1, infile);
  chkfwrite(&baryepoch, sizeof(double), 1, infile);
  chkfwrite(&f, sizeof(double), 1, infile);
  chkfwrite(&fd, sizeof(double), 1, infile);
  chkfwrite(&fdd, sizeof(double), 1, infile);
  chkfwrite(&foldf, sizeof(double), 1, infile);
  chkfwrite(&foldfd, sizeof(double), 1, infile);
  chkfwrite(&foldfdd, sizeof(double), 1, infile);
  dtmp = cmd->npart;
  chkfwrite(&dtmp, sizeof(double), 1, infile);
  dtmp = cmd->nsub;
  chkfwrite(&dtmp, sizeof(double), 1, infile);
  dtmp = proflen;
  chkfwrite(&dtmp, sizeof(double), 1, infile);
  chkfwrite(profs, sizeof(double), cmd->nsub * cmd->npart * proflen, 
	    infile);
  chkfwrite(stats, sizeof(foldstats), cmd->nsub * cmd->npart, infile);
  fclose(infile);

  /*
   *   Perform the candidate optimization search
   */

  printf("\nOptimizing...\n\n");
  {
    int numtrials, pdelay, pddelay, profindex;
    double dphase, dp, dpd, po, pdo, pddo, lop, lopd, pofact;
    double *pdprofs, currentfd, *currentprof;
    foldstats currentstats;

    /* The number of trials for the P-dot and P searches */

    numtrials = 2 * proflen + 1;

    /* Initialize a bunch of variables */

    periods = gen_dvect(numtrials);
    pdots = gen_dvect(numtrials);
    pdprofs = gen_dvect(cmd->npart * proflen);
    currentprof = gen_dvect(proflen);
    bestprof = gen_dvect(proflen);
    initialize_foldstats(&beststats);

    /* Convert the folding freqs and derivs to periods */
    
    switch_f_and_p(foldf, foldfd, foldfdd, &po, &pdo, &pddo);
    
    /* Our P and P-dot steps are the changes that cause the pulse */
    /* to be delayed 1 bin between the first and last bins.       */
    
    dphase = po / proflen;
    pofact = 1.0 / (po * po);
    dtmp = po / T;
    dp = dphase / T;
    dpd = 2.0 * dphase * dtmp * dtmp;
    lop = po - (numtrials - 1) / 2 * dp;
    lopd = pdo - (numtrials - 1) / 2 * dpd;
    for (ii = 0; ii < numtrials; ii++){
      periods[ii] = lop + ii * dp;
      pdots[ii] = lopd + ii * dpd;
    }

    /* If we are searching through DM space */
    
    if (cmd->nsub > 1){
      int numdmtrials, *dmdelays;
      double lodm, ddm, hifdelay, *ddprofs, *subbanddelays;
      foldstats *ddstats;
      
      dmdelays = gen_ivect(cmd->nsub);
      ddprofs = gen_dvect(cmd->npart * proflen);
      ddstats = (foldstats *)malloc(cmd->npart * sizeof(foldstats));
      
      /* Our DM step is the change in DM that would cause the pulse */
      /* to be delayed 1 phasebin at the lowest frequency.          */
      
      ddm = dm_from_delay(dphase, obsf[0]);
      
      /* Insure that we don't try a dm < 0.0 */
      
      if (cmd->dm - proflen * ddm < 0.0){
	numdmtrials = (int)(cmd->dm / ddm) + proflen + 1;
	lodm = cmd->dm - (int)(cmd->dm / ddm);
      } else {
	numdmtrials = 4 * proflen + 1;
	lodm = cmd->dm - (numdmtrials - 1) / 2 * ddm;
      }
      dms = gen_dvect(numdmtrials);
      
      /* De-disperse and combine the subbands */
      
      for (ii = 0; ii < numdmtrials; ii++){  /* Loop over DMs */
	dms[ii] = lodm + ii * ddm;
	hifdelay = delay_from_dm(dms[ii], obsf[numchan - 1]);
	subbanddelays = subband_delays(numchan, cmd->nsub, dms[ii], 
				       idata.freq, idata.chan_wid,
				       avgvoverc);
	for (jj = 0; jj < cmd->nsub; jj++)
	  dmdelays[jj] = ((int) ((subbanddelays[jj] - hifdelay) / 
				 dphase + 0.5)) % proflen;
	free(subbanddelays);
	combine_subbands(profs, stats, cmd->npart, cmd->nsub, proflen, 
			 dmdelays, ddprofs, ddstats);
	
	/* Perform the P-dot and Period searches */

	for (jj = 0; jj < numtrials; jj++){
	  currentfd = -pdots[jj] * pofact;

	  /* Correct each part for the current pdot */

	  for (kk = 0; kk < cmd->npart; kk++){
	    profindex = kk * proflen;
	    pddelay = (int) (((currentfd - foldfd) * parttimes[kk] * 
			     parttimes[kk]) / (2.0 * dphase) + 0.5);
	    shift_prof(ddprofs + profindex, proflen, pddelay, 
		       pdprofs + profindex);
	  }

	  /* Search over the periods */

	  for (kk = 0; kk < numtrials; kk++){
	    pdelay = kk - (numtrials - 1) / 2;
	    combine_profs(pdprofs, ddstats, cmd->npart, proflen, pdelay, 
			  currentprof, &currentstats);
	    if (currentstats.redchi > beststats.redchi){
	      bestdm = dms[ii];
	      bestp = periods[kk];
	      bestpd = pdots[jj];
	      beststats = currentstats;
	      memcpy(bestprof, currentprof, sizeof(double) * proflen);

printf("%ld %ld:  dm = %f  p = %17.15f   pd = %12.6e  reduced chi = %f\n", 
       kk, jj, dms[ii], bestp, bestpd, beststats.redchi);

	    }
	  }
	}
      }
      free(dmdelays);
      free(ddprofs);
      free(ddstats);
      
      /* We are not searching through DM space */
      
    } else {

      /* Perform the P-dot and Period searches */

      for (jj = 0; jj < numtrials; jj++){
	currentfd = -pdots[jj] * pofact;
	
	/* Correct each part for the current pdot */
	
	for (kk = 0; kk < cmd->npart; kk++){
	  profindex = kk * proflen;
	  pddelay = (int) (((currentfd - foldfd) * parttimes[kk] * 
			    parttimes[kk]) / (2.0 * dphase) + 0.5);
	  shift_prof(profs + profindex, proflen, pddelay, 
		     pdprofs + profindex);
	}

	/* Search over the periods */
	
	for (kk = 0; kk < numtrials; kk++){
	  pdelay = kk - (numtrials - 1) / 2;
	  combine_profs(pdprofs, stats, cmd->npart, proflen, pdelay, 
			currentprof, &currentstats);
	  if (currentstats.redchi > beststats.redchi){
	    bestp = periods[kk];
	    bestpd = pdots[jj];
	    beststats = currentstats;
	    memcpy(bestprof, currentprof, sizeof(double) * proflen);

printf("%ld %ld:  p = %17.15f   pd = %12.6e  reduced chi = %f\n", 
       kk, jj, bestp, bestpd, beststats.redchi);

	  }
	}
      }
    }

quick_plot(bestprof, proflen);

    free(pdprofs);
    free(currentprof);
  }

  /*
   *   Plot our results
   */

  printf("Done.\n\n");

  /* Free our memory  */

  free(data);
  free(profs);
  free(stats);
  free(rootfilenm);
  free(outfilenm);
  free(parttimes);
  free(periods);
  free(pdots);
  free(bestprof);
  if (cmd->nsub > 1) free(dms);
  if (binary){
    free(Ep);
    free(tp);
  }
  if (obsf) free(obsf);
  if (dispdts) free(dispdts);
  if (idata.onoff) free(idata.onoff);
  return (0);
}

    
int read_floats(FILE *file, float *data, int numpts,
		double *dispdelays, int numsubbands, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
  /* The following 2 lines just get rid of some compiler warnings */

  *dispdelays = *dispdelays;
  numsubbands = numsubbands;

  /* Read the raw data and return numbar read */

  return chkfread(data, sizeof(float),
		  (unsigned long) (numpts * numchan), file) / numchan;
}




