#include "presto.h"
#include "prepfold_cmd.h"
#include "multibeam.h"

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

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infile=NULL, *filemarker;
  float *data=NULL;
  double p=0.0, pd=0.0, pdd=0.0, f=0.0, fd=0.0, fdd=0.0;
  double difft, tt, nc, pl, recdt, *dispdts=NULL;
  double orb_baryepoch=0.0, topoepoch=0.0, baryepoch=0.0, barydispdt;
  double dtmp, *Ep=NULL, *tp=NULL, startE=0.0, orbdt=1.0;
  double tdf=0.0, N=0.0, dt=0.0, T, endtime=0.0, dtdays, avg_voverc;
  double *voverc=NULL, *tobsf=NULL;
  double *profs=NULL, *barytimes=NULL, *topotimes=NULL;
  char obs[3], ephem[10], *outfilenm, *rootfilenm;
  char pname[30], rastring[50], decstring[50], *cptr;
  int numchan=1, binary=0, np, pnum, numdelays, slen, ptsperrec=1;
  long ii, jj, kk, numbarypts=0, worklen=0, numread=0, reads_per_part;
  long numfolded=0, totnumfolded=0;
  long lorec=0, hirec=0, numrecs=0, totnumrecs=0;
  long numbinpoints=0, proflen, currentrec;
  unsigned long numrec=0, arrayoffset=0;
  multibeam_tapehdr hdr;
  fourierprops rzwcand;
  orbitparams orb;
  psrparams psr;
  psrdatabase pdata;
  infodata idata, rzwidata;
  Cmdline *cmd;
  foldstats *stats, beststats, currentstats;

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
    T = numrec * recdt;

    /* Determine the number of records to use from the command line */

    lorec = (long) (cmd->startT * numrec + DBLCORRECT);
    hirec = (long) (cmd->endT * numrec + DBLCORRECT);
    numrec = hirec - lorec;
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

    /* Read the first header file and generate an infofile from it */

    /* OBS code for TEMPO */

    strcpy(obs, "EF");

    /* The data collection routine to use */

    /* readrec_ptr = read_ebpp; */

    /* The number of data points to work with at a time */

    worklen = 1024;
  }

  /* Raw floating point data (already de-dispersed if radio data) */
  /* and already barycentered.                                    */

  if (!cmd->ebppP && !cmd->pkmbP){

    if (!cmd->nobaryP){
      printf("\nIf you are trying to fold single channel data, \n");
      printf("the data must be barycentered and you must specify\n");
      printf("'-nobary' on the command line.  Exiting.\n\n");
      exit(0);
    }

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
    T = N * dt;
    numrec = N / worklen;
    endtime = T + 2 * TDT;
    if (idata.mjd_i && idata.mjd_f)
      baryepoch = (double) idata.mjd_i + 
	idata.mjd_f + lorec * dt / SECPERDAY;

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
      strncpy(cptr, cmd->argv[0], ii);
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
  
  if (cmd->proflenP) {
    proflen = cmd->proflen;
  } else {
    proflen = (long) (p / dt + 0.5);
  }

  /* Determine the frequency shifts caused by the orbit if needed */

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

  if (!strcmp(idata.band, "Radio")) {
    
    /* The topocentric spacing between channels */
    
    tdf = idata.chan_wid;
    
    /* The topocentric observation frequencies */
    
    tobsf = gen_dvect(numchan);
    tobsf[0] = idata.freq;
    for (ii = 0; ii < numchan; ii++)
      tobsf[ii] = tobsf[0] + ii * tdf;
    dispdts = subband_search_delays(numchan, cmd->nsub, cmd->dm,
				    tobsf[0], tdf); 
  }
  
  printf("Starting work on '%s'...\n\n", cmd->argv[0]);
    
  /* Allocate and initialize some arrays and other information */
  
  data = gen_fvect(cmd->nsub * worklen);
  profs = gen_dvect(cmd->nsub * cmd->npart * proflen);
  stats = (foldstats *)malloc(sizeof(stats) * cmd->nsub * cmd->npart);
  for (ii = 0 ; ii < cmd->nsub * cmd->npart; ii++){
    for (jj = 0 ; jj < proflen; jj++)
      profs[jj] = 0.0;
    stats[ii].numdata = 0.0;
    stats[ii].data_avg = 0.0;
    stats[ii].data_var = 0.0;
  }
  currentrec = 0;
    
  /* Move to the correct starting record */
  
  printf("Folded %ld points of %ld", totnumfolded, N);
  if (cmd->pkmbP)
    currentrec = skip_to_multibeam_rec(infile, lorec);
  else
    currentrec = chkfileseek(infile, sizeof(float) * lorec, 
			     SEEKSET) / sizeof(float);

  /* The number of reads from the file we need for */
  /* each sub-integration.                         */
  
  reads_per_part = numrec / cmd->npart;

  /* Main loop if we are not barycentering... */
  
  if (cmd->nobaryP) {
    
    /* Step through the sub-integrations of time */

    for (ii = 0; ii < cmd->npart; ii++){

      /* Step through the records in the sub-integration */
      
      for (jj = 0; jj < reads_per_part; ii++){
	
	printf("\rFolded %ld points of %ld", totnumfolded, N);
	fflush(stdout);
	
	/* Read the next record (or records) */
	
	numread = readrec_ptr(infile, data, worklen, dispdts, 
			      cmd->nsub, numchan);

	/* tt is (topocentric) seconds from first point */
      
	tt = (ii * reads_per_part + jj) * worklen * dt;
      
	/* Step through channels */
      
	for (kk = 0; kk < numchan ; kk++)
	  fold(data + kk * worklen, numread, dt, tt, 
	       profs + kk * proflen, proflen, cmd->phs, 
	       f, fd, fdd, 1, Ep, tp, numdelays, NULL, &stats);
      totnumfolded += numfolded;
      numrecs++;
    }


    /*****  Need to do this  *****/
    
  /* Main loop if we are barycentering... */

  } else {
    
    /* The number of topo to bary points to generate with TEMPO */
    
    numbarypts = endtime / TDT + 1;

    /* Allocate some arrays */
    
    barytimes = gen_dvect(numbarypts);
    topotimes = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);

    /* topocentric times in days from data start */

    for (ii = 0 ; ii < numbarypts ; ii++)
      topotimes[ii] = topoepoch + (double) ii * TDT / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("\nGenerating barycentric corrections...\n");
    barycenter(topotimes, barytimes, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    printf("   Insure you check the file %s.tempo_out for\n", \
	   cmd->argv[0]);
    printf("   errors from TEMPO when complete.\n\n");

    /* Determine the average V/c */

    avg_voverc = 0.0;
    for (ii = 0 ; ii < numbarypts - 1 ; ii++)
      avg_voverc += voverc[ii];
    avg_voverc /= (numbarypts - 1.0);
    free(voverc);

    if (binary){

      /* Modify the binary times so that they refer to topocentric   */
      /* reference times.  This way we can barycenter while we fold. */
      
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
    } else {
      tp = topotimes;
      Ep = barytimes;
      numdelays = numbarypts;
    }

    /* Step through records */
    
    for (ii = lorec; ii <= hirec; ii++){
      
      printf("\rCompleted record #%ld of %ld", numrecs, totnumrecs);
      fflush(stdout);
      
      /* Read the next record (or records) */
      
      readrec_ptr(infile, data, worklen, dispdts, cmd->nsub, 
		  numchan);
      
      /* tt is topocentric seconds from data start */
      
      tt = ii * worklen * dt;
      
      /* Step through channels */
      
      for (jj = 0 ; jj < numchan ; jj++){
	numfolded = 0;
	fold(data + jj * worklen, worklen, dt, tt, 
	     profs + jj * proflen, proflen, cmd->phs, 
	     f, fd, fdd, 1, Ep, tp, numdelays, NULL, &stats);
      }
      totnumfolded += numfolded;
      numrecs++;
    }
    free(barytimes);
    free(topotimes);
  }
  
  printf("\rCompleted record #%ld of %ld\n\n", 
	 totnumrecs, totnumrecs);
  printf("Folded %d profiles with %ld points each.\n", 
	 numchan, totnumfolded);

  fclose(infile);

  /* Write our results. */

  printf("\nWriting %s.\n", outfilenm);
  infile = chkfopen(outfilenm,"wb");
  nc = numchan;
  pl = proflen;
  chkfwrite(&nc, sizeof(double), 1, infile);
  chkfwrite(&pl, sizeof(double), 1, infile);
  chkfwrite(&p, sizeof(double), 1, infile);
  chkfwrite(&topoepoch, sizeof(double), 1, infile);
  chkfwrite(&baryepoch, sizeof(double), 1, infile);
  chkfwrite(&idata.freq, sizeof(double), 1, infile);
  chkfwrite(&idata.chan_wid, sizeof(double), 1, infile);
  chkfwrite(profs, sizeof(double), 
	    (unsigned long) (numchan * proflen), infile);
  fclose(infile);
  printf("Done.\n\n");

  /* Free our memory  */

  free(data);
  free(profs);
  free(rootfilenm);
  free(outfilenm);
  if (binary){
    free(Ep);
    free(tp);
  }
  if (tobsf) free(tobsf);
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




