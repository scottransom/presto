#include "prepfold.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* 
 * The main program 
 */

int main(int argc, char *argv[])
{
  FILE *infile=NULL, *filemarker, *binproffile;
  float *data=NULL;
  double f=0.0, fd=0.0, fdd=0.0, foldf=0.0, foldfd=0.0, foldfdd=0.0;
  double recdt=0.0, barydispdt, N=0.0, T=0.0, proftime;
  double *obsf=NULL, *dispdts=NULL, *parttimes, *Ep=NULL, *tp=NULL;
  double *barytimes=NULL, *topotimes=NULL, *bestprof, dtmp;
  char *plotfilenm, *outfilenm, *rootnm, *binproffilenm;
  char obs[3], ephem[6], pname[30], rastring[50], decstring[50];
  int numchan=1, binary=0, numdelays=0, numbarypts=0;
  int info, ptsperrec=1, flags=1;
  long ii, jj, kk, worklen=0, numread=0, reads_per_part=0;
  long totnumfolded=0, lorec=0, hirec=0, numbinpoints=0, currentrec=0;
  unsigned long numrec=0, arrayoffset=0;
  infodata idata;
  foldstats beststats;
  prepfoldinfo search;
  Cmdline *cmd;

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
  printf(" Used for DM, Period, and P-dot tweaking of PSR candidates.\n");
  printf("                 by Scott M. Ransom\n");
  printf("                    14 Jan 2000\n\n");

  init_prepfoldinfo(&search);

  /* Manipulate the file names we will use  */

  {
    char *path, *root, *suffix;
    int slen;

    /* Split argv[0] into a path and a file name */

    split_path_file(cmd->argv[0], &path, &search.filenm);
    free(path);

    /* Split the filename into a rootname and a suffix */

    if (split_root_suffix(search.filenm, &root, &suffix)==0){
      printf("\nThe input filename (%s) must have a suffix!\n\n", 
	     search.filenm);
      free(search.filenm);
      exit(1);
    }
    free(suffix);

    /* Split argv[0] into a rootname and a suffix */

    split_root_suffix(cmd->argv[0], &rootnm, &suffix);
    free(suffix);

    if (!cmd->pkmbP && !cmd->ebppP){
      printf("Reading input data from '%s'\n", search.filenm);
      printf("Reading information from '%s.inf'\n\n", root);
      
      /* Read the info file if available */
      
      readinf(&idata, rootnm);
    } else if (cmd->pkmbP){
      printf("Reading Parkes Multibeam data from '%s'\n\n", 
	     search.filenm);
    } else if (cmd->ebppP){
      printf("Reading New Effelsberg PSR data from '%s'\n\n", 
	     search.filenm);
    }
    free(root);

    /* Determine the candidate name */
    
    if (cmd->psrnameP){
      slen = strlen(cmd->psrname) + 5;
      search.candnm = (char *)calloc(slen, sizeof(char));
      sprintf(search.candnm, "PSR_%s", cmd->psrname);
    } else if (cmd->rzwcandP) {						
      slen = 20;
      search.candnm = (char *)calloc(slen, sizeof(char));
      sprintf(search.candnm, "RZW_Cand_%d", cmd->rzwcand);
    } else {
      slen = 20;
      search.candnm = (char *)calloc(slen, sizeof(char));
      if (cmd->pP)
	sprintf(search.candnm, "%.2fms_Cand", cmd->p * 1000.0);
      else if (cmd->fP)
	sprintf(search.candnm, "%.2fHz_Cand", cmd->f);
      else {
	printf("\nYou must specify candidate parameters (i.e period).\n\n");
	exit(1);
      }
    }

    /* Determine the output and plot file names */
    
    slen = strlen(rootnm) + strlen(search.candnm) + 6;
    outfilenm = (char *)calloc(slen, sizeof(char));
    sprintf(outfilenm, "%s_%s.pfd", rootnm, search.candnm);
    plotfilenm = (char *)calloc(slen + 3, sizeof(char));
    sprintf(plotfilenm, "%s_%s.pfd.ps", rootnm, search.candnm);
    binproffilenm = (char *)calloc(slen + 9, sizeof(char));
    sprintf(binproffilenm, "%s_%s.pfd.binprofs", rootnm, search.candnm);
    search.pgdev = (char *)calloc(slen + 7, sizeof(char));
    sprintf(search.pgdev, "%s/CPS", plotfilenm);
    printf("Folding a %s candidate.\n\n", search.candnm);
    printf("Output data file is '%s'.\n", outfilenm);
    printf("Output plot file is '%s'.\n", plotfilenm);
    printf("Raw profile file is '%s'.\n", binproffilenm);
  }

  /* Open the raw data file */
  
  infile = chkfopen(cmd->argv[0], "rb");
  binproffile = chkfopen(binproffilenm, "wb");
    
  /* What ephemeris will we use?  (Default is DE200) */
  
  if (cmd->de405P)
    strcpy(ephem, "DE405");
  else
    strcpy(ephem, "DE200");
    
  /* Set-up values if we are using the Parkes Multibeam System */

  if (cmd->pkmbP){
    multibeam_tapehdr hdr;

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    rewind(infile);
    multibeam_hdr_to_inf(&hdr, &idata);
    idata.dm = cmd->dm;

    /* Some information about the size of the records */

    numrec = chkfilelen(infile, RECLEN);
    ptsperrec = DATLEN * 8 / idata.num_chan;

    /* OBS code for TEMPO */

    strcpy(obs, "PK");
    search.telescope = (char *)calloc(20, sizeof(char));
    strcpy(search.telescope, "Parkes Multibeam");

    /* Define the RA and DEC of the observation */
  
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Define some variables */

    search.dt = idata.dt;
    recdt = search.dt * ptsperrec;

    /* Determine the number of records to use from the command line */

    search.startT = cmd->startT;
    search.endT = cmd->endT;
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

    /* Topocentric and barycentric times of folding epoch data */

    if (idata.mjd_i) {
      search.tepoch = (double) idata.mjd_i + idata.mjd_f + 
	lorec * recdt / SECPERDAY;
      barycenter(&search.tepoch, &search.bepoch, &dtmp, 1, rastring,
		 decstring, obs, ephem);

      /* Correct the barycentric time for the dispersion delay.     */
      /* This converts the barycentric time to infinite frequency.  */

      barydispdt = delay_from_dm(cmd->dm, idata.freq + 
				 (idata.num_chan - 1) * idata.chan_wid);
      search.bepoch -= (barydispdt / SECPERDAY);
    }

    /* The data collection routine to use */

    readrec_ptr = read_multibeam_subbands;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    worklen = ptsperrec;
  }

  /* Using the Effelsberg-Berkeley Pulsar Processor routines   */
  /*   NOTE:  This code is not yet implemented.                */

  if (cmd->ebppP) {

    /* OBS code for TEMPO */

    strcpy(obs, "EF");
    search.telescope = (char *)calloc(20, sizeof(char));
    strcpy(search.telescope, "Effelsberg");
  }

  /* Raw floating point data (already de-dispersed if radio data) */

  if (!cmd->ebppP && !cmd->pkmbP){

    /* Some information about the size of the records */

    cmd->nsub = 1;
    numchan = 1;
    worklen = 1024;
    search.dt = idata.dt;
    N = chkfilelen(infile, sizeof(float));

    /* Determine the number of records to use from the command line */

    search.startT = cmd->startT;
    search.endT = cmd->endT;
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
    T = N * search.dt;

    /* Until I figure out a better way to do this... */

    search.telescope = (char *)calloc(strlen(idata.telescope)+1, 
				      sizeof(char));
    strcpy(search.telescope, idata.telescope);

    if (cmd->nobaryP){
      if (idata.mjd_i){
	if (idata.bary)
	  search.bepoch = (double) idata.mjd_i + 
	    idata.mjd_f + lorec * search.dt / SECPERDAY;
	else
	  search.tepoch = (double) idata.mjd_i + 
	    idata.mjd_f + lorec * search.dt / SECPERDAY;
      }
    } else {
      
      /* OBS code for TEMPO */
      
      if (cmd->obscodeP)
	strcpy(obs, cmd->obscode);
      else {
	printf("\nIf you want to barycenter you must specify an \n");
	printf("observatory found in TEMPO using the '-obs' argument.\n\n");
	exit(1);
      }
      
      /* Define the RA and DEC of the observation */
      
      ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
      ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);
      
      /* Topocentric and barycentric times of folding epoch */

      if (idata.mjd_i) {
	search.tepoch = (double) idata.mjd_i + 
	  idata.mjd_f + lorec * search.dt / SECPERDAY;
	barycenter(&search.tepoch, &search.bepoch, &dtmp, 1, rastring,
		   decstring, obs, ephem);

	if (!strcmp(idata.band, "Radio")) {

	  /* Correct the barycentric time for the dispersion delay.    */
	  /* This converts the barycentric time to infinite frequency. */

	  barydispdt = delay_from_dm(idata.dm, idata.freq + 
				     (idata.num_chan - 1) * 
				     idata.chan_wid);
	  search.bepoch -= (barydispdt / SECPERDAY);
	}
      }
    }

    /* The data collection routine to use */

    readrec_ptr = read_floats;
  }

  /* Read the pulsar database if needed */

  if (cmd->psrnameP) {
    int np, pnum;
    psrparams psr;
    psrdatabase pdata;

    np = read_database(&pdata);
    if (search.bepoch==0.0){
      printf("\nYou must not use the '-nobary' flag if you want to\n");
      printf("access the pulsar database.  Exiting.\n\n");
      exit(1);
    }
    pnum = get_psr_at_epoch(cmd->psrname, search.bepoch, &pdata, &psr);
    if (!pnum) {
      printf("The pulsar is not in the database.  Exiting.\n\n");
      exit(1);
    }
    if (psr.ntype & 8){  /* Checks if the pulsar is in a binary */
      binary = 1;
      search.orb = psr.orb;
    }
    search.bary.p1 = psr.p;
    search.bary.p2 = psr.pd;
    search.bary.p3 = psr.pdd;
    f = psr.f;
    fd = psr.fd;
    fdd = psr.fdd;
    strcpy(pname, psr.jname);
    
    /* If the user specifies all of the binaries parameters */	
    
  } else if (cmd->binaryP) {				
							
    /* Assume that the psr characteristics were measured at the time */
    /* of periastron passage (cmd->To)                               */
    
    if (search.bepoch==0.0)
      dtmp = SECPERDAY * (search.tepoch - cmd->To);
    else 
      dtmp = SECPERDAY * (search.bepoch - cmd->To);
    search.orb.p = cmd->pb;
    search.orb.x = cmd->asinic;
    search.orb.e = cmd->e;
    search.orb.t = fmod(dtmp, search.orb.p);
    if (search.orb.t < 0.0)
      search.orb.t += search.orb.p;
    search.orb.w = (cmd->w + dtmp * cmd->wdot / SECPERJULYR);
    binary = 1;

  } else if (cmd->rzwcandP) {
    fourierprops rzwcand;
    infodata rzwidata;
    char *cptr;

    if (!cmd->rzwfileP) {					
      printf("\nYou must enter a name for the rzw candidate ");
      printf("file (-rzwfile filename)\n");
      printf("Exiting.\n\n");
      exit(1);
    } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))){
      ii = (long) (cptr - cmd->rzwfile);
      cptr = (char *)calloc(ii + 1, sizeof(char));
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
	f += lorec * search.dt * fd;
      if (rzwidata.bary)
	switch_f_and_p(f, fd, fdd, &search.bary.p1, \
		       &search.bary.p2, &search.bary.p3);
      else
	switch_f_and_p(f, fd, fdd, &search.topo.p1, \
		       &search.topo.p2, &search.topo.p3);
    } else {
      printf("\nCould not read the rzwfile.\nExiting.\n\n");
      exit(1);
    }
    if (rzwidata.onoff) free(idata.onoff);
  }
  
  /* Determine the pulsar parameters to fold if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */
  
  if (!cmd->rzwcandP && !cmd->psrnameP) {
    double p=0.0, pd=0.0, pdd=0.0;

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
    if (idata.bary){
      search.bary.p1 = p;
      search.bary.p2 = pd;
      search.bary.p3 = pdd;
    } else {
      search.topo.p1 = p;
      search.topo.p2 = pd;
      search.topo.p3 = pdd;
    }
  }
  
  /* Determine the length of the profile */				
  
  if (cmd->proflenP)
    search.proflen = cmd->proflen;
  else
    if (idata.bary || cmd->psrnameP)
      search.proflen = (long) (search.bary.p1 / search.dt + 0.5);
    else
      search.proflen = (long) (search.topo.p1 / search.dt + 0.5);

  /* Determine the phase delays caused by the orbit if needed */

  if (binary) {
    double orbdt = 1.0, startE=0.0;

    /* Save the orbital solution every half second               */
    /* The times in *tp are now calculated as barycentric times. */
    /* Later, we will change them to topocentric times after     */
    /* applying corrections to Ep using TEMPO.                   */
    
    startE = keplars_eqn(search.orb.t, search.orb.p, 
			 search.orb.e, 1.0E-15);
    if (T > 2048) orbdt = 0.5;
    else orbdt = T / 4096.0;
    numbinpoints = (long) floor(T/orbdt + 0.5) + 1;
    Ep = dorbint(startE, numbinpoints, orbdt, &search.orb);
    tp = gen_dvect(numbinpoints);
    for (ii = 0; ii < numbinpoints; ii++) tp[ii] = ii * orbdt;

    /* Convert Eccentric anomaly to time delays */			
    
    search.orb.w *= DEGTORAD;
    E_to_phib(Ep, numbinpoints, &search.orb);
    numdelays = numbinpoints;
    if (search.bepoch == 0.0)
      search.orb.t = -search.orb.t / SECPERDAY + search.tepoch;
    else 
      search.orb.t = -search.orb.t / SECPERDAY + search.bepoch;
  }

  /* Output some informational data on the screen and to the */
  /* output file.                                            */
  
  fprintf(stdout, "\n");
  filemarker = stdout;
  for (ii = 0 ; ii < 1 ; ii++){
    double p, pd, pdd;

    if (cmd->psrnameP)
      fprintf(filemarker, 
	      "Pulsar                       =  %s\n", pname);
    if (search.tepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (topo) epoch  (MJD)  =  %-17.11f\n", search.tepoch);
    if (search.bepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (bary) epoch  (MJD)  =  %-17.11f\n", search.bepoch);
    fprintf(filemarker, 
	    "Data pt duration (dt)   (s)  =  %-.12f\n", search.dt);
    fprintf(filemarker, 
	    "Total number of data points  =  %-.0f\n", N);
    fprintf(filemarker, 
	    "Number of profile bins       =  %d\n", search.proflen);
    switch_f_and_p(f, fd, fdd, &p, &pd, &pdd);
    fprintf(filemarker, "Folding period          (s)  =  %-.15f\n", p);
    if (pd != 0.0)
      fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10e\n", pd);
    if (pdd != 0.0)
      fprintf(filemarker, "Folding p-dotdot    (s/s^2)  =  %-.10e\n", pdd);
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
	      "Orbital period          (s)  =  %-.10f\n", search.orb.p);
      fprintf(filemarker, 
	      "a*sin(i)/c (x)     (lt-sec)  =  %-.10f\n", search.orb.x);
      fprintf(filemarker, 
	      "Eccentricity                 =  %-.10f\n", search.orb.e);
      fprintf(filemarker, 
	      "Longitude of peri (w) (deg)  =  %-.10f\n", 
	      search.orb.w / DEGTORAD); 
      fprintf(filemarker, 
	      "Epoch of periapsis    (MJD)  =  %-17.11f\n", search.orb.t);
    }
  }

  printf("\nStarting work on '%s'...\n\n", search.filenm);
    
  /* Allocate and initialize some arrays and other information */
  
  search.nsub = cmd->nsub;
  search.npart = cmd->npart;
  data = gen_fvect(cmd->nsub * worklen);
  search.rawfolds = gen_dvect(cmd->nsub * cmd->npart * search.proflen);
  search.stats = (foldstats *)malloc(sizeof(foldstats) * 
				     cmd->nsub * cmd->npart);
  for (ii = 0 ; ii < cmd->npart * cmd->nsub * search.proflen; ii++)
    search.rawfolds[ii] = 0.0;
  for (ii = 0 ; ii < cmd->npart * cmd->nsub; ii++){
    search.stats[ii].numdata = 0.0;
    search.stats[ii].data_avg = 0.0;
    search.stats[ii].data_var = 0.0;
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
    double *voverc;

    /* The number of topo to bary points to generate with TEMPO */
    
    numbarypts = T / TDT + 10;
    barytimes = gen_dvect(numbarypts);
    topotimes = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);

    /* topocentric times in days from data start */

    for (ii = 0 ; ii < numbarypts ; ii++)
      topotimes[ii] = search.tepoch + (double) ii * TDT / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(topotimes, barytimes, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    printf("   Insure you check the file %s.tempo_out for\n", \
	   cmd->argv[0]);
    printf("   errors from TEMPO when complete.\n\n");

    /* Determine the avg v/c of the Earth's motion during the obs */

    for (ii = 0 ; ii < numbarypts - 1 ; ii++)
      search.avgvoverc += voverc[ii];
    search.avgvoverc /= (numbarypts - 1.0);
    free(voverc);
    printf("The average topocentric velocity is %.3g (units of c).\n\n", 
	   search.avgvoverc);
    printf("Barycentric folding frequency    (hz)  =  %-.12f\n", f);
    printf("Barycentric folding f-dot      (hz/s)  =  %-.8e\n", fd);
    printf("Barycentric folding f-dotdot (hz/s^2)  =  %-.8e\n", fdd);

    /* Convert the barycentric folding parameters into topocentric */

    if ((info=bary2topo(topotimes, barytimes, numbarypts, 
			f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
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
	dtmp = search.bepoch + tp[ii] / SECPERDAY;
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
  }

  /* Record the f, fd, and fdd we used to do the raw folding */

  if (idata.bary)
    search.fold.pow = 1.0;
  search.fold.p1 = foldf;
  search.fold.p2 = foldfd;
  search.fold.p3 = foldfdd;

  /* If this is 'raw' radio data, determine the dispersion delays */

  if (!strcmp(idata.band, "Radio")) {
    
    /* The observation frequencies */

    obsf = gen_dvect(numchan);
    obsf[0] = idata.freq;
    search.numchan = numchan;
    search.lofreq = idata.freq;
    search.chan_wid = idata.chan_wid;
    for (ii = 0; ii < numchan; ii++)
      obsf[ii] = obsf[0] + ii * idata.chan_wid;
    if (!cmd->nobaryP){
      for (ii = 0; ii < numchan; ii++)
	obsf[ii] = doppler(obsf[ii], search.avgvoverc);
    } 
    dispdts = subband_search_delays(numchan, cmd->nsub, cmd->dm,
				    idata.freq, idata.chan_wid, 
				    search.avgvoverc); 

    /* Convert the delays in seconds to delays in bins */

    for (ii = 0; ii < numchan; ii++)
      dispdts[ii] /= search.dt;
  }
  
  /* 
   *   Perform the actual folding of the data
   */
  
  proftime = worklen * search.dt;
  parttimes = gen_dvect(cmd->npart);
  printf("Folded %ld points of %.0f", totnumfolded, N);
  
  /* sub-integrations in time  */
  
  dtmp = (double) cmd->npart;
  chkfwrite(&dtmp, sizeof(double), 1, binproffile);
  dtmp = (double) cmd->nsub;
  chkfwrite(&dtmp, sizeof(double), 1, binproffile);
  dtmp = (double) search.proflen;
  chkfwrite(&dtmp, sizeof(double), 1, binproffile);
  for (ii = 0; ii < cmd->npart; ii++){
    parttimes[ii] = ii * reads_per_part * proftime;
    
    /* reads per sub-integration */
    
    for (jj = 0; jj < reads_per_part; jj++){
      numread = readrec_ptr(infile, data, worklen, dispdts, 
			    cmd->nsub, numchan);
     
      /* frequency sub-bands */
      
      for (kk = 0; kk < cmd->nsub; kk++)
	fold(data + kk * worklen, numread, search.dt, 
	     parttimes[ii] + jj * proftime, 
	     search.rawfolds + (ii * cmd->nsub + kk) * search.proflen, 
	     search.proflen, cmd->phs, foldf, foldfd, foldfdd, 
	     flags, Ep, tp, numdelays, NULL, 
	     &(search.stats[ii * cmd->nsub + kk]));
      totnumfolded += numread;
    }

    /* Write the binary profiles */

    for (kk = 0; kk < cmd->nsub; kk++){
      chkfwrite(&(search.stats[ii * cmd->nsub + kk]), 
		sizeof(foldstats), 1, binproffile);
      chkfwrite(search.rawfolds + (ii * cmd->nsub + kk) * 
		search.proflen, sizeof(double), search.proflen, 
		binproffile);
    }
    printf("\rFolded %ld points of %.0f", totnumfolded, N);
  }
  fclose(infile);
  fclose(binproffile);

  /*
   *   Perform the candidate optimization search
   */

  printf("\nOptimizing...\n\n");
  bestprof = gen_dvect(search.proflen);
  {
    int numtrials, pdelay, pddelay, profindex;
    double dphase, po, pdo, pddo, pofact;
    double *pdprofs, *currentprof;
    foldstats currentstats;

    search.ndmfact = cmd->ndmfact;
    search.npfact = cmd->npfact;
    search.pstep = cmd->pstep;
    search.pdstep = cmd->pdstep;
    search.dmstep = cmd->dmstep;

    /* The number of trials for the P-dot and P searches */

    numtrials = 2 * search.npfact * search.proflen + 1;

    /* Initialize a bunch of variables */

    search.periods = gen_dvect(numtrials);
    search.numperiods = numtrials;
    search.pdots = gen_dvect(numtrials);
    search.numpdots = numtrials;
    pdprofs = gen_dvect(cmd->npart * search.proflen);
    currentprof = gen_dvect(search.proflen);
    initialize_foldstats(&beststats);

    /* Convert the folding freqs and derivs to periods */
    
    switch_f_and_p(foldf, foldfd, foldfdd, &po, &pdo, &pddo);
    
    /* Our P and P-dot steps are the changes that cause the pulse      */
    /* to be delayed a number of bins between the first and last time. */
    
    dphase = po / search.proflen;
    pofact = foldf * foldf;
    for (ii = 0; ii < numtrials; ii++){
      pdelay = ii - (numtrials - 1) / 2;
      dtmp = pdelay * dphase * search.pstep;
      search.periods[ii] = 1.0 / (foldf + dtmp / T);
      dtmp = pdelay * dphase * search.pdstep;
      search.pdots[ii] = -((2.0 * dtmp / (T * T) + foldfd) / pofact);
      if (search.pdots[ii]==-0) search.pdots[ii] = 0.0;
    }

    /* If we are searching through DM space */
    
    if (cmd->nsub > 1){
      int numdmtrials, *dmdelays;
      double lodm, ddm, hifdelay, *ddprofs, *subbanddelays;
      foldstats *ddstats;
      
      dmdelays = gen_ivect(cmd->nsub);
      ddprofs = gen_dvect(cmd->npart * search.proflen);
      ddstats = (foldstats *)malloc(cmd->npart * sizeof(foldstats));
      
      /* Our DM step is the change in DM that would cause the pulse   */
      /* to be delayed a number of phasebins at the lowest frequency. */
      
      ddm = dm_from_delay(dphase * search.dmstep, obsf[0]);
    
      /* Insure that we don't try a dm < 0.0 */
      
      numdmtrials = 2 * search.ndmfact * search.proflen + 1;
      lodm = cmd->dm - (numdmtrials - 1) / 2 * ddm;
      if (lodm < 0.0) lodm = 0.0;
      search.dms = gen_dvect(numdmtrials);
      search.numdms = numdmtrials;
      
      printf("Searching %d DMs, %d periods, and %d p-dots...\n\n", 
	     search.numdms, search.numperiods, search.numpdots);

      /* De-disperse and combine the subbands */
      
      for (ii = 0; ii < numdmtrials; ii++){  /* Loop over DMs */
	search.dms[ii] = lodm + ii * ddm;
	hifdelay = delay_from_dm(search.dms[ii], obsf[numchan - 1]);
	subbanddelays = subband_delays(numchan, cmd->nsub, 
				       search.dms[ii], idata.freq, 
				       idata.chan_wid, search.avgvoverc);
	for (jj = 0; jj < cmd->nsub; jj++)
	  dmdelays[jj] = ((int) ((subbanddelays[jj] - hifdelay) / 
				 dphase + 0.5) % search.proflen);
	free(subbanddelays);
	combine_subbands(search.rawfolds, search.stats, cmd->npart, 
			 cmd->nsub, search.proflen, dmdelays, 
			 ddprofs, ddstats);
	
	/* Perform the P-dot and Period searches */

	for (jj = 0; jj < numtrials; jj++){

	  /* Correct each part for the current pdot */

	  for (kk = 0; kk < cmd->npart; kk++){
	    profindex = kk * search.proflen;
	    pddelay = (int) ((-0.5 * parttimes[kk] * parttimes[kk] * 
			      (search.pdots[jj] * foldf * foldf + 
			       foldfd) / dphase) + 0.5);
	    if (pddelay==-0) pddelay = 0.0;
	    shift_prof(ddprofs + profindex, search.proflen, pddelay, 
		       pdprofs + profindex);
	  }

	  /* Search over the periods */

	  for (kk = 0; kk < numtrials; kk++){
	    pdelay = search.pstep * (kk - (numtrials - 1) / 2);
	    combine_profs(pdprofs, ddstats, cmd->npart, search.proflen, 
			  pdelay, currentprof, &currentstats);
	    if (currentstats.redchi > beststats.redchi){
	      search.bestdm = search.dms[ii];
	      if (idata.bary){
		search.bary.p1 = search.periods[kk];
		search.bary.p2 = search.pdots[jj];
	      } else {
		search.topo.p1 = search.periods[kk];
		search.topo.p2 = search.pdots[jj];
	      }
	      beststats = currentstats;
	      memcpy(bestprof, currentprof, sizeof(double) * 
		     search.proflen);
	    }
	  }
	}
      }
      free(dmdelays);
      free(ddprofs);
      free(ddstats);
      
      /* We are not searching through DM space */
      
    } else {

      printf("Searching %d periods, and %d p-dots...\n\n", 
	     search.numperiods, search.numpdots);

      /* Perform the P-dot and Period searches */

      for (jj = 0; jj < numtrials; jj++){
	
	/* Correct each part for the current pdot */
	
	for (kk = 0; kk < cmd->npart; kk++){
	  profindex = kk * search.proflen;
	  pddelay = (int) ((-0.5 * parttimes[kk] * parttimes[kk] * 
			    (search.pdots[jj] * foldf * foldf + 
			     foldfd) / dphase) + 0.5);
	  if (pddelay==-0) pddelay = 0.0;
	  shift_prof(search.rawfolds + profindex, search.proflen, pddelay, 
		     pdprofs + profindex);
	}

	/* Search over the periods */
	
	for (kk = 0; kk < numtrials; kk++){
	  pdelay = search.pstep * (kk - (numtrials - 1) / 2);
	  combine_profs(pdprofs, search.stats, cmd->npart, search.proflen, 
			pdelay, currentprof, &currentstats);
	  if (currentstats.redchi > beststats.redchi){
	    if (idata.bary){
	      search.bary.p1 = search.periods[kk];
	      search.bary.p2 = search.pdots[jj];
	    } else {
	      search.topo.p1 = search.periods[kk];
	      search.topo.p2 = search.pdots[jj];
	    }
	    beststats = currentstats;
	    memcpy(bestprof, currentprof, sizeof(double) * 
		   search.proflen);
	  }
	}
      }
    }
    free(pdprofs);
    free(currentprof);
  }
  printf("Done searching.\n\n");

  {
    double perr, pderr, pdderr;
    char out[100];
  
    printf("Maximum reduced chi-squared found  =  %-.5f\n", 
	   beststats.redchi);
    if (cmd->nsub > 1)
      printf("Best DM  =  %-.4f\n", search.bestdm);
    
    /* Convert best params from/to barycentric to/from topocentric */
    
    if (cmd->nobaryP){

      /* Data was barycentered */

      if (idata.bary){

	/* Calculate the errors in our new pulsation quantities */
	
	fold_errors(bestprof, search.proflen, search.dt, N, 
		    beststats.data_var, search.bary.p1, search.bary.p2, 
		    search.bary.p3, &perr, &pderr, &pdderr);
	
	nice_output_2(out, search.bary.p1, perr, 0);
	printf("Best period        (s)  =  %s\n", out);
	nice_output_2(out, search.bary.p2, pderr, 0);
	printf("Best p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.bary.p3, pdderr, 0);
	printf("Best p-dotdot  (s/s^2)  =  %s\n", out);

      /* Data was topocentric */

      } else {

	/* Calculate the errors in our new pulsation quantities */
	
	fold_errors(bestprof, search.proflen, search.dt, N, 
		    beststats.data_var, search.topo.p1, search.topo.p2, 
		    search.topo.p3, &perr, &pderr, &pdderr);
	
	nice_output_2(out, search.topo.p1, perr, 0);
	printf("Best period        (s)  =  %s\n", out);
	nice_output_2(out, search.topo.p2, pderr, 0);
	printf("Best p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.topo.p3, pdderr, 0);
	printf("Best p-dotdot  (s/s^2)  =  %s\n", out);

      }

    } else {

      if (idata.bary){

	/* Calculate the errors in our new pulsation quantities */
      
	fold_errors(bestprof, search.proflen, search.dt, N, 
		    beststats.data_var, search.bary.p1, search.bary.p2, 
		    search.bary.p3, &perr, &pderr, &pdderr);

	nice_output_2(out, search.bary.p1, perr, 0);
	printf("Best barycentric period        (s)  =  %s\n", out);
	nice_output_2(out, search.bary.p2, pderr, 0);
	printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.bary.p3, pdderr, 0);
	printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);
	
	/* Convert the barycentric folding parameters into topocentric */
	
	switch_f_and_p(search.bary.p1, search.bary.p2, search.bary.p3,
		       &f, &fd, &fdd);
	if ((info=bary2topo(topotimes, barytimes, numbarypts, 
			    f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
	  printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
	switch_f_and_p(foldf, foldfd, foldfdd, 
		       &search.topo.p1, &search.topo.p2, &search.topo.p3);
	
	nice_output_2(out, search.topo.p1, perr, 0);
	printf("Best topocentric period        (s)  =  %s\n", out);
	nice_output_2(out, search.topo.p2, pderr, 0);
	printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.topo.p3, pdderr, 0);
	printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

      } else {
	
	/* Calculate the errors in our new pulsation quantities */
      
	fold_errors(bestprof, search.proflen, search.dt, N, 
		    beststats.data_var, search.topo.p1, search.topo.p2, 
		    search.topo.p3, &perr, &pderr, &pdderr);

	nice_output_2(out, search.topo.p1, perr, 0);
	printf("Best topocentric period        (s)  =  %s\n", out);
	nice_output_2(out, search.topo.p2, pderr, 0);
	printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.topo.p3, pdderr, 0);
	printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

	/* Convert the barycentric folding parameters into topocentric */
	
	switch_f_and_p(search.topo.p1, search.topo.p2, search.topo.p3,
		       &f, &fd, &fdd);
	if ((info=bary2topo(barytimes, topotimes, numbarypts, 
			    f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
	  printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
	switch_f_and_p(foldf, foldfd, foldfdd, 
		       &search.bary.p1, &search.bary.p2, &search.bary.p3);
	
	nice_output_2(out, search.bary.p1, perr, 0);
	printf("Best barycentric period        (s)  =  %s\n", out);
	nice_output_2(out, search.bary.p2, pderr, 0);
	printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
	nice_output_2(out, search.bary.p3, pdderr, 0);
	printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);
      }
    }
  }
  
  printf("\nMaking plots.\n\n");

  /*
   *   Write the raw prepfoldinfo structure
   */

  write_prepfoldinfo(&search, outfilenm);

  /*
   *   Plot our results
   */

  if (cmd->xwinP)
    prepfold_plot(&search, 1);
  else 
    prepfold_plot(&search, 0);

  /* Free our memory  */

  delete_prepfoldinfo(&search);
  free(data);
  free(rootnm);
  free(outfilenm);
  free(plotfilenm);
  free(binproffilenm);
  free(parttimes);
  free(bestprof);
  if (binary){
    free(Ep);
    free(tp);
  }
  if (!cmd->nobaryP){
    free(barytimes);
    free(topotimes);
  }
  if (!strcmp(idata.band, "Radio")){
    free(obsf);
    free(dispdts);
  }
  printf("Done.\n\n");
  return (0);
}
  
