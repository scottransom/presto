#include <ctype.h>
#include "prepfold.h"
#include "prepfold_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP)

extern int getpoly(double mjd, double duration, double *dm, FILE *fp, char *pname);
extern void phcalc(double mjd0, double mjd1, 
		   double *phase, double *psrfreq);
extern int get_psr_from_parfile(char *parfilenm, double epoch, 
				psrparams *psr);
extern char *make_polycos(char *parfilenm, infodata *idata);

/* 
 * The main program 
 */

int main(int argc, char *argv[])
{
  FILE *infiles[MAXPATCHFILES], *filemarker, *binproffile;
  float *data=NULL;
  double f=0.0, fd=0.0, fdd=0.0, foldf=0.0, foldfd=0.0, foldfdd=0.0;
  double recdt=0.0, barydispdt=0.0, N=0.0, T=0.0, proftime, startTday=0.0;
  double polyco_phase=0.0, polyco_phase0=0.0;
  double *obsf=NULL, *dispdts=NULL, *parttimes=NULL, *Ep=NULL, *tp=NULL;
  double *barytimes=NULL, *topotimes=NULL, *bestprof, dtmp;
  double *buffers, *phasesadded, *events=NULL;
  char *plotfilenm, *outfilenm, *rootnm, *binproffilenm, *root, *suffix;
  char obs[3], ephem[6], pname[30], rastring[50], decstring[50];
  int numfiles, numevents, numchan=1, binary=0, numdelays=0, numbarypts=0;
  int info, ptsperrec=1, flags=1, padding=0, arrayoffset=0, useshorts=0;
  int *maskchans=NULL, nummasked=0;
  long ii=0, jj, kk, worklen=0, numread=0, reads_per_part=0;
  long totnumfolded=0, lorec=0, hirec=0, numbinpoints=0, currentrec=0;
  unsigned long numrec=0;
  BPP_ifs bppifs=SUMIFS;
  infodata idata;
  foldstats beststats;
  prepfoldinfo search;
  Cmdline *cmd;
  plotflags pflags;
  mask obsmask;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
    Program = argv[0];
    usage();
    exit(0);
 }
  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);
  obsmask.numchan = obsmask.numint = 0;
  if (cmd->timingP){
    cmd->nopsearchP = 1;
    cmd->nopdsearchP = 1;
    cmd->nodmsearchP = 1;
    cmd->npart = 60;
    cmd->fineP = 1;
  }
  if (cmd->slowP){
    cmd->fineP = 1;
    cmd->proflenP = 1;
    cmd->proflen = 100;
    cmd->nsub = 16;
  }
  if (cmd->fineP){
    cmd->ndmfact = 1;
    cmd->dmstep = 1;
    cmd->npfact = 1;
    cmd->pstep = 1;
    cmd->pdstep = 2;
  }
  if (cmd->coarseP){
    cmd->npfact = 4;
    if (cmd->pstep==1)
      cmd->pstep = 2;
    else
      cmd->pstep = 3;
    if (cmd->pstep==2)
      cmd->pstep = 4;
    else
      cmd->pstep = 6;
  }
  pflags.events = cmd->eventsP;
  pflags.nosearch = cmd->nosearchP;
  pflags.scaleparts = cmd->scalepartsP;
  pflags.justprofs = cmd->justprofsP;
  pflags.allgrey = cmd->allgreyP;

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("        Pulsar Raw-Data Folding Search Routine\n");
  printf(" Used for DM, Period, and P-dot tweaking of PSR candidates.\n");
  printf("                 by Scott M. Ransom\n\n");

  init_prepfoldinfo(&search);

  numfiles = cmd->argc;
  {
    char *path;

    split_path_file(cmd->argv[0], &path, &search.filenm);
    free(path);

    if (!cmd->outfileP){
      char *tmprootnm;

      split_root_suffix(cmd->argv[0], &tmprootnm, &suffix);
      if ((cmd->startT != 0.0) || (cmd->endT != 1.0)){
	rootnm = (char *)calloc(strlen(tmprootnm)+11, sizeof(char));
	sprintf(rootnm, "%s_%4.2f-%4.2f", tmprootnm, cmd->startT, cmd->endT);
      } else {
	rootnm = tmprootnm;
      }
      free(suffix);
    } else {
      rootnm = cmd->outfile;
    }
  }

  if (!RAWDATA){
    /* Split the filename into a rootname and a suffix */
    if (split_root_suffix(cmd->argv[0], &root, &suffix)==0){
      printf("\nThe input filename (%s) must have a suffix!\n\n", 
	     search.filenm);
      exit(1);
    } else {
      if (strcmp(suffix, "sdat")==0){
	useshorts = 1;
      } else if (strcmp(suffix, "bcpm1")==0 || 
		 strcmp(suffix, "bcpm2")==0){
	printf("Assuming the data is from a GBT BCPM...\n");
	cmd->bcpmP = 1;
      } else if (strcmp(suffix, "pkmb")==0){
	printf("Assuming the data is from the Parkes Multibeam system...\n");
	cmd->pkmbP = 1;
      } else if (isdigit(suffix[0]) &&
		 isdigit(suffix[1]) &&
		 isdigit(suffix[2])){
	printf("Assuming the data is from the Arecibo WAPP system...\n");
	cmd->wappP = 1;
      }
    }
    if (RAWDATA){ /* Clean-up a bit */
      free(root);
      free(suffix);
    }
  }

  if (!RAWDATA){
    printf("Reading input data from '%s'.\n", cmd->argv[0]);
    printf("Reading information from '%s.inf'.\n\n", root);
    /* Read the info file if available */
    readinf(&idata, root);
    free(root);
    free(suffix);
    /* Use events instead of a time series */
    if (cmd->eventsP){
      int eventtype=0; /* 0=sec since .inf, 1=days since .inf, 2=MJDs */

      /* The following allows using inf files from searches of a subset */
      /* of events from an event file.                                  */
      if (cmd->rzwcandP || cmd->accelcandP) {
	infodata rzwidata;
	char *cptr=NULL;
	
	if (cmd->rzwcandP) {
	  if (!cmd->rzwfileP) {
	    printf("\nYou must enter a name for the rzw candidate ");
	    printf("file (-rzwfile filename)\n");
	    printf("Exiting.\n\n");
	    exit(1);
	  } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))){
	    ii = (long) (cptr - cmd->rzwfile);
	  } else if (NULL != (cptr = strstr(cmd->rzwfile, "_ACCEL"))){
	    ii = (long) (cptr - cmd->rzwfile);
	  }
	  cptr = (char *)calloc(ii + 1, sizeof(char));
	  strncpy(cptr, cmd->rzwfile, ii);
	}
	if (cmd->accelcandP) {
	  if (!cmd->accelfileP) {
	    printf("\nYou must enter a name for the ACCEL candidate ");
	    printf("file (-accelfile filename)\n");
	    printf("Exiting.\n\n");
	    exit(1);
	  } else if (NULL != (cptr = strstr(cmd->accelfile, "_rzw"))){
	    ii = (long) (cptr - cmd->accelfile);
	  } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))){
	    ii = (long) (cptr - cmd->accelfile);
	  }
	  cptr = (char *)calloc(ii + 1, sizeof(char));
	  strncpy(cptr, cmd->accelfile, ii);
	}
	readinf(&rzwidata, cptr);
	free(cptr);
	idata.mjd_i = rzwidata.mjd_i;
	idata.mjd_f = rzwidata.mjd_f;
	idata.N = rzwidata.N;
	idata.dt = rzwidata.dt;
      }
      printf("Assuming the events are barycentered.\n");
      if (!cmd->proflenP){
	cmd->proflenP = 1;
	cmd->proflen = 20;
	printf("Using %d bins in the profile since not specified.\n", 
	       cmd->proflen);
      }
      if (cmd->doubleP)
	infiles[0] = chkfopen(cmd->argv[0], "rb");
      else
	infiles[0] = chkfopen(cmd->argv[0], "r");
      if (cmd->daysP)
	eventtype = 1;
      else if (cmd->mjdsP)
	eventtype = 2;
      events = read_events(infiles[0], cmd->doubleP, eventtype, &numevents,
			   idata.mjd_i+idata.mjd_f, idata.N*idata.dt, 
			   cmd->startT, cmd->endT, cmd->offset);
      if (cmd->rzwcandP || cmd->accelcandP)
	T = idata.N*idata.dt;
      else
	T = events[numevents-1];
    } else {
      infiles[0] = chkfopen(cmd->argv[0], "rb");
    }
  } else if (cmd->pkmbP){
    if (numfiles > 1)
      printf("Reading Parkes PKMB data from %d files:\n", numfiles);
    else
      printf("Reading Parkes PKMB data from 1 file:\n");
  } else if (cmd->bcpmP){
    if (numfiles > 1)
      printf("Reading Green Bank BCPM data from %d files:\n", numfiles);
    else
      printf("Reading Green Bank BCPM data from 1 file:\n");
  } else if (cmd->wappP){
    if (numfiles > 1)
      printf("Reading Arecibo WAPP data from %d files:\n", numfiles);
    else
      printf("Reading Arecibo WAPP data from 1 file:\n");
  }

  /* Open the raw data files */

  if (RAWDATA){
    for (ii=0; ii<numfiles; ii++){
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }
    printf("\n");
    if (cmd->maskfileP){
      read_mask(cmd->maskfile, &obsmask);
      printf("Read mask information from '%s'\n", cmd->maskfile);
    } else {
      obsmask.numchan = obsmask.numint = 0;
    }
  }

  /* Manipulate the file names we will use  */

  {
    int slen;

    /* Determine the candidate name */
    
    if (cmd->psrnameP){
      search.candnm = (char *)calloc(strlen(cmd->psrname)+5, sizeof(char));
      sprintf(search.candnm, "PSR_%s", cmd->psrname);
    } else if (cmd->parnameP || cmd->timingP){
      int retval;
      psrparams psr;
      if (cmd->timingP){
	cmd->parnameP = 1;
	cmd->parname = cmd->timing;
      }
      /* Read the par file just to get the PSR name */
      retval = get_psr_from_parfile(cmd->parname, 51000.0, &psr);
      search.candnm = (char *)calloc(strlen(psr.jname)+5, sizeof(char));
      sprintf(search.candnm, "PSR_%s", psr.jname);
    } else if (cmd->rzwcandP) {						
      slen = 20;
      search.candnm = (char *)calloc(slen, sizeof(char));
      sprintf(search.candnm, "RZW_Cand_%d", cmd->rzwcand);
    } else if (cmd->accelcandP) {					       
      slen = 22;
      search.candnm = (char *)calloc(slen, sizeof(char));
      sprintf(search.candnm, "ACCEL_Cand_%d", cmd->accelcand);
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
  }

  /* What ephemeris will we use?  (Default is DE200) */
  
  if (cmd->de405P)
    strcpy(ephem, "DE405");
  else
    strcpy(ephem, "DE200");
    
  /* Set-up values if we are using the Parkes Multibeam System, */
  /* a Berkeley-Caltech Pulsar Machine, or the Arecibo WAPP.    */

  if (RAWDATA){
    double local_dt, local_T;
    long long local_N;

    if (cmd->pkmbP){
      PKMB_tapehdr hdr;
      
      printf("PKMB input file information:\n");
      get_PKMB_file_info(infiles, numfiles, &local_N, &ptsperrec, 
			 &numchan, &local_dt, &local_T, 1);

      /* Read the first header file and generate an infofile from it */
      
      chkfread(&hdr, 1, HDRLEN, infiles[0]);
      rewind(infiles[0]);
      PKMB_hdr_to_inf(&hdr, &idata);
      PKMB_update_infodata(numfiles, &idata);

      /* OBS code for TEMPO */
      strcpy(obs, "PK");
      search.telescope = (char *)calloc(20, sizeof(char));
      strcpy(search.telescope, "Parkes Multibeam");

    } else if (cmd->bcpmP){

      printf("BCPM input file information:\n");
      get_BPP_file_info(infiles, numfiles, &local_N, &ptsperrec, &numchan, 
			&local_dt, &local_T, &idata, 1);
      BPP_update_infodata(numfiles, &idata);

      /* Which IFs will we use? */
      
      if (cmd->ifsP){
	if (cmd->ifs==0)
	  bppifs = IF0;
	else if (cmd->ifs==1)
	  bppifs = IF1;
	else
	  bppifs = SUMIFS;
      }
      
      /* OBS code for TEMPO */
      
      /* The following is for the Green Bank 85-3
      strcpy(obs, "G8");
      search.telescope = (char *)calloc(20, sizeof(char));
      strcpy(search.telescope, "GBT");
      */

      /* The following is for the Green Bank Telescope */
      strcpy(obs, "GB");
      search.telescope = (char *)calloc(20, sizeof(char));
      strcpy(search.telescope, "GBT");

    } else if (cmd->wappP){

      printf("WAPP input file information:\n");
      get_WAPP_file_info(infiles, cmd->numwapps, numfiles, cmd->clip,
			 &local_N, &ptsperrec, &numchan, 
			 &local_dt, &local_T, &idata, 1);
      WAPP_update_infodata(numfiles, &idata);

      /* OBS code for TEMPO */
      
      /* The following is for Arecibo */
      strcpy(obs, "AO");
      search.telescope = (char *)calloc(20, sizeof(char));
      strcpy(search.telescope, "Arecibo");
    }

    idata.dm = cmd->dm;
    numrec = local_N / ptsperrec;
    if (cmd->maskfileP)
      maskchans = gen_ivect(idata.num_chan);

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
    startTday = lorec*recdt/SECPERDAY;
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
      search.tepoch = idata.mjd_i + idata.mjd_f + startTday;

      if (!cmd->polycofileP && !cmd->timingP && !cmd->parnameP) {
	barycenter(&search.tepoch, &search.bepoch, &dtmp, 1, rastring,
		   decstring, obs, ephem);
	
	/* Correct the barycentric time for the dispersion delay.     */
	/* This converts the barycentric time to infinite frequency.  */
	if (cmd->dm > 0.0){
	  barydispdt = delay_from_dm(cmd->dm, idata.freq + 
				     (idata.num_chan - 1) * idata.chan_wid);
	  search.bepoch -= (barydispdt / SECPERDAY);
	}
      }
    }
    worklen = ptsperrec;

  } else { /* Raw floating point or event data (already de-dispersed if radio data) */

    cmd->nsub = 1;
    search.startT = cmd->startT;
    search.endT = cmd->endT;
    
    if (!cmd->eventsP){
      
      /* Some information about the size of the records */
      
      numchan = 1;
      worklen = 1024;
      search.dt = idata.dt;
      if (useshorts)
	N = chkfilelen(infiles[0], sizeof(short));
      else
	N = chkfilelen(infiles[0], sizeof(float));

      /* Determine the number of records to use from the command line */
      
      lorec = (long) (cmd->startT * N + DBLCORRECT);
      hirec = (long) (cmd->endT * N + DBLCORRECT);
      startTday = lorec * search.dt / SECPERDAY;
      numrec = (hirec - lorec) / worklen;
      recdt = worklen * search.dt;

      /* The number of reads from the file we need for */
      /* each sub-integration.                         */ 
      
      reads_per_part = numrec / cmd->npart;
      
      /* Correct numrec so that each part will contain */
      /* the same number of records.                   */
      
      numrec = reads_per_part * cmd->npart;
      N = numrec * worklen;
      T = N * search.dt;
    }

    /* Until I figure out a better way to do this... */
      
    search.telescope = (char *)calloc(strlen(idata.telescope)+1, 
				      sizeof(char));
    strcpy(search.telescope, idata.telescope);

    if (idata.mjd_i){
      if (idata.bary)
	search.bepoch = idata.mjd_i + idata.mjd_f + startTday;
      else
	search.tepoch = idata.mjd_i + idata.mjd_f + startTday;
    }
  }

  printf("Folding a %s candidate.\n\n", search.candnm);
  printf("Output data file is '%s'.\n", outfilenm);
  printf("Output plot file is '%s'.\n", plotfilenm);
  printf("Raw profile file is '%s'.\n", binproffilenm);
  printf("Best profile is in  '%s.bestprof'.\n", outfilenm);
  
  /* Generate polycos if required and set the pulsar name */
  if ((cmd->timingP || cmd->parnameP) && !idata.bary){
    char *polycofilenm;
    cmd->psrnameP = 1;
    if (cmd->timingP)
      cmd->psrname = make_polycos(cmd->timing, &idata);
    else
      cmd->psrname = make_polycos(cmd->parname, &idata);
    polycofilenm = (char *)calloc(strlen(outfilenm)+9, sizeof(char));
    sprintf(polycofilenm, "%s.polycos", outfilenm);
    printf("Polycos used are in '%s'.\n", polycofilenm);
    rename("polyco.dat", polycofilenm);
    cmd->polycofileP = 1;
    cmd->polycofile = (char *)calloc(strlen(polycofilenm)+1, sizeof(char));
    strcpy(cmd->polycofile, polycofilenm);
    free(polycofilenm);
  }

  /* Read the pulsar database if needed */
  if (cmd->psrnameP) {
    if (cmd->polycofileP) {
      FILE *polycofileptr;
      int numsets;
      double polyco_dm;
      
      polycofileptr = chkfopen(cmd->polycofile, "r");
      numsets = getpoly(search.tepoch, T/SECPERDAY, &polyco_dm, 
			polycofileptr, cmd->psrname);
      fclose(polycofileptr);
      if (cmd->dm > 0.0){
	printf("\nRead %d set(s) of polycos for PSR %s at %18.12f\n", 
	       numsets, cmd->psrname, search.tepoch);
	printf("Overriding polyco DM = %f with %f\n", 
	       polyco_dm, cmd->dm);
      } else {
	printf("\nRead %d set(s) of polycos for PSR %s at %18.12f (DM = %.5g)\n", 
	       numsets, cmd->psrname, search.tepoch, polyco_dm);
	cmd->dm = polyco_dm;
      }
      phcalc(idata.mjd_i, idata.mjd_f+startTday, &polyco_phase0, &f);
      search.topo.p1 = 1.0/f;
      search.topo.p2 = fd = 0.0;
      search.topo.p3 = fdd = 0.0;
      strcpy(pname, cmd->psrname);
    } else {  /* Use the database */
      int pnum;
      psrparams psr;
      
      if (search.bepoch==0.0){
	printf("\nYou cannot fold topocentric data with the pulsar database.\n");
	printf("Use '-timing' or polycos instead.  Exiting.\n\n");
	exit(1);
      }
      pnum = get_psr_at_epoch(cmd->psrname, search.bepoch, &psr);
      if (!pnum) {
	printf("The pulsar is not in the database.  Exiting.\n\n");
	exit(1);
      }
      if (psr.orb.p != 0.0){  /* Checks if the pulsar is in a binary */
	binary = 1;
	search.orb = psr.orb;
      }
      search.bary.p1 = psr.p;
      search.bary.p2 = psr.pd;
      search.bary.p3 = psr.pdd;
      if (cmd->dm==0.0)
	cmd->dm = psr.dm;
      f = psr.f;
      fd = psr.fd;
      fdd = psr.fdd;
      strcpy(pname, psr.jname);
    }
    
  } else if (cmd->parnameP) { /* Read ephemeris from a par file */
    int retval;
    psrparams psr;
    
    if (search.bepoch==0.0){
      printf("\nYou cannot fold topocentric data with a par file.\n");
      printf("Use '-timing' or polycos instead.  Exiting.\n\n");
      exit(1);
    }

    /* Read the par file*/
    retval = get_psr_from_parfile(cmd->parname, search.bepoch, &psr);

    if (psr.orb.p != 0.0){  /* Checks if the pulsar is in a binary */
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
    if (cmd->dm==0.0)
      cmd->dm = psr.dm;

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
    search.orb.w = (cmd->w + dtmp*cmd->wdot/SECPERJULYR);
    binary = 1;

  } else if (cmd->rzwcandP || cmd->accelcandP) {
    fourierprops rzwcand;
    infodata rzwidata;
    char *cptr=NULL;

    if (cmd->rzwcandP) {
      if (!cmd->rzwfileP) {
	printf("\nYou must enter a name for the rzw candidate ");
	printf("file (-rzwfile filename)\n");
	printf("Exiting.\n\n");
	exit(1);
      } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))){
	ii = (long) (cptr - cmd->rzwfile);
      } else if (NULL != (cptr = strstr(cmd->rzwfile, "_ACCEL"))){
	ii = (long) (cptr - cmd->rzwfile);
      }
      cptr = (char *)calloc(ii + 1, sizeof(char));
      strncpy(cptr, cmd->rzwfile, ii);
    }
    if (cmd->accelcandP) {
      if (!cmd->accelfileP) {
	printf("\nYou must enter a name for the ACCEL candidate ");
	printf("file (-accelfile filename)\n");
	printf("Exiting.\n\n");
	exit(1);
      } else if (NULL != (cptr = strstr(cmd->accelfile, "_rzw"))){
	ii = (long) (cptr - cmd->accelfile);
      } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))){
	ii = (long) (cptr - cmd->accelfile);
      }
      cptr = (char *)calloc(ii + 1, sizeof(char));
      strncpy(cptr, cmd->accelfile, ii);
    }
    fprintf(stderr, "\nAttempting to read '%s.inf'.  ", cptr);
    readinf(&rzwidata, cptr);
    free(cptr);
    fprintf(stderr, "Successful.\n");
    if (cmd->rzwfileP)
      get_rzw_cand(cmd->rzwfile, cmd->rzwcand, &rzwcand);	
    if (cmd->accelfileP)
      get_rzw_cand(cmd->accelfile, cmd->accelcand, &rzwcand);	
    f = (rzwcand.r - 0.5 * rzwcand.z) / 
      (rzwidata.dt * rzwidata.N);
    fd = rzwcand.z / ((rzwidata.dt * rzwidata.N) * 
		      (rzwidata.dt * rzwidata.N));
    
    /* Now correct for the fact that we may not be starting */
    /* to fold at the same start time as the rzw search.    */
    
    if (RAWDATA)
      f += lorec * recdt * fd;
    else
      f += lorec * search.dt * fd;
    if (rzwidata.bary)
      switch_f_and_p(f, fd, fdd, &search.bary.p1, \
		     &search.bary.p2, &search.bary.p3);
    else
      switch_f_and_p(f, fd, fdd, &search.topo.p1, \
		     &search.topo.p2, &search.topo.p3);
  }
  
  /* Determine the pulsar parameters to fold if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */
  
  if (!cmd->rzwcandP && !cmd->accelcandP && 
      !cmd->psrnameP && !cmd->parnameP) {
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
  
  /* Modify the fold period or frequency */

  if (cmd->pfact != 1.0 ||
      cmd->ffact != 1.0){
    if (cmd->pfact==0.0 || cmd->ffact==0.0){
      printf("\nFoldint factors (-pfact or -ffact) cannot be 0!  Exiting\n");
      exit(1);
    }
    if (cmd->pfact != 1.0)
      cmd->ffact = 1.0 / cmd->pfact;
    else if (cmd->ffact != 1.0)
      cmd->pfact = 1.0 / cmd->ffact;
    search.bary.p1 *= cmd->pfact;
    search.bary.p2 /= cmd->pfact;
    search.topo.p1 *= cmd->pfact;
    search.topo.p2 /= cmd->pfact;
    f *= cmd->ffact;
    fd *= cmd->ffact;
    fdd *= cmd->ffact;
  }

  /* Determine the length of the profile */				
  
  if (cmd->proflenP) {
    search.proflen = cmd->proflen;
  } else {
    if (search.topo.p1 == 0.0)
      search.proflen = (long) (search.bary.p1 / search.dt + 0.5);
    else
      search.proflen = (long) (search.topo.p1 / search.dt + 0.5);
    if (cmd->timingP)
      search.proflen = next2_to_n(search.proflen);
    if (search.proflen > 64)
      search.proflen = 64;
  }

  /* Determine the phase delays caused by the orbit if needed */

  if (binary && !cmd->eventsP) {
    double orbdt = 1.0, startE=0.0;

    /* Save the orbital solution every half second               */
    /* The times in *tp are now calculated as barycentric times. */
    /* Later, we will change them to topocentric times after     */
    /* applying corrections to Ep using TEMPO.                   */
    
    startE = keplars_eqn(search.orb.t, search.orb.p, 
			 search.orb.e, 1.0E-15);
    if (T > 2048) orbdt = 0.5;
    else orbdt = T/4096.0;
    numbinpoints = (long) floor(T/orbdt + 0.5) + 1;
    Ep = dorbint(startE, numbinpoints, orbdt, &search.orb);
    tp = gen_dvect(numbinpoints);
    for (ii=0; ii<numbinpoints; ii++) tp[ii] = ii*orbdt;

    /* Convert Eccentric anomaly to time delays */			
    
    search.orb.w *= DEGTORAD;
    E_to_phib(Ep, numbinpoints, &search.orb);
    numdelays = numbinpoints;
    if (search.bepoch == 0.0)
      search.orb.t = -search.orb.t/SECPERDAY + search.tepoch;
    else 
      search.orb.t = -search.orb.t/SECPERDAY + search.bepoch;
  }

  if (RAWDATA && cmd->dm==0.0 && !cmd->polycofileP){
    /* Correct the barycentric time for the dispersion delay.     */
    /* This converts the barycentric time to infinite frequency.  */
    barydispdt = delay_from_dm(cmd->dm, idata.freq + 
			       (idata.num_chan - 1) * idata.chan_wid);
    search.bepoch -= (barydispdt / SECPERDAY);
  }

  if (cmd->eventsP){
    search.dt = (search.bary.p1 + 0.5*T*search.bary.p2)/search.proflen;
    N = ceil(T/search.dt);
  }

  /* Output some informational data on the screen and to the */
  /* output file.                                            */

  fprintf(stdout, "\n");
  filemarker = stdout;
  for (ii=0; ii<1; ii++){
    double p, pd, pdd;

    if (cmd->psrnameP)
      fprintf(filemarker, 
	      "Pulsar                       =  %s\n", pname);
    if (search.tepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (topo) epoch  (MJD)  =  %-17.12f\n", search.tepoch);
    if (search.bepoch != 0.0)
      fprintf(filemarker, 
	      "Folding (bary) epoch  (MJD)  =  %-17.12f\n", search.bepoch);
    fprintf(filemarker, 
	    "Data pt duration (dt)   (s)  =  %-.12g\n", search.dt);
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
      double tmpTo;
      fprintf(filemarker, 
	      "Orbital period          (s)  =  %-.10g\n", search.orb.p);
      fprintf(filemarker, 
	      "a*sin(i)/c (x)     (lt-sec)  =  %-.10g\n", search.orb.x);
      fprintf(filemarker, 
	      "Eccentricity                 =  %-.10g\n", search.orb.e);
      fprintf(filemarker, 
	      "Longitude of peri (w) (deg)  =  %-.10g\n", 
	      search.orb.w / DEGTORAD); 
      tmpTo = search.orb.t;
      if (cmd->eventsP){
	if (search.bepoch == 0.0)
	  tmpTo = -search.orb.t/SECPERDAY + search.tepoch;
	else 
	  tmpTo = -search.orb.t/SECPERDAY + search.bepoch;
      }
      fprintf(filemarker, 
	      "Epoch of periapsis    (MJD)  =  %-17.11f\n", tmpTo);
    }
  }

  /* Allocate and initialize some arrays and other information */
  
  search.nsub = cmd->nsub;
  search.npart = cmd->npart;
  search.rawfolds = gen_dvect(cmd->nsub * cmd->npart * search.proflen);
  search.stats = (foldstats *)malloc(sizeof(foldstats) * 
				     cmd->nsub * cmd->npart);
  for (ii = 0; ii < cmd->npart * cmd->nsub * search.proflen; ii++)
    search.rawfolds[ii] = 0.0;
  for (ii = 0; ii < cmd->npart * cmd->nsub; ii++){
    search.stats[ii].numdata = 0.0;
    search.stats[ii].data_avg = 0.0;
    search.stats[ii].data_var = 0.0;
  }
  if (numdelays == 0) flags = 0;

  if (cmd->eventsP){  /* Fold events instead of a time series */
    double event, dtmp, cts, phase;
    double tf, tfd, tfdd;
    int partnum, binnum;
    
    binproffile = chkfopen(binproffilenm, "wb");
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    dtmp = (double) cmd->nsub;
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    dtmp = (double) search.proflen;
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    foldf = f;  foldfd = fd;  foldfdd = fdd;
    search.fold.pow = 1.0;
    search.fold.p1 = f;
    search.fold.p2 = fd;
    search.fold.p3 = fdd;
    tf = f;
    tfd = fd/2.0;
    tfdd = fdd/6.0;
    dtmp = cmd->npart/T;
    parttimes = gen_dvect(cmd->npart);
    if (binary) search.orb.w *= DEGTORAD;
    for (ii=0; ii<numevents; ii++){
      event = events[ii];
      if (binary){
	double tt, delay;
	tt = fmod(search.orb.t+event, search.orb.p);
	delay = keplars_eqn(tt, search.orb.p, search.orb.e, 1.0E-15);
	E_to_phib(&delay, 1, &search.orb);
	event -= delay;
      }
      partnum = (int) floor(event*dtmp);
      phase = event*(event*(event*tfdd+tfd)+tf);
      binnum = (int)((phase-(int)phase)*search.proflen);
      search.rawfolds[partnum*search.proflen + binnum] += 1.0;
    }
    if (binary){
      if (search.bepoch == 0.0)
	search.orb.t = -search.orb.t/SECPERDAY + search.tepoch;
      else 
	search.orb.t = -search.orb.t/SECPERDAY + search.bepoch;
    }
    for (ii=0; ii<cmd->npart; ii++){
      parttimes[ii] = (T*ii)/cmd->npart;
      cts = 0.0;
      for (jj=ii*search.proflen; jj<(ii+1)*search.proflen; jj++)
	cts += search.rawfolds[jj];
      search.stats[ii].numdata = ceil((T/cmd->npart)/search.dt);
      search.stats[ii].numprof = search.proflen;
      search.stats[ii].prof_avg = search.stats[ii].prof_var = \
	cts/search.proflen;
      search.stats[ii].data_avg = search.stats[ii].data_var = \
	search.stats[ii].prof_avg/search.stats[ii].numdata;
      /* Compute the Chi-Squared probability that there is a signal */
      /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */
      search.stats[ii].redchi = 0.0;
      for (jj=ii*search.proflen; jj<(ii+1)*search.proflen; jj++){
	dtmp = search.rawfolds[jj] - search.stats[ii].prof_avg;
	search.stats[ii].redchi += dtmp * dtmp;
      }
      search.stats[ii].redchi /= (search.stats[ii].prof_var * 
				  (search.proflen - 1));
      chkfwrite(search.stats + ii, sizeof(foldstats), 1, binproffile);
      chkfwrite(search.rawfolds + ii*search.proflen, sizeof(double), 
		search.proflen, binproffile);
    }
    printf("\r  Folded %d events.", numevents);
    fflush(NULL);

  } else {
   
    buffers = gen_dvect(cmd->nsub*search.proflen);
    phasesadded = gen_dvect(cmd->nsub);
    for (ii=0; ii<cmd->nsub*search.proflen; ii++)
      buffers[ii] = 0.0;
    for (ii=0; ii<cmd->nsub; ii++)
      phasesadded[ii] = 0.0;

    /* Move to the correct starting record */

    data = gen_fvect(cmd->nsub * worklen);
    currentrec = lorec;
    if (cmd->pkmbP)
      skip_to_PKMB_rec(infiles, numfiles, lorec+1);
    else if (cmd->bcpmP)
      skip_to_BPP_rec(infiles, numfiles, lorec+1);
    else if (cmd->wappP)
      skip_to_WAPP_rec(infiles, numfiles, lorec+1);
    else {
      if (useshorts)
	chkfileseek(infiles[0], lorec, sizeof(short), SEEK_SET);
      else
	chkfileseek(infiles[0], lorec, sizeof(float), SEEK_SET);
    }
    
    /* Data is already barycentered or the polycos will take care of the barycentering*/
    if (!RAWDATA || cmd->polycofileP){
      foldf = f;  foldfd = fd;  foldfdd = fdd;
    } else { /* Correct our fold parameters if we are barycentering */
      double *voverc;
    
      /* The number of topo to bary points to generate with TEMPO */
    
      numbarypts = T/TDT + 10;
      barytimes = gen_dvect(numbarypts);
      topotimes = gen_dvect(numbarypts);
      voverc = gen_dvect(numbarypts);
    
      /* topocentric times in days from data start */
    
      for (ii=0; ii<numbarypts; ii++)
	topotimes[ii] = search.tepoch + (double) ii*TDT/SECPERDAY;
    
      /* Call TEMPO for the barycentering */
    
      printf("Generating barycentric corrections...\n");
      barycenter(topotimes, barytimes, voverc, numbarypts, \
		 rastring, decstring, obs, ephem);
    
      /* Determine the avg v/c of the Earth's motion during the obs */
    
      for (ii=0; ii<numbarypts-1; ii++)
	search.avgvoverc += voverc[ii];
      search.avgvoverc /= (numbarypts-1.0);
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
	for (ii=0; ii<numbinpoints; ii++){
	  arrayoffset++;  /* Beware nasty NR zero-offset kludges! */
	  dtmp = search.bepoch + tp[ii]/SECPERDAY;
	  hunt(barytimes, numbarypts, dtmp, &arrayoffset);
	  arrayoffset--;
	  tp[ii] = LININTERP(dtmp, barytimes[arrayoffset], 
			     barytimes[arrayoffset+1], 
			     topotimes[arrayoffset], 
			     topotimes[arrayoffset+1]);
	}
	numdelays = numbinpoints;
	dtmp = tp[0];
	for (ii=0 ; ii<numdelays ; ii++)
	  tp[ii] = (tp[ii]-dtmp)*SECPERDAY;
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
      if (RAWDATA){
	for (ii = 0; ii < numchan; ii++)
	  obsf[ii] = doppler(obsf[ii], search.avgvoverc);
      } 
      dispdts = subband_search_delays(numchan, cmd->nsub, cmd->dm,
				      idata.freq, idata.chan_wid, 
				      search.avgvoverc); 
    
      /* Convert the delays in seconds to delays in bins */
    
      for (ii = 0; ii < numchan; ii++)
	dispdts[ii] /= search.dt;

      if (cmd->nsub > 1 && RAWDATA){
	int numdmtrials;
	double dphase, lodm, hidm, ddm;
      
	dphase = 1/(foldf*search.proflen);
	ddm = dm_from_delay(dphase*cmd->dmstep, obsf[0]);
	numdmtrials = 2*cmd->ndmfact*search.proflen + 1;
	lodm = cmd->dm - (numdmtrials-1)/2*ddm;
	if (lodm < 0.0) lodm = 0.0;
	hidm = lodm + numdmtrials * ddm;
	printf("\nWill search %d DMs from %.3f to %.3f (ddm = %.4f)\n", 
	       numdmtrials, lodm, hidm, ddm);
      }
    }
  
    /* 
     *   Perform the actual folding of the data
     */
  
    printf("\nStarting work on '%s'...\n\n", search.filenm);
    proftime = worklen * search.dt;
    parttimes = gen_dvect(cmd->npart);
    printf("  Folded %ld points of %.0f", totnumfolded, N);
  
    /* sub-integrations in time  */
  
    dtmp = (double) cmd->npart;
    binproffile = chkfopen(binproffilenm, "wb");
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    dtmp = (double) cmd->nsub;
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    dtmp = (double) search.proflen;
    chkfwrite(&dtmp, sizeof(double), 1, binproffile);
    for (ii = 0; ii < cmd->npart; ii++){
      parttimes[ii] = ii * reads_per_part * proftime;
    
      /* reads per sub-integration */
    
      for (jj = 0; jj < reads_per_part; jj++){
	double fold_time0;

	if (cmd->pkmbP)
	  numread = read_PKMB_subbands(infiles, numfiles, data, 
				       dispdts, cmd->nsub, 1, &padding,
				       maskchans, &nummasked, &obsmask);
	else if (cmd->bcpmP)
	  numread = read_BPP_subbands(infiles, numfiles, data, 
				      dispdts, cmd->nsub, 1, &padding,
				      maskchans, &nummasked, &obsmask, 
				      bppifs);
	else if (cmd->wappP)
	  numread = read_WAPP_subbands(infiles, numfiles, data, 
				       dispdts, cmd->nsub, 1, &padding,
				       maskchans, &nummasked, &obsmask);
	else {
	  int mm;
	  float runavg=0.0;
	
	  if (useshorts)
	    numread = read_shorts(infiles[0], data, worklen, numchan);
	  else
	    numread = read_floats(infiles[0], data, worklen, numchan);
	  if (cmd->runavgP){
	    for (mm=0; mm<numread; mm++)
	      runavg += data[mm];
	    runavg /= numread;
	    for (mm=0; mm<numread; mm++)
	      data[mm] -= runavg;
	  }
	}
      
	if (cmd->polycofileP){  /* Update the period/phase */
	  double mjdf, currentsec, currentday, offsetphase, orig_cmd_phs=0.0;

	  if (ii==0 && jj==0) orig_cmd_phs = cmd->phs;
	  currentsec = parttimes[ii] + jj*proftime;
	  currentday = currentsec/SECPERDAY;
	  mjdf = idata.mjd_f + startTday + currentday;
	  /* Calculate the pulse phase at the start of the current block */
	  phcalc(idata.mjd_i, mjdf, &polyco_phase, &foldf);
	  polyco_phase -= polyco_phase0;
	  if (polyco_phase < 0.0) polyco_phase += 1.0;
	  /* Calculate the folding frequency at the middle of the current block */
	  phcalc(idata.mjd_i, mjdf+0.5*proftime/SECPERDAY, &offsetphase, &foldf);
	  cmd->phs = orig_cmd_phs + polyco_phase;
	  fold_time0 = 0.0;
	} else {
	  fold_time0 = parttimes[ii]+jj*proftime;
	}

	/* Fold the frequency sub-bands */
	for (kk=0; kk<cmd->nsub; kk++)
	  fold(data+kk*worklen, numread, search.dt, 
	       fold_time0, search.rawfolds+(ii*cmd->nsub+kk)*search.proflen, 
	       search.proflen, cmd->phs, buffers+kk*search.proflen, 
	       phasesadded+kk, foldf, foldfd, foldfdd, 
	       flags, Ep, tp, numdelays, NULL, 
	       &(search.stats[ii*cmd->nsub+kk]));
	totnumfolded += numread;

      }

      /* Write the binary profiles */
      
      for (kk=0; kk<cmd->nsub; kk++){
	chkfwrite(&(search.stats[ii * cmd->nsub + kk]), 
		  sizeof(foldstats), 1, binproffile);
	chkfwrite(search.rawfolds + (ii*cmd->nsub + kk) * 
		  search.proflen, sizeof(double), search.proflen, 
		  binproffile);
      }
      printf("\r  Folded %ld points of %.0f", totnumfolded, N);
      fflush(NULL);
    }
    free(buffers);
    free(phasesadded);
  }
  for (ii=0; ii<numfiles; ii++)
    fclose(infiles[ii]);
  fclose(binproffile);
  
  /*
   *   Perform the candidate optimization search
   */
  
  printf("\n\nOptimizing...\n\n");
  bestprof = gen_dvect(search.proflen);
  {
    int ll, numtrials, pdelay, pddelay, profindex;
    int good_dm_index=0, good_p_index=0, good_pd_index=0, good_pdd_index=0;
    double dphase, po, pdo, pddo;
    double *pdprofs, *pddprofs=NULL, *currentprof, *fdots, *fdotdots=NULL;
    foldstats currentstats;

    search.ndmfact = cmd->ndmfact;
    search.npfact = cmd->npfact;
    search.pstep = cmd->pstep;
    search.pdstep = cmd->pdstep;
    search.dmstep = cmd->dmstep;

    /* The number of trials for the P-dot and P searches */

    numtrials = 2*search.npfact*search.proflen + 1;

    /* Initialize a bunch of variables */

    search.numperiods = numtrials;
    search.periods = gen_dvect(numtrials);
    search.pdots = gen_dvect(numtrials);
    fdots = gen_dvect(numtrials);
    if (cmd->searchfddP)
      cmd->searchpddP = 1;
    if (cmd->nopsearchP && cmd->nopdsearchP){
      if (cmd->nsub > 1 && cmd->nodmsearchP)
	cmd->nosearchP = 1;
      else if (cmd->nsub==1)
	cmd->nosearchP = 1;
      pflags.nosearch = cmd->nosearchP;
    }
    if (cmd->nosearchP){
      cmd->nopsearchP = cmd->nopdsearchP = 1;
      if (cmd->nsub > 1)
	cmd->nodmsearchP = 1;
    }
    if (cmd->searchpddP)
      fdotdots = gen_dvect(numtrials);
    search.numpdots = numtrials;
    pdprofs = gen_dvect(cmd->npart * search.proflen);
    if (cmd->searchpddP)
      pddprofs = gen_dvect(cmd->npart * search.proflen);
    currentprof = gen_dvect(search.proflen);
    initialize_foldstats(&beststats);
    if (cmd->nopsearchP)
      good_p_index = search.npfact*search.proflen;
    if (cmd->nopdsearchP)
      good_pd_index = search.npfact*search.proflen;
    if (cmd->nosearchP && cmd->searchpddP)
      good_pdd_index = search.npfact*search.proflen;

    /* Convert the folding freqs and derivs to periods */
    
    switch_f_and_p(foldf, foldfd, foldfdd, &po, &pdo, &pddo);
    
    /* Our P and P-dot steps are the changes that cause the pulse      */
    /* to be delayed a number of bins between the first and last time. */
    
    dphase = po / search.proflen;
    for (ii = 0; ii < numtrials; ii++){
      pdelay = ii - (numtrials - 1) / 2;
      dtmp = (double) (pdelay * search.pstep) / search.proflen;
      search.periods[ii] = 1.0 / (foldf + dtmp / T);
      dtmp = (double) (pdelay * search.pdstep) / search.proflen;
      fdots[ii] = phasedelay2fdot(dtmp, T);
      search.pdots[ii] = switch_pfdot(foldf, foldfd + fdots[ii]);
      if (cmd->searchpddP)
	fdotdots[ii] = phasedelay2fdotdot(dtmp, T);
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
      
      numdmtrials = 2*search.ndmfact*search.proflen + 1;
      lodm = cmd->dm - (numdmtrials-1)/2*ddm;
      if (cmd->nodmsearchP)
	good_dm_index = search.ndmfact*search.proflen;
      if (lodm < 0.0){
	lodm = 0.0;
	/* Find the closest DM to the requested DM */
	if (cmd->nodmsearchP){
	  double mindmerr=1000.0, dmerr, trialdm;
	  for (ii=0; ii<numdmtrials; ii++){  /* Loop over DMs */
	    trialdm = lodm + ii*ddm;
	    dmerr = fabs(trialdm - cmd->dm);
	    if (dmerr < mindmerr){
	      good_dm_index = ii;
	      mindmerr = dmerr;
	    }
	  }
	}
      }

      search.dms = gen_dvect(numdmtrials);
      search.numdms = numdmtrials;
      
      if (cmd->searchpddP)
	printf("  Searching %d DMs, %d periods, %d p-dots, and %d p-dotdots...\n", 
	       search.numdms, search.numperiods, search.numpdots, search.numpdots);
      else
	printf("  Searching %d DMs, %d periods, and %d p-dots...\n", 
	       search.numdms, search.numperiods, search.numpdots);

      /* De-disperse and combine the subbands */
      
      for (ii=0; ii<numdmtrials; ii++){  /* Loop over DMs */
	int numpdds=1;
	if (!cmd->nodmsearchP)
	  good_dm_index = ii;
	search.dms[ii] = lodm + ii*ddm;
	hifdelay = delay_from_dm(search.dms[ii], obsf[numchan - 1]);
	subbanddelays = subband_delays(numchan, cmd->nsub, 
				       search.dms[ii], idata.freq, 
				       idata.chan_wid, search.avgvoverc);
	for (jj=0; jj<cmd->nsub; jj++)
	  dmdelays[jj] = NEAREST_INT((subbanddelays[jj] - hifdelay) / 
				     dphase) % search.proflen;
	free(subbanddelays);
	combine_subbands(search.rawfolds, search.stats, cmd->npart, 
			 cmd->nsub, search.proflen, dmdelays, 
			 ddprofs, ddstats);
	
	/* Perform the Period and P-dot (and possibly P-dotdot) searches */

	if (cmd->searchpddP)
	  numpdds = numtrials;
	else
	  pddprofs = ddprofs;

	for (ll=0; ll<numpdds; ll++){
	  if (!(cmd->nosearchP && cmd->searchpddP))
	    good_pdd_index = ll;

	  /* Correct each part for the current pdotdot (if required) */

	  if (cmd->searchpddP){
	    for (kk=0; kk<cmd->npart; kk++){
	      profindex = kk*search.proflen;
	      pddelay = NEAREST_INT(fdotdot2phasedelay(fdotdots[ll], parttimes[kk]) * 
				    search.proflen);
	      shift_prof(ddprofs+profindex, search.proflen, pddelay, 
			 pddprofs+profindex);
	    }
	  }
	    
	  for (jj=0; jj<numtrials; jj++){
	    if (!cmd->nopdsearchP)
	      good_pd_index = jj;

	    /* Correct each part for the current pdot */

	    for (kk=0; kk<cmd->npart; kk++){
	      profindex = kk*search.proflen;
	      pddelay = NEAREST_INT(fdot2phasedelay(fdots[jj], parttimes[kk]) * 
				    search.proflen);
	      shift_prof(pddprofs+profindex, search.proflen, pddelay, 
			 pdprofs+profindex);
	    }
	    
	    /* Search over the periods */
	    
	    for (kk=0; kk<numtrials; kk++){
	      if (!cmd->nopsearchP)
		good_p_index = kk;
	      pdelay = search.pstep*(kk-(numtrials-1)/2);
	      combine_profs(pdprofs, ddstats, cmd->npart, search.proflen, 
			  pdelay, currentprof, &currentstats);
	      if (ii==good_dm_index &&
		  ll==good_pdd_index &&
		  jj==good_pd_index &&
		  kk==good_p_index){
		if (cmd->nosearchP ||
		    (currentstats.redchi > beststats.redchi && !cmd->nosearchP)){
		  search.bestdm = search.dms[ii];
		  if (idata.bary){
		    search.bary.p1 = search.periods[kk];
		    search.bary.p2 = search.pdots[jj];
		    if (cmd->searchpddP)
		      search.bary.p3 = switch_pfdotdot(1.0/search.periods[kk], 
						       foldfd+fdots[jj], foldfdd+fdotdots[ll]);
		  } else {
		    search.topo.p1 = search.periods[kk];
		    search.topo.p2 = search.pdots[jj];
		    if (cmd->searchpddP)
		      search.topo.p3 = switch_pfdotdot(1.0/search.periods[kk], 
						       foldfd+fdots[jj], foldfdd+fdotdots[ll]);
		  }
		  beststats = currentstats;
		  memcpy(bestprof, currentprof, sizeof(double) * 
			 search.proflen);
		}
	      }
	    }
	  }
	}
      }
      free(dmdelays);
      free(ddprofs);
      free(ddstats);
      
      /* We are not searching through DM space */
      
    } else {
      int numpdds=1;

      if (cmd->searchpddP)
	printf("  Searching %d periods, %d p-dots, and %d p-dotdots...\n", 
	       search.numperiods, search.numpdots, search.numpdots);
      else
	printf("  Searching %d periods, and %d p-dots...\n", 
	       search.numperiods, search.numpdots);

      /* Perform the Period and P-dot (and possibly P-dotdot) searches */

      if (cmd->searchpddP)
	numpdds = numtrials;
      else
	pddprofs = search.rawfolds;

      for (ll=0; ll<numpdds; ll++){
	if (!(cmd->nosearchP && cmd->searchpddP))
	  good_pdd_index = ll;

	/* Correct each part for the current pdotdot (if required) */
	
	if (cmd->searchpddP){
	  for (kk=0; kk<cmd->npart; kk++){
	    profindex = kk*search.proflen;
	    pddelay = NEAREST_INT(fdotdot2phasedelay(fdotdots[ll], 
						     parttimes[kk]) * search.proflen);
	    shift_prof(search.rawfolds+profindex, search.proflen, pddelay, 
		       pddprofs+profindex);
	  }
	}
	
	/* Perform the P-dot and Period searches */
	
	for (jj=0; jj<numtrials; jj++){
	  if (!cmd->nopdsearchP)
	    good_pd_index = jj;
	  
	  /* Correct each part for the current pdot */
	  
	  for (kk=0; kk<cmd->npart; kk++){
	    profindex = kk*search.proflen;
	    pddelay = NEAREST_INT(fdot2phasedelay(fdots[jj], parttimes[kk]) * 
				  search.proflen);
	    shift_prof(pddprofs+profindex, search.proflen, pddelay, 
		       pdprofs+profindex);
	  }
	  /* Search over the periods */
	  
	  for (kk=0; kk<numtrials; kk++){
	    if (!cmd->nopsearchP)
	      good_p_index = kk;
	    pdelay = search.pstep*(kk-(numtrials-1)/2);
	    combine_profs(pdprofs, search.stats, cmd->npart, search.proflen, 
			  pdelay, currentprof, &currentstats);
	    if (ll==good_pdd_index &&
		jj==good_pd_index &&
		kk==good_p_index){
	      if (cmd->nosearchP ||
		  (currentstats.redchi > beststats.redchi && !cmd->nosearchP)){
		if (idata.bary){
		  search.bary.p1 = search.periods[kk];
		  search.bary.p2 = search.pdots[jj];
		  if (cmd->searchpddP)
		    search.bary.p3 = switch_pfdotdot(1.0/search.periods[kk], 
						     foldfd+fdots[jj], foldfdd+fdotdots[ll]);
		} else {
		  search.topo.p1 = search.periods[kk];
		  search.topo.p2 = search.pdots[jj];
		  if (cmd->searchpddP)
		    search.topo.p3 = switch_pfdotdot(1.0/search.periods[kk], 
						     foldfd+fdots[jj], foldfdd+fdotdots[ll]);
		}
		beststats = currentstats;
		memcpy(bestprof, currentprof, sizeof(double) * 
		       search.proflen);
	      }
	    }
	  }
	}
      }
    }
    free(pdprofs);
    free(currentprof);
    free(fdots);
    if (cmd->searchpddP){
      free(fdotdots);
      free(pddprofs);
    }
  }
  printf("  Done searching.\n\n");

  {
    double perr, pderr, pdderr;
    char out[100];
  
    printf("Maximum reduced chi-squared found  =  %-.5f\n", 
	   beststats.redchi);
    if (cmd->nsub > 1)
      printf("Best DM     (pc cm^-3)  =  %-.4f\n", search.bestdm);
    
    /* Convert best params from/to barycentric to/from topocentric */
    
    if (!RAWDATA || cmd->polycofileP){

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

  if (cmd->noxwinP)
    prepfold_plot(&search, &pflags, 0);
  else 
    prepfold_plot(&search, &pflags, 1);

  /* Free our memory  */

  delete_prepfoldinfo(&search);
  free(data);
  if (!cmd->outfileP)
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
  if (RAWDATA){
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
  
