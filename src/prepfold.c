#include "presto.h"
#include "prepfold_cmd.h"
#include "multibeam.h"

/* Some function definitions */

int (*readrec_ptr)(FILE * file, float *data, int numpts, \
		   double *dispdelays, int numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, int numpts, \
		double *dispdelays, int numchan);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE *infile = NULL, *filemarker;
  float *data = NULL;
  double p_psr = 0.0, pdot_psr = 0.0, freq = 0.0, dfdt = 0.0;
  double difft, tt, nc, pl, diff_epoch = 0.0, recdt = 0.0;
  double orb_baryepoch = 0.0, topoepoch = 0.0, baryepoch = 0.0;
  double dtmp = 0.0, dtmp2 = 0.0, tmptopoepoch = 0.0, tmpbaryepoch = 0.0;
  double *Ep = NULL, startE = 0.0, orbdt = 1.0, orbdtdays = 0.0;
  double tdf = 0.0, N = 0.0, dt = 0.0, T, endtime = 0.0, dtdays;
  double *btoa = NULL, *voverc = NULL, *bobsf = NULL, *tobsf = NULL;
  double fakeonoffpair[2], *profs = NULL;
  double *delta_ts = NULL, *topo_ts = NULL;
  char obs[3], ephem[10], *outfilenm, *rootfilenm;
  char pname[30], rastring[50], decstring[50], *cptr;
  int numchan = 1, binary = 0, np, pnum, lorec = 0, hirec = 0;
  int slen, numonoffpairs = 1, ptsperrec = 1;
  long ii, jj, kk, numbarypts = 0, numpts = 0;
  long numfolded = 0, totnumfolded = 0;
  long numrecs = 0, totnumrecs = 0;
  long numbinpoints, proflen, currentrec;
  unsigned long numrec = 0;
  multibeam_tapehdr hdr;
  fourierprops rzwcand;
  orbitparams orb;
  psrparams psr;
  psrdatabase pdata;
  infodata idata, rzwidata;
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
  printf("  Used for DM and period determination of a PSR candidate.\n");
  printf("                 by Scott M. Ransom\n");
  printf("                    25 Nov, 1999\n\n");

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
    printf("Reading input  data from '%s'.\n", cmd->argv[0]);
    printf("Reading information from '%s.inf'.\n\n", rootfilenm);

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
    
  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    chkfileseek(infile, 0L, sizeof(char), SEEK_SET);
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

    lorec = (int) (cmd->startT / T * numrec + DBLCORRECT);
    hirec = (int) (cmd->endT / T * numrec + DBLCORRECT);
    numrec = hirec - lorec;
    T = numrec * recdt;
    N = numrec * ptsperrec;
    
    /* 1000.0 extra secs allows for worst case barycentric correction */

    endtime = T + 1000;

    /* Topocentric and barycentric times of folding epoch data */

    if (idata.mjd_i && idata.mjd_f) {
      topoepoch = (double) idata.mjd_i + idata.mjd_f + 
	lorec * recdt / SECPERDAY;
      barycenter(&topoepoch, &baryepoch, &dtmp, 1, rastring,
		 decstring, obs, ephem);
    }

    /* The data collection routine to use */

    readrec_ptr = read_multibeam_subbands;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    numpts = ptsperrec;
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

    numpts = 1024;
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
    p_psr = psr.p;
    pdot_psr = psr.pd;
    freq = psr.f;
    dfdt = psr.fd;
    strcpy(pname, psr.jname);
    
    /* If the user specifies all of the binaries parameters */	
    
  } else if (cmd->binaryP) {				
							
    /* Assume that the psr characteristics were measured at the time */
    /* of periastron passage (cmd->To)                               */
    
    difft = SECPERDAY * (topoepoch - cmd->To);				
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
      freq = (rzwcand.r - 0.5 * rzwcand.z) / 
	(rzwidata.dt * rzwidata.N);
      p_psr = 1.0 / freq;
      dfdt = rzwcand.z / ((rzwidata.dt * rzwidata.N) * 
			  (rzwidata.dt * rzwidata.N));
      pdot_psr = -dfdt / (freq * freq);
    } else {
      printf("\nCould not read the rzwfile.\nExiting.\n\n");
      exit(1);
    }

  }
  
  /* Determine the pulsar parameters to fold if we are not getting   */
  /* the data from a .cand file, the pulsar database, or a makefile. */
  
  if (!cmd->rzwcandP && !cmd->psrnameP) {
      
    if (cmd->pP) {							
      p_psr = cmd->p;							
      freq = 1.0 / p_psr;						
    }									
    if (cmd->freqP) {							
      freq = cmd->freq;						
      p_psr = 1.0 / freq;						
    }									
    if (cmd->pdot != 0.0) {						
      pdot_psr = cmd->pdot;						
      dfdt = -pdot_psr / (p_psr * p_psr);				
    }									
    if (cmd->dfdt != 0.0) {						
      dfdt = cmd->dfdt;						
      pdot_psr = -dfdt / (freq * freq);				
    }									
  }
  
  /* Determine the length of the profile */				
  
  if (cmd->proflenP) {						
    proflen = cmd->proflen;						
  } else {								
    proflen = (long) (p_psr / dt + 0.5);				
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

    /* Convert Eccentric anomaly to time delays */			
    
    orb.w *= DEGTORAD;
    E_to_phib(Ep, numbinpoints, &orb);		
  }

  /* Determine the spacing to use when generating */
  /* TEMPO barycentric corrections                */

  orbdtdays = orbdt / SECPERDAY;

  /* Output some informational data on the screen and to the */
  /* output file.                                            */
  
  fprintf(stdout, "\n");
  filemarker = stdout;
  for (ii = 0 ; ii < 1 ; ii++){
    if (cmd->psrnameP) {
      fprintf(filemarker, "Pulsar                       =  %s\n", pname);
    }
    if (tmptopoepoch != 0.0) {
      fprintf(filemarker, "Folding (topo) epoch  (MJD)  =  %-17.11f\n",
	      tmptopoepoch);
    }
    if (tmpbaryepoch != 0.0) {
      fprintf(filemarker, "Folding (bary) epoch  (MJD)  =  %-17.11f\n",
	      tmpbaryepoch);
    }
    fprintf(filemarker, "Data pt duration (dt)   (s)  =  %-.12f\n", dt);
    fprintf(filemarker, "Number of data points        =  %ld\n", (long) N);
    fprintf(filemarker, "Number of profile bins       =  %ld\n", proflen);
    fprintf(filemarker, "Folding period          (s)  =  %-.15f\n", p_psr);
    if (pdot_psr != 0.0) {
      fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10e\n",
	      pdot_psr);
    }
    fprintf(filemarker, "Folding frequency      (hz)  =  %-.12f\n", freq);
    if (pdot_psr != 0.0) {
	fprintf(filemarker, "Folding f-dot        (hz/s)  =  %-.8e\n", dfdt);
    }
    if (binary) {
      fprintf(filemarker, "Orbital period          (s)  =  %-.10f\n", orb.p);
      fprintf(filemarker, "a*sin(i)/c (x)     (lt-sec)  =  %-.10f\n", orb.x);
      fprintf(filemarker, "Eccentricity                 =  %-.10f\n", orb.e);
      fprintf(filemarker, "Longitude of peri (w) (deg)  =  %-.10f\n", \
	      orb.w/DEGTORAD);
      if (cmd->psrnameP){
	fprintf(filemarker, "Epoch of periapsis    (MJD)  =  %-17.11f\n",
		baryepoch + orb_baryepoch);
      }
    }
  }

  /* The number of topo to bary time points to generate with TEMPO */

  numbarypts = (long) floor((endtime - 400.0)/ orbdt);

  if (!strcmp(idata.band, "Radio")) {
    
    /* The topocentric spacing between channels */
    
    tdf = idata.chan_wid;
    
    /* The topocentric observation frequencies */
    
    tobsf = gen_dvect(numchan);
    
    tobsf[0] = idata.freq;
    for (ii = 0; ii < numchan; ii++) {
      tobsf[ii] = tobsf[0] + ii * tdf;
    }
    
  } else {
    
    printf("\n  This routine is only for use with multi-channel radio\n");
    printf("  data.  Exiting.\n\n");
    exit(1);
    
  }
  
  /* Main loop if we are not barycentering... */
  
  if (cmd->nobaryP) {
    
    printf("Massaging the data ...\n\n");
    
    /*****  Need to do this  *****/
    
  /* Main loop if we are barycentering... */

  } else {
    
    /* Allocate some arrays */
    
    bobsf = gen_dvect(numchan);
    btoa = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);
    delta_ts = gen_dvect(numbarypts);
    topo_ts = gen_dvect(numbarypts);
    for (ii = 0 ; ii < numbarypts ; ii++){
      /* topocentric times in days from data start */
      topo_ts[ii] = topoepoch + (double) ii * orbdtdays;
    }

    /* Call TEMPO for the barycentering */

    printf("\nGenerating barycentric corrections...\n");
    barycenter(topo_ts, btoa, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    free(voverc);

    /* printf("Given... (topoepoch = %15.9f, baryepoch = %15.9f)\n", \ */
    /*        topoepoch, baryepoch); */
    /* for (i = 0 ; i < 20 ; i++){ */
    /*   printf("Tt = %15.9f  Tb = %15.9f  orbcorr = %15.9f\n", \ */
    /* 	 (topo_ts[i] - topoepoch) * SECPERDAY, \ */
    /* 	 (btoa[i] - baryepoch) * SECPERDAY, \ */
    /* 	 (btoa[i] - baryepoch) * SECPERDAY - Ep[i]); */
    /* } */

    printf("   Insure you check the file %s.tempo_out for\n", \
	   cmd->argv[0]);
    printf("   errors from TEMPO when complete.\n\n");

    printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);
    
    /* Modify the binary times so that they refer to topocentric   */
    /* reference times.  This way we can barycenter while we fold. */

    diff_epoch = (baryepoch - topoepoch) * SECPERDAY;
    if (binary){
      for (ii = 0 ; ii < numbarypts ; ii++){
	dtmp = (btoa[ii] - baryepoch) * SECPERDAY;
	dtmp2 = (topo_ts[ii] - btoa[ii]) * SECPERDAY + diff_epoch;
	/* delta_ts are offsets from the appropriate topocentric times */
	delta_ts[ii] = lin_interp_E(Ep, dtmp, 0.0, orbdt, \
				    endtime - orbdt) + dtmp2;
      }
      free(Ep);
    } else {
      for (ii = 0 ; ii < numbarypts ; ii++){
	/* delta_ts are offsets from the appropriate topocentric times */
	delta_ts[ii] = diff_epoch - (btoa[ii] - topo_ts[ii]) * SECPERDAY;
      }
    }

    /* printf("Calculated...\n"); */
    /* for (i = 0 ; i < 20 ; i++){ */
    /*   printf("Tt = %15.9f  Tb = %15.9f\n", i * orbdt, \ */
    /* 	 i * orbdt - lin_interp_E(delta_ts, i * orbdt, \ */
    /* 				  0.0, orbdt, endtime - orbdt)); */
    /* } */

    free(btoa);
    free(topo_ts);
    
    /* Allocate and initialize some arrays and other information */
    
    data = gen_fvect(numchan * numpts);
    profs = gen_dvect(numchan * proflen);
    tt = 0;
    for (ii = 0 ; ii < numchan * proflen ; ii++) profs[ii] = 0.0;
    currentrec = 0;
    
    /* Move to the correct starting record */
    
    currentrec = skip_to_multibeam_rec(infile, lorec);
    
    /* Count the total number of records to fold */

    printf("Completed record 0 of %ld", totnumrecs);
    fakeonoffpair[0] = 0;
    fakeonoffpair[1] = numpts-1;

    /* Now, fold the data for each channel */
    
    /* Step through onoffpairs */

    for (ii = 0 ; ii < numonoffpairs ; ii++){
      
      /* Step through records */

      for (jj = lorec; jj <= hirec; jj++){

 	printf("\rCompleted record #%ld of %ld", numrecs, totnumrecs);
	fflush(stdout);

	read_mb_chan_to_vecs(infile, data, numpts, numchan);

	/* tt is topocentric seconds from data start */
	tt = jj * numpts * dt;

	/* Step through channels */
	for (kk = 0 ; kk < numchan ; kk++){
	  numfolded = 0;
	  fold(data + kk * numpts, numpts, dt, tt, profs + kk * proflen, \
	       proflen, freq, dfdt, 0.0, 0, delta_ts, 0.0, orbdt, \
	       numbarypts, fakeonoffpair, &numfolded);
	}
	totnumfolded += numfolded;
	numrecs++;
      }
    } 
  }

  printf("\rCompleted record #%ld of %ld\n\n", totnumrecs, totnumrecs);
  printf("Folded %d profiles with %ld points each.\n", numchan, totnumfolded);

  fclose(infile);

  /* Write our results. */

  printf("\nWriting %s.\n", outfilenm);
  infile = chkfopen(outfilenm,"wb");
  nc = numchan;
  pl = proflen;
  chkfwrite(&nc, sizeof(double), 1, infile);
  chkfwrite(&pl, sizeof(double), 1, infile);
  chkfwrite(&p_psr, sizeof(double), 1, infile);
  chkfwrite(&topoepoch, sizeof(double), 1, infile);
  chkfwrite(&baryepoch, sizeof(double), 1, infile);
  chkfwrite(&idata.freq, sizeof(double), 1, infile);
  chkfwrite(&idata.chan_wid, sizeof(double), 1, infile);
  chkfwrite(profs, sizeof(double), (unsigned long) (numchan * proflen), \
	    infile);
  fclose(infile);
  printf("Done.\n\n");

  /* Close the files and cleanup */
  free(data);
  free(tobsf);
  free(bobsf);
  free(profs);
  free(delta_ts);
  if (idata.onoff) free(idata.onoff);
  return (0);
}

    
int read_floats(FILE * file, float *data, int numpts, \
		double *dispdelays, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
  /* Read the raw data and return numbar read */

  *dispdelays = *dispdelays;
  return chkfread(data, sizeof(float), (unsigned long) (numpts * numchan), \
		  file) / numchan;
}




