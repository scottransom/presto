#include "presto.h"
#include "prepfold_cmd.h"
#include "multibeam.h"

/* Data points per record */

#define PTSPERREC 1536

/* Some function definitions */

int (*readrec_ptr) (FILE * infile, float *data,
		    long numpts, long numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, long numpts, \
		double *dispdelays, long numchan);

/* The main program */

int main(int argc, char *argv[])
/* This routine generates barycentric info for a multichannel radio */
/* data file and then folds each channel at a specific frequency.   */
/* The resulting profiles are plotted and used to solve for the     */
/* dispersion mearure to the pulsar.                                */
/* It uses either the DE200 or the DE405 ephemeris.                 */
{

  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */

  FILE *infile = NULL, *filemarker;
  float *data = NULL;
  double p_psr = 0.0, pdot_psr = 0.0, freq = 0.0, dfdt = 0.0;
  double difft, tt, nc, pl, diff_epoch = 0.0, recdt = 0.0;
  double orb_baryepoch = 0.0, topoepoch = 0.0, baryepoch = 0.0;
  double dtmp = 0.0, dtmp2 = 0.0, tmptopoepoch = 0.0, tmpbaryepoch = 0.0;
  double *Ep = NULL, *tp = NULL, startE = 0.0, orbdt = 1.0, orbdtdays = 0.0;
  double tdf = 0.0, N = 0.0, dt = 0.0, T, endtime = 0.0, dtdays;
  double *btoa = NULL, *voverc = NULL, *bobsf = NULL, *tobsf = NULL;
  double onoffpairs[40], fakeonoffpair[2], *onptr, *offptr, *profs = NULL;
  double *delta_ts = NULL, *topo_ts = NULL;
  char obs[3], ephem[10], datafilenm[200];
  char pname[30], rastring[50], decstring[50], *ctmp;
  int numchan = 1, binary = 0, np, pnum, numonoffpairs = 1;
  long i, j, k, numbarypts = 0, numpts = 0, numfolded = 0, totnumfolded = 0;
  long numrecs = 0, totnumrecs = 0;
  long numbinpoints, proflen, currentrec;
  unsigned long numrec = 0;
  multibeam_tapehdr hdr;
  fourierprops rzwcand;
  orbitparams orb;
  psrparams psr;
  psrdatabase pdata;
  infodata idata;
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
  printf("           Pulsar Raw-Data Folding Routine\n");
  printf("  Used primarily for DM determination and RFI identification.\n");
  printf("                 by Scott M. Ransom\n");
  printf("                    30 Jun, 1998\n\n");

  /* Determine the output file name */

  sprintf(datafilenm, "%s.msv", cmd->argv[0]);

  /* What ephemeris will we use?  (Default is DE200) */
  
  if (cmd->de405P) {
    strcpy(ephem, "DE405");
  } else {
    strcpy(ephem, "DE200");
  }
    
  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* Open the raw input file */

    infile = chkfopen(cmd->argv[0], "rb");
    numrec = chkfilelen(infile, RECLEN);

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    rewind(infile);
    /* print_multibeam_hdr(&hdr); */
    multibeam_hdr_to_inf(&hdr, &idata);
    if (idata.object) {						
      printf("Folding a %s candidate from '%s'.\n", \
	     idata.object, cmd->argv[0]);
    } else {								
      printf("Folding a candidate from '%s'.\n", cmd->argv[0]);
    }									

    /* Define the RA and DEC of the observation */
  
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Define some variables */

    dt = idata.dt;
    dtdays = idata.dt / SECPERDAY;
    recdt = dt * PTSPERREC;
    N = numrec * PTSPERREC;
    T = dt * N;
    /* 1000.0 extra secs allows for worst case barycentric correction */
    endtime = T + 1000;
    if (idata.mjd_i && idata.mjd_f) {
      /* Topocentric time of the beginning of the data */
      topoepoch = (double) idata.mjd_i + idata.mjd_f;			
      barycenter(&topoepoch, &baryepoch, &dtmp, 1, rastring, \
		 decstring, obs, ephem);
    }

    /* The data collection routine to use */

    readrec_ptr = read_mb_chan_to_vecs;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    numpts = PTSPERREC;

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
    if (psr.ntype & 8){
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
    }									
    
    get_rzw_cand(cmd->rzwfile, cmd->rzwcand, &rzwcand);		
    freq = (rzwcand.r - 0.5 * rzwcand.z) / (dt * idata.N);	
    p_psr = 1.0 / freq;						
    dfdt = rzwcand.z / ((dt * idata.N) * (dt * idata.N));	
    pdot_psr = -dfdt / (freq * freq);					
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
    
    startE = keplars_eqn(orb.t, 1.0 / orb.p, orb.e, 1.0E-15);
    if (endtime > 2048) orbdt = 0.5;
    else orbdt = endtime / 4096.0;
    numbinpoints = (long) floor(endtime/orbdt + 0.5) + 1;
    Ep = gen_dvect(numbinpoints);
    tp = gen_dvect(numbinpoints);
    dorbint(Ep, startE, tp, 0.0, numbinpoints, orbdt, &orb);
    free(tp);

    /* Convert Eccentric anomaly to time delays */			
    
    orb.w *= DEGTORAD;
    E_to_phib(Ep, numbinpoints, &orb);		
    
  }
  
  /* Determine the spacing to use when generating */
  /* TEMPO barycentric corrections                */
  
  orbdtdays = orbdt / SECPERDAY;
  
  /* Determine if we will take any breaks during the folding */
  
  for (i = 0 ; i < 40 ; i++) onoffpairs[i] = 0.0;
  i = 0;
  if (cmd->onoffP){
    numonoffpairs = 0;
    ctmp = strtok(cmd->onoff, " \t\n");
    while (NULL != ctmp){
      onoffpairs[i] = strtod(ctmp, NULL);
      if (onoffpairs[i] < 0.0 || onoffpairs[i] > 1.0){
	printf("\nonoff pairs must be between 0.0 and 1.0 inclusive.\n\n");
	exit(1);
      }
      if (onoffpairs[i] == 1.0) onoffpairs[i] = numrec-1.0;
      else onoffpairs[i] = floor(onoffpairs[i] * (numrec-1));
      if (i >= 1 && onoffpairs[i] < onoffpairs[i-1]){
	printf("\nonoff values must increase from left to right.\n\n");
	exit(1);
      }
      if (onoffpairs[i] >= numrec) ctmp = NULL;
      if (!(i & 1)) numonoffpairs++;
      i++;
      ctmp = strtok(NULL, " \t\n");
    }
    if (i & 1) onoffpairs[i] = numrec - 1;
    
    /* Adjust various values to the beginning of the first bin */
    
    tmptopoepoch = topoepoch + onoffpairs[0] * dt * PTSPERREC / SECPERDAY;
    if (cmd->pkmbP) {
      barycenter(&tmptopoepoch, &tmpbaryepoch, &dtmp2, 1, rastring, \
		 decstring, obs, ephem);
      dtmp2 = (tmpbaryepoch - baryepoch) * SECPERDAY;
    }
  } else {
    onoffpairs[0] = 0.0;
    onoffpairs[1] = numrec-1.0;
  }
  onptr = onoffpairs;
  offptr = onoffpairs + 1;
  
  /* Output some informational data on the screen and to the */	
  /* output file.                                            */	
  
  fprintf(stdout, "\n");						
  filemarker = stdout;
  dtmp = p_psr + dtmp2 * pdot_psr;
  for (i = 0 ; i < 1 ; i++){
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
    fprintf(filemarker, "Folding period          (s)  =  %-.15f\n", dtmp);
    if (pdot_psr != 0.0) {						
      fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10e\n",	
	      pdot_psr);						
    }									
    fprintf(filemarker, "Folding frequency      (hz)  =  %-.12f\n",	
	    1.0 / dtmp);
    if (pdot_psr != 0.0) {						
	fprintf(filemarker, "Folding f-dot        (hz/s)  =  %-.8e\n",	
		-pdot_psr / (dtmp * dtmp));				
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
    for (i = 0; i < numchan; i++) {
      tobsf[i] = tobsf[0] + i * tdf;
    }
    
  } else {
    
    printf("\n  This routine is only for use with multi-channel radio\n");
    printf("  data.  Exiting.\n\n");
    exit(1);
    
  }
  
  /* Main loop if we are not barycentering... */
  
  if (cmd->nobaryP) {
    
    printf("Massaging the data ...\n\n");
    
    /* Main loop if we are barycentering... */

    /*****  Need to do this  *****/
    
  } else {
    
    /* Allocate some arrays */
    
    bobsf = gen_dvect(numchan);
    btoa = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);
    delta_ts = gen_dvect(numbarypts);
    topo_ts = gen_dvect(numbarypts);
    for (i = 0 ; i < numbarypts ; i++){
      /* topocentric times in days from data start */
      topo_ts[i] = topoepoch + (double) i * orbdtdays;
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
      for(i = 0 ; i < numbarypts ; i++){
	dtmp = (btoa[i] - baryepoch) * SECPERDAY;
	dtmp2 = (topo_ts[i] - btoa[i]) * SECPERDAY + diff_epoch;
	/* delta_ts are offsets from the appropriate topocentric times */
	delta_ts[i] = lin_interp_E(Ep, dtmp, 0.0, orbdt, \
				   endtime - orbdt) + dtmp2;
      }
      free(Ep);
    } else {
      for(i = 0 ; i < numbarypts ; i++){
	/* delta_ts are offsets from the appropriate topocentric times */
	delta_ts[i] = diff_epoch - (btoa[i] - topo_ts[i]) * SECPERDAY;
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
    for (i = 0 ; i < numchan * proflen ; i++) profs[i] = 0.0;
    currentrec = 0;
    
    /* Move to the correct starting record */
    
    currentrec = skip_to_multibeam_rec(infile, \
				       (long) (*onptr + DBLCORRECT));
    
    /* Count the total number of records to fold */

    for (i = 0 ; i < numonoffpairs ; i++){
      /* onoffpairs are records -- not bins or fractional times  */
      totnumrecs += onoffpairs[2*i+1] - onoffpairs[2*i] + 1;
    }
    printf("Completed record 0 of %ld", totnumrecs);
    fakeonoffpair[0] = 0;
    fakeonoffpair[1] = numpts-1;

    /* Now, fold the data for each channel */
    
    /* Step through onoffpairs */

    for (i = 0 ; i < numonoffpairs ; i++){
      
      /* Step through records */

      for (j = (long) (*onptr + DBLCORRECT) ; \
	   j <= (long) (*offptr + DBLCORRECT); \
	   j++){
 	printf("\rCompleted record #%ld of %ld", numrecs, totnumrecs);
	fflush(stdout);

	read_mb_chan_to_vecs(infile, data, numpts, numchan);

	/* tt is topocentric seconds from data start */
	tt = j * numpts * dt;

	/* Step through channels */
	for (k = 0 ; k < numchan ; k++){
	  numfolded = 0;
	  fold(data + k * numpts, numpts, dt, tt, profs + k * proflen, \
	       proflen, freq, dfdt, 0.0, 0, delta_ts, 0.0, orbdt, \
	       numbarypts, fakeonoffpair, &numfolded);
	}
	totnumfolded += numfolded;
	numrecs++;
      }
      if (numonoffpairs > 1){
	onptr += 2;
	offptr += 2;
	currentrec = skip_to_multibeam_rec(infile, \
					   (long) (*onptr + DBLCORRECT));
      }
    } 
  }
  printf("\rCompleted record #%ld of %ld\n\n", totnumrecs, totnumrecs);
  printf("Folded %d profiles with %ld points each.\n", numchan, totnumfolded);

  fclose(infile);

  /* Write our results. */

  printf("\nWriting %s.\n", datafilenm);
  infile = chkfopen(datafilenm,"wb");
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

    
int read_floats(FILE * file, float *data, long numpts, \
		double *dispdelays, long numchan)
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




