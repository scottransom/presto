#include <limits.h>
#include "presto.h"
#include "prepdata_cmd.h"
#include "multibeam.h"

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

/* Some function definitions */

int (*readrec_ptr)(FILE * file, float *data, int numpts,
		   double *dispdelays, int numchan);
static int read_floats(FILE * file, float *data, int numpts,
		       double *dispdelays, int numchan);
static void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE *infile = NULL, *outfile;
  float *outdata=NULL;
  double tdf=0.0, dtmp=0.0, barydispdt=0.0;
  double *dispdt, *tobsf=NULL, tlotoa=0.0, blotoa=0.0;
  double max=-9.9E30, min=9.9E30, var=0.0, avg=0.0;
  char obs[3], ephem[10], *datafilenm, *rootfilenm, *outinfonm;
  char rastring[50], decstring[50], *cptr;
  int numchan=1, newper=0, oldper=0, slen;
  long ii, numbarypts=0, worklen=65536, next_pow_of_two=0;
  long numread=0, wrote=0, totwrote=0;
  multibeam_tapehdr hdr;
  infodata idata;
  Cmdline *cmd;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
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
  printf("           Pulsar Data Preparation Routine\n");
  printf("    Type conversion, de-dispersion, barycentering.\n");
  printf("                 by Scott M. Ransom\n");
  printf("           Last Modification:  5 July, 2000\n\n");

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

  /* Determine the other file names and open the output data file */

  slen = strlen(cmd->outfile)+6;
  datafilenm = (char *)malloc(slen+1);
  datafilenm[slen] = '\0';
  sprintf(datafilenm, "%s.dat", cmd->outfile);
  outfile = chkfopen(datafilenm, "wb");
  sprintf(idata.name, "%s", cmd->outfile);
  outinfonm = (char *)malloc(slen+1);
  outinfonm[slen] = '\0';
  sprintf(outinfonm, "%s.inf", cmd->outfile);
  printf("Writing output data to '%s'.\n", datafilenm);
  printf("Writing information to '%s'.\n\n", outinfonm);

  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    rewind(infile);
    multibeam_hdr_to_inf(&hdr, &idata);
    idata.dm = cmd->dm;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    worklen = DATLEN * 8 / numchan;

    /* Set idata.N to the actual size of the raw data file */

    idata.N = (double) (chkfilelen(infile, RECLEN)) * worklen;

    /* Compare the size of the data to the size of output we request */

    if (cmd->numoutP && cmd->numout != idata.N) {
      dtmp = idata.N;
      idata.N = cmd->numout;
    }

    /* Write a new info-file */

    writeinf(&idata);
    if (cmd->numoutP && cmd->numout != idata.N) {
      idata.N = dtmp;
    }

    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* The data collection routine to use */

    readrec_ptr = read_multibeam;

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;
  }

  /* Set-up values if we are using the Effelsberg-Berkeley Pulsar Processor */
  /*   NOTE:  This code is not yet implemented.                             */

  if (cmd->ebppP) {

    /* Read the first header file and generate an infofile from it */

    /* OBS code for TEMPO */

    strcpy(obs, "EF");

    /* The data collection routine to use */

    /* readrec_ptr = read_ebpp; */

    /* The number of data points to work with at a time */

    worklen = 1024;
  }

  /* Determine our initialization data if we do _not_ have Parkes */
  /* or Effelesberg data sets.                                    */

  if (!cmd->pkmbP && !cmd->ebppP) {

    /* The data collection routine to use */

    readrec_ptr = read_floats;

    /* If we will be barycentering... */

    if (!cmd->nobaryP) {

      /* The number of topo to bary time points to generate with TEMPO */

      numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;

      /* The number of data points to work with at a time */

      if (worklen > idata.N)
	worklen = idata.N;
      worklen = (long) (worklen / 1024) * 1024;

      /* What telescope are we using? */

      if (!strcmp(idata.telescope, "Arecibo")) {
	strcpy(obs, "AO");
      } else if (!strcmp(idata.telescope, "Parkes")) {
	strcpy(obs, "PK");
      } else if (!strcmp(idata.telescope, "Effelsberg")) {
	strcpy(obs, "EF");
      } else if (!strcmp(idata.telescope, "MMT")) {
	strcpy(obs, "MT");
      } else {
	printf("\nYou need to choose a telescope whose data is in\n");
	printf("$TEMPO/obsys.dat.  Exiting.\n\n");
	exit(1);
      }

    } else {			/* No barycentering... */

      /* The number of data points to work with at a time */

      if (worklen > idata.N)
	worklen = idata.N;
      worklen = (long) (worklen / 1024) * 1024;

    }

  }

  /* The topocentric epoch of the start of the data */

  tlotoa = (double) idata.mjd_i + idata.mjd_f;

  if (!strcmp(idata.band, "Radio") && (cmd->pkmbP || cmd->ebppP)) {

    /* The topocentric spacing between channels */

    tdf = idata.chan_wid;
    numchan = idata.num_chan;

    /* The topocentric observation frequencies */

    tobsf = gen_dvect(numchan);
    tobsf[0] = idata.freq;
    for (ii = 0; ii < numchan; ii++)
      tobsf[ii] = tobsf[0] + ii * tdf;

    /* The dispersion delays (in time bins) */

    dispdt = gen_dvect(numchan);

    if (cmd->nobaryP) {

      /* Determine our dispersion time delays for each channel */

      for (ii = 0; ii < numchan; ii++)
	dispdt[ii] = delay_from_dm(cmd->dm, tobsf[ii]);

      /* The highest frequency channel gets no delay                 */
      /* All other delays are positive fractions of bin length (dt)  */

      dtmp = dispdt[numchan - 1];
      for (ii = 0; ii < numchan; ii++)
	dispdt[ii] = (dispdt[ii] - dtmp) / idata.dt;
      worklen *= ((int)(fabs(dispdt[0])) / worklen) + 1;
    }

  } else {     /* For non Parkes or Effelsberg raw data */
    tobsf = gen_dvect(numchan);
    dispdt = gen_dvect(numchan);
    dispdt[0] = 0.0;
    if (!strcmp(idata.band, "Radio")){
      tobsf[0] = idata.freq + (idata.num_chan - 1) * idata.chan_wid;
      cmd->dm = idata.dm;
    } else {
      tobsf[0] = 0.0;
      cmd->dm = 0.0;
    }
  }

  /* Allocate our data array */
    
  outdata = gen_fvect(worklen);
    
  /* Main loop if we are not barycentering... */

  if (cmd->nobaryP) {

    /* Open our new output data file */
    
    outfile = chkfopen(datafilenm, "wb");
    
    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");
    
    while ((numread = 
	    (*readrec_ptr)(infile, outdata, worklen, dispdt, numchan))){
      
      /* Print percent complete */
      
      if (cmd->numoutP)
	newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
      else
	newper = (int) ((float) totwrote / idata.N * 100.0) + 1;
      if (newper > oldper) {
	printf("\rAmount Complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }
      
      /* Write the latest chunk of data, but don't   */
      /* write more than cmd->numout points.         */
      
      wrote = numread;
      if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	wrote = cmd->numout - totwrote;
      chkfwrite(outdata, sizeof(float), (unsigned long) wrote, outfile);
      
      /* Update the statistics */
      
      for (ii = 0; ii < numread; ii++)
	update_stats(totwrote + ii, outdata[ii], &min, &max, &avg, &var);
      totwrote += wrote;
      
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote >= cmd->numout))
	break;
     }
      
    /* Main loop if we are barycentering... */

  } else {

    double avgvoverc=0.0, *voverc=NULL;
    double *bobsf=NULL, *btoa=NULL, *ttoa=NULL;
    int *diffbins, *diffbinptr;

    /* What ephemeris will we use?  (Default is DE200) */

    if (cmd->de405P)
      strcpy(ephem, "DE405");
    else
      strcpy(ephem, "DE200");

    /* Define the RA and DEC of the observation */
    
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Allocate some arrays */

    bobsf = gen_dvect(numchan);
    btoa = gen_dvect(numbarypts);
    ttoa = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);
    for (ii = 0 ; ii < numbarypts ; ii++)
      ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(ttoa, btoa, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    for (ii = 0 ; ii < numbarypts ; ii++)
      avgvoverc =+ voverc[ii];
    avgvoverc /= numbarypts;
    free(voverc);
    blotoa = btoa[0];

    printf("   Insure you check the files tempoout_times.tmp and\n");
    printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n\n");

    printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

    /* Determine the initial dispersion time delays for each channel */

    for (ii = 0; ii < numchan; ii++) {
      bobsf[ii] = doppler(tobsf[ii], avgvoverc);
      dispdt[ii] = delay_from_dm(cmd->dm, bobsf[ii]);
    }

    /* The highest frequency channel gets no delay                   */
    /* All other delays are positive fractions of bin length (dt)    */

    barydispdt = dispdt[numchan - 1];
    for (ii = 0; ii < numchan; ii++)
      dispdt[ii] = (dispdt[ii] - barydispdt) / idata.dt;
    worklen *= ((int)(dispdt[0]) / worklen) + 1;

    /* If the data is de-dispersed radio data...*/

    if (!strcmp(idata.band, "Radio")) {
      printf("The DM of %.2f at the barycentric observing freq of %.3f MHz\n",
	     idata.dm, bobsf[numchan-1]);
      printf("   causes a delay of %f seconds compared to infinite freq.\n",
	     barydispdt);
      printf("   This delay is removed from the barycented times.\n\n");
    }
    printf("Topocentric epoch (at data start) is:\n");
    printf("   %17.11f\n\n", tlotoa);
    printf("Barycentric epoch (infinite obs freq at data start) is:\n");
    printf("   %17.11f\n\n", blotoa - (barydispdt / SECPERDAY));

    /* Convert the bary TOAs to differences from the topo TOAs in */
    /* units of bin length (dt) rounded to the nearest integer.   */

    dtmp = (btoa[0] - ttoa[0]);
    for (ii = 0 ; ii < numbarypts ; ii++)
      btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / idata.dt;

    /* Find the points where we need to add or remove bins */

    {
      int oldbin=0, currentbin, numdiffbins;
      double lobin, hibin, calcpt;

      numdiffbins = abs(NEAREST_INT(btoa[numbarypts-1])) + 1;
      diffbins = gen_ivect(numdiffbins);
      diffbinptr = diffbins;
      for (ii = 0 ; ii < numbarypts ; ii++){
	currentbin = NEAREST_INT(btoa[ii]);
	if (currentbin != oldbin){
	  calcpt = (currentbin > 0) ? oldbin + 0.5 : oldbin - 0.5;
	  while(fabs(calcpt) < fabs(btoa[ii])){
	    lobin = (ii-1) * TDT / idata.dt;
	    hibin = ii * TDT / idata.dt;
	    /* Negative bin number means remove that bin */
	    /* Positive bin number means add a bin there */
	    *diffbinptr = NEAREST_INT(LININTERP(calcpt, btoa[ii-1],
						btoa[ii], lobin, hibin));
            diffbinptr++;
	    calcpt = (currentbin > 0) ? calcpt + 1.0 : calcpt - 1.0;
	  }
	  oldbin = currentbin;
	}
      }
      *diffbinptr = INT_MAX;  /* Used as a marker */
    }
    diffbinptr = diffbins;

    /* Now perform the barycentering */

    printf("Massaging the data...\n\n");
    printf("Amount Complete = 0%%");
    
    while ((numread = 
	    (*readrec_ptr)(infile, outdata, worklen, 
			   dispdt, numchan))){
      
      /* Print percent complete */
      
      if (cmd->numoutP)
	newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
      else 
	newper = (int) ((float) totwrote / idata.N * 100.0) + 1;
      if (newper > oldper) {
 	printf("\rAmount Complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }
      
      /* Simply write the data if we don't have to add or */
      /* remove any bins from this batch.                 */
      
      if (totwrote + numread < abs(*diffbinptr)){

	/* Write the latest chunk of data, but don't   */
	/* write more than cmd->numout points.         */
	
	wrote = numread;
	if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	  wrote = cmd->numout - totwrote;
	chkfwrite(outdata, sizeof(float), wrote, outfile);
	
	/* Update the statistics */
	
	for (ii = 0; ii < numread; ii++)
	  update_stats(totwrote + ii, outdata[ii], &min, &max, 
		       &avg, &var);
	totwrote += wrote;
	
      } else {  /* Add or remove bins */
	float favg;
	int skip;

	/* Write the first part */

	wrote = abs(*diffbinptr) - totwrote;
	if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	  wrote = cmd->numout - totwrote;
	chkfwrite(outdata, sizeof(float), wrote, outfile);

	/* Update the statistics */
	
	for (ii = 0; ii < wrote; ii++)
	  update_stats(totwrote + ii, outdata[ii], &min, &max, 
		       &avg, &var);
	totwrote += wrote;

	/* Add a bin with the average value */

	skip = wrote;
	wrote = numread - skip;
	if (*diffbinptr > 0){
	  if (!(cmd->numoutP && (totwrote + 1) > cmd->numout)){
	    favg = (float) avg;
	    chkfwrite(&favg, sizeof(float), 1, outfile);
	    update_stats(totwrote, avg, &min, &max, 
			 &avg, &var);
	    totwrote++;
	  }
	} else {  /* Remove a bin from the data */
	  skip++;
	  wrote--;
	}

	/* Write the rest */

	if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	  wrote = cmd->numout - totwrote;
	chkfwrite(outdata+skip, sizeof(float), wrote, outfile);
	for (ii = 0; ii < wrote; ii++)
	  update_stats(totwrote + ii, outdata[skip+ii], &min, &max, 
		       &avg, &var);
	totwrote += wrote;

	/* Increment the diffbinptr */

	diffbinptr++;
      }

      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote >= cmd->numout))
	break;
    }

    /* Free the arrays used in barycentering */

    free(diffbins);
    free(bobsf);
    free(btoa);
    free(ttoa);
  }

  /* Print simple stats and results */

  var /= (totwrote - 1);
  printf("\n\nDone.\n\nSimple statistics of the output data:\n");
  printf("  Number of points written:    %ld\n", totwrote);
  printf("  Maximum value:               %.2f\n", max);
  printf("  Minimum value:               %.2f\n", min);
  printf("  Mean value:                  %.2f\n", avg);
  printf("  Standard deviation:          %.2f\n", sqrt(var));
  printf("\n");

  /* Calculate what the next power of two number of points would be */

  next_pow_of_two = next2_to_n(totwrote);

  /* Write the new info file for the output data */

  if (cmd->pkmbP) {
    idata.numonoff = 0;
  }
  idata.N = totwrote;
  if (!cmd->nobaryP) {
    idata.bary = 1;
    idata.mjd_i = (int) floor(blotoa - (barydispdt / SECPERDAY));
    idata.mjd_f = blotoa - (barydispdt / SECPERDAY) - idata.mjd_i;
  }
  if (cmd->numoutP && (cmd->numout > totwrote)) {
    idata.N = cmd->numout;
    free(idata.onoff);
    idata.numonoff = 2;
    idata.onoff = gen_dvect(2 * idata.numonoff);
    idata.onoff[0] = 0;
    idata.onoff[1] = totwrote - 1;
    idata.onoff[2] = cmd->numout - 1;
    idata.onoff[3] = cmd->numout - 1;
    printf("NOTE:  If there were on-off pairs in the input info file,\n");
    printf("       the output info file pairs are now incorrect!\n\n");
  }
  if ((cmd->pad0P || cmd->padavgP) && (next_pow_of_two != totwrote)) {
    idata.N = next_pow_of_two;
    free(idata.onoff);
    idata.numonoff = 2;
    idata.onoff = gen_dvect(2 * idata.numonoff);
    idata.onoff[0] = 0;
    idata.onoff[1] = totwrote - 1;
    idata.onoff[2] = next_pow_of_two - 1;
    idata.onoff[3] = next_pow_of_two - 1;
    printf("NOTE:  If there were on-off pairs in the input info file,\n");
    printf("       the output info file pairs are now incorrect!\n\n");
  }
  writeinf(&idata);

  /* Pad the file if needed to pow_of_two_length */

  if (cmd->pad0P || cmd->padavgP) {
    if (cmd->pad0P)
      for (ii = 0; ii < worklen; ii++)
	outdata[ii] = 0.0;
    else
      for (ii = 0; ii < worklen; ii++)
	outdata[ii] = avg;
    for (ii = 0; ii < (next_pow_of_two - totwrote) / worklen; ii++)
      chkfwrite(outdata, sizeof(float), worklen, outfile);
    chkfwrite(outdata, sizeof(float), \
	      (next_pow_of_two - totwrote) % worklen, outfile);
  }

  /* Pad the file if needed to the requested length */

  if (cmd->numoutP && (cmd->numout > totwrote)) {
    for (ii = 0; ii < worklen; ii++)
      outdata[ii] = avg;
    for (ii = 0; ii < (cmd->numout - totwrote) / worklen; ii++)
      chkfwrite(outdata, sizeof(float), worklen, outfile);
    chkfwrite(outdata, sizeof(float), 
	      (cmd->numout - totwrote) % worklen, outfile);
  }
  free(outdata);

  /* Close the files and cleanup */

  fclose(infile);
  fclose(outfile);
  free(tobsf);
  free(dispdt);
  free(rootfilenm);
  free(outinfonm);
  free(datafilenm);
  if (idata.onoff)
    free(idata.onoff);
  return (0);
}


static int read_floats(FILE * file, float *data, int numpts, \
		       double *dispdelays, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
  /* Read the raw data and return numbar read */

  *dispdelays = *dispdelays;
  return chkfread(data, sizeof(float), \
		  (unsigned long) (numpts * numchan), file) / numchan;
}


static void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var)
/* Update time series statistics using one-pass technique */
{
  double dev;

  /* Check the max and min values */
  
  if (x > *max) *max = x;
  if (x < *min) *min = x;
  
  /* Use clever single pass mean and variance calculation */
  
  dev = x - *avg;
  *avg += dev / (N + 1.0);
  *var += dev * (x - *avg);
}
