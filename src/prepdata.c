#include "presto.h"
#include "prepdata_cmd.h"
#include "multibeam.h"

/* This causes the barycentric motion to be calculated once per second */

#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Some function definitions */

int (*readrec_ptr)(FILE * file, float *data, int numpts,
		   double *dispdelays, int numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, int numpts,
		double *dispdelays, int numchan);
void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE *infile = NULL, *outfile;
  float *data, *outdata;
  double tdf=0.0, *dispdt, rdt=0.0, fact;
  double *tobsf=NULL, *bobsf=NULL, tlotoa=0.0, blotoa=0.0;
  double intime=0.0, intimelo=0.0, outtime=0.0, dtmp=0.0, lopoint, hipoint;
  double intimenext=0.0, outtimenext=0.0, deltapoints, rdeltapoints;
  double *delayptr=NULL, *delaytimeptr=NULL;
  double delaytlo=0.0, delaythi=0.0, delaylo=0.0, delayhi=0.0;
  double max=-9.9E30, min=9.9E30, var=0.0, avg=0.0, dev=0.0;
  double *btoa=NULL, *ttoa=NULL, *voverc=NULL, barydispdt=0.0;
  char obs[3], ephem[10], *datafilenm, *rootfilenm, *outinfonm;
  char rastring[50], decstring[50], *cptr;
  int numchan=1, newper=0, oldper=0, slen;
  long ii, numbarypts=0, worklen=65536, next_pow_of_two=0;
  long numread=0, wlen2, wrote=0;
  long totread=0, totwrote=0;
  unsigned long lobin, hibin, index=0, numbins, arrayoffset=0;
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
  printf("                    16 Nov, 1999\n\n");

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

  /* Main loop if we are not barycentering... */

  if (cmd->nobaryP) {

    /* Allocate and initialize our data array and others */
    
    data = gen_fvect(worklen);
    wlen2 = worklen * 2;
    outdata = gen_fvect(wlen2);
    for (ii = 0; ii < wlen2; ii++)
      outdata[ii] = 0.0;
    
    /* Open our new output data file */
    
    outfile = chkfopen(datafilenm, "wb");
    
    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");
    
    wrote = worklen;
    
    while ((numread = 
	   (*readrec_ptr)(infile, data, worklen, dispdt, numchan))){

      /* Insure we don't write out too much data */

      if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	wrote = cmd->numout - totwrote;

      chkfwrite(data, sizeof(float), (unsigned long) wrote, outfile);
      totwrote += wrote;

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

      /* Determine some data statistics */

      for (ii = 0; ii < wrote; ii++) {

	/* Find the max and min values */
	
	if (data[ii] > max)
	  max = data[ii];
	if (data[ii] < min)
	  min = data[ii];
	
	/* Use clever single pass mean and variance calculation */
	
	dev = data[ii] - avg;
	avg += dev / (totwrote + ii + 1);
	var += dev * (data[ii] - avg);
      }

      if (cmd->numoutP && (totwrote >= cmd->numout))
	break;
	
     }
      
    /* Main loop if we are barycentering... */

  } else {

    /* What ephemeris will we use?  (Default is DE200) */

    if (cmd->de405P) {
      strcpy(ephem, "DE405");
    } else {
      strcpy(ephem, "DE200");
    }

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
    blotoa = btoa[0];

    printf("   Insure you check the files tempoout_times.tmp and\n");
    printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n\n");

    printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

    /* Determine the initial dispersion time delays for each channel */

    for (ii = 0; ii < numchan; ii++) {
      bobsf[ii] = doppler(tobsf[ii], voverc[0]);
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

    /* Allocate and initialize some arrays and other information */
    
    data = gen_fvect(worklen);
    wlen2 = worklen * 2;
    outdata = gen_fvect(wlen2);
    for (ii = 0; ii < wlen2; ii++)
      outdata[ii] = 0.0;
    rdt = 1.0 / idata.dt;

    /* Convert the topo TOAs to seconds from start */
    /* Convert the bary TOAs to delays from the topo TOAs */

    dtmp = (ttoa[0] - btoa[0]);
    for (ii = 0 ; ii < numbarypts ; ii++){
      btoa[ii] = ((ttoa[ii] - btoa[ii]) - dtmp) * SECPERDAY;
      ttoa[ii] = TDT * ii;
    }

    /* Initiate some variables */
    
    intime = 0.0;
    arrayoffset = 1;  /* Beware nasty NR zero-offset kludges! */
    hunt(ttoa-1, numbarypts, intime, &arrayoffset);
    arrayoffset--;
    delaytimeptr = ttoa + arrayoffset;
    delayptr = btoa + arrayoffset;
    delaytlo = *delaytimeptr;
    delaythi = *(delaytimeptr + 1);
    delaylo = *delayptr;
    delayhi = *(delayptr + 1);

    /* Adjust the start time for the barycentric delay */
    
    outtime = intime - 
      LININTERP(outtime, delaytlo, delaythi, delaylo, delayhi);
    
    /* Now perform the barycentering */

    printf("Massaging the data...\n\n");
    printf("Amount Complete = 0%%");

    while ((numread = 
	   (*readrec_ptr)(infile, data, worklen, dispdt, numchan))){

      intimelo = totread * idata.dt;
      totread += numread;

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
      
      /* Work through each batch of topocentric data */
      
      for (ii = 0; ii < numread; ii++) {

	/* Calculate the topocentric time for the next point. */
      
	intimenext = intimelo + (ii + 1) * idata.dt;

	/* Set the delay pointers and variables */
	
	if (intimenext > delaythi){

	  /* Guess that the next delay we want is the next available */
	  
	  arrayoffset += 2;  /* Beware nasty NR zero-offset kludges! */
	  hunt(ttoa-1, numbarypts, intimenext, &arrayoffset);
	  arrayoffset--;
	  delaytimeptr = ttoa + arrayoffset;
	  delayptr = btoa + arrayoffset;
	  delaytlo = *delaytimeptr;
	  delaythi = *(delaytimeptr + 1);
	  delaylo = *delayptr;
	  delayhi = *(delayptr + 1);
	}

	/* Delay the topocentric time by the appropriate amount */
	
	outtimenext = intimenext - 
	  LININTERP(intimenext, delaytlo, delaythi, delaylo, delayhi);
	
	/* How much time does the data point cover? */
      
      
	/* Find the output points we will add data to.   */
	
	lopoint = outtime * rdt;
	hipoint = outtimenext * rdt;
	deltapoints = (outtimenext - outtime) * rdt;
	rdeltapoints = 1.0 / deltapoints;
	lobin = (unsigned long) lopoint;
	hibin = (unsigned long) hipoint;

	/* How many output points we will spread the data over? */
	
	numbins = hibin - lobin + 1;
	index = lobin - totwrote;
      
	/* Barycentric bin fits in one output bin */

	if (numbins == 1) {
	  outdata[index] += data[ii];

	/* Barycentric bin covers 2 output bins */

	} else if (numbins == 2) {
	  index = lobin - totwrote;
	  fact = (hibin - lopoint) * rdeltapoints;
	  outdata[index] += data[ii] * fact;
	  outdata[index+1] += data[ii] * (1.0 - fact);
	  index++;

	/* Barycentric bin covers 3 output bins */

	} else {
	  fact = ((lobin + 1) - lopoint) * rdeltapoints;
	  dtmp = idata.dt * rdeltapoints;
	  outdata[index] += data[ii] * fact;
	  outdata[index+1] += data[ii] * dtmp;
	  outdata[index+2] += data[ii] * (1.0 - (fact + dtmp));
	  index += 2;

	}

	/* Increment our time steps */

	intime = intimenext;
	outtime = outtimenext;
      }

      /* Write the latest chunk of data, but don't   */
      /* write more than cmd->numout points.         */

      wrote = index;
      if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	wrote = cmd->numout - totwrote;
      chkfwrite(outdata, sizeof(float), wrote, outfile);

      /* The following accounts for overlap effects between reads */
      
      dtmp = outdata[index];
      for (ii = 0; ii < wrote; ii++) {

	/* Find the max and min values */
	
	if (outdata[ii] > max)
	  max = outdata[ii];
	if (outdata[ii] < min)
	  min = outdata[ii];
	
	/* Use clever single pass mean and variance calculation */
	
	dev = outdata[ii] - avg;
	avg += dev / (totwrote + ii + 1);
	var += dev * (outdata[ii] - avg);
	
	/* Reset the output array */
	
	outdata[ii] = 0.0;
      }

      /* Set the lowest output  point to the last fractional point */

      outdata[0] = dtmp;

      /* Update the total number of points written */

      totwrote += wrote;

      /* Stop if we have written out all the data we need to */

      if (cmd->numoutP && (totwrote >= cmd->numout))
	break;
    }
  }

  /* Print simple stats and results */

  var /= (totwrote - 1);
  printf("\n\nDone.\n\nSimple statistics of the output data:\n");
  printf("  Number of points written:    %ld\n", totwrote);
  printf("  Maximum value:               %f\n", max);
  printf("  Minimum value:               %f\n", min);
  printf("  Mean value:                  %f\n", avg);
  printf("  Standard deviation:          %f\n", sqrt(var));
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
      for (ii = 0; ii < wlen2; ii++)
	outdata[ii] = 0.0;
    else
      for (ii = 0; ii < wlen2; ii++)
	outdata[ii] = avg;
    for (ii = 0; ii < (next_pow_of_two - totwrote) / wlen2; ii++)
      chkfwrite(outdata, sizeof(float), wlen2, outfile);
    chkfwrite(outdata, sizeof(float), \
	      (next_pow_of_two - totwrote) % wlen2, outfile);
  }

  /* Pad the file if needed to the requested length */

  if (cmd->numoutP && (cmd->numout > totwrote)) {
    for (ii = 0; ii < wlen2; ii++)
      outdata[ii] = avg;
    for (ii = 0; ii < (cmd->numout - totwrote) / wlen2; ii++)
      chkfwrite(outdata, sizeof(float), wlen2, outfile);
    chkfwrite(outdata, sizeof(float), 
	      (cmd->numout - totwrote) % wlen2, outfile);
  }

  /* Close the files and cleanup */

  fclose(infile);
  fclose(outfile);
  free(data);
  free(outdata);
  free(tobsf);
  free(dispdt);
  free(rootfilenm);
  free(outinfonm);
  free(datafilenm);
  if (idata.onoff)
    free(idata.onoff);
  if (!cmd->nobaryP) {
    free(bobsf);
    free(btoa);
    free(ttoa);
    free(voverc);
  }
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
  return chkfread(data, sizeof(float), \
		  (unsigned long) (numpts * numchan), file) / numchan;
}




