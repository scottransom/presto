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

static int read_floats(FILE *file, float *data, int numpts, int numchan);
static void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var);
static void update_infodata(infodata *idata, long datawrote, long padwrote, 
			    int *barybins, int numbarybins);
/* Update our infodata for barycentering and padding */

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE *infiles[MAXPATCHFILES], *outfile;
  float *outdata=NULL;
  double tdf=0.0, dtmp=0.0, barydispdt=0.0;
  double *dispdt, *tobsf=NULL, tlotoa=0.0, blotoa=0.0;
  double max=-9.9E30, min=9.9E30, var=0.0, avg=0.0;
  char obs[3], ephem[10], *datafilenm, *outinfonm;
  char rastring[50], decstring[50];
  int numfiles, numchan=1, newper=0, oldper=0, nummasked=0;
  int slen, numadded=0, numremoved=0, padding=0, *maskchans=NULL;
  long ii, numbarypts=0, worklen=65536;
  long numread=0, numtowrite=0, totwrote=0, datawrote=0;
  long padwrote=0, padtowrite=0, statnum=0;
  int numdiffbins=0, *diffbins=NULL, *diffbinptr=NULL;
  PKMB_tapehdr hdr;
  infodata idata;
  Cmdline *cmd;
  mask obsmask;

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
  printf("           Last Modification:  16 Dec, 2000\n\n");

  numfiles = cmd->argc;
  if (!cmd->pkmbP && !cmd->ebppP){
    char *root, *suffix;

    /* Split the filename into a rootname and a suffix */

    if (split_root_suffix(cmd->argv[0], &root, &suffix)==0){
      printf("\nThe input filename (%s) must have a suffix!\n\n", 
	     cmd->argv[0]);
      exit(1);
    }
    free(suffix);

    printf("Reading input data from '%s'.\n", cmd->argv[0]);
    printf("Reading information from '%s.inf'.\n\n", root);

    /* Read the info file if available */

    readinf(&idata, root);
    infiles[0] = chkfopen(cmd->argv[0], "rb");
    free(root);
  } else if (cmd->pkmbP){
    if (numfiles > 1)
      printf("Reading Parkes PKMB data from %d files:\n", numfiles);
    else
      printf("Reading Parkes PKMB data from 1 file:\n");
  } else if (cmd->ebppP){
    if (numfiles > 1)
      printf("Reading Effelsberg RBPP data from %d files:\n", numfiles);
    else
      printf("Reading Effelsberg RBPP data from 1 file:\n");
  }

  /* Open the raw data files */

  if (cmd->pkmbP || cmd->ebppP){
    for (ii=0; ii<numfiles; ii++){
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }
    printf("\n");
  }

  /* Determine the other file names and open the output data file */

  slen = strlen(cmd->outfile)+8;
  datafilenm = (char *)calloc(slen, 1);
  sprintf(datafilenm, "%s.dat", cmd->outfile);
  outfile = chkfopen(datafilenm, "wb");
  sprintf(idata.name, "%s", cmd->outfile);
  outinfonm = (char *)calloc(slen, 1);
  sprintf(outinfonm, "%s.inf", cmd->outfile);

  /* Read an input mask if wanted */

  if (cmd->maskfileP)
    read_mask(cmd->maskfile, &obsmask);
  else
    obsmask.numchan = obsmask.numint = 0;
  
  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {
    double dt, T;
    int ptsperblock;
    long long N;

    get_PKMB_file_info(infiles, numfiles, &N, &ptsperblock, &numchan, 
		       &dt, &T, 1);

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infiles[0]);
    rewind(infiles[0]);
    PKMB_hdr_to_inf(&hdr, &idata);
    PKMB_update_infodata(numfiles, &idata);
    idata.dm = cmd->dm;
    worklen = ptsperblock;
    if (cmd->maskfileP)
      maskchans = gen_ivect(idata.num_chan);

    /* Compare the size of the data to the size of output we request */

    if (cmd->numoutP) {
      dtmp = idata.N;
      idata.N = cmd->numout;
      writeinf(&idata);
      idata.N = dtmp;
    } else {
      cmd->numout = INT_MAX;
      writeinf(&idata);
    }
     
    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;
  }

  /* Set-up values if we are using the Effelsberg-Berkeley Pulsar Processor */
  /*   NOTE:  This code is not yet implemented.                             */

  if (cmd->ebppP) {

    /* Read the first header file and generate an infofile from it */

    /* OBS code for TEMPO */

    strcpy(obs, "EF");

    /* The number of data points to work with at a time */

    worklen = 1024;
  }

  /* Determine our initialization data if we do _not_ have Parkes */
  /* or Effelesberg data sets.                                    */

  if (!cmd->pkmbP && !cmd->ebppP) {

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

  printf("Writing output data to '%s'.\n", datafilenm);
  printf("Writing information to '%s'.\n\n", outinfonm);

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

  if (cmd->nobaryP) { /* Main loop if we are not barycentering... */

    /* Allocate our data array */
    
    outdata = gen_fvect(worklen);
    
    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");
    
    if (cmd->pkmbP)
      numread = read_PKMB(infiles, numfiles, outdata, worklen, 
			  dispdt, &padding, maskchans, &nummasked,
			  &obsmask);
    else
      numread = read_floats(infiles[0], outdata, worklen, numchan);

    while (numread){

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
      
      numtowrite = numread;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      chkfwrite(outdata, sizeof(float), numtowrite, outfile);
      totwrote += numtowrite;
      
      /* Update the statistics */
      
      if (!padding){
	for (ii = 0; ii < numtowrite; ii++)
	  update_stats(statnum + ii, outdata[ii], &min, &max, &avg, &var);
	statnum += numtowrite;
      }
      
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote == cmd->numout))
	break;

      if (cmd->pkmbP)
	numread = read_PKMB(infiles, numfiles, outdata, worklen, 
			    dispdt, &padding, maskchans, &nummasked,
			    &obsmask);
      else
	numread = read_floats(infiles[0], outdata, worklen, numchan);
     }
      
    datawrote = totwrote;

  } else { /* Main loop if we are barycentering... */

    double avgvoverc=0.0, maxvoverc=-1.0, minvoverc=1.0, *voverc=NULL;
    double *bobsf=NULL, *btoa=NULL, *ttoa=NULL;

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
    for (ii = 0 ; ii < numbarypts ; ii++){
      if (voverc[ii] > maxvoverc) maxvoverc = voverc[ii];
      if (voverc[ii] < minvoverc) minvoverc = voverc[ii];
      avgvoverc += voverc[ii];
    }
    avgvoverc /= numbarypts;
    free(voverc);
    blotoa = btoa[0];

    printf("   Insure you check the files tempoout_times.tmp and\n");
    printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n");
    printf("   Average topocentric velocity (c) = %.5g.\n", avgvoverc);
    printf("   Maximum topocentric velocity (c) = %.5g.\n", maxvoverc);
    printf("   Minimum topocentric velocity (c) = %.5g.\n\n", minvoverc);
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

    { /* Find the points where we need to add or remove bins */

      int oldbin=0, currentbin;
      double lobin, hibin, calcpt;

      numdiffbins = abs(NEAREST_INT(btoa[numbarypts-1])) + 1;
      diffbins = gen_ivect(numdiffbins);
      diffbinptr = diffbins;
      for (ii = 1 ; ii < numbarypts ; ii++){
	currentbin = NEAREST_INT(btoa[ii]);
	if (currentbin != oldbin){
	  if (currentbin > 0){
	    calcpt = oldbin + 0.5;
	    lobin = (ii-1) * TDT / idata.dt;
	    hibin = ii * TDT / idata.dt;
	  } else {
	    calcpt = oldbin - 0.5;
	    lobin = -((ii-1) * TDT / idata.dt);
	    hibin = -(ii * TDT / idata.dt);
	  }
	  while(fabs(calcpt) < fabs(btoa[ii])){
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
    
    /* Allocate our data array */
    
    outdata = gen_fvect(worklen);
    
    if (cmd->pkmbP)
      numread = read_PKMB(infiles, numfiles, outdata, worklen, 
			  dispdt, &padding, maskchans, &nummasked,
			  &obsmask);
    else
      numread = read_floats(infiles[0], outdata, worklen, numchan);
    
    while (numread){ /* Loop to read and write the data */
      int numwritten=0;
      
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
      /* OR write the amount of data up to cmd->numout or */
      /* the next bin that will be added or removed.      */
      
      numtowrite = abs(*diffbinptr) - datawrote;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      if (numtowrite > numread)
	numtowrite = numread;
      chkfwrite(outdata, sizeof(float), numtowrite, outfile);
      datawrote += numtowrite;
      totwrote += numtowrite;
      numwritten += numtowrite;
      
      /* Update the statistics */
      
      if (!padding){
	for (ii = 0; ii < numtowrite; ii++)
	  update_stats(statnum + ii, outdata[ii], &min, &max, 
		       &avg, &var);
	statnum += numtowrite;
      }
      
      if ((datawrote == abs(*diffbinptr)) &&
	  (totwrote < cmd->numout)){ /* Add/remove a bin */
	float favg;
	int skip, nextdiffbin;
	
	skip = numtowrite;
	
	do { /* Write the rest of the data after adding/removing a bin  */
	  
	  if (*diffbinptr > 0){
	    
	    /* Add a bin */
	    
	    favg = (float) avg;
	    chkfwrite(&favg, sizeof(float), 1, outfile);
	    numadded++;
	    totwrote++;
	  } else {
	    
	    /* Remove a bin */
	    
	    numremoved++;
	    datawrote++;
	    numwritten++;
	    skip++;
	  }
	  diffbinptr++;
	  
	  /* Write the part after the diffbin */
	  
	  numtowrite = numread - numwritten;
	  if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	    numtowrite = cmd->numout - totwrote;
	  nextdiffbin = abs(*diffbinptr) - datawrote;
	  if (numtowrite > nextdiffbin)
	    numtowrite = nextdiffbin;
	  chkfwrite(outdata + skip, sizeof(float), numtowrite, outfile);
	  numwritten += numtowrite;
	  datawrote += numtowrite;
	  totwrote += numtowrite;
	   
	  /* Update the statistics and counters */
	  
	  if (!padding){
	    for (ii = 0; ii < numtowrite; ii++)
	      update_stats(statnum + ii, outdata[skip + ii], 
			   &min, &max, &avg, &var);
	    statnum += numtowrite;
	  }
	  skip += numtowrite;

	   /* Stop if we have written out all the data we need to */
      
	  if (cmd->numoutP && (totwrote == cmd->numout)) 
	    break;
	} while (numwritten < numread);
      }
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote == cmd->numout))
	break;

      if (cmd->pkmbP)
	numread = read_PKMB(infiles, numfiles, outdata, worklen, 
			    dispdt, &padding, maskchans, &nummasked,
			    &obsmask);
      else
	numread = read_floats(infiles[0], outdata, worklen, numchan);
    }

    /* Free the arrays used in barycentering */

    free(bobsf);
    free(btoa);
    free(ttoa);
  }

  /* Calculate what the amount of padding we need  */

  if (cmd->numoutP && (cmd->numout > totwrote))
    padwrote = padtowrite = cmd->numout - totwrote;
  
  
  /* Write the new info file for the output data */

  if (!cmd->nobaryP) {
    idata.bary = 1;
    idata.mjd_i = (int) floor(blotoa - (barydispdt / SECPERDAY));
    idata.mjd_f = blotoa - (barydispdt / SECPERDAY) - idata.mjd_i;
  }
  update_infodata(&idata, totwrote, padtowrite, diffbins, numdiffbins);
  writeinf(&idata);

  /* Set the padded points equal to the average data point */

  if (idata.numonoff > 1){
    int jj, index, startpad, endpad;

    for (ii = 0; ii < worklen; ii++)
      outdata[ii] = avg;
    fclose(outfile);
    outfile = chkfopen(datafilenm, "rb+");
    for (ii=0; ii<idata.numonoff; ii++){
      index = 2 * ii;
      startpad = idata.onoff[index+1];
      if (ii==idata.numonoff-1)
	endpad = idata.N - 1;
      else
	endpad = idata.onoff[index+2];
      chkfseek(outfile, (startpad+1)*sizeof(float), SEEK_SET);
      padtowrite = endpad - startpad;
      for (jj=0; jj<padtowrite/worklen; jj++)
	chkfwrite(outdata, sizeof(float), worklen, outfile);
      chkfwrite(outdata, sizeof(float), padtowrite%worklen, outfile);
    }
  }
  free(outdata);

  /* Print simple stats and results */

  var /= (datawrote - 1);
  printf("\n\nDone.\n\nSimple statistics of the output data:\n");
  printf("             Data points written:  %ld\n", totwrote);
  if (padwrote)
    printf("          Padding points written:  %ld\n", padwrote);
  if (!cmd->nobaryP) {
    if (numadded)
      printf("    Bins added for barycentering:  %d\n", numadded);
    if (numremoved)
      printf("  Bins removed for barycentering:  %d\n", numremoved);
  }
  printf("           Maximum value of data:  %.2f\n", max);
  printf("           Minimum value of data:  %.2f\n", min);
  printf("              Data average value:  %.2f\n", avg);
  printf("         Data standard deviation:  %.2f\n", sqrt(var));
  printf("\n");

  /* Close the files and cleanup */

  if (cmd->maskfileP){
    free_mask(obsmask);
    free(maskchans);
  }
  for (ii=0; ii<numfiles; ii++)
    fclose(infiles[ii]);
  fclose(outfile);
  free(tobsf);
  free(dispdt);
  free(outinfonm);
  free(datafilenm);
  if (!cmd->nobaryP)
    free(diffbins);
  return (0);
}


static int read_floats(FILE *file, float *data, int numpts, \
		       int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
  /* Read the raw data and return numbar read */

  return chkfread(data, sizeof(float), (numpts * numchan), 
		  file) / numchan;
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


static void update_infodata(infodata *idata, long datawrote, long padwrote, 
			    int *barybins, int numbarybins)
/* Update our infodata for barycentering and padding */
{
  int ii, jj, index;

  idata->N = datawrote + padwrote;
  if (idata->numonoff==0){
    if (padwrote){
      idata->numonoff = 2;
      idata->onoff[0] = 0.0;
      idata->onoff[1] = datawrote-1;
      idata->onoff[2] = idata->N-1;
      idata->onoff[3] = idata->N-1;
    }
    return;
  }
  
  /* Determine the barycentric onoff bins (approximate) */

  if (numbarybins){
    int numadded=0, numremoved=0;

    ii = 1; /* onoff index    */
    jj = 0; /* barybins index */
    while (ii < idata->numonoff * 2){
      while (abs(barybins[jj]) <= idata->onoff[ii] &&
	     jj < numbarybins){
	if (barybins[jj] < 0)
	  numremoved++;
	else
	  numadded++;
	jj++;
      }
      idata->onoff[ii] += numadded - numremoved;
      ii++;
    }
  }

  /* Now cut off the extra onoff bins */

  for (ii=1, index=1; ii<=idata->numonoff; ii++, index+=2){
    if (idata->onoff[index-1] > idata->N - 1){
      idata->onoff[index-1] = idata->N - 1;
      idata->onoff[index] = idata->N - 1;
      break;
    }
    if (idata->onoff[index] > datawrote - 1){
      idata->onoff[index] = datawrote - 1;
      idata->numonoff = ii;
      if (padwrote){
	idata->numonoff++;
	idata->onoff[index+1] = idata->N - 1;
	idata->onoff[index+2] = idata->N - 1;
      }
      break;
    }
  }
}


