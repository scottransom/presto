#include "presto.h"
#include "prepdata_cmd.h"
#include "multibeam.h"

/* This causes the barycentric motion to be calculated once per second */

#define TDT 10.0

/* Some function definitions */

int (*readrec_ptr)(FILE * file, float *data, int numpts,
		   double *dispdelays, int numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, int numpts,
		double *dispdelays, int numchan);

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
  double tdf = 0.0, tdtbins, tt1, tt2, ttdt, *dispdt, rdt = 0.0, rbdt = 0.0;
  double *tobsf = NULL, *bobsf = NULL, tlotoa = 0.0, blotoa = 0.0;
  double *btl = NULL, *bth = NULL, bt, bdt, fact, dtmp = 0.0;
  double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0, dev = 0.0;
  double *btoa = NULL, *ttoa = NULL, *voverc = NULL, barydispdt = 0.0;
  char obs[3], ephem[10], *datafilenm, *rootfilenm, *outinfonm;
  char rastring[50], decstring[50], *cptr;
  int numchan = 1, newper = 0, oldper = 0, slen, hiindex;
  long i, j, k, numbarypts = 0, worklen = 65536, next_pow_of_two = 0;
  long numread = 0, totread = 0, wlen2, wrote = 0, totwrote = 0, bindex = 0;
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

  slen = strlen(cmd->outfile)+5;
  datafilenm = (char *)malloc(slen);
  sprintf(datafilenm, "%s.dat", cmd->outfile);
  outfile = chkfopen(datafilenm, "wb");
  sprintf(idata.name, "%s", cmd->outfile);
  outinfonm = (char *)malloc(slen);
  sprintf(outinfonm, "%s.inf", cmd->outfile);
  printf("Writing output data to '%s'.\n", datafilenm);
  printf("Writing information to '%s'.\n\n", outinfonm);

  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    chkfileseek(infile, 0L, sizeof(char), SEEK_SET);
    multibeam_hdr_to_inf(&hdr, &idata);
    idata.dm = cmd->dm;

    /* Compare the size of the data to the size of output we request */

    if (cmd->numoutP && cmd->numout != idata.N) {
      dtmp = idata.N;
      idata.N = cmd->numout;
    }

    /* Write a new info-file */

    writeinf(&idata);
    if (cmd->numoutP && cmd->numout == idata.N) {
      idata.N = dtmp;
    }

    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* The data collection routine to use */

    readrec_ptr = read_multibeam;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    worklen = DATLEN * 8 / numchan;

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5);

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

      numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5);

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
    for (i = 0; i < numchan; i++)
      tobsf[i] = tobsf[0] + i * tdf;

    /* The dispersion delays (in time bins) */

    dispdt = gen_dvect(numchan);

    if (cmd->nobaryP) {

      /* Determine our dispersion time delays for each channel */

      for (i = 0; i < numchan; i++)
	dispdt[i] = delay_from_dm(cmd->dm, tobsf[i]);

      /* The highest frequency channel gets no delay                 */
      /* All other delays are positive fractions of bin length (dt)  */

      dtmp = dispdt[numchan - 1];
      for (i = 0; i < numchan; i++)
	dispdt[i] = (dispdt[i] - dtmp) / idata.dt;
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
    for (i = 0; i < wlen2; i++)
      outdata[i] = 0.0;
    
    /* Open our new output data file */
    
    outfile = chkfopen(datafilenm, "wb");
    
    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");
    
    wrote = worklen;
    
    do {

      numread = (*readrec_ptr) (infile, data, worklen, dispdt, numchan);

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

      for (j = 0; j < wrote; j++) {

	/* Find the max and min values */
	
	if (data[j] > max)
	  max = data[j];
	if (data[j] < min)
	  min = data[j];
	
	/* Use clever single pass mean and variance calculation */
	
	dev = data[j] - avg;
	avg += dev / (totwrote + j + 1);
	var += dev * (data[j] - avg);
      }

      if (cmd->numoutP && (totwrote >= cmd->numout))
	break;
	
     } while (numread == worklen);
      
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
    btoa = gen_dvect(numbarypts + 1);
    ttoa = gen_dvect(numbarypts + 1);
    voverc = gen_dvect(numbarypts + 1);
    for (i = 0 ; i <= numbarypts ; i++)
      ttoa[i] = tlotoa + TDT * i / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(ttoa, btoa, voverc, numbarypts + 1, \
	       rastring, decstring, obs, ephem);
    blotoa = btoa[0];

    printf("   Insure you check the files tempoout_times.tmp and\n");
    printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n\n");

    printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

    /* Determine the initial dispersion time delays for each channel */

    for (i = 0; i < numchan; i++) {
      bobsf[i] = doppler(tobsf[i], voverc[0]);
      dispdt[i] = delay_from_dm(cmd->dm, bobsf[i]);
    }

    /* The highest frequency channel gets no delay                   */
    /* All other delays are positive fractions of bin length (dt)    */

    barydispdt = dispdt[numchan - 1];
    for (i = 0; i < numchan; i++)
      dispdt[i] = (dispdt[i] - barydispdt) / idata.dt;
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
    for (i = 0; i < wlen2; i++)
      outdata[i] = 0.0;
    rdt = 1.0 / idata.dt;

    /* Convert the barycentric TOAs to seconds from the beginning */

    for (i = 0 ; i <= numbarypts ; i++)
      btoa[i] = (btoa[i] - blotoa) * SECPERDAY;

    /* The length of a time bin in terms of the topocentric spacing   */
    /* of the points where we determined the barycentric corrections. */

    tdtbins = idata.dt / TDT;

    /* Initialize some data */

    bt = 0.0;
    bdt = 0.0;
    ttdt = idata.dt * worklen;
    tt1 = 0.0;
    tt2 = ttdt;
    wrote = worklen;
    hiindex = 0;

    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");

    /* Now perform the barycentering */

    /* Blatantly stolen then subsequently modified */
    /* from Steve Eikenberry's bstretch.c          */

    do {

      numread = (*readrec_ptr) (infile, data, worklen, dispdt, numchan);

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

      /* Determine the barycentric time for the data  */
      /* segment using linear interpolation.          */

      dtmp = totread * tdtbins;
      bindex = (long)(dtmp + DBLCORRECT); 
      fact = dtmp - bindex;
      totread += numread;
      btl = btoa + bindex;
      bth = btl + 1;
      bdt = (*bth - *btl);
      bt = *btl + fact * bdt;
      bdt *= tdtbins;
      rbdt = 1.0 / bdt;

      /* Loop over the newly read data segment */

      for (i = 0; i < worklen; i++, bt += bdt) {

	/* Get new barycentered dt if needed */

	if (bt > *bth) {
	  btl++;
	  bth++;
	  bdt = (*bth - *btl) * tdtbins;
	  rbdt = 1.0 / bdt;
	}

	/* Write a chunk of the output data set */

	if (bt >= tt2) {

	  /* If the barycentering compresses the data, don't */
	  /* write a bunch of zeros.  Also, don't write more */
	  /* than cmd->numout points.                        */

	  wrote = hiindex;
	  if (cmd->numoutP && (totwrote + wrote) > cmd->numout)
	    wrote = cmd->numout - totwrote;

	  chkfwrite(outdata, sizeof(float), wrote, outfile);

	  /* The following accounts for overlap effects between reads */

	  for (j = 0; j < wrote; j++) {

	    /* Find the max and min values */

	    if (outdata[j] > max)
	      max = outdata[j];
	    if (outdata[j] < min)
	      min = outdata[j];

	    /* Use clever single pass mean and variance calculation */

	    dev = outdata[j] - avg;
	    avg += dev / (totwrote + j + 1);
	    var += dev * (outdata[j] - avg);

	    /* Move the output array around for the next pass */

	    outdata[j] = outdata[j + wrote];
	  }

	  /* Re-set the high end of the buffer */

	  for (j = wrote; j < wlen2; j++) 
	    outdata[j] = 0.0;

	  /* Increment our topocentric time steppers */

	  tt1 = tt2;
	  totwrote += wrote;
	  tt2 = (totread + worklen) * idata.dt;
	  hiindex = 0;
	}

	/* Determine if barycentric interval is less than, equal to, */
	/* or greater than the topocentric interval.                 */

	j = (long) ((bt - tt1) * rdt + DBLCORRECT);       /* Bin start */
	k = (long) ((bt - tt1 + bdt) * rdt + DBLCORRECT); /* Bin end   */
	hiindex = k;

	/* Barycentric bin fits in one topocentric bin */

	if (k == j) {
	  outdata[j] += data[i];

	/* Barycentric bin covers 2 topocentric bins */

	} else if (k - j == 1) {
	  fact = ((k * idata.dt + tt1) - bt) * rbdt;
	  outdata[j] += data[i] * fact;
	  outdata[k] += data[i] * (1.0 - fact);

	/* Barycentric bin covers 3 topocentric bins */

	} else {
	  fact = (((j+1) * idata.dt + tt1) - bt) * rbdt;
	  dtmp = idata.dt * rbdt;
	  outdata[j] += data[i] * fact;
	  outdata[j+1] += data[i] * dtmp;
	  outdata[k] += data[i] * (1.0 - (fact + dtmp));
	}

	if (cmd->numoutP && totwrote >= cmd->numout)
	  break;

      }

      if (cmd->numoutP && totwrote >= cmd->numout)
	break;

    } while (numread == worklen);

    /* Write the final few points if necessary */

    if (!(cmd->numoutP && totwrote >= cmd->numout)) {
      i = wlen2 - 1;
      while (outdata[i] == 0.0)
	i--;
      if (i < 0) i = 0;

      /* Insure we don't write out too much data */

      if (cmd->numoutP && totwrote + i > cmd->numout)
	i = cmd->numout - totwrote;

      chkfwrite(outdata, sizeof(float), i, outfile);
      totwrote += i;
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
    idata.numonoff = 2;
    idata.onoff[0] = 0;
    idata.onoff[1] = totwrote - 1;
    idata.onoff[2] = cmd->numout - 1;
    idata.onoff[3] = cmd->numout - 1;
    printf("NOTE:  If there were on-off pairs in the input info file,\n");
    printf("       the output info file pairs are now incorrect!\n\n");
  }
  if ((cmd->pad0P || cmd->padavgP) && (next_pow_of_two != totwrote)) {
    idata.N = next_pow_of_two;
    idata.numonoff = 2;
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
    if (cmd->pad0P) {
      for (i = 0; i < wlen2; i++)
	outdata[i] = 0.0;
    } else {
      for (i = 0; i < wlen2; i++)
	outdata[i] = avg;
    }
    for (i = 0; i < (next_pow_of_two - totwrote) / wlen2; i++) {
      chkfwrite(outdata, sizeof(float), (unsigned long) wlen2, outfile);
    }
    totwrote += i * wlen2;
    chkfwrite(outdata, sizeof(float), \
	      (unsigned long) (next_pow_of_two - totwrote), outfile);
  }
  /* Pad the file if needed to the requested length */

  if (cmd->numoutP && cmd->numout > totwrote) {
    for (i = 0; i < wlen2; i++)
      outdata[i] = avg;
    for (i = 0; i < (cmd->numout - totwrote) / wlen2; i++) {
      chkfwrite(outdata, sizeof(float), (unsigned long) wlen2, outfile);
    }
    totwrote += i * wlen2;
    chkfwrite(outdata, sizeof(float), \
	      (unsigned long) (cmd->numout - totwrote), outfile);
  }
  /* Close the files and cleanup */

  fclose(infile);
  fclose(outfile);
  free(data);
  free(outdata);
  free(tobsf);
  free(dispdt);
  free(rootfilenm);
  free(datafilenm);
  free(outinfonm);
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




