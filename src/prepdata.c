#include "presto.h"
#include "prepdata_cmd.h"
#include "multibeam.h"

/* This causes the binary orbit to be calculated about once every second */

#define TDT 0.00001525878906250

/* Some function definitions */

int (*readrec_ptr) (FILE * file, float *data, long numpts, \
		    double *dispdelays, long numchan);
int read_resid_rec(FILE * file, double *toa, double *obsf);
int read_floats(FILE * file, float *data, long numpts, \
		double *dispdelays, long numchan);


/* The main program */

int main(int argc, char *argv[])
/* This routine generates barycentric info for a data   */
/* file.  The file can contain many spectral channels.  */
/* It uses either the DE200 or the DE405 ephemeris.     */
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE *infile = NULL, *outfile;
  float *data, *outdata;
  double tdf = 0.0, dt = 0.0, tdtbins, tt1, tt2, ttdt, *dispdt;
  double *tobsf = NULL, *bobsf = NULL, tlotoa = 0.0, blotoa = 0.0;
  double *btl = NULL, *bth = NULL, bt, bdt, fact, dtmp;
  double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0, dev = 0.0;
  double *btoa = NULL, *ttoa = NULL, *voverc = NULL;
  char obs[3], ephem[10], datafilenm[120], filenm[120];
  char rastring[50], decstring[50];
  int numchan = 1, newper = 0, oldper = 0;
  long i, j, k, numbarypts = 0, worklen = 65536, next_pow_of_two = 0;
  long numread = 0, totread = 0, totwrote = 0, wrote = 0, wlen2;
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
  printf("                    13 Feb, 1998\n\n");

  /* Determine the output file name */

  sprintf(datafilenm, "%s.dat", cmd->outfile);

  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    if (!cmd->outfileP) {
      printf("You must enter an output file name if you are using the \n");
      printf("Parkes Multibeam format.  Exiting.\n\n");
      exit(1);
    }
    /* Open the raw input file */

    infile = chkfopen(cmd->argv[0], "rb");

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    chkfileseek(infile, 0L, sizeof(char), SEEK_SET);
    multibeam_hdr_to_inf(&hdr, &idata);
    sprintf(filenm, "%s.inf", cmd->outfile);
    sprintf(idata.name, "%s", cmd->outfile);
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
    /* The topocentric spacing of the observed data points (days) */

    dt = idata.dt / SECPERDAY;

    /* OBS code for TEMPO */

    strcpy(obs, "PK");

    /* The data collection routine to use */

    readrec_ptr = read_multibeam;

    /* The number of data points to work with at a time */

    numchan = idata.num_chan;
    worklen = DATLEN * 8 / numchan;

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (long) (dt * idata.N * 1.1 / TDT + 5.5);

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

    /* Determine the input file name */

    sprintf(filenm, "%s.raw", cmd->argv[0]);

    if (!cmd->outfileP){

      /* Determine the output file name */

      sprintf(datafilenm, "%s.dat", cmd->argv[0]);
    }

    /* Open the raw input file */

    infile = chkfopen(filenm, "rb");

    /* Open and read the infofile */

    readinf(&idata, cmd->argv[0]);

    /* The topocentric spacing of the observed data points (days) */

    dt = idata.dt / SECPERDAY;

    /* The data collection routine to use */

    readrec_ptr = read_floats;

    /* If we will be barycentering... */

    if (!cmd->nobaryP) {

      /* The number of topo to bary time points to generate with TEMPO */

      numbarypts = (long) (dt * idata.N * 1.1/ TDT + 5.5);

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
  /* Allocate and initialize some arrays and other information */

  data = gen_fvect(worklen);
  wlen2 = worklen * 2;
  outdata = gen_fvect(wlen2);
  for (i = 0; i < wlen2; i++)
    outdata[i] = 0.0;

  /* The topocentric epoch of the start of the data */

  tlotoa = (double) idata.mjd_i + idata.mjd_f;

  if (!strcmp(idata.band, "Radio") && (cmd->pkmbP || cmd->ebppP)) {

    /* The topocentric spacing between channels */

    tdf = idata.chan_wid;
    numchan = idata.num_chan;

    /* The topocentric observation frequencies */

    tobsf = gen_dvect(numchan);

    tobsf[0] = idata.freq;
    for (i = 0; i < numchan; i++) {
      tobsf[i] = tobsf[0] + i * tdf;
    }

    /* The dispersion delays (in days) */

    dispdt = gen_dvect(numchan);

    if (cmd->nobaryP) {

      /* Determine our dispersion time delays for each channel */

      for (i = 0; i < numchan; i++) {
	dispdt[i] = delay_from_dm(cmd->dm, tobsf[i]) / SECPERDAY;
      }

      /* The highest frequency channel gets no delay                   */
      /* All other delays are positive fractions of bin length (dt)    */

      dtmp = dispdt[numchan - 1];
      for (i = 0; i < numchan; i++) {
	dispdt[i] -= dtmp;
	dispdt[i] /= dt;
      }

    }

  } else {			/* For non-radio data */

    /* The topocentric observation frequencies */

    tobsf = gen_dvect(numchan);

    /* The dispersion delays (in days) */

    dispdt = gen_dvect(numchan);

    /* Set the observing frequency to infinity */

    tobsf[0] = 0.0;

    /* Our effective DM */

    cmd->dm = 0.0;

    /* Our dispersion time delay... */

    dispdt[0] = 0.0;

  }

  /* Main loop if we are not barycentering... */

  if (cmd->nobaryP) {

    /* Open our new output data file */

    outfile = chkfopen(datafilenm, "wb");

    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");

    wrote = worklen;

    do {

      numread = (*readrec_ptr) (infile, data, worklen, dispdt, numchan);

      /* Insure we don't write out too much data */

      if (cmd->numoutP && totwrote + wrote > cmd->numout)
	wrote = cmd->numout - totwrote;

      chkfwrite(data, sizeof(float), (unsigned long) wrote, outfile);
      totwrote += wrote;

      /* Print percent complete */

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
    for (i = 0 ; i <= numbarypts ; i++){
      ttoa[i] = tlotoa + (double) i * TDT;
    }

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(ttoa, btoa, voverc, numbarypts + 1, \
	       rastring, decstring, obs, ephem);

    printf("   Insure you check the files tempoout_times.tmp and\n");
    printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n\n");

    printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

    tlotoa = ttoa[0];
    blotoa = btoa[0];
    printf("Topocentric epoch (at data start) is:\n");
    printf("   %17.11f\n\n", tlotoa);
    printf("Barycentric epoch (infinite obs freq at data start) is:\n");
    printf("   %17.11f\n\n", blotoa);

    /* Determine the initial dispersion time delays for each channel */

    for (i = 0; i < numchan; i++) {
      bobsf[i] = doppler(tobsf[i], voverc[0]);
      dispdt[i] = delay_from_dm(cmd->dm, bobsf[i]) / SECPERDAY;
    }

    /* The highest frequency channel gets no delay                   */
    /* All other delays are positive fractions of bin length (dt)    */

    dtmp = dispdt[numchan - 1];
    for (i = 0; i < numchan; i++) {
      dispdt[i] -= dtmp;
      dispdt[i] /= dt;
    }

    /* Convert the rest of the barycentric TOAs */

    for (i = 0 ; i <= numbarypts ; i++) btoa[i] -= blotoa;

    /* The topocentric spacing of the points where we determined  */
    /* the barycentric corrections (days).                        */

    tdtbins = dt / TDT;

    /* Initialize some data */

    bt = 0.0;
    bdt = 0.0;
    tt1 = 0.0;
    ttdt = dt * worklen;
    tt2 = ttdt;
    wrote = worklen;

    /* Open our new output data file */

    outfile = chkfopen(datafilenm, "wb");

    printf("Massaging the data ...\n\n");
    printf("Amount Complete = 0%%");

    /* Now perform the barycentering */

    /* Blatantly stolen then subsequently modified */
    /* from Steve Eikenberry's bstretch.c          */

    do {

      numread = (*readrec_ptr) (infile, data, worklen, dispdt, numchan);

      /* Print percent complete */

      newper = (int) ((float) totwrote / idata.N * 100.0) + 1;
      if (newper > oldper) {
	printf("\rAmount Complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }
      /* Determine the barycentric time for the data segment */

      fact = modf(totread * tdtbins, &dtmp);
      totread += numread;
      btl = btoa + (long) (dtmp + 0.0000000001);
      bth = btl + 1;
      bdt = (*bth - *btl);
      bt = *btl + fact * bdt;
      bdt *= tdtbins;

      /* Loop over the newly read data segment */

      for (i = 0; i < worklen; i++, bt += bdt) {

	/* Get new barycentered dt if needed */

	if (bt > *bth) {
	  btl++;
	  bth++;
	  bdt = (*bth - *btl) * tdtbins;
	}
	/* Write a chunk of the output data set */

	if (bt >= tt2) {

	  /* If the barycentering compresses the data, don't */
	  /* write a bunch of zeros.                         */

	  if (totwrote && outdata[worklen - 1] == 0.0) {
	    wrote = 0;
	    while (outdata[wrote] != 0.0)
	      wrote++;
	  } else {
	    wrote = worklen;
	  }

	  /* Insure we don't write out too much data */

	  if (cmd->numoutP && totwrote + wrote > cmd->numout)
	    wrote = cmd->numout - totwrote;

	  chkfwrite(outdata, sizeof(float), (unsigned long) wrote, \
		    outfile);

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

	    outdata[j] = outdata[j + worklen];
	    outdata[j + worklen] = 0.0;
	  }

	  /* Increment our topocentric time steppers */

	  tt1 = tt2;
	  totwrote += wrote;
	  tt2 = (totwrote + worklen) * dt;
	}
	/* Determine if barycentric interval is less than, equal to, */
	/* or greater than the topocentric interval.                 */

	j = (long) floor((bt - tt1) / dt);
	k = (long) floor((bt - tt1 + bdt) / dt);

	if (k == j) {

	  /* Barycentric and topocentric bins are the same size */

	  outdata[j] += data[i];

	} else if (k - j == 1) {

	  /* Barycentric bins cover 2 topocentric bins */

	  fact = k * dt + tt1 - bt;
	  outdata[j] += data[i] * fact / bdt;
	  outdata[k] += data[i] * (1.0 - fact / bdt);

	} else {

	  /* Barycentric bins cover 3 topocentric bins */

	  fact = (j + 1) * dt + tt1 - bt;
	  outdata[j] += data[i] * fact / bdt;
	  outdata[j + 1] += data[i] * dt / bdt;
	  outdata[k] += data[i] * (1.0 - (fact + dt) / bdt);

	}

	if ((cmd->numoutP && totwrote >= cmd->numout) \
	    ||(wrote != worklen))
	  break;

      }

      if ((cmd->numoutP && totwrote >= cmd->numout) \
	  ||(wrote != worklen))
	break;

    } while (numread == worklen);


    /* Write the final few points if necessary */

    if (!(cmd->numoutP && totwrote >= cmd->numout)) {
      i = 0;
      while (outdata[i] != 0.0)
	i++;

      /* Insure we don't write out too much data */

      if (cmd->numoutP && totwrote + i > cmd->numout)
	i = cmd->numout - totwrote;

      chkfwrite(outdata, sizeof(float), (unsigned long) i, outfile);
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

  /* Update the info in our '.inf' file if this routine generated it */

  if (cmd->pkmbP) {
    idata.N = totwrote;
    idata.numonoff = 1;
    idata.onoff[0] = 0;
    idata.onoff[1] = totwrote - 1;
    if (!cmd->nobaryP) {
      idata.bary = 1;
      idata.mjd_i = (int) floor(blotoa);
      idata.mjd_f = blotoa - idata.mjd_i;
    }
    if (cmd->numoutP && cmd->numout > totwrote) {
      idata.N = cmd->numout;
      idata.onoff[2] = cmd->numout - 1;
      idata.onoff[3] = cmd->numout - 1;
    }
    if ((cmd->pad0P || cmd->padavgP) && next_pow_of_two != totwrote) {
      idata.N = next_pow_of_two;
      idata.onoff[2] = next_pow_of_two - 1;
      idata.onoff[3] = next_pow_of_two - 1;
    }
    writeinf(&idata);
  }

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
  if (idata.onoff) free(idata.onoff);
  if (!cmd->nobaryP) {
    free(bobsf);
    free(btoa);
    free(ttoa);
    free(voverc);
  }
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
  return chkfread(data, sizeof(float), \
		  (unsigned long) (numpts * numchan), file) / numchan;
}




