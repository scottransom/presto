#include <limits.h>
#include "presto.h"
#include "prepsubband_cmd.h"
#include "multibeam.h"

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

static void write_data(FILE *outfiles[], int numfiles, float **outdata, 
		       int startpoint, int numtowrite);
static void write_padding(FILE *outfiles[], int numfiles, float value, 
			  int numtowrite);
static int get_data(FILE *infiles[], int numfiles, float **outdata, 
		    int numchan, int blocklen, int blocksperread, 
		    mask *obsmask, double *dispdts, int **offsets, 
		    int *padding);
static void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var);
static void update_infodata(infodata *idata, int datawrote, int padwrote, 
			    int *barybins, int numbarybins);
static void print_percent_complete(int current, int number);

/* From CLIG */
static Cmdline *cmd;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE **infiles, **outfiles;
  float **outdata;
  double dtmp, *dms, avgdm=0.0, maxdm;
  double *dispdt, tlotoa=0.0, blotoa=0.0;
  double max=-9.9E30, min=9.9E30, var=0.0, avg=0.0;
  double *btoa=NULL, *ttoa=NULL, avgvoverc=0.0;
  char obs[3], ephem[10], rastring[50], decstring[50];
  int numinfiles, numchan=1, totnumtowrite, **offsets;
  int ii, jj, numadded=0, numremoved=0, padding=0;
  int numbarypts=0, blocklen=0, blocksperread=0, worklen=0;
  int numread=0, numtowrite=0, totwrote=0, datawrote=0;
  int padwrote=0, padtowrite=0, statnum=0;
  int numdiffbins=0, *diffbins=NULL, *diffbinptr=NULL;
  char *datafilenm;
  PKMB_tapehdr hdr;
  infodata idata;
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
  printf("          Pulsar Subband De-dispersion Routine\n");
  printf("                 by Scott M. Ransom\n");
  printf("            Last Modification:  18 Mar, 2001\n\n");

  numinfiles = cmd->argc;
  if (cmd->pkmbP){
    if (numinfiles > 1)
      printf("Reading Parkes PKMB data from %d files:\n", numinfiles);
    else
      printf("Reading Parkes PKMB data from 1 file:\n");
  } else if (cmd->ebppP){
    if (numinfiles > 1)
      printf("Reading Effelsberg RBPP data from %d files:\n", numinfiles);
    else
      printf("Reading Effelsberg RBPP data from 1 file:\n");
  }

  /* Open the raw data files */

  infiles = (FILE **)malloc(numinfiles * sizeof(FILE *));
  for (ii=0; ii<numinfiles; ii++){
    printf("  '%s'\n", cmd->argv[ii]);
    infiles[ii] = chkfopen(cmd->argv[ii], "rb");
  }
  if (!cmd->numoutP)
    cmd->numout = INT_MAX;

  /* Determine the output file names and open them */

  datafilenm = (char *)calloc(strlen(cmd->outfile)+20, 1);
  dms = gen_dvect(cmd->numdms);
  outfiles = (FILE **)malloc(cmd->numdms * sizeof(FILE *));
  printf("\nWriting output data to:\n");
  for (ii=0; ii<cmd->numdms; ii++){
    dms[ii] = cmd->lodm + ii * cmd->dmstep;
    avgdm += dms[ii];
    sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
    outfiles[ii] = chkfopen(datafilenm, "wb");
    printf("   '%s'\n", datafilenm);
  }
  avgdm /= cmd->numdms;
  maxdm = dms[cmd->numdms-1];

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

    printf("\nPKMB input file information:\n");
    get_PKMB_file_info(infiles, numinfiles, &N, &ptsperblock, &numchan, 
		       &dt, &T, 1);
    chkfread(&hdr, 1, HDRLEN, infiles[0]);
    rewind(infiles[0]);
    PKMB_hdr_to_inf(&hdr, &idata);
    PKMB_update_infodata(numinfiles, &idata);
    idata.dm = avgdm;
    blocklen = ptsperblock;
    blocksperread = ((int)(delay_from_dm(maxdm, idata.freq)/dt) 
		     / ptsperblock + 1);
    worklen = blocklen * blocksperread;
    strcpy(obs, "PK");  /* OBS code for TEMPO */

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
  }

  /* Set-up values if we are using the Effelsberg-Berkeley Pulsar Processor */
  /*   NOTE:  This code is not yet implemented.                             */

  if (cmd->ebppP) {
    strcpy(obs, "EF");  /* OBS code for TEMPO */
  }

  tlotoa = idata.mjd_i + idata.mjd_f;  /* Topocentric epoch */

  if (cmd->numoutP)
    totnumtowrite = cmd->numout;
  else
    totnumtowrite = (int) idata.N;

  if (cmd->nobaryP) { /* Main loop if we are not barycentering... */
    
    /* Dispersion delays (in bins).  The high freq gets no delay   */
    /* All other delays are positive fractions of bin length (dt)  */
    
    dispdt = subband_search_delays(numchan, cmd->numsub, avgdm, 
				   idata.freq, idata.chan_wid, 0.0);
    for (ii=0; ii<numchan; ii++)
      dispdt[ii] /= idata.dt;
    
    /* The subband dispersion delays (see note above) */

    offsets = gen_imatrix(cmd->numdms, cmd->numsub);
    for (ii=0; ii<cmd->numdms; ii++){
      double *subdispdt;

      subdispdt = subband_delays(numchan, cmd->numsub, dms[ii], 
				 idata.freq, idata.chan_wid, 0.0);
      dtmp = subdispdt[cmd->numsub-1];
      for (jj=0; jj<cmd->numsub; jj++)
	offsets[ii][jj] = (int)((subdispdt[jj] - dtmp) / idata.dt + 0.5);
      free(subdispdt);
    }

    /* Allocate our data array and start getting data */
    
    outdata = gen_fmatrix(cmd->numsub, worklen);
    numread = get_data(infiles, numinfiles, outdata, 
		       numchan, blocklen, blocksperread, 
		       &obsmask, dispdt, offsets, &padding);
    printf("De-dispersing using:\n");
    printf("     Subbands = %d\n", cmd->numsub);
    printf("   Average DM = %.7g\n\n", avgdm);
    
    while (numread==worklen){

      print_percent_complete(totwrote, totnumtowrite);

      /* Write the latest chunk of data, but don't   */
      /* write more than cmd->numout points.         */
      
      numtowrite = numread;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
      totwrote += numtowrite;
      
      /* Update the statistics */
      
      if (!padding){
	for (ii=0; ii<numtowrite; ii++)
	  update_stats(statnum+ii, outdata[0][ii], &min, &max, &avg, &var);
	statnum += numtowrite;
      }
      
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote == cmd->numout))
	break;
      
      numread = get_data(infiles, numinfiles, outdata, 
			 numchan, blocklen, blocksperread, 
			 &obsmask, dispdt, offsets, &padding);
    }
    datawrote = totwrote;

  } else { /* Main loop if we are barycentering... */
    double maxvoverc=-1.0, minvoverc=1.0, *voverc=NULL;

    /* What ephemeris will we use?  (Default is DE200) */

    if (cmd->de405P)
      strcpy(ephem, "DE405");
    else
      strcpy(ephem, "DE200");

    /* Define the RA and DEC of the observation */
    
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Allocate some arrays */

    btoa = gen_dvect(numbarypts);
    ttoa = gen_dvect(numbarypts);
    voverc = gen_dvect(numbarypts);
    for (ii=0; ii<numbarypts; ii++)
      ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

    /* Call TEMPO for the barycentering */

    printf("Generating barycentric corrections...\n");
    barycenter(ttoa, btoa, voverc, numbarypts, \
	       rastring, decstring, obs, ephem);
    for (ii=0; ii<numbarypts; ii++){
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
    printf("De-dispersing and barycentering using:\n");
    printf("     Subbands = %d\n", cmd->numsub);
    printf("   Average DM = %.7g\n\n", avgdm);

    /* Dispersion delays (in bins).  The high freq gets no delay   */
    /* All other delays are positive fractions of bin length (dt)  */
    
    dispdt = subband_search_delays(numchan, cmd->numsub, avgdm, 
				   idata.freq, idata.chan_wid, avgvoverc);
    for (ii=0; ii<numchan; ii++)
      dispdt[ii] /= idata.dt;
    
    /* The subband dispersion delays (see note above) */

    offsets = gen_imatrix(cmd->numdms, cmd->numsub);
    for (ii=0; ii<cmd->numdms; ii++){
      double *subdispdt;

      subdispdt = subband_delays(numchan, cmd->numsub, dms[ii], 
				 idata.freq, idata.chan_wid, avgvoverc);
      dtmp = subdispdt[cmd->numsub-1];
      for (jj=0; jj<cmd->numsub; jj++)
	offsets[ii][jj] = (int)((subdispdt[jj] - dtmp) / idata.dt + 0.5);
      free(subdispdt);
    }

    /* Convert the bary TOAs to differences from the topo TOAs in */
    /* units of bin length (dt) rounded to the nearest integer.   */

    dtmp = (btoa[0] - ttoa[0]);
    for (ii=0; ii<numbarypts; ii++)
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

    outdata = gen_fmatrix(cmd->numsub, worklen);
    numread = get_data(infiles, numinfiles, outdata, 
		       numchan, blocklen, blocksperread, 
		       &obsmask, dispdt, offsets, &padding);
    
    while (numread==worklen){ /* Loop to read and write the data */
      int numwritten=0;

      print_percent_complete(totwrote, totnumtowrite);

      /* Simply write the data if we don't have to add or */
      /* remove any bins from this batch.                 */
      /* OR write the amount of data up to cmd->numout or */
      /* the next bin that will be added or removed.      */
      
      numtowrite = abs(*diffbinptr) - datawrote;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      if (numtowrite > numread)
	numtowrite = numread;
      write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
      datawrote += numtowrite;
      totwrote += numtowrite;
      numwritten += numtowrite;
      
      /* Update the statistics */
      
      if (!padding){
	for (ii = 0; ii < numtowrite; ii++)
	  update_stats(statnum + ii, outdata[0][ii], &min, &max, 
		       &avg, &var);
	statnum += numtowrite;
      }
      
      if ((datawrote == abs(*diffbinptr)) &&
	  (numwritten != numread) &&
	  (totwrote < cmd->numout)){ /* Add/remove a bin */
	int skip, nextdiffbin;
	
	skip = numtowrite;
	
	do { /* Write the rest of the data after adding/removing a bin  */
	  
	  if (*diffbinptr > 0){
	    /* Add a bin */
	    write_padding(outfiles, cmd->numdms, avg, 1);
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
	  write_data(outfiles, cmd->numdms, outdata, skip, numtowrite);
	  numwritten += numtowrite;
	  datawrote += numtowrite;
	  totwrote += numtowrite;
	   
	  /* Update the statistics and counters */
	  
	  if (!padding){
	    for (ii=0; ii<numtowrite; ii++)
	      update_stats(statnum+ii, outdata[0][skip+ii], 
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

      numread = get_data(infiles, numinfiles, outdata, 
			 numchan, blocklen, blocksperread, 
			 &obsmask, dispdt, offsets, &padding);
    }
  }

  /* Calculate what the amount of padding we need  */

  if (cmd->numoutP && (cmd->numout > totwrote))
    padwrote = padtowrite = cmd->numout - totwrote;
  
  /* Write the new info file for the output data */

  update_infodata(&idata, totwrote, padtowrite, diffbins, numdiffbins);
  for (ii=0; ii<cmd->numdms; ii++){
    idata.dm = dms[ii];
    if (!cmd->nobaryP) {
      double baryepoch, barydispdt, baryhifreq;

      baryhifreq = idata.freq+(numchan-1)*idata.chan_wid;
      barydispdt = delay_from_dm(dms[ii], doppler(baryhifreq, avgvoverc));
      baryepoch = blotoa - (barydispdt / SECPERDAY);
      idata.bary = 1;
      idata.mjd_i = (int) floor(baryepoch);
      idata.mjd_f = baryepoch - idata.mjd_i;
    }
    sprintf(idata.name, "%s_DM%.2f", cmd->outfile, dms[ii]);
    writeinf(&idata);
  }

  /* Set the padded points equal to the average data point */

  if (idata.numonoff > 1){
    int index, startpad, endpad;

    for (ii=0; ii<cmd->numdms; ii++){
      fclose(outfiles[ii]);
      sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
      outfiles[ii] = chkfopen(datafilenm, "rb+");
    }
    for (ii=0; ii<idata.numonoff; ii++){
      index = 2 * ii;
      startpad = idata.onoff[index+1];
      if (ii==idata.numonoff-1)
	endpad = idata.N - 1;
      else
	endpad = idata.onoff[index+2];
      for (jj=0; jj<cmd->numdms; jj++)
	chkfseek(outfiles[jj], (startpad+1)*sizeof(float), SEEK_SET);
      padtowrite = endpad - startpad;
      write_padding(outfiles, cmd->numdms, avg, padtowrite);
    }
  }

  /* Print simple stats and results */

  var /= (datawrote - 1);
  print_percent_complete(1, 1);
  printf("\n\nDone.\n\nSimple statistics of the output data:\n");
  printf("             Data points written:  %d\n", totwrote);
  if (padwrote)
    printf("          Padding points written:  %d\n", padwrote);
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

  if (cmd->maskfileP)
    free_mask(obsmask);
  for (ii=0; ii<numinfiles; ii++)
    fclose(infiles[ii]);
  free(infiles);
  for (ii=0; ii<cmd->numdms; ii++)
    fclose(outfiles[ii]);
  free(outdata[0]);
  free(outdata);
  free(outfiles);
  free(dms);
  free(dispdt);
  free(offsets[0]);
  free(offsets);
  free(datafilenm);
  if (!cmd->nobaryP){
    free(btoa);
    free(ttoa);
    free(diffbins);
  }
  return (0);
}

static void write_data(FILE *outfiles[], int numfiles, float **outdata, 
		       int startpoint, int numtowrite)
{
  int ii;

  for (ii=0; ii<numfiles; ii++)
    chkfwrite(outdata[ii]+startpoint, sizeof(float), 
	      numtowrite, outfiles[ii]);
}


static void write_padding(FILE *outfiles[], int numfiles, float value, 
			  int numtowrite)
{
  int ii;

  if (numtowrite<=0){
    return;
  } else if (numtowrite==1){
    for (ii=0; ii<numfiles; ii++)
      chkfwrite(&value, sizeof(float), 1, outfiles[ii]);
  } else {
    int maxatonce=8192, veclen, jj;
    float *buffer;   
    veclen = (numtowrite > maxatonce) ? maxatonce : numtowrite;
    buffer = gen_fvect(veclen);
    for (ii=0; ii<veclen; ii++)
      buffer[ii] = value;
    if (veclen==numtowrite){
      for (ii=0; ii<numfiles; ii++)
	chkfwrite(buffer, sizeof(float), veclen, outfiles[ii]);
    } else {
      for (ii=0; ii<=numtowrite/veclen; ii++){
	for (jj=0; jj<numfiles; jj++)
	  chkfwrite(buffer, sizeof(float), veclen, outfiles[jj]);
      }
      for (jj=0; jj<numfiles; jj++)
	chkfwrite(buffer, sizeof(float), numtowrite%veclen, outfiles[jj]);
    }
    free(buffer);
  }
}


static int get_data(FILE *infiles[], int numfiles, float **outdata, 
		    int numchan, int blocklen, int blocksperread, 
		    mask *obsmask, double *dispdts, int **offsets, 
		    int *padding)
{
  static int firsttime=1, worklen, *maskchans=NULL, blocksize;
  static float *tempzz, *data1, *data2; 
  static float *currentdata, *lastdata;
  int totnumread=0, numread, ii, jj, tmppad=0, nummasked;
  
  if (firsttime){
    if (cmd->maskfileP)
      maskchans = gen_ivect(numchan);
    worklen = blocklen * blocksperread;
    blocksize = blocklen * cmd->numsub;
    data1 = gen_fvect(cmd->numsub * worklen);
    data2 = gen_fvect(cmd->numsub * worklen);
    currentdata = data1;
    lastdata = data2;
    if (cmd->pkmbP){
      for (ii=0; ii<blocksperread; ii++){
	numread = read_PKMB_subbands(infiles, numfiles, 
				     currentdata+ii*blocksize, 
				     dispdts, cmd->numsub, 0, &tmppad, 
				     maskchans, &nummasked, obsmask);
	if (numread!=blocklen){
	  for (jj=ii*blocksize; jj<(ii+1)*blocksize; jj++)
	    currentdata[jj] = 0.0;
	}
	if (tmppad) 
	  *padding = 1;
      }
    }
    SWAP(currentdata, lastdata);
    firsttime = 0;
  }
  if (cmd->pkmbP){
    for (ii=0; ii<blocksperread; ii++){
      numread = read_PKMB_subbands(infiles, numfiles, 
				   currentdata+ii*blocksize, 
				   dispdts, cmd->numsub, 0, &tmppad, 
				   maskchans, &nummasked, obsmask);
      totnumread += numread;
      if (numread!=blocklen){
	for (jj=ii*blocksize; jj<(ii+1)*blocksize; jj++)
	  currentdata[jj] = 0.0;
      }
      if (tmppad) 
	*padding = 1;
    }
  }
  for (ii=0; ii<cmd->numdms; ii++){
    float_dedisp(currentdata, lastdata, worklen, cmd->numsub, 
		 offsets[ii], 0.0, outdata[ii]);
  }
  SWAP(currentdata, lastdata);
  if (totnumread != worklen){
    if (cmd->maskfileP)
      free(maskchans);
    free(data1);
    free(data2);
  }
  return totnumread;
}


static void print_percent_complete(int current, int number)
{
  static int newper=0, oldper=-1;
 
  newper = (int) (current / (float)(number) * 100.0);
  if (newper < 0) newper = 0;
  if (newper > 100) newper = 100;
  if (newper > oldper) {
    printf("\rAmount complete = %3d%%", newper);
    fflush(stdout);
    oldper = newper;
  }
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


static void update_infodata(infodata *idata, int datawrote, int padwrote, 
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


