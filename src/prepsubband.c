#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "prepsubband_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP)

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
			    int *barybins, int numbarybins, int downsamp);
static void print_percent_complete(int current, int number);

/* From CLIG */
static Cmdline *cmd;
static BPP_ifs bppifs=SUMIFS;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE **infiles, **outfiles;
  float **outdata;
  double dtmp, *dms, avgdm=0.0, maxdm, dsdt=0;
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
  numinfiles = cmd->argc;

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("          Pulsar Subband De-dispersion Routine\n");
  printf("                 by Scott M. Ransom\n\n");

  if (!RAWDATA){
    char *root, *suffix;
    /* Split the filename into a rootname and a suffix */
    if (split_root_suffix(cmd->argv[0], &root, &suffix)==0){
      printf("\nThe input filename (%s) must have a suffix!\n\n", 
	     cmd->argv[0]);
      exit(1);
    } else {
      if (strcmp(suffix, "bcpm1")==0 || 
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
      } else {
	printf("\nCannot determine the format of the input files '%s'...\n\n", 
	       cmd->argv[0]);
	exit(1);
      }
      free(root);
      free(suffix);
    }
  }

  if (cmd->pkmbP){
    if (numinfiles > 1)
      printf("Reading Parkes PKMB data from %d files:\n", numinfiles);
    else
      printf("Reading Parkes PKMB data from 1 file:\n");
  } else if (cmd->bcpmP){
    if (numinfiles > 1)
      printf("Reading Green Bank BCPM data from %d files:\n", numinfiles);
    else
      printf("Reading Green Bank BCPM data from 1 file:\n");
  } else if (cmd->wappP){
    if (numinfiles > 1)
      printf("Reading Arecibo WAPP data from %d files:\n", numinfiles);
    else
      printf("Reading Arecibo WAPP data from 1 file:\n");
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

  if (cmd->maskfileP){
    read_mask(cmd->maskfile, &obsmask);
    printf("Read mask information from '%s'\n", cmd->maskfile);
  } else {
    obsmask.numchan = obsmask.numint = 0;
  }
  
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
    dsdt = cmd->downsamp * idata.dt;
    idata.dm = avgdm;
    blocklen = ptsperblock;
    blocksperread = ((int)(delay_from_dm(maxdm, idata.freq)/dsdt) 
		     / ptsperblock + 1);
    worklen = blocklen * blocksperread;
    strcpy(obs, "PK");  /* OBS code for TEMPO */

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
  }

  /* Set-up values if we are using the Berkeley-Caltech */
  /* Pulsar Machine (or BPP) format.                    */

  if (cmd->bcpmP) {
    double dt, T;
    int ptsperblock;
    long long N;

    printf("\nBCPM input file information:\n");
    get_BPP_file_info(infiles, numinfiles, &N, &ptsperblock, &numchan, 
		      &dt, &T, &idata, 1);
    BPP_update_infodata(numinfiles, &idata);
    dsdt = cmd->downsamp * idata.dt;
    idata.dm = avgdm;
    blocklen = ptsperblock;
    blocksperread = ((int)(delay_from_dm(maxdm, idata.freq)/dsdt) 
		     / ptsperblock + 1);
    worklen = blocklen * blocksperread;

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
    */
    /* The following is for the Green Bank Telescope */
    strcpy(obs, "GB");

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
  }

  /* Set-up values if we are using the Arecobo WAPP */

  if (cmd->wappP) {
    double dt, T;
    int ptsperblock;
    long long N;

    printf("\nWAPP input file information:\n");
    get_WAPP_file_info(infiles, cmd->numwapps, numinfiles, cmd->clip,
		       &N, &ptsperblock, &numchan, 
		       &dt, &T, &idata, 1);
    WAPP_update_infodata(numinfiles, &idata);
    dsdt = cmd->downsamp * idata.dt;
    idata.dm = avgdm;
    blocklen = ptsperblock;
    blocksperread = ((int)(delay_from_dm(maxdm, idata.freq)/dsdt) 
		     / ptsperblock + 1);
    worklen = blocklen * blocksperread;

    /* OBS code for TEMPO */
    /* The following is for Arecibo */
    strcpy(obs, "AO");

    /* The number of topo to bary time points to generate with TEMPO */

    numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
  }

  if (cmd->numsub > numchan){
    printf("Warning:  The number of requested subbands (%d) is larger than the number of channels (%d).\n",
	   cmd->numsub, numchan);
    printf("          Re-setting the number of subbands to %d.\n\n", numchan);
    cmd->numsub = numchan;
  }

  if (blocklen % cmd->downsamp){
    printf("Error:  The downsample factor (%d) must be a factor of the\n",
	   cmd->downsamp);
    printf("        blocklength (%d).  Exiting.\n\n", 
	   blocklen);
    exit(1);
  }

  tlotoa = idata.mjd_i + idata.mjd_f;  /* Topocentric epoch */

  if (cmd->numoutP)
    totnumtowrite = cmd->numout;
  else
    totnumtowrite = (int) idata.N/cmd->downsamp;

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
	offsets[ii][jj] = NEAREST_INT((subdispdt[jj]-dtmp)/dsdt);
      free(subdispdt);
    }

    /* Allocate our data array and start getting data */
    
    printf("De-dispersing using:\n");
    printf("       Subbands = %d\n", cmd->numsub);
    printf("     Average DM = %.7g\n", avgdm);
    if (cmd->downsamp > 1){
      printf("     Downsample = %d\n", cmd->downsamp);
      printf("  New sample dt = %.10g\n", dsdt);
    }
    printf("\n");
    
    outdata = gen_fmatrix(cmd->numsub, worklen/cmd->downsamp);
    numread = get_data(infiles, numinfiles, outdata, 
		       numchan, blocklen, blocksperread, 
		       &obsmask, dispdt, offsets, &padding);

    while (numread==worklen){

      numread /= cmd->downsamp;
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

    printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
    printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
    printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
    printf("De-dispersing and barycentering using:\n");
    printf("       Subbands = %d\n", cmd->numsub);
    printf("     Average DM = %.7g\n", avgdm);
    if (cmd->downsamp > 1){
      printf("     Downsample = %d\n", cmd->downsamp);
      printf("  New sample dt = %.10g\n", dsdt);
    }
    printf("\n");

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
	offsets[ii][jj] = NEAREST_INT((subdispdt[jj]-dtmp)/dsdt);
      free(subdispdt);
    }

    /* Convert the bary TOAs to differences from the topo TOAs in */
    /* units of bin length (dt) rounded to the nearest integer.   */

    dtmp = (btoa[0] - ttoa[0]);
    for (ii=0; ii<numbarypts; ii++)
      btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

    { /* Find the points where we need to add or remove bins */

      int oldbin=0, currentbin;
      double lobin, hibin, calcpt;

      numdiffbins = abs(NEAREST_INT(btoa[numbarypts-1])) + 1;
      diffbins = gen_ivect(numdiffbins);
      diffbinptr = diffbins;
      for (ii=1; ii<numbarypts; ii++){
	currentbin = NEAREST_INT(btoa[ii]);
	if (currentbin != oldbin){
	  if (currentbin > 0){
	    calcpt = oldbin + 0.5;
	    lobin = (ii-1) * TDT / dsdt;
	    hibin = ii * TDT / dsdt;
	  } else {
	    calcpt = oldbin - 0.5;
	    lobin = -((ii-1) * TDT / dsdt);
	    hibin = -(ii * TDT / dsdt);
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

    outdata = gen_fmatrix(cmd->numsub, worklen/cmd->downsamp);
    numread = get_data(infiles, numinfiles, outdata, 
		       numchan, blocklen, blocksperread, 
		       &obsmask, dispdt, offsets, &padding);
    
    while (numread==worklen){ /* Loop to read and write the data */
      int numwritten=0;

      numread /= cmd->downsamp;
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

  /* Calculate the amount of padding we need  */

  if (cmd->numoutP && (cmd->numout > totwrote))
    padwrote = padtowrite = cmd->numout - totwrote;

  /* Write the new info file for the output data */

  idata.dt = dsdt;
  update_infodata(&idata, totwrote, padtowrite, diffbins, 
		  numdiffbins, cmd->downsamp);
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
      for (ii=0; ii<numtowrite/veclen; ii++){
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
  static int dsworklen;
  static float *tempzz, *data1, *data2, *dsdata1=NULL, *dsdata2=NULL; 
  static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
  int totnumread=0, numread=0, ii, jj, tmppad=0, nummasked=0;
  
  if (firsttime){
    if (cmd->maskfileP)
      maskchans = gen_ivect(numchan);
    worklen = blocklen * blocksperread;
    dsworklen = worklen / cmd->downsamp;
    blocksize = blocklen * cmd->numsub;
    data1 = gen_fvect(cmd->numsub * worklen);
    data2 = gen_fvect(cmd->numsub * worklen);
    currentdata = data1;
    lastdata = data2;
    if (cmd->downsamp > 1){
      dsdata1 = gen_fvect(cmd->numsub * dsworklen);
      dsdata2 = gen_fvect(cmd->numsub * dsworklen);
      currentdsdata = dsdata1;
      lastdsdata = dsdata2;
    } else {
      currentdsdata = data1;
      lastdsdata = data2;
    }
  }
  while (firsttime >= 0){
    if (cmd->pkmbP || cmd->bcpmP || cmd->wappP){
      for (ii=0; ii<blocksperread; ii++){
	if (cmd->pkmbP)
	  numread = read_PKMB_subbands(infiles, numfiles, 
				       currentdata+ii*blocksize, 
				       dispdts, cmd->numsub, 0, &tmppad, 
				       maskchans, &nummasked, obsmask);
	else if (cmd->bcpmP)
	  numread = read_BPP_subbands(infiles, numfiles, 
				      currentdata+ii*blocksize, 
				      dispdts, cmd->numsub, 0, &tmppad, 
				      maskchans, &nummasked, obsmask, bppifs);
	else if (cmd->wappP)
	  numread = read_WAPP_subbands(infiles, numfiles, 
				       currentdata+ii*blocksize, 
				       dispdts, cmd->numsub, 0, &tmppad, 
				       maskchans, &nummasked, obsmask);
	if (firsttime==0) totnumread += numread;
	if (numread!=blocklen){
	  for (jj=ii*blocksize; jj<(ii+1)*blocksize; jj++)
	    currentdata[jj] = 0.0;
	}
	if (tmppad) 
	  *padding = 1;
      }
    }
    /* Downsample the subband data if needed */
    if (cmd->downsamp > 1){
      int kk, offset, dsoffset, index, dsindex;
      for (ii=0; ii<dsworklen; ii++){
	dsoffset = ii * cmd->numsub;
	offset = dsoffset * cmd->downsamp;
	for (jj=0; jj<cmd->numsub; jj++){
	  dsindex = dsoffset + jj;
	  index = offset + jj;
	  currentdsdata[dsindex] = 0.0;
	  for (kk=0; kk<cmd->downsamp; kk++){
	    currentdsdata[dsindex] += currentdata[index];
	    index += cmd->numsub;
	  }
	}
      }
    }
    if (firsttime){
      SWAP(currentdata, lastdata);
      SWAP(currentdsdata, lastdsdata);
    }
    firsttime--;
  }
  firsttime = 0;
  for (ii=0; ii<cmd->numdms; ii++)
    float_dedisp(currentdsdata, lastdsdata, dsworklen, 
		 cmd->numsub, offsets[ii], 0.0, outdata[ii]);
  SWAP(currentdata, lastdata);
  SWAP(currentdsdata, lastdsdata);
  if (totnumread != worklen){
    if (cmd->maskfileP)
      free(maskchans);
    free(data1);
    free(data2);
    if (cmd->downsamp > 1){
      free(dsdata1);
      free(dsdata2);
    }
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
			    int *barybins, int numbarybins, int downsamp)
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
  } else {
    for (ii=0; ii<idata->numonoff; ii++){
      idata->onoff[ii*2] /= downsamp;
      idata->onoff[ii*2+1] /= downsamp;
    }
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


