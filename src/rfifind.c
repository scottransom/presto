#include <limits.h>
#include "presto.h"
#include "rfifind_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "rfifind.h"

/* Some function definitions */

void rfifind_plot(int numchan, int numint, int ptsperint, 
		  float timesigma, float freqsigma, 
		  float **dataavg, float **datastd, float **datapow,
		  int *userchan, int numuserchan, 
		  int *userints, int numuserints, 
		  infodata *idata, unsigned char **bytemask, 
		  mask *oldmask, mask *newmask, 
		  rfi *rfivect, int numrfi, 
		  int rfixwin, int rfips, int xwin);
static void write_rfifile(char *rfifilenm, rfi *rfivect, int numrfi,
			  int numchan, int numint, int ptsperint, 
			  int lobin, int numbetween, int harmsum,
			  float fracterror, float freqsigma);
static void write_statsfile(char *statsfilenm, float *datapow,
			    float *dataavg, float *datastd,
			    int numchan, int numint, int ptsperint, 
			    int lobin, int numbetween);
static void read_rfifile(char *rfifilenm, rfi **rfivect, int *numrfi,
			 int *numchan, int *numint, int *ptsperint, 
			 int *lobin, int *numbetween, int *harmsum,
			 float *fracterror, float *freqsigma);
static void read_statsfile(char *statsfilenm, float ***datapow,
			   float ***dataavg, float ***datastd,
			   int *numchan, int *numint, int *ptsperint, 
			   int *lobin, int *numbetween);
int compare_rfi_sigma(const void *ca, const void *cb);
int compare_rfi_numobs(const void *ca, const void *cb);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infiles[MAXPATCHFILES], *bytemaskfile;
  float **dataavg=NULL, **datastd=NULL, **datapow=NULL;
  float *chandata=NULL, powavg, powstd, powmax;
  float inttime, norm, fracterror=RFI_FRACTERROR, freqsigma, timesigma;
  unsigned char *rawdata=NULL, **bytemask=NULL;
  char *outfilenm, *statsfilenm, *maskfilenm;
  char *bytemaskfilenm, *rfifilenm;
  int numchan=0, numint=0, newper=0, oldper=0, numfiles;
  int blocksperint, ptsperint=0, ptsperblock=0, padding=0;
  int numcands, candnum, numrfi=0, numrfivect=NUM_RFI_VECT;
  int ii, jj, kk, slen, numread=0;
  int harmsum=RFI_NUMHARMSUM, lobin=RFI_LOBIN, numbetween=RFI_NUMBETWEEN;
  double davg, dvar, freq, dt, T;
  long long N;
  BPP_ifs bppifs=SUMIFS;
  presto_interptype interptype;
  rfi *rfivect=NULL;
  mask oldmask, newmask;
  fftcand *cands;
  PKMB_tapehdr hdr;
  infodata idata;
  Cmdline *cmd;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
    Program = argv[0];
    printf("\n");
    usage();
    exit(1);
  }

  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);
  freqsigma = cmd->freqsigma;
  timesigma = cmd->timesigma;
  slen = strlen(cmd->outfile)+20;
  numfiles = cmd->argc;

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("               Pulsar Data RFI Finder\n");
  printf("                 by Scott M. Ransom\n");
  printf("           Last Modification:  4 Jan, 2001\n\n");

  /* The following is the root of all the output files */

  outfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(outfilenm, "%s_rfifind", cmd->outfile);

  /* And here are the output file names */

  maskfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(maskfilenm, "%s.mask", outfilenm);
  bytemaskfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(bytemaskfilenm, "%s.bytemask", outfilenm);
  rfifilenm = (char *)calloc(slen, sizeof(char));
  sprintf(rfifilenm, "%s.rfi", outfilenm);
  statsfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(statsfilenm, "%s.stats", outfilenm);
  sprintf(idata.name, "%s", outfilenm);

  /* Read an input mask if wanted */

  if (cmd->maskfileP)
    read_mask(cmd->maskfile, &oldmask);
  else
    oldmask.numchan = oldmask.numint = 0;
  
  if (!cmd->nocomputeP){

    if (!cmd->pkmbP && !cmd->bcpmP && !cmd->wappP){
      printf("\nYou must specify the data format.  Legal types are:\n");
      printf("   -pkmb : Parkes Multibeam\n");
      printf("   -bcpm : Berkeley-Caltech Pulsar Machine\n");
      printf("   -wapp : Arecibo WAPP\n");
      printf("\n");
      exit(0);
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

    for (ii=0; ii<numfiles; ii++){
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }
    printf("\n");

    if (cmd->pkmbP){

      /* Set-up values if we are using the Parkes multibeam */
    
      printf("PKMB input file information:\n");
      get_PKMB_file_info(infiles, numfiles, &N, &ptsperblock, &numchan, 
			 &dt, &T, 1);

      /* Read the first header file and generate an infofile from it */
      
      chkfread(&hdr, 1, HDRLEN, infiles[0]);
      rewind(infiles[0]);
      PKMB_hdr_to_inf(&hdr, &idata);
      PKMB_update_infodata(numfiles, &idata);
      idata.dm = 0.0;
      writeinf(&idata);

    } else if (cmd->bcpmP){

      /* Set-up for the BCPM machines at Green Bank  */

      printf("BCPM input file information:\n");
      get_BPP_file_info(infiles, numfiles, &N, &ptsperblock, &numchan, 
			&dt, &T, &idata, 1);
      BPP_update_infodata(numfiles, &idata);
      idata.dm = 0.0;
      writeinf(&idata);

      /* Which IFs will we use? */
      
      if (cmd->ifsP){
	if (cmd->ifs==0)
	  bppifs = IF0;
	else if (cmd->ifs==1)
	  bppifs = IF1;
	else
	  bppifs = SUMIFS;
      }

    } else if (cmd->wappP){

      /* Set-up for the WAPP machine at Arecibo */

      printf("WAPP input file information:\n");
      get_WAPP_file_info(infiles, cmd->numwapps, numfiles, cmd->clip,
			 &N, &ptsperblock, &numchan, 
			 &dt, &T, &idata, 1);
      WAPP_update_infodata(numfiles, &idata);
      idata.dm = 0.0;
      writeinf(&idata);
    }
    
    /* The number of data points and blocks to work with at a time */
    
    blocksperint = (int) (cmd->time / 
			  (ptsperblock * idata.dt) + 0.5);
    ptsperint = blocksperint * ptsperblock;
    numint = (long long) idata.N / ptsperint;
    if ((long long) idata.N % ptsperint) numint++;
    inttime = ptsperint * idata.dt;
    
    /* Allocate our workarrays */
    
    if (cmd->pkmbP)
      rawdata = gen_bvect(DATLEN * blocksperint);
    else if (cmd->bcpmP || cmd->wappP)
      /* This allocates extra incase both IFs were stored */
      rawdata = gen_bvect(idata.num_chan * ptsperblock * blocksperint);
    dataavg = gen_fmatrix(numint, numchan);
    datastd = gen_fmatrix(numint, numchan);
    datapow = gen_fmatrix(numint, numchan);
    chandata = gen_fvect(ptsperint);
    bytemask = gen_bmatrix(numint, numchan);
    for (ii=0; ii<numint; ii++)
      for (jj=0; jj<numchan; jj++)
	bytemask[ii][jj] = GOODDATA;
    rfivect = rfi_vector(rfivect, numchan, numint, 0, numrfivect);
    if (numbetween==2)
      interptype = INTERBIN;
    else
      interptype = INTERPOLATE;
    
    /* Main loop */

    printf("Writing mask data  to '%s'.\n", maskfilenm);
    printf("Writing  RFI data  to '%s'.\n", rfifilenm);
    printf("Writing statistics to '%s'.\n\n", statsfilenm);
    printf("Massaging the data ...\n\n");
    printf("Amount Complete = %3d%%", oldper);
    fflush(stdout);
    
    for (ii=0; ii<numint; ii++){  /* Loop over the intervals */
      newper = (int) ((float) ii / numint * 100.0 + 0.5);
      if (newper > oldper) {
	printf("\rAmount Complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }
      
      /* Read a chunk of data */
      
      if (cmd->pkmbP)
	numread = read_PKMB_rawblocks(infiles, numfiles, 
				      rawdata, blocksperint, &padding);
      else if (cmd->bcpmP)
	numread = read_BPP_rawblocks(infiles, numfiles, 
				     rawdata, blocksperint, &padding);
      else if (cmd->wappP)
	numread = read_WAPP_rawblocks(infiles, numfiles, 
				      rawdata, blocksperint, &padding);
      if (padding)
	for (jj=0; jj<numchan; jj++)
	  bytemask[ii][jj] |= PADDING;

      for (jj=0; jj<numchan; jj++){  /* Loop over the channels */

	if (cmd->pkmbP)
	  get_PKMB_channel(jj, chandata, rawdata, blocksperint);
	else if (cmd->bcpmP)
	  get_BPP_channel(jj, chandata, rawdata, blocksperint, bppifs);
	else if (cmd->wappP)
	  get_WAPP_channel(jj, chandata, rawdata, blocksperint);

	/* Calculate the averages and standard deviations */
	/* for each point in time.                        */
      
	if (padding){
	  if (cmd->pkmbP){
	    dataavg[ii][jj] = 0.5;
	    datastd[ii][jj] = 0.5;
	    datapow[ii][jj] = 1.0;
	  } else {
	    dataavg[ii][jj] = 0.0;
	    datastd[ii][jj] = 0.0;
	    datapow[ii][jj] = 1.0;
	  }
	} else {
	  avg_var(chandata, ptsperint, &davg, &dvar);
	  dataavg[ii][jj] = davg;
	  datastd[ii][jj] = sqrt(dvar);
	  realfft(chandata, ptsperint, -1);
	  numcands=0;
	  norm = datastd[ii][jj] * datastd[ii][jj] * ptsperint;
	  if (norm==0.0) 
	    norm = (chandata[0]==0.0) ? 1.0 : chandata[0];
	  cands = search_fft((fcomplex *)chandata, ptsperint / 2, 
			     lobin, ptsperint / 2, harmsum, 
			     numbetween, interptype, norm, freqsigma,
			     &numcands, &powavg, &powstd, &powmax);
	  datapow[ii][jj] = powmax;
	  
	  /* Record the birdies */
	  
	  if (numcands){
	    for (kk=0; kk<numcands; kk++){
	      freq = cands[kk].r/inttime;
	      candnum = find_rfi(rfivect, numrfi, freq, RFI_FRACTERROR);
	      if (candnum >= 0){
		update_rfi(rfivect+candnum, freq, cands[kk].sig, jj, ii);
	      } else {
		update_rfi(rfivect+numrfi, freq, cands[kk].sig, jj, ii);
		numrfi++;
		if (numrfi==numrfivect){
		  numrfivect *= 2;
		  rfivect = rfi_vector(rfivect, numchan, numint, 
				       numrfivect/2, numrfivect);
		}
	      }
	    }
	    free(cands);
	  }
	}
      }
    }

    /* Write the data to the output files */
    
    write_rfifile(rfifilenm, rfivect, numrfi, numchan, numint, 
		  ptsperint, lobin, numbetween, harmsum,
		  fracterror, freqsigma);
    write_statsfile(statsfilenm, datapow[0], dataavg[0], datastd[0],
		    numchan, numint, ptsperint, lobin, 
		    numbetween);

  } else { /* If "-nocompute" */

    /* Read the data from the output files */
    
    printf("Reading  RFI data  from '%s'.\n", rfifilenm);
    printf("Reading statistics from '%s'.\n", statsfilenm);
    readinf(&idata, outfilenm);
    read_rfifile(rfifilenm, &rfivect, &numrfi, &numchan, &numint, 
		 &ptsperint, &lobin, &numbetween, &harmsum,
		 &fracterror, &freqsigma);
    numrfivect = numrfi;
    read_statsfile(statsfilenm, &datapow, &dataavg, &datastd,
		   &numchan, &numint, &ptsperint, &lobin, 
		   &numbetween);
    bytemask = gen_bmatrix(numint, numchan);
    printf("Reading  bytemask  from '%s'.\n\n", bytemaskfilenm);
    bytemaskfile = chkfopen(bytemaskfilenm, "rb");
    chkfread(bytemask[0], numint*numchan, 1, bytemaskfile);
    fclose(bytemaskfile);
    for (ii=0; ii<numint; ii++)
      for (jj=0; jj<numchan; jj++)
	bytemask[ii][jj] &= PADDING; /* Clear all but the PADDING bits */
    inttime = ptsperint * idata.dt;
  }

  /* Make the plots */

  rfifind_plot(numchan, numint, ptsperint, timesigma, freqsigma, 
	       dataavg, datastd, datapow, cmd->zapchan, cmd->zapchanC,
	       cmd->zapints, cmd->zapintsC, &idata, bytemask, 
	       &oldmask, &newmask, rfivect, numrfi, 
	       cmd->rfixwinP, cmd->rfipsP, cmd->xwinP);

  /* Write the new mask and bytemask to the file */

  write_mask(maskfilenm, &newmask);
  bytemaskfile = chkfopen(bytemaskfilenm, "wb");
  chkfwrite(bytemask[0], numint*numchan, 1, bytemaskfile);
  fclose(bytemaskfile);
  
  /* Determine the percent of good and bad data */

  {
    int numpad=0, numbad=0, numgood=0;

    for (ii=0; ii<numint; ii++){
      for (jj=0; jj<numchan; jj++){
	if (bytemask[ii][jj]==GOODDATA){
	  numgood++;
	} else {
	  if (bytemask[ii][jj] & PADDING)
	    numpad++;
	  else 
	    numbad++;
	}
      }
    }
    printf("\nTotal number of intervals in the data:  %d\n\n", 
	   numint*numchan);
    printf("  Number of padded intervals:  %7d  (%6.3f%%)\n", 
	   numpad, (float) numpad / (float)(numint*numchan) * 100.0);
    printf("  Number of  good  intervals:  %7d  (%6.3f%%)\n", 
	   numgood, (float) numgood / (float)(numint*numchan) * 100.0);
    printf("  Number of  bad   intervals:  %7d  (%6.3f%%)\n\n", 
	   numbad, (float) numbad / (float)(numint*numchan) * 100.0);
    qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_sigma);
    printf("  Ten most significant birdies:\n");
    printf("#  Sigma     Period(ms)      Freq(Hz)       Number \n");
    printf("----------------------------------------------------\n");
    for(ii=0; ii<10; ii++){
      double pperr;
      char temp1[40], temp2[40];

      if (rfivect[ii].freq_var==0.0){
	pperr = 0.0;
	sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
	sprintf(temp2, " %-14g", 1000.0/rfivect[ii].freq_avg);
      } else {
	pperr = 1000.0 * sqrt(rfivect[ii].freq_var) / 
	  (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
	nice_output_2(temp1, rfivect[ii].freq_avg, sqrt(rfivect[ii].freq_var), -15);
	nice_output_2(temp2, 1000.0/rfivect[ii].freq_avg, pperr, -15);
      }
      printf("%-2d %-8.2f %13s %13s %-8d\n", ii+1, rfivect[ii].sigma_avg, 
	     temp2, temp1, rfivect[ii].numobs);
    } 
    qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_numobs);
    printf("\n  Ten most numerous birdies:\n");
    printf("#  Number    Period(ms)      Freq(Hz)       Sigma \n");
    printf("----------------------------------------------------\n");
    for(ii=0; ii<10; ii++){
      double pperr;
      char temp1[40], temp2[40];

      if (rfivect[ii].freq_var==0.0){
	pperr = 0.0;
	sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
	sprintf(temp2, " %-14g", 1000.0/rfivect[ii].freq_avg);
      } else {
	pperr = 1000.0 * sqrt(rfivect[ii].freq_var) / 
	  (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
	nice_output_2(temp1, rfivect[ii].freq_avg, sqrt(rfivect[ii].freq_var), -15);
	nice_output_2(temp2, 1000.0/rfivect[ii].freq_avg, pperr, -15);
      }
      printf("%-2d %-8d %13s %13s %-8.2f\n", ii+1, rfivect[ii].numobs, 
	     temp2, temp1, rfivect[ii].sigma_avg);
    } 
    printf("\nDone.\n\n");
  }

  /* Close the files and cleanup */

  free_rfi_vector(rfivect, numrfivect);
  free_mask(newmask);
  if (cmd->maskfileP)
    free_mask(oldmask);
  free(outfilenm);
  free(statsfilenm);
  free(bytemaskfilenm);
  free(maskfilenm);
  free(rfifilenm);
  free(dataavg[0]); free(dataavg);
  free(datastd[0]); free(datastd);
  free(datapow[0]); free(datapow);
  free(bytemask[0]); free(bytemask);
  if (!cmd->nocomputeP){
    for (ii=0; ii<numfiles; ii++)
      fclose(infiles[ii]);
    free(chandata);
    free(rawdata);
  }
  return (0);
}

static void write_rfifile(char *rfifilenm, rfi *rfivect, int numrfi,
			  int numchan, int numint, int ptsperint, 
			  int lobin, int numbetween, int harmsum,
			  float fracterror, float freqsigma)
{
  FILE *outfile;
  int ii;
  
  outfile = chkfopen(rfifilenm, "wb");
  chkfwrite(&numchan, sizeof(int), 1, outfile);
  chkfwrite(&numint, sizeof(int), 1, outfile);
  chkfwrite(&ptsperint, sizeof(int), 1, outfile);
  chkfwrite(&lobin, sizeof(int), 1, outfile);
  chkfwrite(&numbetween, sizeof(int), 1, outfile);
  chkfwrite(&harmsum, sizeof(int), 1, outfile);
  chkfwrite(&numrfi, sizeof(int), 1, outfile);
  chkfwrite(&fracterror, sizeof(float), 1, outfile);
  chkfwrite(&freqsigma, sizeof(float), 1, outfile);
  for (ii=0; ii<numrfi; ii++)
    write_rfi(outfile, rfivect+ii, numchan, numint);
  fclose(outfile);
}

static void write_statsfile(char *statsfilenm, float *datapow,
			    float *dataavg, float *datastd,
			    int numchan, int numint, int ptsperint, 
			    int lobin, int numbetween)
{
  FILE *outfile;
  
  outfile = chkfopen(statsfilenm, "wb");
  chkfwrite(&numchan, sizeof(int), 1, outfile);
  chkfwrite(&numint, sizeof(int), 1, outfile);
  chkfwrite(&ptsperint, sizeof(int), 1, outfile);
  chkfwrite(&lobin, sizeof(int), 1, outfile);
  chkfwrite(&numbetween, sizeof(int), 1, outfile);
  chkfwrite(datapow, sizeof(float), numchan * numint, outfile);
  chkfwrite(dataavg, sizeof(float), numchan * numint, outfile);
  chkfwrite(datastd, sizeof(float), numchan * numint, outfile);
  fclose(outfile);
}

static void read_rfifile(char *rfifilenm, rfi **rfivect, int *numrfi,
			 int *numchan, int *numint, int *ptsperint, 
			 int *lobin, int *numbetween, int *harmsum,
			 float *fracterror, float *freqsigma)
{
  FILE *outfile;
  int ii;
  
  outfile = chkfopen(rfifilenm, "rb");
  chkfread(numchan, sizeof(int), 1, outfile);
  chkfread(numint, sizeof(int), 1, outfile);
  chkfread(ptsperint, sizeof(int), 1, outfile);
  chkfread(lobin, sizeof(int), 1, outfile);
  chkfread(numbetween, sizeof(int), 1, outfile);
  chkfread(harmsum, sizeof(int), 1, outfile);
  chkfread(numrfi, sizeof(int), 1, outfile);
  chkfread(fracterror, sizeof(float), 1, outfile);
  chkfread(freqsigma, sizeof(float), 1, outfile);
  *rfivect = rfi_vector(*rfivect, *numchan, *numint, 0, *numrfi);
  for (ii=0; ii<*numrfi; ii++)
    read_rfi(outfile, *rfivect+ii, *numchan, *numint);
  fclose(outfile);
}

static void read_statsfile(char *statsfilenm, float ***datapow,
			   float ***dataavg, float ***datastd,
			   int *numchan, int *numint, int *ptsperint, 
			   int *lobin, int *numbetween)
{
  FILE *outfile;
  
  outfile = chkfopen(statsfilenm, "rb");
  chkfread(numchan, sizeof(int), 1, outfile);
  chkfread(numint, sizeof(int), 1, outfile);
  chkfread(ptsperint, sizeof(int), 1, outfile);
  chkfread(lobin, sizeof(int), 1, outfile);
  chkfread(numbetween, sizeof(int), 1, outfile);
  *dataavg = gen_fmatrix(*numint, *numchan);
  *datastd = gen_fmatrix(*numint, *numchan);
  *datapow = gen_fmatrix(*numint, *numchan);
  chkfread(*(datapow[0]), sizeof(float), *numchan * *numint, outfile);
  chkfread(*(dataavg[0]), sizeof(float), *numchan * *numint, outfile);
  chkfread(*(datastd[0]), sizeof(float), *numchan * *numint, outfile);
  fclose(outfile);
}
