#include <limits.h>
#include "presto.h"
#include "rfifind_cmd.h"
#include "multibeam.h"
#include "rfifind.h"

/* Some function definitions */

int (*readrec_ptr)(FILE *file, float *data, int numchan, int numblocks);
void rfifind_plot(int numchan, int numint, int ptsperint, float sigma, 
		  float **dataavg, float **datastd, float **datapow, 
		  unsigned char **datamask, infodata *idata, int xwin);
static void write_rfifile(char *rfifilenm, rfi *rfivect, int numrfi,
			  int numchan, int numint, int ptsperint, 
			  int lobin, int numbetween, int harmsum,
			  float fracterror, float sigma);
static void write_statsfile(char *statsfilenm, float *datapow,
			    float *dataavg, float *datastd,
			    int numchan, int numint, int ptsperint, 
			    int lobin, int numbetween, float sigma);
static void read_rfifile(char *rfifilenm, rfi **rfivect, int *numrfi,
			 int *numchan, int *numint, int *ptsperint, 
			 int *lobin, int *numbetween, int *harmsum,
			 float *fracterror, float *sigma);
static void read_statsfile(char *statsfilenm, float ***datapow,
			   float ***dataavg, float ***datastd,
			   int *numchan, int *numint, int *ptsperint, 
			   int *lobin, int *numbetween, float *sigma);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infiles[MAXPATCHFILES];
  float **dataavg=NULL, **datastd=NULL, **datapow=NULL;
  float *outdata=NULL, **indata=NULL, powavg, powstd, powmax;
  float inttime, norm, fracterror=RFI_FRACTERROR, sigma;
  unsigned char **datamask=NULL, *tbuffer=NULL;
  char *outfilenm, *statsfilenm, *maskfilenm, *rfifilenm;
  int numchan=0, numint=0, newper=0, oldper=0, numfiles;
  int blocksperint, ptsperint=0, ptsperblock=0, tbuffer_size;
  int numcands, candnum, numrfi=0, numrfivect=NUM_RFI_VECT;
  int ii, jj, kk, slen, numread=0, compute=1;
  int harmsum=RFI_NUMHARMSUM, lobin=RFI_LOBIN, numbetween=RFI_NUMBETWEEN;
  double davg, dvar, freq, dt;
  long long numpoints[MAXPATCHFILES], padpoints[MAXPATCHFILES];
  presto_interptype interptype;
  rfi *rfivect=NULL;
  fftcand *cands;
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
  sigma = cmd->sigma;
  slen = strlen(cmd->outfile)+20;
  numfiles = cmd->argc;
   
#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("               Pulsar Data RFI Finder\n");
  printf("                 by Scott M. Ransom\n");
  printf("           Last Modification:  5 Dec, 2000\n\n");

  /* The following is the root of all the output files */

  outfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(outfilenm, "%s_rfifind", cmd->outfile);

  /* And here are the output file names */

  maskfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(maskfilenm, "%s.mask", outfilenm);
  rfifilenm = (char *)calloc(slen, sizeof(char));
  sprintf(rfifilenm, "%s.rfi", outfilenm);
  statsfilenm = (char *)calloc(slen, sizeof(char));
  sprintf(statsfilenm, "%s.stats", outfilenm);
  sprintf(idata.name, "%s", outfilenm);

  if (compute){

    if (!cmd->pkmbP && !cmd->ebppP && !cmd->gbppP){
      printf("\nYou must specify the data format.  Legal types are:\n");
      printf("   -pkmb : Parkes Multibeam\n");
      printf("   -ebpp : Effelsberg-Berkeley Pulsar Processor\n");
      printf("   -gbpp : Green Bank-Berkeley Pulsar Processor\n");
      printf("\n");
      exit(0);
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
    } else if (cmd->gbppP){
      if (numfiles > 1)
	printf("Reading Green Bank GBPP data from %d files:\n", numfiles);
      else
	printf("Reading Green Bank GBPP data from 1 file:\n");
    }

    /* Open the raw data file */

    for (ii=0; ii<numfiles; ii++){
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }
    printf("\nWriting mask  data to '%s'.\n", maskfilenm);
    printf("Writing RFI   data to '%s'.\n", rfifilenm);
    printf("Writing stats data to '%s'.\n\n", statsfilenm);
    
    /* Set-up values if we are using the Parkes multibeam */
    
    if (cmd->pkmbP) {
      
      get_pkmb_file_info(infiles, numfiles, numpoints, padpoints, 
			 &dt, &numchan, 1);

      /* Read the first header file and generate an infofile from it */
      
      chkfread(&hdr, 1, HDRLEN, infiles[0]);
      rewind(infiles[0]);
      multibeam_hdr_to_inf(&hdr, &idata);
      idata.dm = 0.0;
      idata.N = 0.0;
      for (ii=0; ii<numfiles; ii++)
	idata.N += numpoints[ii] + padpoints[ii];
      writeinf(&idata);
      ptsperblock = DATLEN * 8 / numchan;
      
      /* The data collection routine to use */
      
      readrec_ptr = read_rawmultibeam;
    }
    
    /* Set-up values if we are using the EBPP     */
    
    if (cmd->ebppP) {
      /* This code is not yet implemented. */
    }
    
    /* The number of data points and blocks to work with at a time */
    
    blocksperint = (int) (cmd->time * 60.0 / 
			  (ptsperblock * idata.dt) + 0.5);
    ptsperint = blocksperint * ptsperblock;
    numint = (long long) idata.N / ptsperint;
    if ((long long) idata.N % ptsperint) numint++;
    inttime = ptsperint * idata.dt;
    
exit(0);

    /* Allocate our workarrays */
    
    dataavg = gen_fmatrix(numint, numchan);
    datastd = gen_fmatrix(numint, numchan);
    datapow = gen_fmatrix(numint, numchan);
    indata = gen_fmatrix(ptsperint, numchan);
    datamask = gen_bmatrix(numint, numchan);
    for (ii=0; ii<numint; ii++)
      for (jj=0; jj<numchan; jj++)
	datamask[ii][jj] = GOOD;
    outdata = gen_fvect(ptsperint);
    tbuffer_size = (numint + numchan) / 2;
    tbuffer = gen_bvect(tbuffer_size);
    rfivect = rfi_vector(rfivect, numchan, numint, 0, numrfivect);
    if (numbetween==2)
      interptype = INTERBIN;
    else
      interptype = INTERPOLATE;
    
    /* Main loop */

    printf("Massaging the data ...\n\n");
    printf("Amount Complete = %3d%%", oldper);
    fflush(stdout);
    
    for (ii=0; ii<numint; ii++){
      newper = (int) ((float) ii / numint * 100.0 + 0.5);
      if (newper > oldper) {
	printf("\rAmount Complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }
      
      /* Read a chunk of data */
      
      numread = (*readrec_ptr)(infiles[0], indata[0], numchan, blocksperint);
      
      /* Transpose the data so we can work on the channels */
      
      transpose_float(indata[0], numint, numchan, tbuffer, tbuffer_size);
      
      /* Calculate the averages and standard deviations */
      /* for each point in time.                        */
      
      for (jj=0; jj<numchan; jj++){
	avg_var(indata[jj], ptsperint, &davg, &dvar);
	dataavg[ii][jj] = davg;
	datastd[ii][jj] = sqrt(dvar);
      }
      
      /* Calculate the FFT of each time interval and search it */
      
      for (jj=0; jj<numchan; jj++){
	for (kk=0; kk<ptsperint; kk++)
	  outdata[kk] = *(indata[jj]+kk);
	realfft(outdata, ptsperint, -1);
	numcands=0;
	norm = datastd[ii][jj] * datastd[ii][jj] * ptsperint;
	cands = search_fft((fcomplex *)outdata, ptsperint / 2, 
			   lobin, harmsum, numbetween,
			   interptype, norm, sigma,
			   &numcands, &powavg, &powstd, &powmax);
	printf("powavg = %f  powstd = %f  powmax = %f\n", 
	       powavg, powstd, powmax);
	datapow[ii][jj] = powmax;
	
	/* Record the birdies */
	
	if (numcands){
	  for (kk=0; kk<numcands; kk++){
	    freq = cands[kk].r/inttime;
	    candnum = find_rfi(rfivect, numrfi, freq, RFI_FRACTERROR);
	    if (candnum >= 0){
	      printf("  Another %.4f Hz birdie (channel = %d sigma = %.2f)\n", 
		     freq, jj, cands[kk].sig);
	      update_rfi(rfivect+candnum, freq, cands[kk].sig, jj, ii);
	    } else {
	      printf("\nNew %.4f Hz birdie (channel = %d, sigma = %.2f)\n", 
		     freq, jj, cands[kk].sig);
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

    /* Write the data to the output files */
    
    write_rfifile(rfifilenm, rfivect, numrfi, numchan, numint, 
		  ptsperint, lobin, numbetween, harmsum,
		  fracterror, sigma);
    write_statsfile(statsfilenm, datapow[0], dataavg[0], datastd[0],
		    numchan, numint, ptsperint, lobin, 
		    numbetween, sigma);

  } else {

    /* Read the data from the output files */
    
    printf("Reading mask  data from '%s'.\n", maskfilenm);
    printf("Reading RFI   data from '%s'.\n", rfifilenm);
    printf("Reading stats data from '%s'.\n\n", statsfilenm);
    readinf(&idata, outfilenm);
    read_rfifile(rfifilenm, &rfivect, &numrfi, &numchan, &numint, 
		 &ptsperint, &lobin, &numbetween, &harmsum,
		 &fracterror, &sigma);
    read_statsfile(statsfilenm, &datapow, &dataavg, &datastd,
		   &numchan, &numint, &ptsperint, &lobin, 
		   &numbetween, &sigma);
    inttime = ptsperint * idata.dt;
  }

  /* Generate the mask */

  /* Make the plots */

  rfifind_plot(numchan, numint, ptsperint, cmd->sigma, 
	       dataavg, datastd, datapow, datamask, &idata, 0);

  /* Close the files and cleanup */

  free_rfi_vector(rfivect, numrfivect);
  free(outfilenm);
  free(statsfilenm);
  free(maskfilenm);
  free(rfifilenm);
  if (compute){
    for (ii=0; ii<numfiles; ii++)
      fclose(infiles[ii]);
    free(tbuffer);
    free(dataavg[0]); free(dataavg);
    free(datastd[0]); free(datastd);
    free(datapow[0]); free(datapow);
    free(indata[0]); free(indata);
    free(datamask[0]); free(datamask);
    free(outdata);
  }
  return (0);
}

static void write_rfifile(char *rfifilenm, rfi *rfivect, int numrfi,
			  int numchan, int numint, int ptsperint, 
			  int lobin, int numbetween, int harmsum,
			  float fracterror, float sigma)
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
  chkfwrite(&sigma, sizeof(float), 1, outfile);
  for (ii=0; ii<numrfi; ii++)
    write_rfi(outfile, rfivect+ii, numchan, numint);
  fclose(outfile);
}

static void write_statsfile(char *statsfilenm, float *datapow,
			    float *dataavg, float *datastd,
			    int numchan, int numint, int ptsperint, 
			    int lobin, int numbetween, float sigma)
{
  FILE *outfile;
  
  outfile = chkfopen(statsfilenm, "wb");
  chkfwrite(&numchan, sizeof(int), 1, outfile);
  chkfwrite(&numint, sizeof(int), 1, outfile);
  chkfwrite(&ptsperint, sizeof(int), 1, outfile);
  chkfwrite(&lobin, sizeof(int), 1, outfile);
  chkfwrite(&numbetween, sizeof(int), 1, outfile);
  chkfwrite(&sigma, sizeof(float), 1, outfile);
  chkfwrite(datapow, sizeof(float), numchan * numint, outfile);
  chkfwrite(dataavg, sizeof(float), numchan * numint, outfile);
  chkfwrite(datastd, sizeof(float), numchan * numint, outfile);
  fclose(outfile);
}

static void read_rfifile(char *rfifilenm, rfi **rfivect, int *numrfi,
			 int *numchan, int *numint, int *ptsperint, 
			 int *lobin, int *numbetween, int *harmsum,
			 float *fracterror, float *sigma)
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
  chkfread(sigma, sizeof(float), 1, outfile);
  *rfivect = rfi_vector(*rfivect, *numchan, *numint, 0, *numrfi);
  for (ii=0; ii<*numrfi; ii++)
    read_rfi(outfile, *rfivect+ii, *numchan, *numint);
  fclose(outfile);
}

static void read_statsfile(char *statsfilenm, float ***datapow,
			   float ***dataavg, float ***datastd,
			   int *numchan, int *numint, int *ptsperint, 
			   int *lobin, int *numbetween, float *sigma)
{
  FILE *outfile;
  
  outfile = chkfopen(statsfilenm, "rb");
  chkfread(numchan, sizeof(int), 1, outfile);
  chkfread(numint, sizeof(int), 1, outfile);
  chkfread(ptsperint, sizeof(int), 1, outfile);
  chkfread(lobin, sizeof(int), 1, outfile);
  chkfread(numbetween, sizeof(int), 1, outfile);
  chkfread(sigma, sizeof(float), 1, outfile);
  *dataavg = gen_fmatrix(*numint, *numchan);
  *datastd = gen_fmatrix(*numint, *numchan);
  *datapow = gen_fmatrix(*numint, *numchan);
  chkfread(*(datapow[0]), sizeof(float), *numchan * *numint, outfile);
  chkfread(*(dataavg[0]), sizeof(float), *numchan * *numint, outfile);
  chkfread(*(datastd[0]), sizeof(float), *numchan * *numint, outfile);
  fclose(outfile);
}
