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

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infile, *outfile, *maskfile, *rfifile;
  float **dataavg, **datastd, **datapow, *outdata, **indata;
  float powavg, powstd, powmax, inttime, norm;
  unsigned char **datamask, *tbuffer;
  char *datafilenm, *rootfilenm, *maskfilenm, *rfifilenm, *cptr;
  int numblocks=0, numchan=0, numint, newper=0, oldper=0;
  int blocksperint, ptsperint, ptsperblock=0, tbuffer_size;
  int bitsperpt, numcands, candnum, numrfi=0, numrfivect=NUM_RFI_VECT;
  int ii, jj, kk, slen, numread=0;
  double davg, dvar, freq;
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

#ifdef DEBUG
  showOptionValues();
#endif

  printf("\n\n");
  printf("               Pulsar Data RFI Finder\n");
  printf("                 by Scott M. Ransom\n");
  printf("           Last Modification:  5 Dec, 2000\n\n");

  /* Determine the root input file name and the input info file name */

  cptr = strrchr(cmd->argv[0], '.');
  if (cptr==NULL){
    printf("\nThe input filename (%s) must have a suffix!\n\n", 
	   cmd->argv[0]);
    exit(1);
  }
  slen = cptr - cmd->argv[0];
  rootfilenm = (char *)calloc(slen+1, sizeof(char));
  strncpy(rootfilenm, cmd->argv[0], slen);
  maskfilenm = (char *)calloc(slen+6, sizeof(char));
  sprintf(maskfilenm, "%s.mask", rootfilenm);
  maskfile = chkfopen(maskfilenm, "wb");
  rfifilenm = (char *)calloc(slen+5, sizeof(char));
  sprintf(rfifilenm, "%s.rfi", rootfilenm);
  rfifile = chkfopen(rfifilenm, "w");
  if (!cmd->pkmbP && !cmd->ebppP && !cmd->gbppP){
    printf("\nYou must specify the data format.  Legal types are:\n");
    printf("   -pkmb : Parkes Multibeam\n");
    printf("   -ebpp : Effelsberg-Berkeley Pulsar Processor\n");
    printf("   -gbpp : Green Bank-Berkeley Pulsar Processor\n");
    printf("\n");
    exit(0);
  } else if (cmd->pkmbP){
    printf("Reading Parkes Multibeam data from '%s'\n\n", 
	   cmd->argv[0]);
  } else if (cmd->ebppP){
    printf("Reading Effelsberg PSR data from '%s'\n\n", 
	   cmd->argv[0]);
  } else if (cmd->gbppP){
    printf("Reading GBPP PSR data from '%s'\n\n", 
	   cmd->argv[0]);
  }

  /* Open the raw data file */

  infile = chkfopen(cmd->argv[0], "rb");

  /* Determine the other file names and open the output data file */

  slen = strlen(cmd->outfile)+20;
  datafilenm = (char *)calloc(slen, sizeof(char));
  sprintf(datafilenm, "%s_DM0.0_masked.dat", cmd->outfile);
  printf("Writing masked DM=0.0 data to '%s'.\n", datafilenm);
  sprintf(idata.name, "%s_DM0.0_masked", cmd->outfile);
  printf("Writing mask data to '%s'.\n", maskfilenm);
  maskfile = chkfopen(maskfilenm, "wb");

  /* Set-up values if we are using the Parkes multibeam */

  if (cmd->pkmbP) {

    /* Read the first header file and generate an infofile from it */

    chkfread(&hdr, 1, HDRLEN, infile);
    rewind(infile);
    multibeam_hdr_to_inf(&hdr, &idata);
    idata.dm = 0.0;
    numchan = idata.num_chan;
    ptsperblock = DATLEN * 8 / numchan;
    numblocks = chkfilelen(infile, RECLEN);
    idata.N = (double) (numblocks * ptsperblock);
    writeinf(&idata);
    bitsperpt = 1;

    /* The data collection routine to use */

    readrec_ptr = read_rawmultibeam;
  }

  /* Set-up values if we are using the Effelsberg-Berkeley Pulsar Processor */
  /*   NOTE:  This code is not yet implemented.                             */

  if (cmd->ebppP) {
    bitsperpt = 4;
  }

  /* The number of data points and blocks to work with at a time */
  
  blocksperint = (int) (cmd->time * 60.0 / 
			(ptsperblock * idata.dt) + 0.5);
  numint = numblocks / blocksperint;
  if (numblocks % blocksperint) numint++;
  ptsperint = blocksperint * ptsperblock;
  inttime = ptsperint * idata.dt;

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

    numread = (*readrec_ptr)(infile, indata[0], numchan, blocksperint);
      
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
			 RFI_LOBIN, RFI_NUMHARMSUM, RFI_NUMBETWEEN,
			 INTERBIN, norm, cmd->sigma,
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

  /* Generate the mask */

  /* Make the plots */

  rfifind_plot(numchan, numint, ptsperint, cmd->sigma, 
	       dataavg, datastd, datapow, datamask, &idata, 1);

  /* Close the files and cleanup */

  free_rfi_vector(rfivect, numrfivect);
  fclose(infile);
  /* fclose(outfile); */
  fclose(maskfile);
  fclose(rfifile);
  free(tbuffer);
  free(rootfilenm);
  free(maskfilenm);
  free(datafilenm);
  free(rfifilenm);
  free(dataavg[0]); free(dataavg);
  free(datastd[0]); free(datastd);
  free(datapow[0]); free(datapow);
  free(indata[0]); free(indata);
  free(datamask[0]); free(datamask);
  free(outdata);
  return (0);
}
