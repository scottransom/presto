#include <limits.h>
#include "presto.h"
#include "rfifind_cmd.h"
#include "multibeam.h"
#include "rfifind.h"

/* Some function definitions */

int (*readrec_ptr)(FILE *file, float *data, int numchan, int numblocks);
static void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infile, *outfile, *maskfile, *rfifile;
  float **dataavg, **datavar, **datapow, *outdata;
  float *avg_chan, *avg_int, *var_chan, *var_int, *pow_chan, *pow_int;
  float powavg, powvar, powmax, inttime;
  unsigned char **datamask, *tbuffer;
  char *datafilenm, *rootfilenm, *maskfilenm, *rfifilenm;
  int numblocks, numchan, numint, newper=0, oldper=0, index;
  int blocksperint, ptsperint, ptsperblock, tbuffer_size;
  int bitsperpt, numcands, candnum, numrfi=0, numrfiobs=NUM_RFI_OBS;
  long ii, jj, kk;
  long numread=0, numtowrite=0, totwrote=0, datawrote=0;
  rfi_obs **rfiobs=NULL;
  rfi_instance newrfi;
  rawbincand cands[10];
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
  rfifilenm = (char *)calloc(slen+5, sizeof(char));
  sprintf(rfifilenm, "%s.rfi", rootfilenm);
  if (!cmd->pkmbP && !cmd->ebppP && !cmd->gbpp){
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
  outfile = chkfopen(datafilenm, "wb");
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
  datavar = gen_fmatrix(numint, numchan);
  datapow = gen_fmatrix(numint, numchan);
  indata = gen_fmatrix(ptsperint, numchan);
  datamask = gen_bvect(numint * numchan);
  for (ii=0; ii<numint * numchan; ii++)
    datamask[ii] = GOOD;
  avg_chan = gen_fvect(numchan);
  avg_int = gen_fvect(numint);
  var_chan = gen_fvect(numchan);
  var_int = gen_fvect(numint);
  pow_chan = gen_fvect(numchan);
  pow_int = gen_fvect(numint);
  outdata = gen_fvect(ptsperint);
  tbuffer_size = (numint + numchan) / 2;
  tbuffer = gen_bvect(tbuffer_size);
  create_rfi_obs_vector(rfiobs, 0, numrfiobs);

  /* Main loop */

  printf("Massaging the data ...\n\n");
  printf("Amount Complete = %3d%%", oldper);

  for (ii=0; ii<numint; ii++){
    newper = (int) ((float) ii / numint * 100.0 + 0.5);
    if (newper > oldper) {
      printf("\rAmount Complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
    }

    /* Read a chunk of data */

    numread = (*readrec_ptr)(infile, indata, numchan, blocksperint);
      
    /* Calculate the averages and standard deviations */
    /* for each point in time.                        */

    for (jj=0; jj<ptsperint; jj++)
      avg_var(indata[jj], numchan, &(dataavg[ii][jj]), 
	      &(datavar[ii][jj]));

    /* Transpose the data so we can work on the channels */

    transpose_float(indata[0], numint, numchan, tbuffer, tbuffer_size);

    /* Calculate the FFT of each time interval and search it */

    for (jj=0; jj<numchan; jj++){
      index = jj * ptsperint;
      for (kk=0; kk<ptsperint; kk++)
	outdata[kk] = data[index+kk];
      realfft(outdata, ptsperint, -1);
      numcands=0;
      cands = search_fft((fcomplex *)outdata, ptsperint>>1, 
			 RFI_LOBIN, RFI_NUMHARMSUM, RFI_NUMBETWEEN,
			 INTERBIN, outdata[0], cmd->sigma,
			 &numcands, &powavg, &powvar, &powmax);
      printf("\npowavg = %.3f powvar = %.3f powmax = %.3f\n", 
	     powavg, powvar, powmax); 
      datapow[ii][jj] = powmax;
      if (numcands){
	for (kk=0; kk<numcands; kk++){
	  /* Store the birdie */
	  candnum = find_rfind(rfiobs, cands[kk], RFI_FRACTERROR);
	  newrfi.freq = cands[kk].r/inttime;
	  newrfi.power = cands[kk].p;
	  newrfi.fftbin = cands[kk].r;
	  newrfi.fftbins = ptsperint;
	  newrfi.inttime = inttime;
	  newrfi.channel = jj;
	  newrfi.intnum = ii;
	  if (candnum >= 0){
	    printf("  Found another %.4f Hz birdie (chan = %d pow = %.2f)\n", 
		   newrfi.freq, newrfi.channel, newrfi.power);
	    add_rfi_instance(rfiobs[candnum], newrfi);
	  } else {
	    printf("New %.4f Hz birdie (chan = %d pow = %.2f)\n", 
		   newrfi.freq, newrfi.channel, newrfi.power);
	    rfiobs[numrfi] = create_rfi_obs(newrfi);
	    numrfi++;
	    if (numrfi==numrifobs-1){
	      create_rfi_obs_vector(rfiobs, numrifobs, numrfiobs * 2);
	      numrfiobs *= 2;
	    }
	  }
	}
	free(cands);
      }
    }
  }

  /* Close the files and cleanup */

  fclose(infile);
  fclose(outfile);
  fclose(maskfile);
  fclose(rfifile);
  free(tbuffer);
  free(rootfilenm);
  free(maskfilenm);
  free(datafilenm);
  free(rfifilenm);
  free(dataavg[0]); free(dataavg);
  free(datavar[0]); free(datavar);
  free(datapow[0]); free(datapow);
  free(indata[0]); free(indata);
  free(datamask);
  free(avg_chan);
  free(avg_int);
  free(var_chan);
  free(var_int);
  free(pow_chan);
  free(pow_int);
  free(outdata);
  return (0);
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
