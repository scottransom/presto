#include "presto.h"
#include "mask.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "fitsfile.h"
#include "fitshead.h"

void spigot2sigprocfb(SPIGOT_INFO *spigot, sigprocfb *fb, char *filenmbase)
{
  int h_or_d, m;
  double s;

  strncpy(fb->inpfile, filenmbase, 40);
  strncpy(fb->source_name, spigot->object, 80);
  fb->nifs = 1;
  if (spigot->num_samplers==1)
    strncpy(fb->ifstream, "YXXX", 8);
  else if (spigot->num_samplers==2)
    strncpy(fb->ifstream, "YYXX", 8);
  fb->tstart = spigot->MJD_obs;
  fb->tsamp = spigot->dt_us/1e6;
  hours2hms(spigot->ra/15.0, &h_or_d, &m, &s);
  fb->src_raj = h_or_d*10000.0+m*100.0+s;
  deg2dms(spigot->dec, &h_or_d, &m, &s);
  fb->src_dej = h_or_d*10000.0+m*100.0+s;
  fb->az_start = 0.0;
  fb->za_start = 0.0;
  fb->nchans = spigot->lags_per_sample;
  fb->foff = spigot->bandwidth/fb->nchans;
  fb->fch1 = spigot->freq_ctr + (fb->nchans/2-0.5)*fb->foff;
  fb->foff = -fb->foff;
  fb->machine_id = 7;
  fb->telescope_id = 6;
  fb->nbits = 8;
  fb->sumifs = spigot->summed_pols;
  if (fb->sumifs) fb->nifs = 1;
  else {
    if (spigot->num_samplers==2) fb->nifs = 2;
    else fb->nifs = 1;
  }
      
  /* The following are not necessary for writing filterbank files */
  fb->headerlen = 0;
  fb->N = 0;
}


int main(int argc, char *argv[])
{
  FILE **infiles, *outfile, *scalingfile;
  int filenum, ii=0, ptsperblock, numchan, numfiles, bytes_per_read, scaling=0, numscalings;
  long long N;
  char outfilenm[200], filenmbase[200], scalingnm[200], rawlags[4096];
  unsigned char output_samples[2048];
  float *scalings;
  double dt, T;
  SPIGOT_INFO *spigots;
  sigprocfb fb;
  infodata idata;

  if (argc==1){
    fprintf(stderr, "Usage: spigot2filterbank SPIGOT_files\n");
    exit(0);
  }
  printf("\nConverting raw SPIGOT FITs data into SIGPROC filterbank format.\n\n");
  strncpy(filenmbase, argv[1], strlen(argv[1])-10);
  filenmbase[strlen(argv[1])-10] = '\0';
  sprintf(outfilenm, "%s.fil", filenmbase);
  printf("Writing data to file '%s'.\n\n", outfilenm);

  /* Attempt to read a file with lag scalings in it */
  if ((scalingfile=fopen(scalingnm, "fb"))){
    /* Determine the length of the file */
    numscalings = (int)chkfilelen(scalingfile, sizeof(float));
    /* Create the array and read 'em */
    scalings = gen_fvect(numscalings);
    chkfread(scalings, sizeof(float), numscalings, scalingfile);
    scaling = 1;
    /* close the scaling file */
    fclose(scalingfile);
    printf("Scaling the lags with the %d values found in '%s'\n\n", 
	   numscalings, scalingnm);
  }


  /* Read and convert the basic SPIGOT file information */
  
  numfiles = argc-1;
  printf("Spigot card input file information:\n");
  spigots = (SPIGOT_INFO *)malloc(sizeof(SPIGOT_INFO)*numfiles);
  infiles = (FILE **)malloc(sizeof(FILE *)*numfiles);
  for (filenum=1; filenum<=numfiles; filenum++){
    printf("  '%s'\n", argv[filenum]);
    infiles[filenum-1] = chkfopen(argv[filenum], "rb");
    read_SPIGOT_header(argv[filenum], spigots+filenum-1);
    rewind(infiles[filenum-1]);
  }
  printf("\n");
  /* The following is necessary in order to initialize all the */
  /* static variables in spigot.c                              */
  get_SPIGOT_file_info(infiles, spigots, numfiles, 0, 0, &N, 
		       &ptsperblock, &numchan, &dt, &T, &idata, 1);
  spigot2sigprocfb(&(spigots[0]), &fb, filenmbase);
  fb.N = N;
  free(spigots);

  /* Write the header */
  outfile = chkfopen(outfilenm, "wb");
  write_filterbank_header(&fb, outfile);

  /* Step throught he SPIGOT files */
  ii = 0;
  for (filenum=1; filenum<argc; filenum++){
    printf("Reading from file '%s'...\n", argv[filenum]);
    chkfseek(infiles[filenum-1], spigots[filenum-1].header_len, SEEK_SET);
    bytes_per_read = spigots[filenum-1].lags_per_sample*\
      spigots[filenum-1].bits_per_lag/8;
    
    /* Loop over the samples in the file */
    while (chkfread(rawlags, bytes_per_read, 1, infiles[filenum-1])){
      if (scaling)
	convert_SPIGOT_point(rawlags, output_samples, SUMIFS, scalings[ii]);
      else
	convert_SPIGOT_point(rawlags, output_samples, SUMIFS, 1.0);
      ii++;
      /* Invert the band so that the high freqs are first */
      /* This is how SIGPROC stores its data.             */
      {
	int jj;
	unsigned char tempzz=0.0, *loptr, *hiptr;
	loptr = output_samples + 0;
	hiptr = output_samples + fb.nchans - 1;
	for (jj=0; jj<fb.nchans/2; jj++, loptr++, hiptr--){
	  SWAP(*loptr, *hiptr);
	}
      }
      chkfwrite(output_samples, sizeof(unsigned char), fb.nchans, outfile);
    }
    fclose(infiles[filenum-1]);
  }
  fclose(outfile);
  if (scaling)
    free(scalings);
  fprintf(stderr, "Converted and wrote %d samples.\n\n", ii);
  free(infiles);
  return 0;
}
