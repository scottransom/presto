/*     Real-Valued Data FFT Program        */
/*          by Scott Ransom                */
/*            Version 3.0                  */

#include <time.h>
#include <sys/times.h>
#include "clk_tck.h"
#include "misc_utils.h"
#include "chkio.h"
#include "ransomfft.h"
#include "vectors.h"
#include "realfft_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

long long llnext2_to_n(long long x)
/* Return the first value of 2^n >= x */
{
  long long i = 1;
  
  while (i < x) i <<= 1;
  return i;
}

/*  This program calculates the FFT of a file containing    */
/*  a number of single-precision floats representing        */
/*  real numbers.  (i.e. a normal time series)              */
/*  The data is assumed to be located in the directory,     */
/*  specified in the input filename(s).                     */
/*  Input filename(s) must include '.dat' or '.fft'         */
/*        suffixes.  The output file(s) will have the       */
/*        appropriate other suffix.                         */
/*  Do not end paths in '/'.                                */
/*  Scratch file(s) are the same size as the input file(s). */
/*  If '-inv' is specified, the file to be transformed      */
/*        should end in '.fft'.  Otherwise, it should end   */
/*        in '.dat'.                                        */

int main(int argc, char *argv[])
{
  multifile *datfile, *tmpfile=NULL, *outfile;
  char *datdir, **datfilenms, **tmpfilenms, **outfilenms;
  float *data;
  int ii, isign=-1, numfiles;
  long long numdata=0, maxfilelen=0;
  struct tms runtimes;
  double ttim, stim, utim, tott;
  Cmdline *cmd;
  char datsuffix[]="dat";
  char outsuffix[]="fft";
  char tmpsuffix[]="tmp";

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

  tott = times(&runtimes) / (double) CLK_TCK;
  printf("\n");
  printf("   Real-Valued Data FFT Program v3.0\n");
  printf("        by Scott M. Ransom\n");
  printf("            3 Oct, 2000\n\n");

  /* Get our file information */

  numfiles = cmd->argc;
  printf("Checking data in %d input file(s):\n", numfiles);
  datfilenms = (char **)malloc(numfiles * sizeof(char *));
  tmpfilenms = (char **)malloc(numfiles * sizeof(char *));
  outfilenms = (char **)malloc(numfiles * sizeof(char *));
  {
    int hassuffix=0, new_hassuffix=0, filenmlen;
    char *dir, *filenm, *root, *suffix;
    
    for (ii=0; ii<numfiles; ii++){
      if (ii==0){
	split_path_file(cmd->argv[0], &datdir, &filenm);
	hassuffix = split_root_suffix(filenm, &root, &suffix);
	if (hassuffix){
	  if (strcmp(suffix, "fft")==0){
	    isign = 1;
	    strcpy(datsuffix, "fft");
	    strcpy(outsuffix, "dat");
	  }
	  free(suffix);
	}
	free(filenm);
      } else {
	split_path_file(cmd->argv[ii], &dir, &filenm);
	new_hassuffix = split_root_suffix(filenm, &root, &suffix);
	if (new_hassuffix && hassuffix){
	  if (strcmp(datsuffix, suffix)){
	    printf("\nAll input files must have the same suffix!\n\n");
	    exit(1);
	  }
	}
	if (strcmp(datdir, dir)){
	  printf("\nAll input files must be in the same directory!\n\n");
	  exit(1);
	}
	if (new_hassuffix) free(suffix);
	free(dir);
	free(filenm);
      }
      filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
      datfilenms[ii] = (char *)calloc(filenmlen, 1);
      sprintf(datfilenms[ii], "%s/%s.%s", datdir, root, datsuffix);
      if (cmd->tmpdirP){
	filenmlen = strlen(cmd->tmpdir) + 1 + strlen(root) + 5;
	tmpfilenms[ii] = (char *)calloc(filenmlen, 1);
	sprintf(tmpfilenms[ii], "%s/%s.%s", cmd->tmpdir, root, tmpsuffix);
      } else {
	filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
	tmpfilenms[ii] = (char *)calloc(filenmlen, 1);
	sprintf(tmpfilenms[ii], "%s/%s.%s", datdir, root, tmpsuffix);
      }
      if (cmd->outdirP){
	filenmlen = strlen(cmd->outdir) + 1 + strlen(root) + 5;
	outfilenms[ii] = (char *)calloc(filenmlen, 1);
	sprintf(outfilenms[ii], "%s/%s.%s", cmd->outdir, root, outsuffix);
      } else {
	filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
	outfilenms[ii] = (char *)calloc(filenmlen, 1);
	sprintf(outfilenms[ii], "%s/%s.%s", datdir, root, outsuffix);
      }
      free(root);
    }
    free(datdir);
  }

  /* Force a forward or inverse transform.   */
  /* Note that the suffixes do _not_ change! */

  if (cmd->forwardP)
    isign = -1;
  if (cmd->inverseP)
    isign = 1;
  if (cmd->diskfftP && cmd->memfftP){
    printf("\nYou cannot take both an in- and out-of-core FFT!\n\n");
    exit(1);
  }

  /* Open and check data and output files */

  datfile = fopen_multifile(numfiles, datfilenms, "r", 0);
  for (ii=0; ii<numfiles; ii++){
    if (datfile->filelens[ii] > maxfilelen) 
      maxfilelen = datfile->filelens[ii];
  }
  outfile = fopen_multifile(numfiles, outfilenms, "w", maxfilelen);
  for (ii=0; ii<numfiles; ii++)
    printf("   %d:  '%s'\n", ii+1, datfile->filenames[ii]);
  numdata = datfile->length / sizeof(float);
  if (isign==-1){
    if (datfile->length % sizeof(float)){
      printf("\nInput file does not contain the correct number of\n");
      printf("   bytes for it to be floating point data!\n\n");
      exit(1);
    }
    printf("\nData OK.  There are %lld floats.\n\n", numdata);
  } else {
    if (datfile->length % sizeof(fcomplex)){
      printf("\nInput file does not contain the correct number of\n");
      printf("   bytes for it to be single precision complex data!\n\n");
      exit(1);
    }
    printf("\nData OK.  There are %lld complex points.\n\n", numdata/2);
  }
  printf("Result will be written to %d output file(s):\n", numfiles);
  for (ii=0; ii<numfiles; ii++)
    printf("   %d:  '%s'\n", ii+1, outfilenms[ii]);

  /*  Start the transform sequence  */
  
  if ((numdata > MAXREALFFT || cmd->diskfftP) && !cmd->memfftP){

    /*  Perform Two-Pass, Out-of-Core, FFT  */

    printf("\nPerforming Out-of-Core Two-Pass FFT on data.\n");

    /* Make sure the number of points is a power-of-two! */

    if (llnext2_to_n(numdata) != numdata){
      printf("\nMust have a power-of-two number of points for\n");
      printf("   an Out-of-Core FFT!  Exiting.\n\n");
      exit(1);
    }

    tmpfile = fopen_multifile(numfiles, tmpfilenms, "w", maxfilelen);
    if (isign==1) {
      realfft_scratch_inv(datfile, tmpfile, numdata);
    } else {
      realfft_scratch_fwd(datfile, tmpfile, numdata);
    }
    fclose_multifile(tmpfile);
    for (ii=0; ii<numfiles; ii++)
      remove(tmpfilenms[ii]);
    
  } else {

    /* Perform standard FFT for real functions  */
    
    printf("\nPerforming In-Core FFT on data:\n");
    printf("   Reading.\n");
    data = gen_fvect(numdata);
    fread_multifile(data, sizeof(float), numdata, datfile);
    printf("   Transforming.\n");
    realfft(data, numdata, isign);
    printf("   Writing.\n");
    fwrite_multifile(data, sizeof(float), numdata, outfile);
  }

  /* Close our input and output files */

  fclose_multifile(datfile);
  fclose_multifile(outfile);
  
  /* Delete the input files if requested */
  
  if (cmd->deleteP){
    for (ii=0; ii<numfiles; ii++)
      remove(datfilenms[ii]);
  }

  /* Output the timing information */
  
  printf("Finished.\n\n");
  printf("Timing summary:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("  CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	 ttim, utim, stim);
  printf("  Total time elapsed:  %.3f sec\n\n", tott);

  /*
    fftw_print_max_memory_usage();
    fftw_check_memory_leaks();
  */

  for (ii=0; ii<numfiles; ii++){
    free(datfilenms[ii]);
    free(tmpfilenms[ii]);
    free(outfilenms[ii]);
  }
  free(datfilenms);
  free(tmpfilenms);
  free(outfilenms);
  exit(0);
}
