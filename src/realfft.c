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

#define DEBUG 1

/*  This program calculates the FFT of a file containing    */
/*  a number of single-precision floats representing        */
/*  real numbers.  (i.e. a normal time series)              */
/*  The data is assumed to be located in the directory,     */
/*  specified in the input filename(s).                     */
/*  Input filename(s) must not include '.dat' or '.fft'     */
/*        suffixes. They will be added by the program.      */
/*  Do not end paths in '/'.                                */
/*  Scratch file(s) are the same size as the input file(s). */
/*  If '-inv' is specified, the file to be transformed      */
/*        should end in '.fft'.  Otherwise, it should end   */
/*        in '.dat'.                                        */

int main(int argc, char *argv[])
{
  FILE **data_files, **scratch_files, **result_files;
  char *data_dir, *scratch_dir, *result_dir, **filenm_roots;
  float *data;
  int ii, status, isign=-1, numfiles;
  long long numdata=0, *files_numdata;
  struct tms runtimes;
  double ttim, stim, utim, tott;
  Cmdline *cmd;
  char data_suffix[]="dat";
  char result_suffix[]="fft";
  char scratch_suffix[]="scratch";

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
  printf("\n\n");
  printf("   Real-Valued Data FFT Program v3.0\n");
  printf("        by Scott M. Ransom\n");
  printf("           27 Sept, 2000\n\n");

  /* Get our file information */

  numfiles = cmd->argc;
  printf("Checking data in %d input file(s):\n", numfiles);
  data_files = malloc(sizeof(FILE *) * numfiles);
  filenm_roots = malloc(sizeof(char *) * numfiles);
  files_numdata = malloc(sizeof(long long) * numfiles);
  {
    int hassuffix=0, new_hassuffix=0;
    char *dir, *filenm, *suffix;
    
    for (ii=0; ii<numfiles; ii++){
      if (ii==0){
	split_path_file(cmd->argv[ii], &data_dir, &filenm);
	hassuffix = split_root_suffix(filenm, &(filenm_roots[ii]), 
				      &suffix);
	if (hassuffix){
	  if (strcmp(suffix, "fft")==0){
	    isign = 1;
	    strcpy(data_suffix, "fft");
	    strcpy(result_suffix, "dat");
	  }
	  free(suffix);
	}
      } else {
	split_path_file(cmd->argv[ii], &dir, &filenm);
	new_hassuffix = split_root_suffix(filenm, &(filenm_roots[ii]), 
					  &suffix);
	if (new_hassuffix && hassuffix){
	  if (strcmp(data_suffix, suffix)){
	    printf("\nAll input files must have the same suffix!\n\n");
	    exit(1);
	  }
	}
	if (new_hassuffix) free(suffix);
	if (strcmp(data_dir, dir)){
	  printf("\nAll input files must be in the same directory!\n\n");
	  exit(1);
	}
	free(dir);
      }
      free(filenm);
      data_files[ii] = chkfopen(cmd->argv[ii], "r");
      files_numdata[ii] = chkfilelen(data_files[ii], sizeof(float));
      numdata += files_numdata[ii];
      if (isign==1)
	printf("   '%s' has %lld complex pts\n", cmd->argv[ii], 
	       files_numdata[ii]/2);
      else
	printf("   '%s' has %lld real pts\n", cmd->argv[ii],
	       files_numdata[ii]);
      fclose(data_files[ii]);
    }
  }
  printf("\nData OK.  There are %lld points.\n\n", numdata);
  exit(0);

  /*  Start the transform sequence  */

  if (numdata > MAXREALFFT) {

    /*  Perform Two-Pass, Out-of-Core, FFT  */

    printf("Performing Out-of-Core Two-Pass FFT on data.\n\n");
    /* printf("Result will be stored in the file \"%s\".\n", resultfilenm); */

    printf("Finished.\n\n");

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);

  } else {

    /* Perform standard FFT for real functions  */

    printf("Performing in-core FFT for real functions on data.\n\n");
    /* printf("FFT will be stored in the file \"%s\".\n\n", resultfilenm); * /

    /* Output the timing information */

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);

  }

  /*
    fftw_print_max_memory_usage();
    fftw_check_memory_leaks();
   */
  for (ii=0; ii<numfiles; ii++)
    free(filenm_roots[ii]);
  free(filenm_roots);
  free(files_numdata);
  free(data_dir);
  free(data_files);
  exit(0);
}
