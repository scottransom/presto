/*     Real-Valued Data FFT Program        */
/*          by Scott Ransom                */

#include <time.h>
#include <sys/times.h>
#include "ransomfft.h"
#include "misc_utils.h"
#include "chkio.h"
#include "vectors.h"

int main(int argc, char *argv[])
{
  float *data;
  int status, isign, numfileops, ct;
  long dataptr, blocksize, nwrite, next2ton, twon;
  unsigned long numdata;
  FILE *datafile, *bigfft[5], *temperr;
  char datafilenm[100], resultfilenm[100], bigfftnm[4][100];
  char inpath[100], outpath[100], command[100];
  struct tms runtimes;
  double ttim, utim, stim, tott;

  tott = times(&runtimes) / (double) CLOCKS_PER_SEC;
  printf("\n\n");
  printf("   Real-Valued Data FFT Program v1.3\n");
  printf("        by Scott M. Ransom\n");
  printf("           23 July, 1997\n\n");

  if ((argc > 1) && (argc < 7)) {

    /*      Open and check data file.    */

    if (argc == 3) {
      sprintf(datafilenm, "%s.", argv[2]);
      sprintf(resultfilenm, "%s.", argv[2]);
      sprintf(inpath, "./");
      sprintf(outpath, "./");
    }
    if (argc == 4) {
      sprintf(datafilenm, "%s/%s.", argv[3], argv[2]);
      sprintf(resultfilenm, "%s/%s.", argv[3], argv[2]);
      sprintf(inpath, "%s/", argv[3]);
      sprintf(outpath, "./");
    }
    if (argc == 5) {
      sprintf(datafilenm, "%s/%s.", argv[3], argv[2]);
      sprintf(resultfilenm, "%s/%s.", argv[3], argv[2]);
      sprintf(inpath, "%s/", argv[3]);
      sprintf(outpath, "%s/", argv[4]);
    }
    isign = atoi(argv[1]);

    /*  Add the appropriate suffix to the filenames. */

    if (isign == 1) {
      strcat(datafilenm, "fft");
      strcat(resultfilenm, "dat");
    } else {
      strcat(datafilenm, "dat");
      strcat(resultfilenm, "fft");
    }

    /*  Check the input data set...  */

    printf("Checking data in \"%s\".\n", datafilenm);
    datafile = chkfopen(datafilenm, "rb");

    /* # of real data points */

    numdata = get_filelen(datafile, sizeof(float));
    next2ton = next2_to_n(numdata);
    if (numdata != next2ton) {
      printf("\nNumber of data pts must be an integer power of two,\n");
      printf("     or data must be single precision floats.\n");
      printf("Exiting.\n\n");
      fclose(datafile);
      exit(1);
    }
    printf("Data OK.  There are %ld points.\n\n", numdata);
    fclose(datafile);
  } else {
    printf("\nUsage:  realfft sign datafilename [ data path] ");
    printf("[scratch path]\n\n");
    printf("  This program calculates the FFT of a file containing\n");
    printf("  a power-of-two number of floats representing\n");
    printf("  real numbers.  (i.e. a normal time series)\n\n");
    printf("  THE INPUT FILE WILL BE OVERWRITTEN.\n\n");
    printf("  \"sign\" is the sign of the exponential during the FFT. \n");
    printf("  (i.e. -1 is the standard forward transform, 1 is the\n");
    printf("  inverse transform times N (the number of floats))\n");
    printf("  If both paths are omitted, the data is assumed to be\n");
    printf("  located in the working directory, and the scratch space\n");
    printf("  if needed will be located there. If only one path is\n");
    printf("  given, both input and scratch files will reside there.\n\n");
    printf("  Notes:\n");
    printf("    - datafilename must not include \".dat\" or \".fft\"\n");
    printf("        suffix. It will be added by the program.\n");
    printf("    - Do not end paths in \"/\".\n");
    printf("    - The scratch file is the same size as the input file.\n");
    printf("    - If \"sign\"=1, the file to be transformed should end\n");
    printf("        with \".fft\".  Otherwise, it should end with\n");
    printf("        \".dat\".\n\n");
    exit(0);
  }

  /*     Start the transform sequence               */

  if (numdata > MAXREALFFT) {

    /*  Perform singleton's Method FFT  */

    printf("Performing Singleton's Method FFT on data.\n");
    printf("FFT will be stored in the file \"%s\".\n", resultfilenm);

    /*  Initialize files. */

    printf("\nPrepping data files...\n");
    twon = (numdata << 1);
    blocksize = (twon >= 1024) ? 1024 : twon;
    nwrite = twon / blocksize;
    temperr = freopen("/dev/null", "r+", stderr);

    sprintf(command,
	    "dd if=%s of=%sbigfft.000 bs=%d count=%d\n",
	    datafilenm, outpath, (int) blocksize, (int) nwrite);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (dd) failed.\n");
      exit(1);
    }
    sprintf(command,
	    "dd if=%s of=%sbigfft.001 bs=%d skip=%d count=%d\n", \
	    datafilenm, outpath, (int) blocksize, (int) nwrite, \
	    (int) nwrite);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (dd) failed.\n");
      exit(1);
    }
    temperr = freopen("/dev/tty", "r+", stderr);
    sprintf(bigfftnm[0], "%sbigfft.000", outpath);
    bigfft[1] = chkfopen(bigfftnm[0], "rb+");
    sprintf(bigfftnm[1], "%sbigfft.001", outpath);
    bigfft[2] = chkfopen(bigfftnm[1], "rb+");
    sprintf(bigfftnm[2], "%sbigfft.002", inpath);
    bigfft[3] = chkfopen(bigfftnm[2], "wb+");
    sprintf(bigfftnm[3], "%sbigfft.003", inpath);
    bigfft[4] = chkfopen(bigfftnm[3], "wb+");
    printf("Transforming...\n");
    sprintf(command, "rm -f %s\n", datafilenm);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (rm) failed.\n");
      exit(1);
    }
    /*     Call singletons routine  */

    realsingfft(bigfft, numdata, isign, outpath, inpath);
    sprintf(command, "mv %ssingresult.fft %s\n", outpath, resultfilenm);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    printf("Finished.\n\n");

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLOCKS_PER_SEC - tott;
    utim = runtimes.tms_utime / (double) CLOCKS_PER_SEC;
    stim = runtimes.tms_stime / (double) CLOCKS_PER_SEC;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);

  } else {

    /*  Perform standard FFT for real functions  */

    printf("Performing in-core FFT for real functions on data.\n\n");
    printf("FFT will be stored in the file \"%s\".\n\n", resultfilenm);

    /*      Open the data and fft results files.    */

    datafile = chkfopen(datafilenm, "rb+");

    /*    Read data from file to data array   */

    printf("Reading data.\n\n");
    data = gen_fvect(numdata);
    if (numdata > FILEBUFFSIZE) {
      numfileops = numdata / FILEBUFFSIZE;
      dataptr = 0L;
      for (ct = 0; ct < numfileops; ct++, dataptr += FILEBUFFSIZE) {
	chkfread(&data[dataptr], sizeof(float), FILEBUFFSIZE,
		 datafile);
      }
    } else {
      chkfread(data, sizeof(float), (unsigned long) numdata, datafile);
    }

    rewind(datafile);

    /*    Start and time the transform   */

    printf("Transforming...\n");
    realfft(data, numdata, isign);
    printf("Finished.\n\n");

    /*    Write data from FFT array to file  */

    printf("Writing FFT data.\n\n");
    if (numdata > FILEBUFFSIZE) {
      numfileops = numdata / FILEBUFFSIZE;
      dataptr = 0L;
      for (ct = 0; ct < numfileops; ct++, dataptr += FILEBUFFSIZE) {
	chkfwrite(&data[dataptr], sizeof(float), FILEBUFFSIZE,
		  datafile);
      }
    } else {
      chkfwrite(data, sizeof(float), (unsigned long) numdata, datafile);
    }

    free(data);

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLOCKS_PER_SEC - tott;
    utim = runtimes.tms_utime / (double) CLOCKS_PER_SEC;
    stim = runtimes.tms_stime / (double) CLOCKS_PER_SEC;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);

    sprintf(command, "mv %s %s\n", datafilenm, resultfilenm);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    fclose(datafile);

  }

  exit(0);

}

#undef NRANSI
