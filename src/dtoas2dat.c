#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#define WORKLEN 65536
#define SECPERDAY 86400
#define MAXREAD 32768

unsigned long getfilelen(FILE *file, size_t size)
{
  int filenum, rt;
  struct stat buf;
  
  filenum = fileno(file);
  rt = fstat(filenum, &buf);
  if (rt == -1){
    perror("\nError in getfilelen()");
    printf("\n");
    exit(-1);
  }
  return (unsigned long) (buf.st_size / size);
}

int main(int argc, char *argv[])
/* Convert a binary file of double precision MJD TOAs       */
/* into a floating point time series.  The Time series will */
/* have 'N' points with each bin of length 'dt' seconds.    */
{
  long ii, jj, N, ntoas, days=1, numwrites, numtowrite;
  double To, dt, toa, *toaptr, *ddata, lotime, hitime, dtfract;
  float *fdata;
  char outfilenm[200];
  FILE *infile, *outfile;

  if (argc < 4 || argc > 5){
    printf("\nUsage:  dtoas2dat filenm dt N [type]\n\n");
    printf("          'filenm' is the TOA filename\n");
    printf("          'dt' is the time interval for the output data (sec)\n");
    printf("          'N' is the number of output data points\n");
    printf("          'type' is the optional format of the TOAs:\n");
    printf("              either 'd' for days (default) or 's' for sec\n\n");
    exit(0);
  }

  /* Open our files */

  sprintf(outfilenm, "%s.dat", argv[1]);
  printf("\nReading TOAs from '%s'.\n", argv[1]);
  infile = fopen(argv[1], "rb");
  outfile = fopen(outfilenm, "wb");

  /* Get the other command line arguments */

  dt = strtod(argv[2], NULL);
  dtfract = 1.0 / dt;
  N = strtol(argv[3], NULL, 10);
  if (argc==5){
    if (argv[4][0]=='d')
      days = 1;
    else if (argv[4][0]=='s')
      days = 0;
    else {
      printf("\nUnrecognized data type '%s'.\n", argv[4]);
      printf("Must be 'd' for days or 's' for seconds.\n\n");
      exit(-1);
    }
  }

  /* Get the number of TOAs in the TOA file */

  ntoas = getfilelen(infile, sizeof(double));
  printf("Found %ld TOAs.\n", ntoas);

  /* Allocate our data arrays */
 
  ddata = (double *)malloc(sizeof(double) * ntoas);
  fdata = (float *)malloc(sizeof(float) * WORKLEN);
  printf("\nWriting time series of %ld points of\n", N); 
  printf("length %f seconds to '%s'.\n\n", dt, outfilenm); 

  /* Read the TOAs */

  jj = fread(ddata, sizeof(double), ntoas, infile);
  if (jj != ntoas){
    printf("\nError reading TOA file.  Only %ld points read.\n\n", jj);
    exit(-1);
  }

  /* Convert the TOAs to seconds offset from the first TOA */

  To = ddata[0];
  if (days)
    for (ii = 0; ii < ntoas; ii++)
      ddata[ii] = (ddata[ii] - To) * SECPERDAY;
  else
    for (ii = 0; ii < ntoas; ii++)
      ddata[ii] = ddata[ii] - To;
  toaptr = ddata;
  toa = *toaptr;

  /* Determine the number of output writes we need */

  numwrites = (N % WORKLEN) ? N / WORKLEN + 1 : N / WORKLEN;

  /* Loop over the number of writes */

  for (ii = 0; ii < numwrites; ii++){

    /* Determine the beginning and ending times of the output array */

    numtowrite = ((N % WORKLEN) && (ii == (numwrites - 1))) ? 
      N % WORKLEN : WORKLEN;
    lotime = ii * WORKLEN * dt;
    hitime = lotime + numtowrite * dt;

    /* Initialize the output data array to all zeros */

    for (jj = 0; jj < WORKLEN; jj++)
      fdata[ii] = 0.0;

    /* Place any TOAs we need to in the current output array */

    while ((toa >= lotime && toa < hitime)
	   && (toaptr - ddata) < ntoas-1){
      fdata[(int)((toa - lotime) * dtfract)] += 1.0;
      toaptr++;
      toa = *toaptr;
    }

    /* Write the output data */

    fwrite(fdata, sizeof(float), numtowrite, outfile);
  }

  /* Cleanup */

  printf("Done.\n\n");
  fclose(infile);
  fclose(outfile);
  free(fdata);
  free(ddata);
  exit(0);
}

    
