#include "presto.h"
#include "plot2d.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *datfile;
  int startbin, numpoints;
  float *data, *times;
  double dt;
  char filenm[100], psfilenm[100];
  infodata idata;

  if (argc != 4) {
    printf("\nUsage: timeseries startbin numpoints filename \n");
    printf("     This routine displays a '.dat' time series. \n");
    printf("     \n");
    printf("     startbin is the first point to read\n");
    printf("     numpoints is the total number of points to read\n");
    printf("       (0=all)\n");
    exit(0);
  }

  startbin = atoi(argv[1]);
  numpoints = atoi(argv[2]);

  /* open the input file */
  sprintf(filenm, "%s.dat", argv[3]);
  datfile = chkfopen(filenm, "r");

  /* open inf file to get the time resolution of the data */
  readinf(&idata, argv[3]);
  dt = idata.dt;
  if (numpoints==0) numpoints = idata.N;
  printf("    Total number of time bins  = %.0f\n",idata.N);
  printf("    Time resolution (sec)      = %-17.15g\n",idata.dt);

  /* get output ps file name */
  sprintf(psfilenm, "%s.dat_%d-%d.ps", argv[3], 
	  startbin, startbin+numpoints);

  /* read from the file */
  printf("    Reading from file '%s'\n",filenm);
  data = read_float_file(datfile, startbin, numpoints);

  /* calculate times */
  times = gen_freqs(numpoints, startbin*dt, dt);

  /* plot data points */
  printf("    Writing to file '%s'\n",psfilenm);
  cpgstart_ps(psfilenm, "landscape");
  xyline(numpoints, times, data, "Time (s)", "Amplitude", 1);
  cpgend();

  fclose(datfile);
  free(times);
  free(data);
  exit(0);
}
