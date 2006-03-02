#include "presto.h"
#include "plot2d.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   FILE *fftfile;
   int startbin, numpoints;
   fcomplex *data;
   float *powers = NULL, *freqs;
   double T, df;
   double powargr, powargi;
   int ii = 0;
   char filenm[100], psfilenm[100];
   infodata idata;


   if (argc != 4) {
      printf("\nUsage: pspec startbin numpoints filename \n");
      printf("     This routine displays the power spectrum\n");
      printf("     as calculated by realfft.\n");
      printf("     \n");
      printf("     startbin is the first point to read\n");
      printf("     numpoints is the total number of points to read\n");
      printf("       (0=all)\n");
      exit(0);
   }

   startbin = atoi(argv[1]);
   numpoints = atoi(argv[2]);

   /* open the input file */
   sprintf(filenm, "%s.fft", argv[3]);
   fftfile = chkfopen(filenm, "r");

   /* open inf file to get the freq resolution of the fft */
   readinf(&idata, argv[3]);
   T = idata.N * idata.dt;
   df = 1.0 / T;
   if (numpoints == 0)
      numpoints = idata.N / 2;
   printf("    Total number of time bins  = %.0f\n", idata.N);
   printf("    Time resolution (sec)      = %-17.15g\n", idata.dt);

   /* get output ps file name */
   sprintf(psfilenm, "%s.fft_%d-%d.ps", argv[3], startbin, startbin + numpoints);
   printf("    Reading from file '%s'\n", filenm);

   powers = gen_fvect(numpoints);

   /* read from the file */
   data = read_fcomplex_file(fftfile, startbin, numpoints);
   for (ii = 0; ii < numpoints; ii++)
      powers[ii] = POWER(data[ii].r, data[ii].i);

   /* calculate frequencies */
   freqs = gen_freqs(numpoints, startbin * df, df);

   /* plot powers */
   printf("    Writing to file '%s'\n", psfilenm);
   cpgstart_ps(psfilenm, "landscape");
   xyline(numpoints, freqs, powers, "Frequency (Hz)", "Raw Spectral Power", 1);
   cpgend();

   fclose(fftfile);
   free(data);
   free(powers);
   free(freqs);
   exit(0);
}
