#include <presto.h>

int read_resid_rec(FILE * file, double *toa, double *obsf);

double doppler(double freq_observed, double voverc)
/* This routine returns the frequency emitted by a pulsar */
/* (in MHz) given that we observe the pulsar at frequency */
/* freq_observed (MHz) while moving with radial velocity  */
/* (in units of v/c) of voverc wrt the pulsar.            */
{
  return freq_observed * (1.0 + voverc);
}


int read_resid_rec(FILE * file, double *toa, double *obsf)
/* This routine reads a single record (i.e. 1 TOA) from */
/* the file resid2.tmp which is written by TEMPO.       */
/* It returns 1 if successful, 0 if unsuccessful.       */
{
  static int l;
  static double d[9];

  fread(&l, sizeof(int), 1, file);
  fread(&d, sizeof(double), 9, file);
  *toa = d[0];
  *obsf = d[4];
  return fread(&l, sizeof(int), 1, file);
}
  

void barycenter(double *topotimes, double *barytimes, \
		double *voverc, long N, char *ra, char *dec, \
		char *obs, char *ephem)
/* This routine uses TEMPO to correct a vector of           */
/* topocentric times (in *topotimes) to barycentric times   */
/* (in *barytimes) assuming an infinite observation         */
/* frequency.  The routine also returns values for the      */
/* radial velocity of the observation site (in units of     */
/* v/c) at the barycentric times.  All three vectors must   */
/* be initialized prior to calling.  The vector length for  */
/* all the vectors is 'N' points.  The RA and DEC (J2000)   */
/* of the observed object are passed as strings in the      */
/* following format: "hh:mm:ss.ssss" for RA and             */
/* "dd:mm:ss.ssss" for DEC.  The observatory site is passed */
/* as a 2 letter ITOA code.  This observatory code must be  */
/* found in obsys.dat (in the TEMPO paths).  The ephemeris  */
/* is either "DE200" or "DE400".                            */
{

  FILE *outfile;
  long i;
  double fobs = 1000.0, femit, dtmp;
  char command[100], temporaryfile[100];

  /* Write the free format TEMPO file to begin barycentering */
  
  strcpy(temporaryfile, "bary.tmp");
  outfile = chkfopen(temporaryfile, "w");
  fprintf(outfile, "C  Header Section         \n");
  fprintf(outfile, "  HEAD                    \n");
  fprintf(outfile, "  PSR                 bary\n");
  fprintf(outfile, "  NPRNT                  2\n");
  fprintf(outfile, "  P0                   1.0\n");
  fprintf(outfile, "  P1                   0.0\n");
  fprintf(outfile, "  CLK            UTC(NIST)\n");
  fprintf(outfile, "  COORD              J2000\n");
  fprintf(outfile, "  RA                    %s\n", ra);
  fprintf(outfile, "  DEC                   %s\n", dec);
  fprintf(outfile, "  DM                   0.0\n");
  fprintf(outfile, "  EPHEM                 %s\n", ephem);
  fprintf(outfile, "C  TOA Section (uses ITAO Format)\n");
  fprintf(outfile, "  TOA\n");

  /* Write the TOAs for infinite frequencies */

  for (i = 0; i < N; i++) {
    fprintf(outfile, "topocent %20.14f 0.00  0.000000000 0.000000 %s\n", \
	    topotimes[i], obs);
  }
  fprintf(outfile, "topocent %20.14f 0.00  0.000000000 0.000000 %s\n", \
	  topotimes[N-1]+0.0001, obs);
  fclose(outfile);

  /* Call TEMPO */

  /* Check the TEMPO *.tmp and *.lis files for errors when done. */

  sprintf(command, "tempo bary.tmp > tempoout_times.tmp");
  system(command);

  /* Now read the TEMPO results */
  
  strcpy(temporaryfile, "resid2.tmp");
  outfile = chkfopen(temporaryfile, "rb");

  /* Read the barycentric TOAs for infinite frequencies */
  
  for (i = 0; i < N; i++) {
    read_resid_rec(outfile, &barytimes[i], &dtmp);
  }
  fclose(outfile);
  
  /* Write the free format TEMPO file to begin barycentering */
  
  strcpy(temporaryfile, "bary.tmp");
  outfile = chkfopen(temporaryfile, "w");
  fprintf(outfile, "C  Header Section         \n");
  fprintf(outfile, "  HEAD                    \n");
  fprintf(outfile, "  PSR                 bary\n");
  fprintf(outfile, "  NPRNT                  2\n");
  fprintf(outfile, "  P0                   1.0\n");
  fprintf(outfile, "  P1                   0.0\n");
  fprintf(outfile, "  CLK            UTC(NIST)\n");
  fprintf(outfile, "  COORD              J2000\n");
  fprintf(outfile, "  RA                    %s\n", ra);
  fprintf(outfile, "  DEC                   %s\n", dec);
  fprintf(outfile, "  DM                   0.0\n");
  fprintf(outfile, "  EPHEM                 %s\n", ephem);
  fprintf(outfile, "C  TOA Section (uses ITAO Format)\n");
  fprintf(outfile, "  TOA\n");

  /* Write the TOAs for finite frequencies */

  for (i = 0; i < N; i++) {
    fprintf(outfile, "topocent %20.14f 0.00  %11.6f 0.000000 %s\n", \
	    topotimes[i], fobs, obs);
  }
  fprintf(outfile, "topocent %20.14f 0.00  %11.6f 0.000000 %s\n", \
	  topotimes[N-1]+0.0001, fobs, obs);
  fclose(outfile);

  /* Call TEMPO */

  /* Insure you check the file tempoout.tmp for  */
  /* errors from TEMPO when complete.            */

  sprintf(command, "tempo bary.tmp > tempoout_vels.tmp");
  system(command);

  /* Now read the TEMPO results */
  
  strcpy(temporaryfile, "resid2.tmp");
  outfile = chkfopen(temporaryfile, "rb");

  /* Determine the radial velocities using the emitted freq */

  for (i = 0; i < N; i++) {
    read_resid_rec(outfile, &dtmp, &femit);
    voverc[i] = femit/fobs - 1.0;
  }
  fclose(outfile);

  /* Cleanup the temp files */

  rename("tempo.lis", "tempoout.tmp");
  remove("resid2.tmp");
  remove("bary.tmp");
  remove("matrix.tmp");
  remove("bary.par");
}

