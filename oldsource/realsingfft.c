#include <time.h>
#include "nrutil.h"
#include "nr.h"
#include "chkio.h"
#include "ransomfft.h"

#define NRANSI

void realsingfft(FILE * bigfft[5], long numdata, int isign, \
		 char *inpath, char *outpath)
/*  This routine performs a forward 1-D FFT on a HUGE real data set  */
/*    numdata is the number of real data pts  */
{
  FILE *fileloin, *filehiin, *fileloout, *filehiout;
  char command[100], bigfftnm[4][100], resultfile[100], tempfile[100];
  int endi, status;
  long nworkflts, nloio;
  unsigned long nn[2], n, hifileptr, numchunks;
  unsigned long i, i1, i2, i3, i4, ct;
  float *datahi, *datalo, *sav, tmp;
  float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  float h2rwr, h2rwi, h2iwr, h2iwi;
  double filetdiff, wr, wi, wpr, wpi, wtemp, theta;
  struct stat values1, values3;

  n = numdata;
  nn[1] = numdata >> 1;
  nworkflts = MAXREALFFT >> 1;

  /*  Set up trig recursion  */

  isign = -isign;
  theta = PI / (double) (n >> 1);
  if (isign == 1) {
    c2 = -0.5;
    fourfs(bigfft, nn, 1, -1);
    theta = -theta;
  } else {
    c2 = 0.5;
  }
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;

  /*  Reset the files. */

  sprintf(bigfftnm[0], "%sbigfft.000", inpath);
  sprintf(bigfftnm[1], "%sbigfft.001", inpath);
  sprintf(bigfftnm[2], "%sbigfft.002", outpath);
  sprintf(bigfftnm[3], "%sbigfft.003", outpath);

  for (ct = 1; ct <= 4; ct++) {
    fclose(bigfft[ct]);
    bigfft[ct] = chkfopen(bigfftnm[ct - 1], "rb+");
  }

  stat(bigfftnm[0], &values1);
  stat(bigfftnm[2], &values3);
  filetdiff = difftime(values3.st_ctime, values1.st_ctime);

  if (filetdiff > 0) {
    fileloin = bigfft[3];
    filehiin = bigfft[4];
    fileloout = bigfft[1];
    filehiout = bigfft[2];
  } else {
    fileloin = bigfft[1];
    filehiin = bigfft[2];
    fileloout = bigfft[3];
    filehiout = bigfft[4];
  }

  /*  Rearrange the output. */

  /*  Initialize variables and files */

  datalo = vector(1, nworkflts);
  datahi = vector(1, nworkflts);
  for (ct = 1; ct <= 4; ct++)
    rewind(bigfft[ct]);
  nloio = nworkflts;
  endi = nworkflts >> 1;

  /*  Get initial values from the data arrays and save to results */

  sav = vector(1, 2);
  chkfread(&sav[1], sizeof(float), 2, fileloin);
  sav[1] = (tmp = sav[1]) + sav[2];
  sav[2] = tmp - sav[2];
  chkfwrite(&sav[1], sizeof(float), 2, fileloout);
  chkfread(&sav[1], sizeof(float), 2, filehiin);

  /*  Start the conversion loop through the files */

  numchunks = ((n >> 1) / nworkflts) + 1;
  for (ct = 1, hifileptr = 2 * (numdata - MAXREALFFT); \
       ct < numchunks; ct++, hifileptr -= 2 * MAXREALFFT) {

    /* Read data from remaining bigfft files  */

    chkfseek(filehiin, (long) hifileptr, SEEK_SET);
    if (ct == numchunks - 1)
      nloio -= 2;
    chkfread(&datalo[1], sizeof(float), (unsigned long) nloio, fileloin);
    chkfread(&datahi[1], sizeof(float), (unsigned long) nworkflts, filehiin);

    /*  Do the conversion */

    if (ct == numchunks - 1)
      endi--;
    for (i = 1; i <= (unsigned long) endi; i++) {
      i2 = i + i;
      i1 = i2 - 1;
      i3 = nworkflts - i1;
      i4 = i3 + 1;
      h1r = c1 * (datalo[i1] + datahi[i3]);
      h1i = c1 * (datalo[i2] - datahi[i4]);
      h2r = -c2 * (datalo[i2] + datahi[i4]);
      h2i = c2 * (datalo[i1] - datahi[i3]);
      h2rwr = h2r * wr;
      h2rwi = h2r * wi;
      h2iwr = h2i * wr;
      h2iwi = h2i * wi;
      datalo[i1] = h1r + h2rwr - h2iwi;
      datalo[i2] = h1i + h2iwr + h2rwi;
      datahi[i3] = h1r - h2rwr + h2iwi;
      datahi[i4] = -h1i + h2iwr + h2rwi;
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    if (ct == numchunks - 1) {
      datahi[1] = sav[1];
      datahi[2] = sav[2];
    }
    /*  Write the corrected data */

    chkfwrite(&datalo[1], sizeof(float), (unsigned long) nloio, fileloout);
    chkfseek(filehiout, (long) hifileptr, SEEK_SET);
    chkfwrite(&datahi[1], sizeof(float), (unsigned long) nworkflts, filehiout);
  }

  /*  Close files, move files, and cleanup */

  for (ct = 1; ct <= 4; ct++)
    fclose(bigfft[ct]);

  /* Note backward use of paths to get paths straight... */

  sprintf(resultfile, "%ssingresult.fft", inpath);
  sprintf(tempfile, "%stempresult.fft", outpath);

  if (filetdiff > 0) {

    if (remove(bigfftnm[2]))
      nrerror("Can't delete file bigfft.002\n");
    if (remove(bigfftnm[3]))
      nrerror("Can't delete file bigfft.003\n");

    sprintf(command, "mv %s %s\n", bigfftnm[0], resultfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    sprintf(command, "mv %s %s\n", bigfftnm[1], tempfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    sprintf(command, "cat %s >> %s\n", tempfile, resultfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (cat) failed.\n");
      exit(1);
    }
    if (remove(tempfile))
      nrerror("Can't delete file tempresult.fft\n");

  } else {

    if (remove(bigfftnm[0]))
      nrerror("Can't delete file bigfft.000\n");
    if (remove(bigfftnm[1]))
      nrerror("Can't delete file bigfft.001\n");

    sprintf(command, "mv %s %s\n", bigfftnm[2], resultfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    sprintf(command, "mv %s %s\n", bigfftnm[3], tempfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (mv) failed.\n");
      exit(1);
    }
    sprintf(command, "cat %s >> %s\n", tempfile, resultfile);
    if ((status = (system(command))) == -1 || status == 127) {
      printf("System call (cat) failed.\n");
      exit(1);
    }
    if (remove(tempfile))
      nrerror("Can't delete file tempresult.fft\n");

  }
  nworkflts = MAXREALFFT >> 1;
  free_vector(datalo, 1, nworkflts);
  free_vector(datahi, 1, nworkflts);
  free_vector(sav, 1, 2);
}

#undef NRANSI
