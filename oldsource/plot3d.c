#include "presto.h"
#include "plot3d.h"

void plot3d_complex(float **array, int nx, int ny, char *palette)
{
  FILE *fp;
  color col;
  int i, j, hii, hij;
  float powargr, powargi;
  float z, maxpower, temp, scale;
  void (*color_ptr) (color *x, double val);

  /*  Determine the color table to use */

  if (strcmp(palette, "heat") == 0)
    color_ptr = heat;
  else if (strcmp(palette, "rainbow") == 0)
    color_ptr = rainbow;
  else if (strcmp(palette, "hue") == 0)
    color_ptr = hue;
  else
    color_ptr = gray;

  /*  Start "Geomview" through "togeomview"  */

  scale = 0.8 * (double) ((nx > (ny / 2)) ? nx : ny / 2);
  system("togeomview -g < /dev/null");
  fp = fopen("/tmp/geomview/OOGL", "w");
  if (!fp) {
    fprintf(stderr, "Unable to open Geomview pipe.\n");
    return;
  }
  /*  Write to the pipe...  */

  fprintf(fp, "appearance { shading smooth}\n");
  fprintf(fp, "CZMESH\n");
  fprintf(fp, "%d %d\n", nx, ny);
  maxpower = 0.0;

  /*  Find the highest z value in the array  */

  for (i = 0; i < ny; i++)
    for (j = 0; j < 2 * nx; j += 2) {
      {
	if ((temp = POWER(array[i][j], array[i][j + 1])) > \
	    maxpower) {
	  maxpower = temp;
	  hii = i;
	  hij = j;
	}
      }
    }

  /*  Write the data to the pipe... */

  for (i = 0; i < ny; i++)
    for (j = 0; j < 2 * nx; j += 2) {
      {
	z = POWER(array[i][j], array[i][j + 1]) / maxpower;
	(*color_ptr) (&col, z);
	z *= scale;
	fprintf(fp, "  %f %5.3f %5.3f %5.3f 1.0", z, \
		col.r, col.g, col.b);
      }
      fprintf(fp, "\n");
    }

  /*  Cleanup and return... */

  fclose(fp);
  return;
}


void quick3d_complex(FILE * file, double fftfreq, char *palette)
{
  long nextbin, startbin;
  int nf = 105, nz = 101, fbetween = 4, fftlen = 1024, offset, zct;
  double bigz = 50.0, df, frac, tmp;
  float **result, *dataptr;
  fcomplex **array;
  long m = 0;

  /*  Get the array to display from the file  */

  df = 1.0 / fbetween;
  frac = modf(fftfreq, &tmp);

  /*  the freq at bin (nf-1)/2 */

  fftfreq = tmp + floor(frac * fbetween + DBLCORRECT) * df;
  startbin = (unsigned long) (fftfreq) - ((nf - 1) / 2) / fbetween;
  offset = (int) ((fftfreq - (unsigned long) (fftfreq)) \
		  * fbetween + DBLCORRECT);
  result = corr_rz_plane_file(file, fbetween, startbin, -bigz, bigz,
			      nz, fftlen, LOWACC, &m, &nextbin);

  /*  Get the part that we need  */

  array = gen_cmatrix(nz, nf);

  for (zct = 0; zct < nz; zct++) {
    dataptr = result[zct] + offset;
    memcpy(array[zct], dataptr, sizeof(fcomplex) * nf);
  }

  /*  Plot it and cleanup */

  plot3d_complex((float) array, nf, nz, palette);

  free(array[0]);
  free(array);
  free(result[0]);
  free(result);
  return;
}
