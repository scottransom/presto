#include "plot3d.h"
#include "presto.h"

void threedplot(float **array, int nx, int ny, char *palette)
{
  FILE *fp;
  color col;
  int i, j, hii, hij;
  float powargr, powargi;
  float z, maxpower, temp, scale;
  color(*color_ptr) (float val);

  /*  Determine the color table to use */
  if (strcmp(palette, "heat") == 0)
    color_ptr = heat;
  else if (strcmp(palette, "rainbow") == 0)
    color_ptr = rainbow;
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
	if ((temp = POWER(array[i][j], array[i][j + 1])) > maxpower) {
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
	col = (*color_ptr) (z);
	z *= scale;
	fprintf(fp, "  %f %5.3f %5.3f %5.3f 1.0", z, col.r, col.g, col.b);
      }
      fprintf(fp, "\n");
    }

  /*  Cleanup and return... */
  fclose(fp);
  return;
}
