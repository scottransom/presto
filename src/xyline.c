#include "plot2d.h"
#include <string.h>

void cpgstart_ps(const char *filenm, const char *orientation)
{
  char tmp[100];
  if (0 == strcmp(orientation, "portrait")){
    sprintf(tmp, "%s/VCPS", filenm);
    if (cpgopen(tmp) <= 0)
      exit(EXIT_FAILURE);
  } else {
    sprintf(tmp, "%s/CPS", filenm);
    if (cpgopen(tmp) <= 0)
      exit(EXIT_FAILURE);
  }
}

void cpgstart_x(const char *orientation)
{
  if (cpgopen("/XWIN") <= 0)
    exit(EXIT_FAILURE);
  if (0 == strcmp(orientation, "portrait")){
    cpgpap(8.5, 11.0/8.5);
  } else {
    cpgpap(11.0, 8.5/11.0);
  }
  cpgask(1);
}

void find_min_max_arr(int npts, float *arr, float *min, float *max)
{
  int i;

  *min = 1.1E38;
  *max = -1.1E38;
  for (i = 0; i < npts; i++) {
    if (arr[i] > *max)
      *max = arr[i];
    if (arr[i] < *min)
      *min = arr[i];
  }
}

void dfind_min_max_arr(int npts, double *arr, double *min, double *max)
{
  int i;

  *min = 1.1E100;
  *max = -1.1E100;
  for (i = 0; i < npts; i++) {
    if (arr[i] > *max)
      *max = arr[i];
    if (arr[i] < *min)
      *min = arr[i];
  }
}

void dxyline(int npts, double *x, double *y, const char *xlab, \
	     const char *ylab, int id)
/* Wrapper to plot double precision vectors */
{
  float *fx, *fy;
  long i;

  /* Copy our double vectors to float vectors */

  fx = (float *) malloc(sizeof(float) * npts);
  fy = (float *) malloc(sizeof(float) * npts);

  for (i = 0; i < npts; i++) {
    fx[i] = (float) x[i];
    fy[i] = (float) y[i];
  }

  /* Call xyline */

  xyline(npts, fx, fy, xlab, ylab, id);

  /* Free our memory */

  free(fx);
  free(fy);
}


void xyline(int npts, float *x, float *y, const char *xlab, \
	    const char *ylab, int id)
{
  float xmin, xmax, ymin, ymax;
  float overy, over = 0.1;

  /* Determine min and max values to plot and scaling: */
  find_min_max_arr(npts, x, &xmin, &xmax);
  find_min_max_arr(npts, y, &ymin, &ymax);
  overy = over * (ymax - ymin);
  ymax += overy;
  ymin -= overy;

  /* Setup the plot screen: */
  cpgenv(xmin, xmax, ymin, ymax, 0, 0);

  /* Choose the font: */
  cpgscf(2);

  /* Label the axes: */
  cpglab(xlab, ylab, "");

  /* Add ID line if required */
  if (id == 1)
    cpgiden();

  /* Plot the points: */
  cpgline(npts, x, y);

}


void dxybinned(int npts, double *x, double *y, const char *xlab, \
	       const char *ylab, int id)
/* Wrapper to plot double precision vectors */
{
  float *fx, *fy;
  long i;

  /* Copy our double vectors to float vectors */

  fx = (float *) malloc(sizeof(float) * npts);
  fy = (float *) malloc(sizeof(float) * npts);

  for (i = 0; i < npts; i++) {
    fx[i] = (float) x[i];
    fy[i] = (float) y[i];
  }

  /* Call xyline */

  xybinned(npts, fx, fy, xlab, ylab, id);

  /* Free our memory */

  free(fx);
  free(fy);
}


void xybinned(int npts, float *x, float *y, const char *xlab, \
	      const char *ylab, int id)
{
  float xmin, xmax, ymin, ymax;
  float overy, over = 0.1;

  /* Determine min and max values to plot and scaling: */
  find_min_max_arr(npts, x, &xmin, &xmax);
  find_min_max_arr(npts, y, &ymin, &ymax);
  overy = over * (ymax - ymin);
  ymax += overy;
  ymin -= overy;

  /* Setup the plot screen: */
  cpgenv(xmin, xmax, ymin, ymax, 0, 0);

  /* Choose the font: */
  cpgscf(2);

  /* Label the axes: */
  cpglab(xlab, ylab, "");

  /* Add ID line if required */
  if (id == 1)
    cpgiden();

  /* Plot the points: */
  cpgbin(npts, x, y, 0);
}


void xyline2lab(int npts, float *x, float *y, float *y2, const char *xlab, \
		const char *ylab, const char *ylab2, int id)
{
  float xmin, xmax, ymin, ymax, ymin2, ymax2;
  float overy, over = 0.1;

  /* Determine min and max values to plot and scaling: */
  find_min_max_arr(npts, x, &xmin, &xmax);
  find_min_max_arr(npts, y, &ymin, &ymax);
  find_min_max_arr(npts, y2, &ymin2, &ymax2);
  overy = over * (ymax - ymin);
  ymax += overy;
  ymin -= overy;
  overy = over * (ymax2 - ymin2);
  ymax2 += overy;
  ymin2 -= overy;

  /* Choose the font: */
  cpgscf(2);

  /* Setup the plot screen for the first set of y's: */
  cpgpage();
  cpgvstd();
  cpgswin(xmin, xmax, ymin, ymax);
  cpgbox("BCNST", 0.0, 0, "BNST", 0.0, 0);
  cpgmtxt("B", 3.0, 0.5, 0.5, xlab);
  cpgmtxt("L", 2.6, 0.5, 0.5, ylab);

  /* Plot the points for the 1st y axis: */
  cpgline(npts, x, y);

  /* Setup the plot screen for the second set of y's: */
  cpgvstd();
  cpgswin(xmin, xmax, ymin2, ymax2);
  cpgbox("", 0.0, 0, "CMST", 0.0, 0);
  cpgmtxt("R", 3.0, 0.5, 0.5, ylab2);

  /* Plot the points for the 2nd y axis: */
  cpgline(npts, x, y2);

  /* Add ID line if required */
  if (id == 1)
    cpgiden();

}
