#include "cpgplot.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

void find_min_max_arr(int npts, float *arr, float *min, float *max);
void dfind_min_max_arr(int npts, double *arr, double *min, double *max);
void xyline(int npts, float *x, float *y, const char *xlab, \
	    const char *ylab, int id);
void dxyline(int npts, double *x, double *y, const char *xlab, \
	     const char *ylab, int id);
void xybinned(int npts, float *x, float *y, const char *xlab, \
	      const char *ylab, int id);
void dxybinned(int npts, double *x, double *y, const char *xlab, \
	       const char *ylab, int id);
void xyline2lab(int npts, float *x, float *y, float *y2, const char *xlab, \
		const char *ylab, const char *ylab2, int id);
void powerplot(int npts, float *freqs, float *amp, float norm, int id);
void cpgstart_ps(const char *filenm, const char *orientation);
void cpgstart_x(const char *orientation);
double plot_power(double rl, double im);
void multi_prof_plot(int proflen, int numprofs, double *profiles, \
		     double *sumprof, const char *xlbl, \
		     double loly, double ldy, const char *lylbl, \
		     double lory, double rdy, const char *rylbl);
void plot_profile(int proflen, float *profile, const char *title, \
		  const char *probtxt, const char *foldtxt, \
		  int showerr, float *errors, int showid);
void plot_spectrum(fcomplex *spect, int numspect, 
		   double lor, double dr, double T, double average);
/* Plot a chunk of the Fourier power spectrum normalized by average  */
/* The viewing area is left defined with the xvals as _bins_.        */
int cpgnice_output_2(char *out, double val, double err, int len);
/* Return a string that has the same formatting as       */
/* nice_output_2(), but for PGPLOT.  This way, exponents */
/* are actually in superscript!  Woo-hoo!                */


