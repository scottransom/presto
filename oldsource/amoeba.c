/* Note:  This version of the Downhill Simplex Algorithm is     */
/*        a modified version on the NR in C 2ed version         */
/*        by Press et el.  Their copyright statement is below.  */
/*        Modification by Scott M. Ransom, 18 June 1999         */
/* (C) Copr. 1986-92 Numerical Recipes Software 3#1y-i.31-.     */

#include <math.h>
#include "vectors.h"

static double amotry(double **p, double *y, double *psum, int ndim,
		     double (*funk)(double []), int ihi, double fac);

void amoeba(double **p, double *y, int ndim, double ftol,
	    double (*funk)(double []), int *nfunk)
{
  int ii, ihi, ilo, inhi, jj, mpts = ndim + 1;
  double rtol, sum, swap, ysave, ytry, *psum;
  
  psum = gen_dvect(ndim);
  *nfunk = 0;
  for (jj = 0; jj < ndim; jj++){
    for (sum = 0.0, ii = 0; ii < mpts; ii++) sum += p[ii][jj];
    psum[jj] = sum;
  }
  for (;;) {
    ilo = 0;
    ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
    for (ii = 0; ii < mpts; ii++) {
      if (y[ii] <= y[ilo]) ilo = ii;
      if (y[ii] > y[ihi]) {
	inhi = ihi;
	ihi = ii;
      } else if (y[ii] > y[inhi] && ii != ihi) inhi = ii;
    }
    rtol = 2.0 * fabs(y[ihi] - y[ilo]) / \
      (fabs(y[ihi]) + fabs(y[ilo]) + 1.0e-10);
    if (rtol < ftol) {
      swap = y[0]; 
      y[0] = y[ilo]; 
      y[ilo] = swap;
      for (ii = 0; ii < ndim; ii++){
	swap = p[0][ii]; 
	p[0][ii] = p[ilo][ii]; 
	p[ilo][ii] = swap;
      }
      break;
    }
    if (*nfunk >= 5000){
      printf("\n Max # of iterations exceeded in amoeba().  Exiting.\n\n");
      exit(1);
    }
    *nfunk += 2;
    ytry = amotry(p, y, psum, ndim, funk, ihi, -1.0);
    if (ytry <= y[ilo])
      ytry = amotry(p, y, psum, ndim, funk, ihi, 2.0);
    else if (ytry >= y[inhi]) {
      ysave = y[ihi];
      ytry = amotry(p, y, psum, ndim, funk, ihi, 0.5);
      if (ytry >= ysave) {
	for (ii = 0; ii < mpts; ii++) {
	  if (ii != ilo) {
	    for (jj = 0; jj < ndim; jj++)
	      p[ii][jj] = psum[jj] = 0.5 * (p[ii][jj] + p[ilo][jj]);
	    y[ii] = (*funk)(psum);
	  }
	}
	*nfunk += ndim;
	for (jj = 0; jj < ndim; jj++){
	  for (sum = 0.0, ii = 0; ii < mpts; ii++) sum += p[ii][jj];
	  psum[jj] = sum;
        }
      }
    } else --(*nfunk);
  }
  free(psum);
}


double amotry(double **p, double *y, double *psum, int ndim,
	      double (*funk)(double []), int ihi, double fac)
{
  int jj;
  double fac1, fac2, ytry, *ptry;
  
  ptry = gen_dvect(ndim);
  fac1 = (1.0 - fac) / ndim;
  fac2 = fac1 - fac;
  for (jj = 0; jj < ndim; jj++) ptry[jj] = psum[jj] * \
			       fac1 - p[ihi][jj] * fac2;
  ytry = (*funk)(ptry);
  if (ytry < y[ihi]) {
    y[ihi] = ytry;
    for (jj = 0; jj < ndim; jj++) {
      psum[jj] += ptry[jj] - p[ihi][jj];
      p[ihi][jj]=ptry[jj];
    }
  }
  free(ptry);
  return ytry;
}

