#include <math.h>
#include "vectors.h"

double *bessjtable(double x, int *N, int tol);

/* This routine will return a table of Integer Order     */
/* Bessel Functions of the First kind, J_n(x).  This     */
/* table will be returned in a dynamically allocated     */
/* double prec vector that includes all orders giving    */
/* an absolute value greater than a tolerance            */
/* determined below.  The table length will be returned  */
/* in the value N, and the max order returned will       */
/* equal (N-1)/2.  (i.e. the jable will go from J_-10(x) */
/* to J_10(x) when the returned N=21.  N is always odd.  */
/* N is determined by finding the lowest magnitude       */
/* value to be returned which is within 'tol' absolute   */
/* magnitudes of the maximum value to be returned.       */

/* The routine will give accuracy out to approx. the     */
/* # of significant figures equal to the square root     */
/* of ACC.  This routine is very fast and accurate       */
/* for virtually all ranges of n and x.                  */

#define ACC   100.0
#define BIGNO 1.0e12
#define BIGNI 1.0e-12
#define MINM  10

double *bessjtable(double x, int *N, int tol)
{
  int i, j, jsum, m, twom, pm, neg = 1, newn;
  double ax, bj, bjm, bjp, sum, tox, norm, max = 0.0, dtol;
  double *bes, *trim;

  ax = fabs(x);
  tol = abs(tol);
  if (x < 0)
    neg = -1;

  /* Do x=0 as a special case */

  if (ax == 0.0) {
    *N = 3;
    bes = gen_dvect(*N);
    for (i = 0; i < *N; i++) {
      bes[i] = 0.0;
    }
    bes[1] = 1.0;
    return bes;

    /* Do all the rest with recursion */

  } else {
    tox = 2.0 / ax;
    m = 2 * ((int) (ax + sqrt(ACC * ax)) / 2);
    if (m < MINM)
      m = MINM;
    twom = 2 * m;
    *N = twom + 1;
    bes = gen_dvect(*N);
    jsum = 0;
    pm = *N - 1;
    bjp = sum = 0.0;
    bj = 1.0;
    pm = *N - 1;
    for (j = m; j > 0; j--, pm--) {
      bjm = j * tox * bj - bjp;
      bjp = bj;
      bj = bjm;

      /* Re-normalize to avoid overflow */

      if (fabs(bj) > BIGNO) {
	bj *= BIGNI;
	bjp *= BIGNI;
	sum *= BIGNI;
	for (i = m; i > j; i--, twom--) {
	  bes[twom] *= BIGNI;
	  bes[*N - 1 - twom] *= BIGNI;
	}
	twom = 2 * m;
      }
      bes[pm] = bjp;
      bes[*N - 1 - pm] = -neg * bjp;
      if (jsum) {
	sum += bj;
	bes[pm] = neg * bjp;
      } else {
	bes[*N - 1 - pm] = bjp;
      }
      jsum = !jsum;
    }
    bes[pm] = bj;
    sum = 2.0 * sum - bj;
    norm = 1.0 / sum;

    /* Normalize and pick the max value */

    for (j = 0; j < *N; j++) {
      bes[j] *= norm;
      if (fabs(bes[j]) > max)
	max = fabs(bes[j]);
    }
  }

  /*  Trim the array to give only the values that fall       */
  /*  within 'tol' order of magnitudes from the max value.   */

  pm = 0;
  max = 1.0 / max;
  dtol = pow(10.0, (double) tol);
  while (max / fabs(bes[pm]) > dtol) {
    pm++;			/* Value of pm is first value to keep */
  }

  /* Create the new array */

  newn = *N - 2 * pm;
  trim = gen_dvect(newn);

  /* Copy the values we need */

  for (i = 0, j = pm; i < newn; i++, j++)
    trim[i] = bes[j];

  /* Erase the old array */

  free(bes);
  *N = newn;

  return trim;
}
#undef ACC
#undef BIGNO
#undef BIGNI
#undef MINM
