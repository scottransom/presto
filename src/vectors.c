#include "vectors.h"

float *gen_fvect(long length)
{
  float *v;

  v = (float *) malloc((size_t) (sizeof(float) * length));
  if (!v) {
    perror("\nError in gen_fvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


double *gen_dvect(long length)
{
  double *v;

  v = (double *) malloc((size_t) (sizeof(double) * length));
  if (!v) {
    perror("\nError in gen_dvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


fcomplex *gen_cvect(long length)
{
  fcomplex *v;

  v = (fcomplex *) malloc((size_t) (sizeof(fcomplex) * length));
  if (!v) {
    perror("\nError in gen_cvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


int *gen_ivect(long length)
{
  int *v;

  v = (int *) malloc((size_t) (sizeof(float) * length));
  if (!v) {
    perror("\nError in gen_ivect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


long *gen_lvect(long length)
{
  long *v;

  v = (long *) malloc((size_t) (sizeof(float) * length));
  if (!v) {
    perror("\nError in gen_lvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


unsigned char *gen_bvect(long length)
{
  unsigned char *v;

  v = (unsigned char *) malloc((size_t) (sizeof(unsigned char) * length));
  if (!v) {
    perror("\nError in gen_bvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


rawtype *gen_rawvect(long length)
/* Allocate a vector of length 'length' rawtype points */
{
  rawtype *v;

  v = (rawtype *) malloc((size_t) (sizeof(rawtype) * length));
  if (!v) {
    perror("\nError in gen_rawvect()");
    printf("\n");
    exit(-1);
  }
  return v;
}


unsigned char **gen_bmatrix(long nrows, long ncols)
{
  /* Note:  To free this matrix, assuming you called it with:    */
  /*             x = gen_bmatrix(10,10);                         */
  /*        all you need to do is the following:                 */
  /*             free(x[0]) ; free(x) ;                          */
  /*        The order is important!                              */

  long i;
  unsigned char **m;

  m = (unsigned char **) malloc((size_t) (nrows * sizeof(unsigned char *)));
  if (!m) {
    perror("\nError in 1st malloc() in gen_bmatrix()");
    printf("\n");
    exit(-1);
  }
  m[0] = (unsigned char *) malloc((size_t) ((nrows * ncols) * sizeof(unsigned char)));
  if (!m[0]) {
    perror("\nError in 2nd malloc() in gen_bmatrix()");
    printf("\n");
    exit(-1);
  }
  for (i = 1; i < nrows; i++)
    m[i] = m[i - 1] + ncols;
  return m;
}


float **gen_fmatrix(long nrows, long ncols)
{
  /* Note:  To free this matrix, assuming you called it with:    */
  /*             x = gen_fmatrix(10,10);                         */
  /*        all you need to do is the following:                 */
  /*             free(x[0]) ; free(x) ;                          */
  /*        The order is important!                              */

  long i;
  float **m;

  m = (float **) malloc((size_t) (nrows * sizeof(float *)));
  if (!m) {
    perror("\nError in 1st malloc() in gen_fmatrix()");
    printf("\n");
    exit(-1);
  }
  m[0] = (float *) malloc((size_t) ((nrows * ncols) * sizeof(float)));
  if (!m[0]) {
    perror("\nError in 2nd malloc() in gen_fmatrix()");
    printf("\n");
    exit(-1);
  }
  for (i = 1; i < nrows; i++)
    m[i] = m[i - 1] + ncols;
  return m;
}


double **gen_dmatrix(long nrows, long ncols)
{
  /* Note:  To free this matrix, assuming you called it with:    */
  /*             x = gen_dmatrix(10,10);                         */
  /*        all you need to do is the following:                 */
  /*             free(x[0]) ; free(x) ;                          */
  /*        The order is important!                              */

  long i;
  double **m;

  m = (double **) malloc((size_t) (nrows * sizeof(double *)));
  if (!m) {
    perror("\nError in 1st malloc() in gen_dmatrix()");
    printf("\n");
    exit(-1);
  }
  m[0] = (double *) malloc((size_t) ((nrows * ncols) * sizeof(double)));
  if (!m[0]) {
    perror("\nError in 2nd malloc() in gen_dmatrix()");
    printf("\n");
    exit(-1);
  }
  for (i = 1; i < nrows; i++)
    m[i] = m[i - 1] + ncols;
  return m;
}


fcomplex **gen_cmatrix(long nrows, long ncols)
{
  /* Note:  To free this matrix, assuming you called it with:    */
  /*             x = gen_cmatrix(10,10);                         */
  /*        all you need to do is the following:                 */
  /*             free(x[0]) ; free(x) ;                          */
  /*        The order is important!                              */

  long i;
  fcomplex **m;

  /* allocate pointers to rows */

  m = (fcomplex **) malloc((size_t) (nrows * sizeof(fcomplex *)));
  if (!m) {
    perror("\nError in 1st malloc() in gen_cmatrix()");
    printf("\n");
    exit(-1);
  }
  /* allocate rows and set pointers to them */

  m[0] = (fcomplex *) malloc((size_t) ((nrows * ncols) * sizeof(fcomplex)));
  if (!m[0]) {
    perror("\nError in 2nd malloc() in gen_cmatrix()");
    printf("\n");
    exit(-1);
  }
  for (i = 1; i < nrows; i++)
    m[i] = m[i - 1] + ncols;

  /* return pointer to array of pointers to rows */

  return m;
}
