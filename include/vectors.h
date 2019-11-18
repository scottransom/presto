#include <stdio.h>
#include <stdlib.h>
#include "fftw3.h"
#include "rawtype.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* Some vector routines by Scott Ransom since he hates the non-zero */
/* offset general vector routines in Numerical Recipes.  ;)         */
/* The routines are in vectors.c                                    */

float *gen_fvect(long length);
/* Generate a floating point vector */

double *gen_dvect(long length);
/* Generate a double precision vector */

fcomplex *gen_cvect(long length);
/* Generate a floating complex number vector */

short *gen_svect(long length);
/* Generate an short integer vector */

int *gen_ivect(long length);
/* Generate an integer vector */

long *gen_lvect(long length);
/* Generate a long integer vector */

unsigned char *gen_bvect(long length);
/* Generate a 'byte' or unsigned character vector */

rawtype *gen_rawvect(long length);
/* Generate a 'rawtype' (can be anything defined above) vector */

unsigned char **gen_bmatrix(long nrows, long ncols);
/* Generate a 'byte' or unsigned char matrix (2 dimensions) */

short **gen_smatrix(long nrows, long ncols);
/* Generate a short int matrix (2 dimensions) */

int **gen_imatrix(long nrows, long ncols);
/* Generate an integer matrix (2 dimensions) */

float **gen_fmatrix(long nrows, long ncols);
/* Generate a floating point matrix (2 dimensions) */

double **gen_dmatrix(long nrows, long ncols);
/* Generate a double precision matrix (2 dimensions) */

fcomplex **gen_cmatrix(long nrows, long ncols);
/* Generate a floating complex number matrix (2 dimensions) */

float ***gen_f3Darr(long nhgts, long nrows, long ncols);
/* Generate a floating point 3D array */

fcomplex ***gen_c3Darr(long nhgts, long nrows, long ncols);
/* Generate a floating complex 3D array */

void vect_free(void *vect);
/* Free a generated vector */ 

/*  Note:  To free memory allocated by these routines simply use  */
/*         the free() function.                                   */
/*                                                                */
/*  Example:                                                      */
/*                                                                */
/*  x = gen_fvect(100);   // Generate a 100 point float vector.   */
/*  free(x);              // Free the vector.                     */
