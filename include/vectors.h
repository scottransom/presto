#include <stdlib.h>
#include <stdio.h>
#include "rawtype.h"

/* Some vector routines by Scott Ransom since he hates the non-zero */
/* offset general vector routines in Numerical Recipes.  ;)         */
/* The routines are in vectors.c                                    */

float *gen_fvect(long length);
/* Generate a floating point vector */

double *gen_dvect(long length);
/* Generate a double precision vector */

fcomplex *gen_cvect(long length);
/* Generate a floating complex number vector */

int *gen_ivect(long length);
/* Generate an integer vector */

long *gen_lvect(long length);
/* Generate a long integer vector */

unsigned char *gen_bvect(long length);
/* Generate a 'byte' or unsigned character vector */

rawtype *gen_rawvect(long length);
/* Generate a 'rawtype' (can be anything defined above) vector */

float **gen_fmatrix(long nrows, long ncols);
/* Generate a floating point matrix (2 dimensions) */

double **gen_dmatrix(long nrows, long ncols);
/* Generate a double precision matrix (2 dimensions) */

fcomplex **gen_cmatrix(long nrows, long ncols);
/* Generate a floating complex number matrix (2 dimensions) */


/*  Note:  To free memory allocated by these routines simply use  */
/*         the free() function.                                   */
/*                                                                */
/*  Example:                                                      */
/*                                                                */
/*  x = gen_fvect(100);   // Generate a 100 point float vector.   */
/*  free(x);              // Free the vector.                     */
