#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vectors.h"

float *complex_arr(long N);
/* Returns a CFLOAT Python array.       */

void mult_arr(float *arr, float val, long N);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */

void dgenrotate_1d(double *data, long numbins, double bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */

void dstats(double *x, int n, double *mean, double *var, 
            double *skew, double *kurt);
/* For a double precision vector, *x, of length n, this routine  */
/* returns the mean, variance, skewness, and kurtosis of *x.     */

long get_filelen(FILE *file);
/* Return the length of the file *file in bytes without */
/* losing your place in the file.                       */
