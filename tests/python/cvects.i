%module cvects
%{
#include "cvects.h"
%}

%include numpy.i

float_complex *complex_arr(long ARRAYLEN);
/* Returns a CFLOAT Python array.       */

void mult_arr(float *IN_1D_CFLOAT, float val, long N);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */

void dgenrotate_1d(double *IN_1D_DOUBLE, long numbins, double bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */

void dstats(double *IN_1D_DOUBLE, int n, double *OUTDOUBLE, double *OUTDOUBLE, 
            double *OUTDOUBLE, double *OUTDOUBLE);
/* For a double precision vector, *x, of length n, this routine  */
/* returns the mean, variance, skewness, and kurtosis of *x.     */

long get_filelen(FILE *file);
/* Return the length of the file *file in bytes without */
/* losing your place in the file.                       */
