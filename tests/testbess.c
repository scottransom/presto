#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *bessjtable(double x, int *N, int tol);
/* Generates the bessel functions needed by fmresp.                    */

#define TOL -15

int main(int argc, char *argv[])
{
	int N, i, offset;
	double *bessj, x;
	
	x = strtod(argv[1], NULL);
	printf("x = %f\n",x);
	printf("-------------------------------\n");
	bessj = bessjtable(x, &N, TOL);
	offset = (N-1)/2;
	for (i = 0; i < N; i++)
		printf("J_%-4d =  % 11.8f\n", i-offset, bessj[i]);
	free(bessj);
}
		