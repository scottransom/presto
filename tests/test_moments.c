#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/vectors.h"

void stats(float *x, int n, double *mean, double *var, 
	   double *skew, double *kurt);

int main(int argc, char *argv[])
{
  long N, i;
  float *x;
  double mean, var, skew, kurt;
  
  N = strtol(argv[1], NULL, 10);
  x = gen_fvect(N);
  printf("Number of points  = %ld\n",N);
  printf("x = {\n");
  for (i = 0 ; i < N ; i++){
    x[i] = (10.0 * rand() / (RAND_MAX + 1.0)) - 5.0;
    printf("%15.12f",x[i]);
    if (i < N-1) printf(",\n");
    else printf("\n");
  }
  printf("}\n");
  stats(x, N, &mean, &var, &skew, &kurt);
  printf("\nmean = %20.10f\n",mean);
  printf("var  = %20.10f\n",var);
  printf("skew = %20.10f\n",skew);
  printf("kurt = %20.10f\n",kurt);
  free(x);
  return 0;
}

