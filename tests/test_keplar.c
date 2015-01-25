#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/orbint.h"

int main (int argc, char **argv)
{
  double u, e, p, t = 0.0;
  long i, npts;

  npts = atoi(argv[1]);
  p = atof(argv[2]);
  e = atof(argv[3]);
  for (i = 0 ; i <= npts ; i++){
    t = (double) i / (double) npts * p;
    u = keplers_eqn(t, 1.0/p, e, 1.0E-14);
    printf(" t = %10.6f   u = %17.15f\n",t,u);
  }
}
