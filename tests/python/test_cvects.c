#include "cvects.h"

int main(void)
{
  long i, N=10;
  float *x, *y;
  double *z;

  x = complex_arr(N);
  y = complex_arr(N);
  mult_arr(y, 3.0, 2*N);
  for(i=0; i<2*N; i++){
    x[i] += 2.0;
    printf("x[%ld] = %f   y[%ld] = %f\n", i, x[i], i, y[i]);
  }
  free(x);
  free(y);
  z = gen_dvect(10);
  for (i = 0; i<10; i++){
    z[i] = (double) i;
    printf("x[%ld] = %f\n", i, z[i]);
  }
  printf("rotate left...\n");
  dgenrotate_1d(z, 10, 3.5);
  for (i = 0; i<10; i++){
    printf("z[%ld] = %f\n", i, z[i]);
  }
  printf("rotate right...\n");
  dgenrotate_1d(z, 10, -3.5);
  for (i = 0; i<10; i++){
    printf("z[%ld] = %f\n", i, z[i]);
  }
  free(z);
  return 0;
}
