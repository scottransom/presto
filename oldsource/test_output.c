#include "presto.h"

int main(void)
{
  int ii, jj, kk;
  double val[7], err[7];
  char out[7][30];

  for (ii = -20 ; ii < 21 ; ii++){

    val[0] = 0.000 * pow(10.0, ii);
    val[1] = 0.400 * pow(10.0, ii);
    val[2] = 0.490 * pow(10.0, ii);
    val[3] = 0.499 * pow(10.0, ii);
    val[4] = 0.500 * pow(10.0, ii);
    val[5] = 0.590 * pow(10.0, ii);
    val[6] = 0.999 * pow(10.0, ii);
    printf("\n                  %14.5g    %14.5g    %14.5g    %14.5g    %14.5g    %14.5g    %14.5g  \n", val[0], val[1], val[2], val[3], val[4], val[5], val[6]);
    printf("--------------------------------------------------------------------------------------------------------------------------------------------------\n");

    for (jj = -20 ; jj < 21 ; jj++){

      err[0] = 0.000 * pow(10.0, jj);
      err[1] = 0.400 * pow(10.0, jj);
      err[2] = 0.490 * pow(10.0, jj);
      err[3] = 0.499 * pow(10.0, jj);
      err[4] = 0.500 * pow(10.0, jj);
      err[5] = 0.590 * pow(10.0, jj);
      err[6] = 0.999 * pow(10.0, jj);
      
      for (kk = 0 ; kk < 7 ; kk++){

	nice_output_2(out[0], val[0], err[kk], 16);
	nice_output_2(out[1], val[1], err[kk], 16);
	nice_output_2(out[2], val[2], err[kk], 16);
	nice_output_2(out[3], val[3], err[kk], 16);
	nice_output_2(out[4], val[4], err[kk], 16);
	nice_output_2(out[5], val[5], err[kk], 16);
	nice_output_2(out[6], val[6], err[kk], 16);

	printf("  (%10.5g)   %16s  %16s  %16s  %16s  %16s  %16s  %16s\n", err[kk], out[0], out[1], out[2], out[3], out[4], out[5], out[6]);
      }

      printf("\n");
    }

  }
  return 0;
}
