#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "meminfo.h"
#include "sfftw.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(void)
{
  FILE *wisdomfile;
  fftw_plan plan;
  int ii, fftlen;
  int padlen[13] = {288, 540, 1080, 2100, 4200, 8232, 16464, 32805, 
		    65610, 131220, 262440, 525000, 1050000};

  fftlen = 2;

  /* Generate the wisdom... */

  printf("\nCreating Wisdom for FFTW.\n");
  printf("This may take a while...\n\n");
  printf("Generating plans for FFTs of length:\n");

  while (fftlen <= BIGFFTWSIZE) {
    printf("   %d\n", fftlen);
    plan = fftw_create_plan(fftlen, -1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
    plan = fftw_create_plan(fftlen, 1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
    fftlen <<= 1;
  }

  for (ii = 0; ii < 13; ii++){
    fftlen = padlen[ii];
    printf("   %d\n", fftlen);
    plan = fftw_create_plan(fftlen, -1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
    plan = fftw_create_plan(fftlen, 1, \
	FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE);
    fftw_destroy_plan(plan);
  }

  printf("Exporting wisdom to 'fftw_wisdom.txt'\n");

  /* Open wisdom file for writing... */

  wisdomfile = fopen("fftw_wisdom.txt", "w");

  /* Write the wisdom... */

  fftw_export_wisdom_to_file(wisdomfile);

  /* Cleanup... */

  fftw_forget_wisdom();
  fclose(wisdomfile);
  printf("Done.\n\n");

  return (0);

}
