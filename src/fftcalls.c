#include "ransomfft.h"

#if defined USEFFTW

#define NUMPLANS 20

void fftwcall(fcomplex *indata, long nn, int isign)
/* This routine calls the FFTW complex-complex FFT using stored wisdom */
/* files.  It is VERY fast.  nn doe _not_ have to be a power of two    */
/* size.  indata is a complex array but stored as floats.              */
{
  int i;
  FILE *wisdomfile;
  fftw_plan plan_forward, plan_inverse;
  static fftw_plan last_plan_forward = NULL, last_plan_inverse = NULL;
  static int firsttime = 1, lastnn = 0;
  static char wisdomfilenm[120];

  /* Call the six-step algorithm if the FFT is too big to */
  /* be efficiently handled by FFTW.                      */

  if (nn > BIGFFTWSIZE){
    tablesixstepfft(indata, nn, isign);
    return;
  }

  /* If calling for the first time, read the wisdom file */

  if (firsttime) {
    sprintf(wisdomfilenm, "%s/fftw_wisdom.txt", DATABASE);
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
      printf("Error opening '%s'.  Run makewisdom again.\n", \
	     wisdomfilenm);
      printf("Exiting.\n");
      exit(1);
    }
    if (FFTW_FAILURE == fftw_import_wisdom_from_file(wisdomfile)) {
      printf("Error importing FFTW wisdom.\n");
      printf("Exiting.\n");
      exit(1);
    }
    fclose(wisdomfile);
  }

  /* If we used the same plan during the last call, use it again */

  if (nn == lastnn){
    plan_forward = last_plan_forward;
    plan_inverse = last_plan_inverse;
  } else {
    if (!firsttime){
      fftw_destroy_plan(last_plan_forward);
      fftw_destroy_plan(last_plan_inverse);
    }
    plan_forward = fftw_create_plan(nn, -1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    plan_inverse = fftw_create_plan(nn, +1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    last_plan_forward = plan_forward;
    last_plan_inverse = plan_inverse;
    lastnn = nn;
  }

  /* Call the transform */

  if (isign == -1){
    fftw_one(plan_forward, (FFTW_COMPLEX *) indata, NULL);
  } else {
    fftw_one(plan_inverse, (FFTW_COMPLEX *) indata, NULL);
  }

  firsttime = 0;
}


#elif defined USESGIFFT


void sgifftcall(fcomplex *indata, long nn, int isign)
{
  int expon;
  double fracpart;
  static complex *coeff[30];

  /* Determine the twoth power of the length of the data */

  fracpart = frexp((double) nn, &expon);
  expon--;

  /* If we are calling using an nn we haven't seen before */

  if (coeff[expon] == NULL) {

    /* Allocate coefficient array */

    coeff[expon] = cfft1di(nn, NULL);

  }
  /* Do the FFT */

  cfft1d(isign, nn, (complex *) indata, 1, coeff[expon]);
}

#endif
