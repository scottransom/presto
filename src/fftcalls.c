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
  static fftw_plan plan_forward[NUMPLANS], plan_inverse[NUMPLANS];
  static fftw_plan *good_plan_forward, *good_plan_inverse;
  static fftw_plan est_plan_forward, est_plan_inverse;
  static long fftsizes[NUMPLANS];
  static int firsttime = 1, currentindex = 0;
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

    firsttime = 0;
  }

  /* Search the stack of plans for the proper size transform */

  for(i=currentindex-1; i>=0; i--){
    if (fftsizes[i] == nn){
      good_plan_forward = &plan_forward[i];
      good_plan_inverse = &plan_inverse[i];
      break;
    }
  }
    
  /* If we are calling using an nn we haven't seen before */

  if (i < 0) {

    if (currentindex < NUMPLANS){

      /* Create new plans if space in the array */
      
      plan_forward[currentindex] = fftw_create_plan(nn, -1, FFTW_MEASURE | \
						    FFTW_USE_WISDOM | \
						    FFTW_IN_PLACE);
      plan_inverse[currentindex] = fftw_create_plan(nn, +1, FFTW_MEASURE | \
						    FFTW_USE_WISDOM | \
						    FFTW_IN_PLACE);
      good_plan_forward = &plan_forward[currentindex];
      good_plan_inverse = &plan_inverse[currentindex];
      currentindex++;

    } else {

      /* Estimate a plan if we are out of space */

      est_plan_forward = fftw_create_plan(nn, -1, FFTW_ESTIMATE | \
					  FFTW_USE_WISDOM | \
					  FFTW_IN_PLACE);
      est_plan_inverse = fftw_create_plan(nn, +1, FFTW_ESTIMATE | \
					  FFTW_USE_WISDOM | \
					  FFTW_IN_PLACE);
      good_plan_forward = &est_plan_forward;
      good_plan_inverse = &est_plan_inverse;

    }
  }

  /* Call the transform */

  if (isign == -1){

    fftw(*good_plan_forward, 1, (FFTW_COMPLEX *) indata, 1, 1, \
	 NULL, 1, 1);

  } else {

    fftw(*good_plan_inverse, 1, (FFTW_COMPLEX *) indata, 1, 1, \
	 NULL, 1, 1);

  }
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
