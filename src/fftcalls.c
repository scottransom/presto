#include "ransomfft.h"

#if defined USEFFTW

void fftwcall(fcomplex *indata, long nn, int isign)
/* This routine calls the FFTW complex-complex FFT using stored wisdom */
/* files.  It is VERY fast.  nn doe _not_ have to be a power of two    */
/* size.  indata is a complex array but stored as floats.              */
{
  FILE *wisdomfile;
  fftw_plan plan_forward, plan_inverse;
  static fftw_plan plancache_forward[4] = {NULL, NULL, NULL, NULL};
  static fftw_plan plancache_inverse[4] = {NULL, NULL, NULL, NULL};
  static int firsttime = 1, nncache[4] = {0, 0, 0, 0};
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

  /* If we used the same plan during the last few calls, use it again */
  /* We keep, in effect, a stack of the 4 most recent plans.          */

  if (nn == nncache[0]){
printf("oldplan0\n");
    plan_forward = plancache_forward[0];
    plan_inverse = plancache_inverse[0];
  } else if (nn == nncache[1]){
printf("oldplan1\n");
    plan_forward = plancache_forward[1];
    plan_inverse = plancache_inverse[1];
    plancache_forward[1] = plancache_forward[0];
    plancache_inverse[1] = plancache_inverse[0];
    nncache[1] = nncache[0];
    plancache_forward[0] = plan_forward;
    plancache_inverse[0] = plan_inverse;
    nncache[0] = nn;
  } else if (nn == nncache[2]){
printf("oldplan2\n");
    plan_forward = plancache_forward[2];
    plan_inverse = plancache_inverse[2];
    plancache_forward[2] = plancache_forward[1];
    plancache_inverse[2] = plancache_inverse[1];
    nncache[2] = nncache[1];
    plancache_forward[1] = plancache_forward[0];
    plancache_inverse[1] = plancache_inverse[0];
    nncache[1] = nncache[0];
    plancache_forward[0] = plan_forward;
    plancache_inverse[0] = plan_inverse;
    nncache[0] = nn;
  } else if (nn == nncache[3]){
printf("oldplan3\n");
    plan_forward = plancache_forward[3];
    plan_inverse = plancache_inverse[3];
    plancache_forward[3] = plancache_forward[2];
    plancache_inverse[3] = plancache_inverse[2];
    nncache[3] = nncache[2];
    plancache_forward[2] = plancache_forward[1];
    plancache_inverse[2] = plancache_inverse[1];
    nncache[2] = nncache[1];
    plancache_forward[1] = plancache_forward[0];
    plancache_inverse[1] = plancache_inverse[0];
    nncache[1] = nncache[0];
    plancache_forward[0] = plan_forward;
    plancache_inverse[0] = plan_inverse;
    nncache[0] = nn;
  } else {
printf("*** newplan ***\n");
    if (!firsttime){
      plancache_forward[3] = plancache_forward[2];
      plancache_inverse[3] = plancache_inverse[2];
      nncache[3] = nncache[2];
      plancache_forward[2] = plancache_forward[1];
      plancache_inverse[2] = plancache_inverse[1];
      nncache[2] = nncache[1];
      plancache_forward[1] = plancache_forward[0];
      plancache_inverse[1] = plancache_inverse[0];
      nncache[1] = nncache[0];
      fftw_destroy_plan(plancache_forward[0]);
      fftw_destroy_plan(plancache_inverse[0]);
    }
    plancache_forward[0] = fftw_create_plan(nn, -1, FFTW_ESTIMATE | \
					    FFTW_USE_WISDOM | \
					    FFTW_IN_PLACE);
    plancache_inverse[0] = fftw_create_plan(nn, +1, FFTW_ESTIMATE | \
					    FFTW_USE_WISDOM | \
					    FFTW_IN_PLACE);
    nncache[0] = nn;
    plan_forward = plancache_forward[0];
    plan_inverse = plancache_inverse[0];
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
