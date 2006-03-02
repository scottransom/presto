#include "ransomfft.h"

#if defined USEFFTW

void fftwcall(fcomplex * indata, long nn, int isign)
/* This routine calls the FFTW complex-complex FFT using stored wisdom */
/* files.  It is VERY fast.  nn does _not_ have to be a power of two   */
/* size.  indata is a complex array but stored as floats.              */
{
   FILE *wisdomfile;
   fftwf_plan *plan_forward, *plan_inverse;
   int ii;
   unsigned long oldestplan = 0;
   static fftwf_plan plancache_forward[4] = { NULL, NULL, NULL, NULL };
   static fftwf_plan plancache_inverse[4] = { NULL, NULL, NULL, NULL };
   static fcomplex *indatacache[4] = { NULL, NULL, NULL, NULL };
   static int firsttime = 1, nncache[4] = { 0, 0, 0, 0 };
   static unsigned long lastuse[4] = { 0, 0, 0, 0 };
   static char wisdomfilenm[120];

   /* Call the six-step algorithm if the FFT is too big to     */
   /* be efficiently handled by FFTW.                          */

   if (nn > BIGFFTWSIZE) {
      tablesixstepfft(indata, nn, isign);
      return;
   }

   /* If calling for the first time, read the wisdom file */

   if (firsttime) {
      /* First try to import the system wisdom if available */
      fftwf_import_system_wisdom();
      sprintf(wisdomfilenm, "%s/fftw_wisdom.txt", DATABASE);
      wisdomfile = fopen(wisdomfilenm, "r");
      if (wisdomfile == NULL) {
         printf("Error opening '%s'.  Run makewisdom again.\n", wisdomfilenm);
         printf("Exiting.\n");
         exit(1);
      }
      if (!fftwf_import_wisdom_from_file(wisdomfile)) {
         printf("Error importing FFTW wisdom.\n");
         printf("Exiting.\n");
         exit(1);
      }
      fclose(wisdomfile);
   }

   /* If we used the same plan during the last few calls, use it again */
   /* We keep, in effect, a stack of the 4 most recent plans.          */

   if (nn == nncache[0] && indata == indatacache[0]) {
      plan_forward = &plancache_forward[0];
      plan_inverse = &plancache_inverse[0];
      lastuse[0] = 0;
      lastuse[1]++;
      lastuse[2]++;
      lastuse[3]++;
      /* printf("Found old plan.  %d  %d\n", nn, ct++); */
   } else if (nn == nncache[1] && indata == indatacache[1]) {
      plan_forward = &plancache_forward[1];
      plan_inverse = &plancache_inverse[1];
      lastuse[1] = 0;
      lastuse[0]++;
      lastuse[2]++;
      lastuse[3]++;
      /* printf("Found old plan.  %d  %d\n", nn, ct++); */
   } else if (nn == nncache[2] && indata == indatacache[2]) {
      plan_forward = &plancache_forward[2];
      plan_inverse = &plancache_inverse[2];
      lastuse[2] = 0;
      lastuse[0]++;
      lastuse[1]++;
      lastuse[3]++;
      /* printf("Found old plan.  %d  %d\n", nn, ct++); */
   } else if (nn == nncache[3] && indata == indatacache[3]) {
      plan_forward = &plancache_forward[3];
      plan_inverse = &plancache_inverse[3];
      lastuse[3] = 0;
      lastuse[0]++;
      lastuse[1]++;
      lastuse[2]++;
      /* printf("Found old plan.  %d  %d\n", nn, ct++); */
   } else {
      if (!firsttime) {
         for (ii = 3; ii >= 0; ii--)
            if (lastuse[ii] >= oldestplan)
               oldestplan = ii;
         if (plancache_forward[oldestplan])
            fftwf_destroy_plan(plancache_forward[oldestplan]);
         if (plancache_inverse[oldestplan])
            fftwf_destroy_plan(plancache_inverse[oldestplan]);
      }
      /* printf("Dammit.  Making a new plan.  %d  %d\n", nn, ct++); */
      plancache_forward[oldestplan] = fftwf_plan_dft_1d(nn,
                                                        (fftwf_complex *) indata,
                                                        (fftwf_complex *) indata,
                                                        -1, FFTW_ESTIMATE);
      plancache_inverse[oldestplan] = fftwf_plan_dft_1d(nn,
                                                        (fftwf_complex *) indata,
                                                        (fftwf_complex *) indata,
                                                        +1, FFTW_ESTIMATE);
      nncache[oldestplan] = nn;
      indatacache[oldestplan] = indata;
      plan_forward = &plancache_forward[oldestplan];
      plan_inverse = &plancache_inverse[oldestplan];
      lastuse[0]++;
      lastuse[1]++;
      lastuse[2]++;
      lastuse[3]++;
      lastuse[oldestplan] = 0;
   }

   /* Call the transform */

   if (isign == -1) {
      fftwf_execute(*plan_forward);
   } else {
      fftwf_execute(*plan_inverse);
   }

   firsttime = 0;
}


#elif defined USESGIFFT


void sgifftcall(fcomplex * indata, long nn, int isign)
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
