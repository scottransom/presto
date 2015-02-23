#include "ransomfft.h"

#if defined USEFFTW

void read_wisdom(void)
{
    FILE *wisdomfile;
    static char wisdomfilenm[120];

    /* First try to import the system wisdom if available */
    fftwf_import_system_wisdom();
    sprintf(wisdomfilenm, "%s/lib/fftw_wisdom.txt", getenv("PRESTO"));
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
        printf("Warning:  Couldn't open '%s'\n"
               "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n", 
               wisdomfilenm);
    } else {
        if (!fftwf_import_wisdom_from_file(wisdomfile))
            printf("Warning:  '%s' is not up-to-date.\n"
                   "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n", 
                   wisdomfilenm);
        fclose(wisdomfile);
    }
}


void fftwcall(fcomplex * indata, long nn, int isign)
/* This routine calls the FFTW complex-complex FFT using stored wisdom */
/* files.  It is VERY fast.  nn does _not_ have to be a power of two   */
/* size.  indata is a complex array but stored as floats.              */
{
    fftwf_plan *plan_forward, *plan_inverse;
    fftwf_complex *dataptr = (fftwf_complex *) indata;
    int ii, indata_align, slot, incache = 0, oldestplan = 0;
    //static int goodct = 0, badct = 0;
    static fftwf_plan plancache_forward[4] = { NULL, NULL, NULL, NULL };
    static fftwf_plan plancache_inverse[4] = { NULL, NULL, NULL, NULL };
    static int aligncache[4] = { -99, -99, -99, -99 };
    static int firsttime = 1;
    static int lastslot = 0, lastused[4] = { 0, 0, 0, 0 };
    static long nncache[4] = { 0, 0, 0, 0 };

    // Call the six-step algorithm if the FFT is too big to be
    // efficiently handled by FFTW.
    if (nn > BIGFFTWSIZE) {
        tablesixstepfft(indata, nn, isign);
        return;
    }

    // If calling for the first time, read the wisdom file
    if (firsttime) read_wisdom();

    // This determines the alignment of the input array.  Allows
    // more flexible calling of FFTW using its plans.
    // A return value of 0 is "properly" aligned.
    indata_align = fftwf_alignment_of((double *) indata);

    // If we used the same plan during the last few calls, use it
    // again.  We keep, in effect, a stack of the 4 most recent plans.
    ii = 0;
    slot = lastslot;
    while (ii < 4) {
        if (nn == nncache[slot] && indata_align == aligncache[slot]) {
            plan_forward = &plancache_forward[slot];
            plan_inverse = &plancache_inverse[slot];
            lastused[slot] = 0;
            lastused[(slot+1)%4]++;
            lastused[(slot+2)%4]++;
            lastused[(slot+3)%4]++;
            //printf("Found plan in slot %d (iter = %d):  nn=%ld  align=%d  number=%d\n",
            //       slot, ii, nn, aligncache[slot], goodct++);
            lastslot = slot;
            incache = 1;
            break;
        }
        slot = (slot + 1) % 4;
        ii++;
    }
    if (!incache) {
        unsigned int planflag;
        if (!firsttime) {
            for (ii = 3; ii >= 0; ii--)
                if (lastused[ii] >= oldestplan)
                    oldestplan = ii;
            // Delete the old plans to prevent memory leaks
            if (plancache_forward[oldestplan])
                fftwf_destroy_plan(plancache_forward[oldestplan]);
            if (plancache_inverse[oldestplan])
                fftwf_destroy_plan(plancache_inverse[oldestplan]);
        }
        //printf("Making a new plan for nn=%ld (dropping nn=%ld) %d\n",
        //       nn, nncache[oldestplan], badct++);
        // We don't want to wait around to measure huge transforms
        planflag = (nn > 90000) ? FFTW_ESTIMATE : FFTW_MEASURE;
        // Actually make the plans
        plancache_forward[oldestplan] = \
            fftwf_plan_dft_1d(nn, dataptr, dataptr, -1, planflag);
        plancache_inverse[oldestplan] = \
            fftwf_plan_dft_1d(nn, dataptr, dataptr, +1, planflag);
        nncache[oldestplan] = nn;
        aligncache[oldestplan] = indata_align;
        plan_forward = &plancache_forward[oldestplan];
        plan_inverse = &plancache_inverse[oldestplan];
        lastused[oldestplan] = 0;
        lastused[(oldestplan+1)%4]++;
        lastused[(oldestplan+2)%4]++;
        lastused[(oldestplan+3)%4]++;
        lastslot = oldestplan;
    }

    // Call the transform using the "new-array" functionality of FFTW
    if (isign == -1) {
        fftwf_execute_dft(*plan_forward, dataptr, dataptr);
    } else {
        fftwf_execute_dft(*plan_inverse, dataptr, dataptr);
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
