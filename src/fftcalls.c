#include "ransomfft.h"
#include "stdint.h"
#include <errno.h>

// Following gives the same as FFTW's fftwf_alignment_of when
// BYTE_COUNT = 16, which is what we need for SSE.
// 0 means that it is aligned on BYTE_COUNT boundaries
#define is_aligned(POINTER, BYTE_COUNT) \
    ((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT)

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
    // The following resets errno if one of the wisdom files was not found
    errno = 0;
}


void fftwcallsimple(fcomplex * data, long nn, int isign)
/* Simple FFTW calling function for testing */
{
    static int firsttime = 1;
    fftwf_plan plan;
    if (firsttime) {
        read_wisdom();
        firsttime = 0;
    }
    // Note: We need to use FFTW_ESTIMATE since other
    // plan-making destroys the input and output arrays
    plan = fftwf_plan_dft_1d(nn, (fftwf_complex *) data,
                             (fftwf_complex *) data, isign, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
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

    // This determines the alignment of the input array.  Allows
    // more flexible calling of FFTW using its plans.
    // A return value of 0 is "properly" aligned.
    indata_align = is_aligned(indata, 16);
    //if (indata_align)
    //    printf("Data not properly aligned (%d)!\n", indata_align);

    // Call the six-step algorithm if the FFT is too big to be
    // efficiently handled by FFTW
    if (nn > BIGFFTWSIZE) {
        tablesixstepfft(indata, nn, isign);
        return;
    }
    // If calling for the first time, read the wisdom file
    if (firsttime)
        read_wisdom();

    // If we used the same plan during the last few calls, use it
    // again.  We keep, in effect, a stack of the 4 most recent plans.
    ii = 0;
    slot = lastslot;
    while (ii < 4) {
        if (nn == nncache[slot] && indata_align == aligncache[slot]) {
            plan_forward = &plancache_forward[slot];
            plan_inverse = &plancache_inverse[slot];
            lastused[slot] = 0;
            lastused[(slot + 1) % 4]++;
            lastused[(slot + 2) % 4]++;
            lastused[(slot + 3) % 4]++;
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
        fcomplex *datacopy;
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
        //printf("Making a new plan for nn=%ld, align=%d (dropping nn=%ld) %d\n",
        //       nn, indata_align, nncache[oldestplan], badct++);
        // We don't want to wait around to measure huge transforms
        planflag = (nn > 16384) ? FFTW_ESTIMATE : FFTW_MEASURE;
        // FFTW_MEASURE will destroy the input/output data, so copy it
        datacopy = gen_cvect(nn);
        memcpy(datacopy, dataptr, nn * sizeof(fcomplex));
        // Actually make the plans
        plancache_forward[oldestplan] =
            fftwf_plan_dft_1d(nn, dataptr, dataptr, -1, planflag);
        plancache_inverse[oldestplan] =
            fftwf_plan_dft_1d(nn, dataptr, dataptr, +1, planflag);
        // Now copy the input data back
        memcpy(dataptr, datacopy, nn * sizeof(fcomplex));
        vect_free(datacopy);
        nncache[oldestplan] = nn;
        aligncache[oldestplan] = indata_align;
        plan_forward = &plancache_forward[oldestplan];
        plan_inverse = &plancache_inverse[oldestplan];
        lastused[oldestplan] = 0;
        lastused[(oldestplan + 1) % 4]++;
        lastused[(oldestplan + 2) % 4]++;
        lastused[(oldestplan + 3) % 4]++;
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
