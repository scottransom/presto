#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "meminfo.h"
#include "fftw3.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(void)
{
    FILE *wisdomfile;
    fftwf_plan plan;
    fftwf_complex *inout;
    int ii, fftlen;
    int padlen[13] = { 288, 540, 1080, 2100, 4200, 8232, 16464, 32805,
        65610, 131220, 262440, 525000, 1050000
    };

    fftlen = 2;

    /* Generate the wisdom... */

    printf("\nAttempting to read the system wisdom file...\n");
    if (!fftwf_import_system_wisdom())
        printf
            ("  failed.  The file probably does not exist.  Tell your sysadmin.\n\n");
    else
        printf("  succeded.  Good.  We'll use it.\n\n");

    printf("Creating Wisdom for FFTW.\n");
    printf("This may take a while...\n\n");
    printf("Generating plans for FFTs of length:\n");

    inout = fftwf_malloc(sizeof(fftwf_complex) * BIGFFTWSIZE + 2);
    while (fftlen <= 1.1e6) {
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        printf("   %d real-to-complex\n", fftlen);
        plan = fftwf_plan_dft_r2c_1d(fftlen, (float *) inout, inout, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        fftlen <<= 1;
    }
    fftwf_free(inout);

    fftlen = 10;

    while (fftlen <= 1.1e6) {
        inout = fftwf_malloc(sizeof(fftwf_complex) * fftlen);
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        fftlen *= 10;
        fftwf_free(inout);
    }

    for (ii = 0; ii < 13; ii++) {
        fftlen = padlen[ii];
        inout = fftwf_malloc(sizeof(fftwf_complex) * fftlen);
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        fftwf_free(inout);
    }

    printf("Exporting wisdom to 'fftw_wisdom.txt'\n");

    /* Open wisdom file for writing... */

    wisdomfile = fopen("fftw_wisdom.txt", "w");

    /* Write the wisdom... */

    fftwf_export_wisdom_to_file(wisdomfile);

    /* Cleanup... */

    fclose(wisdomfile);
    printf("Done.\n\n");

    return (0);

}
