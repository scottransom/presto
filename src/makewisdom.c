#include <stdio.h>
#include <stdlib.h>
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
    int padlen[20] = { 192, 288, 384, 540, 768, 1080, 1280, 2100, 4200, 5120,
                       7680, 8232, 10240, 12288, 15360, 16464, 25600, 32805, 65610, 131220
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
    while (fftlen <= 1.1e5) {
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

    while (fftlen <= 1.1e5) {
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
