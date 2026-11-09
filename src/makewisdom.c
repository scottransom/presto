#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "meminfo.h"
#include "misc_utils.h"
#include "fftw3.h"
#include "cs_fftplan.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *wisdomfile;
    char *wisdomfilenm = NULL;
    fftwf_plan plan;
    fftwf_complex *inout;
    int ii, fftlen, patient = 0;
    unsigned planflag;
    int max_pow2, max_pow10, num_padlen;
    int padlen[20] = { 192, 288, 384, 540, 768, 1080, 1280, 2100, 4200, 5120,
                       7680, 8232, 10240, 12288, 15360, 16464, 25600, 32805, 65610, 131220
    };

    /* Parse the command line:  an optional '-patient' flag and an optional  */
    /* explicit output path (in either order).                               */
    for (ii = 1; ii < argc; ii++) {
        if (strcmp(argv[ii], "-patient") == 0)
            patient = 1;
        else
            wisdomfilenm = strdup(argv[ii]);
    }

    /* By default use FFTW_MEASURE and skip the two largest sizes in each of  */
    /* the loops below -- this is much faster and covers the sizes the tools  */
    /* actually use most.  '-patient' uses FFTW_PATIENT over the full set of  */
    /* sizes for higher-quality (but slow to generate) wisdom.               */
    if (patient) {
        planflag = FFTW_PATIENT;
        max_pow2 = 65536;       /* powers of 2 up to 2^16          */
        max_pow10 = 100000;     /* powers of 10 up to 1e5          */
        num_padlen = 13;        /* padlen[] indices 0..12          */
    } else {
        planflag = FFTW_MEASURE;
        max_pow2 = 16384;       /* skip 32768, 65536               */
        max_pow10 = 1000;       /* skip 10000, 100000              */
        num_padlen = 11;        /* skip padlen 8232, 10240         */
    }

    /* Generate the wisdom... */

    printf("\nAttempting to read the system wisdom file...\n");
    if (!fftwf_import_system_wisdom())
        printf
            ("  failed.  The file probably does not exist.  Tell your sysadmin.\n\n");
    else
        printf("  succeded.  Good.  We'll use it.\n\n");

    if (patient)
        printf("Creating Wisdom for FFTW (FFTW_PATIENT, all sizes).\n");
    else
        printf("Creating Wisdom for FFTW (FFTW_MEASURE, skipping the largest\n"
               "sizes).  Use '-patient' for the fuller, slower, higher-quality\n"
               "wisdom.\n");
    printf("This may take a while...\n\n");
    printf("Generating plans for FFTs of length:\n");

    fftlen = 2;
    inout = fftwf_malloc(sizeof(fftwf_complex) * BIGFFTWSIZE + 2);
    while (fftlen <= max_pow2) {
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, planflag);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, planflag);
        fftwf_destroy_plan(plan);
        printf("   %d real-to-complex\n", fftlen);
        plan = fftwf_plan_dft_r2c_1d(fftlen, (float *) inout, inout, planflag);
        fftwf_destroy_plan(plan);
        fftlen <<= 1;
    }
    fftwf_free(inout);

    fftlen = 10;

    while (fftlen <= max_pow10) {
        inout = fftwf_malloc(sizeof(fftwf_complex) * fftlen);
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, planflag);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, planflag);
        fftwf_destroy_plan(plan);
        fftlen *= 10;
        fftwf_free(inout);
    }

    for (ii = 0; ii < num_padlen; ii++) {
        fftlen = padlen[ii];
        inout = fftwf_malloc(sizeof(fftwf_complex) * fftlen);
        printf("   %d forward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_FORWARD, planflag);
        fftwf_destroy_plan(plan);
        printf("   %d backward\n", fftlen);
        plan = fftwf_plan_dft_1d(fftlen, inout, inout, FFTW_BACKWARD, planflag);
        fftwf_destroy_plan(plan);
        fftwf_free(inout);
    }

    /* coherentsearch's batched complex-to-real transforms, for its default
     * number of harmonics and block size and every decimation factor. */
    {
        const int nh = CS_DEFAULT_NHARMS, nprof = CS_DEFAULT_BLOCKSIZE;
        fftwf_complex *in = fftwf_malloc(sizeof(fftwf_complex) * (nh + 1) * nprof);
        float *out = fftwf_malloc(sizeof(float) * 2 * nh * nprof);
        for (int k = 1; k <= CS_DEFAULT_MAXDECIM; k++) {
            printf("   %d x %d-point complex-to-real (coherentsearch, k=%d)\n", nprof,
                   2 * (nh / k), k);
            plan = cs_plan_c2r(nh, nprof, k, in, out, planflag);
            fftwf_destroy_plan(plan);
        }
        fftwf_free(in);
        fftwf_free(out);
    }

    /* Write to an explicit path if given, else to the PRESTO data directory
     * ($PRESTO/lib if set, otherwise <prefix>/share/presto) where the tools
     * look for it at runtime. */
    if (wisdomfilenm == NULL)
        wisdomfilenm = presto_data_writepath("fftw_wisdom.txt");

    printf("Exporting wisdom to '%s'\n", wisdomfilenm);

    /* Open wisdom file for writing... */

    wisdomfile = fopen(wisdomfilenm, "w");
    if (wisdomfile == NULL) {
        printf("\nError:  Could not open '%s' for writing.\n", wisdomfilenm);
        printf("        Pass an explicit output path as an argument, or make\n");
        printf("        sure the directory exists and is writable.\n\n");
        free(wisdomfilenm);
        return (1);
    }

    /* Write the wisdom... */

    fftwf_export_wisdom_to_file(wisdomfile);

    /* Cleanup... */

    fclose(wisdomfile);
    free(wisdomfilenm);
    printf("Done.\n\n");

    return (0);

}
