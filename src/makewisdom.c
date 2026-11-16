#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "meminfo.h"
#include "misc_utils.h"
#include "fftw3.h"
#include "cs_fftplan.h"

int next_good_fftlen(int N);   /* corr_prep.c, in libpresto */

/* accelsearch plans its correlation FFTs with FFTW_MEASURE up to this */
/* length (and with FFTW_ESTIMATE, which needs no wisdom, above it).   */
#define ACCEL_MAXWISELEN 16384

/* The plans accelsearch makes for one correlation length: in place in    */
/* both directions (fftwcall(), for the kernels and the forward data FFT, */
/* and add_accel_plans()), and the out-of-place, input-destroying inverse */
/* of add_accel_plans().  Returns how many of the three could be made     */
/* (with FFTW_WISDOM_ONLY in planflag, how many the wisdom covers).       */
static int accel_plans(int fftlen, fftwf_complex *in, fftwf_complex *out,
                       unsigned planflag)
{
    fftwf_plan plan[3];
    int ii, num = 0;

    plan[0] = fftwf_plan_dft_1d(fftlen, in, in, FFTW_FORWARD, planflag);
    plan[1] = fftwf_plan_dft_1d(fftlen, in, in, FFTW_BACKWARD, planflag);
    plan[2] = fftwf_plan_dft_1d(fftlen, in, out, FFTW_BACKWARD,
                                planflag | FFTW_DESTROY_INPUT);
    for (ii = 0; ii < 3; ii++) {
        if (plan[ii]) {
            fftwf_destroy_plan(plan[ii]);
            num++;
        }
    }
    return num;
}

/* makewisdom -check:  does the wisdom (in wisdomfilenm, or where the     */
/* tools read it) cover every FFT accelsearch plans with FFTW_MEASURE?    */
/* Asks FFTW itself, with FFTW_WISDOM_ONLY, so it also catches wisdom     */
/* that FFTW rejects or that was made for another FFT-length list.       */
/* Returns the exit status:  0 if covered, 1 if not, 2 if unreadable.    */
static int check_accel_wisdom(char *wisdomfilenm)
{
    FILE *wisdomfile;
    fftwf_complex *in, *out;
    int fftlen, numlens = 0, nummissing = 0;

    if (wisdomfilenm == NULL)
        wisdomfilenm = presto_data_path("fftw_wisdom.txt");
    printf("FFTW:     %s\n", fftwf_version);
    printf("Wisdom:   %s\n", wisdomfilenm);
    fftwf_import_system_wisdom();
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
        printf("MISSING:  cannot open it.  Run 'makewisdom'.\n");
        return 2;
    }
    if (!fftwf_import_wisdom_from_file(wisdomfile)) {
        fclose(wisdomfile);
        printf("REJECTED: FFTW cannot use it (made by another FFTW build?).\n"
               "          Rerun 'makewisdom' with this FFTW.\n");
        return 2;
    }
    fclose(wisdomfile);
    in = fftwf_malloc(sizeof(fftwf_complex) * ACCEL_MAXWISELEN);
    out = fftwf_malloc(sizeof(fftwf_complex) * ACCEL_MAXWISELEN);
    printf("Missing accelsearch lengths:");
    for (fftlen = next_good_fftlen(1); fftlen <= ACCEL_MAXWISELEN;
         fftlen = next_good_fftlen(fftlen + 1)) {
        numlens++;
        if (accel_plans(fftlen, in, out, FFTW_MEASURE | FFTW_WISDOM_ONLY) < 3) {
            printf(" %d", fftlen);
            nummissing++;
        }
    }
    printf("%s\n", nummissing ? "" : " none");
    fftwf_free(in);
    fftwf_free(out);
    if (nummissing) {
        printf("STALE:    %d of %d accelsearch FFT lengths have no wisdom, so\n"
               "          accelsearch would re-measure them on every run.\n"
               "          Rerun 'makewisdom'.\n", nummissing, numlens);
        return 1;
    }
    printf("OK:       all %d accelsearch FFT lengths (up to %d) are covered.\n",
           numlens, ACCEL_MAXWISELEN);
    return 0;
}

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *wisdomfile;
    char *wisdomfilenm = NULL;
    fftwf_plan plan;
    fftwf_complex *inout;
    int ii, fftlen, patient = 0, check = 0;
    unsigned planflag;
    int max_pow2, max_pow10, num_padlen;
    int padlen[20] = { 192, 288, 384, 540, 768, 1080, 1280, 2100, 4200, 5120,
                       7680, 8232, 10240, 12288, 15360, 16464, 25600, 32805, 65610, 131220
    };

    /* Parse the command line:  an optional '-patient' or '-check' flag and */
    /* an optional explicit wisdom-file path (in either order).             */
    for (ii = 1; ii < argc; ii++) {
        if (strcmp(argv[ii], "-patient") == 0)
            patient = 1;
        else if (strcmp(argv[ii], "-check") == 0)
            check = 1;
        else if (argv[ii][0] == '-') {
            printf("\nUsage:  makewisdom [-patient] [wisdomfile]\n"
                   "        makewisdom -check [wisdomfile]\n\n"
                   "Writes FFTW wisdom (by default to where the PRESTO tools read it).\n"
                   "-check reports whether the existing wisdom covers every FFT that\n"
                   "accelsearch plans with FFTW_MEASURE (exit status 0 if so).\n\n");
            return (argv[ii][1] == 'h' ? 0 : 1);
        } else
            wisdomfilenm = strdup(argv[ii]);
    }
    if (check)
        return check_accel_wisdom(wisdomfilenm);

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
     * and -deep numbers of harmonics at the default block size, and every
     * decimation factor of each.  They
     * are small, so they always get FFTW_PATIENT (under a second in all; it
     * makes them ~2x faster than FFTW_ESTIMATE, and FFTW_EXHAUSTIVE is no
     * better for 20x the planning time). */
    {
        const int nhs[2] = { CS_DEFAULT_NHARMS, CS_DEEP_NHARMS };
        const int maxk[2] = { CS_DEFAULT_MAXDECIM, CS_DEEP_MAXDECIM };
        const int nprof = CS_DEFAULT_BLOCKSIZE;
        for (int ii = 0; ii < 2; ii++) {
            const int nh = nhs[ii];
            fftwf_complex *in = fftwf_malloc(sizeof(fftwf_complex) * (nh + 1) * nprof);
            float *out = fftwf_malloc(sizeof(float) * 2 * nh * nprof);
            for (int k = 1; k <= maxk[ii]; k++) {
                printf("   %d x %d-point complex-to-real (coherentsearch%s, k=%d)\n", nprof,
                       2 * (nh / k), ii ? " -deep" : "", k);
                plan = cs_plan_c2r(nh, nprof, k, in, out, FFTW_PATIENT);
                fftwf_destroy_plan(plan);
            }
            fftwf_free(in);
            fftwf_free(out);
        }
    }

    /* accelsearch's correlation lengths: every length next_good_fftlen()
     * can return, up to the 16384 above which accelsearch plans with
     * FFTW_ESTIMATE.  Without wisdom, each one costs a FFTW_MEASURE planning
     * (repeatedly, as fftwcall() caches only 4 plans) on every run.  Plan the
     * same transforms accelsearch does: in place in both directions (as
     * fftwcall() does, for the kernels and the forward data FFT), and the
     * out-of-place, input-destroying inverse of add_accel_plans(). */
    {
        fftwf_complex *in = fftwf_malloc(sizeof(fftwf_complex) * ACCEL_MAXWISELEN);
        fftwf_complex *out = fftwf_malloc(sizeof(fftwf_complex) * ACCEL_MAXWISELEN);
        for (fftlen = next_good_fftlen(1); fftlen <= ACCEL_MAXWISELEN;
             fftlen = next_good_fftlen(fftlen + 1)) {
            printf("   %d forward, backward, out-of-place backward (accelsearch)\n",
                   fftlen);
            accel_plans(fftlen, in, out, planflag);
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
