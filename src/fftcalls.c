#include "ransomfft.h"
#include "misc_utils.h"
#include "stdint.h"
#include "string.h"
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>

// Following gives the same as FFTW's fftwf_alignment_of when
// BYTE_COUNT = 16, which is what we need for SSE.
// 0 means that it is aligned on BYTE_COUNT boundaries
#define is_aligned(POINTER, BYTE_COUNT) \
    ((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT)

/* What read_wisdom() found, for print_perfinfo() */
static enum { WISDOM_UNREAD, WISDOM_OK, WISDOM_MISSING, WISDOM_REJECTED } wisdom_state;
static char wisdom_path[4096];

void read_wisdom(void)
{
    FILE *wisdomfile;
    char *wisdomfilenm;

    /* First try to import the system wisdom if available */
    fftwf_import_system_wisdom();
    wisdomfilenm = presto_data_path("fftw_wisdom.txt");
    snprintf(wisdom_path, sizeof(wisdom_path), "%s", wisdomfilenm);
    wisdomfile = fopen(wisdomfilenm, "r");
    wisdom_state = WISDOM_OK;
    if (wisdomfile == NULL) {
        wisdom_state = WISDOM_MISSING;
        printf("Warning:  Couldn't open '%s'\n"
               "          You should run 'makewisdom'.  See INSTALL.md.\n",
               wisdomfilenm);
    } else {
        if (!fftwf_import_wisdom_from_file(wisdomfile)) {
            wisdom_state = WISDOM_REJECTED;
            printf("Warning:  FFTW rejected the wisdom in '%s'.\n"
                   "          It was likely written by a different FFTW build or\n"
                   "          version than the one in use now (the wisdom cookie is\n"
                   "          per-build).  Re-run 'makewisdom' with the same FFTW.\n"
                   "          See INSTALL.md.\n",
                   wisdomfilenm);
        }
        fclose(wisdomfile);
    }
    free(wisdomfilenm);
    // The following resets errno if one of the wisdom files was not found
    errno = 0;
}



/* The CPU's model name (Linux), or NULL */
static char *cpu_model_name(char *buf, int len)
{
    FILE *f = fopen("/proc/cpuinfo", "r");
    char line[512];
    char *res = NULL;

    if (f == NULL)
        return NULL;
    while (fgets(line, sizeof(line), f)) {
        char *colon = strchr(line, ':');
        if (colon && strncmp(line, "model name", 10) == 0) {
            colon++;
            while (*colon == ' ' || *colon == '\t')
                colon++;
            colon[strcspn(colon, "\n")] = '\0';
            snprintf(buf, len, "%s", colon);
            res = buf;
            break;
        }
    }
    fclose(f);
    return res;
}


void print_perfinfo(const char *simd, int nthreads)
/* Print the build and runtime facts that decide a tool's speed (for its   */
/* -perfinfo flag): version, compiler and target, CPU, the SIMD kernels    */
/* used (simd, from the tool; NULL if it has no dispatched kernels), its   */
/* OpenMP threads (nthreads; < 0 if built without OpenMP), the             */
/* OPENBLAS_NUM_THREADS setting, and the FFTW version and wisdom.          */
{
    char cpu[256];
    const char *ob = getenv("OPENBLAS_NUM_THREADS");
    const char *target =
#if defined(__AVX512F__)
        "AVX-512";
#elif defined(__AVX2__)
        "AVX2";
#elif defined(__AVX__)
        "AVX";
#elif defined(__x86_64__)
        "baseline x86-64 (SSE2-SSE3)";
#elif defined(__aarch64__) || defined(__ARM_NEON)
        "ARM (NEON)";
#else
        "generic";
#endif

    if (wisdom_state == WISDOM_UNREAD)
        read_wisdom();
    printf("\nPerformance information (-perfinfo):\n");
#ifdef PRESTO_VERSION
    printf("  PRESTO version        %s\n", PRESTO_VERSION);
#endif
#if defined(__clang__)
    printf("  Compiler              clang %s\n", __clang_version__);
#elif defined(__GNUC__)
    printf("  Compiler              gcc %s\n", __VERSION__);
#endif
    printf("  Compiled for          %s%s\n", target,
#if defined(__FMA__)
           ", with FMA"
#else
           ""
#endif
        );
    printf("  CPU                   %s\n",
           cpu_model_name(cpu, sizeof(cpu)) ? cpu : "(unknown)");
    printf("  SIMD kernels          %s\n",
           simd ? simd : "as compiled (no runtime-dispatched kernels in this program)");
    if (nthreads < 0)
        printf("  OpenMP threads        (not built with OpenMP)\n");
    else
        printf("  OpenMP threads        %d (of %ld online CPUs)\n", nthreads,
               sysconf(_SC_NPROCESSORS_ONLN));
    printf("  OPENBLAS_NUM_THREADS  %s\n",
           ob ? ob : "unset (set it to 1 to stop idle OpenBLAS threads spinning)");
    printf("  FFTW                  %s\n", fftwf_version);
    printf("  FFTW wisdom           %s (%s)\n\n", wisdom_path,
           wisdom_state == WISDOM_OK ? "imported" :
           wisdom_state == WISDOM_MISSING ? "MISSING: run makewisdom" :
           "REJECTED: rerun makewisdom with this FFTW");
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
