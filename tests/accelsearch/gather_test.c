/*
 * gather_test: check every version of accelsearch's add_gathered_row() that
 * this CPU can run (scalar, AVX2, AVX-512, NEON and the dispatched one)
 * against the scalar loop, bit for bit, on the r index tables accelsearch
 * builds (ffdot_geometry() in accel_utils.c) for every harmonic of
 * numharm = 2 ... 32.  Each input row is allocated with exactly
 * ACCEL_ROWPAD floats after the last one used, so a sanitizer build also
 * checks that the SIMD versions stay inside that padding.
 *
 * Run by tests/test_accel_gather.py.  Prints a line per kernel; exits 1 on
 * any mismatch.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "accel.h"

#define MAXKERN 8

static unsigned int rng_state = 12345u;
static float rand_float(void)
{
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float) (rng_state >> 8) / 16777216.0f * 50.0f;
}

/* The fundamental r bins fullrlo + k * ACCEL_DR (k < n) and the matching */
/* bins of harmonic harmnum/numharm, as ffdot_geometry() computes them.   */
/* Returns the subharmonic row length (numrs).                            */
static int make_rinds(double fullrlo, int n, int numharm, int harmnum, int *rinds)
{
    const double harm_fract = (double) harmnum / (double) numharm;
    const double fullrhi = fullrlo + (n - 1) * ACCEL_DR;
    const double drlo = rint(ACCEL_RDR * fullrlo * harm_fract) * ACCEL_DR;
    const double drhi = rint(ACCEL_RDR * fullrhi * harm_fract) * ACCEL_DR;
    const long long rlo = (long long) floor(drlo);

    for (int ii = 0; ii < n; ii++) {
        const double rr = fullrlo + ii * ACCEL_DR;
        const double subr = rint(ACCEL_RDR * rr * harm_fract) * ACCEL_DR;
        rinds[ii] = (int) ((subr - rlo) * ACCEL_RDR + DBLCORRECT);
    }
    int numrs = (int) ((ceil(drhi) - floor(drlo)) * ACCEL_RDR + DBLCORRECT) + 1;
    if (numrs % ACCEL_RDR)
        numrs = (numrs / ACCEL_RDR + 1) * ACCEL_RDR;
    return numrs;
}

int main(void)
{
    const char *names[MAXKERN];
    accel_gather_fn fns[MAXKERN];
    const int numkern = accel_gather_kernels(names, fns, MAXKERN);
    const double rlos[] = { 1000.0, 1000.5, 12345.0, 777.25, 3.3, 98765.75 };
    const int lens[] = { 1544, 1804, 1, 7, 8, 15, 16, 17, 31, 33, 100 };
    long ntests[MAXKERN] = { 0 }, nfail[MAXKERN] = { 0 };
    int bad = 0;

    printf("add_gathered_row() in this run: %s\n", accel_simd_name());
    for (int numharm = 2; numharm <= 32; numharm *= 2) {
        for (int harmnum = 1; harmnum < numharm; harmnum += 2) {
            for (size_t ir = 0; ir < sizeof(rlos) / sizeof(rlos[0]); ir++) {
                for (size_t il = 0; il < sizeof(lens) / sizeof(lens[0]); il++) {
                    const int n = lens[il];
                    int *rinds = (int *) malloc(n * sizeof(int));
                    const int numrs = make_rinds(rlos[ir], n, numharm, harmnum, rinds);

                    /* The property the SIMD versions rely on */
                    for (int kk = 0; kk < n; kk++) {
                        const int step = kk ? rinds[kk] - rinds[kk - 1] : 0;
                        if (rinds[kk] < 0 || rinds[kk] >= numrs || step < 0 || step > 1) {
                            printf("bad r index table: numharm %d harmnum %d rlo %g "
                                   "k %d: rinds %d, numrs %d\n", numharm, harmnum,
                                   rlos[ir], kk, rinds[kk], numrs);
                            bad = 1;
                            break;
                        }
                    }
                    float *in = (float *) malloc((numrs + ACCEL_ROWPAD) * sizeof(float));
                    float *out0 = (float *) malloc(n * sizeof(float));
                    float *ref = (float *) malloc(n * sizeof(float));
                    float *out = (float *) malloc(n * sizeof(float));
                    for (int kk = 0; kk < numrs + ACCEL_ROWPAD; kk++)
                        in[kk] = rand_float();
                    for (int kk = 0; kk < n; kk++)
                        out0[kk] = rand_float();
                    memcpy(ref, out0, n * sizeof(float));
                    fns[0](ref, in, rinds, n);
                    for (int kern = 0; kern < numkern; kern++) {
                        memcpy(out, out0, n * sizeof(float));
                        fns[kern](out, in, rinds, n);
                        ntests[kern]++;
                        if (memcmp(out, ref, n * sizeof(float))) {
                            if (!nfail[kern])
                                printf("%s differs: numharm %d harmnum %d rlo %g n %d\n",
                                       names[kern], numharm, harmnum, rlos[ir], n);
                            nfail[kern]++;
                            bad = 1;
                        }
                    }
                    free(in);
                    free(out0);
                    free(ref);
                    free(out);
                    free(rinds);
                }
            }
        }
    }
    for (int kern = 0; kern < numkern; kern++)
        printf("%-10s %ld cases, %ld differ\n", names[kern], ntests[kern], nfail[kern]);
    return bad;
}
