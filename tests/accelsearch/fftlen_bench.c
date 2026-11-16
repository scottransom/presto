/*
 * fftlen_bench: time accelsearch's per-kernel correlation step (complex
 * multiply by a kernel + inverse FFT, as in subharm_fderivs_vol()) for every
 * candidate FFT length, single-threaded, with the same FFTW plans accelsearch
 * makes (FFTW_MEASURE, out-of-place inverse, PRESTO wisdom loaded).
 *
 * It answers two questions for this machine:
 *   1. Subharmonics: their useful output is fixed, so for a needed length n
 *      the best choice is simply the cheapest transform of length >= n.
 *      Compared against next_good_fftlen().
 *   2. The fundamental: for kernel width K, the best length L maximizes the
 *      good points per second, i.e. minimizes (t(L) + t_pow*(L-K)) / (L-K).
 *      Compared against fftlen_from_kernwidth() with create_accelobs()'s
 *      2048 floor.
 *
 * Usage: accel_fftlen_bench [-lo N] [-hi N] [-maxprime P] [-nkern K]
 *                           [-rounds R] [-mintime S] [-need n1,n2,...]
 *                           [-kernwidth k1,k2,...] [-patient] [-quiet]
 *
 * Not installed.  Each length's time is the minimum over R interleaved rounds
 * (to shrug off frequency drift and other processes).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "accel.h"

#define MAXLENS 4096
#define MAXLIST 256

typedef struct {
    int len;
    double t;                   /* best seconds per correlation */
} lenres;

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int is_smooth(int n, int maxprime)
{
    static const int primes[] = { 2, 3, 5, 7, 11, 13 };
    for (int ii = 0; ii < 6 && primes[ii] <= maxprime; ii++)
        while (n % primes[ii] == 0)
            n /= primes[ii];
    return n == 1;
}

static void factor_string(int n, char *out, size_t sz)
{
    static const int primes[] = { 2, 3, 5, 7, 11, 13 };
    out[0] = '\0';
    for (int ii = 0; ii < 6; ii++) {
        int e = 0;
        while (n % primes[ii] == 0) {
            n /= primes[ii];
            e++;
        }
        if (e) {
            char tmp[32];
            if (e > 1)
                snprintf(tmp, sizeof(tmp), "%s%d^%d", out[0] ? "*" : "", primes[ii], e);
            else
                snprintf(tmp, sizeof(tmp), "%s%d", out[0] ? "*" : "", primes[ii]);
            strncat(out, tmp, sz - strlen(out) - 1);
        }
    }
}

static int parse_list(const char *s, int *out)
{
    int n = 0;
    char *copy = strdup(s), *tok, *save = NULL;
    for (tok = strtok_r(copy, ",", &save); tok && n < MAXLIST;
         tok = strtok_r(NULL, ",", &save))
        out[n++] = atoi(tok);
    free(copy);
    return n;
}

static unsigned int rng = 2463534242u;
static float frand(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    return (float) (rng >> 8) / 16777216.0f - 0.5f;
}

typedef struct {
    int len;
    fftwf_plan inv;
    float *pdata, *kerns, *tmpdat, *tmpout;
} lenwork;

static void setup(lenwork *w, int len, int nkern, unsigned flag)
{
    w->len = len;
    w->pdata = fftwf_malloc(sizeof(float) * 2 * len);
    w->kerns = fftwf_malloc(sizeof(float) * 2 * len * (size_t) nkern);
    w->tmpdat = fftwf_malloc(sizeof(float) * 2 * len);
    w->tmpout = fftwf_malloc(sizeof(float) * 2 * len);
    /* The same inverse plan as add_accel_plans() in accel_utils.c */
    w->inv = fftwf_plan_dft_1d(len, (fftwf_complex *) w->tmpdat,
                               (fftwf_complex *) w->tmpout, +1,
                               flag | FFTW_DESTROY_INPUT);
    for (int ii = 0; ii < 2 * len; ii++)
        w->pdata[ii] = frand();
    for (size_t ii = 0; ii < 2 * (size_t) len * nkern; ii++)
        w->kerns[ii] = frand();
}

static void teardown(lenwork *w)
{
    fftwf_destroy_plan(w->inv);
    fftwf_free(w->pdata);
    fftwf_free(w->kerns);
    fftwf_free(w->tmpdat);
    fftwf_free(w->tmpout);
}

static volatile float sink;

/* One pass over all nkern kernels: multiply + inverse FFT, as in */
/* subharm_fderivs_vol().  Returns seconds per correlation.       */
static double run_pass(lenwork *w, int nkern, double mintime)
{
    const int len = w->len;
    long ncorr = 0;
    double t0 = now(), t1;
    do {
        for (int jj = 0; jj < nkern; jj++) {
            const float *fkern = w->kerns + 2 * (size_t) len * jj;
            const float *fpdata = w->pdata;
            float *fdata = w->tmpdat;
#if (defined(__GNUC__) || defined(__GNUG__)) &&         \
    !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC ivdep
#endif
            for (int kk = 0; kk < len * 2; kk += 2) {
                const float dr = fpdata[kk], di = fpdata[kk + 1];
                const float kr = fkern[kk], ki = fkern[kk + 1];
                fdata[kk] = dr * kr + di * ki;
                fdata[kk + 1] = di * kr - dr * ki;
            }
            fftwf_execute_dft(w->inv, (fftwf_complex *) w->tmpdat,
                              (fftwf_complex *) w->tmpout);
            sink = w->tmpout[jj % len];
        }
        ncorr += nkern;
        t1 = now();
    } while (t1 - t0 < mintime);
    return (t1 - t0) / ncorr;
}

/* Seconds per output point of the power loop in subharm_fderivs_vol() */
static double time_power(double mintime)
{
    const int n = 2048;
    float *in = fftwf_malloc(sizeof(float) * 2 * n);
    float *out = fftwf_malloc(sizeof(float) * n);
    for (int ii = 0; ii < 2 * n; ii++)
        in[ii] = frand();
    const float norm = 1.0f / (n * (float) n);
    double best = 1e30;
    for (int rr = 0; rr < 5; rr++) {
        long npts = 0;
        double t0 = now(), t1;
        do {
            for (int rep = 0; rep < 100; rep++) {
#if (defined(__GNUC__) || defined(__GNUG__)) &&         \
    !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC ivdep
#endif
                for (int kk = 0; kk < n; kk++) {
                    const int ind = 2 * kk;
                    out[kk] = (in[ind] * in[ind] + in[ind + 1] * in[ind + 1]) * norm;
                }
                sink = out[rep];
            }
            npts += 100L * n;
            t1 = now();
        } while (t1 - t0 < mintime);
        if ((t1 - t0) / npts < best)
            best = (t1 - t0) / npts;
    }
    fftwf_free(in);
    fftwf_free(out);
    return best;
}

static const lenres *cheapest_at_least(const lenres *res, int nres, int n)
{
    const lenres *best = NULL;
    for (int ii = 0; ii < nres; ii++)
        if (res[ii].len >= n && (!best || res[ii].t < best->t))
            best = &res[ii];
    return best;
}

static const lenres *find_len(const lenres *res, int nres, int len)
{
    for (int ii = 0; ii < nres; ii++)
        if (res[ii].len == len)
            return &res[ii];
    return NULL;
}

int main(int argc, char *argv[])
{
    int lo = 128, hi = 4096, maxprime = 7, nkern = 128, rounds = 5, quiet = 0;
    double mintime = 0.02;
    unsigned flag = FFTW_MEASURE;
    int needs[MAXLIST], nneeds = 0, kws[MAXLIST], nkws = 0;

    for (int ii = 1; ii < argc; ii++) {
        if (!strcmp(argv[ii], "-lo") && ii + 1 < argc)
            lo = atoi(argv[++ii]);
        else if (!strcmp(argv[ii], "-hi") && ii + 1 < argc)
            hi = atoi(argv[++ii]);
        else if (!strcmp(argv[ii], "-maxprime") && ii + 1 < argc)
            maxprime = atoi(argv[++ii]);
        else if (!strcmp(argv[ii], "-nkern") && ii + 1 < argc)
            nkern = atoi(argv[++ii]);
        else if (!strcmp(argv[ii], "-rounds") && ii + 1 < argc)
            rounds = atoi(argv[++ii]);
        else if (!strcmp(argv[ii], "-mintime") && ii + 1 < argc)
            mintime = atof(argv[++ii]);
        else if (!strcmp(argv[ii], "-need") && ii + 1 < argc)
            nneeds = parse_list(argv[++ii], needs);
        else if (!strcmp(argv[ii], "-kernwidth") && ii + 1 < argc)
            nkws = parse_list(argv[++ii], kws);
        else if (!strcmp(argv[ii], "-patient"))
            flag = FFTW_PATIENT;
        else if (!strcmp(argv[ii], "-quiet"))
            quiet = 1;
        else {
            fprintf(stderr, "Usage: %s [-lo N] [-hi N] [-maxprime P] [-nkern K] "
                    "[-rounds R] [-mintime S]\n"
                    "         [-need n1,n2,...] [-kernwidth k1,k2,...] [-patient] [-quiet]\n",
                    argv[0]);
            return 1;
        }
    }

    read_wisdom();

    /* Candidate lengths: even (ACCEL_NUMBETWEEN = 2) and maxprime-smooth */
    static lenwork work[MAXLENS];
    static lenres res[MAXLENS];
    int nres = 0;
    for (int len = lo + (lo & 1); len <= hi && nres < MAXLENS; len += 2)
        if (is_smooth(len, maxprime))
            res[nres++].len = len;
    printf("# %d candidate lengths %d..%d (%d-smooth, even), %d kernels each, "
           "%d rounds, %s\n", nres, lo, hi, maxprime, nkern, rounds,
           flag == FFTW_PATIENT ? "FFTW_PATIENT" : "FFTW_MEASURE");
    fflush(stdout);
    for (int ii = 0; ii < nres; ii++) {
        setup(&work[ii], res[ii].len, nkern, flag);
        res[ii].t = 1e30;
    }

    /* Interleaved rounds: every length once per round, keep the minimum */
    for (int rr = 0; rr < rounds; rr++) {
        for (int ii = 0; ii < nres; ii++) {
            double t = run_pass(&work[ii], nkern, mintime);
            if (t < res[ii].t)
                res[ii].t = t;
        }
        if (!quiet) {
            fprintf(stderr, "  round %d of %d done\n", rr + 1, rounds);
        }
    }
    for (int ii = 0; ii < nres; ii++)
        teardown(&work[ii]);
    const double tpow = time_power(mintime);

    /* 1. The raw table */
    printf("#\n# len  factors          ns/corr   ns/point  ns/(N log2 N)  "
           "cheapest>=len\n");
    for (int ii = 0; ii < nres; ii++) {
        char fs[64];
        factor_string(res[ii].len, fs, sizeof(fs));
        const lenres *c = cheapest_at_least(res, nres, res[ii].len);
        printf("%5d  %-15s %9.1f %10.3f %13.4f  %s\n", res[ii].len, fs,
               1e9 * res[ii].t, 1e9 * res[ii].t / res[ii].len,
               1e9 * res[ii].t / (res[ii].len * log2(res[ii].len)),
               c == &res[ii] ? "*" : "");
    }
    printf("# power loop: %.3f ns per output point\n", 1e9 * tpow);

    /* 2. The useful set for subharmonics: lengths cheaper than every longer one */
    printf("#\n# Useful subharmonic lengths (the lower envelope):\n#");
    int nuse = 0;
    for (int ii = 0; ii < nres; ii++)
        if (cheapest_at_least(res, nres, res[ii].len) == &res[ii]) {
            printf(" %d", res[ii].len);
            if (++nuse % 16 == 0)
                printf("\n#");
        }
    printf("\n");

    /* 3. Subharmonic needs: today's length vs. the cheapest */
    if (nneeds) {
        double tnow = 0.0, tbest = 0.0;
        printf("#\n# need   now  ns/corr   best  ns/corr  saving\n");
        for (int ii = 0; ii < nneeds; ii++) {
            const int nowlen = next_good_fftlen(needs[ii]);
            const lenres *n = find_len(res, nres, nowlen);
            const lenres *b = cheapest_at_least(res, nres, needs[ii]);
            if (!n || !b) {
                printf("%6d  (outside the benchmarked range)\n", needs[ii]);
                continue;
            }
            printf("%6d %5d %8.1f  %5d %8.1f  %5.1f%%\n", needs[ii], nowlen,
                   1e9 * n->t, b->len, 1e9 * b->t, 100.0 * (1.0 - b->t / n->t));
            tnow += n->t;
            tbest += b->t;
        }
        if (tnow > 0.0)
            printf("# unweighted total: %.1f%% less correlation time\n",
                   100.0 * (1.0 - tbest / tnow));
    }

    /* 4. Fundamental: maximize good points per second */
    if (nkws) {
        printf("#\n# kernwidth   now  ns/goodpt   best  ns/goodpt  saving\n");
        for (int ii = 0; ii < nkws; ii++) {
            const int K = kws[ii];
            int nowlen = fftlen_from_kernwidth(K);
            if (nowlen < 2048)
                nowlen = 2048;
            const lenres *n = find_len(res, nres, nowlen);
            const lenres *b = NULL;
            double bcost = 1e30, ncost = 1e30;
            for (int jj = 0; jj < nres; jj++) {
                const int good = (res[jj].len - K) / ACCEL_RDR * ACCEL_RDR;
                if (good < 2 * ACCEL_RDR)
                    continue;
                const double cost = (res[jj].t + tpow * good) / good;
                if (cost < bcost) {
                    bcost = cost;
                    b = &res[jj];
                }
                if (&res[jj] == n)
                    ncost = cost;
            }
            if (!n || !b) {
                printf("%10d  (outside the benchmarked range)\n", K);
                continue;
            }
            printf("%10d %5d %10.3f  %5d %10.3f  %5.1f%%\n", K, nowlen, 1e9 * ncost,
                   b->len, 1e9 * bcost, 100.0 * (1.0 - bcost / ncost));
        }
    }
    return 0;
}
