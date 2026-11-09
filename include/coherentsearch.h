/*
 * coherentsearch: a coherent harmonic-summing pulsar search over PRESTO .fft
 * files.  A C port of the production search in CoherentSearch.jl v1.0.0.
 *
 * The algorithm, stage by stage and with every convention spelled out, is
 * tests/coherentsearch/cs_oracle.py; tests/test_coherentsearch.py pins both it
 * and this code to numbers written by the Julia package itself.  In brief, for
 * a chunk of n trial fundamentals r_t = r_lo + t*lodr (global trial index t):
 *
 *   1. interpolate each harmonic h*r_t exactly (Eqn. 30, m bins) into a
 *      harmonic stack, leaving a harmonic's row zero for the whole chunk unless
 *      it fits at both the chunk's first and last trial;
 *   2. inverse-FFT (unnormalised c2r) each stack into a 2*nharms-bin profile,
 *      and every k-th harmonic of it into a 2*(nharms/k)-bin profile for each
 *      decimation factor k;
 *   3. score each profile with the peak boxcar matched filter over a pruned
 *      width bank, normalised by a per-chunk noise scale (analytic or measured)
 *      and a per-width factor for band-limited noise.
 *
 * Precision follows the Julia production defaults: float32 interpolation
 * weights, harmonic stacks, profiles and boxcar scans; the per-width shape
 * factors and the analytic noise scale are computed in double.
 */
#ifndef COHERENTSEARCH_H
#define COHERENTSEARCH_H

#include "presto.h"

#define CS_MAXDECIM 64          /* most decimation factors in one search */
#define CS_SIGMA_SAMPLES 8192   /* bins pooled for the measured noise scale */
#define CS_LADDER_WMAX 6        /* widest boxcar kept on any fold but the shallowest */
#define CS_DIRECT_QMAX 65536    /* largest trial-grid denominator tabulated */

typedef struct {
    int nharms;              /* harmonics coherently summed (default 60) */
    int m;                   /* Fourier bins in the interpolation kernel (even, default 16) */
    double hidr;             /* Fourier-bin step at the highest harmonic (default 0.5) */
    double threshold;        /* candidate S/N cutoff (default 8) */
    double boxcar_fsp;       /* width-bank growth factor (default 1.5) */
    double boxcar_maxfrac;   /* widest boxcar as a fraction of the profile (default 0.3) */
    int ndecim;              /* number of decimation factors, decims[0] == 1 */
    int decims[CS_MAXDECIM];
    int sigma_measured;      /* 0: analytic noise scale (default), 1: measured (MAD) */
} cs_params;

/* The coherent_search CLI defaults, with decimations 1..maxdecim. */
void cs_default_params(cs_params *p, int maxdecim);

/* The trial step hidr/nharms in Fourier bins. */
static inline double cs_lodr(const cs_params *p)
{
    return p->hidr / p->nharms;
}

/* A PRESTO .fft file held in memory. */
typedef struct {
    fcomplex *amps;          /* the packed amplitudes (amps[0] = DC.r + i Nyquist) */
    long long numamps;
    long long N;             /* time-series length, from the .inf */
    double T;                /* observation length (s) */
} cs_fft;

/* Read basename.fft and basename.inf (filenm may end in .fft). */
void cs_read_fft(const char *filenm, cs_fft *fft);
void cs_free_fft(cs_fft *fft);

/*
 * Per-harmonic interpolation plan.  The trial grid is the exact rational
 * r_t = r_lo + t*pnum/q, so harmonic h advances by h*pnum/q bins per trial:
 * an integer part base_adv plus a residue step s (mod q).  Only P of the q
 * residues are ever visited; row[] maps a residue to its row in W/A/carry.
 * P == 0 marks the fallback for a trial step that is not a small rational,
 * which recomputes the weights per trial.
 */
typedef struct {
    int h, m, P;
    long long q, s, base_adv, pnum;
    long long rfloor0;       /* floor(h * r_lo) */
    double dr0;              /* h * r_lo - rfloor0 */
    double dh;               /* h * lodr (fallback path only) */
    int *row;                /* (q) residue -> row, -1 if never visited */
    float *W;                /* (P, m): 1/(dr - j), j = 1-m/2 .. m/2 */
    fcomplex *A;             /* (P): conj(sin(pi dr) exp(i pi dr)) / pi */
    unsigned char *carry;    /* (P): frac(h r_lo) + res/q wrapped past 1 */
} cs_hplan;

/* Build the nharms plans for a search whose global trial 0 is at r_lo. */
cs_hplan *cs_build_plans(const cs_params *p, double r_lo);
void cs_free_plans(cs_hplan *plans, int nharms);

/* Whether harmonic plan hp has a full m-bin window below Nyquist at trial t. */
int cs_harmonic_fits(const cs_hplan *hp, const cs_fft *fft, long long t);

/* Interpolate harmonic hp->h at trials t0..t0+n-1 into out[j*stride]
 * (j = 0..n-1).  Returns 0, writing nothing, if the harmonic does not fit at
 * both ends of the range. */
int cs_fill_row(const cs_hplan *hp, const cs_fft *fft, long long t0, int n,
                fcomplex *out, int stride);

/* --- the metric -------------------------------------------------------- */

/* Geometric boxcar width bank for nbins-bin profiles; returns the count. */
int cs_boxcar_widths(int nbins, double fsp, double maxfrac, int *widths);
/* The bank for decimation k, pruned against the rest of the ladder. */
int cs_ladder_widths(int nbins, int k, const cs_params *p, int *widths);
/* Noise scale of a decimation-k profile from the filled harmonics. */
double cs_analytic_sigma(const unsigned char *filled, int k, int Hk);
/* Per-width sigma / sd(boxcar sum) for band-limited noise. */
void cs_boxcar_shape(const unsigned char *filled, int k, int Hk,
                     const int *widths, int nwidths, double *shape);
/* 1.4826 * median(|x|) pooled over (a subsample of) n nbins-bin profiles.
 * buf must hold CS_SIGMA_SAMPLES floats. */
float cs_block_sigma(const float *profs, int nbins, int n, float *buf);
/* Peak boxcar S/N of one profile.  psum needs nbins + widths[nw-1] + 1. */
float cs_profile_boxcar(const float *prof, int nbins, const int *widths, int nw,
                        float invsigma, const double *shape, float *psum);

/* --- one chunk ------------------------------------------------------------ */

typedef struct {
    int k, Hk, nbins, nwidths;
    int widths[64];
    double shape[64];
    fftwf_plan plan;         /* c2r of the stride-k rows of ftprofs */
    float *profs;            /* (nprof, nbins) */
    float *psum;
    float *mvals;            /* (nprof) metrics of the last scored chunk */
    int nvalid;              /* trials scored (k*r below Nyquist) */
    double sigma_analytic;
    float sigma_measured;
    unsigned char *shape_filled;  /* the filled set shape[] was computed for */
    int shape_valid;
} cs_decim;

typedef struct {
    int nharms, nprof, ndecim;
    fcomplex *ftprofs;       /* (nprof, nharms+1); row 0 (DC) always zero */
    unsigned char *filled;   /* (nharms+1), by harmonic number */
    float *sigbuf;           /* CS_SIGMA_SAMPLES */
    cs_decim *d;             /* ndecim, d[0].k == 1 */
} cs_workspace;

cs_workspace *cs_workspace_new(const cs_params *p, int nprof);
void cs_workspace_free(cs_workspace *ws);

/* Fill ws->ftprofs and ws->filled for trials t0..t0+n-1 (n <= nprof). */
void cs_fill_chunk(cs_workspace *ws, const cs_hplan *plans, const cs_fft *fft,
                   const cs_params *p, long long t0, int n);

/* Inverse-FFT decimation d's profiles for the current chunk, find its valid
 * trials (rstart is the chunk's first fundamental) and both noise scales, and
 * refresh its shape factors. */
void cs_prepare_decim(cs_workspace *ws, cs_decim *d, const cs_fft *fft,
                      const cs_params *p, double rstart, int n);

/* Score decimation d's nvalid profiles into d->mvals with noise scale sigma. */
void cs_score_decim(cs_decim *d, double sigma);

#endif
