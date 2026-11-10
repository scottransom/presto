/*
 * coherentsearch: reading .fft files and exact Fourier interpolation of the
 * harmonic stacks (CoherentSearch.jl's directinterp.jl).
 *
 * Harmonic h of trial t sits at h*r_lo + h*t*pnum/q Fourier bins.  The integer
 * and fractional parts are tracked with exact integer arithmetic, so the m
 * weights 1/(dr - j) and the scale conj(sin(pi dr) e^{i pi dr})/pi depend only
 * on the residue (h*t*pnum) mod q, and are tabulated once per harmonic.
 */
#include "coherentsearch.h"

void cs_default_params(cs_params *p, int maxdecim)
{
    p->nharms = 60;
    p->m = 16;
    p->hidr = 0.5;
    p->threshold = 8.0;
    p->boxcar_fsp = 1.5;
    p->boxcar_maxfrac = 0.3;
    p->sigma_measured = 0;
    cs_set_decimations(p, maxdecim);
}

void cs_set_decimations(cs_params *p, int maxdecim)
{
    /* Keep only the factors that still leave >= 2 harmonics to fold. */
    p->ndecim = 0;
    for (int k = 1; k <= (maxdecim < 1 ? 1 : maxdecim) && p->ndecim < CS_MAXDECIM; k++)
        if (p->nharms / k >= 2)
            p->decims[p->ndecim++] = k;
}

void cs_read_fft(const char *filenm, cs_fft *fft)
{
    cs_read_fft_info(filenm, fft, 1);
}

void cs_read_fft_info(const char *filenm, cs_fft *fft, int read_amps)
{
    char *base = calloc(strlen(filenm) + 8, 1);
    strcpy(base, filenm);
    size_t len = strlen(base);
    if (len > 4 && strcmp(base + len - 4, ".fft") == 0)
        base[len - 4] = '\0';
    infodata idata;
    memset(&idata, 0, sizeof(idata));
    readinf(&idata, base);
    char *fftnm = calloc(strlen(base) + 8, 1);
    sprintf(fftnm, "%s.fft", base);
    fft->amps = NULL;
    fft->numamps = 0;
    if (read_amps) {
        FILE *f = chkfopen(fftnm, "rb");
        fft->numamps = chkfilelen(f, sizeof(fcomplex));
        fft->amps = gen_cvect(fft->numamps);
        chkfread(fft->amps, sizeof(fcomplex), fft->numamps, f);
        fclose(f);
    }
    fft->N = (long long) idata.N;
    fft->T = idata.N * idata.dt;
    fft->dm = idata.dm;
    fft->idata = idata;
    snprintf(fft->object, sizeof(fft->object), "%s", idata.object);
    free(fftnm);
    free(base);
}

void cs_free_fft(cs_fft *fft)
{
    vect_free(fft->amps);
    fft->amps = NULL;
}

/* The trial step as an exact small rational pnum/q (Julia's
 * rationalize(Int, lodr; tol=1e-13) followed by the same sanity checks), via
 * continued-fraction convergents.  Returns 0 if there is none. */
static int trial_grid_rational(double lodr, long long *pnum, long long *q)
{
    if (!(lodr > 0.0))
        return 0;
    long long p0 = 0, q0 = 1, p1 = 1, q1 = 0;   /* convergents h_{-2}/k_{-2}, h_{-1}/k_{-1} */
    double x = lodr;
    for (int it = 0; it < 64; it++) {
        double a = floor(x);
        if (a > 1e15)
            break;
        long long ai = (long long) a;
        long long p2 = ai * p1 + p0, q2 = ai * q1 + q0;
        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;
        if (fabs((double) p1 / (double) q1 - lodr) <= 1e-13)
            break;
        double frac = x - a;
        if (frac <= 0.0)
            break;
        x = 1.0 / frac;
    }
    if (q1 <= 0 || q1 > CS_DIRECT_QMAX || p1 < 1)
        return 0;
    if (fabs((double) p1 / (double) q1 - lodr) > 8 * (nextafter(lodr, INFINITY) - lodr))
        return 0;
    *pnum = p1;
    *q = q1;
    return 1;
}

static long long gcdll(long long a, long long b)
{
    while (b) {
        long long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

/* Tabulate hp's trials-axis group blocks (CoherentSearch.jl's
 * build_group_table).  Group g starts at residue g*V*s mod q; within it the
 * k-th trial's window starts delta_k = b_k - b_0 bins after the first one
 * (0/1 steps plus base_adv per trial), so lane k's m weights sit at rows
 * delta_k .. delta_k+m-1 of the block and it reads one contiguous slice of
 * the bin planes. */
static void build_group_table(cs_hplan *hp)
{
    enum { V = CS_GROUP_V };
    const int m = hp->m;
    const long long q = hp->q, s = hp->s;
    const long long gstep = (V * s) % q;
    hp->ngrp = gstep == 0 ? 1 : (int) (q / gcdll(gstep, q));
    hp->grow = malloc(q * sizeof(int));
    for (long long i = 0; i < q; i++)
        hp->grow[i] = -1;
    hp->goff = malloc((hp->ngrp + 1) * sizeof(long));
    hp->gnj = malloc(hp->ngrp * sizeof(int));
    hp->gA = malloc((size_t) hp->ngrp * V * sizeof(fcomplex));
    int ps[V], ds[V];
    /* First pass: sizes. */
    long total = 0;
    long long res0 = 0;
    for (int g = 0; g < hp->ngrp; g++) {
        long long res = res0, B = 0, B0 = 0;
        for (int k = 0; k < V; k++) {
            const int pp = hp->row[res];
            const long long Bk = B + hp->carry[pp];
            if (k == 0)
                B0 = Bk;
            ps[k] = pp;
            ds[k] = (int) (Bk - B0);
            hp->gA[(size_t) g * V + k] = hp->A[pp];
            res += s;
            B += hp->base_adv;
            if (res >= q) {
                res -= q;
                B++;
            }
        }
        hp->grow[res0] = g;
        hp->gnj[g] = m + ds[V - 1];
        hp->goff[g] = total;
        total += (long) hp->gnj[g] * V;
        res0 = (res0 + gstep) % q;
    }
    hp->goff[hp->ngrp] = total;
    hp->gW = calloc(total, sizeof(float));
    /* Second pass: fill the blocks. */
    res0 = 0;
    for (int g = 0; g < hp->ngrp; g++) {
        long long res = res0, B = 0, B0 = 0;
        for (int k = 0; k < V; k++) {
            const int pp = hp->row[res];
            const long long Bk = B + hp->carry[pp];
            if (k == 0)
                B0 = Bk;
            const int dk = (int) (Bk - B0);
            float *blk = hp->gW + hp->goff[g];
            for (int i = 0; i < m; i++)
                blk[(long) (dk + i) * V + k] = hp->W[(size_t) pp * m + i];
            res += s;
            B += hp->base_adv;
            if (res >= q) {
                res -= q;
                B++;
            }
        }
        res0 = (res0 + gstep) % q;
    }
}

cs_hplan *cs_build_plans(const cs_params *p, double r_lo)
{
    const int m = p->m, m2 = p->m / 2;
    const double lodr = cs_lodr(p);
    long long pnum = 0, q = 0;
    int rational = trial_grid_rational(lodr, &pnum, &q);
    cs_hplan *plans = calloc(p->nharms, sizeof(cs_hplan));
    for (int h = 1; h <= p->nharms; h++) {
        cs_hplan *hp = plans + h - 1;
        double hr = h * r_lo;
        hp->h = h;
        hp->m = m;
        hp->rfloor0 = (long long) floor(hr);
        hp->dr0 = hr - (double) hp->rfloor0;
        hp->dh = h * lodr;
        if (!rational) {
            hp->P = 0;
            hp->q = 1;
            hp->pnum = 1;
            continue;
        }
        long long step = h * pnum;
        hp->q = q;
        hp->pnum = pnum;
        hp->s = step % q;
        hp->base_adv = step / q;
        hp->P = (int) (q / gcdll(step, q));
        hp->row = malloc(q * sizeof(int));
        for (long long i = 0; i < q; i++)
            hp->row[i] = -1;
        hp->W = malloc((size_t) hp->P * m * sizeof(float));
        hp->A = malloc((size_t) hp->P * sizeof(fcomplex));
        hp->carry = malloc((size_t) hp->P);
        long long res = 0;
        for (int pp = 0; pp < hp->P; pp++) {
            hp->row[res] = pp;
            /* Computed in double and rounded once on storage. */
            double u = hp->dr0 + (double) res / (double) q;
            hp->carry[pp] = (u >= 1.0);
            double dr = hp->carry[pp] ? u - 1.0 : u;
            float *w = hp->W + (size_t) pp * m;
            if (dr == 0.0) {
                /* r is exactly a Fourier bin: the kernel is a delta on it. */
                hp->A[pp].r = 1.0f;
                hp->A[pp].i = 0.0f;
                for (int i = 0; i < m; i++)
                    w[i] = (i == m2 - 1) ? 1.0f : 0.0f;
            } else {
                /* conj(sin(pi dr) * exp(i pi dr) / pi) */
                double s = sin(M_PI * dr) / M_PI;
                hp->A[pp].r = (float) (s * cos(M_PI * dr));
                hp->A[pp].i = (float) (-s * sin(M_PI * dr));
                for (int i = 0; i < m; i++)
                    w[i] = (float) (1.0 / (dr - (i + 1 - m2)));
            }
            res += hp->s;
            if (res >= q)
                res -= q;
        }
        build_group_table(hp);
    }
    return plans;
}

void cs_free_plans(cs_hplan *plans, int nharms)
{
    if (!plans)
        return;
    for (int h = 0; h < nharms; h++) {
        free(plans[h].row);
        free(plans[h].W);
        free(plans[h].A);
        free(plans[h].carry);
        free(plans[h].grow);
        free(plans[h].goff);
        free(plans[h].gnj);
        free(plans[h].gW);
        free(plans[h].gA);
    }
    free(plans);
}

/* Residue and accumulated integer-bin count of harmonic hp->h at trial t:
 * floor(h r_t) = rfloor0 + qint + carry[row[res]]. */
static inline void hstate(const cs_hplan *hp, long long t, long long *res, long long *qint)
{
    /* t * h * pnum fits easily in 64 bits for any real search (t < 1e12). */
    long long T = t * hp->h * hp->pnum;
    *qint = T / hp->q;
    *res = T % hp->q;
}

long long cs_hbin(const cs_hplan *hp, long long t, long long *nyqbin)
{
    if (hp->P == 0) {
        /* The fallback's arithmetic, reproduced exactly -- including that its
         * Nyquist test floors without the 1e-15 nudge the window uses. */
        double r = (hp->rfloor0 + hp->dr0) + t * hp->dh;
        *nyqbin = (long long) floor(r);
        return (long long) floor(r + 1e-15);
    }
    long long res, qint;
    hstate(hp, t, &res, &qint);
    *nyqbin = hp->rfloor0 + qint;
    return hp->rfloor0 + qint + hp->carry[hp->row[res]];
}

int cs_harmonic_fits(const cs_hplan *hp, const cs_fft *fft, long long t)
{
    const int m2 = hp->m / 2;
    long long nyqbin, fb = cs_hbin(hp, t, &nyqbin);
    return fb + 1 - m2 >= 0 && fb + m2 <= fft->numamps - 1 && nyqbin < fft->N / 2;
}

/* The per-trial fallback: weights recomputed in double for every trial. */
static void fill_row_slow(const cs_hplan *hp, const cs_fft *fft, long long t0, int n,
                          fcomplex *out, int stride)
{
    const int m = hp->m, m2 = hp->m / 2;
    const double hr_lo = hp->rfloor0 + hp->dr0;
    for (int j = 0; j < n; j++) {
        double r = hr_lo + (t0 + j) * hp->dh;
        long long fb = (long long) floor(r + 1e-15);
        double dr = fmod(r, 1.0);
        fcomplex *o = out + (size_t) j * stride;
        if (dr == 0.0) {
            *o = fft->amps[fb];
            continue;
        }
        double sre = 0.0, sim = 0.0;
        for (int i = 0; i < m; i++) {
            int jj = i + 1 - m2;
            double w = 1.0 / (dr - jj);
            sre += w * fft->amps[fb + jj].r;
            sim += w * fft->amps[fb + jj].i;
        }
        double s = sin(M_PI * dr) / M_PI, ar = s * cos(M_PI * dr), ai = -s * sin(M_PI * dr);
        o->r = (float) (ar * sre - ai * sim);
        o->i = (float) (ar * sim + ai * sre);
    }
}

/* One group: the (nj x V) weight block times the nj-long slice of the bin
 * planes, as V independent accumulators per plane -- no gather and no
 * horizontal reduce.  (Measured: this plain form beats both a single
 * V-float GCC vector type, which gcc lowers badly onto 16-byte registers,
 * and an array of 4-float vectors.) */
static inline void group_lanes(const float *restrict w, int nj, const float *restrict re,
                               const float *restrict im, float *restrict ar,
                               float *restrict ai)
{
    enum { V = CS_GROUP_V };
    float sr[V], si[V];
    for (int k = 0; k < V; k++)
        sr[k] = si[k] = 0.0f;
    for (int j = 0; j < nj; j++) {
        const float rv = re[j], iv = im[j];
        const float *restrict wj = w + (long) j * V;
        for (int k = 0; k < V; k++) {
            sr[k] += wj[k] * rv;
            si[k] += wj[k] * iv;
        }
    }
    for (int k = 0; k < V; k++) {
        ar[k] = sr[k];
        ai[k] = si[k];
    }
}

int cs_fill_row(const cs_hplan *hp, const cs_fft *fft, long long t0, int n,
                fcomplex *out, int stride, float *re, float *im, long nplane)
{
    enum { V = CS_GROUP_V };
    if (n < 1)
        return 0;
    if (!(cs_harmonic_fits(hp, fft, t0) && cs_harmonic_fits(hp, fft, t0 + n - 1)))
        return 0;
    if (hp->P == 0) {
        fill_row_slow(hp, fft, t0, n, out, stride);
        return 1;
    }
    const int m2 = hp->m / 2;
    /* Groups are anchored to the GLOBAL trial grid (multiples of V), so a
     * trial's arithmetic does not depend on where chunk boundaries fall.  The
     * end groups may hang off the chunk: they are computed in full and stored
     * only where they land inside it, reading zeros past the ends of the file. */
    const long long T_first = (t0 / V) * V, T_last = ((t0 + n - 1) / V) * V + V - 1;
    long long res, qint, resx, qintx;
    hstate(hp, T_first, &res, &qint);
    hstate(hp, T_last, &resx, &qintx);
    const long long lo = hp->rfloor0 + qint + hp->carry[hp->row[res]] + 1 - m2;
    const long long hi = hp->rfloor0 + qintx + hp->carry[hp->row[resx]] + m2;
    const long nw = (long) (hi - lo + 1);
    if (nw > nplane) {
        fprintf(stderr, "cs_fill_row: harmonic %d spans %ld bins but the planes hold %ld\n",
                hp->h, nw, nplane);
        exit(1);
    }
    /* De-interleave the window, zero outside the file. */
    for (long i = 0; i < nw; i++) {
        const long long b = lo + i;
        if (b < 0 || b >= fft->numamps) {
            re[i] = im[i] = 0.0f;
        } else {
            re[i] = fft->amps[b].r;
            im[i] = fft->amps[b].i;
        }
    }
    float ar[V], ai[V];
    for (long long T = T_first; T <= T_last; T += V) {
        const int g = hp->grow[res];
        const long b0 = (long) (hp->rfloor0 + qint + hp->carry[hp->row[res]] + 1 - m2 - lo);
        group_lanes(hp->gW + hp->goff[g], hp->gnj[g], re + b0, im + b0, ar, ai);
        const fcomplex *A = hp->gA + (size_t) g * V;
        const long long j0 = T - t0;
        for (int k = 0; k < V; k++) {
            const long long j = j0 + k;
            if (j < 0 || j >= n)
                continue;
            fcomplex *o = out + (size_t) j * stride;
            o->r = A[k].r * ar[k] - A[k].i * ai[k];
            o->i = A[k].r * ai[k] + A[k].i * ar[k];
        }
        /* Advance V trials: exactly V steps of the per-trial recurrence. */
        const long long tot = res + V * hp->s;
        qint += V * hp->base_adv + tot / hp->q;
        res = tot % hp->q;
    }
    return 1;
}
