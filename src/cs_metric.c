/*
 * coherentsearch: the detection metric (CoherentSearch.jl's search.jl).
 *
 * Each folded profile is scored by the peak boxcar matched filter over a
 * geometric bank of widths, wrapping around the profile:
 *
 *   max_{w,p} (sum_{i=p}^{p+w-1} P_i - (w/nbins) sum_i P_i) * shape_w / sigma
 *
 * sigma is the per-bin noise scale of the chunk (computed from which
 * harmonics are filled, or measured with a MAD) and shape_w converts it to the
 * noise scale of a width-w sum of band-limited (correlated) profile bins.
 */
#include "coherentsearch.h"

int cs_boxcar_widths(int nbins, double fsp, double maxfrac, int *widths)
{
    if (nbins < 2) {
        fprintf(stderr, "cs_boxcar_widths: need nbins >= 2, got %d\n", nbins);
        exit(1);
    }
    int wmax = (int) floor(maxfrac * nbins);
    if (wmax < 1)
        wmax = 1;
    if (wmax > nbins - 1)
        wmax = nbins - 1;
    int nw = 0, w = 1;
    while (w <= wmax) {
        if (nw == 64) {
            fprintf(stderr, "cs_boxcar_widths: more than 64 widths\n");
            exit(1);
        }
        widths[nw++] = w;
        int next = (int) floor(fsp * w);
        w = next > w + 1 ? next : w + 1;
    }
    return nw;
}

/* Widths that the rest of the decimation ladder already covers are dropped:
 * w = 1 on all but the deepest fold (w = 2 one rung deeper is the same
 * filter at twice the harmonics), and w > CS_LADDER_WMAX on all but the
 * shallowest.  Only sound when the rungs are close enough together (each
 * k' <= (CS_LADDER_WMAX/2) k), and never when only one fold is searched. */
int cs_ladder_widths(int nbins, int k, const cs_params *p, int *widths)
{
    int all[64];
    int nall = cs_boxcar_widths(nbins, p->boxcar_fsp, p->boxcar_maxfrac, all);
    int ks[CS_MAXDECIM], nk = p->ndecim;
    memcpy(ks, p->decims, nk * sizeof(int));
    for (int i = 1; i < nk; i++)        /* insertion sort; nk is tiny */
        for (int j = i; j > 0 && ks[j - 1] > ks[j]; j--) {
            int t = ks[j];
            ks[j] = ks[j - 1];
            ks[j - 1] = t;
        }
    int prune = nk > 1;
    for (int i = 0; prune && i < nk - 1; i++)
        if (ks[i + 1] > (CS_LADDER_WMAX / 2) * ks[i])
            prune = 0;
    int nw = 0;
    if (prune) {
        int kmin = ks[0], kmax = ks[nk - 1];
        for (int i = 0; i < nall; i++)
            if ((all[i] >= 2 || k == kmin) && (all[i] <= CS_LADDER_WMAX || k == kmax))
                widths[nw++] = all[i];
    }
    if (nw == 0) {              /* no ladder, or nothing would be left */
        memcpy(widths, all, nall * sizeof(int));
        nw = nall;
    }
    return nw;
}

/* sqrt(2 nlow + 0.5 nnyq): each filled harmonic below the profile's own
 * Nyquist contributes real and imaginary parts of variance 1/2 to every
 * unnormalised c2r output bin (times 2 for the Hermitian pair); the Nyquist
 * harmonic contributes only its real part. */
double cs_analytic_sigma(const unsigned char *filled, int k, int Hk)
{
    int nlow = 0;
    for (int j = 1; j < Hk; j++)
        nlow += filled[j * k] ? 1 : 0;
    int nnyq = filled[Hk * k] ? 1 : 0;
    if (nlow + nnyq < 1)
        return 0.0;
    return sqrt(2.0 * nlow + 0.5 * nnyq);
}

/* var(S_w) = 2 sum_{j filled, j < Hk} |D_j|^2 + (1/2) |D_Hk|^2 [if filled],
 * with the Dirichlet kernel D_j = sin(pi j w / nbins) / sin(pi j / nbins). */
void cs_boxcar_shape(const unsigned char *filled, int k, int Hk,
                     const int *widths, int nwidths, double *shape)
{
    const int nbins = 2 * Hk;
    const double sigma = cs_analytic_sigma(filled, k, Hk);
    for (int i = 0; i < nwidths; i++) {
        const int w = widths[i];
        double v = 0.0;
        for (int j = 1; j < Hk; j++) {
            if (!filled[j * k])
                continue;
            double d = sin(M_PI * j * w / nbins) / sin(M_PI * j / nbins);
            v += 2.0 * d * d;
        }
        if (filled[Hk * k]) {
            double d = sin(M_PI * w / 2.0);
            v += 0.5 * d * d;
        }
        shape[i] = (sigma > 0.0 && v > 0.0) ? sigma / sqrt(v) : 0.0;
    }
}

/* Quickselect: put the kth smallest (0-based) of v[0..n-1] at v[k], with
 * everything before it no larger. */
static void select_kth(float *v, int n, int k)
{
    int lo = 0, hi = n - 1;
    while (hi > lo) {
        float pivot = v[lo + (hi - lo) / 2];
        int i = lo, j = hi;
        while (i <= j) {
            while (v[i] < pivot)
                i++;
            while (v[j] > pivot)
                j--;
            if (i <= j) {
                float t = v[i];
                v[i] = v[j];
                v[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j)
            hi = j;
        else if (k >= i)
            lo = i;
        else
            return;
    }
}

/* Median of v[0..n-1] (destroys the order): the mean of the two central
 * order statistics for even n. */
static float median_destroy(float *v, int n)
{
    int half = n / 2;
    select_kth(v, n, half);
    float upper = v[half];
    if (n % 2)
        return upper;
    float lower = v[0];
    for (int i = 1; i < half; i++)
        if (v[i] > lower)
            lower = v[i];
    return 0.5f * (lower + upper);
}

/* Pooled about zero, not the sample median: every profile's DC is held at
 * zero, so its bins have mean exactly 0.  Above CS_SIGMA_SAMPLES bins, whole
 * profiles are taken, evenly spaced through the chunk, and the last one
 * truncated to fit. */
float cs_block_sigma(const float *profs, int nbins, int n, float *buf)
{
    const long long N = (long long) nbins * n;
    if (N == 0)
        return 0.0f;
    int ns = 0;
    if (N <= CS_SIGMA_SAMPLES) {
        for (long long t = 0; t < N; t++)
            buf[t] = fabsf(profs[t]);
        ns = (int) N;
    } else {
        int ncol = CS_SIGMA_SAMPLES / nbins;
        if (ncol < 1)
            ncol = 1;
        if (ncol > n)
            ncol = n;
        for (int c = 0; c < ncol; c++) {
            const float *col = profs + (size_t) (((long long) c * n) / ncol) * nbins;
            int take = CS_SIGMA_SAMPLES - ns < nbins ? CS_SIGMA_SAMPLES - ns : nbins;
            if (take == 0)
                break;
            for (int i = 0; i < take; i++)
                buf[ns++] = fabsf(col[i]);
        }
    }
    return 1.4826f * median_destroy(buf, ns);
}

float cs_profile_boxcar(const float *prof, int nbins, const int *widths, int nw,
                        float invsigma, const double *shape, float *psum)
{
    if (!(invsigma > 0.0f))
        return 0.0f;            /* degenerate (flat or empty) chunk */
    const int wmax = widths[nw - 1];
    /* Prefix sums, tiled by wmax extra bins so a boxcar can wrap. */
    psum[0] = 0.0f;
    for (int i = 0; i < nbins + wmax; i++)
        psum[i + 1] = psum[i] + prof[i < nbins ? i : i - nbins];
    const float stot = psum[nbins] - psum[0];
    float best = -INFINITY;
    for (int iw = 0; iw < nw; iw++) {
        const int w = widths[iw];
        const float duty = (float) w / (float) nbins;
        const float invsw = (float) (invsigma * shape[iw]);
        float mx = psum[w] - psum[0];
        for (int p = 1; p < nbins; p++) {
            float s = psum[p + w] - psum[p];
            mx = s > mx ? s : mx;
        }
        float cand = (mx - duty * stot) * invsw;
        if (cand > best)
            best = cand;
    }
    return best;
}
