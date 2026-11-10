/*
 * coherentsearch: the full search over a band, and candidate handling
 * (CoherentSearch.jl's search, chunk_starts, _sigma_sanity_check,
 * remove_duplicates, remove_harmonics, measure_ducy and write_candidates).
 */
#include "coherentsearch.h"
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void cs_candlist_push(cs_candlist *l, cs_cand c)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? 2 * l->cap : 256;
        l->c = realloc(l->c, l->cap * sizeof(cs_cand));
        if (!l->c) {
            perror("coherentsearch: candidate list");
            exit(1);
        }
    }
    l->c[l->n++] = c;
}

void cs_candlist_free(cs_candlist *l)
{
    free(l->c);
    l->c = NULL;
    l->n = l->cap = 0;
}

cs_band cs_make_band(const cs_fft *fft, const cs_params *p, double lofreq,
                     double hifreq, int lobin)
{
    cs_band b;
    b.r_lo = lobin >= 0 ? (double) lobin : lofreq * fft->T;
    b.r_hi = hifreq * fft->T;
    long long total = (long long) floor((b.r_hi - b.r_lo) / cs_lodr(p)) + 1;
    b.total = total > 0 ? total : 0;
    return b;
}

/* --- chunk boundaries ------------------------------------------------------ */

/* First t in 0..total-1 with pred(t) == want, for a pred that is monotone in t
 * (total if there is none). */
static long long first_where(int (*pred)(const cs_hplan *, const cs_fft *, long long),
                             const cs_hplan *hp, const cs_fft *fft, long long total, int want)
{
    if (total < 1)
        return total;
    if ((pred(hp, fft, 0) != 0) == want)
        return 0;
    if ((pred(hp, fft, total - 1) != 0) != want)
        return total;
    long long lo = 0, hi = total - 1;
    while (hi - lo > 1) {
        long long mid = lo + (hi - lo) / 2;
        if ((pred(hp, fft, mid) != 0) == want)
            hi = mid;
        else
            lo = mid;
    }
    return hi;
}

/* The two monotone halves of cs_harmonic_fits. */
static int fits_lo(const cs_hplan *hp, const cs_fft *fft, long long t)
{
    long long nyq;
    (void) fft;
    return cs_hbin(hp, t, &nyq) + 1 - hp->m / 2 >= 0;
}

static int fits_hi(const cs_hplan *hp, const cs_fft *fft, long long t)
{
    long long nyq, fb = cs_hbin(hp, t, &nyq);
    return fb + hp->m / 2 <= fft->numamps - 1 && nyq < fft->N / 2;
}

static int cmp_ll(const void *a, const void *b)
{
    long long x = *(const long long *) a, y = *(const long long *) b;
    return (x > y) - (x < y);
}

/* One flag per harmonic per chunk is only right if every trial in a chunk
 * agrees on which harmonics fit, so chunks are also split wherever that
 * changes.  Each harmonic is usable on one contiguous interval, found by
 * bisection; that adds at most 2*nharms boundaries. */
long long *cs_chunk_starts(const cs_hplan *plans, const cs_fft *fft, const cs_params *p,
                           long long total, int nprof, long long *nchunks)
{
    if (total < 1) {
        *nchunks = 0;
        return NULL;
    }
    long long nreg = (total + nprof - 1) / nprof;
    long long *s = malloc((nreg + 2 * p->nharms) * sizeof(long long));
    long long n = 0;
    for (long long t = 0; t < total; t += nprof)
        s[n++] = t;
    if (total > 1) {
        for (int h = 0; h < p->nharms; h++) {
            long long a = first_where(fits_lo, plans + h, fft, total, 1);
            long long b = first_where(fits_hi, plans + h, fft, total, 0);
            if (a > 0 && a < total)
                s[n++] = a;
            if (b > 0 && b < total)
                s[n++] = b;
        }
    }
    qsort(s, n, sizeof(long long), cmp_ll);
    long long u = 0;
    for (long long i = 0; i < n; i++)
        if (u == 0 || s[i] != s[u - 1])
            s[u++] = s[i];
    *nchunks = u;
    return s;
}

/* --- the normalisation check ---------------------------------------------- */

#define CS_NORM_CHECK_LO 10000     /* first bin checked */
#define CS_NORM_CHECK_N 10000      /* bins checked */
#define CS_NORM_CHECK_TOL 0.10     /* tolerated error in the noise scale */

static int cmp_float(const void *a, const void *b)
{
    float x = *(const float *) a, y = *(const float *) b;
    return (x > y) - (x < y);
}

/* The analytic noise scale assumes the input is normalised (Fourier powers of
 * mean 1), and fails silently -- inflating every S/N -- when it is not.  Check
 * that on a block of bins away from the lowest frequencies, where leftover red
 * noise (rednoise's running median rests on few bins there) would raise false
 * alarms.  Normalised powers are exponential with mean 1, so median/ln 2 should
 * be 1; the median shrugs off a pulsar or RFI in the block.  10000 bins give
 * a scatter of about 1.4% in the power, 0.7% in the noise scale. */
int cs_check_normalisation(const cs_fft *fft, const char *fname)
{
    long long lo = CS_NORM_CHECK_LO, n = CS_NORM_CHECK_N;
    if (fft->numamps < lo + n) {    /* short FFT: the top half, at most n bins */
        lo = fft->numamps / 2;
        n = fft->numamps - lo < n ? fft->numamps - lo : n;
    }
    if (n < 100)
        return 1;
    float *pow = malloc(n * sizeof(float));
    for (long long i = 0; i < n; i++) {
        const fcomplex a = fft->amps[lo + i];
        pow[i] = a.r * a.r + a.i * a.i;
    }
    qsort(pow, n, sizeof(float), cmp_float);
    const double med = n % 2 ? pow[n / 2] : 0.5 * (pow[n / 2 - 1] + pow[n / 2]);
    free(pow);
    const double ratio = sqrt(med / M_LN2);   /* (actual / assumed) noise scale */
    if (fabs(ratio - 1.0) <= CS_NORM_CHECK_TOL)
        return 1;
    fprintf(stderr, "\nWARNING (%s): the Fourier powers do not look normalised: the noise\n"
            "  scale is ", fname);
    /* Quote a ratio, not a percentage, when it is far off: a 1000x error
     * would otherwise read like a 99.9% one. */
    if (ratio > 2.0 || ratio < 0.5)
        fprintf(stderr, "off by a factor of %.3g", ratio > 1.0 ? ratio : 1.0 / ratio);
    else
        fprintf(stderr, "%.1f%% too %s", 100.0 * fabs(ratio - 1.0), ratio > 1.0 ? "high" : "low");
    fprintf(stderr, ", and every reported S/N is off by that factor.\n"
            "  Normalised powers have mean 1 (median ln 2 = 0.693); bins %lld-%lld have\n"
            "  median %.4g.  Run PRESTO's rednoise on the FFT, or use -measuredsigma.\n\n",
            lo, lo + n - 1, med);
    return 0;
}

/* --- the search ------------------------------------------------------------ */

cs_searcher *cs_searcher_new(const cs_params *p, int nprof, int nthreads)
{
    cs_searcher *s = calloc(1, sizeof(cs_searcher));
    s->p = *p;
    s->nprof = nprof;
    s->nthreads = nthreads < 1 ? 1 : nthreads;
#ifndef _OPENMP
    s->nthreads = 1;
#endif
    s->ws = malloc(s->nthreads * sizeof(cs_workspace *));
    for (int t = 0; t < s->nthreads; t++)
        s->ws[t] = cs_workspace_new(p, nprof);
    return s;
}

void cs_searcher_free(cs_searcher *s)
{
    if (!s)
        return;
    for (int t = 0; t < s->nthreads; t++)
        cs_workspace_free(s->ws[t]);
    free(s->ws);
    free(s);
}

/* One chunk: fill, transform and score every decimation, appending the
 * trials above threshold to out. */
static void search_chunk(cs_workspace *ws, const cs_hplan *plans, const cs_fft *fft,
                         const cs_params *p, double r_lo, long long i0, int n,
                         double threshold, cs_candlist *out)
{
    const double lodr = cs_lodr(p);
    const double rstart = r_lo + i0 * lodr;
    cs_fill_chunk(ws, plans, fft, p, i0, n);
    for (int di = 0; di < ws->ndecim; di++) {
        cs_decim *d = ws->d + di;
        cs_prepare_decim(ws, d, fft, p, rstart, n, p->sigma_measured);
        cs_score_decim(d, p->sigma_measured ? d->sigma_measured : d->sigma_analytic);
        for (int j = 0; j < d->nvalid; j++) {
            if (d->mvals[j] > threshold) {
                cs_cand cand;
                cand.r = d->k == 1 ? rstart + j * lodr : d->k * (rstart + j * lodr);
                cand.freq = cand.r / fft->T;
                cand.metric = d->mvals[j];
                cand.nharm = d->Hk;
                cand.ducy = NAN;
                cand.note[0] = '\0';
                cs_candlist_push(out, cand);
            }
        }
    }
}

void cs_search(cs_searcher *s, const cs_fft *fft, const cs_band *band, double threshold,
               int progress, cs_candlist *out)
{
    if (band->total < 1)
        return;
    const cs_params *p = &s->p;
    /* The phase tables key off the global trial index, so they are built once
     * against r_lo and shared, read-only, by every thread. */
    cs_hplan *plans = cs_build_plans(p, band->r_lo);
    long long nchunks;
    long long *starts = cs_chunk_starts(plans, fft, p, band->total, s->nprof, &nchunks);
    /* Each chunk keeps its own candidates and they are joined in chunk order,
     * so the output does not depend on the number of threads. */
    cs_candlist *found = calloc(nchunks, sizeof(cs_candlist));
    long long done = 0;
    int lastpct = -1;
    const int nt = (int) (s->nthreads < nchunks ? s->nthreads : nchunks);
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        cs_workspace *ws = s->ws[tid];
        /* Round-robin: whole chunks, one at a time per thread. */
        for (long long c = tid; c < nchunks; c += nt) {
            const long long i0 = starts[c];
            const int n = (int) ((c < nchunks - 1 ? starts[c + 1] : band->total) - i0);
            search_chunk(ws, plans, fft, p, band->r_lo, i0, n, threshold, found + c);
            long long ndone;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
            ndone = ++done;
            /* One thread owns the meter, so lines do not interleave. */
            if (progress && tid == 0) {
                int pct = (int) (100 * ndone / nchunks);
                if (pct != lastpct) {
                    printf("\r%3d%% searched", pct);
                    fflush(stdout);
                    lastpct = pct;
                }
            }
        }
    }
    if (progress) {
        printf("\r100%% searched\n");
        fflush(stdout);
    }
    for (long long c = 0; c < nchunks; c++) {
        for (long i = 0; i < found[c].n; i++)
            cs_candlist_push(out, found[c].c[i]);
        cs_candlist_free(found + c);
    }
    free(found);
    free(starts);
    cs_free_plans(plans, p->nharms);
}

/* --- candidate handling ---------------------------------------------------- */

/* Stable sorts: qsort is not stable, so the original position breaks ties
 * (Julia's sorts are stable, and ties decide which of two equal candidates
 * survives a collapse). */
typedef struct {
    cs_cand c;
    long idx;
} cs_idxcand;

static int cmp_r(const void *a, const void *b)
{
    const cs_idxcand *x = a, *y = b;
    if (x->c.r != y->c.r)
        return x->c.r < y->c.r ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

static int cmp_freq(const void *a, const void *b)
{
    const cs_idxcand *x = a, *y = b;
    if (x->c.freq != y->c.freq)
        return x->c.freq < y->c.freq ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

static int cmp_metric_desc(const void *a, const void *b)
{
    const cs_idxcand *x = a, *y = b;
    if (x->c.metric != y->c.metric)
        return x->c.metric > y->c.metric ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

static void stable_sort(cs_candlist *l, int (*cmp)(const void *, const void *))
{
    if (l->n < 2)
        return;
    cs_idxcand *t = malloc(l->n * sizeof(cs_idxcand));
    for (long i = 0; i < l->n; i++) {
        t[i].c = l->c[i];
        t[i].idx = i;
    }
    qsort(t, l->n, sizeof(cs_idxcand), cmp);
    for (long i = 0; i < l->n; i++)
        l->c[i] = t[i].c;
    free(t);
}

void cs_sort_by_metric(cs_candlist *l)
{
    stable_sort(l, cmp_metric_desc);
}

void cs_remove_duplicates(cs_candlist *l, double dr_tol)
{
    if (l->n == 0)
        return;
    stable_sort(l, cmp_r);
    long nk = 0;
    cs_cand best = l->c[0];
    double prev_r = best.r;
    for (long i = 1; i < l->n; i++) {
        cs_cand c = l->c[i];
        /* A chain: each candidate within dr_tol of the previous one joins the
         * current cluster, however far that drifts from its start. */
        if (c.r - prev_r <= dr_tol) {
            if (c.metric > best.metric)
                best = c;
        } else {
            l->c[nk++] = best;
            best = c;
        }
        prev_r = c.r;
    }
    l->c[nk++] = best;
    l->n = nk;
    stable_sort(l, cmp_freq);
}

/* Whether r1 and r2 are harmonics of a common fundamental: hi/lo ~ n/m with
 * 1 <= m, n <= numharm, holding hi to tol bins on the shared comb. */
static int harmonically_related(double r1, double r2, int numharm, double tol)
{
    double lo = r1 < r2 ? r1 : r2, hi = r1 < r2 ? r2 : r1;
    if (!(lo > 0.0))
        return 0;
    for (int m = 1; m <= numharm; m++) {
        double n = nearbyint(m * hi / lo);      /* round half to even, as Julia */
        if (n < 1 || n > numharm)
            continue;
        if (fabs(m * hi - n * lo) <= tol * m)
            return 1;
    }
    return 0;
}

void cs_remove_harmonics(cs_candlist *l, int numharm, double tol)
{
    if (l->n == 0)
        return;
    stable_sort(l, cmp_metric_desc);    /* strongest first */
    long nk = 0;
    for (long i = 0; i < l->n; i++) {
        int related = 0;
        for (long j = 0; j < nk && !related; j++)
            related = harmonically_related(l->c[i].r, l->c[j].r, numharm, tol);
        if (!related)
            l->c[nk++] = l->c[i];
    }
    l->n = nk;
    stable_sort(l, cmp_freq);
}

long cs_remove_harmonics_of(cs_candlist *l, const cs_candlist *ref, int numharm, double tol)
{
    long nk = 0;
    for (long i = 0; i < l->n; i++) {
        int related = 0;
        for (long j = 0; j < ref->n && !related; j++)
            related = ref->c[j].metric > l->c[i].metric &&
                harmonically_related(l->c[i].r, ref->c[j].r, numharm, tol);
        if (!related)
            l->c[nk++] = l->c[i];
    }
    const long nremoved = l->n - nk;
    l->n = nk;
    return nremoved;
}

/* Exact Fourier interpolation at r with an m-bin kernel, in double
 * (CoherentSearch.jl's fourier_interp). */
static void interp_exact(const cs_fft *fft, double r, int m, double *re, double *im)
{
    const int m2 = m / 2;
    const long long fb = (long long) floor(r + 1e-15);
    const double dr = fmod(r, 1.0);
    double sr = 0.0, si = 0.0;
    for (int j = 1 - m2; j <= m2; j++) {
        double x = dr - j;
        double sinc = x == 0.0 ? 1.0 : sin(M_PI * x) / (M_PI * x);
        double cr = sinc * cos(M_PI * x), ci = -sinc * sin(M_PI * x);   /* conj */
        fcomplex a = fft->amps[fb + j];
        sr += cr * a.r - ci * a.i;
        si += cr * a.i + ci * a.r;
    }
    *re = sr;
    *im = si;
}

/* The exactly interpolated (m = 64) complex amplitudes of harmonics
 * 1..nharm of r into xr[h], xi[h], stopping at the first whose kernel
 * crosses Nyquist or the end of the data.  Returns the last one filled. */
static int candidate_harmonics(const cs_fft *fft, double r, int nharm, double *xr, double *xi)
{
    const int m = 64, m2 = m / 2;
    int H = 0;
    for (int h = 1; h <= nharm; h++) {
        double rh = r * h;
        long long rint = (long long) floor(rh + 1e-15) + 1;
        if (!(rint - m2 >= 0 && rint + m2 <= fft->numamps && rh < (double) (fft->N / 2)))
            break;
        interp_exact(fft, rh, m, xr + h, xi + h);
        H = h;
    }
    return H;
}

int cs_candidate_profile(const cs_fft *fft, double r, int nharm, double *prof)
{
    double *xr = calloc(nharm + 1, sizeof(double)), *xi = calloc(nharm + 1, sizeof(double));
    const int H = candidate_harmonics(fft, r, nharm, xr, xi);
    const int nbins = H >= 1 ? 2 * H : 2;
    /* Inverse real DFT of the stack, DC held at zero (normalised like irfft). */
    for (int i = 0; i < nbins; i++) {
        double s = 0.0;
        for (int h = 1; h < H; h++) {
            double ph = 2.0 * M_PI * h * i / nbins;
            s += 2.0 * (xr[h] * cos(ph) - xi[h] * sin(ph));
        }
        if (H >= 1)
            s += xr[H] * ((i % 2) ? -1.0 : 1.0);    /* the profile's own Nyquist term */
        prof[i] = s / nbins;
    }
    free(xr);
    free(xi);
    return nbins;
}

long cs_remove_periodic(const cs_fft *fft, cs_candlist *l, int nlow, int minharm,
                        double ratio)
{
    long nk = 0;
    for (long i = 0; i < l->n; i++) {
        const cs_cand *c = l->c + i;
        int bad = 0;
        if (c->nharm >= minharm && c->r > 0.0) {
            double *xr = calloc(c->nharm + 1, sizeof(double));
            double *xi = calloc(c->nharm + 1, sizeof(double));
            const int H = candidate_harmonics(fft, c->r, c->nharm, xr, xi);
            double lowmax = 0.0, topmax = 0.0;
            int toph = 0;
            for (int h = 1; h <= H; h++) {
                const double pw = xr[h] * xr[h] + xi[h] * xi[h];
                if (h <= nlow && pw > lowmax)
                    lowmax = pw;
                if (pw > topmax) {
                    topmax = pw;
                    toph = h;
                }
            }
            bad = toph >= minharm && topmax > ratio * lowmax;
            free(xr);
            free(xi);
        }
        if (!bad)
            l->c[nk++] = *c;
    }
    const long nremoved = l->n - nk;
    l->n = nk;
    return nremoved;
}

/* The width of the best boxcar on a candidate's profile, refolded at its own
 * r with its own nharm harmonics.  Its noise scale does not matter: it scales
 * every width's score alike. */
static double candidate_ducy(const cs_fft *fft, const cs_params *p, double r, int nharm)
{
    double *prof = malloc((2 * nharm + 2) * sizeof(double));
    const int nbins = cs_candidate_profile(fft, r, nharm, prof);
    int widths[64];
    int nw = cs_boxcar_widths(nbins, p->boxcar_fsp, p->boxcar_maxfrac, widths);
    unsigned char *filled = malloc(nbins / 2 + 1);
    memset(filled, 1, nbins / 2 + 1);
    double shape[64];
    cs_boxcar_shape(filled, 1, nbins / 2, widths, nw, shape);
    /* Prefix sums about the mean (conditioning only), tiled to wrap. */
    double mean = 0.0;
    for (int i = 0; i < nbins; i++)
        mean += prof[i];
    mean /= nbins;
    double *psum = malloc((nbins + widths[nw - 1] + 1) * sizeof(double));
    psum[0] = 0.0;
    for (int i = 0; i < nbins + widths[nw - 1]; i++)
        psum[i + 1] = psum[i] + (prof[i < nbins ? i : i - nbins] - mean);
    const double stot = psum[nbins] - psum[0];
    double best = -INFINITY;
    int bestw = widths[0];
    for (int iw = 0; iw < nw; iw++) {
        int w = widths[iw];
        double mx = psum[w] - psum[0];
        for (int q = 1; q < nbins; q++)
            if (psum[q + w] - psum[q] > mx)
                mx = psum[q + w] - psum[q];
        double cand = (mx - ((double) w / nbins) * stot) * shape[iw];
        if (cand > best) {      /* strict: the narrowest wins a tie */
            best = cand;
            bestw = w;
        }
    }
    free(prof);
    free(filled);
    free(psum);
    return (double) bestw / nbins;
}

void cs_measure_ducy(const cs_fft *fft, const cs_params *p, cs_candlist *l)
{
    for (long i = 0; i < l->n; i++) {
        cs_cand *c = l->c + i;
        c->ducy = (c->nharm >= 1 && c->r > 0.0) ? candidate_ducy(fft, p, c->r, c->nharm) : NAN;
    }
}

int cs_identify_known(const cs_fft *fft, cs_candlist *l)
{
    const infodata *id = &fft->idata;
    /* Same guards as accelsearch: real data with a position. */
    if (strncmp(id->telescope, "None", 4) == 0 || id->telescope[0] == '\0')
        return 0;
    if (hms2rad(id->ra_h, id->ra_m, id->ra_s) == 0.0 ||
        dms2rad(id->dec_d, id->dec_m, id->dec_s) == 0.0)
        return 0;
    char *dbnm = presto_data_path("pulsars.cat");
    FILE *db = fopen(dbnm, "rb");
    if (!db) {
        fprintf(stderr, "coherentsearch: no pulsar database (%s), so no known-pulsar "
                "labels\n", dbnm);
        free(dbnm);
        return 0;
    }
    fclose(db);
    free(dbnm);
    infodata idcopy = *id;      /* comp_psr_to_cand takes a non-const pointer */
    int nfound = 0;
    for (long i = 0; i < l->n; i++) {
        fourierprops props;
        memset(&props, 0, sizeof(props));
        /* Only r matters for the match (z always matches); give the rest
         * harmless values for the error estimates it also computes. */
        props.r = l->c[i].r;
        props.rerr = props.zerr = props.werr = 1.0;
        props.pow = props.rawpow = props.locpow = 1.0;
        char out[256] = "";
        if (comp_psr_to_cand(&props, &idcopy, out, 0)) {
            snprintf(l->c[i].note, sizeof(l->c[i].note), "%s", out);
            nfound++;
        } else {
            l->c[i].note[0] = '\0';
        }
    }
    return nfound;
}

#define CS_KNOWN_HDR "   Known PSR"

void cs_write_cohout(const cs_candlist *l, FILE *f)
{
    int anynote = 0;
    for (long i = 0; i < l->n; i++)
        anynote |= l->c[i].note[0] != '\0';
    /* The leading column is the rank, and the header says so. */
    fprintf(f, "#Num    'S/N'      Frequency (Hz)        Period (ms)    #Harm  Ducy(%%)%s\n",
            anynote ? CS_KNOWN_HDR : "");
    for (long i = 0; i < l->n; i++) {
        const cs_cand *c = l->c + i;
        char ducy[32];
        if (isnan(c->ducy))
            strcpy(ducy, "-");
        else
            snprintf(ducy, sizeof(ducy), "%.2f", 100.0 * c->ducy);
        fprintf(f, "%-4ld  %8.2f  %18.12f  %18.12f   %3d   %6s", i + 1, c->metric, c->freq,
                1000.0 / c->freq, c->nharm, ducy);
        if (c->note[0])
            fprintf(f, "   %s", c->note);
        fprintf(f, "\n");
    }
}

void cs_read_cohout(const char *filenm, double T, long maxcands, cs_candlist *l)
{
    FILE *f = chkfopen((char *) filenm, "r");
    char line[1024];
    while (fgets(line, sizeof(line), f) && l->n < maxcands) {
        if (line[0] == '#')
            continue;
        long rank;
        cs_cand c;
        double period;
        char ducy[64];
        int nf = sscanf(line, "%ld %lf %lf %lf %d %63s", &rank, &c.metric, &c.freq, &period,
                        &c.nharm, ducy);
        if (nf < 5 || !(c.freq > 0.0))
            continue;
        c.r = c.freq * T;
        c.ducy = (nf == 6 && strcmp(ducy, "-") != 0) ? atof(ducy) / 100.0 : NAN;
        c.note[0] = '\0';      /* recomputed by cs_identify_known */
        cs_candlist_push(l, c);
    }
    fclose(f);
}

/* End of the nf'th whitespace-separated field of line. */
static char *field_end(char *line, int nf)
{
    char *q = line;
    for (int i = 0; i < nf; i++) {
        while (*q == ' ' || *q == '\t')
            q++;
        while (*q && !isspace((unsigned char) *q))
            q++;
    }
    return q;
}

long cs_relabel_cohout(const char *filenm, const cs_fft *fft, int *nknown)
{
    FILE *f = chkfopen((char *) filenm, "r");
    char **lines = NULL;
    long *row = NULL;               /* each line's candidate, or -1 */
    long nlines = 0, cap = 0;
    cs_candlist l = { 0 };
    char buf[1024];
    while (fgets(buf, sizeof(buf), f)) {
        if (nlines == cap) {
            cap = cap ? 2 * cap : 256;
            lines = realloc(lines, cap * sizeof(char *));
            row = realloc(row, cap * sizeof(long));
        }
        buf[strcspn(buf, "\n")] = '\0';
        row[nlines] = -1;
        lines[nlines] = strdup(buf);
        if (buf[0] != '#') {
            long rank;
            cs_cand c;
            double period;
            char ducy[64];
            int nf = sscanf(buf, "%ld %lf %lf %lf %d %63s", &rank, &c.metric, &c.freq,
                            &period, &c.nharm, ducy);
            if (nf >= 5 && c.freq > 0.0) {
                c.r = c.freq * fft->T;
                c.ducy = NAN;
                c.note[0] = '\0';
                /* Drop any old label: keep the row's text up to its last numeric field. */
                *field_end(lines[nlines], nf) = '\0';
                row[nlines] = l.n;
                cs_candlist_push(&l, c);
            }
        } else if (strncmp(buf, "#Num", 4) == 0) {
            char *k = strstr(lines[nlines], CS_KNOWN_HDR);
            if (k)
                *k = '\0';
        }
        nlines++;
    }
    fclose(f);
    *nknown = cs_identify_known(fft, &l);

    /* Write beside it and rename, so a failure cannot leave it half-written. */
    char *tmpnm = malloc(strlen(filenm) + 8);
    sprintf(tmpnm, "%s.tmp", filenm);
    FILE *out = chkfopen(tmpnm, "w");
    for (long i = 0; i < nlines; i++) {
        fputs(lines[i], out);
        if (row[i] >= 0 && l.c[row[i]].note[0])
            fprintf(out, "   %s", l.c[row[i]].note);
        else if (*nknown && strncmp(lines[i], "#Num", 4) == 0)
            fputs(CS_KNOWN_HDR, out);
        fputc('\n', out);
        free(lines[i]);
    }
    if (fclose(out) != 0 || rename(tmpnm, filenm) != 0) {
        perror(filenm);
        exit(1);
    }
    free(tmpnm);
    free(lines);
    free(row);
    const long n = l.n;
    cs_candlist_free(&l);
    return n;
}
