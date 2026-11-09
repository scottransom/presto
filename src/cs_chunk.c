/*
 * coherentsearch: one chunk of trial fundamentals -- interpolate every
 * harmonic, inverse-FFT the stacks for each decimation factor k, and score
 * them (CoherentSearch.jl's fill_chunk_profiles!, decim_pass! and
 * _search_region!).
 */
#include "coherentsearch.h"

cs_workspace *cs_workspace_new(const cs_params *p, int nprof)
{
    cs_workspace *ws = calloc(1, sizeof(cs_workspace));
    ws->nharms = p->nharms;
    ws->nprof = nprof;
    ws->ndecim = p->ndecim;
    ws->ftprofs = gen_cvect((long) nprof * (p->nharms + 1));
    ws->filled = calloc(p->nharms + 1, 1);
    ws->sigbuf = gen_fvect(CS_SIGMA_SAMPLES);
    ws->nplane = cs_plane_len(p, nprof);
    ws->re = gen_fvect(ws->nplane);
    ws->im = gen_fvect(ws->nplane);
    ws->d = calloc(p->ndecim, sizeof(cs_decim));
    ws->wisdom = 1;
    for (int i = 0; i < p->ndecim; i++) {
        cs_decim *d = ws->d + i;
        d->k = p->decims[i];
        d->Hk = p->nharms / d->k;
        d->nbins = 2 * d->Hk;
        d->nwidths = cs_ladder_widths(d->nbins, d->k, p, d->widths);
        d->profs = gen_fvect((long) nprof * d->nbins);
        d->psum = gen_fvect(d->nbins + d->widths[d->nwidths - 1] + 1);
        d->mvals = gen_fvect(nprof);
        d->bcX = gen_fvect((long) d->nbins * CS_BC_BATCH);
        d->bcT = gen_fvect((long) (d->nbins + d->widths[d->nwidths - 1] + 1) * CS_BC_BATCH);
        d->bcwork = gen_fvect(3 * CS_BC_BATCH);
        d->shape_filled = calloc(p->nharms + 1, 1);
        d->shape_valid = 0;
        /* One batched c2r for all nprof stacks, preserving its input: the
         * decimated passes read ftprofs after the base pass.  A plan from
         * wisdom (makewisdom stores the default geometry) is ~1.8x faster
         * than FFTW_ESTIMATE's.  Measuring here instead would make results
         * depend on timings, i.e. differ in the last bits from run to run,
         * so without wisdom fall back to ESTIMATE, which is deterministic. */
        d->plan = cs_plan_c2r(p->nharms, nprof, d->k, (fftwf_complex *) ws->ftprofs,
                              d->profs, FFTW_MEASURE | FFTW_WISDOM_ONLY);
        if (!d->plan) {
            ws->wisdom = 0;
            d->plan = cs_plan_c2r(p->nharms, nprof, d->k, (fftwf_complex *) ws->ftprofs,
                                  d->profs, FFTW_ESTIMATE);
        }
        if (!d->plan) {
            fprintf(stderr, "coherentsearch: FFTW could not plan the k=%d transform\n", d->k);
            exit(1);
        }
    }
    /* Planning may write into the arrays it plans on (it does for anything
     * above FFTW_ESTIMATE), so clear them only after the last plan. */
    memset(ws->ftprofs, 0, sizeof(fcomplex) * (size_t) nprof * (p->nharms + 1));
    return ws;
}

void cs_workspace_free(cs_workspace *ws)
{
    if (!ws)
        return;
    for (int i = 0; i < ws->ndecim; i++) {
        cs_decim *d = ws->d + i;
        fftwf_destroy_plan(d->plan);
        vect_free(d->profs);
        vect_free(d->psum);
        vect_free(d->mvals);
        vect_free(d->bcX);
        vect_free(d->bcT);
        vect_free(d->bcwork);
        free(d->shape_filled);
    }
    free(ws->d);
    vect_free(ws->ftprofs);
    vect_free(ws->sigbuf);
    vect_free(ws->re);
    vect_free(ws->im);
    free(ws->filled);
    free(ws);
}

void cs_fill_chunk(cs_workspace *ws, const cs_hplan *plans, const cs_fft *fft,
                   const cs_params *p, long long t0, int n)
{
    if (n > ws->nprof) {
        fprintf(stderr, "cs_fill_chunk: chunk of %d trials > workspace's %d\n", n, ws->nprof);
        exit(1);
    }
    const int stride = p->nharms + 1;
    /* Every stack, not just the first n: the transforms run on all nprof. */
    memset(ws->ftprofs, 0, sizeof(fcomplex) * (size_t) ws->nprof * stride);
    ws->filled[0] = 0;
    for (int h = 1; h <= p->nharms; h++)
        ws->filled[h] = cs_fill_row(plans + h - 1, fft, t0, n, ws->ftprofs + h, stride,
                                    ws->re, ws->im, ws->nplane);
}

void cs_prepare_decim(cs_workspace *ws, cs_decim *d, const cs_fft *fft,
                      const cs_params *p, double rstart, int n, int want_measured)
{
    fftwf_execute(d->plan);
    if (d->k == 1) {
        d->nvalid = n;
    } else {
        /* Decimated trials whose fundamental k*r reaches Nyquist are skipped;
         * k*r grows with the trial, so the valid ones are a prefix. */
        const double lodr = cs_lodr(p), nyq = fft->N / 2.0;
        int nv = 0;
        while (nv < n && d->k * (rstart + nv * lodr) < nyq)
            nv++;
        d->nvalid = nv;
    }
    d->sigma_analytic = cs_analytic_sigma(ws->filled, d->k, d->Hk);
    d->sigma_measured = want_measured ? cs_block_sigma(d->profs, d->nbins, d->nvalid,
                                                       ws->sigbuf) : 0.0f;
    /* The shape factors depend only on which harmonics are filled, which
     * changes a handful of times in a whole search. */
    if (!d->shape_valid || memcmp(d->shape_filled, ws->filled, ws->nharms + 1) != 0) {
        cs_boxcar_shape(ws->filled, d->k, d->Hk, d->widths, d->nwidths, d->shape);
        memcpy(d->shape_filled, ws->filled, ws->nharms + 1);
        d->shape_valid = 1;
    }
}

void cs_score_decim(cs_decim *d, double sigma)
{
    const float sig = (float) sigma;
    const float invsigma = sig > 0.0f ? 1.0f / sig : 0.0f;
    int j = 0;
    for (; j + CS_BC_BATCH <= d->nvalid; j += CS_BC_BATCH)
        cs_boxcar_batch(d->profs + (size_t) j * d->nbins, d->nbins, d->widths, d->nwidths,
                        invsigma, d->shape, d->bcX, d->bcT, d->bcwork, d->mvals + j);
    for (; j < d->nvalid; j++)      /* the leftover partial batch */
        d->mvals[j] = cs_profile_boxcar(d->profs + (size_t) j * d->nbins, d->nbins,
                                        d->widths, d->nwidths, invsigma, d->shape, d->psum);
}
