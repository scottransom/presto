/*
 * coherentsearch's FFTW plan, shared with makewisdom (which must plan exactly
 * the same transforms for its wisdom to apply).  Depends only on FFTW.
 */
#ifndef CS_FFTPLAN_H
#define CS_FFTPLAN_H

#include "fftw3.h"

/* The batched c2r for decimation k: all nprof harmonic stacks (element j of
 * stack b at in[b*(nharms+1) + j*k]) into nprof 2*(nharms/k)-bin profiles.
 * makewisdom plans exactly this for the default nharms and blocksize, so the
 * two must stay identical. */
static inline fftwf_plan cs_plan_c2r(int nharms, int nprof, int k, fftwf_complex *in,
                                     float *out, unsigned flags)
{
    int nbins = 2 * (nharms / k);
    return fftwf_plan_many_dft_c2r(1, &nbins, nprof, in, NULL, k, nharms + 1,
                                   out, NULL, 1, nbins, flags | FFTW_PRESERVE_INPUT);
}

#define CS_DEFAULT_NHARMS 60
#define CS_DEFAULT_BLOCKSIZE 2048
#define CS_DEFAULT_MAXDECIM 6

#endif
