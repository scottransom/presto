#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "presto.h"
#include "backend_common.h"

// Bandpass flattening ("-flatten") for prepdata, prepsubband and
// mpiprepsubband.
//
// Each raw channel is turned into (x - mean_c(t)) * w_c, where
// w_c = 1/sigma_c for good channels and 0 for bad ones.  sigma_c is
// fixed for the whole observation (relative channel variances barely
// change in time, and weights that change in time would put steps, and
// so birdies, into the data).  The means follow slow drifts by linear
// interpolation between estimates spread through the observation.
//
// The statistics come from the rfifind .stats file that goes with the
// -mask file, if there is one.  If not, a short pre-pass reads blocks
// spread evenly through the observation and measures them, keeping only
// the spectra that pass the same zero-DM cut that clip_times() uses.
//
// Bad channels get zero weight: dead channels (sigma ~ 0, e.g. band
// edges and gaps), channels whose sigma/mean departs from the radiometer
// value (used only if the data look like total power), channels whose
// sigma varies a lot in time, channels with sigma far above (or far
// below) their neighbours, channels that the mask zaps most of the time,
// and -ignorechan channels.

#define FLATTEN_NSAMPLES     32     // pre-pass blocks spread through the obs
#define FLATTEN_MINSPECT     4096   // minimum spectra per pre-pass sample
#define FLATTEN_DEADFRAC     1e-3   // sigma < this * median sigma -> dead
#define FLATTEN_LOWFRAC      0.1    // sigma < this * local sigma -> low
#define FLATTEN_NSIG         5.0    // outlier threshold in robust sigmas
#define FLATTEN_MASKFRAC     0.5    // masked more than this -> zero weight
#define FLATTEN_MINTIMES     8      // need this many samples for the time test

static const char *flag_names[] = {
    "good", "dead", "masked", "ignored", "ratio", "variable", "high", "low"
};

static float fmedian_copy(const float *arr, int n, float *work)
// Median of n floats without changing arr (work must hold n floats).
{
    if (n <= 0)
        return 0.0;
    memcpy(work, arr, n * sizeof(float));
    return median(work, n);
}

static float robust_spread(const float *x, int n, float *work)
// 1.4826 * MAD of x (a robust standard deviation).
{
    float med;
    int ii;
    if (n < 2)
        return 0.0;
    med = fmedian_copy(x, n, work);
    for (ii = 0; ii < n; ii++)
        work[ii] = fabsf(x[ii] - med);
    return 1.4826 * median(work, n);
}

void flatten_alloc(struct spectra_info *s, int numtimes)
// Allocate s->flatten for numtimes mean estimates.
{
    struct flatten_info *fi;
    int numchan = s->num_channels;

    fi = (struct flatten_info *) calloc(1, sizeof(struct flatten_info));
    fi->numchan = numchan;
    fi->numtimes = numtimes;
    fi->specnums = gen_dvect(numtimes);
    fi->avgs = gen_fvect((long) numtimes * numchan);
    fi->slopes = gen_fvect((long) numtimes * numchan);
    fi->stds = gen_fvect(numchan);
    fi->wts = gen_fvect(numchan);
    fi->flags = (unsigned char *) calloc(numchan, 1);
    fi->curmean = gen_fvect(numchan);
    s->flatten = fi;
}

void flatten_free(struct spectra_info *s)
{
    struct flatten_info *fi = s->flatten;
    if (fi == NULL)
        return;
    vect_free(fi->specnums);
    vect_free(fi->avgs);
    vect_free(fi->slopes);
    vect_free(fi->stds);
    vect_free(fi->wts);
    free(fi->flags);
    vect_free(fi->curmean);
    free(fi);
    s->flatten = NULL;
}

void flatten_finish(struct spectra_info *s)
// Compute the interpolation slopes between the mean estimates.  Call
// after avgs and specnums are filled (also on MPI nodes after a
// broadcast).
{
    struct flatten_info *fi = s->flatten;
    int ii, jj, nc = fi->numchan;
    for (ii = 0; ii < fi->numtimes; ii++) {
        float *a0 = fi->avgs + (long) ii * nc, *sl = fi->slopes + (long) ii * nc;
        if (ii == fi->numtimes - 1) {
            for (jj = 0; jj < nc; jj++)
                sl[jj] = 0.0;
        } else {
            float *a1 = a0 + nc;
            double dspec = fi->specnums[ii + 1] - fi->specnums[ii];
            for (jj = 0; jj < nc; jj++)
                sl[jj] = (dspec > 0.0) ? (a1[jj] - a0[jj]) / dspec : 0.0;
        }
    }
}


/* ---------------------------------------------------------------- */
/* Gathering the statistics: a (numtimes x numchan) table of channel */
/* means and standard deviations, plus a "bad cell" table.           */
/* ---------------------------------------------------------------- */

typedef struct {
    int numtimes, numchan;
    long long ptsper;           // spectra behind each estimate
    double *specnums;           // centre spectrum of each estimate
    float *avg, *std;           // numtimes x numchan
    unsigned char *bad;         // numtimes x numchan (1 = don't use)
} chanstats;

static void free_chanstats(chanstats * cs)
{
    vect_free(cs->specnums);
    vect_free(cs->avg);
    vect_free(cs->std);
    free(cs->bad);
}

static void alloc_chanstats(chanstats * cs, int numtimes, int numchan)
{
    cs->numtimes = numtimes;
    cs->numchan = numchan;
    cs->specnums = gen_dvect(numtimes);
    cs->avg = gen_fvect((long) numtimes * numchan);
    cs->std = gen_fvect((long) numtimes * numchan);
    cs->bad = (unsigned char *) calloc((size_t) numtimes * numchan, 1);
}

static void mark_masked_cells(chanstats * cs, double dt, mask * obsmask)
// Flag the cells that the rfifind mask zaps.
{
    int ii, jj, nummasked, *maskchans;
    if (obsmask == NULL || obsmask->numchan != cs->numchan)
        return;
    maskchans = gen_ivect(cs->numchan);
    for (ii = 0; ii < cs->numtimes; ii++) {
        double t0 = (cs->specnums[ii] - 0.5 * cs->ptsper) * dt;
        nummasked = check_mask(t0, cs->ptsper * dt, obsmask, maskchans);
        if (nummasked == -1) {
            memset(cs->bad + (long) ii * cs->numchan, 1, cs->numchan);
        } else {
            for (jj = 0; jj < nummasked; jj++)
                cs->bad[(long) ii * cs->numchan + maskchans[jj]] = 1;
        }
    }
    vect_free(maskchans);
}

static int stats_from_rfifind(char *maskfilenm, struct spectra_info *s,
                              mask * obsmask, chanstats * cs)
// Read the rfifind .stats file that goes with maskfilenm.  Returns 1 on
// success, 0 if there is no usable file.
{
    FILE *statsfile;
    char *statsfilenm, *root, *suffix;
    int ii, jj, numchan, numint, ptsperint, lobin, numbetween;
    float *pow;

    if (maskfilenm == NULL)
        return 0;
    if (split_root_suffix(maskfilenm, &root, &suffix) == 0)
        return 0;
    statsfilenm = (char *) calloc(strlen(root) + 10, sizeof(char));
    sprintf(statsfilenm, "%s.stats", root);
    free(root);
    free(suffix);
    statsfile = fopen(statsfilenm, "rb");
    if (statsfile == NULL) {
        printf("   No rfifind statistics file '%s'.\n", statsfilenm);
        free(statsfilenm);
        return 0;
    }
    chkfread(&numchan, sizeof(int), 1, statsfile);
    chkfread(&numint, sizeof(int), 1, statsfile);
    chkfread(&ptsperint, sizeof(int), 1, statsfile);
    chkfread(&lobin, sizeof(int), 1, statsfile);
    chkfread(&numbetween, sizeof(int), 1, statsfile);
    if (numchan != s->num_channels || numint < 1 || ptsperint < 1) {
        printf("   '%s' has %d channels (the data have %d).  Not using it.\n",
               statsfilenm, numchan, s->num_channels);
        fclose(statsfile);
        free(statsfilenm);
        return 0;
    }
    alloc_chanstats(cs, numint, numchan);
    cs->ptsper = ptsperint;
    pow = gen_fvect((long) numint * numchan);
    chkfread(pow, sizeof(float), (long) numint * numchan, statsfile);
    chkfread(cs->avg, sizeof(float), (long) numint * numchan, statsfile);
    chkfread(cs->std, sizeof(float), (long) numint * numchan, statsfile);
    vect_free(pow);
    fclose(statsfile);
    for (ii = 0; ii < numint; ii++) {
        int allzero = 1;
        float *row = cs->std + (long) ii * numchan;
        cs->specnums[ii] = (ii + 0.5) * ptsperint;
        // rfifind writes all zeros for intervals that hold padding
        for (jj = 0; jj < numchan; jj++)
            if (row[jj] != 0.0) {
                allzero = 0;
                break;
            }
        if (allzero)
            memset(cs->bad + (long) ii * numchan, 1, numchan);
    }
    mark_masked_cells(cs, s->dt, obsmask);
    printf("   Read channel statistics for %d intervals from '%s'.\n",
           numint, statsfilenm);
    free(statsfilenm);
    return 1;
}

static int good_spectra(const float *data, int numspect, int numchan,
                        unsigned char *good, float *work)
// Mark the spectra whose zero-DM sum lies within 3 sigma of the block
// median: the same cut that clip_times() uses to pick the spectra that
// set its running channel averages.  Returns the number of good ones.
{
    int ii, jj, numgood = 0;
    double zavg, zvar;
    float zmed, zstd;
    float *zdm = work + numspect;
    for (ii = 0; ii < numspect; ii++) {
        const float *row = data + (long) ii * numchan;
        float sum = 0.0;
        for (jj = 0; jj < numchan; jj++)
            sum += row[jj];
        zdm[ii] = sum;
    }
    avg_var(zdm, numspect, &zavg, &zvar);
    zstd = sqrt(zvar);
    zmed = fmedian_copy(zdm, numspect, work);
    for (ii = 0; ii < numspect; ii++) {
        good[ii] = (fabsf(zdm[ii] - zmed) < 3.0 * zstd) ? 1 : 0;
        numgood += good[ii];
    }
    return numgood;
}

static int stats_from_data(struct spectra_info *s, mask * obsmask,
                           chanstats * cs)
// Pre-pass: read FLATTEN_NSAMPLES stretches of data spread through the
// observation and measure the channel means and standard deviations.
// Leaves the input positioned at spectrum 0.  Returns 1 on success.
{
    int ii, jj, kk, numchan = s->num_channels, spb = s->spectra_per_subint;
    int blocks_per_sample, numsamples, nfile_blocks, isamp, padding;
    long long totblocks = 0, numspect;
    float *rawdata, *work;
    unsigned char *good;
    double *sum, *sum2;
    long long *blockstart;

    blocks_per_sample = (FLATTEN_MINSPECT + spb - 1) / spb;
    numspect = (long long) blocks_per_sample * spb;
    // Only sample whole blocks that start on a block boundary of a file
    for (ii = 0; ii < s->num_files; ii++)
        totblocks += s->num_spec[ii] / spb;
    numsamples = totblocks / blocks_per_sample;
    if (numsamples > FLATTEN_NSAMPLES)
        numsamples = FLATTEN_NSAMPLES;
    if (numsamples < 1) {
        printf("   Too little data for the -flatten pre-pass.\n");
        return 0;
    }
    // Start spectrum of each sample, spread evenly over the files
    blockstart = (long long *) malloc(numsamples * sizeof(long long));
    for (isamp = 0; isamp < numsamples; isamp++) {
        long long target = (numsamples == 1) ? 0 :
            (long long) ((double) isamp / (numsamples - 1) *
                         (totblocks - blocks_per_sample));
        for (ii = 0; ii < s->num_files; ii++) {
            nfile_blocks = s->num_spec[ii] / spb;
            if (target < nfile_blocks || ii == s->num_files - 1) {
                if (target > nfile_blocks - blocks_per_sample)
                    target = nfile_blocks - blocks_per_sample;
                if (target < 0)
                    target = 0;
                blockstart[isamp] = s->start_spec[ii] + target * spb;
                break;
            }
            target -= nfile_blocks;
        }
    }
    alloc_chanstats(cs, numsamples, numchan);
    cs->ptsper = numspect;
    rawdata = gen_fvect(2L * spb * numchan + numspect * numchan);
    work = gen_fvect(2 * numspect);
    good = (unsigned char *) malloc(numspect);
    sum = gen_dvect(numchan);
    sum2 = gen_dvect(numchan);
    printf("   Measuring channel statistics from %d stretches of %lld spectra...\n",
           numsamples, numspect);
    for (isamp = 0; isamp < numsamples; isamp++) {
        float *block = rawdata + 2L * spb * numchan;
        int anypad = 0, numgood;
        offset_to_spectra(blockstart[isamp], s);
        for (kk = 0; kk < blocks_per_sample; kk++) {
            padding = 0;
            if (!s->get_rawblock(rawdata, s, &padding)) {
                anypad = 1;
                break;
            }
            if (padding)
                anypad = 1;
            memcpy(block + (long) kk * spb * numchan, rawdata,
                   (size_t) spb * numchan * sizeof(float));
        }
        cs->specnums[isamp] = blockstart[isamp] + 0.5 * numspect;
        if (anypad) {
            memset(cs->bad + (long) isamp * numchan, 1, numchan);
            continue;
        }
        numgood = good_spectra(block, numspect, numchan, good, work);
        if (numgood < numspect / 4) {
            memset(cs->bad + (long) isamp * numchan, 1, numchan);
            continue;
        }
        for (jj = 0; jj < numchan; jj++)
            sum[jj] = sum2[jj] = 0.0;
        for (ii = 0; ii < numspect; ii++) {
            const float *row = block + (long) ii * numchan;
            if (!good[ii])
                continue;
            for (jj = 0; jj < numchan; jj++) {
                sum[jj] += row[jj];
                sum2[jj] += (double) row[jj] * row[jj];
            }
        }
        for (jj = 0; jj < numchan; jj++) {
            double avg = sum[jj] / numgood;
            double var = sum2[jj] / numgood - avg * avg;
            cs->avg[(long) isamp * numchan + jj] = avg;
            cs->std[(long) isamp * numchan + jj] = (var > 0.0) ? sqrt(var) : 0.0;
        }
    }
    mark_masked_cells(cs, s->dt, obsmask);
    offset_to_spectra(0, s);
    free(blockstart);
    free(good);
    vect_free(rawdata);
    vect_free(work);
    vect_free(sum);
    vect_free(sum2);
    return 1;
}


/* ---------------------------------------------------------------- */
/* Deciding the weights                                              */
/* ---------------------------------------------------------------- */

static void choose_weights(chanstats * cs, struct spectra_info *s, int *use_ratio)
{
    struct flatten_info *fi = s->flatten;
    int ii, jj, nt = cs->numtimes, nc = cs->numchan, nlive, ngoodtimes;
    float *col, *work, *vals, *mu, *tv, *ls;
    float medsig, spread;

    col = gen_fvect(nt);
    work = gen_fvect(nt > nc ? nt : nc);
    vals = gen_fvect(nc);
    mu = gen_fvect(nc);
    tv = gen_fvect(nc);
    ls = gen_fvect(nc);

    // Each channel's median sigma and mean over its good intervals
    for (jj = 0; jj < nc; jj++) {
        int ngood = 0;
        for (ii = 0; ii < nt; ii++)
            if (!cs->bad[(long) ii * nc + jj])
                col[ngood++] = cs->std[(long) ii * nc + jj];
        fi->flags[jj] = FLATTEN_GOOD;
        if (ngood == 0 || ngood < (1.0 - FLATTEN_MASKFRAC) * nt) {
            fi->flags[jj] = FLATTEN_MASKED;
            fi->stds[jj] = (ngood) ? median(col, ngood) : 0.0;
            mu[jj] = 0.0;
            continue;
        }
        fi->stds[jj] = median(col, ngood);
        ngood = 0;
        for (ii = 0; ii < nt; ii++)
            if (!cs->bad[(long) ii * nc + jj])
                col[ngood++] = cs->avg[(long) ii * nc + jj];
        mu[jj] = median(col, ngood);
    }
    for (jj = 0; jj < s->num_ignorechans; jj++)
        if (s->ignorechans[jj] >= 0 && s->ignorechans[jj] < nc)
            fi->flags[s->ignorechans[jj]] = FLATTEN_IGNORED;

    // Dead channels: sigma tiny compared with the typical channel
    nlive = 0;
    for (jj = 0; jj < nc; jj++)
        if (fi->flags[jj] == FLATTEN_GOOD && isfinite(fi->stds[jj])
            && fi->stds[jj] > 0.0)
            vals[nlive++] = fi->stds[jj];
    medsig = (nlive) ? median(vals, nlive) : 0.0;
    for (jj = 0; jj < nc; jj++)
        if (fi->flags[jj] == FLATTEN_GOOD &&
            (!isfinite(fi->stds[jj]) || fi->stds[jj] <= FLATTEN_DEADFRAC * medsig))
            fi->flags[jj] = FLATTEN_DEAD;

    // Radiometer test: for total-power data sigma/mean is nearly the same
    // in every channel (1/sqrt(npol*df*dt)).  Use it only if the measured
    // ratio is close to that value, which tells us the means are real
    // powers (not 2-bit levels or offset-subtracted values).
    *use_ratio = 0;
    {
        int nq = 0, npos = 0, ngood = 0;
        double r1 = 1.0 / sqrt(fabs(s->df) * 1e6 * s->dt), qmed = 0.0;
        for (jj = 0; jj < nc; jj++) {
            if (fi->flags[jj] != FLATTEN_GOOD)
                continue;
            ngood++;
            if (mu[jj] > 0.0) {
                npos++;
                vals[nq++] = logf(fi->stds[jj] / mu[jj]);
            }
        }
        if (nq > 0)
            qmed = fmedian_copy(vals, nq, work);
        if (ngood > 0 && npos > 0.9 * ngood && r1 > 0.0 &&
            exp(qmed) / r1 > 0.35 && exp(qmed) / r1 < 2.0) {
            *use_ratio = 1;
            spread = robust_spread(vals, nq, work);
            if (spread < 0.05)
                spread = 0.05;
            for (jj = 0; jj < nc; jj++) {
                if (fi->flags[jj] != FLATTEN_GOOD)
                    continue;
                if (mu[jj] <= 0.0 ||
                    fabs(logf(fi->stds[jj] / mu[jj]) - qmed) > FLATTEN_NSIG * spread)
                    fi->flags[jj] = FLATTEN_RATIO;
            }
            printf("   Median sigma/mean = %.4g (radiometer value %.4g-%.4g): "
                   "using the sigma/mean test.\n", exp(qmed), r1 / sqrt(2.0), r1);
        } else if (nq > 0) {
            printf("   Median sigma/mean = %.4g is not the radiometer value "
                   "(%.4g-%.4g): skipping the sigma/mean test.\n",
                   exp(qmed), r1 / sqrt(2.0), r1);
        }
    }

    // Time test: a channel whose sigma (relative to the other channels at
    // the same time) wanders a lot holds intermittent RFI.
    ngoodtimes = 0;
    for (ii = 0; ii < nt; ii++) {
        int nb = 0;
        for (jj = 0; jj < nc; jj++)
            nb += cs->bad[(long) ii * nc + jj];
        if (nb < nc)
            ngoodtimes++;
    }
    if (ngoodtimes >= FLATTEN_MINTIMES) {
        float *rowmed = gen_fvect(nt);
        int ntv = 0;
        for (ii = 0; ii < nt; ii++) {
            int nv = 0;
            for (jj = 0; jj < nc; jj++)
                if (fi->flags[jj] == FLATTEN_GOOD && !cs->bad[(long) ii * nc + jj])
                    vals[nv++] = cs->std[(long) ii * nc + jj];
            rowmed[ii] = (nv) ? median(vals, nv) : 0.0;
        }
        for (jj = 0; jj < nc; jj++) {
            int ngood = 0;
            float rmed;
            tv[jj] = 0.0;
            if (fi->flags[jj] != FLATTEN_GOOD)
                continue;
            for (ii = 0; ii < nt; ii++)
                if (!cs->bad[(long) ii * nc + jj] && rowmed[ii] > 0.0)
                    col[ngood++] = cs->std[(long) ii * nc + jj] / rowmed[ii];
            if (ngood < FLATTEN_MINTIMES)
                continue;
            rmed = fmedian_copy(col, ngood, work);
            tv[jj] = (rmed > 0.0) ? robust_spread(col, ngood, work) / rmed : 0.0;
            vals[ntv++] = tv[jj];
        }
        if (ntv > 0) {
            // Clean channels sit at the statistical floor ~1/sqrt(2*ptsper)
            float thresh = 5.0 * median(vals, ntv);
            if (thresh < 0.03)
                thresh = 0.03;
            for (jj = 0; jj < nc; jj++)
                if (fi->flags[jj] == FLATTEN_GOOD && tv[jj] > thresh)
                    fi->flags[jj] = FLATTEN_VARIABLE;
        }
        vect_free(rowmed);
    }

    // Local test: log(sigma) against a running median of its good
    // neighbours.  Catches narrow RFI and roll-off channels; wide RFI
    // bands are left to the tests above and to the mask.
    {
        int win = nc / 64, half, nr = 0;
        float *loc = gen_fvect(nc), *r = gen_fvect(nc);
        if (win < 33)
            win = 33;
        half = win / 2;
        for (jj = 0; jj < nc; jj++)
            ls[jj] = (fi->flags[jj] == FLATTEN_GOOD) ? logf(fi->stds[jj]) : 0.0;
        for (jj = 0; jj < nc; jj++) {
            int lo = jj - half, hi = jj + half, nv = 0;
            if (lo < 0)
                lo = 0;
            if (hi > nc - 1)
                hi = nc - 1;
            for (ii = lo; ii <= hi; ii++)
                if (fi->flags[ii] == FLATTEN_GOOD)
                    vals[nv++] = ls[ii];
            loc[jj] = (nv) ? median(vals, nv) : ls[jj];
            if (fi->flags[jj] == FLATTEN_GOOD)
                r[nr++] = ls[jj] - loc[jj];
        }
        spread = robust_spread(r, nr, work);
        if (spread < 0.1)
            spread = 0.1;
        for (jj = 0; jj < nc; jj++) {
            float rr;
            if (fi->flags[jj] != FLATTEN_GOOD)
                continue;
            rr = ls[jj] - loc[jj];
            if (rr > FLATTEN_NSIG * spread)
                fi->flags[jj] = FLATTEN_HIGH;
            else if (rr < logf(FLATTEN_LOWFRAC))
                fi->flags[jj] = FLATTEN_LOW;
        }
        vect_free(loc);
        vect_free(r);
    }

    for (jj = 0; jj < nc; jj++)
        fi->wts[jj] = (fi->flags[jj] == FLATTEN_GOOD) ? 1.0 / fi->stds[jj] : 0.0;

    // The means to subtract.  Use each interval's mean unless it is
    // masked or far (> 3 sigma) from the channel's median mean; fill
    // those cells from the nearest good ones in time.
    for (jj = 0; jj < nc; jj++) {
        int last = -1;
        for (ii = 0; ii < nt; ii++) {
            long idx = (long) ii * nc + jj;
            int ok = !cs->bad[idx] && fabsf(cs->avg[idx] - mu[jj]) < 3.0 * fi->stds[jj];
            if (!ok)
                continue;
            fi->avgs[idx] = cs->avg[idx];
            if (last < ii - 1) {        // fill the gap since the last good one
                int kk;
                for (kk = last + 1; kk < ii; kk++) {
                    long kidx = (long) kk * nc + jj;
                    if (last < 0)
                        fi->avgs[kidx] = cs->avg[idx];
                    else {
                        float f = (float) (kk - last) / (ii - last);
                        fi->avgs[kidx] = (1.0 - f) * cs->avg[(long) last * nc + jj] +
                            f * cs->avg[idx];
                    }
                }
            }
            last = ii;
        }
        if (last < 0) {
            for (ii = 0; ii < nt; ii++)
                fi->avgs[(long) ii * nc + jj] = mu[jj];
        } else {
            for (ii = last + 1; ii < nt; ii++)
                fi->avgs[(long) ii * nc + jj] = cs->avg[(long) last * nc + jj];
        }
    }
    for (ii = 0; ii < nt; ii++)
        fi->specnums[ii] = cs->specnums[ii];

    // Padding (from the readers, before flattening) should land near
    // zero after flattening, so set unset padding values to the means.
    {
        int allzero = 1;
        for (jj = 0; jj < nc; jj++)
            if (s->padvals[jj] != 0.0)
                allzero = 0;
        if (allzero)
            for (jj = 0; jj < nc; jj++)
                s->padvals[jj] = mu[jj];
    }

    vect_free(col);
    vect_free(work);
    vect_free(vals);
    vect_free(mu);
    vect_free(tv);
    vect_free(ls);
}

static void report_and_write(struct spectra_info *s, char *outroot, int from_stats,
                             int use_ratio)
{
    struct flatten_info *fi = s->flatten;
    int ii, jj, counts[FLATTEN_NUMFLAGS], nzero;
    char *filenm;
    FILE *outfile;

    for (ii = 0; ii < FLATTEN_NUMFLAGS; ii++)
        counts[ii] = 0;
    for (jj = 0; jj < fi->numchan; jj++)
        counts[fi->flags[jj]]++;
    nzero = fi->numchan - counts[FLATTEN_GOOD];
    printf("   Zero-weighted %d of %d channels (%.1f%%):", nzero, fi->numchan,
           100.0 * nzero / fi->numchan);
    for (ii = 1; ii < FLATTEN_NUMFLAGS; ii++)
        if (counts[ii])
            printf(" %d %s", counts[ii], flag_names[ii]);
    printf("\n");
    printf("   Note:  -flatten weights each channel by 1/sigma, which equalises\n"
           "   the noise across the band.  That is S/N-optimal only if Tsys is\n"
           "   flat once the gain is removed (the data cannot tell gain from Tsys).\n"
           "   Channel weights are fixed for the whole observation; the channel\n"
           "   means follow slow drifts.\n");
    if (outroot == NULL)
        return;
    filenm = (char *) calloc(strlen(outroot) + 10, sizeof(char));
    sprintf(filenm, "%s.flatten", outroot);
    outfile = fopen(filenm, "w");
    if (outfile == NULL) {
        printf("   Warning:  could not write '%s'.\n", filenm);
        free(filenm);
        return;
    }
    fprintf(outfile, "# Bandpass flattening (-flatten) from %s\n",
            from_stats ? "the rfifind .stats file" : "a pre-pass through the data");
    fprintf(outfile, "# sigma/mean (radiometer) test: %s\n",
            use_ratio ? "used" : "not used");
    fprintf(outfile, "# Zero-weighted %d of %d channels\n", nzero, fi->numchan);
    fprintf(outfile, "# Zero-weighted channels (for -ignorechan): ");
    {
        int first = 1, start = -1;
        for (jj = 0; jj <= fi->numchan; jj++) {
            int bad = (jj < fi->numchan) && fi->flags[jj] != FLATTEN_GOOD;
            if (bad && start < 0)
                start = jj;
            if (!bad && start >= 0) {
                fprintf(outfile, "%s%d", first ? "" : ",", start);
                if (jj - 1 > start)
                    fprintf(outfile, ":%d", jj - 1);
                first = 0;
                start = -1;
            }
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "#%6s %14s %14s %14s %14s  %s\n", "chan", "freq(MHz)",
            "mean", "sigma", "weight", "flag");
    for (jj = 0; jj < fi->numchan; jj++) {
        double freq = s->lo_freq + jj * s->df;
        float meanval = 0.0;
        for (ii = 0; ii < fi->numtimes; ii++)
            meanval += fi->avgs[(long) ii * fi->numchan + jj];
        meanval /= fi->numtimes;
        fprintf(outfile, "%7d %14.6f %14.6g %14.6g %14.6g  %s\n", jj, freq,
                meanval, fi->stds[jj], fi->wts[jj], flag_names[fi->flags[jj]]);
    }
    fclose(outfile);
    printf("   Wrote the channel statistics and weights to '%s'.\n", filenm);
    free(filenm);
}


int flatten_init(struct spectra_info *s, char *maskfilenm, mask * obsmask,
                 char *outroot)
// Set up bandpass flattening for raw data in s.  Uses the rfifind .stats
// file that goes with maskfilenm (may be NULL) if there is one, and a
// pre-pass through the data if not.  Writes '<outroot>.flatten' (if
// outroot is not NULL).  Call after the mask and padding values are set
// and before any offset into the data.  Returns 1 if it used the stats
// file, 0 if it used the pre-pass.
{
    chanstats cs;
    int from_stats, use_ratio, zerodm = s->remove_zerodm;

    printf("\nBandpass flattening (-flatten):\n");
    // Zero-DM is done after flattening (on the flattened data), and the
    // pre-pass must see the data without it.
    s->remove_zerodm = 0;
    from_stats = stats_from_rfifind(maskfilenm, s, obsmask, &cs);
    if (!from_stats) {
        if (!stats_from_data(s, obsmask, &cs)) {
            printf("   Can't measure the channel statistics.  Exiting.\n\n");
            exit(1);
        }
    }
    flatten_alloc(s, cs.numtimes);
    s->flatten->zerodm = zerodm;
    choose_weights(&cs, s, &use_ratio);
    flatten_finish(s);
    free_chanstats(&cs);
    report_and_write(s, outroot, from_stats, use_ratio);
    printf("\n");
    return from_stats;
}


void flatten_apply(float *data, int numspect, long long startspec,
                   struct spectra_info *s)
// Flatten numspect spectra (time-major, num_channels per spectrum) that
// start at spectrum startspec: x -> (x - mean(t)) * weight.  Then
// remove the zero-DM signal (mean over the weighted channels) if -zerodm
// was requested.
{
    struct flatten_info *fi = s->flatten;
    int ii, jj, nc = fi->numchan, seg = 0;

    for (ii = 0; ii < numspect; ii++) {
        double spec = (double) startspec + ii;
        float *row = data + (long) ii * nc;
        const float *a, *sl;
        float ds;
        // Find the segment of the mean interpolation we are in
        while (seg < fi->numtimes - 1 && spec >= fi->specnums[seg + 1])
            seg++;
        a = fi->avgs + (long) seg * nc;
        sl = fi->slopes + (long) seg * nc;
        ds = spec - fi->specnums[seg];
        if (ds < 0.0 || seg == fi->numtimes - 1)
            ds = 0.0;
        for (jj = 0; jj < nc; jj++)
            row[jj] = (row[jj] - (a[jj] + ds * sl[jj])) * fi->wts[jj];
        if (fi->zerodm) {
            float sum = 0.0;
            int ngood = 0;
            for (jj = 0; jj < nc; jj++)
                if (fi->wts[jj] > 0.0) {
                    sum += row[jj];
                    ngood++;
                }
            if (ngood) {
                sum /= ngood;
                for (jj = 0; jj < nc; jj++)
                    if (fi->wts[jj] > 0.0)
                        row[jj] -= sum;
            }
        }
    }
}
