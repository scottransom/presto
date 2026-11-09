/*
 * cs_dump: run the C coherentsearch stages on one chunk and write them in the
 * layout of tests/coherentsearch/fixtures/<case>/ (see README.md), so that
 * tests/test_coherentsearch.py can compare them with CoherentSearch.jl's.
 * Built but not installed.
 *
 *   cs_dump file.fft outdir t0 n [lofreq [maxdecim]]
 *
 * Uses the coherent_search defaults otherwise; the global trial grid starts at
 * r_lo = lofreq * T (default lofreq 0.1).
 */
#include <errno.h>
#include <sys/stat.h>
#include "coherentsearch.h"

static void writebin(const char *dir, const char *name, const void *data, size_t nbytes)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    FILE *f = chkfopen(path, "wb");
    chkfwrite((void *) data, 1, nbytes, f);
    fclose(f);
}

static void print_intlist(FILE *f, const char *key, const int *v, int n)
{
    fprintf(f, "%s = [", key);
    for (int i = 0; i < n; i++)
        fprintf(f, "%s%d", i ? ", " : "", v[i]);
    fprintf(f, "]\n");
}

int main(int argc, char *argv[])
{
    if (argc < 5) {
        fprintf(stderr, "usage: %s file.fft outdir t0 n [lofreq [maxdecim]]\n", argv[0]);
        return 1;
    }
    const char *outdir = argv[2];
    long long t0 = atoll(argv[3]);
    int n = atoi(argv[4]);
    double lofreq = argc > 5 ? atof(argv[5]) : 0.1;
    int maxdecim = argc > 6 ? atoi(argv[6]) : 6;
    if (mkdir(outdir, 0755) != 0 && errno != EEXIST) {
        perror(outdir);
        return 1;
    }

    cs_params p;
    cs_default_params(&p, maxdecim);
    cs_fft fft;
    cs_read_fft(argv[1], &fft);
    const double r_lo = lofreq * fft.T;
    const double rstart = r_lo + t0 * cs_lodr(&p);

    cs_hplan *plans = cs_build_plans(&p, r_lo);
    cs_workspace *ws = cs_workspace_new(&p, n);
    cs_fill_chunk(ws, plans, &fft, &p, t0, n);
    writebin(outdir, "ftprofs.c64", ws->ftprofs, sizeof(fcomplex) * (size_t) n * (p.nharms + 1));

    char path[4096];
    snprintf(path, sizeof(path), "%s/meta.toml", outdir);
    FILE *meta = chkfopen(path, "w");
    int filled[CS_MAXDECIM * 64], nf = 0;
    for (int h = 1; h <= p.nharms; h++)
        if (ws->filled[h])
            filled[nf++] = h;
    print_intlist(meta, "filled", filled, nf);
    fprintf(meta, "n = %d\nt0 = %lld\nr_lo = %.17g\nrstart = %.17g\nlodr = %.17g\n",
            n, t0, r_lo, rstart, cs_lodr(&p));

    double *mv = malloc(sizeof(double) * n);
    for (int i = 0; i < ws->ndecim; i++) {
        cs_decim *d = ws->d + i;
        cs_prepare_decim(ws, d, &fft, &p, rstart, n);
        char name[64];
        if (d->k == 1)
            writebin(outdir, "profs_k1.f32", d->profs, sizeof(float) * (size_t) n * d->nbins);
        fprintf(meta, "\n[k%d]\nHk = %d\nnbins = %d\nnvalid = %d\n", d->k, d->Hk, d->nbins,
                d->nvalid);
        fprintf(meta, "sigma_analytic = %.17g\nsigma_measured = %.9g\n", d->sigma_analytic,
                (double) d->sigma_measured);
        print_intlist(meta, "widths", d->widths, d->nwidths);
        fprintf(meta, "shape = [");
        for (int w = 0; w < d->nwidths; w++)
            fprintf(meta, "%s%.17g", w ? ", " : "", d->shape[w]);
        fprintf(meta, "]\n");
        for (int s = 0; s < 2; s++) {
            cs_score_decim(d, s == 0 ? d->sigma_analytic : d->sigma_measured);
            for (int j = 0; j < d->nvalid; j++)
                mv[j] = d->mvals[j];
            snprintf(name, sizeof(name), "metric_k%d_%s.f64", d->k, s == 0 ? "analytic" : "measured");
            writebin(outdir, name, mv, sizeof(double) * d->nvalid);
        }
    }
    fclose(meta);
    free(mv);
    cs_workspace_free(ws);
    cs_free_plans(plans, p.nharms);
    cs_free_fft(&fft);
    return 0;
}
