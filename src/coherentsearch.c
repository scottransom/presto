/*
 * coherentsearch: a coherent harmonic-summing pulsar search of normalised
 * PRESTO .fft files.  See include/coherentsearch.h and docs/coherentsearch.1.
 */
#include <errno.h>
#include <sys/stat.h>
#include <time.h>
#include "coherentsearch.h"
#include "coherentsearch_cmd.h"

/* <input without .fft>.cohout, beside the input or in outdir. */
static char *output_name(const char *infile, const char *outdir)
{
    const char *base = infile;
    if (outdir) {
        const char *slash = strrchr(infile, '/');
        base = slash ? slash + 1 : infile;
    }
    size_t len = strlen(base);
    if (len > 4 && strcmp(base + len - 4, ".fft") == 0)
        len -= 4;
    char *out = malloc((outdir ? strlen(outdir) + 1 : 0) + len + 8);
    out[0] = '\0';
    if (outdir) {
        strcpy(out, outdir);
        strcat(out, "/");
    }
    strncat(out, base, len);
    strcat(out, ".cohout");
    return out;
}

static double seconds_since(const struct timespec *t0)
{
    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    return (t1.tv_sec - t0->tv_sec) + 1e-9 * (t1.tv_nsec - t0->tv_nsec);
}

int main(int argc, char *argv[])
{
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }
    Cmdline *cmd = parseCmdline(argc, argv);

    if (cmd->m % 2) {
        fprintf(stderr, "coherentsearch: -m must be even (got %d)\n", cmd->m);
        exit(1);
    }
    if (cmd->argc > 1 && cmd->outfileP) {
        fprintf(stderr, "coherentsearch: -o names one output file but %d input files "
                "were given; omit it or use -outdir\n", cmd->argc);
        exit(1);
    }
    if (cmd->outdirP && mkdir(cmd->outdir, 0755) != 0 && errno != EEXIST) {
        perror(cmd->outdir);
        exit(1);
    }

    cs_params p;
    cs_default_params(&p, cmd->maxdecim);
    p.nharms = cmd->nharms;
    p.m = cmd->m;
    p.hidr = cmd->hidr;
    p.threshold = cmd->threshold;
    p.sigma_measured = cmd->measuredsigmaP;
    cs_set_decimations(&p, cmd->maxdecim);

    for (int fi = 0; fi < cmd->argc; fi++) {
        const char *infile = cmd->argv[fi];
        struct timespec t0;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        cs_fft fft;
        cs_read_fft(infile, &fft);
        cs_band band = cs_make_band(&fft, &p, cmd->lofreq, cmd->hifreq,
                                    cmd->lobinP ? cmd->lobin : -1);
        const int nprof = cmd->blocksize;
        fprintf(stderr, "coherentsearch: %s (T = %.3f s), %d harmonics, k = 1..%d, "
                "%.4f-%.4f Hz\n", infile, fft.T, p.nharms, p.decims[p.ndecim - 1],
                band.r_lo / fft.T, band.r_hi / fft.T);

        cs_hplan *plans = cs_build_plans(&p, band.r_lo);
        if (cmd->verboseP) {
            long long nchunks;
            long long *starts = cs_chunk_starts(plans, &fft, &p, band.total, nprof, &nchunks);
            fprintf(stderr, "  %lld trials in %lld chunks of <= %d over r = %.1f .. %.1f bins\n",
                    band.total, nchunks, nprof, band.r_lo, band.r_hi);
            if (plans[0].P == 0)
                fprintf(stderr, "  NOTE: -hidr %g gives no small rational trial step; "
                        "weights are computed per trial (slower)\n", p.hidr);
            else
                fprintf(stderr, "  trial step %lld/%lld bins; m = %d; %s noise scale\n",
                        plans[0].pnum, plans[0].q, p.m,
                        p.sigma_measured ? "measured" : "analytic");
            free(starts);
        }
        if (!p.sigma_measured && band.total > 0) {
            cs_workspace *ws = cs_workspace_new(&p, nprof);
            cs_sigma_sanity_check(ws, plans, &fft, &p, &band, nprof, infile);
            cs_workspace_free(ws);
        }
        cs_free_plans(plans, p.nharms);

        cs_candlist cands = { 0 };
        cs_search(&fft, &p, &band, nprof, p.threshold, !cmd->noprogressP, &cands);
        const long nraw = cands.n;
        if (!cmd->noremoveP)
            cs_remove_duplicates(&cands, cmd->drtol);
        if (!cmd->noharmremoveP)
            cs_remove_harmonics(&cands, cmd->numharm, 1.0);
        cs_sort_by_metric(&cands);
        if (cands.n > cmd->ncands)
            cands.n = cmd->ncands;
        cs_measure_ducy(&fft, &p, &cands);

        char *outnm = cmd->outfileP ? strdup(cmd->outfile)
                                    : output_name(infile, cmd->outdirP ? cmd->outdir : NULL);
        FILE *out = chkfopen(outnm, "w");
        cs_write_cohout(&cands, out);
        fclose(out);
        fprintf(stderr, "  %ld trials above S/N %g -> %ld candidates in %s (%.1f s)\n",
                nraw, p.threshold, cands.n, outnm, seconds_since(&t0));
        free(outnm);
        cs_candlist_free(&cands);
        cs_free_fft(&fft);
    }
    return 0;
}
