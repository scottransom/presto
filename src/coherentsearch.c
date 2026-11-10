/*
 * coherentsearch: a coherent harmonic-summing pulsar search of normalised
 * PRESTO .fft files.  See include/coherentsearch.h and docs/coherentsearch.1.
 */
#include <errno.h>
#include <limits.h>
#include <sys/stat.h>
#include <time.h>
#include "coherentsearch.h"
#include "coherentsearch_cmd.h"

/* <input without .fft><suffix>, beside the input or in outdir. */
static char *output_name(const char *infile, const char *outdir, const char *suffix)
{
    const char *base = infile;
    if (outdir) {
        const char *slash = strrchr(infile, '/');
        base = slash ? slash + 1 : infile;
    }
    size_t len = strlen(base);
    if (len > 4 && strcmp(base + len - 4, ".fft") == 0)
        len -= 4;
    char *out = malloc((outdir ? strlen(outdir) + 1 : 0) + len + strlen(suffix) + 1);
    out[0] = '\0';
    if (outdir) {
        strcpy(out, outdir);
        strcat(out, "/");
    }
    strncat(out, base, len);
    strcat(out, suffix);
    return out;
}

/* The candidate file for infile: -o if given, else <input>.cohout (or
 * <input>_deep.cohout with -deep). */
static char *cohout_name(const Cmdline *cmd, const char *infile)
{
    return cmd->outfileP ? strdup(cmd->outfile)
                         : output_name(infile, cmd->outdirP ? cmd->outdir : NULL,
                                       cmd->deepP ? "_deep.cohout" : ".cohout");
}

/* Whether option opt was given on the command line (the parser marks options
 * with defaults as present either way). */
static int option_given(int argc, char *argv[], const char *opt)
{
    for (int i = 1; i < argc && strcmp(argv[i], "--") != 0; i++)
        if (strcmp(argv[i], opt) == 0)
            return 1;
    return 0;
}

/* Plot l to <candidate file stem>_cohout.ps (and _NN.png, given gs). */
static void plot(const Cmdline *cmd, const cs_fft *fft, const cs_candlist *l,
                 const cs_params *p, const char *infile)
{
    char *cohnm = cohout_name(cmd, infile);
    size_t len = strlen(cohnm);
    if (len > 7 && strcmp(cohnm + len - 7, ".cohout") == 0)
        cohnm[len - 7] = '\0';
    char *stem = malloc(strlen(cohnm) + 16);
    sprintf(stem, "%s_cohout", cohnm);
    const int ret = cs_plot_candidates(fft, l, p->nharms, p->threshold, infile, stem);
    if (ret >= 0 && l->n > 0)
        printf("\nPlotted %ld candidates to '%s.ps'%s\n", l->n, stem,
               ret ? " (and _NN.png)" : "");
    free(stem);
    free(cohnm);
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
    if (cmd->deepP) {
        /* The deep pass's defaults; explicit options still win. */
        if (!option_given(argc, argv, "-nharms"))
            cmd->nharms = CS_DEEP_NHARMS;
        if (!option_given(argc, argv, "-maxdecim"))
            cmd->maxdecim = CS_DEEP_MAXDECIM;
        if (!option_given(argc, argv, "-hifreq"))
            cmd->hifreq = CS_DEEP_HIFREQ;
    }

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

    /* Plans and workspaces do not depend on the file, so build them once.
     * FFTW wisdom (from makewisdom) makes the small batched FFTs ~1.8x faster. */
    read_wisdom();
    cs_searcher *srch = cs_searcher_new(&p, cmd->blocksize, cmd->ncpus);

    printf("\n    Coherent Harmonic Summing Pulsar Search Routine\n");
    printf("                by Scott M. Ransom\n");

    for (int fi = 0; fi < cmd->argc; fi++) {
        const char *infile = cmd->argv[fi];
        struct timespec t0;
        cs_fft fft;
        /* -rewrite alone needs only the .inf (T and the position). */
        cs_read_fft_info(infile, &fft, !cmd->rewriteP || cmd->replotP);
        printf("\nAnalyzing %s data from '%s'.\n", fft.object[0] ? fft.object : "unknown",
               infile);
        printf("  Data are of duration %.3f seconds.\n\n", fft.T);
        if (cmd->replotP || cmd->rewriteP) {
            char *cohnm = cohout_name(cmd, infile);
            if (cmd->rewriteP) {
                int nknown;
                long n = cs_relabel_cohout(cohnm, &fft, &nknown);
                printf("Rewrote %ld candidates in '%s'\n", n, cohnm);
                printf("  %d match known pulsars\n", nknown);
            }
            cs_candlist cands = { 0 };
            if (cmd->replotP) {
                cs_read_cohout(cohnm, fft.T, cmd->ncands, &cands);
                const int nknown = cs_identify_known(&fft, &cands);
                printf("%sReplotting %ld candidates from '%s'\n", cmd->rewriteP ? "\n" : "",
                       cands.n, cohnm);
                printf("  %d match known pulsars\n", nknown);
                plot(cmd, &fft, &cands, &p, infile);
            }
            free(cohnm);
            cs_candlist_free(&cands);
            cs_free_fft(&fft);
            continue;
        }
        if (!p.sigma_measured) {
            fflush(stdout);
            cs_check_normalisation(&fft, infile);
        }
        cs_band band = cs_make_band(&fft, &p, cmd->lofreq, cmd->hifreq,
                                    cmd->lobinP ? cmd->lobin : -1);
        const int nprof = cmd->blocksize;
        printf("Using:\n");
        printf("  %d CPU core(s)\n", srch->nthreads);
        printf("  S/N threshold of %g\n", p.threshold);
        printf("  %d harmonics\n", p.nharms);
        printf("  f = %.3f to %.3f Hz for fundamental\n", band.r_lo / fft.T, band.r_hi / fft.T);
        printf("  k = 1..%d decimation factors\n", p.decims[p.ndecim - 1]);
        if (p.sigma_measured)
            printf("  Measuring chunk sigmas\n");
        if (cmd->deepP)
            printf("  Deep pass (-deep)\n");
        cs_hplan *plans = cs_build_plans(&p, band.r_lo);
        if (cmd->verboseP) {
            long long nchunks;
            long long *starts = cs_chunk_starts(plans, &fft, &p, band.total, nprof, &nchunks);
            printf("  FFTs planned %s\n", srch->ws[0]->wisdom ? "from FFTW wisdom"
                   : "with FFTW_ESTIMATE (no matching wisdom: run makewisdom)");
            printf("  %lld trials in %lld chunks of <= %d over r = %.1f .. %.1f bins\n",
                   band.total, nchunks, nprof, band.r_lo, band.r_hi);
            if (plans[0].P == 0)
                printf("  NOTE: -hidr %g gives no small rational trial step; "
                       "weights are computed per trial (slower)\n", p.hidr);
            else
                printf("  trial step %lld/%lld bins; m = %d\n", plans[0].pnum, plans[0].q,
                       p.m);
            free(starts);
        }
        printf("\n");
        fflush(stdout);
        cs_free_plans(plans, p.nharms);

        clock_gettime(CLOCK_MONOTONIC, &t0);
        cs_candlist cands = { 0 };
        cs_search(srch, &fft, &band, p.threshold, !cmd->noprogressP, &cands);
        long nleft = cands.n;
        long ndups = -1, nharm = -1, nperiodic = -1, ndeepharm = -1;
        if (!cmd->noremoveP) {
            cs_remove_duplicates(&cands, cmd->drtol);
            ndups = nleft - cands.n;
            nleft = cands.n;
        }
        if (!cmd->noharmremoveP) {
            cs_remove_harmonics(&cands, cmd->numharm, 1.0);
            nharm = nleft - cands.n;
            nleft = cands.n;
        }
        if (!cmd->noperiodicP)
            nperiodic = cs_remove_periodic(&fft, &cands, CS_PERIODIC_NLOW,
                                           CS_PERIODIC_MINHARM, CS_PERIODIC_RATIO);
        /* A deep pass drops what the default search of the same file already
         * found more strongly (e.g. a bright pulsar's subharmonics, whose
         * fundamental is above the deep band). */
        char *refnm = NULL;
        if (cmd->deepP && !cmd->noharmremoveP) {
            refnm = output_name(infile, cmd->outdirP ? cmd->outdir : NULL, ".cohout");
            char *outnm = cohout_name(cmd, infile);
            if (strcmp(outnm, refnm) != 0 && access(refnm, R_OK) == 0) {
                cs_candlist ref = { 0 };
                cs_read_cohout(refnm, fft.T, LONG_MAX, &ref);
                ndeepharm = cs_remove_harmonics_of(&cands, &ref, cmd->numharm, 1.0);
                cs_candlist_free(&ref);
            }
            free(outnm);
        }
        cs_sort_by_metric(&cands);
        const long nbeyond = cands.n > cmd->ncands ? cands.n - cmd->ncands : 0;
        if (nbeyond)
            cands.n = cmd->ncands;
        cs_measure_ducy(&fft, &p, &cands);
        const int nknown = cs_identify_known(&fft, &cands);

        char *outnm = cohout_name(cmd, infile);
        FILE *out = chkfopen(outnm, "w");
        cs_write_cohout(&cands, out);
        fclose(out);

        printf("Search took %.2f seconds\n\n", seconds_since(&t0));
        printf("Wrote %ld candidates to '%s'\n\n", cands.n, outnm);
        if (ndups >= 0)
            printf("  %ld removed as duplicates\n", ndups);
        if (nharm >= 0)
            printf("  %ld removed as harmonics\n", nharm);
        if (refnm) {
            if (ndeepharm >= 0)
                printf("  %ld removed as harmonics of stronger candidates in '%s'\n",
                       ndeepharm, refnm);
            else
                printf("  NOTE: no default search ('%s') to remove harmonics of "
                       "(run one before -deep)\n", refnm);
            free(refnm);
        }
        if (nperiodic >= 0)
            printf("  %ld removed as periodic\n", nperiodic);
        if (nbeyond)
            printf("  %ld beyond the %d strongest (-ncands)\n", nbeyond, cmd->ncands);
        printf("  %d match known pulsars\n", nknown);
        if (cmd->plotP)
            plot(cmd, &fft, &cands, &p, infile);
        fflush(stdout);
        free(outnm);
        cs_candlist_free(&cands);
        cs_free_fft(&fft);
    }
    cs_searcher_free(srch);
    return 0;
}
