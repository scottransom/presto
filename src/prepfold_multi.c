#include <ctype.h>
#include <unistd.h>
#include "prepfold.h"
#include "prepfold_multi_cmd.h"
#include "mask.h"
#include "backend_common.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

/* Fold-memory safety net.  prepfold_multi keeps every candidate's folds        */
/* resident at once (Σ_c nsub·npart·proflen doubles) plus per-thread dedispersion */
/* scratch, so a careless run with thousands of candidates could exhaust RAM.    */
/* We estimate the total up front, print it, warn past a soft threshold, and     */
/* fail gracefully (naming the offending candidate) past a hard cap -- rather    */
/* than letting gen_dvect() abort mid-allocation or the OOM killer step in.      */
#define BYTES_PER_GB              (1024.0 * 1024.0 * 1024.0)
#define FOLD_MEM_WARN_GB_FALLBACK 4.0
#define FOLD_MEM_CAP_GB_FALLBACK  64.0

/* malloc that aborts with a clear message instead of returning NULL, so an     */
/* allocation failure surfaces as a diagnostic here rather than a later NULL     */
/* dereference.  The large fold arrays go through gen_dvect/gen_fvect, which      */
/* already abort on failure; this covers the small bookkeeping allocations.      */
static void *chkmalloc(size_t n, const char *what)
{
    void *p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "\nprepfold_multi: out of memory allocating %s.\n\n", what);
        exit(1);
    }
    return p;
}

/* Resolve the soft-warn and hard-cap fold-memory thresholds from the host's    */
/* physical RAM so the guard scales with the machine instead of a fixed guess:  */
/* cap at total RAM (a run needing more cannot possibly fit), warn at half.  If  */
/* sysconf() cannot report RAM, fall back to the fixed GB thresholds above.      */
static void fold_mem_limits(size_t *warn, size_t *cap)
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long pgsize = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && pgsize > 0) {
        size_t ram = (size_t) pages * (size_t) pgsize;
        *cap = ram;
        *warn = ram / 2;
    } else {
        *cap = (size_t) (FOLD_MEM_CAP_GB_FALLBACK * BYTES_PER_GB);
        *warn = (size_t) (FOLD_MEM_WARN_GB_FALLBACK * BYTES_PER_GB);
    }
}

/* set_posn lives in prepfold_utils.c (no public header); declare it here as   */
/* prepfold.c does.                                                            */
void set_posn(prepfoldinfo * in, infodata * idata);

/*
 * prepfold_multi: fold many candidates -- each with its own period/frequency,
 * derivatives, and DM -- in a single pass over the raw data.  The expensive
 * resource (raw-data disk I/O) is paid once: each raw block is read and cleaned
 * once via the shared rawblock_reader, then dedispersed and folded per
 * candidate.  All of the fold math is reused from prepfold_pipeline.c (the same
 * code prepfold uses); this driver only adds the candidate-file parser and the
 * read-once/dedisperse-many control flow.
 */

/* ------------------------------------------------------------------------- */
/* Candidate file parsing                                                    */
/* ------------------------------------------------------------------------- */

#define CANDLINE_LEN 2048
#define MAX_CANDCOLS 8

/* One candidate row.  x0/x1/x2 hold p,pd,pdd (period mode) or f,fd,fdd        */
/* (frequency mode); which is decided file-wide by the header.                 */
typedef struct {
    char *name;                 /* NULL => autogenerate like prepfold */
    double dm;
    double x0, x1, x2;
} multicand;

enum { COL_NAME = 1, COL_DM, COL_X0, COL_X1, COL_X2 };

static void candfile_error(const char *fn, int lineno, const char *msg, const char *arg)
{
    fprintf(stderr, "\nprepfold_multi: error in candidate file '%s'", fn);
    if (lineno > 0)
        fprintf(stderr, " (line %d)", lineno);
    fprintf(stderr, ":\n  %s", msg);
    if (arg)
        fprintf(stderr, " '%s'", arg);
    fprintf(stderr, "\n\n");
    exit(1);
}

static char *dupstr(const char *s)
{
    char *d = (char *) chkmalloc(strlen(s) + 1, "candidate string");
    strcpy(d, s);
    return d;
}

static void lowercase(char *s)
{
    for (; *s; s++)
        *s = (char) tolower((unsigned char) *s);
}

static double parse_candval(const char *fn, int lineno, const char *col, char *v)
{
    char *end;
    double x = strtod(v, &end);
    if (end == v || *end != '\0') {
        char msg[80];
        snprintf(msg, sizeof(msg), "expected a number in column '%s', got", col);
        candfile_error(fn, lineno, msg, v);
    }
    return x;
}

/* Parse the candidate file.  Returns a malloc'd array of `ncand` candidates    */
/* and sets *is_period_out to 1 for a period (p/pd/pdd) file, 0 for frequency.  */
static multicand *parse_candfile(const char *fn, int *ncand_out, int *is_period_out)
{
    FILE *f = chkfopen((char *) fn, "r");
    char line[CANDLINE_LEN];
    int col_field[MAX_CANDCOLS];
    const char *col_name[MAX_CANDCOLS];
    int ncols = 0, lineno = 0, header_done = 0, is_period = 0;
    int seen_p = 0, seen_f = 0, seen_pd = 0, seen_fd = 0, seen_pdd = 0, seen_fdd = 0;
    int seen_dm = 0;
    multicand *cands = NULL;
    int ncand = 0, cap = 0;

    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        lineno++;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p == '\0' || *p == '\n' || *p == '\r')
            continue;           /* blank line */

        if (!header_done) {     /* the first non-blank line is the header */
            char *tok;
            if (*p != '#')
                candfile_error(fn, lineno,
                               "the first non-blank line must be a header starting with '#'",
                               NULL);
            p++;                /* skip the leading '#' */
            for (tok = strtok(p, " \t\n\r"); tok; tok = strtok(NULL, " \t\n\r")) {
                if (ncols >= MAX_CANDCOLS)
                    candfile_error(fn, lineno, "too many header columns", NULL);
                lowercase(tok);
                if (!strcmp(tok, "name")) {
                    col_field[ncols] = COL_NAME; col_name[ncols] = "name";
                } else if (!strcmp(tok, "dm")) {
                    col_field[ncols] = COL_DM;   col_name[ncols] = "DM"; seen_dm = 1;
                } else if (!strcmp(tok, "p")) {
                    col_field[ncols] = COL_X0;   col_name[ncols] = "p"; seen_p = 1;
                } else if (!strcmp(tok, "f")) {
                    col_field[ncols] = COL_X0;   col_name[ncols] = "f"; seen_f = 1;
                } else if (!strcmp(tok, "pd")) {
                    col_field[ncols] = COL_X1;   col_name[ncols] = "pd"; seen_pd = 1;
                } else if (!strcmp(tok, "fd")) {
                    col_field[ncols] = COL_X1;   col_name[ncols] = "fd"; seen_fd = 1;
                } else if (!strcmp(tok, "pdd")) {
                    col_field[ncols] = COL_X2;   col_name[ncols] = "pdd"; seen_pdd = 1;
                } else if (!strcmp(tok, "fdd")) {
                    col_field[ncols] = COL_X2;   col_name[ncols] = "fdd"; seen_fdd = 1;
                } else {
                    candfile_error(fn, lineno, "unknown column key", tok);
                }
                ncols++;
            }
            if (!seen_dm)
                candfile_error(fn, lineno, "header is missing the required 'DM' column", NULL);
            if (seen_p == seen_f)
                candfile_error(fn, lineno,
                               "header must contain exactly one of 'p' or 'f'", NULL);
            is_period = seen_p;
            if (is_period && (seen_fd || seen_fdd))
                candfile_error(fn, lineno,
                               "cannot mix period 'p' with frequency derivatives 'fd'/'fdd'",
                               NULL);
            if (!is_period && (seen_pd || seen_pdd))
                candfile_error(fn, lineno,
                               "cannot mix frequency 'f' with period derivatives 'pd'/'pdd'",
                               NULL);
            header_done = 1;
            continue;
        }

        if (*p == '#')
            continue;           /* comment line in the body */

        {                       /* a data row */
            char *toks[MAX_CANDCOLS + 1];
            char *tok;
            int nt = 0, ci;
            multicand *cd;

            for (tok = strtok(p, " \t\n\r"); tok && nt <= MAX_CANDCOLS;
                 tok = strtok(NULL, " \t\n\r"))
                toks[nt++] = tok;
            if (nt != ncols) {
                char msg[96];
                snprintf(msg, sizeof(msg),
                         "found %d columns but the header defines %d", nt, ncols);
                candfile_error(fn, lineno, msg, NULL);
            }

            if (ncand == cap) {
                cap = cap ? cap * 2 : 8;
                cands = (multicand *) realloc(cands, cap * sizeof(multicand));
            }
            cd = &cands[ncand];
            cd->name = NULL;
            cd->dm = cd->x0 = cd->x1 = cd->x2 = 0.0;
            for (ci = 0; ci < ncols; ci++) {
                switch (col_field[ci]) {
                case COL_NAME: cd->name = dupstr(toks[ci]); break;
                case COL_DM:   cd->dm  = parse_candval(fn, lineno, col_name[ci], toks[ci]); break;
                case COL_X0:   cd->x0  = parse_candval(fn, lineno, col_name[ci], toks[ci]); break;
                case COL_X1:   cd->x1  = parse_candval(fn, lineno, col_name[ci], toks[ci]); break;
                case COL_X2:   cd->x2  = parse_candval(fn, lineno, col_name[ci], toks[ci]); break;
                }
            }
            if (cd->dm < 0.0)
                candfile_error(fn, lineno, "DM must be >= 0", NULL);
            if (cd->x0 <= 0.0)
                candfile_error(fn, lineno,
                               is_period ? "period 'p' must be > 0" : "frequency 'f' must be > 0",
                               NULL);
            ncand++;
        }
    }
    fclose(f);

    if (!header_done)
        candfile_error(fn, 0, "no header line found", NULL);
    if (ncand == 0)
        candfile_error(fn, 0, "no candidates found", NULL);
    *ncand_out = ncand;
    *is_period_out = is_period;
    return cands;
}

/* ------------------------------------------------------------------------- */
/* Per-candidate setup helpers (the explicit-p/f subset of prepfold's         */
/* resolve_output_names / resolve_fold_params, driven by a candidate row).    */
/* ------------------------------------------------------------------------- */

/* Has candidate name `nm` already been assigned to one of fc[0..c-1]?          */
static int candnm_taken(foldcand *fc, int c, const char *nm)
{
    int j;
    for (j = 0; j < c; j++)
        if (fc[j].search.candnm && !strcmp(fc[j].search.candnm, nm))
            return 1;
    return 0;
}

/* Build candidate c's name (given, or autogenerated exactly as prepfold does)   */
/* and its .pfd output / plot / pgplot-device names.  Duplicate names (e.g. two  */
/* unnamed candidates with the same period, or two identical 'name' columns)     */
/* would map to the same output files and silently clobber each other, so any    */
/* collision against an earlier candidate is disambiguated with a '_N' suffix.   */
static void multi_output_names(foldcand *fc_array, int c, char *rootnm,
                               multicand *cand, int is_period)
{
    foldcand *fc = &fc_array[c];
    prepfoldinfo *search = &fc->search;
    char *base = cand->name ? dupstr(cand->name)
                            : make_autogen_candnm(is_period, cand->x0);

    if (candnm_taken(fc_array, c, base)) {
        char *uniq = (char *) chkmalloc(strlen(base) + 16, "candidate name");
        int k = 2;
        do {
            sprintf(uniq, "%s_%d", base, k++);
        } while (candnm_taken(fc_array, c, uniq));
        fprintf(stderr,
                "prepfold_multi: warning: duplicate candidate name '%s'; using"
                " '%s' to avoid overwriting output.\n", base, uniq);
        free(base);
        base = uniq;
    }
    search->candnm = base;

    build_output_filenames(fc, rootnm);
}

/* Convert the candidate's explicit p/f (+ derivatives) into fold frequencies   */
/* and seed search->topo/bary, then pick the profile length -- byte-for-byte    */
/* the explicit-p/f path of prepfold's resolve_fold_params.                     */
static void multi_resolve_fold_params(Cmdline *cmd, infodata *idata, foldcand *fc,
                                      multicand *cand, int is_period)
{
    prepfoldinfo *search = &fc->search;
    double f, fd, fdd, p, pd, pdd;

    /* Convert the candidate's explicit period or frequency (+ derivatives) into  */
    /* both representations with the shared converter, so the period<->frequency  */
    /* math lives in exactly one place (characteristics.c switch_f_and_p).        */
    if (is_period) {
        p = cand->x0; pd = cand->x1; pdd = cand->x2;
        switch_f_and_p(p, pd, pdd, &f, &fd, &fdd);
    } else {
        f = cand->x0; fd = cand->x1; fdd = cand->x2;
        switch_f_and_p(f, fd, fdd, &p, &pd, &pdd);
    }
    if (idata->bary) {
        search->bary.p1 = p; search->bary.p2 = pd; search->bary.p3 = pdd;
    } else {
        search->topo.p1 = p; search->topo.p2 = pd; search->topo.p3 = pdd;
    }
    fc->f = f;
    fc->fd = fd;
    fc->fdd = fdd;

    resolve_default_proflen(cmd, search);
}

/* Barycentric folding epoch for one candidate's DM (matches the DM correction  */
/* in compute_obs_timing: barycenter the topo epoch, then subtract the          */
/* dispersion delay to infinite frequency).                                     */
static double multi_bary_epoch(infodata *idata, double tepoch, double dm,
                               char *rastring, char *decstring, char *obs, char *ephem)
{
    double bepoch = 0.0, dtmp;
    barycenter(&tepoch, &bepoch, &dtmp, 1, rastring, decstring, obs, ephem);
    /* Shared epoch-to-infinite-frequency correction (same as compute_obs_timing). */
    return bary_epoch_to_infinite_freq(idata, bepoch, dm);
}

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    double recdt = 0.0, N = 0.0, T = 0.0, startTday = 0.0, proftime;
    double *events = NULL;
    char *rootnm, *shared_filenm, *shared_telescope;
    char obs[3], ephem[6], pname[30], rastring[50], decstring[50];
    int numevents = 0, numchan = 1, useshorts = 0;
    int ptsperrec = 1, insubs = 0, shared_nsub, ncand = 0, is_period = 0;
    int *maskchans = NULL;
    long worklen = 0, reads_per_part = 0, blockfloats;
    long long lorec = 0, numrec = 0, totnumfolded = 0, totnumblocks = 0;
    long ii, jj;
    int c;
    double tepoch, shared_dt;
    struct spectra_info s;
    infodata idata;
    foldcand fcs, *fc;
    multicand *cands;
    Cmdline *cmd;
    plotflags pflags;
    mask obsmask;
    /* Per-candidate fold scratch, plus the SHARED barycentric tables: the topo<->  */
    /* bary times and average v/c depend only on the epoch/obs/sky, so they are     */
    /* computed once (one TEMPO call) and reused for every candidate.               */
    double **buffers, **phasesadded;
    int shared_numbarypts = 0;
    double *shared_barytimes = NULL, *shared_topotimes = NULL, shared_avgvoverc = 0.0;
    rawblock_reader *reader = NULL;
    float **subdatas = NULL, *insub_data = NULL;
    int padding = 0, nummasked = 0;
    int nthreads = 1;
    size_t fold_mem_bytes = 0, fold_mem_warn = 0, fold_mem_cap = 0;

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(0);
    }
    pname[0] = '\0';

    cmd = parseCmdline(argc, argv);

    /* Parse the candidate list up front so we fail fast on a bad file. */
    cands = parse_candfile(cmd->candfile, &ncand, &is_period);

    /* Shared spectra_info setup (mirrors prepfold's main()). */
    init_foldcand(&fcs);
    fcs.dm = cands[0].dm;
    fcs.nsub = cmd->nsub;
    fcs.npart = cmd->npart;
    fcs.phs = cmd->phs;
    fcs.flags = 1;
    /* Placeholders so the shared-setup functions' progress prints (which expect */
    /* per-candidate names) have valid strings; the real names are per-candidate.*/
    fcs.search.candnm = dupstr("shared");
    fcs.outfilenm = dupstr("(prepfold_multi)");
    fcs.plotfilenm = dupstr("(prepfold_multi)");

    spectra_info_set_defaults(&s);
    s.filenames = cmd->argv;
    s.num_files = cmd->argc;
    if (cmd->zerodmP)
        cmd->noclipP = 1;
    s.clip_sigma = cmd->clip;
    s.apply_flipband = (cmd->invertP) ? 1 : -1;
    s.apply_weight = (cmd->noweightsP) ? 0 : -1;
    s.apply_scale = (cmd->noscalesP) ? 0 : -1;
    s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
    s.remove_zerodm = (cmd->zerodmP) ? 1 : 0;
    if (cmd->ncpus > 1) {
#ifdef _OPENMP
        int maxcpus = omp_get_num_procs();
        int openmp_numthreads = (cmd->ncpus <= maxcpus) ? cmd->ncpus : maxcpus;
        omp_set_dynamic(0);
        omp_set_num_threads(openmp_numthreads);
        printf("Using %d threads with OpenMP\n\n", openmp_numthreads);
#endif
    } else {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
    }
    if (cmd->noclipP) {
        cmd->clip = 0.0;
        s.clip_sigma = 0.0;
    }
    if (cmd->ifsP)
        s.use_poln = cmd->ifs + 1;
    obsmask.numchan = obsmask.numint = 0;
    if (cmd->slowP) {
        cmd->fineP = 1;
        if (!cmd->proflenP)
            fcs.search.proflen = 100;
    }
    if (cmd->fineP) {
        cmd->ndmfact = 1;
        cmd->dmstep = 1;
        cmd->npfact = 1;
        cmd->pstep = 1;
        cmd->pdstep = 2;
    }
    if (cmd->coarseP) {
        cmd->npfact = 4;
        cmd->pstep = (cmd->pstep == 1) ? 2 : 3;
        cmd->pdstep = (cmd->pdstep == 2) ? 4 : 6;
    }
    pflags.events = 0;
    pflags.nosearch = cmd->nosearchP;
    pflags.scaleparts = cmd->scalepartsP;
    pflags.justprofs = cmd->justprofsP;
    pflags.allgrey = cmd->allgreyP;
    pflags.fixchi = cmd->fixchiP;
    pflags.samples = cmd->samplesP;
    pflags.showfold = 0;

    printf("\n\n");
    printf("     Multi-Candidate Pulsar Raw-Data Folding Routine\n");
    printf("    Folds many candidates in one pass over the raw data.\n");
    printf("                 by Scott M. Ransom\n\n");

    /* Identify/open the input once.  This also resolves the auto-nsub and reads */
    /* the mask -- all shared across candidates. */
    identify_and_open_input(cmd, &s, &idata, &fcs, &obsmask, &pflags,
                            &rootnm, &ptsperrec, &numrec, &numchan, &useshorts,
                            &insubs, &events, &numevents, &T);

    /* Scope: per-candidate DM only makes sense for raw data or PRESTO subbands. */
    if (cmd->eventsP || !(RAWDATA || insubs)) {
        fprintf(stderr,
                "\nprepfold_multi: input must be raw radio data (PSRFITS/filterbank)\n"
                "  or PRESTO subbands (*.subNN).  Already-dedispersed .dat time series\n"
                "  and event files have no per-candidate DM; use prepfold instead.\n\n");
        exit(1);
    }
    shared_nsub = fcs.nsub;     /* resolved auto-nsub (or cmd->nsub / num subbands) */

    /* Compute the (candidate-independent) observation timing once. */
    compute_obs_timing(cmd, &s, &idata, &fcs, &obsmask, insubs, useshorts,
                       ptsperrec, obs, ephem, rastring, decstring, &numchan,
                       &maskchans, &recdt, &lorec, &startTday, &numrec,
                       &reads_per_part, &T, &N, &worklen);
    tepoch = fcs.search.tepoch;
    shared_dt = fcs.search.dt;
    shared_filenm = fcs.search.filenm;
    shared_telescope = fcs.search.telescope;
    proftime = worklen * shared_dt;
    blockfloats = (long) shared_nsub * worklen;

    printf("\nFolding %d candidate%s (%s) from '%s'.\n", ncand,
           ncand == 1 ? "" : "s", is_period ? "period" : "frequency", cmd->candfile);
    printf("Reading the raw data once; dedispersing and folding per candidate.\n\n");

    /* The barycentric time tables are candidate-INDEPENDENT, so build them once  */
    /* here (one TEMPO call) and reuse for every candidate.  Skipped when folding  */
    /* topocentrically (-topo), exactly as the single-candidate path decides.      */
    if (!cmd->topoP)
        compute_bary_table(tepoch, T, rastring, decstring, obs, ephem,
                           &shared_numbarypts, &shared_barytimes,
                           &shared_topotimes, &shared_avgvoverc);

    /* Per-candidate context and scratch. */
    fc = (foldcand *) chkmalloc(ncand * sizeof(foldcand), "candidate array");
    buffers = (double **) chkmalloc(ncand * sizeof(double *), "fold buffers");
    phasesadded = (double **) chkmalloc(ncand * sizeof(double *), "phase buffers");
    fold_mem_limits(&fold_mem_warn, &fold_mem_cap);

    for (c = 0; c < ncand; c++) {
        prepfoldinfo *search;
        init_foldcand(&fc[c]);
        fc[c].dm = cands[c].dm;
        fc[c].nsub = shared_nsub;
        fc[c].npart = fcs.npart;
        fc[c].phs = cmd->phs;
        fc[c].flags = 1;
        search = &fc[c].search;

        /* Shared observation fields, copied/recomputed per candidate.  This setup */
        /* loop is SERIAL, so setting idata.dm on the one shared idata struct just  */
        /* before the code that reads it (set_posn, compute_dispersion_delays) is   */
        /* safe here.  The parallel optimize loop later takes a per-candidate idata */
        /* copy, since idata.dm cannot be shared across concurrent candidates.      */
        idata.dm = fc[c].dm;
        search->dt = shared_dt;
        search->startT = cmd->startT;
        search->endT = cmd->endT;
        search->tepoch = tepoch;
        search->filenm = dupstr(shared_filenm);
        search->telescope = dupstr(shared_telescope);
        set_posn(search, &idata);
        if (idata.mjd_i && !cmd->topoP)
            search->bepoch = multi_bary_epoch(&idata, tepoch, fc[c].dm,
                                              rastring, decstring, obs, ephem);

        multi_output_names(fc, c, rootnm, &cands[c], is_period);
        printf("\n==> Candidate %d/%d: '%s'\n", c + 1, ncand, search->candnm);
        printf("Output data file is '%s'.\n", fc[c].outfilenm);
        printf("Output plot file is '%s'.\n", fc[c].plotfilenm);

        multi_resolve_fold_params(cmd, &idata, &fc[c], &cands[c], is_period);

        /* Memory guard: tally this candidate's fold arrays (rawfolds dominate)  */
        /* before allocating them, and bail with a clear message if the running  */
        /* total crosses the hard cap -- so a huge candidate list fails here     */
        /* rather than deep inside gen_dvect() or via the OOM killer.            */
        {
            size_t cand_bytes =
                (size_t) fc[c].nsub * fc[c].npart * search->proflen * sizeof(double)
                + (size_t) fc[c].nsub * fc[c].npart * sizeof(foldstats)
                + (size_t) fc[c].nsub * search->proflen * sizeof(double)
                + (size_t) fc[c].nsub * sizeof(double);
            fold_mem_bytes += cand_bytes;
            if (fold_mem_bytes > fold_mem_cap) {
                fprintf(stderr,
                        "\nprepfold_multi: estimated fold memory (%.1f GB by candidate"
                        " %d/%d, '%s')\n  exceeds the %.1f GB physical-RAM safety cap."
                        "  Fold fewer candidates per\n  run, or fold individual"
                        " candidates with prepfold.\n\n",
                        fold_mem_bytes / BYTES_PER_GB, c + 1, ncand, search->candnm,
                        fold_mem_cap / BYTES_PER_GB);
                exit(1);
            }
        }

        print_fold_info(cmd, &fc[c], N, T, pname);
        allocate_fold_arrays(cmd, &fc[c]);

        /* Barycentric fold-frequency correction (or fold topocentrically).  The  */
        /* shared tables built above feed each candidate's cheap per-candidate     */
        /* bary2topo correction -- no extra TEMPO call here.                       */
        if (cmd->topoP) {
            fc[c].foldf = fc[c].f;
            fc[c].foldfd = fc[c].fd;
            fc[c].foldfdd = fc[c].fdd;
            fc[c].orig_foldf = fc[c].foldf;
        } else {
            search->avgvoverc = shared_avgvoverc;
            apply_bary_fold_correction(&fc[c], shared_numbarypts,
                                       shared_barytimes, shared_topotimes);
        }
        if (idata.bary)
            search->fold.pow = 1.0;
        search->fold.p1 = fc[c].foldf;
        search->fold.p2 = fc[c].foldfd;
        search->fold.p3 = fc[c].foldfdd;

        /* Per-candidate dispersion delays at this candidate's DM. */
        compute_dispersion_delays(cmd, &idata, &fc[c], numchan, insubs);

        /* Over-dispersion guard (raw data only; mirrors read_subbands' guard). */
        if (RAWDATA &&
            !rawblock_dispersion_fits(fc[c].idispdts[0], s.spectra_per_subint)) {
            fprintf(stderr,
                    "\nprepfold_multi: candidate '%s' (DM=%.3f) over-disperses a raw\n"
                    "  block: the lowest-channel delay (%d) exceeds the %d spectra per\n"
                    "  subint.  Lower the DM or fold this candidate with prepfold.\n\n",
                    search->candnm, fc[c].dm, fc[c].idispdts[0], s.spectra_per_subint);
            exit(1);
        }

        /* Sub-integration time grid for the optimizer. */
        fc[c].parttimes = gen_dvect(fc[c].npart);
        for (ii = 0; ii < fc[c].npart; ii++)
            fc[c].parttimes[ii] = ii * reads_per_part * proftime;

        /* Fold scratch (per candidate: proflen varies with the spin period). */
        buffers[c] = gen_dvect(fc[c].nsub * search->proflen);
        phasesadded[c] = gen_dvect(fc[c].nsub);
        for (ii = 0; ii < fc[c].nsub * search->proflen; ii++)
            buffers[c][ii] = 0.0;
        for (ii = 0; ii < fc[c].nsub; ii++)
            phasesadded[c][ii] = 0.0;
    }

    /*
     *   The single read-once/dedisperse-many fold pass.
     */

    /* Per-thread (not per-candidate) subband scratch: the fold pass dedisperses */
    /* each candidate's block into a private buffer so the candidate loop can run */
    /* in parallel without sharing the dedispersion output.  One buffer per       */
    /* worker thread keeps the scratch footprint independent of the candidate     */
    /* count.                                                                     */
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    subdatas = (float **) chkmalloc(nthreads * sizeof(float *), "per-thread scratch");
    for (ii = 0; ii < nthreads; ii++)
        subdatas[ii] = gen_fvect(blockfloats);

    fold_mem_bytes += (size_t) nthreads * blockfloats * sizeof(float);
    if (fold_mem_bytes > fold_mem_warn)
        printf("Estimated fold memory: %.2f GB for %d candidate%s "
               "(+ %d-thread scratch).\n",
               fold_mem_bytes / BYTES_PER_GB, ncand, ncand == 1 ? "" : "s", nthreads);
    else
        printf("Estimated fold memory: %.1f MB for %d candidate%s "
               "(+ %d-thread scratch).\n",
               fold_mem_bytes / (1024.0 * 1024.0), ncand, ncand == 1 ? "" : "s",
               nthreads);

    if (RAWDATA) {
        printf("\rTrue starting fraction       =  %g\n",
               (double) (lorec * ptsperrec) / s.N);
        offset_to_spectra(lorec * (long long) ptsperrec, &s);
        reader = rawblock_reader_init(&s, &obsmask);
        /* Seed the masking start time for startT != 0 (no-op when lorec == 0). */
        reader->currentspectra = lorec * (long long) ptsperrec;
        /* Prime the reader (fills the "last" slot; produces no foldable data). */
        (void) read_clean_rawblock(reader, &s, maskchans, &nummasked, &obsmask,
                                   &padding);
    } else {                    /* insubs */
        insub_data = gen_fvect(blockfloats);
        for (ii = 0; ii < s.num_files; ii++)
            chkfileseek(s.files[ii], lorec * SUBSBLOCKLEN, sizeof(short), SEEK_SET);
    }

    printf("Folding...\n");
    printf("  Folded %lld points of %.0f", totnumfolded, N);
    for (ii = 0; ii < fcs.npart; ii++) {
        for (jj = 0; jj < reads_per_part; jj++) {
            long numread;
            double fold_time0 = ii * reads_per_part * proftime + jj * proftime;

            if (RAWDATA) {
                numread = read_clean_rawblock(reader, &s, maskchans, &nummasked,
                                              &obsmask, &padding);
                /* Parallelize over candidates: the cleaned block (current/last)  */
                /* is read-only here, and each candidate writes only its own       */
                /* dedispersion scratch and rawfolds/buffers, so there is no       */
                /* cross-candidate contention.  The reader advance/swap happens    */
                /* AFTER this region, so no candidate reads last_clean post-swap.  */
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
                for (c = 0; c < ncand; c++) {
                    int tid = 0;
#ifdef _OPENMP
                    tid = omp_get_thread_num();
#endif
                    float *sub_c = subdatas[tid];
                    dedisp_subbands(reader->current_clean, reader->last_clean,
                                    worklen, numchan, fc[c].idispdts, fc[c].nsub,
                                    sub_c);
                    fold_subband_block(cmd, &fc[c], sub_c, worklen, numread, ii,
                                       fold_time0, fc[c].foldf, fc[c].foldfd,
                                       fc[c].foldfdd, buffers[c], phasesadded[c]);
                }
                rawblock_reader_advance(reader);
            } else {            /* insubs: subbands already dedispersed at the */
                                /* nominal DM; the candidate DM enters at the   */
                                /* optimize stage (correct_subbands_for_DM).    */
                numread = read_PRESTO_subbands(s.files, s.num_files, insub_data,
                                               recdt, maskchans, &nummasked,
                                               &obsmask, s.padvals);
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
                for (c = 0; c < ncand; c++) {
                    int tid = 0;
                    float *src = insub_data;
#ifdef _OPENMP
                    tid = omp_get_thread_num();
#endif
                    if (cmd->runavgP) {     /* runavg mutates data: fold a copy */
                        src = subdatas[tid];
                        memcpy(src, insub_data,
                               (size_t) blockfloats * sizeof(float));
                    }
                    fold_subband_block(cmd, &fc[c], src, worklen, numread, ii,
                                       fold_time0, fc[c].foldf, fc[c].foldfd,
                                       fc[c].foldfdd, buffers[c], phasesadded[c]);
                }
            }
            totnumfolded += numread;
            totnumblocks++;
        }
        printf("\r  Folded %lld points of %.0f", totnumfolded, N);
        fflush(NULL);
    }
    printf("\n");

    /* Every candidate is folded inside the same per-block iteration, sharing the */
    /* one reader and the one `numread` per block, so all candidates necessarily  */
    /* see an identical block count and identical EOF/padding handling.  Verify    */
    /* the loop processed the full npart x reads_per_part block grid -- a hard      */
    /* check, not assert(), so it still guards builds compiled with -DNDEBUG.       */
    if (totnumblocks != (long long) fcs.npart * reads_per_part) {
        fprintf(stderr,
                "\nprepfold_multi: internal error: folded %lld raw blocks but expected"
                " %lld\n  (%d parts x %ld reads per part).  Please report this.\n\n",
                totnumblocks, (long long) fcs.npart * reads_per_part,
                fcs.npart, reads_per_part);
        exit(1);
    }

    if (RAWDATA)
        free_rawblock_reader(reader);
    close_rawfiles(&s);

    /*
     *   Per-candidate optimization and output.
     */

    /* The per-candidate optimize/output stage is embarrassingly parallel -- each  */
    /* candidate searches and writes only its own foldcand.  Two shared-state       */
    /* hazards are removed with per-candidate copies: (1) optimize_candidate()       */
    /* settles a few search flags in the cmd and pflags structs (identical for       */
    /* every candidate since they share nsub) and would otherwise race; (2) idata.dm */
    /* must hold THIS candidate's DM.  infodata/Cmdline/plotflags own no heap that    */
    /* these functions mutate, so a by-value struct copy fully isolates each thread.  */
    /* PGPLOT and the .pfd writer are not thread-safe, so only the output half runs   */
    /* under a critical section; the expensive search runs fully in parallel.         */
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for (c = 0; c < ncand; c++) {
        Cmdline cmd_c = *cmd;
        infodata idata_c = idata;
        plotflags pflags_c = pflags;

        idata_c.dm = fc[c].dm;  /* per-candidate DM (was a shared-idata mutation) */
        optimize_candidate(&cmd_c, &idata_c, &fc[c], &pflags_c, T);

        /* The shared bary tables are read-only here, so every candidate is handed */
        /* the same tables.  Serialize the plot/write (PGPLOT + file + stdout).    */
#ifdef _OPENMP
#pragma omp critical (prepfold_multi_output)
#endif
        {
            printf("\n\nResults for candidate %d/%d ('%s'):\n\n",
                   c + 1, ncand, fc[c].search.candnm);
            write_results_and_plot(&cmd_c, &idata_c, &fc[c], &pflags_c, N, insubs,
                                   shared_barytimes, shared_topotimes,
                                   shared_numbarypts);
        }

        free_foldcand(&fc[c], &cmd_c, &idata_c);
        vect_free(buffers[c]);
        vect_free(phasesadded[c]);
    }

    /* Free shared state. */
    for (ii = 0; ii < nthreads; ii++)
        vect_free(subdatas[ii]);
    free(subdatas);
    if (insub_data)
        vect_free(insub_data);
    if (!cmd->outfileP)
        free(rootnm);
    if (cmd->maskfileP) {
        free_mask(obsmask);
        if (maskchans)
            vect_free(maskchans);
    }
    free(shared_filenm);
    free(shared_telescope);
    free(fcs.search.candnm);
    free(fcs.outfilenm);
    free(fcs.plotfilenm);
    for (c = 0; c < ncand; c++)
        free(cands[c].name);
    free(cands);
    free(fc);
    if (shared_barytimes)
        vect_free(shared_barytimes);
    if (shared_topotimes)
        vect_free(shared_topotimes);
    free(buffers);
    free(phasesadded);
    free(events);

    printf("Done.\n\n");
    return 0;
}
