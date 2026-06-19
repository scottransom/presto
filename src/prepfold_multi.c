#include <ctype.h>
#include "prepfold.h"
#include "prepfold_multi_cmd.h"
#include "mask.h"
#include "backend_common.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

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
    char *d = (char *) malloc(strlen(s) + 1);
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

/* Build the candidate name (given, or autogenerated exactly as prepfold does) */
/* and the .pfd output / plot / pgplot-device names.                          */
static void multi_output_names(foldcand *fc, char *rootnm, multicand *cand,
                               int is_period)
{
    prepfoldinfo *search = &fc->search;

    /* Candidate name: the given one, or prepfold's autogenerated format (shared */
    /* helper).  Output/plot/device names via the same builder prepfold uses.    */
    if (cand->name)
        search->candnm = dupstr(cand->name);
    else
        search->candnm = make_autogen_candnm(is_period, cand->x0);

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
    long long lorec = 0, numrec = 0, totnumfolded = 0;
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
    float *subdata = NULL, *insub_data = NULL;
    int padding = 0, nummasked = 0;

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
    fc = (foldcand *) malloc(ncand * sizeof(foldcand));
    buffers = (double **) malloc(ncand * sizeof(double *));
    phasesadded = (double **) malloc(ncand * sizeof(double *));

    for (c = 0; c < ncand; c++) {
        prepfoldinfo *search;
        init_foldcand(&fc[c]);
        fc[c].dm = cands[c].dm;
        fc[c].nsub = shared_nsub;
        fc[c].npart = fcs.npart;
        fc[c].phs = cmd->phs;
        fc[c].flags = 1;
        search = &fc[c].search;

        /* Shared observation fields, copied/recomputed per candidate.  NOTE:     */
        /* idata is a SINGLE shared struct; idata.dm is set here for each          */
        /* candidate immediately before the code that reads it (set_posn,          */
        /* compute_dispersion_delays).  Stage 5 must give each thread its own       */
        /* idata before parallelizing over candidates -- idata.dm cannot be         */
        /* shared across concurrent candidates.                                     */
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

        multi_output_names(&fc[c], rootnm, &cands[c], is_period);
        printf("\n==> Candidate %d/%d: '%s'\n", c + 1, ncand, search->candnm);
        printf("Output data file is '%s'.\n", fc[c].outfilenm);
        printf("Output plot file is '%s'.\n", fc[c].plotfilenm);

        multi_resolve_fold_params(cmd, &idata, &fc[c], &cands[c], is_period);
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

    subdata = gen_fvect(blockfloats);
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
                for (c = 0; c < ncand; c++) {
                    dedisp_subbands(reader->current_clean, reader->last_clean,
                                    worklen, numchan, fc[c].idispdts, fc[c].nsub,
                                    subdata);
                    fold_subband_block(cmd, &fc[c], subdata, worklen, numread, ii,
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
                for (c = 0; c < ncand; c++) {
                    float *src = insub_data;
                    if (cmd->runavgP) {     /* runavg mutates data: fold a copy */
                        memcpy(subdata, insub_data,
                               (size_t) blockfloats * sizeof(float));
                        src = subdata;
                    }
                    fold_subband_block(cmd, &fc[c], src, worklen, numread, ii,
                                       fold_time0, fc[c].foldf, fc[c].foldfd,
                                       fc[c].foldfdd, buffers[c], phasesadded[c]);
                }
            }
            totnumfolded += numread;
        }
        printf("\r  Folded %lld points of %.0f", totnumfolded, N);
        fflush(NULL);
    }
    printf("\n");

    if (RAWDATA)
        free_rawblock_reader(reader);
    close_rawfiles(&s);

    /*
     *   Per-candidate optimization and output.
     */

    for (c = 0; c < ncand; c++) {
        printf("\n\nOptimizing candidate %d/%d ('%s')...\n\n",
               c + 1, ncand, fc[c].search.candnm);
        idata.dm = fc[c].dm;    /* shared idata: set this candidate's DM before use */
        optimize_candidate(cmd, &idata, &fc[c], &pflags, T);
        /* The shared bary tables are read-only in the output stage, so every     */
        /* candidate is handed the same tables.                                   */
        write_results_and_plot(cmd, &idata, &fc[c], &pflags, N, insubs,
                               shared_barytimes, shared_topotimes,
                               shared_numbarypts);
        free_foldcand(&fc[c], cmd, &idata);
        vect_free(buffers[c]);
        vect_free(phasesadded[c]);
    }

    /* Free shared state. */
    vect_free(subdata);
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
