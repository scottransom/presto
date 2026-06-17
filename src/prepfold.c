#include <ctype.h>
#include "prepfold.h"
#include "prepfold_cmd.h"
#include "mask.h"
#include "backend_common.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

extern int getpoly(double mjd, double duration, double *dm, FILE * fp, char *pname);
extern int phcalc(double mjd0, double mjd1, int last_index,
                  double *phase, double *psrfreq);
extern int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);
extern char *make_polycos(char *parfilenm, infodata * idata, char *polycofilenm,
                          int debug_tempo);
extern int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals);
void set_posn(prepfoldinfo * in, infodata * idata);

/*
 * The main program
 */

int main(int argc, char *argv[])
{
    float *data = NULL;
    double recdt = 0.0, N = 0.0, T = 0.0, startTday = 0.0;
    double *barytimes = NULL, *topotimes = NULL;
    double *events = NULL;
    char *rootnm;
    char obs[3], ephem[6], pname[30], rastring[50], decstring[50];
    int numevents = 0, numchan = 1, numbarypts = 0;
    int ptsperrec = 1, useshorts = 0;
    int *maskchans = NULL, insubs = 0;
    long worklen = 0, reads_per_part = 0;
    long long lorec = 0;
    long long numrec = 0;
    struct spectra_info s;
    infodata idata;
    foldcand fc;
    Cmdline *cmd;
    plotflags pflags;
    mask obsmask;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(0);
    }
    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

    /* Set up the per-candidate fold context from the parsed defaults.  From  */
    /* here on, cmd is read-only parsed input; computed per-candidate state    */
    /* (dm, nsub, npart, phs, proflen, spin params, ...) lives in fc.          */
    init_foldcand(&fc);
    fc.dm = cmd->dm;
    fc.nsub = cmd->nsub;
    fc.npart = cmd->npart;
    fc.phs = cmd->phs;
    fc.flags = 1;

    spectra_info_set_defaults(&s);
    s.filenames = cmd->argv;
    s.num_files = cmd->argc;
    // If we are zeroDMing, make sure that clipping is off.
    if (cmd->zerodmP)
        cmd->noclipP = 1;
    s.clip_sigma = cmd->clip;
    // -1 causes the data to determine if we use weights, scales, &
    // offsets for PSRFITS or flip the band for any data type where
    // we can figure that out with the data
    s.apply_flipband = (cmd->invertP) ? 1 : -1;
    s.apply_weight = (cmd->noweightsP) ? 0 : -1;
    s.apply_scale = (cmd->noscalesP) ? 0 : -1;
    s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
    s.remove_zerodm = (cmd->zerodmP) ? 1 : 0;
    if (cmd->ncpus > 1) {
#ifdef _OPENMP
        int maxcpus = omp_get_num_procs();
        int openmp_numthreads = (cmd->ncpus <= maxcpus) ? cmd->ncpus : maxcpus;
        // Make sure we are not dynamically setting the number of threads
        omp_set_dynamic(0);
        omp_set_num_threads(openmp_numthreads);
        printf("Using %d threads with OpenMP\n\n", openmp_numthreads);
#endif
    } else {
#ifdef _OPENMP
        omp_set_num_threads(1); // Explicitly turn off OpenMP
#endif
    }
    if (cmd->noclipP) {
        cmd->clip = 0.0;
        s.clip_sigma = 0.0;
    }
    if (cmd->ifsP) {
        // 0 = default or summed, 1-4 are possible also
        s.use_poln = cmd->ifs + 1;
    }
    obsmask.numchan = obsmask.numint = 0;
    if (cmd->timingP) {
        cmd->nosearchP = 1;
        cmd->nopsearchP = 1;
        cmd->nopdsearchP = 1;
        cmd->nodmsearchP = 1;
        if (cmd->npart == 64)   /* The default value */
            fc.npart = 60;
        cmd->fineP = 1;
    }
    if (cmd->slowP) {
        cmd->fineP = 1;
        if (!cmd->proflenP)
            fc.search.proflen = 100;
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
        if (cmd->pstep == 1)
            cmd->pstep = 2;
        else
            cmd->pstep = 3;
        if (cmd->pdstep == 2)
            cmd->pdstep = 4;
        else
            cmd->pdstep = 6;
    }
    pflags.events = cmd->eventsP;
    pflags.nosearch = cmd->nosearchP;
    pflags.scaleparts = cmd->scalepartsP;
    pflags.justprofs = cmd->justprofsP;
    pflags.allgrey = cmd->allgreyP;
    pflags.fixchi = cmd->fixchiP;
    pflags.samples = cmd->samplesP;
    pflags.showfold = 0;

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("        Pulsar Raw-Data Folding Search Routine\n");
    printf(" Used for DM, Period, and P-dot tweaking of PSR candidates.\n");
    printf("                 by Scott M. Ransom\n\n");

    /* Identify the data type, open the input file(s), read any mask, and    */
    /* (for events) read the event list.                                     */
    identify_and_open_input(cmd, &s, &idata, &fc, &obsmask, &pflags,
                            &rootnm, &ptsperrec, &numrec, &numchan, &useshorts,
                            &insubs, &events, &numevents, &T);

    /* Manipulate the file names we will use  */
    resolve_output_names(cmd, &fc, rootnm);

    /* Compute the observation timing and set the source position. */
    compute_obs_timing(cmd, &s, &idata, &fc, &obsmask, insubs, useshorts,
                       ptsperrec, obs, ephem, rastring, decstring, &numchan,
                       &maskchans, &recdt, &lorec, &startTday, &numrec,
                       &reads_per_part, &T, &N, &worklen);

    /* Resolve the fold parameters and determine the profile length. */
    resolve_fold_params(cmd, &idata, &fc, insubs, T, startTday, recdt, lorec,
                        pname);

    /* Determine the phase delays caused by the orbit if needed. */
    setup_orbit_delays(cmd, &idata, &fc, insubs, T, &N);

    /* Output some informational data on the screen. */
    print_fold_info(cmd, &fc, N, T, pname);

    /* Allocate and initialize the fold arrays. */
    allocate_fold_arrays(cmd, &fc);

    if (cmd->eventsP) {         /* Fold events instead of a time series */
        fold_events(cmd, &idata, &fc, T, numevents, events, startTday);
    } else {                    /* Fold a time series */
        fold_timeseries(cmd, &s, &idata, &fc, &obsmask, T, N, lorec, ptsperrec,
                        recdt, worklen, reads_per_part, numchan, startTday,
                        useshorts, insubs, maskchans, obs, ephem, rastring,
                        decstring, &data, &numbarypts, &barytimes, &topotimes);
    }
    // This resets foldf (which is used below) to the original value
    if (cmd->polycofileP)
        fc.foldf = fc.orig_foldf;

    //  Close all the raw files and free their vectors
    close_rawfiles(&s);

    /*
     *   Perform the candidate optimization search
     */

    printf("\n\nOptimizing...\n\n");

    optimize_candidate(cmd, &idata, &fc, &pflags, T);

    printf("  Done searching.\n\n");

    /* Write and plot the results. */
    write_results_and_plot(cmd, &idata, &fc, &pflags, N, insubs, barytimes,
                           topotimes, numbarypts);

    /* Free our memory  */
    free_foldcand(&fc, cmd, &idata);
    cleanup_fold(cmd, &obsmask, insubs, data, rootnm, barytimes, topotimes);

    printf("Done.\n\n");
    return (0);
}
