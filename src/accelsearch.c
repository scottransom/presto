#include "accel.h"

/*#undef USEMMAP*/

#ifdef USEMMAP
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
extern void set_openmp_numthreads(int numthreads);
#endif

extern float calc_median_powers(fcomplex * amplitudes, int numamps);
extern void zapbirds(double lobin, double hibin, long numbins, FILE * fftfile,
                     fcomplex * fft, int constamp);

static void print_percent_complete(int current, int number, char *what, int reset)
{
    static int newper = 0, oldper = -1;

    if (reset) {
        oldper = -1;
        newper = 0;
    } else {
        newper = (int) (current / (float) (number) * 100.0);
        if (newper < 0)
            newper = 0;
        if (newper > 100)
            newper = 100;
        if (newper > oldper) {
            printf("\rAmount of %s complete = %3d%%", what, newper);
            fflush(stdout);
            oldper = newper;
        }
    }
}

/* -perfinfo:  what decides the speed (see print_perfinfo()) */
static void show_perfinfo(void)
{
#ifdef _OPENMP
    print_perfinfo(accel_simd_name(), omp_get_max_threads());
#else
    print_perfinfo(accel_simd_name(), -1);
#endif
}

int main(int argc, char *argv[])
{
    int ii, rstep;
    double ttim, utim, stim, tott;
    double clk_tck = (double) sysconf(_SC_CLK_TCK);
    struct tms runtimes;
    subharminfo **subharminfs;
    accelobs obs;
    infodata idata;
    GSList *cands = NULL;
    Cmdline *cmd;

    /* Prep the timer */

    tott = times(&runtimes) / clk_tck;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

    /* With no input file, -perfinfo just reports and exits */
    if (cmd->argc == 0) {
        if (!cmd->perfinfoP) {
            printf("\n");
            usage();
            exit(1);
        }
#ifdef _OPENMP
        omp_set_num_threads(cmd->ncpus > 1 ? cmd->ncpus : 1);
#endif
        show_perfinfo();
        exit(0);
    }

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("    Fourier-Domain Acceleration and Jerk Search Routine\n");
    printf("                    by Scott M. Ransom\n\n");

    /* Create the accelobs structure */
    create_accelobs(&obs, &idata, cmd, 1);

    /* The step-size of blocks to walk through the input data */
    rstep = obs.corr_uselen * ACCEL_DR;

    /* Zap birdies if requested and if in memory */
    if (cmd->zaplistP && !obs.mmap_file && obs.fft) {
        int numbirds;
        double *bird_lobins, *bird_hibins, hibin;

        /* Read the Standard bird list */
        numbirds = get_birdies(cmd->zaplist, obs.T, cmd->baryv,
                               &bird_lobins, &bird_hibins);

        /* Zap the birdies */
        printf("Zapping them using a barycentric velocity of %.5gc.\n\n",
               cmd->baryv);
        hibin = obs.N / 2;
        for (ii = 0; ii < numbirds; ii++) {
            if (bird_lobins[ii] >= hibin)
                break;
            if (bird_hibins[ii] >= hibin)
                bird_hibins[ii] = hibin - 1;
            zapbirds(bird_lobins[ii], bird_hibins[ii], obs.numbins, NULL, obs.fft, 0);
        }

        vect_free(bird_lobins);
        vect_free(bird_hibins);
    }

    printf("Searching with up to %d harmonics summed:\n",
           1 << (obs.numharmstages - 1));
    printf("  f = %.1f to %.1f Hz\n", obs.rlo / obs.T, obs.rhi / obs.T);
    printf("  r = %.1f to %.1f Fourier bins\n", obs.rlo, obs.rhi);
    printf("  z = %.1f to %.1f Fourier bins drifted\n", obs.zlo, obs.zhi);
    if (obs.numw)
        printf("  w = %.1f to %.1f Fourier-derivative bins drifted\n", obs.wlo, obs.whi);

    /* Generate the correlation kernels */

    printf("\nGenerating correlation kernels:\n");
    subharminfs = create_subharminfos(&obs);
    printf("Done generating kernels.\n\n");
    if (cmd->ncpus > 1) {
#ifdef _OPENMP
        set_openmp_numthreads(cmd->ncpus);
#endif
    } else {
#ifdef _OPENMP
        omp_set_num_threads(1); // Explicitly turn off OpenMP
#endif
        printf("Starting the search.\n\n");
    }
    if (cmd->perfinfoP)
        show_perfinfo();
    /* Don't use the *.txtcand files on short in-memory searches */
    if (!obs.dat_input) {
        printf("  Working candidates in a test format are in '%s'.\n\n",
               obs.workfilenm);
    }

    /* Function pointers to make code a bit cleaner */
    void (*fund_to_ffdot)(ffdotpows *, accelobs *) = NULL;
    void (*inmem_add_subharm)(ffdotpows *, accelobs *, int, int) = NULL;
    if (obs.inmem) {
        if (cmd->otheroptP) {
            fund_to_ffdot = &fund_to_ffdotplane_trans;
            inmem_add_subharm = &inmem_add_ffdotpows_trans;
        } else {
            fund_to_ffdot = &fund_to_ffdotplane;
            inmem_add_subharm = &inmem_add_ffdotpows;
        }
    }

    /* Start the main search loop.  Blocks of the spectrum are searched in
     * parallel, one block per thread, and each block's hits are added to
     * the candidate list in block order, so the results do not depend on
     * the number of threads. */
    {
        long numblocks = 0, numpreblocks = 0;
        double *starts = NULL, *prestarts = NULL;

        /* The blocks' starting r values, from the same sums as a serial
         * walk through the spectrum (obs.rlo need not be an integer) */
        {
            double startr;
            long maxblocks = (long) (fmax(obs.highestbin - obs.rlo, 0.0) / rstep) + 2;
            starts = (double *) malloc(maxblocks * sizeof(double));
            for (startr = obs.rlo; startr + rstep < obs.highestbin; startr += rstep)
                starts[numblocks++] = startr;
            if (obs.inmem) {
                /* The saved F-Fdot plane at low freqs, for in-memory
                 * searches of harmonics that are below obs.rlo */
                long maxpre = (long) (fmax(obs.rlo, 0.0) / rstep) + 2;
                prestarts = (double *) malloc(maxpre * sizeof(double));
                for (startr = 8; startr < obs.rlo; startr += rstep)  // a very low Fourier bin
                    prestarts[numpreblocks++] = startr;
            }
        }

        accel_make_plans(subharminfs, &obs);

        if (obs.inmem) {
            /* First fill the full in-core F-Fdot plane with the fundamentals.
             * The low-frequency blocks go first: the first search block
             * overwrites where they overlap it. */
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                accelwork *ws = create_accelwork(subharminfs, &obs);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (long blk = 0; blk < numpreblocks; blk++) {
                    double lastr = (prestarts[blk] + rstep) - ACCEL_DR;
                    fund_to_ffdot(subharm_fderivs_vol(1, 1, prestarts[blk], lastr,
                                                      &subharminfs[0][0], &obs, ws),
                                  &obs);
                }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (long blk = 0; blk < numblocks; blk++) {
                    double lastr = (starts[blk] + rstep) - ACCEL_DR;
                    fund_to_ffdot(subharm_fderivs_vol(1, 1, starts[blk], lastr,
                                                      &subharminfs[0][0], &obs, ws),
                                  &obs);
                }
                free_accelwork(ws);
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            accelwork *ws = create_accelwork(subharminfs, &obs);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1) ordered
#endif
            for (long blk = 0; blk < numblocks; blk++) {
                double startr = starts[blk];
                double lastr = (startr + rstep) - ACCEL_DR;
                ffdotpows *fundamental;

                /* Search the fundamental */
                if (obs.inmem)
                    fundamental = ffdotplane_to_fund(startr, lastr, &subharminfs[0][0],
                                                     &obs, ws, cmd->otheroptP);
                else
                    fundamental = subharm_fderivs_vol(1, 1, startr, lastr,
                                                      &subharminfs[0][0], &obs, ws);
                search_ffdotpows(fundamental, 1, &obs, ws);

                /* Search the subharmonics */
                for (int stage = 1; stage < obs.numharmstages; stage++) {
                    int harmtosum = 1 << stage;
                    if (obs.inmem) {
                        for (int harm = 1; harm < harmtosum; harm += 2)
                            inmem_add_subharm(fundamental, &obs, harmtosum, harm);
                        search_ffdotpows(fundamental, harmtosum, &obs, ws);
                    } else {
                        /* Compute all of the stage's subharmonics, then add
                         * them and search, row by row */
                        ffdotpows *subharmonics[ACCEL_MAXSTAGESUBS];
                        for (int harm = 1; harm < harmtosum; harm += 2)
                            subharmonics[harm / 2] =
                                subharm_fderivs_vol(harmtosum, harm, startr, lastr,
                                                    &subharminfs[stage][harm - 1],
                                                    &obs, ws);
                        add_search_ffdotpows(fundamental, subharmonics, harmtosum,
                                             &obs, ws);
                    }
                }
#ifdef _OPENMP
#pragma omp ordered
#endif
                {
                    cands = insert_accelhits(cands, ws, &obs);
                    print_percent_complete(startr - obs.rlo,
                                           obs.highestbin - obs.rlo, "search", 0);
                }
            }
            free_accelwork(ws);
        }
        accel_free_plans();
        free(starts);
        free(prestarts);
        print_percent_complete(obs.highestbin - obs.rlo,
                               obs.highestbin - obs.rlo, "search", 0);
    }

    printf("\n\nDone searching.  Now optimizing each candidate.\n\n");
    free_subharminfos(&obs, subharminfs);

    {                           /* Candidate list trimming and optimization */
        int numcands = g_slist_length(cands);
        GSList *listptr;
        accelcand *cand;
        fourierprops *props;

        if (numcands) {

            /* Sort the candidates according to the optimized sigmas */
            cands = sort_accelcands(cands);

            /* Eliminate (most of) the harmonically related candidates */
            if ((cmd->numharm > 1) && !(cmd->noharmremoveP))
                eliminate_harmonics(cands, &numcands);

            /* Now optimize each candidate and its harmonics */
            print_percent_complete(0, 0, NULL, 1);
            listptr = cands;
            for (ii = 0; ii < numcands; ii++) {
                print_percent_complete(ii, numcands, "optimization", 0);
                cand = (accelcand *) (listptr->data);
                optimize_accelcand(cand, &obs);
                listptr = listptr->next;
            }
            print_percent_complete(ii, numcands, "optimization", 0);

            /* Calculate the properties of the fundamentals */
            props = (fourierprops *) malloc(sizeof(fourierprops) * numcands);
            listptr = cands;
            for (ii = 0; ii < numcands; ii++) {
                cand = (accelcand *) (listptr->data);
                /* In case the fundamental harmonic is not significant,  */
                /* send the originally determined r and z from the       */
                /* harmonic sum in the search.  Note that the derivs are */
                /* not used for the computations with the fundamental.   */
                calc_props(cand->derivs[0], cand->r, cand->z, cand->w, props + ii);
                /* Override the error estimates based on power */
                props[ii].rerr = (float) (ACCEL_DR) / cand->numharm;
                props[ii].zerr = (float) (ACCEL_DZ) / cand->numharm;
                props[ii].werr = (float) (ACCEL_DW) / cand->numharm;
                listptr = listptr->next;
            }
            /* Write the fundamentals to the output text file */
            output_fundamentals(props, cands, &obs, &idata);

            /* Write the harmonics to the output text file */
            output_harmonics(cands, &obs, &idata);
            
            /* Write the fundamental fourierprops to the cand file */
            obs.workfile = chkfopen(obs.candnm, "wb");
            chkfwrite(props, sizeof(fourierprops), numcands, obs.workfile);
            fclose(obs.workfile);
            free(props);
            printf("\n\n");
        } else {
            printf("No candidates above sigma = %.2f were found.\n\n", obs.sigma);
        }
    }

    /* Finish up */

    printf("Searched the following approx numbers of independent points:\n");
    printf("  %d harmonic:   %9lld\n", 1, obs.numindep[0]);
    for (ii = 1; ii < obs.numharmstages; ii++)
        printf("  %d harmonics:  %9lld\n", 1 << ii, obs.numindep[ii]);

    printf("\nTiming summary:\n");
    tott = times(&runtimes) / clk_tck - tott;
    utim = runtimes.tms_utime / clk_tck;
    stim = runtimes.tms_stime / clk_tck;
    ttim = utim + stim;
    printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n",
           ttim, utim, stim);
    printf("  Total time: %.3f sec\n\n", tott);

    printf("Final candidates in binary format are in '%s'.\n", obs.candnm);
    printf("Final Candidates in a text format are in '%s'.\n\n", obs.accelnm);

    free_accelobs(&obs);
    g_slist_foreach(cands, free_accelcand, NULL);
    g_slist_free(cands);
    
    return (0);
}
