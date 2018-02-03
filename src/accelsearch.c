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
extern void zapbirds(double lobin, double hibin, FILE * fftfile, fcomplex * fft);

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

int main(int argc, char *argv[])
{
    int ii, rstep;
    double ttim, utim, stim, tott;
    struct tms runtimes;
    subharminfo **subharminfs;
    accelobs obs;
    infodata idata;
    GSList *cands = NULL;
    Cmdline *cmd;

    /* Prep the timer */

    tott = times(&runtimes) / (double) CLK_TCK;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

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
            zapbirds(bird_lobins[ii], bird_hibins[ii], NULL, obs.fft);
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
    /* Don't use the *.txtcand files on short in-memory searches */
    if (!obs.dat_input) {
        printf("  Working candidates in a test format are in '%s'.\n\n",
               obs.workfilenm);
    }

    /* Function pointers to make code a bit cleaner */
    void (*fund_to_ffdot)() = NULL;
    void (*add_subharm)() = NULL;
    void (*inmem_add_subharm)() = NULL;
    if (obs.inmem) {
        if (cmd->otheroptP) {
            fund_to_ffdot = &fund_to_ffdotplane_trans;
            inmem_add_subharm = &inmem_add_ffdotpows_trans;
        } else {
            fund_to_ffdot = &fund_to_ffdotplane;
            inmem_add_subharm = &inmem_add_ffdotpows;
        }
    } else {
        if (cmd->otheroptP) {
            add_subharm = &add_ffdotpows_ptrs;
        } else {
            add_subharm = &add_ffdotpows;
        }
    }

    /* Start the main search loop */
    {
        double startr, lastr, nextr;
        ffdotpows *fundamental;

        /* Populate the saved F-Fdot plane at low freqs for in-memory
         * searches of harmonics that are below obs.rlo */
        if (obs.inmem) {
            startr = 8;  // Choose a very low Fourier bin
            lastr = 0;
            nextr = 0;
            while (startr < obs.rlo) {
                nextr = startr + rstep;
                lastr = nextr - ACCEL_DR;
                // Compute the F-Fdot plane
                fundamental = subharm_fderivs_vol(1, 1, startr, lastr,
                                                  &subharminfs[0][0], &obs);
                // Copy it into the full in-core one
                fund_to_ffdot(fundamental, &obs);
                free_ffdotpows(fundamental);
                startr = nextr;
            }
        }    
        
        /* Reset indices if needed and search for real */
        startr = obs.rlo;
        lastr = 0;
        nextr = 0;
        while (startr + rstep < obs.highestbin) {
            /* Search the fundamental */
            print_percent_complete(startr - obs.rlo,
                                   obs.highestbin - obs.rlo, "search", 0);
            nextr = startr + rstep;
            lastr = nextr - ACCEL_DR;
            fundamental = subharm_fderivs_vol(1, 1, startr, lastr,
                                              &subharminfs[0][0], &obs);
            cands = search_ffdotpows(fundamental, 1, &obs, cands);

            if (obs.numharmstages > 1) {        /* Search the subharmonics */
                int stage, harmtosum, harm;
                ffdotpows *subharmonic;

                // Copy the fundamental's ffdot plane to the full in-core one
                if (obs.inmem)
                    fund_to_ffdot(fundamental, &obs);
                for (stage = 1; stage < obs.numharmstages; stage++) {
                    harmtosum = 1 << stage;
                    for (harm = 1; harm < harmtosum; harm += 2) {
                        if (obs.inmem) {
                            inmem_add_subharm(fundamental, &obs, harmtosum, harm);
                        } else {
                            subharmonic =
                                subharm_fderivs_vol(harmtosum, harm, startr, lastr,
                                                    &subharminfs[stage][harm - 1],
                                                    &obs);
                            add_subharm(fundamental, subharmonic, harmtosum, harm);
                            free_ffdotpows(subharmonic);
                        }
                    }
                    cands = search_ffdotpows(fundamental, harmtosum, &obs, cands);
                }
            }
            free_ffdotpows(fundamental);
            startr = nextr;
        }
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
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
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
