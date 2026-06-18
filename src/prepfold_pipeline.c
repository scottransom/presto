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


/* Phase 1: identify the data type, open the input file(s), read any */
/* mask, and (for events) read the event list.                        */
void identify_and_open_input(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                             foldcand *fc, mask *obsmask, plotflags *pflags,
                             char **rootnm_out, int *ptsperrec_out,
                             long long *numrec_out, int *numchan_out,
                             int *useshorts_out, int *insubs_out,
                             double **events_out, int *numevents_out, double *T_out)
{
    prepfoldinfo *search = &fc->search;
    char *rootnm = NULL;
    int ptsperrec = *ptsperrec_out;
    long long numrec = *numrec_out;
    int numchan = *numchan_out;
    int useshorts = *useshorts_out;
    int insubs = *insubs_out;
    int good_padvals = 0;
    double *events = NULL;
    int numevents = 0;
    double T = *T_out;
    long ii = 0;

    // Determine a output filename if necessary
    {
        char *path;

        split_path_file(cmd->argv[0], &path, &search->filenm);
        free(path);

        if (!cmd->outfileP) {
            char *tmprootnm, *suffix;
            split_root_suffix(cmd->argv[0], &tmprootnm, &suffix);
            if ((cmd->startT != 0.0) || (cmd->endT != 1.0)) {
                rootnm = (char *) calloc(strlen(tmprootnm) + 11, sizeof(char));
                sprintf(rootnm, "%s_%4.2f-%4.2f", tmprootnm, cmd->startT, cmd->endT);
                free(tmprootnm);   /* contents copied into rootnm above */
            } else {
                rootnm = tmprootnm;
            }
            free(suffix);
        } else {
            rootnm = cmd->outfile;
        }
    }

    if (RAWDATA) {
        if (cmd->filterbankP)
            s->datatype = SIGPROCFB;
        else if (cmd->psrfitsP)
            s->datatype = PSRFITS;
    } else {                    // Attempt to auto-identify the data
        identify_psrdatatype(s, 1);
        if (s->datatype == SIGPROCFB)
            cmd->filterbankP = 1;
        else if (s->datatype == PSRFITS)
            cmd->psrfitsP = 1;
        else if (s->datatype == EVENTS)
            cmd->eventsP = pflags->events = 1;
        else if (s->datatype == SDAT)
            useshorts = 1;
        else if (s->datatype == DAT)
            useshorts = 0;
        else if (s->datatype == SUBBAND) {
            useshorts = 1;
            insubs = 1;
        } else {
            printf
                ("Error:  Unable to identify input data files.  Please specify type.\n\n");
            exit(1);
        }
    }

    if (!RAWDATA)
        s->files = (FILE **) malloc(sizeof(FILE *) * s->num_files);
    if (RAWDATA || insubs) {
        char description[40];
        psrdatatype_description(description, s->datatype);
        if (s->num_files > 1)
            printf("Reading %s data from %d files:\n", description, s->num_files);
        else
            printf("Reading %s data from 1 file:\n", description);
        for (ii = 0; ii < s->num_files; ii++) {
            printf("  '%s'\n", cmd->argv[ii]);
            if (insubs)
                s->files[ii] = chkfopen(s->filenames[ii], "rb");
        }
        printf("\n");
        if (RAWDATA) {
            read_rawdata_files(s);
            if (cmd->ignorechanstrP) {
                s->ignorechans = get_ignorechans(cmd->ignorechanstr, 0, s->num_channels-1,
                                                &s->num_ignorechans, &s->ignorechans_str);
                if (s->ignorechans_str==NULL) {
                    s->ignorechans_str = (char *)malloc(strlen(cmd->ignorechanstr)+1);
                    strcpy(s->ignorechans_str, cmd->ignorechanstr);
                }
            }
            print_spectra_info_summary(s);
            spectra_info_to_inf(s, idata);
            ptsperrec = s->spectra_per_subint;
            numrec = s->N / ptsperrec;
            numchan = s->num_channels;
            if (!cmd->nsubP) {
                fc->nsub = 1;  // flag
                if (numchan <= 256) {
                    if (numchan % 32 == 0) {
                        fc->nsub = 32;
                    } else if (numchan % 30 == 0) {
                        fc->nsub = 30;
                    } else if (numchan % 25 == 0) {
                        fc->nsub = 25;
                    } else if (numchan % 20 == 0) {
                        fc->nsub = 20;
                    }
                } else if (numchan <= 1024) {
                    if (numchan % 8 == 0) {
                        fc->nsub = numchan / 8;
                    } else if (numchan % 10 == 0) {
                        fc->nsub = numchan / 10;
                    }
                } else {
                    if (numchan % 128 == 0) {
                        fc->nsub = 128;
                    } else if (numchan % 100 == 0) {
                        fc->nsub = 100;
                    }
                }
                if (fc->nsub == 1) {
                    perror("Cannot automatically determine a good value for -nsub");
                    printf("\n");
                    exit(1);
                }
            }
        } else {                // insubs
            fc->nsub = s->num_files;
            s->N = chkfilelen(s->files[0], sizeof(short));
            s->spectra_per_subint = ptsperrec = SUBSBLOCKLEN;
            numrec = s->N / ptsperrec;
            s->padvals = gen_fvect(s->num_files);
            for (ii = 0; ii < s->num_files; ii++)
                s->padvals[ii] = 0.0;
            s->start_MJD = (long double *) malloc(sizeof(long double));
            s->start_spec = (long long *) malloc(sizeof(long long));
            s->num_spec = (long long *) malloc(sizeof(long long));
            s->num_pad = (long long *) malloc(sizeof(long long));
            s->start_spec[0] = 0L;
            s->num_spec[0] = s->N;
            s->num_pad[0] = 0L;
        }
        /* Read an input mask if wanted */
        if (cmd->maskfileP) {
            read_mask(cmd->maskfile, obsmask);
            printf("Read mask information from '%s'\n\n", cmd->maskfile);
            if ((obsmask->numchan != idata->num_chan) ||
                (fabs(obsmask->mjd - (idata->mjd_i + idata->mjd_f)) > 1e-9)) {
                    printf("WARNING!: maskfile has different number of channels or start MJD than raw data! Exiting.\n\n");
                    exit(1);
            }
            good_padvals = determine_padvals(cmd->maskfile, obsmask, s->padvals);
        } else {
            obsmask->numchan = obsmask->numint = 0;
        }
    }

    if (!RAWDATA) {
        char *root, *suffix;
        if (split_root_suffix(s->filenames[0], &root, &suffix) == 0) {
            printf("Error:  The input filename (%s) must have a suffix!\n\n",
                   s->filenames[0]);
            exit(1);
        }
        if (insubs) {
            char *tmpname;
            if (strncmp(suffix, "sub", 3) == 0) {
                tmpname = calloc(strlen(root) + 10, 1);
                sprintf(tmpname, "%s.sub", root);
                readinf(idata, tmpname);
                free(tmpname);
                s->num_channels = numchan = idata->num_chan;
                s->start_MJD[0] = idata->mjd_i + idata->mjd_f;
                s->dt = idata->dt;
                s->T = s->N * s->dt;
                s->lo_freq = idata->freq;
                s->df = idata->chan_wid;
                s->hi_freq = s->lo_freq + (s->num_channels - 1.0) * s->df;
                s->BW = s->num_channels * s->df;
                s->fctr = s->lo_freq - 0.5 * s->df + 0.5 * s->BW;
                s->padvals = gen_fvect(s->num_channels);
                for (ii = 0; ii < s->num_channels; ii++)
                    s->padvals[ii] = 0.0;
                print_spectra_info_summary(s);
            } else {
                printf
                    ("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                     cmd->argv[0]);
                exit(1);
            }
            /* Read an input mask if wanted */
            if (cmd->maskfileP) {
                read_mask(cmd->maskfile, obsmask);
                printf("Read mask information from '%s'\n\n", cmd->maskfile);
                if ((obsmask->numchan != idata->num_chan) ||
                    (fabs(obsmask->mjd - (idata->mjd_i + idata->mjd_f)) > 1e-9)) {
                        printf("WARNING!: maskfile has different number of channels or start MJD than raw data! Exiting.\n\n");
                        exit(1);
                }
                good_padvals = determine_padvals(cmd->maskfile, obsmask, s->padvals);
            } else {
                obsmask->numchan = obsmask->numint = 0;
            }
        } else {
            printf("Reading input data from '%s'.\n", cmd->argv[0]);
            printf("Reading information from '%s.inf'.\n\n", root);
            /* Read the info file if available */
            readinf(idata, root);
            fc->nsub = 1;
        }
        free(root);
        free(suffix);
        /* Use events instead of a time series */
        if (cmd->eventsP) {
            int eventtype = 0;  /* 0=sec since .inf, 1=days since .inf, 2=MJDs */

            /* The following allows using inf files from searches of a subset */
            /* of events from an event file.                                  */
            if (cmd->accelcandP) {
                infodata rzwidata;
                char *cptr = NULL;

                if (!cmd->accelfileP) {
                    printf("\nYou must enter a name for the ACCEL candidate ");
                    printf("file (-accelfile filename)\n");
                    printf("Exiting.\n\n");
                    exit(1);
                } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))) {
                    ii = (long) (cptr - cmd->accelfile);
                }
                cptr = (char *) calloc(ii + 1, sizeof(char));
                strncpy(cptr, cmd->accelfile, ii);
                readinf(&rzwidata, cptr);
                free(cptr);
                idata->mjd_i = rzwidata.mjd_i;
                idata->mjd_f = rzwidata.mjd_f;
                idata->N = rzwidata.N;
                idata->dt = rzwidata.dt;
            }
            if (!cmd->proflenP && search->proflen == 0) {
                search->proflen = 20;
                printf("Using %d bins in the profile since not specified.\n",
                       search->proflen);
            }
            if (cmd->doubleP)
                s->files[0] = chkfopen(s->filenames[0], "rb");
            else
                s->files[0] = chkfopen(s->filenames[0], "r");
            if (cmd->daysP)
                eventtype = 1;
            else if (cmd->mjdsP)
                eventtype = 2;
            events = read_events(s->files[0], cmd->doubleP, eventtype, &numevents,
                                 idata->mjd_i + idata->mjd_f, idata->N * idata->dt,
                                 cmd->startT, cmd->endT, cmd->offset);
            if (cmd->accelcandP)
                T = idata->N * idata->dt;
            else {
                /* The 1e-8 prevents floating point rounding issues from
                   causing the last event to go into slice npart+1.
                   Thanks to Paul Ray for finding this. */
                T = events[numevents - 1] + 1e-8;
            }
            printf("Found %d events starting at MJD %.15f over %.2f seconds\n",
                numevents, idata->mjd_i + idata->mjd_f, T);
        } else {
            if (!insubs)
                s->files[0] = chkfopen(s->filenames[0], "rb");
        }
    }
    (void) good_padvals;

    *rootnm_out = rootnm;
    *ptsperrec_out = ptsperrec;
    *numrec_out = numrec;
    *numchan_out = numchan;
    *useshorts_out = useshorts;
    *insubs_out = insubs;
    *events_out = events;
    *numevents_out = numevents;
    *T_out = T;
}


/* Phase 3: build the candidate name and the output/plot file names. */
void resolve_output_names(Cmdline *cmd, foldcand *fc, char *rootnm)
{
    prepfoldinfo *search = &fc->search;
    char *outfilenm, *plotfilenm;
    int slen;

    /* Determine the candidate name */

    if (cmd->psrnameP) {
        search->candnm = (char *) calloc(strlen(cmd->psrname) + 5, sizeof(char));
        sprintf(search->candnm, "PSR_%s", cmd->psrname);
    } else if (cmd->parnameP || cmd->timingP) {
        psrparams psr;
        if (cmd->timingP) {
            cmd->parnameP = 1;
            cmd->parname = cmd->timing;
        }
        /* Read the par file just to get the PSR name */
        get_psr_from_parfile(cmd->parname, 51000.0, &psr);
        search->candnm = (char *) calloc(strlen(psr.jname) + 5, sizeof(char));
        sprintf(search->candnm, "PSR_%s", psr.jname);
    } else if (cmd->accelcandP) {
        char *cptr = NULL;
        slen = 22;
        search->candnm = (char *) calloc(slen, sizeof(char));
        if (NULL != (cptr = strstr(cmd->accelfile, "_JERK")))
            sprintf(search->candnm, "JERK_Cand_%d", cmd->accelcand);
        else
            sprintf(search->candnm, "ACCEL_Cand_%d", cmd->accelcand);
    } else {
        slen = 20;
        search->candnm = (char *) calloc(slen, sizeof(char));
        if (cmd->pP)
            sprintf(search->candnm, "%.2fms_Cand", cmd->p * 1000.0);
        else if (cmd->fP)
            sprintf(search->candnm, "%.2fHz_Cand", cmd->f);
        else {
            printf("\nYou must specify candidate parameters (i.e period).\n\n");
            exit(1);
        }
    }

    /* Determine the output and plot file names */

    slen = strlen(rootnm) + strlen(search->candnm) + 6;
    outfilenm = (char *) calloc(slen, sizeof(char));
    sprintf(outfilenm, "%s_%s.pfd", rootnm, search->candnm);
    plotfilenm = (char *) calloc(slen + 3, sizeof(char));
    sprintf(plotfilenm, "%s_%s.pfd.ps", rootnm, search->candnm);
    search->pgdev = (char *) calloc(slen + 7, sizeof(char));
    sprintf(search->pgdev, "%s/CPS", plotfilenm);

    fc->outfilenm = outfilenm;
    fc->plotfilenm = plotfilenm;
}


/* Phase 2: compute the observation timing (T, N, records, epochs), check */
/* the nsub divisibility, and set the source position.                    */
void compute_obs_timing(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                        foldcand *fc, mask *obsmask, int insubs,
                        int useshorts, int ptsperrec, char *obs, char *ephem,
                        char *rastring, char *decstring, int *numchan_out,
                        int **maskchans_out, double *recdt_out, long long *lorec_out,
                        double *startTday_out, long long *numrec_out,
                        long *reads_per_part_out, double *T_out, double *N_out,
                        long *worklen_out)
{
    prepfoldinfo *search = &fc->search;
    int numchan = *numchan_out;
    int *maskchans = *maskchans_out;
    double recdt = *recdt_out;
    long long lorec = *lorec_out, hirec = 0;
    double startTday = *startTday_out;
    long long numrec = *numrec_out;
    long reads_per_part = *reads_per_part_out;
    double T = *T_out, N = *N_out;
    long worklen = *worklen_out;
    double barydispdt, dtmp;

    /* What ephemeris will we use?  (Default is DE405) */
    strcpy(ephem, "DE405");

    // Set-up values if we are using raw radio pulsar data

    if (RAWDATA || insubs) {

        // Identify the TEMPO observatory code
        search->telescope = (char *) calloc(40, sizeof(char));
        telescope_to_tempocode(idata->telescope, search->telescope, obs);

        idata->dm = fc->dm;
        if (cmd->maskfileP)
            maskchans = gen_ivect(obsmask->numchan);

        /* Define the RA and DEC of the observation */

        ra_dec_to_string(rastring, idata->ra_h, idata->ra_m, idata->ra_s);
        ra_dec_to_string(decstring, idata->dec_d, idata->dec_m, idata->dec_s);

        /* Define some variables */

        search->dt = idata->dt;
        recdt = search->dt * ptsperrec;

        /* Determine the number of records to use from the command line */

        search->startT = cmd->startT;
        search->endT = cmd->endT;
        lorec = (long) (cmd->startT * numrec + DBLCORRECT);
        hirec = (long) (cmd->endT * numrec + DBLCORRECT);
        startTday = lorec * recdt / SECPERDAY;
        numrec = hirec - lorec;

        /* The number of reads from the file we need for */
        /* each sub-integration.                         */

        reads_per_part = numrec / fc->npart;

        /* If the number of records is less than the number of parts    */
        /* then set the number of parts equat to the number of records. */
        if (numrec < fc->npart) {
            reads_per_part = 1;
            fc->npart = numrec;
            printf
                ("Overriding -npart to be %lld, the number of raw (requested) records.\n",
                 numrec);
        }

        /* Correct numrec so that each part will contain */
        /* the same number of records.                   */

        numrec = reads_per_part * fc->npart;
        T = numrec * recdt;
        N = numrec * ptsperrec;

        /* Topocentric and barycentric times of folding epoch data */

        if (idata->mjd_i) {
            search->tepoch = idata->mjd_i + idata->mjd_f + startTday;

            if (!cmd->polycofileP && !cmd->timingP && !cmd->topoP && !cmd->parnameP) {
                barycenter(&search->tepoch, &search->bepoch, &dtmp, 1, rastring,
                           decstring, obs, ephem);

                /* Correct the barycentric time for the dispersion delay.     */
                /* This converts the barycentric time to infinite frequency.  */
                if (fc->dm > 0.0) {
                    barydispdt = delay_from_dm(fc->dm, idata->freq +
                                               (idata->num_chan -
                                                1) * idata->chan_wid);
                    search->bepoch -= (barydispdt / SECPERDAY);
                }
            }
        }
        worklen = ptsperrec;

    } else {                    /* Raw floating point or event data (already de-dispersed if radio data) */

        fc->nsub = 1;
        search->startT = cmd->startT;
        search->endT = cmd->endT;

        if (!cmd->eventsP) {

            /* Some information about the size of the records */

            numchan = 1;
            worklen = SUBSBLOCKLEN;
            search->dt = idata->dt;
            if (useshorts)
                N = chkfilelen(s->files[0], sizeof(short));
            else
                N = chkfilelen(s->files[0], sizeof(float));

            /* Determine the number of records to use from the command line */

            lorec = (long) (cmd->startT * N + DBLCORRECT);
            hirec = (long) (cmd->endT * N + DBLCORRECT);
            startTday = lorec * search->dt / SECPERDAY;
            numrec = (hirec - lorec) / worklen;
            recdt = worklen * search->dt;

            /* The number of reads from the file we need for */
            /* each sub-integration.                         */

            reads_per_part = numrec / fc->npart;

            /* Correct numrec so that each part will contain */
            /* the same number of records.                   */

            numrec = reads_per_part * fc->npart;
            N = numrec * worklen;
            T = N * search->dt;
        }

        /* Until I figure out a better way to do this... */

        search->telescope =
            (char *) calloc(strlen(idata->telescope) + 1, sizeof(char));
        strcpy(search->telescope, idata->telescope);

        if (idata->mjd_i) {
            if (idata->bary)
                search->bepoch = idata->mjd_i + idata->mjd_f + startTday;
            else
                search->tepoch = idata->mjd_i + idata->mjd_f + startTday;
        }
    }

    /* Make sure that the number of subbands evenly divides the number of channels */
    if (numchan % fc->nsub != 0) {
        printf("Error:  # of channels (%d) not divisible by # of subbands (%d)!\n",
               numchan, fc->nsub);
        exit(1);
    }

    set_posn(search, idata);
    printf("Folding a %s candidate.\n\n", search->candnm);
    printf("Output data file is '%s'.\n", fc->outfilenm);
    printf("Output plot file is '%s'.\n", fc->plotfilenm);
    printf("Best profile is in  '%s.bestprof'.\n", fc->outfilenm);

    *numchan_out = numchan;
    *maskchans_out = maskchans;
    *recdt_out = recdt;
    *lorec_out = lorec;
    *startTday_out = startTday;
    *numrec_out = numrec;
    *reads_per_part_out = reads_per_part;
    *T_out = T;
    *N_out = N;
    *worklen_out = worklen;
}


/* Phase 4+5: resolve the fold parameters (polycos / psr database / par file / */
/* accelcand / binary / explicit p-f), apply -pfact/-ffact, and determine the  */
/* profile length.                                                             */
void resolve_fold_params(Cmdline *cmd, infodata *idata, foldcand *fc,
                         int insubs, double T, double startTday, double recdt,
                         long long lorec, char *pname)
{
    prepfoldinfo *search = &fc->search;
    double f = fc->f, fd = fc->fd, fdd = fc->fdd;
    int binary = fc->binary;
    int polyco_index = fc->polyco_index;
    double polyco_phase0 = fc->polyco_phase0;
    double dtmp;
    long ii = 0;

    /* Generate polycos if required and set the pulsar name */
    if (((cmd->timingP || cmd->parnameP) && (!idata->bary)) ||
        (idata->bary && cmd->barypolycosP)) {
        char *polycofilenm;
        polycofilenm = (char *) calloc(strlen(fc->outfilenm) + 9, sizeof(char));
        sprintf(polycofilenm, "%s.polycos", fc->outfilenm);
        cmd->psrnameP = 1;
        if (cmd->timingP)
            cmd->psrname = make_polycos(cmd->timing, idata, polycofilenm, cmd->debugP);
        else
            cmd->psrname = make_polycos(cmd->parname, idata, polycofilenm, cmd->debugP);
        printf("Polycos used are in '%s'.\n", polycofilenm);
        cmd->polycofileP = 1;
        cmd->polycofile = (char *) calloc(strlen(polycofilenm) + 1, sizeof(char));
        strcpy(cmd->polycofile, polycofilenm);
        free(polycofilenm);
    }

    /* Read the pulsar database if needed */
    if (cmd->psrnameP) {
        if (cmd->polycofileP) {
            FILE *polycofileptr;
            int numsets;
            double polyco_dm, epoch = search->tepoch;

            if (idata->bary)
                epoch = search->bepoch;
            polycofileptr = chkfopen(cmd->polycofile, "r");
            numsets = getpoly(epoch, T / SECPERDAY, &polyco_dm,
                              polycofileptr, cmd->psrname);
            fclose(polycofileptr);
            if (fc->dm > 0.0) {
                printf("\nRead %d set(s) of polycos for PSR %s at %18.12f\n",
                       numsets, cmd->psrname, epoch);
                printf("Overriding polyco DM = %f with %f\n", polyco_dm, fc->dm);
            } else {
                printf
                    ("\nRead %d set(s) of polycos for PSR %s at %18.12f (DM = %.5g)\n",
                     numsets, cmd->psrname, epoch, polyco_dm);
                fc->dm = polyco_dm;
            }
            polyco_index = phcalc(idata->mjd_i, idata->mjd_f + startTday,
                                  polyco_index, &polyco_phase0, &f);
            search->topo.p1 = 1.0 / f;
            search->topo.p2 = fd = 0.0;
            search->topo.p3 = fdd = 0.0;
            if (idata->bary) {
                search->bary.p1 = search->topo.p1;
                search->bary.p2 = search->topo.p2;
                search->bary.p3 = search->topo.p3;
            }
            strcpy(pname, cmd->psrname);
        } else {                /* Use the database */
            int pnum;
            psrparams psr;

            if (search->bepoch == 0.0) {
                printf
                    ("\nYou cannot fold topocentric data with the pulsar database.\n");
                printf("Use '-timing' or polycos instead.  Exiting.\n\n");
                exit(1);
            }
            pnum = get_psr_at_epoch(cmd->psrname, search->bepoch, &psr);
            if (!pnum) {
                printf("The pulsar is not in the database.  Exiting.\n\n");
                exit(1);
            }
            if (psr.orb.p != 0.0) {     /* Checks if the pulsar is in a binary */
                binary = 1;
                search->orb = psr.orb;
            }
            search->bary.p1 = psr.p;
            search->bary.p2 = psr.pd;
            search->bary.p3 = psr.pdd;
            if (fc->dm == 0.0)
                fc->dm = psr.dm;
            f = psr.f;
            fd = psr.fd;
            fdd = psr.fdd;
            strcpy(pname, psr.jname);
        }

    } else if (cmd->parnameP) { /* Read ephemeris from a par file */
        psrparams psr;

        if (search->bepoch == 0.0) {
            printf("\nYou cannot fold topocentric data with a par file.\n");
            printf("Use '-timing' or polycos instead.  Exiting.\n\n");
            exit(1);
        }

        /* Read the par file */
        get_psr_from_parfile(cmd->parname, search->bepoch, &psr);

        if (psr.orb.p != 0.0) { /* Checks if the pulsar is in a binary */
            binary = 1;
            search->orb = psr.orb;
        }
        search->bary.p1 = psr.p;
        search->bary.p2 = psr.pd;
        search->bary.p3 = psr.pdd;
        f = psr.f;
        fd = psr.fd;
        fdd = psr.fdd;
        strcpy(pname, psr.jname);
        if (fc->dm == 0.0)
            fc->dm = psr.dm;

        /* If the user specifies all of the binaries parameters */

    } else if (cmd->binaryP) {

        /* Assume that the psr characteristics were measured at the time */
        /* of periastron passage (cmd->To)                               */

        if (search->bepoch == 0.0)
            dtmp = SECPERDAY * (search->tepoch - cmd->To);
        else
            dtmp = SECPERDAY * (search->bepoch - cmd->To);
        search->orb.p = cmd->pb;
        search->orb.x = cmd->asinic;
        search->orb.e = cmd->e;
        search->orb.t = fmod(dtmp, search->orb.p);
        if (search->orb.t < 0.0)
            search->orb.t += search->orb.p;
        search->orb.w = (cmd->w + dtmp * cmd->wdot / SECPERJULYR);
        binary = 1;

    } else if (cmd->accelcandP) {
        fourierprops rzwcand;
        infodata rzwidata;
        char *cptr = NULL;
        double T, r0, z0;

        if (!cmd->accelfileP) {
            printf("\nYou must enter a name for the ACCEL candidate ");
            printf("file (-accelfile filename)\n");
            printf("Exiting.\n\n");
            exit(1);
        } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))) {
            ii = (long) (cptr - cmd->accelfile);
        }
        cptr = (char *) calloc(ii + 1, sizeof(char));
        strncpy(cptr, cmd->accelfile, ii);
        printf("\nAttempting to read '%s.inf'.  ", cptr);
        readinf(&rzwidata, cptr);
        free(cptr);
        printf("Successful.\n");
        get_rzw_cand(cmd->accelfile, cmd->accelcand, &rzwcand);
        T = rzwidata.dt * rzwidata.N;
        // fourier props file reports average r and average z.
        // We need the starting values.
        z0 = rzwcand.z - 0.5 * rzwcand.w;
        r0 = rzwcand.r - 0.5 * z0 - rzwcand.w / 6.0;
        f = r0 / T;
        fd = z0 / (T * T);
        fdd = rzwcand.w / (T * T * T);

        /* Now correct for the fact that we may not be starting */
        /* to fold at the same start time as the rzw search.    */

        if (RAWDATA || insubs)
            f += lorec * recdt * fd;
        else
            f += lorec * search->dt * fd;
        if (rzwidata.bary)
            switch_f_and_p(f, fd, fdd, &search->bary.p1,
                           &search->bary.p2, &search->bary.p3);
        else
            switch_f_and_p(f, fd, fdd, &search->topo.p1,
                           &search->topo.p2, &search->topo.p3);
    }

    /* Determine the pulsar parameters to fold if we are not getting   */
    /* the data from a .cand file, the pulsar database, or a makefile. */

    if (!cmd->accelcandP && !cmd->psrnameP && !cmd->parnameP) {
        double p = 0.0, pd = 0.0, pdd = 0.0;

        if (cmd->pP) {
            p = cmd->p;
            f = 1.0 / p;
        }
        if (cmd->fP) {
            f = cmd->f;
            p = 1.0 / f;
        }
        if (cmd->pd != 0.0) {
            pd = cmd->pd;
            fd = -pd / (p * p);
        }
        if (cmd->fd != 0.0) {
            fd = cmd->fd;
            pd = -fd / (f * f);
        }
        if (cmd->pdd != 0.0) {
            pdd = cmd->pdd;
            fdd = 2 * pd * pd / (p * p * p) - pdd / (p * p);
        }
        if (cmd->fdd != 0.0) {
            fdd = cmd->fdd;
            pdd = 2 * fd * fd / (f * f * f) - fdd / (f * f);
        }
        if (idata->bary) {
            search->bary.p1 = p;
            search->bary.p2 = pd;
            search->bary.p3 = pdd;
        } else {
            search->topo.p1 = p;
            search->topo.p2 = pd;
            search->topo.p3 = pdd;
        }
    }

    /* Modify the fold period or frequency */

    if (cmd->pfact != 1.0 || cmd->ffact != 1.0) {
        if (cmd->pfact == 0.0 || cmd->ffact == 0.0) {
            printf("\nFolding factors (-pfact or -ffact) cannot be 0!  Exiting\n");
            exit(1);
        }
        if (cmd->pfact != 1.0)
            cmd->ffact = 1.0 / cmd->pfact;
        else if (cmd->ffact != 1.0)
            cmd->pfact = 1.0 / cmd->ffact;
        f *= cmd->ffact;
        fd *= cmd->ffact;
        fdd *= cmd->ffact;
        if (idata->bary) {
            switch_f_and_p(f, fd, fdd,\
                &search->bary.p1, &search->bary.p2, &search->bary.p3);
        } else {
            switch_f_and_p(f, fd, fdd,\
                &search->topo.p1, &search->topo.p2, &search->topo.p3);
        }
    }

    /* Determine the length of the profile */

    if (cmd->proflenP) {
        search->proflen = cmd->proflen;
    } else if (search->proflen == 0) {
        if (search->topo.p1 == 0.0)
            search->proflen = (long) (search->bary.p1 / search->dt + 0.5);
        else
            search->proflen = (long) (search->topo.p1 / search->dt + 0.5);
        if (cmd->timingP)
            search->proflen = next2_to_n(search->proflen);
        if (search->proflen > 64)
            search->proflen = 64;
    }

    fc->f = f;
    fc->fd = fd;
    fc->fdd = fdd;
    fc->binary = binary;
    fc->polyco_index = polyco_index;
    fc->polyco_phase0 = polyco_phase0;
}


/* Phase 6: compute the orbital phase-delay tables (binary), apply the      */
/* barycentric-epoch dispersion correction, and (for events) set the dt/N.  */
void setup_orbit_delays(Cmdline *cmd, infodata *idata, foldcand *fc,
                        int insubs, double T, double *N_out)
{
    prepfoldinfo *search = &fc->search;
    double *Ep = fc->Ep, *tp = fc->tp;
    int binary = fc->binary;
    long long numbinpoints = fc->numbinpoints;
    int numdelays = fc->numdelays;
    double N = *N_out;
    double barydispdt;
    long ii;

    /* Determine the phase delays caused by the orbit if needed */

    if (binary && !cmd->eventsP) {
        double orbdt = 1.0, startE = 0.0;

        /* Save the orbital solution every half second               */
        /* The times in *tp are now calculated as barycentric times. */
        /* Later, we will change them to topocentric times after     */
        /* applying corrections to Ep using TEMPO.                   */

        startE = keplers_eqn(search->orb.t, search->orb.p, search->orb.e, 1.0E-15);
        if (T > 2048)
            orbdt = 0.5;
        else
            orbdt = T / 4096.0;
        numbinpoints = (long) floor(T / orbdt + 0.5) + 1;
        Ep = dorbint(startE, numbinpoints, orbdt, &search->orb);
        tp = gen_dvect(numbinpoints);
        for (ii = 0; ii < numbinpoints; ii++)
            tp[ii] = ii * orbdt;

        /* Convert Eccentric anomaly to time delays */

        E_to_phib(Ep, numbinpoints, &search->orb);
        numdelays = numbinpoints;
        if (search->bepoch == 0.0)
            search->orb.t = -search->orb.t / SECPERDAY + search->tepoch;
        else
            search->orb.t = -search->orb.t / SECPERDAY + search->bepoch;
    }

    if (((RAWDATA || insubs) && !cmd->topoP)
        && fc->dm == 0.0 && !cmd->polycofileP) {
        /* Correct the barycentric time for the dispersion delay.     */
        /* This converts the barycentric time to infinite frequency.  */
        barydispdt = delay_from_dm(fc->dm, idata->freq +
                                   (idata->num_chan - 1) * idata->chan_wid);
        search->bepoch -= (barydispdt / SECPERDAY);
    }

    if (cmd->eventsP) {
        if (search->topo.p1 == 0.0)
            search->dt = (search->bary.p1 + 0.5 * T * search->bary.p2) / search->proflen;
        else
            search->dt = (search->topo.p1 + 0.5 * T * search->topo.p2) / search->proflen;
        N = ceil(T / search->dt);
    }

    fc->Ep = Ep;
    fc->tp = tp;
    fc->numbinpoints = numbinpoints;
    fc->numdelays = numdelays;
    *N_out = N;
}


/* Phase 7: print the informational fold parameters and the sub-interval */
/* warning (which may override npart for events).                        */
void print_fold_info(Cmdline *cmd, foldcand *fc, double N, double T, char *pname)
{
    prepfoldinfo *search = &fc->search;
    double f = fc->f, fd = fc->fd, fdd = fc->fdd;
    int binary = fc->binary;
    FILE *filemarker;
    long ii;

    /* Output some informational data on the screen and to the */
    /* output file.                                            */

    fprintf(stdout, "\n");
    filemarker = stdout;
    for (ii = 0; ii < 1; ii++) {
        double p, pd, pdd;

        if (cmd->psrnameP)
            fprintf(filemarker, "Pulsar                       =  %s\n", pname);
        if (search->tepoch != 0.0)
            fprintf(filemarker,
                    "Folding (topo) epoch  (MJD)  =  %-17.12f\n", search->tepoch);
        if (search->bepoch != 0.0)
            fprintf(filemarker,
                    "Folding (bary) epoch  (MJD)  =  %-17.12f\n", search->bepoch);
        fprintf(filemarker, "Data pt duration (dt)   (s)  =  %-.12g\n", search->dt);
        fprintf(filemarker, "Total number of data points  =  %-.0f\n", N);
        fprintf(filemarker, "Number of profile bins       =  %d\n", search->proflen);
        switch_f_and_p(f, fd, fdd, &p, &pd, &pdd);
        fprintf(filemarker, "Folding period          (s)  =  %-.15g\n", p);
        if (pd != 0.0)
            fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10g\n", pd);
        if (pdd != 0.0)
            fprintf(filemarker, "Folding p-dotdot    (s/s^2)  =  %-.10g\n", pdd);
        fprintf(filemarker, "Folding frequency      (hz)  =  %-.12g\n", f);
        if (fd != 0.0)
            fprintf(filemarker, "Folding f-dot        (hz/s)  =  %-.8g\n", fd);
        if (fdd != 0.0)
            fprintf(filemarker, "Folding f-dotdot   (hz/s^2)  =  %-.8g\n", fdd);
        if (binary) {
            double tmpTo;
            fprintf(filemarker,
                    "Orbital period          (s)  =  %-.10g\n", search->orb.p);
            fprintf(filemarker,
                    "a*sin(i)/c (x)     (lt-sec)  =  %-.10g\n", search->orb.x);
            fprintf(filemarker,
                    "Eccentricity                 =  %-.10g\n", search->orb.e);
            fprintf(filemarker,
                    "Longitude of peri (w) (deg)  =  %-.10g\n",
                    search->orb.w);
            tmpTo = search->orb.t;
            if (cmd->eventsP) {
                if (search->bepoch == 0.0)
                    tmpTo = -search->orb.t / SECPERDAY + search->tepoch;
                else
                    tmpTo = -search->orb.t / SECPERDAY + search->bepoch;
            }
            fprintf(filemarker, "Epoch of periapsis    (MJD)  =  %-17.11f\n", tmpTo);
        }
    }

    /* Check to see if the folding period is close to the time in a sub-interval */
    if (!cmd->eventsP) {
        if (T / fc->npart < 5.0 / f)
            printf
                ("\nWARNING:  The pulse period is close to the duration of a folding\n"
                 "  sub-interval.  This may cause artifacts in the plot and/or excessive\n"
                 "  loss of data during the fold.  I recommend re-running with -npart set\n"
                 "  to a significantly smaller value than the current value of %d.\n",
                 fc->npart);
    } else {
        if (T / fc->npart < 5.0 / f) {
            fc->npart = (int) (T * f / 5.0);
            printf("\nOverriding default number of sub-intervals due to low number\n"
                   "  of pulses during the observation.  Current npart = %d\n",
                   fc->npart);
        }
    }
}


/* Phase 8: allocate and zero the rawfolds/stats arrays. */
void allocate_fold_arrays(Cmdline *cmd, foldcand *fc)
{
    prepfoldinfo *search = &fc->search;
    int numdelays = fc->numdelays;
    int flags = fc->flags;
    long ii;

    /* Allocate and initialize some arrays and other information */

    search->nsub = fc->nsub;
    search->npart = fc->npart;
    search->rawfolds = gen_dvect(fc->nsub * fc->npart * search->proflen);
    search->stats = (foldstats *) malloc(sizeof(foldstats) * fc->nsub * fc->npart);
    for (ii = 0; ii < fc->npart * fc->nsub * search->proflen; ii++)
        search->rawfolds[ii] = 0.0;
    for (ii = 0; ii < fc->npart * fc->nsub; ii++) {
        search->stats[ii].numdata = 0.0;
        search->stats[ii].data_avg = 0.0;
        search->stats[ii].data_var = 0.0;
    }
    if (numdelays == 0)
        flags = 0;

    fc->flags = flags;
}


/* Phase 9: fold an event list (instead of a time series). */
void fold_events(Cmdline *cmd, infodata *idata, foldcand *fc, double T,
                 int numevents, double *events, double startTday)
{
    prepfoldinfo *search = &fc->search;
    double f = fc->f, fd = fc->fd, fdd = fc->fdd;
    int binary = fc->binary;
    double foldf = fc->foldf, foldfd = fc->foldfd, foldfdd = fc->foldfdd;
    double orig_foldf = fc->orig_foldf;
    int polyco_index = fc->polyco_index;
    double *parttimes = fc->parttimes;
    long ii, jj;

    {
        double event, dtmp, cts, phase, begphs, endphs, dphs, lphs, rphs;
        double tf, tfd, tfdd, totalphs, calctotalphs, numwraps;
        int partnum, binnum;

        foldf = f;
        foldfd = fd;
        foldfdd = fdd;
        search->fold.pow = 1.0;  // Data are barycentric
        search->fold.p1 = f;
        search->fold.p2 = fd;
        search->fold.p3 = fdd;
        tf = f;
        tfd = fd / 2.0;
        tfdd = fdd / 6.0;
        dtmp = fc->npart / T;
        parttimes = gen_dvect(fc->npart);
        for (ii = 0; ii < numevents; ii++) {
            event = events[ii];
            if (binary) {
                double tt, delay;
                tt = fmod(search->orb.t + event, search->orb.p);
                delay = keplers_eqn(tt, search->orb.p, search->orb.e, 1.0E-15);
                E_to_phib(&delay, 1, &search->orb);
                event -= delay;
            }
            partnum = (int) floor(event * dtmp);
            if (!cmd->polycofileP)
                phase = event * (event * (event * tfdd + tfd) + tf);
            else {
                double mjdf = idata->mjd_f + startTday + event / SECPERDAY;
                /* Calculate the pulse phase for the event  */
                polyco_index =
                    phcalc(idata->mjd_i, mjdf, polyco_index, &phase, &foldf);
                if (ii == 0)
                    orig_foldf = foldf;
            }
            binnum = (int) ((phase - (long long) phase) * search->proflen);
            search->rawfolds[partnum * search->proflen + binnum] += 1.0;
        }
        if (binary) {
            if (search->bepoch == 0.0)
                search->orb.t = -search->orb.t / SECPERDAY + search->tepoch;
            else
                search->orb.t = -search->orb.t / SECPERDAY + search->bepoch;
        }
        for (ii = 0; ii < fc->npart; ii++) {
            parttimes[ii] = (T * ii) / (double) (fc->npart);
            /* Correct each part for the "exposure".  This gives us a count rate. */
            event = parttimes[ii];
            begphs = event * (event * (event * tfdd + tfd) + tf);
            event = (T * (ii + 1)) / (double) (fc->npart);
            endphs = event * (event * (event * tfdd + tfd) + tf);
            totalphs = endphs - begphs;
            begphs = begphs < 0.0 ? fmod(begphs, 1.0) + 1.0 : fmod(begphs, 1.0);
            endphs = endphs < 0.0 ? fmod(endphs, 1.0) + 1.0 : fmod(endphs, 1.0);
            dphs = 1.0 / search->proflen;
            /* printf("%.4f %.4f %5.4f\n", begphs, endphs, totalphs); */
            calctotalphs = 0.0;
            for (jj = 0; jj < search->proflen; jj++) {
                numwraps = floor(totalphs);
                lphs = jj * dphs;
                rphs = (jj + 1) * dphs;
                if (begphs <= lphs) {   /* BLR */
                    if (rphs <= endphs) {       /* LRE */
                        numwraps += 1.0;
                    } else if (endphs <= lphs) {        /* ELR */
                        if (endphs <= begphs)
                            numwraps += 1.0;
                    } else {    /* LER */
                        numwraps += (endphs - lphs) * search->proflen;
                    }
                } else if (rphs <= begphs) {    /* LRB */
                    if (rphs <= endphs) {       /* LRE */
                        if (endphs <= begphs)
                            numwraps += 1.0;
                    } else if (lphs <= endphs && endphs <= rphs) {      /* LER */
                        numwraps += (endphs - lphs) * search->proflen;
                    }
                } else {        /* LBR */
                    numwraps += (rphs - begphs) * search->proflen;       /* All E's */
                    if (lphs <= endphs && endphs <= rphs) {     /* LER */
                        numwraps += (endphs - lphs) * search->proflen;
                        if (begphs <= endphs)
                            numwraps -= 1.0;
                    }
                }
                /* printf("%.2f ", numwraps); */
                calctotalphs += numwraps;
                if (numwraps > 0)
                    search->rawfolds[ii * search->proflen + jj] *=
                        (totalphs / numwraps);
            }
            /* printf("\n"); */
            calctotalphs /= search->proflen;
            if (fabs(totalphs - calctotalphs) > 0.00001)
                printf
                    ("\nThere seems to be a problem in the \"exposure\" calculation\n"
                     "  in prepfold (npart = %ld):  totalphs = %.6f but calctotalphs = %0.6f\n",
                     ii, totalphs, calctotalphs);
            cts = 0.0;
            for (jj = ii * search->proflen; jj < (ii + 1) * search->proflen; jj++)
                cts += search->rawfolds[jj];
            search->stats[ii].numdata = ceil((T / fc->npart) / search->dt);
            search->stats[ii].numprof = search->proflen;
            search->stats[ii].prof_avg = search->stats[ii].prof_var =
                cts / search->proflen;
            search->stats[ii].data_avg = search->stats[ii].data_var =
                search->stats[ii].prof_avg / search->stats[ii].numdata;
            /* Compute the Chi-Squared probability that there is a signal */
            /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */
            search->stats[ii].redchi = 0.0;
            for (jj = ii * search->proflen; jj < (ii + 1) * search->proflen; jj++) {
                dtmp = search->rawfolds[jj] - search->stats[ii].prof_avg;
                search->stats[ii].redchi += dtmp * dtmp;
            }
            search->stats[ii].redchi /= (search->stats[ii].prof_var *
                                        (search->proflen - 1));
        }
        printf("\r  Folded %d events.", numevents);
        fflush(NULL);
    }

    fc->foldf = foldf;
    fc->foldfd = foldfd;
    fc->foldfdd = foldfdd;
    fc->orig_foldf = orig_foldf;
    fc->polyco_index = polyco_index;
    fc->parttimes = parttimes;
}


/* Phase 10: generate the topocentric<->barycentric corrections and convert */
/* the barycentric fold params to topocentric ones (called by fold_timeseries). */
void compute_bary_corrections(Cmdline *cmd, foldcand *fc, double T,
                              char *rastring, char *decstring, char *obs,
                              char *ephem, int *numbarypts_out,
                              double **barytimes_out, double **topotimes_out)
{
    prepfoldinfo *search = &fc->search;
    int binary = fc->binary;
    long long numbinpoints = fc->numbinpoints;
    double *tp = fc->tp;
    double f = fc->f, fd = fc->fd, fdd = fc->fdd;
    double foldf = fc->foldf, foldfd = fc->foldfd, foldfdd = fc->foldfdd;
    int numbarypts = *numbarypts_out;
    double *barytimes = *barytimes_out, *topotimes = *topotimes_out;
    int numdelays = fc->numdelays;
    int info, arrayoffset = 0;
    double dtmp;
    long ii;

    {                           /* Correct our fold parameters if we are barycentering */
        double *voverc;

        /* The number of topo to bary points to generate with TEMPO */

        numbarypts = T / TDT + 10;
        barytimes = gen_dvect(numbarypts);
        topotimes = gen_dvect(numbarypts);
        voverc = gen_dvect(numbarypts);

        /* topocentric times in days from data start */

        for (ii = 0; ii < numbarypts; ii++)
            topotimes[ii] = search->tepoch + (double) ii *TDT / SECPERDAY;

        /* Call TEMPO for the barycentering */

        printf("\nGenerating barycentric corrections...\n");
        barycenter(topotimes, barytimes, voverc, numbarypts,
                   rastring, decstring, obs, ephem);

        /* Determine the avg v/c of the Earth's motion during the obs */

        for (ii = 0; ii < numbarypts - 1; ii++)
            search->avgvoverc += voverc[ii];
        search->avgvoverc /= (numbarypts - 1.0);
        vect_free(voverc);
        printf("The average topocentric velocity is %.6g (units of c).\n\n",
               search->avgvoverc);
        printf("Barycentric folding frequency    (hz)  =  %-.12g\n", f);
        printf("Barycentric folding f-dot      (hz/s)  =  %-.8g\n", fd);
        printf("Barycentric folding f-dotdot (hz/s^2)  =  %-.8g\n", fdd);

        /* Convert the barycentric folding parameters into topocentric */

        if ((info = bary2topo(topotimes, barytimes, numbarypts,
                              f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
            printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
        printf("Topocentric folding frequency    (hz)  =  %-.12g\n", foldf);
        if (foldfd != 0.0)
            printf("Topocentric folding f-dot      (hz/s)  =  %-.8g\n", foldfd);
        if (foldfdd != 0.0)
            printf("Topocentric folding f-dotdot (hz/s^2)  =  %-.8g\n", foldfdd);
        printf("\n");

        /* Modify the binary delay times so they refer to */
        /* topocentric reference times.                   */

        if (binary) {
            for (ii = 0; ii < numbinpoints; ii++) {
                arrayoffset++;      /* Beware nasty NR zero-offset kludges! */
                dtmp = search->bepoch + tp[ii] / SECPERDAY;
                hunt(barytimes, numbarypts, dtmp, &arrayoffset);
                arrayoffset--;
                tp[ii] = LININTERP(dtmp, barytimes[arrayoffset],
                                   barytimes[arrayoffset + 1],
                                   topotimes[arrayoffset],
                                   topotimes[arrayoffset + 1]);
            }
            numdelays = numbinpoints;
            dtmp = tp[0];
            for (ii = 0; ii < numdelays; ii++)
                tp[ii] = (tp[ii] - dtmp) * SECPERDAY;
        }
    }

    fc->foldf = foldf;
    fc->foldfd = foldfd;
    fc->foldfdd = foldfdd;
    *numbarypts_out = numbarypts;
    *barytimes_out = barytimes;
    *topotimes_out = topotimes;
    fc->numdelays = numdelays;
}


/* Phase 11: for raw radio data, determine the per-channel dispersion delays */
/* in bins and print the DM-trial info (called by fold_timeseries).          */
void compute_dispersion_delays(Cmdline *cmd, infodata *idata, foldcand *fc,
                               int numchan, int insubs)
{
    prepfoldinfo *search = &fc->search;
    double foldf = fc->foldf;
    double *obsf = fc->obsf;
    int *idispdts = fc->idispdts;
    long ii;

    /* If this is 'raw' radio data, determine the dispersion delays */

    if (!strcmp(idata->band, "Radio")) {

        /* The observation frequencies */

        obsf = gen_dvect(numchan);
        obsf[0] = idata->freq;
        search->numchan = numchan;
        search->lofreq = idata->freq;
        search->bestdm = idata->dm;
        search->chan_wid = idata->chan_wid;
        for (ii = 1; ii < numchan; ii++)
            obsf[ii] = obsf[0] + ii * idata->chan_wid;
        if (RAWDATA || insubs) {
            for (ii = 0; ii < numchan; ii++)
                obsf[ii] = doppler(obsf[ii], search->avgvoverc);
        }
        {
            double *dispdts;
            dispdts = subband_search_delays(numchan, fc->nsub, fc->dm,
                                            idata->freq, idata->chan_wid,
                                            search->avgvoverc);
            idispdts = gen_ivect(numchan);
            /* Convert the delays in seconds to delays in bins */
            for (ii = 0; ii < numchan; ii++)
                idispdts[ii] = (int) (dispdts[ii] / search->dt + 0.5);
            vect_free(dispdts);
        }

        if (fc->nsub > 1 && (RAWDATA || insubs)) {
            int numdmtrials;
            double dphase, lodm, hidm, ddm;

            dphase = 1 / (foldf * search->proflen);
            ddm = dm_from_delay(dphase * cmd->dmstep, obsf[0]);
            numdmtrials = 2 * cmd->ndmfact * search->proflen + 1;
            lodm = fc->dm - (numdmtrials - 1) / 2 * ddm;
            if (lodm < 0.0)
                lodm = 0.0;
            hidm = lodm + numdmtrials * ddm;
            printf("Will search %d DMs from %.3f to %.3f (ddm = %.4f)\n",
                   cmd->nodmsearchP ? 1 : numdmtrials, lodm, hidm, ddm);
        }
    }

    fc->obsf = obsf;
    fc->idispdts = idispdts;
}


/* Phase 12: fold a time series.  Seeks to the start record, computes the    */
/* barycentric/dispersion corrections, then loops over npart sub-integrations */
/* x reads_per_part reads, folding each frequency subband.                   */
void fold_timeseries(Cmdline *cmd, struct spectra_info *s, infodata *idata,
                     foldcand *fc, mask *obsmask, double T, double N, long long lorec,
                     int ptsperrec, double recdt, long worklen, long reads_per_part,
                     int numchan, double startTday, int useshorts, int insubs,
                     int *maskchans, char *obs, char *ephem, char *rastring,
                     char *decstring, float **data_out, int *numbarypts_out,
                     double **barytimes_out, double **topotimes_out)
{
    prepfoldinfo *search = &fc->search;
    double f = fc->f, fd = fc->fd, fdd = fc->fdd;
    double *Ep = fc->Ep, *tp = fc->tp;
    int flags = fc->flags;
    double polyco_phase0 = fc->polyco_phase0;
    float *data = *data_out;
    double foldf = fc->foldf, foldfd = fc->foldfd, foldfdd = fc->foldfdd;
    double orig_foldf = fc->orig_foldf;
    int numbarypts = *numbarypts_out;
    double *barytimes = *barytimes_out, *topotimes = *topotimes_out;
    double *obsf = fc->obsf;
    int *idispdts = fc->idispdts;
    double *parttimes = fc->parttimes;
    int numdelays = fc->numdelays;
    int polyco_index = fc->polyco_index;
    double *buffers, *phasesadded;
    double proftime, polyco_phase = 0.0, dtmp;
    long long totnumfolded = 0;
    long numread = 0;
    int padding = 0, nummasked = 0;
    long ii, jj, kk;

    buffers = gen_dvect(fc->nsub * search->proflen);
    phasesadded = gen_dvect(fc->nsub);
    for (ii = 0; ii < fc->nsub * search->proflen; ii++)
        buffers[ii] = 0.0;
    for (ii = 0; ii < fc->nsub; ii++)
        phasesadded[ii] = 0.0;

    /* Move to the correct starting record */

    data = gen_fvect(fc->nsub * worklen);
    if (RAWDATA) {
        printf("\rTrue starting fraction       =  %g\n",
               (double) (lorec * ptsperrec) / s->N);
        offset_to_spectra(lorec * ptsperrec, s);
    } else {
        if (useshorts) {
            int reclen = 1;

            if (insubs)
                reclen = SUBSBLOCKLEN;
            /* Use a loop to accommodate subband data */
            for (ii = 0; ii < s->num_files; ii++)
                chkfileseek(s->files[ii], lorec * reclen, sizeof(short),
                            SEEK_SET);
        } else {
            chkfileseek(s->files[0], lorec, sizeof(float), SEEK_SET);
        }
    }

    /* Data is already barycentered or
       the polycos will take care of the barycentering or
       we are folding topocentrically.  */
    if (!(RAWDATA || insubs) || cmd->polycofileP || cmd->topoP) {
        foldf = f;
        foldfd = fd;
        foldfdd = fdd;
        orig_foldf = foldf;
    } else {                /* Correct our fold parameters if we are barycentering */
        compute_bary_corrections(cmd, fc, T, rastring, decstring, obs, ephem,
                                 &numbarypts, &barytimes, &topotimes);
        foldf = fc->foldf;
        foldfd = fc->foldfd;
        foldfdd = fc->foldfdd;
        numdelays = fc->numdelays;
    }

    /* Record the f, fd, and fdd we used to do the raw folding */

    if (idata->bary)
        search->fold.pow = 1.0;
    search->fold.p1 = foldf;
    search->fold.p2 = foldfd;
    search->fold.p3 = foldfdd;

    /* If this is 'raw' radio data, determine the dispersion delays */

    fc->foldf = foldf;
    compute_dispersion_delays(cmd, idata, fc, numchan, insubs);
    obsf = fc->obsf;
    idispdts = fc->idispdts;

    /*
     *   Perform the actual folding of the data
     */

    printf("\nStarting work on '%s'...\n\n", search->filenm);
    proftime = worklen * search->dt;
    parttimes = gen_dvect(fc->npart);
    printf("  Folded %lld points of %.0f", totnumfolded, N);

    /* sub-integrations in time  */

    dtmp = (double) fc->npart;
    for (ii = 0; ii < fc->npart; ii++) {
        parttimes[ii] = ii * reads_per_part * proftime;

        /* reads per sub-integration */

        for (jj = 0; jj < reads_per_part; jj++) {
            double fold_time0;

            if (RAWDATA) {
                numread =
                    read_subbands(data, idispdts, fc->nsub, s, 1, &padding,
                                  maskchans, &nummasked, obsmask);
            } else if (insubs) {
                numread = read_PRESTO_subbands(s->files, s->num_files, data, recdt,
                                               maskchans, &nummasked, obsmask,
                                               s->padvals);
            } else {
                int mm;
                float runavg = 0.0;
                static float oldrunavg = 0.0;
                static int firsttime = 1;

                if (useshorts)
                    numread = read_shorts(s->files[0], data, worklen, numchan);
                else
                    numread = read_floats(s->files[0], data, worklen, numchan);
                if (cmd->runavgP) {
                    for (mm = 0; mm < numread; mm++)
                        runavg += data[mm];
                    runavg /= numread;
                    if (firsttime) {
                        firsttime = 0;
                    } else {
                        // Use a running average of the block averages to subtract...
                        runavg = 0.95 * oldrunavg + 0.05 * runavg;
                    }
                    oldrunavg = runavg;
                    for (mm = 0; mm < numread; mm++)
                        data[mm] -= runavg;
                }
            }

            if (cmd->polycofileP) { /* Update the period/phase */
                double mjdf, currentsec, currentday, offsetphase, orig_cmd_phs =
                    0.0;

                if (ii == 0 && jj == 0)
                    orig_cmd_phs = fc->phs;
                currentsec = parttimes[ii] + jj * proftime;
                currentday = currentsec / SECPERDAY;
                mjdf = idata->mjd_f + startTday + currentday;
                /* Calculate the pulse phase at the start of the current block */
                polyco_index =
                    phcalc(idata->mjd_i, mjdf, polyco_index, &polyco_phase,
                           &foldf);
                if (!cmd->absphaseP)
                    polyco_phase -= polyco_phase0;
                if (polyco_phase < 0.0)
                    polyco_phase += 1.0;
                /* Calculate the folding frequency at the middle of the current block */
                polyco_index =
                    phcalc(idata->mjd_i, mjdf + 0.5 * proftime / SECPERDAY,
                           polyco_index, &offsetphase, &foldf);
                fc->phs = orig_cmd_phs + polyco_phase;
                fold_time0 = 0.0;
            } else {
                fold_time0 = parttimes[ii] + jj * proftime;
            }

            /* Fold the frequency sub-bands */

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
            for (kk = 0; kk < fc->nsub; kk++) {
                /* This is a quick hack to see if it will remove power drifts */
                if (cmd->runavgP && (numread > 0)) {
                    int dataptr;
                    double avg, var;
                    avg_var(data + kk * worklen, numread, &avg, &var);
                    for (dataptr = 0; dataptr < worklen; dataptr++)
                        data[kk * worklen + dataptr] -= avg;
                }
                fold(data + kk * worklen, numread, search->dt,
                     fold_time0,
                     search->rawfolds + (ii * fc->nsub + kk) * search->proflen,
                     search->proflen, fc->phs, buffers + kk * search->proflen,
                     phasesadded + kk, foldf, foldfd, foldfdd, flags, Ep, tp,
                     numdelays, NULL, &(search->stats[ii * fc->nsub + kk]),
                     !cmd->samplesP);
            }
            totnumfolded += numread;
        }

        printf("\r  Folded %lld points of %.0f", totnumfolded, N);
        fflush(NULL);
    }
    vect_free(buffers);
    vect_free(phasesadded);

    *data_out = data;
    fc->foldf = foldf;
    fc->foldfd = foldfd;
    fc->foldfdd = foldfdd;
    fc->orig_foldf = orig_foldf;
    *numbarypts_out = numbarypts;
    *barytimes_out = barytimes;
    *topotimes_out = topotimes;
    fc->obsf = obsf;
    fc->idispdts = idispdts;
    fc->parttimes = parttimes;
    fc->numdelays = numdelays;
    fc->polyco_index = polyco_index;
}


/* Phase 13: search over DM / period / p-dot / p-dotdot for the best profile, */
/* recording the best stats, the p-pdot plane, and the best p, pd, pdd, DM.   */
void optimize_candidate(Cmdline *cmd, infodata *idata, foldcand *fc,
                        plotflags *pflags, double T)
{
    prepfoldinfo *search = &fc->search;
    foldstats *beststats = &fc->beststats;
    double foldf = fc->foldf, foldfd = fc->foldfd, foldfdd = fc->foldfdd;
    double *obsf = fc->obsf;
    double *parttimes = fc->parttimes;
    float *ppdot = fc->ppdot;
    double *bestprof;
    double dtmp;
    long ii;

    // Normalize the profiles if requested (this changes the stats in
    // the .ofd file, unlike if you use this option for show_pfd)
    if (cmd->normalizeP)
        normalize_stats(search->rawfolds, search->stats, \
            search->npart, search->nsub, search->proflen);

    bestprof = gen_dvect(search->proflen);
    {
        int numtrials, totpdelay;
        int good_idm = 0, good_ip = 0, good_ipd = 0, good_ipdd = 0;
        double dphase, po, pdo, pddo;
        double *currentprof, *fdots, *fdotdots = NULL;
        foldstats currentstats;

        search->ndmfact = cmd->ndmfact;
        search->npfact = cmd->npfact;
        search->pstep = cmd->pstep;
        search->pdstep = cmd->pdstep;
        search->dmstep = cmd->dmstep;

        /* The number of trials for the P-dot and P searches */

        numtrials = 2 * search->npfact * search->proflen + 1;

        /* If we really don't need to search for the pulsations, */
        /* Don't make us search (and display) a hugh p/pd plane  */

        if (cmd->nosearchP && search->proflen > 100)
            numtrials = 201;

        /* Initialize a bunch of variables */

        search->numperiods = numtrials;
        search->periods = gen_dvect(numtrials);
        search->pdots = gen_dvect(numtrials);
        fdots = gen_dvect(numtrials);
        fdotdots = gen_dvect(numtrials);
        if (cmd->searchfddP)
            cmd->searchpddP = 1;
        if (cmd->nopsearchP && cmd->nopdsearchP) {
            if (fc->nsub > 1 && cmd->nodmsearchP)
                cmd->nosearchP = 1;
            else if (fc->nsub == 1)
                cmd->nosearchP = 1;
            pflags->nosearch = cmd->nosearchP;
        }
        if (cmd->nosearchP) {
            cmd->nopsearchP = cmd->nopdsearchP = 1;
            if (fc->nsub > 1)
                cmd->nodmsearchP = 1;
        }
        search->numpdots = numtrials;
        currentprof = gen_dvect(search->proflen);
        initialize_foldstats(beststats);
        if (cmd->nopsearchP)
            good_ip = (numtrials - 1) / 2;
        if (cmd->nopdsearchP)
            good_ipd = (numtrials - 1) / 2;
        if (cmd->nosearchP && cmd->searchpddP)
            good_ipdd = (numtrials - 1) / 2;

        /* Convert the folding freqs and derivs to periods */

        switch_f_and_p(foldf, foldfd, foldfdd, &po, &pdo, &pddo);

        /* Our P and P-dot steps are the changes that cause the pulse      */
        /* to be delayed a number of bins between the first and last time. */

        dphase = po / search->proflen;
        for (ii = 0; ii < numtrials; ii++) {
            totpdelay = ii - (numtrials - 1) / 2;
            dtmp = (double) (totpdelay * search->pstep) / search->proflen;
            search->periods[ii] = 1.0 / (foldf + dtmp / T);
            dtmp = (double) (totpdelay * search->pdstep) / search->proflen;
            fdots[ii] = phasedelay2fdot(dtmp, T);
            search->pdots[ii] = switch_pfdot(foldf, foldfd + fdots[ii]);
            fdotdots[ii] = phasedelay2fdotdot(dtmp, T);
        }

        {                       /* Do the optimization */
            int numdmtrials = 1, numpdds = 1, lodmnum = 0;
            long long totnumtrials, currtrial = 0;
            int oldper = -1, newper = 0;
            int idm, ip, ipd, ipdd, bestidm = 0, bestip = 0, bestipd = 0, bestipdd =
                0;
            double lodm = 0.0, ddm = 0.0;
            double *delays, *pd_delays, *pdd_delays, *ddprofs = search->rawfolds;
            foldstats *ddstats = search->stats;

            if (cmd->searchpddP)
                numpdds = numtrials;
            delays = gen_dvect(fc->npart);
            pd_delays = gen_dvect(fc->npart);
            pdd_delays = gen_dvect(fc->npart);

            if (fc->nsub > 1) {        /* This is only for doing DM searches */
                ddprofs = gen_dvect(fc->npart * search->proflen);
                ddstats = (foldstats *) malloc(fc->npart * sizeof(foldstats));
                numdmtrials = 2 * search->ndmfact * search->proflen + 1;
                search->numdms = numdmtrials;
                search->dms = gen_dvect(numdmtrials);

                /* Our DM step is the change in DM that would cause the pulse   */
                /* to be delayed a number of phasebins at the lowest frequency. */

                ddm = dm_from_delay(dphase * search->dmstep, obsf[0]);

                /* Insure that we don't try a dm < 0.0 */

                lodm = fc->dm - (numdmtrials - 1) / 2 * ddm;
                if (cmd->nodmsearchP)
                    good_idm = (numdmtrials - 1) / 2;
                if (lodm < 0.0) {
                    lodm = 0.0;
                    /* Find the closest DM to the requested DM */
                    if (cmd->nodmsearchP) {
                        double mindmerr = 1000.0, dmerr, trialdm;
                        for (idm = 0; idm < numdmtrials; idm++) {       /* Loop over DMs */
                            trialdm = lodm + idm * ddm;
                            dmerr = fabs(trialdm - fc->dm);
                            if (dmerr < mindmerr) {
                                good_idm = idm;
                                mindmerr = dmerr;
                            }
                        }
                    }
                }
                for (idm = 0; idm < numdmtrials; idm++)
                    search->dms[idm] = lodm + idm * ddm;

                if (cmd->searchpddP)
                    printf
                        ("  Searching %d DMs, %d periods, %d p-dots, and %d p-dotdots...\n",
                         cmd->nodmsearchP ? 1 : search->numdms,
                         cmd->nopsearchP ? 1 : search->numperiods,
                         cmd->nopdsearchP ? 1 : search->numpdots,
                         cmd->searchpddP ? search->numpdots : 1);
                else
                    printf("  Searching %d DMs, %d periods, and %d p-dots...\n",
                           cmd->nodmsearchP ? 1 : search->numdms,
                           cmd->nopsearchP ? 1 : search->numperiods,
                           cmd->nopdsearchP ? 1 : search->numpdots);
            } else {            /* No DM searches are to be done */
                /* Allocate our p-pdot plane to speed up plotting */
                if (!cmd->searchpddP)
                    ppdot = gen_fvect(search->numpdots * search->numperiods);
                if (cmd->searchpddP)
                    printf
                        ("  Searching %d periods, %d p-dots, and %d p-dotdots...\n",
                         cmd->nopsearchP ? 1 : search->numperiods,
                         cmd->nopdsearchP ? 1 : search->numpdots,
                         cmd->searchpddP ? search->numpdots : 1);
                else
                    printf("  Searching %d periods, and %d p-dots...\n",
                           cmd->nopsearchP ? 1 : search->numperiods,
                           cmd->nopdsearchP ? 1 : search->numpdots);
            }

            /* Skip searching over all the DMs if we know that we're */
            /* only going to look at one of them.                    */
            if (cmd->nodmsearchP) {
                lodmnum = good_idm;
                numdmtrials = lodmnum + 1;
            }
            /* The total number of search trials */
            totnumtrials = (cmd->nodmsearchP ? 1l : numdmtrials) *
                (cmd->nopsearchP ? 1l : numtrials) *
                (cmd->nopdsearchP ? 1l : numtrials) *
                (cmd->searchpddP ? numpdds : 1);
            if (totnumtrials < 100000000l) {
                printf("     (%lld total trials)\n", totnumtrials);
            } else if (totnumtrials < 1000000000l) {
                printf("     Warning!:  This is %lld trials!  It will take a while!\n",
                       totnumtrials);
            } else {
                printf("     Warning!:  This is %lld trials!  This will take forever!\n"
                       "                You almost certainly want some type of -nosearch!!\n",
                       totnumtrials);
            }

            for (idm = lodmnum; idm < numdmtrials; idm++) {     /* Loop over DMs */
                if (fc->nsub > 1) {    /* This is only for doing DM searches */
                    if (!cmd->nodmsearchP)
                        good_idm = idm;
                    correct_subbands_for_DM(search->dms[idm], search, ddprofs,
                                            ddstats);
                }

                for (ipdd = 0; ipdd < numpdds; ipdd++) {        /* Loop over pdds */
                    if (!cmd->nosearchP)
                        good_ipdd = ipdd;
                    for (ii = 0; ii < fc->npart; ii++)
                        pdd_delays[ii] = cmd->searchpddP ?
                            fdotdot2phasedelay(fdotdots[ipdd],
                                               parttimes[ii]) * search->proflen : 0.0;

                    for (ipd = 0; ipd < numtrials; ipd++) {     /* Loop over the pds */
                        if (!cmd->nopdsearchP)
                            good_ipd = ipd;
                        else if (cmd->nopdsearchP && fc->nsub > 1
                                 && ipd != good_ipd)
                            continue;
                        for (ii = 0; ii < fc->npart; ii++)
                            pd_delays[ii] = pdd_delays[ii] +
                                fdot2phasedelay(fdots[ipd],
                                                parttimes[ii]) * search->proflen;

                        for (ip = 0; ip < numtrials; ip++) {    /* Loop over the ps */
                            if (!cmd->nopsearchP)
                                good_ip = ip;
                            else if (cmd->nopsearchP && fc->nsub > 1
                                     && ip != good_ip)
                                continue;
                            totpdelay = search->pstep * (ip - (numtrials - 1) / 2);
                            for (ii = 0; ii < fc->npart; ii++)
                                delays[ii] =
                                    pd_delays[ii] +
                                    (double) (ii * totpdelay) / fc->npart;

                            /* Combine the profiles usingthe above computed delays */
                            combine_profs(ddprofs, ddstats, fc->npart,
                                          search->proflen, delays, currentprof,
                                          &currentstats);

                            /* If this is a simple fold, create the chi-square p-pdot plane */
                            if (fc->nsub == 1 && !cmd->searchpddP)
                                ppdot[ipd * search->numpdots + ip] =
                                    currentstats.redchi;

                            /* If this is the best profile or if it is the specific */
                            /* profile that we are looking for, save it.            */
                            if (idm == good_idm && ipdd == good_ipdd
                                && ipd == good_ipd && ip == good_ip) {
                                if (cmd->nosearchP
                                    || (currentstats.redchi > beststats->redchi
                                        && !cmd->nosearchP)) {
                                    bestidm = idm;
                                    bestipdd = ipdd;
                                    bestipd = ipd;
                                    bestip = ip;
                                    *beststats = currentstats;
                                    memcpy(bestprof, currentprof,
                                           sizeof(double) * search->proflen);
                                }
                            }
                            currtrial += 1;
                            newper =
                                (int) (((double) currtrial) / totnumtrials * 100.0 +
                                       0.5);
                            if (newper > oldper) {
                                printf("\r  Amount Complete = %3d%%", newper);
                                fflush(stdout);
                                oldper = newper;
                            }
                        }
                    }
                }
            }

            /* Convert the indices of the folds into p, pd, pdd and DM values */
            if (fc->nsub > 1)
                search->bestdm = search->dms[bestidm];
            {
                double p1 = search->periods[bestip];
                double p2 = search->pdots[bestipd];
                double p3;
                if (cmd->searchpddP)
                    p3 = switch_pfdotdot(1.0 / p1, foldfd + fdots[bestipd],
                                         foldfdd + fdotdots[bestipdd]);
                else
                    p3 = switch_pfdotdot(1.0 / p1, foldfd + fdots[bestipd], foldfdd);
                if (idata->bary) {
                    search->bary.p1 = p1; search->bary.p2 = p2; search->bary.p3 = p3;
                } else {
                    search->topo.p1 = p1; search->topo.p2 = p2; search->topo.p3 = p3;
                }
            }
            vect_free(delays);
            vect_free(pd_delays);
            vect_free(pdd_delays);
            if (fc->nsub > 1) {
                vect_free(ddprofs);
                free(ddstats);
            }
        }
        vect_free(currentprof);
        vect_free(fdots);
        vect_free(fdotdots);
    }

    fc->ppdot = ppdot;
    fc->bestprof = bestprof;
}


/* Phase 14: convert the best parameters between bary/topo, print them, write */
/* the prepfoldinfo structure, and make the plot.                             */
void write_results_and_plot(Cmdline *cmd, infodata *idata, foldcand *fc,
                            plotflags *pflags, double N, int insubs,
                            double *barytimes, double *topotimes, int numbarypts)
{
    prepfoldinfo *search = &fc->search;
    foldstats *beststats = &fc->beststats;
    double *bestprof = fc->bestprof;
    float *ppdot = fc->ppdot;
    char *outfilenm = fc->outfilenm;
    /* Write and plot the results and cleanup */

    {
        double perr, pderr, pdderr;
        char out[100];
        double f, fd, fdd, foldf, foldfd, foldfdd;
        int info;

        printf("Maximum reduced chi-squared found  =  %-.5f\n", beststats->redchi);
        if (fc->nsub > 1)
            printf("Best DM     (pc cm^-3)  =  %-.4f\n", search->bestdm);

        /* Convert best params from/to barycentric to/from topocentric */

        if (!(RAWDATA || insubs) || cmd->polycofileP || cmd->topoP) {

            /* Data was barycentered */

            if (idata->bary && !cmd->polycofileP) {

                /* Calculate the errors in our new pulsation quantities */

                fold_errors(bestprof, search->proflen, search->dt, N,
                            beststats->data_var, search->bary.p1, search->bary.p2,
                            search->bary.p3, &perr, &pderr, &pdderr);

                nice_output_2(out, search->bary.p1, perr, 0);
                printf("Best period        (s)  =  %s\n", out);
                nice_output_2(out, search->bary.p2, pderr, 0);
                printf("Best p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->bary.p3, pdderr, 0);
                printf("Best p-dotdot  (s/s^2)  =  %s\n", out);

                /* Data was topocentric */

            } else {

                /* Calculate the errors in our new pulsation quantities */

                fold_errors(bestprof, search->proflen, search->dt, N,
                            beststats->data_var, search->topo.p1, search->topo.p2,
                            search->topo.p3, &perr, &pderr, &pdderr);

                nice_output_2(out, search->topo.p1, perr, 0);
                printf("Best period        (s)  =  %s\n", out);
                nice_output_2(out, search->topo.p2, pderr, 0);
                printf("Best p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->topo.p3, pdderr, 0);
                printf("Best p-dotdot  (s/s^2)  =  %s\n", out);
                if (idata->bary)
                    search->topo.p1 = 0.0;
            }

        } else {

            if (idata->bary && !cmd->polycofileP) {

                /* Calculate the errors in our new pulsation quantities */

                fold_errors(bestprof, search->proflen, search->dt, N,
                            beststats->data_var, search->bary.p1, search->bary.p2,
                            search->bary.p3, &perr, &pderr, &pdderr);

                nice_output_2(out, search->bary.p1, perr, 0);
                printf("Best barycentric period        (s)  =  %s\n", out);
                nice_output_2(out, search->bary.p2, pderr, 0);
                printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->bary.p3, pdderr, 0);
                printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);

                /* Convert the barycentric folding parameters into topocentric */

                switch_f_and_p(search->bary.p1, search->bary.p2, search->bary.p3,
                               &f, &fd, &fdd);
                if ((info = bary2topo(topotimes, barytimes, numbarypts,
                                      f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
                    printf("\nError in bary2topo().  Argument %d was bad.\n\n",
                           -info);
                switch_f_and_p(foldf, foldfd, foldfdd, &search->topo.p1,
                               &search->topo.p2, &search->topo.p3);

                nice_output_2(out, search->topo.p1, perr, 0);
                printf("Best topocentric period        (s)  =  %s\n", out);
                nice_output_2(out, search->topo.p2, pderr, 0);
                printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->topo.p3, pdderr, 0);
                printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

            } else {

                /* Calculate the errors in our new pulsation quantities */

                fold_errors(bestprof, search->proflen, search->dt, N,
                            beststats->data_var, search->topo.p1, search->topo.p2,
                            search->topo.p3, &perr, &pderr, &pdderr);

                nice_output_2(out, search->topo.p1, perr, 0);
                printf("Best topocentric period        (s)  =  %s\n", out);
                nice_output_2(out, search->topo.p2, pderr, 0);
                printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->topo.p3, pdderr, 0);
                printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

                /* Convert the barycentric folding parameters into topocentric */

                switch_f_and_p(search->topo.p1, search->topo.p2, search->topo.p3,
                               &f, &fd, &fdd);
                if ((info = bary2topo(barytimes, topotimes, numbarypts,
                                      f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
                    printf("\nError in bary2topo().  Argument %d was bad.\n\n",
                           -info);
                switch_f_and_p(foldf, foldfd, foldfdd, &search->bary.p1,
                               &search->bary.p2, &search->bary.p3);

                nice_output_2(out, search->bary.p1, perr, 0);
                printf("Best barycentric period        (s)  =  %s\n", out);
                nice_output_2(out, search->bary.p2, pderr, 0);
                printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
                nice_output_2(out, search->bary.p3, pdderr, 0);
                printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);
            }
        }
    }

    printf("\nMaking plots.\n\n");

    /*
     *   Write the raw prepfoldinfo structure
     */

    write_prepfoldinfo(search, outfilenm);

    /*
     *   Plot our results
     */

    prepfold_plot(search, pflags, !cmd->noxwinP, ppdot);
}


/* Allocate and zero a foldcand, initializing its embedded prepfoldinfo. */
void init_foldcand(foldcand *fc)
{
    memset(fc, 0, sizeof(foldcand));
    init_prepfoldinfo(&fc->search);
}


/* Free the per-candidate (foldcand) allocations.  Mirrors the per-candidate  */
/* portion of the original cleanup_fold, guarding each free exactly as before. */
void free_foldcand(foldcand *fc, Cmdline *cmd, infodata *idata)
{
    if (fc->nsub == 1 && !cmd->searchpddP)
        vect_free(fc->ppdot);
    delete_prepfoldinfo(&fc->search);
    free(fc->outfilenm);
    free(fc->plotfilenm);
    vect_free(fc->parttimes);
    vect_free(fc->bestprof);
    if (fc->binary) {
        vect_free(fc->Ep);
        vect_free(fc->tp);
    }
    if (!strcmp(idata->band, "Radio")) {
        vect_free(fc->obsf);
        vect_free(fc->idispdts);
    }
}


/* Phase 15: free the shared, per-observation allocations.  The per-candidate */
/* allocations are released separately by free_foldcand().                    */
void cleanup_fold(Cmdline *cmd, mask *obsmask, int insubs, float *data,
                  char *rootnm, double *barytimes, double *topotimes)
{
    /* Free our memory  */

    vect_free(data);
    if (!cmd->outfileP)
        free(rootnm);
    if (cmd->maskfileP)
        free_mask(*obsmask);
    if (RAWDATA || insubs) {
        vect_free(barytimes);
        vect_free(topotimes);
    }
}






