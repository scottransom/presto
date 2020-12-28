#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "prepsubband_cmd.h"
#include "mask.h"
#include "backend_common.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_LONG(x) (long) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

static void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite);
static void write_subs(FILE * outfiles[], int numfiles, short **subsdata,
                       int startpoint, int numtowrite);
static void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite);
static int read_PRESTO_subbands(FILE * infiles[], int numfiles,
                                float *subbanddata, double timeperblk,
                                int *maskchans, int *nummasked, mask * obsmask,
                                float clip_sigma, float *padvals);
static int get_data(float **outdata, int blocksperread,
                    struct spectra_info *s,
                    mask * obsmask, int *idispdts, int **offsets,
                    int *padding, short **subsdata);
static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins, int downsamp);
static void print_percent_complete(int current, int number);

/* From CLIG */
static int insubs = 0;
static Cmdline *cmd;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    /* Any variable that begins with 't' means topocentric */
    /* Any variable that begins with 'b' means barycentric */
    FILE **outfiles;
    float **outdata = NULL;
    short **subsdata = NULL;
    double dtmp, *dms = NULL, avgdm = 0.0, maxdm, dsdt = 0;
    double tlotoa = 0.0, blotoa = 0.0, BW_ddelay = 0.0;
    double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
    double *btoa = NULL, *ttoa = NULL, avgvoverc = 0.0;
    char obs[3], ephem[10], rastring[50], decstring[50];
    long totnumtowrite, totwrote = 0, padwrote = 0, datawrote = 0;
    int **offsets;
    int ii, jj, numadded = 0, numremoved = 0, padding = 0;
    int numbarypts = 0, blocksperread = 0, worklen = 0;
    int numread = 0, numtowrite = 0;
    int padtowrite = 0, statnum = 0, good_padvals = 0;
    int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL;
    int *idispdt;
    char *datafilenm;
    int dmprecision = 2;
    struct spectra_info s;
    infodata idata;
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
    spectra_info_set_defaults(&s);
    dmprecision = cmd->dmprec;
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
    if (cmd->noclipP) {
        cmd->clip = 0.0;
        s.clip_sigma = 0.0;
    }
    if (cmd->ifsP) {
        // 0 = default or summed, 1-4 are possible also
        s.use_poln = cmd->ifs + 1;
    }
    if (!cmd->numoutP)
        cmd->numout = LONG_MAX;

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

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("          Pulsar Subband De-dispersion Routine\n");
    printf("                 by Scott M. Ransom\n\n");

    if (RAWDATA) {
        if (cmd->filterbankP)
            s.datatype = SIGPROCFB;
        else if (cmd->psrfitsP)
            s.datatype = PSRFITS;
    } else {                    // Attempt to auto-identify the data
        identify_psrdatatype(&s, 1);
        if (s.datatype == SIGPROCFB)
            cmd->filterbankP = 1;
        else if (s.datatype == PSRFITS)
            cmd->psrfitsP = 1;
        else if (s.datatype == SUBBAND)
            insubs = 1;
        else {
            printf
                ("Error:  Unable to identify input data files.  Please specify type.\n\n");
            exit(1);
        }
    }

    if (!RAWDATA)
        s.files = (FILE **) malloc(sizeof(FILE *) * s.num_files);
    if (RAWDATA || insubs) {
        char description[40];
        psrdatatype_description(description, s.datatype);
        if (s.num_files > 1)
            printf("Reading %s data from %d files:\n", description, s.num_files);
        else
            printf("Reading %s data from 1 file:\n", description);
        for (ii = 0; ii < s.num_files; ii++) {
            printf("  '%s'\n", cmd->argv[ii]);
            if (insubs)
                s.files[ii] = chkfopen(s.filenames[ii], "rb");
        }
        printf("\n");
        if (RAWDATA) {
            read_rawdata_files(&s);
            // Make sure that the requested number of subbands divides into the
            // the raw number of channels.
            if (s.num_channels % cmd->nsub) {
                printf("Error:  The number of subbands (-nsub %d) must divide into the\n"
                       "        number of channels (%d)\n\n",
                       cmd->nsub, s.num_channels);
                exit(1);
            }
            if (cmd->ignorechanstrP) {
                s.ignorechans = get_ignorechans(cmd->ignorechanstr, 0, s.num_channels-1,
                                                &s.num_ignorechans, &s.ignorechans_str);
                if (s.ignorechans_str==NULL) {
                    s.ignorechans_str = (char *)malloc(strlen(cmd->ignorechanstr)+1);
                    strcpy(s.ignorechans_str, cmd->ignorechanstr);
                }
            }
            print_spectra_info_summary(&s);
            spectra_info_to_inf(&s, &idata);
        } else {                // insubs
            cmd->nsub = s.num_files;
            s.N = chkfilelen(s.files[0], sizeof(short));
            s.padvals = gen_fvect(s.num_files);
            for (ii = 0; ii < s.num_files; ii++)
                s.padvals[ii] = 0.0;
            s.start_MJD = (long double *) malloc(sizeof(long double));
            s.start_spec = (long long *) malloc(sizeof(long long));
            s.num_spec = (long long *) malloc(sizeof(long long));
            s.num_pad = (long long *) malloc(sizeof(long long));
            s.start_spec[0] = 0L;
            s.num_spec[0] = s.N;
            s.num_pad[0] = 0L;
        }
        /* Read an input mask if wanted */
        if (cmd->maskfileP) {
            read_mask(cmd->maskfile, &obsmask);
            printf("Read mask information from '%s'\n\n", cmd->maskfile);
            good_padvals = determine_padvals(cmd->maskfile, &obsmask, s.padvals);
        } else {
            obsmask.numchan = obsmask.numint = 0;
        }
    }

    if (insubs) {
        char *root, *suffix;
        if (split_root_suffix(s.filenames[0], &root, &suffix) == 0) {
            printf("Error:  The input filename (%s) must have a suffix!\n\n",
                   s.filenames[0]);
            exit(1);
        }
        if (strncmp(suffix, "sub", 3) == 0) {
            char *tmpname;
            tmpname = calloc(strlen(root) + 10, 1);
            sprintf(tmpname, "%s.sub", root);
            readinf(&idata, tmpname);
            free(tmpname);
            s.num_channels = idata.num_chan;
            s.start_MJD[0] = idata.mjd_i + idata.mjd_f;
            s.dt = idata.dt;
            s.T = s.N * s.dt;
            s.lo_freq = idata.freq;
            s.df = idata.chan_wid;
            s.hi_freq = s.lo_freq + (s.num_channels - 1.0) * s.df;
            s.BW = s.num_channels * s.df;
            s.fctr = s.lo_freq - 0.5 * s.df + 0.5 * s.BW;
            s.spectra_per_subint = SUBSBLOCKLEN;
            print_spectra_info_summary(&s);
        } else {
            printf("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                   cmd->argv[0]);
            exit(1);
        }
        free(root);
        free(suffix);
    }

    /* Determine the output file names and open them */

    datafilenm = (char *) calloc(strlen(cmd->outfile) + 20, 1);
    if (!cmd->subP) {
        printf("Writing output data to:\n");
        outfiles = (FILE **) malloc(cmd->numdms * sizeof(FILE *));
        dms = gen_dvect(cmd->numdms);
        for (ii = 0; ii < cmd->numdms; ii++) {
            dms[ii] = cmd->lodm + ii * cmd->dmstep;
            avgdm += dms[ii];
            sprintf(datafilenm, "%s_DM%.*f.dat", cmd->outfile, dmprecision, dms[ii]);
            outfiles[ii] = chkfopen(datafilenm, "wb");
            printf("   '%s'\n", datafilenm);
        }
        avgdm /= cmd->numdms;
        maxdm = dms[cmd->numdms - 1];
    } else {
        char format_str[30];
        int num_places;

        if (!cmd->nobaryP) {
            printf("\nWarning:  You cannot (currently) barycenter subbands.\n"
                   "          Setting the '-nobary' flag automatically.\n");
            cmd->nobaryP = 1;
        }
        printf("Writing subbands to:\n");
        cmd->numdms = 1;
        dms = gen_dvect(cmd->numdms);
        dms[0] = cmd->subdm;
        cmd->lodm = cmd->subdm;
        avgdm = cmd->subdm;
        maxdm = cmd->subdm;
        outfiles = (FILE **) malloc(cmd->nsub * sizeof(FILE *));
        num_places = (int) ceil(log10(cmd->nsub));
        sprintf(format_str, "%%s_DM%%.*f.sub%%0%dd", num_places);
        for (ii = 0; ii < cmd->nsub; ii++) {
            sprintf(datafilenm, format_str, cmd->outfile, dmprecision, avgdm, ii);
            outfiles[ii] = chkfopen(datafilenm, "wb");
            printf("   '%s'\n", datafilenm);
        }
    }

    /* Set a few other key values */
    if (insubs)
        avgdm = idata.dm;
    if (RAWDATA)
        idata.dm = avgdm;
    dsdt = cmd->downsamp * idata.dt;
    BW_ddelay = delay_from_dm(maxdm, idata.freq) -
        delay_from_dm(maxdm, idata.freq + (idata.num_chan - 1) * idata.chan_wid);
    blocksperread = ((int) (BW_ddelay / idata.dt) / s.spectra_per_subint + 1);
    worklen = s.spectra_per_subint * blocksperread;
    /* The number of topo to bary time points to generate with TEMPO */
    numbarypts = (int) (s.T * 1.1 / TDT + 5.5) + 1;

    // Identify the TEMPO observatory code
    {
        char *outscope = (char *) calloc(40, sizeof(char));
        telescope_to_tempocode(idata.telescope, outscope, obs);
        free(outscope);
    }

    /* If we are offsetting into the file, change inf file start time */
    if (cmd->start > 0.0 || cmd->offset > 0) {
        if (cmd->start > 0.0) /* Offset in units of worklen */
            cmd->offset = (long) (cmd->start *
                                  idata.N / worklen) * worklen;
        add_to_inf_epoch(&idata, cmd->offset * idata.dt);
        printf("Offsetting into the input files by %ld spectra (%.6g sec)\n",
               cmd->offset, cmd->offset * idata.dt);
        if (RAWDATA)
            offset_to_spectra(cmd->offset, &s);
        else { // subbands
            for (ii = 0; ii < s.num_files; ii++)
                chkfileseek(s.files[ii], cmd->offset, sizeof(short), SEEK_SET);
            if (cmd->maskfileP)
                printf("WARNING!:  masking does not work with old-style subbands and -start or -offset!\n");
        }
    }

    if (cmd->nsub > s.num_channels) {
        printf
            ("Warning:  The number of requested subbands (%d) is larger than the number of channels (%d).\n",
             cmd->nsub, s.num_channels);
        printf("          Re-setting the number of subbands to %d.\n\n",
               s.num_channels);
        cmd->nsub = s.num_channels;
    }

    if (s.spectra_per_subint % cmd->downsamp) {
        printf
            ("Error:  The downsample factor (%d) must be a factor of the\n",
             cmd->downsamp);
        printf("        blocklength (%d).  Exiting.\n\n", s.spectra_per_subint);
        exit(1);
    }

    tlotoa = idata.mjd_i + idata.mjd_f; /* Topocentric epoch */

    /* Set the output length to a good number if it wasn't requested */
    if (!cmd->numoutP && !cmd->subP) {
        cmd->numoutP = 1;
        cmd->numout = choose_good_N((long long)(idata.N/cmd->downsamp));
        printf("Setting a 'good' output length of %ld samples\n", cmd->numout);
    }
    if (cmd->subP && (cmd->numout > idata.N/cmd->downsamp))
        cmd->numout = (long long)(idata.N/cmd->downsamp); // Don't pad subbands
    totnumtowrite = cmd->numout;

    if (cmd->nobaryP) {         /* Main loop if we are not barycentering... */
        double *dispdt;

        /* Dispersion delays (in bins).  The high freq gets no delay   */
        /* All other delays are positive fractions of bin length (dt)  */

        dispdt = subband_search_delays(s.num_channels, cmd->nsub, avgdm,
                                       idata.freq, idata.chan_wid, 0.0);
        idispdt = gen_ivect(s.num_channels);
        for (ii = 0; ii < s.num_channels; ii++)
            idispdt[ii] = NEAREST_LONG(dispdt[ii] / idata.dt);
        vect_free(dispdt);

        /* The subband dispersion delays (see note above) */

        offsets = gen_imatrix(cmd->numdms, cmd->nsub);
        for (ii = 0; ii < cmd->numdms; ii++) {
            double *subdispdt;

            subdispdt = subband_delays(s.num_channels, cmd->nsub, dms[ii],
                                       idata.freq, idata.chan_wid, 0.0);
            dtmp = subdispdt[cmd->nsub - 1];
            for (jj = 0; jj < cmd->nsub; jj++)
                offsets[ii][jj] = NEAREST_LONG((subdispdt[jj] - dtmp) / dsdt);
            vect_free(subdispdt);
        }

        /* Allocate our data array and start getting data */

        printf("\nDe-dispersing using:\n");
        printf("       Subbands = %d\n", cmd->nsub);
        printf("     Average DM = %.7g\n", avgdm);
        if (cmd->downsamp > 1) {
            printf("     Downsample = %d\n", cmd->downsamp);
            printf("  New sample dt = %.10g\n", dsdt);
        }
        printf("\n");

        if (cmd->subP)
            subsdata = gen_smatrix(cmd->nsub, worklen / cmd->downsamp);
        else
            outdata = gen_fmatrix(cmd->numdms, worklen / cmd->downsamp);
        numread = get_data(outdata, blocksperread, &s,
                           &obsmask, idispdt, offsets, &padding, subsdata);

        while (numread == worklen) {

            numread /= cmd->downsamp;
            print_percent_complete(totwrote, totnumtowrite);

            /* Write the latest chunk of data, but don't   */
            /* write more than cmd->numout points.         */

            numtowrite = numread;
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            if (cmd->subP)
                write_subs(outfiles, cmd->nsub, subsdata, 0, numtowrite);
            else
                write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
            totwrote += numtowrite;

            /* Update the statistics */

            if (!padding && !cmd->subP) {
                for (ii = 0; ii < numtowrite; ii++)
                    update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg,
                                 &var);
                statnum += numtowrite;
            }

            /* Stop if we have written out all the data we need to */

            if (totwrote == cmd->numout)
                break;

            numread = get_data(outdata, blocksperread, &s,
                               &obsmask, idispdt, offsets, &padding, subsdata);
        }
        datawrote = totwrote;

    } else {                    /* Main loop if we are barycentering... */
        double maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;
        double *dispdt;

        /* What ephemeris will we use?  (Default is DE405) */
        strcpy(ephem, "DE405");

        /* Define the RA and DEC of the observation */

        ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
        ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

        /* Allocate some arrays */

        btoa = gen_dvect(numbarypts);
        ttoa = gen_dvect(numbarypts);
        voverc = gen_dvect(numbarypts);
        for (ii = 0; ii < numbarypts; ii++)
            ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

        /* Call TEMPO for the barycentering */

        printf("\nGenerating barycentric corrections...\n");
        barycenter(ttoa, btoa, voverc, numbarypts, rastring, decstring, obs, ephem);
        for (ii = 0; ii < numbarypts; ii++) {
            if (voverc[ii] > maxvoverc)
                maxvoverc = voverc[ii];
            if (voverc[ii] < minvoverc)
                minvoverc = voverc[ii];
            avgvoverc += voverc[ii];
        }
        avgvoverc /= numbarypts;
        vect_free(voverc);
        blotoa = btoa[0];

        printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
        printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
        printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
        printf("De-dispersing and barycentering using:\n");
        printf("       Subbands = %d\n", cmd->nsub);
        printf("     Average DM = %.7g\n", avgdm);
        if (cmd->downsamp > 1) {
            printf("     Downsample = %d\n", cmd->downsamp);
            printf("  New sample dt = %.10g\n", dsdt);
        }
        printf("\n");

        /* Dispersion delays (in bins).  The high freq gets no delay   */
        /* All other delays are positive fractions of bin length (dt)  */

        dispdt = subband_search_delays(s.num_channels, cmd->nsub, avgdm,
                                       idata.freq, idata.chan_wid, avgvoverc);
        idispdt = gen_ivect(s.num_channels);
        for (ii = 0; ii < s.num_channels; ii++)
            idispdt[ii] = NEAREST_LONG(dispdt[ii] / idata.dt);
        vect_free(dispdt);

        /* The subband dispersion delays (see note above) */

        offsets = gen_imatrix(cmd->numdms, cmd->nsub);
        for (ii = 0; ii < cmd->numdms; ii++) {
            double *subdispdt;

            subdispdt = subband_delays(s.num_channels, cmd->nsub, dms[ii],
                                       idata.freq, idata.chan_wid, avgvoverc);
            dtmp = subdispdt[cmd->nsub - 1];
            for (jj = 0; jj < cmd->nsub; jj++)
                offsets[ii][jj] = NEAREST_LONG((subdispdt[jj] - dtmp) / dsdt);
            vect_free(subdispdt);
        }

        /* Convert the bary TOAs to differences from the topo TOAs in */
        /* units of bin length (dt) rounded to the nearest integer.   */

        dtmp = (btoa[0] - ttoa[0]);
        for (ii = 0; ii < numbarypts; ii++)
            btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

        {                       /* Find the points where we need to add or remove bins */

            int oldbin = 0, currentbin;
            double lobin, hibin, calcpt;

            numdiffbins = labs(NEAREST_LONG(btoa[numbarypts - 1])) + 1;
            diffbins = gen_ivect(numdiffbins);
            diffbinptr = diffbins;
            for (ii = 1; ii < numbarypts; ii++) {
                currentbin = NEAREST_LONG(btoa[ii]);
                if (currentbin != oldbin) {
                    if (currentbin > 0) {
                        calcpt = oldbin + 0.5;
                        lobin = (ii - 1) * TDT / dsdt;
                        hibin = ii * TDT / dsdt;
                    } else {
                        calcpt = oldbin - 0.5;
                        lobin = -((ii - 1) * TDT / dsdt);
                        hibin = -(ii * TDT / dsdt);
                    }
                    while (fabs(calcpt) < fabs(btoa[ii])) {
                        /* Negative bin number means remove that bin */
                        /* Positive bin number means add a bin there */
                        *diffbinptr = NEAREST_LONG(LININTERP(calcpt, btoa[ii - 1],
                                                             btoa[ii], lobin,
                                                             hibin));
                        diffbinptr++;
                        calcpt = (currentbin > 0) ? calcpt + 1.0 : calcpt - 1.0;
                    }
                    oldbin = currentbin;
                }
            }
            *diffbinptr = cmd->numout;  /* Used as a marker */
        }
        diffbinptr = diffbins;

        /* Now perform the barycentering */

        if (cmd->subP)
            subsdata = gen_smatrix(cmd->nsub, worklen / cmd->downsamp);
        else
            outdata = gen_fmatrix(cmd->numdms, worklen / cmd->downsamp);
        numread = get_data(outdata, blocksperread, &s,
                           &obsmask, idispdt, offsets, &padding, subsdata);

        while (numread == worklen) {    /* Loop to read and write the data */
            int numwritten = 0;
            double block_avg, block_var;

            numread /= cmd->downsamp;
            /* Determine the approximate local average */
            avg_var(outdata[0], numread, &block_avg, &block_var);
            print_percent_complete(totwrote, totnumtowrite);

            /* Simply write the data if we don't have to add or */
            /* remove any bins from this batch.                 */
            /* OR write the amount of data up to cmd->numout or */
            /* the next bin that will be added or removed.      */

            numtowrite = abs(*diffbinptr) - datawrote;
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            if (numtowrite > numread)
                numtowrite = numread;
            if (cmd->subP)
                write_subs(outfiles, cmd->nsub, subsdata, 0, numtowrite);
            else
                write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
            datawrote += numtowrite;
            totwrote += numtowrite;
            numwritten += numtowrite;

            /* Update the statistics */

            if (!padding && !cmd->subP) {
                for (ii = 0; ii < numtowrite; ii++)
                    update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg,
                                 &var);
                statnum += numtowrite;
            }

            if ((datawrote == abs(*diffbinptr)) && (numwritten != numread) && (totwrote < cmd->numout)) {       /* Add/remove a bin */
                int skip, nextdiffbin;

                skip = numtowrite;

                do {            /* Write the rest of the data after adding/removing a bin  */

                    if (*diffbinptr > 0) {
                        /* Add a bin */
                        write_padding(outfiles, cmd->numdms, block_avg, 1);
                        numadded++;
                        totwrote++;
                    } else {
                        /* Remove a bin */
                        numremoved++;
                        datawrote++;
                        numwritten++;
                        skip++;
                    }
                    diffbinptr++;

                    /* Write the part after the diffbin */

                    numtowrite = numread - numwritten;
                    if ((totwrote + numtowrite) > cmd->numout)
                        numtowrite = cmd->numout - totwrote;
                    nextdiffbin = abs(*diffbinptr) - datawrote;
                    if (numtowrite > nextdiffbin)
                        numtowrite = nextdiffbin;
                    if (cmd->subP)
                        write_subs(outfiles, cmd->nsub, subsdata, skip, numtowrite);
                    else
                        write_data(outfiles, cmd->numdms, outdata, skip, numtowrite);
                    numwritten += numtowrite;
                    datawrote += numtowrite;
                    totwrote += numtowrite;

                    /* Update the statistics and counters */

                    if (!padding && !cmd->subP) {
                        for (ii = 0; ii < numtowrite; ii++)
                            update_stats(statnum + ii, outdata[0][skip + ii],
                                         &min, &max, &avg, &var);
                        statnum += numtowrite;
                    }
                    skip += numtowrite;

                    /* Stop if we have written out all the data we need to */

                    if (totwrote == cmd->numout)
                        break;
                } while (numwritten < numread);
            }
            /* Stop if we have written out all the data we need to */

            if (totwrote == cmd->numout)
                break;

            numread = get_data(outdata, blocksperread, &s,
                               &obsmask, idispdt, offsets, &padding, subsdata);
        }
    }

    /* Calculate the amount of padding we need (don't pad subbands) */

    if (!cmd->subP && (cmd->numout > totwrote))
        padwrote = padtowrite = cmd->numout - totwrote;

    /* Write the new info file for the output data */

    idata.dt = dsdt;
    update_infodata(&idata, totwrote, padtowrite, diffbins,
                    numdiffbins, cmd->downsamp);
    for (ii = 0; ii < cmd->numdms; ii++) {
        idata.dm = dms[ii];
        if (!cmd->nobaryP) {
            double baryepoch, barydispdt, baryhifreq;

            baryhifreq = idata.freq + (s.num_channels - 1) * idata.chan_wid;
            barydispdt = delay_from_dm(dms[ii], doppler(baryhifreq, avgvoverc));
            baryepoch = blotoa - (barydispdt / SECPERDAY);
            idata.bary = 1;
            idata.mjd_i = (int) floor(baryepoch);
            idata.mjd_f = baryepoch - idata.mjd_i;
        }
        if (cmd->subP)
            sprintf(idata.name, "%s_DM%.*f.sub", cmd->outfile, dmprecision, dms[ii]);
        else
            sprintf(idata.name, "%s_DM%.*f", cmd->outfile, dmprecision, dms[ii]);
        writeinf(&idata);
    }

    /* Set the padded points equal to the average data point */

    if (idata.numonoff >= 1) {
        int index, startpad, endpad;

        for (ii = 0; ii < cmd->numdms; ii++) {
            fclose(outfiles[ii]);
            sprintf(datafilenm, "%s_DM%.*f.dat", cmd->outfile, dmprecision, dms[ii]);
            outfiles[ii] = chkfopen(datafilenm, "rb+");
        }
        for (ii = 0; ii < idata.numonoff; ii++) {
            index = 2 * ii;
            startpad = idata.onoff[index + 1];
            if (ii == idata.numonoff - 1)
                endpad = idata.N - 1;
            else
                endpad = idata.onoff[index + 2];
            for (jj = 0; jj < cmd->numdms; jj++)
                chkfseek(outfiles[jj], (startpad + 1) * sizeof(float), SEEK_SET);
            padtowrite = endpad - startpad;
            write_padding(outfiles, cmd->numdms, avg, padtowrite);
        }
    }

    /* Print simple stats and results */

    if (!cmd->subP) {
        var /= (datawrote - 1);
        print_percent_complete(1, 1);
        printf("\n\nDone.\n\nSimple statistics of the output data:\n");
        printf("             Data points written:  %ld\n", totwrote);
        if (padwrote)
            printf("          Padding points written:  %ld\n", padwrote);
        if (!cmd->nobaryP) {
            if (numadded)
                printf("    Bins added for barycentering:  %d\n", numadded);
            if (numremoved)
                printf("  Bins removed for barycentering:  %d\n", numremoved);
        }
        printf("           Maximum value of data:  %.2f\n", max);
        printf("           Minimum value of data:  %.2f\n", min);
        printf("              Data average value:  %.2f\n", avg);
        printf("         Data standard deviation:  %.2f\n", sqrt(var));
        printf("\n");
    } else {
        printf("\n\nDone.\n");
        printf("             Data points written:  %ld\n", totwrote);
        if (padwrote)
            printf("          Padding points written:  %ld\n", padwrote);
        if (!cmd->nobaryP) {
            if (numadded)
                printf("    Bins added for barycentering:  %d\n", numadded);
            if (numremoved)
                printf("  Bins removed for barycentering:  %d\n", numremoved);
        }
        printf("\n");
    }

    /* Close the files and cleanup */

    if (cmd->maskfileP) {
        free_mask(obsmask);
    }
    //  Close all the raw files and free their vectors
    close_rawfiles(&s);
    for (ii = 0; ii < cmd->numdms; ii++)
        fclose(outfiles[ii]);
    if (cmd->subP) {
        vect_free(subsdata[0]);
        vect_free(subsdata);
    } else {
        vect_free(outdata[0]);
        vect_free(outdata);
    }
    free(outfiles);
    vect_free(dms);
    vect_free(idispdt);
    vect_free(offsets[0]);
    vect_free(offsets);
    free(datafilenm);
    if (!cmd->nobaryP) {
        vect_free(btoa);
        vect_free(ttoa);
        vect_free(diffbins);
    }
    return (0);
}

static void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite)
{
    int ii;

    for (ii = 0; ii < numfiles; ii++)
        chkfwrite(outdata[ii] + startpoint, sizeof(float), numtowrite, outfiles[ii]);
}


static void write_subs(FILE * outfiles[], int numfiles, short **subsdata,
                       int startpoint, int numtowrite)
{
    int ii;

    for (ii = 0; ii < numfiles; ii++)
        chkfwrite(subsdata[ii] + startpoint, sizeof(short), numtowrite,
                  outfiles[ii]);
}


static void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite)
{
    int ii;

    if (numtowrite <= 0) {
        return;
    } else if (numtowrite == 1) {
        for (ii = 0; ii < numfiles; ii++)
            chkfwrite(&value, sizeof(float), 1, outfiles[ii]);
    } else {
        int maxatonce = 8192, veclen, jj;
        float *buffer;
        veclen = (numtowrite > maxatonce) ? maxatonce : numtowrite;
        buffer = gen_fvect(veclen);
        for (ii = 0; ii < veclen; ii++)
            buffer[ii] = value;
        if (veclen == numtowrite) {
            for (ii = 0; ii < numfiles; ii++)
                chkfwrite(buffer, sizeof(float), veclen, outfiles[ii]);
        } else {
            for (ii = 0; ii < numtowrite / veclen; ii++) {
                for (jj = 0; jj < numfiles; jj++)
                    chkfwrite(buffer, sizeof(float), veclen, outfiles[jj]);
            }
            for (jj = 0; jj < numfiles; jj++)
                chkfwrite(buffer, sizeof(float), numtowrite % veclen, outfiles[jj]);
        }
        vect_free(buffer);
    }
}


static int read_PRESTO_subbands(FILE * infiles[], int numfiles,
                                float *subbanddata, double timeperblk,
                                int *maskchans, int *nummasked, mask * obsmask,
                                float clip_sigma, float *padvals)
/* Read short int subband data written by prepsubband */
{
    int ii, jj, index, numread = 0, mask = 0, offset;
    short subsdata[SUBSBLOCKLEN];
    double starttime, run_avg;
    float subband_sum;
    static int currentblock = 0;

    if (obsmask->numchan)
        mask = 1;

    /* Read the data */
    for (ii = 0; ii < numfiles; ii++) {
        numread = chkfread(subsdata, sizeof(short), SUBSBLOCKLEN, infiles[ii]);
        run_avg = 0.0;
        if (cmd->runavgP == 1) {
            for (jj = 0; jj < numread; jj++)
                run_avg += (float) subsdata[jj];
            run_avg /= numread;
        }
        for (jj = 0, index = ii; jj < numread; jj++, index += numfiles)
            subbanddata[index] = (float) subsdata[jj] - run_avg;
        for (jj = numread; jj < SUBSBLOCKLEN; jj++, index += numfiles)
            subbanddata[index] = 0.0;
    }

    if (mask) {
        starttime = currentblock * timeperblk;
        *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
    }

    /* Clip nasty RFI if requested and we're not masking all the channels */
    if ((clip_sigma > 0.0) && !(mask && (*nummasked == -1))) {
        clip_times(subbanddata, SUBSBLOCKLEN, numfiles, clip_sigma, padvals);
    }

    /* Mask it if required */
    if (mask && numread) {
        if (*nummasked == -1) { /* If all channels are masked */
            for (ii = 0; ii < SUBSBLOCKLEN; ii++)
                memcpy(subbanddata + ii * numfiles, padvals,
                       sizeof(float) * numfiles);
        } else if (*nummasked > 0) {    /* Only some of the channels are masked */
            int channum;
            for (ii = 0; ii < SUBSBLOCKLEN; ii++) {
                offset = ii * numfiles;
                for (jj = 0; jj < *nummasked; jj++) {
                    channum = maskchans[jj];
                    subbanddata[offset + channum] = padvals[channum];
                }
            }
        }
    }

    /* Zero-DM removal if required */
    if (cmd->zerodmP == 1) {
        for (ii = 0; ii < SUBSBLOCKLEN; ii++) {
            offset = ii * numfiles;
            subband_sum = 0.0;
            for (jj = offset; jj < offset + numfiles; jj++) {
                subband_sum += subbanddata[jj];
            }
            subband_sum /= (float) numfiles;
            /* Remove the channel average */
            for (jj = offset; jj < offset + numfiles; jj++) {
                subbanddata[jj] -= subband_sum;
            }
        }
    }

    currentblock += 1;
    return numread;
}



static int get_data(float **outdata, int blocksperread,
                    struct spectra_info *s,
                    mask * obsmask, int *idispdts, int **offsets,
                    int *padding, short **subsdata)
{
    static int firsttime = 1, *maskchans = NULL, blocksize;
    static int worklen, dsworklen;
    static float *tempzz, *data1, *data2, *dsdata1 = NULL, *dsdata2 = NULL;
    static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
    static double blockdt;
    int totnumread = 0, numread = 0, ii, jj, tmppad = 0, nummasked = 0;

    if (firsttime) {
        if (cmd->maskfileP)
            maskchans = gen_ivect(s->num_channels);
        worklen = s->spectra_per_subint * blocksperread;
        dsworklen = worklen / cmd->downsamp;
        // Make sure that our working blocks are long enough...
        for (ii = 0; ii < cmd->numdms; ii++) {
            for (jj = 0; jj < cmd->nsub; jj++) {
                if (offsets[ii][jj] > dsworklen)
                    printf
                        ("WARNING!:  (offsets[%d][%d] = %d) > (dsworklen = %d)\n",
                         ii, jj, offsets[ii][jj], dsworklen);
            }
        }

        blocksize = s->spectra_per_subint * cmd->nsub;
        blockdt = s->spectra_per_subint * s->dt;
        data1 = gen_fvect(cmd->nsub * worklen);
        data2 = gen_fvect(cmd->nsub * worklen);
        currentdata = data1;
        lastdata = data2;
        if (cmd->downsamp > 1) {
            dsdata1 = gen_fvect(cmd->nsub * dsworklen);
            dsdata2 = gen_fvect(cmd->nsub * dsworklen);
            currentdsdata = dsdata1;
            lastdsdata = dsdata2;
        } else {
            currentdsdata = data1;
            lastdsdata = data2;
        }
    }
    while (1) {
        if (RAWDATA || insubs) {
            for (ii = 0; ii < blocksperread; ii++) {
                if (RAWDATA)
                    numread = read_subbands(currentdata + ii * blocksize, idispdts,
                                            cmd->nsub, s, 0, &tmppad,
                                            maskchans, &nummasked, obsmask);
                else if (insubs)
                    numread = read_PRESTO_subbands(s->files, s->num_files,
                                                   currentdata + ii * blocksize,
                                                   blockdt, maskchans, &nummasked,
                                                   obsmask, cmd->clip, s->padvals);
                if (!firsttime)
                    totnumread += numread;
                if (numread != s->spectra_per_subint) {
                    for (jj = ii * blocksize; jj < (ii + 1) * blocksize; jj++)
                        currentdata[jj] = 0.0;
                }
                if (tmppad)
                    *padding = 1;
            }
        }
        /* Downsample the subband data if needed */
        if (cmd->downsamp > 1) {
            int kk, index;
            float ftmp;
            for (ii = 0; ii < dsworklen; ii++) {
                const int dsoffset = ii * cmd->nsub;
                const int offset = dsoffset * cmd->downsamp;
                for (jj = 0; jj < cmd->nsub; jj++) {
                    const int dsindex = dsoffset + jj;
                    index = offset + jj;
                    currentdsdata[dsindex] = ftmp = 0.0;
                    for (kk = 0; kk < cmd->downsamp; kk++) {
                        ftmp += currentdata[index];
                        index += cmd->nsub;
                    }
                    currentdsdata[dsindex] += ftmp / cmd->downsamp;
                }
            }
        }
        if (firsttime) {
            SWAP(currentdata, lastdata);
            SWAP(currentdsdata, lastdsdata);
            firsttime = 0;
        } else
            break;
    }
    if (!cmd->subP) {
        for (ii = 0; ii < cmd->numdms; ii++)
            float_dedisp(currentdsdata, lastdsdata, dsworklen,
                         cmd->nsub, offsets[ii], 0.0, outdata[ii]);
    } else {
        /* Input format is sub1[0], sub2[0], sub3[0], ..., sub1[1], sub2[1], sub3[1], ... */
        float infloat;
        for (ii = 0; ii < cmd->nsub; ii++) {
            for (jj = 0; jj < dsworklen; jj++) {
                infloat = lastdsdata[ii + (cmd->nsub * jj)];
                subsdata[ii][jj] = (short) (infloat + 0.5);
                //if ((float) subsdata[ii][jj] != infloat)
                //   printf
                //       ("Warning:  We are incorrectly converting subband data! float = %f  short = %d\n",
                //         infloat, subsdata[ii][jj]);
            }
        }
    }
    SWAP(currentdata, lastdata);
    SWAP(currentdsdata, lastdsdata);
    if (totnumread != worklen) {
        if (cmd->maskfileP)
            vect_free(maskchans);
        vect_free(data1);
        vect_free(data2);
        if (cmd->downsamp > 1) {
            vect_free(dsdata1);
            vect_free(dsdata2);
        }
    }
    return totnumread;
}


static void print_percent_complete(int current, int number)
{
    static int newper = 0, oldper = -1;

    newper = (int) (current / (float) (number) * 100.0);
    if (newper < 0)
        newper = 0;
    if (newper > 100)
        newper = 100;
    if (newper > oldper) {
        printf("\rAmount complete = %3d%%", newper);
        fflush(stdout);
        oldper = newper;
    }
}


static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins, int downsamp)
/* Update our infodata for barycentering and padding */
{
    int ii, jj, index;

    idata->N = datawrote + padwrote;
    if (idata->numonoff == 0) {
        if (padwrote) {
            idata->numonoff = 2;
            idata->onoff[0] = 0.0;
            idata->onoff[1] = datawrote - 1;
            idata->onoff[2] = idata->N - 1;
            idata->onoff[3] = idata->N - 1;
        }
        return;
    } else {
        for (ii = 0; ii < idata->numonoff; ii++) {
            idata->onoff[ii * 2] /= downsamp;
            idata->onoff[ii * 2 + 1] /= downsamp;
        }
    }

    /* Determine the barycentric onoff bins (approximate) */

    if (numbarybins) {
        int numadded = 0, numremoved = 0;

        ii = 1;                 /* onoff index    */
        jj = 0;                 /* barybins index */
        while (ii < idata->numonoff * 2) {
            while (abs(barybins[jj]) <= idata->onoff[ii] && jj < numbarybins) {
                if (barybins[jj] < 0)
                    numremoved++;
                else
                    numadded++;
                jj++;
            }
            idata->onoff[ii] += numadded - numremoved;
            ii++;
        }
    }

    /* Now cut off the extra onoff bins */

    for (ii = 1, index = 1; ii <= idata->numonoff; ii++, index += 2) {
        if (idata->onoff[index - 1] > idata->N - 1) {
            idata->onoff[index - 1] = idata->N - 1;
            idata->onoff[index] = idata->N - 1;
            break;
        }
        if (idata->onoff[index] > datawrote - 1) {
            idata->onoff[index] = datawrote - 1;
            idata->numonoff = ii;
            if (padwrote) {
                idata->numonoff++;
                idata->onoff[index + 1] = idata->N - 1;
                idata->onoff[index + 2] = idata->N - 1;
            }
            break;
        }
    }
}
