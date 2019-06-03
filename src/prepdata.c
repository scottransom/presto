#include <limits.h>
#include <ctype.h>
#include "prepdata_cmd.h"
#include "presto.h"
#include "mask.h"
#include "backend_common.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_LONG(x) (long) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

/* Some function definitions */
static int read_floats(FILE * file, float *data, int numpts, int numchan);
static int read_shorts(FILE * file, float *data, int numpts, int numchan);
static int downsample(float outdata[], int numread, int downsampfact);
static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins);

/* The main program */

int main(int argc, char *argv[])
{
    /* Any variable that begins with 't' means topocentric */
    /* Any variable that begins with 'b' means barycentric */
    FILE *outfile;
    float *outdata = NULL;
    double tdf = 0.0, dtmp = 0.0, barydispdt = 0.0, dsdt = 0.0;
    double *dispdt, *tobsf = NULL, tlotoa = 0.0, blotoa = 0.0;
    double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
    char obs[3], ephem[10], *datafilenm, *outinfonm;
    char rastring[50], decstring[50];
    int numchan = 1, newper = 0, oldper = 0, nummasked = 0, useshorts = 0;
    int numadded = 0, numremoved = 0, padding = 0, *maskchans = NULL, offset = 0;
    long slen, ii, numbarypts = 0, worklen = 65536;
    long numread = 0, numtowrite = 0, totwrote = 0, datawrote = 0;
    long padwrote = 0, padtowrite = 0, statnum = 0;
    int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL, good_padvals = 0;
    int *idispdt;
    struct spectra_info s;
    infodata idata;
    Cmdline *cmd;
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
    printf("           Pulsar Data Preparation Routine\n");
    printf("    Type conversion, de-dispersion, barycentering.\n");
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
        else if (s.datatype == SDAT)
            useshorts = 1;
        else if (s.datatype != DAT) {
            printf
                ("Error:  Unable to identify input data files.  Please specify type.\n\n");
            exit(1);
        }
    }

    if (!RAWDATA) {
        char *root, *suffix;
        /* Split the filename into a rootname and a suffix */
        if (split_root_suffix(s.filenames[0], &root, &suffix) == 0) {
            printf("\nThe input filename (%s) must have a suffix!\n\n",
                   s.filenames[0]);
            exit(1);
        }
        printf("Reading input data from '%s'.\n", s.filenames[0]);
        printf("Reading information from '%s.inf'.\n\n", root);
        /* Read the info file if available */
        readinf(&idata, root);
        free(root);
        free(suffix);
        s.files = (FILE **) malloc(sizeof(FILE *));
        s.files[0] = chkfopen(s.filenames[0], "rb");
    } else {
        char description[40];
        psrdatatype_description(description, s.datatype);
        if (s.num_files > 1)
            printf("Reading %s data from %d files:\n", description, s.num_files);
        else
            printf("Reading %s data from 1 file:\n", description);
        for (ii = 0; ii < s.num_files; ii++) {
            printf("  '%s'\n", cmd->argv[ii]);
        }
        printf("\n");
    }

    /* Determine the other file names and open the output data file */
    slen = strlen(cmd->outfile) + 8;
    datafilenm = (char *) calloc(slen, 1);
    sprintf(datafilenm, "%s.dat", cmd->outfile);
    outfile = chkfopen(datafilenm, "wb");
    sprintf(idata.name, "%s", cmd->outfile);
    outinfonm = (char *) calloc(slen, 1);
    sprintf(outinfonm, "%s.inf", cmd->outfile);

    if (RAWDATA) {
        read_rawdata_files(&s);
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
        /* Finish setting up stuff common to all raw formats */
        idata.dm = cmd->dm;
        worklen = s.spectra_per_subint;

        /* If we are offsetting into the file, change inf file start time */
        if (cmd->start > 0.0 || cmd->offset > 0) {
            if (cmd->start > 0.0) /* Offset in units of worklen */
                cmd->offset = (long) (cmd->start *
                                      idata.N / worklen) * worklen;
            add_to_inf_epoch(&idata, cmd->offset * idata.dt);
            offset_to_spectra(cmd->offset, &s);
            printf("Offsetting into the input files by %ld spectra (%.6g sec)\n",
                   cmd->offset, cmd->offset * idata.dt);
        }
        if (cmd->maskfileP)
            maskchans = gen_ivect(idata.num_chan);

        /* Compare the size of the data to the size of output we request */
        if (cmd->numoutP) {
            dtmp = idata.N;
            idata.N = cmd->numout;
            writeinf(&idata);
            idata.N = dtmp;
        } else {
        /* Set the output length to a good number if it wasn't requested */
            cmd->numoutP = 1;
            cmd->numout = choose_good_N((long long)(idata.N/cmd->downsamp));
            writeinf(&idata);
            printf("Setting a 'good' output length of %ld samples\n", cmd->numout);
        }

        /* The number of topo to bary time points to generate with TEMPO */
        numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;

        // Identify the TEMPO observatory code
        {
            char *outscope = (char *) calloc(40, sizeof(char));
            telescope_to_tempocode(idata.telescope, outscope, obs);
            free(outscope);
        }
    }

    /* Read an input mask if wanted */
    if (cmd->maskfileP) {
        read_mask(cmd->maskfile, &obsmask);
        printf("Read mask information from '%s'\n\n", cmd->maskfile);
        good_padvals = determine_padvals(cmd->maskfile, &obsmask, s.padvals);
    } else {
        obsmask.numchan = obsmask.numint = 0;
    }

    /* Determine our initialization data if we do _not_ have Parkes, */
    /* Green Bank BCPM, or Arecibo WAPP data sets.                   */
    if (!RAWDATA) {

        /* If we will be barycentering... */
        if (!cmd->nobaryP) {
            /* The number of topo to bary time points to generate with TEMPO */
            numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;
            // Identify the TEMPO observatory code
            {
                char *outscope = (char *) calloc(40, sizeof(char));
                telescope_to_tempocode(idata.telescope, outscope, obs);
                free(outscope);
            }
        }

        /* The number of data points to work with at a time */
        if (worklen > idata.N)
            worklen = idata.N;
        worklen = (long) (worklen / 1024) * 1024;

        /* If we are offsetting into the file, change inf file start time */
        if (cmd->start > 0.0 || cmd->offset > 0) {
            if (cmd->start > 0.0) /* Offset in units of worklen */
                cmd->offset = (long) (cmd->start *
                                      idata.N / worklen) * worklen;
            add_to_inf_epoch(&idata, cmd->offset * idata.dt);
            printf("Offsetting into the input files by %ld samples (%.6g sec)\n",
                   cmd->offset, cmd->offset * idata.dt);
            if (useshorts) {
                chkfileseek(s.files[0], cmd->offset, sizeof(short), SEEK_SET);
            } else {
                chkfileseek(s.files[0], cmd->offset, sizeof(float), SEEK_SET);
            }
        }

        /* Set the output length to a good number if it wasn't requested */
        if (!cmd->numoutP) {
            cmd->numoutP = 1;
            cmd->numout = choose_good_N((long long)(idata.N/cmd->downsamp));
            printf("Setting a 'good' output length of %ld samples\n", cmd->numout);
        }
    }

    /* Check if we are downsampling */
    dsdt = idata.dt * cmd->downsamp;
    if (cmd->downsamp > 1) {
        printf("Downsampling by a factor of %d\n", cmd->downsamp);
        printf("New sample dt = %.10g\n\n", dsdt);
        if (worklen % cmd->downsamp) {
            printf("Error:  The downsample factor (%d) must be a factor of the\n",
                   cmd->downsamp);
            printf("        worklength (%ld).  Exiting.\n\n", worklen);
            exit(1);
        }
    }
    printf("Writing output data to '%s'.\n", datafilenm);
    printf("Writing information to '%s'.\n\n", outinfonm);

    /* The topocentric epoch of the start of the data */
    tlotoa = (double) idata.mjd_i + idata.mjd_f;

    if (!strcmp(idata.band, "Radio") && RAWDATA) {

        /* The topocentric spacing between channels */
        tdf = idata.chan_wid;
        numchan = idata.num_chan;

        /* The topocentric observation frequencies */
        tobsf = gen_dvect(numchan);
        tobsf[0] = idata.freq;
        for (ii = 0; ii < numchan; ii++)
            tobsf[ii] = tobsf[0] + ii * tdf;

        /* The dispersion delays (in time bins) */
        dispdt = gen_dvect(numchan);    // full float bins
        idispdt = gen_ivect(numchan);   // nearest integer bins

        if (cmd->nobaryP) {

            /* Determine our dispersion time delays for each channel */
            for (ii = 0; ii < numchan; ii++)
                dispdt[ii] = delay_from_dm(cmd->dm, tobsf[ii]);

            /* The highest frequency channel gets no delay                 */
            /* All other delays are positive fractions of bin length (dt)  */
            dtmp = dispdt[numchan - 1];
            for (ii = 0; ii < numchan; ii++) {
                dispdt[ii] = (dispdt[ii] - dtmp) / idata.dt;
                idispdt[ii] = (int) (dispdt[ii] + 0.5);
            }
            worklen *= ((int) (fabs(dispdt[0])) / worklen) + 1;
        }

    } else {                    /* For unknown radio raw data (Why is this here?) */
        tobsf = gen_dvect(numchan);
        dispdt = gen_dvect(numchan);
        idispdt = gen_ivect(numchan);
        dispdt[0] = 0.0;
        idispdt[0] = 0;
        if (!strcmp(idata.band, "Radio")) {
            tobsf[0] = idata.freq + (idata.num_chan - 1) * idata.chan_wid;
            cmd->dm = idata.dm;
        } else {
            tobsf[0] = 0.0;
            cmd->dm = 0.0;
        }
    }

    if (cmd->nobaryP) {         /* Main loop if we are not barycentering... */

        /* Allocate our data array */
        outdata = gen_fvect(worklen);

        printf("Massaging the data ...\n\n");
        printf("Amount Complete = 0%%");

        do {
            if (RAWDATA)
                numread = read_psrdata(outdata, worklen, &s, idispdt, &padding,
                                       maskchans, &nummasked, &obsmask);
            else if (useshorts)
                numread = read_shorts(s.files[0], outdata, worklen, numchan);
            else
                numread = read_floats(s.files[0], outdata, worklen, numchan);
            if (numread == 0)
                break;

            /* Downsample if requested */
            if (cmd->downsamp > 1)
                numread = downsample(outdata, numread, cmd->downsamp);

            /* Print percent complete */
            newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
            if (newper > oldper) {
                printf("\rAmount Complete = %3d%%", newper);
                fflush(stdout);
                oldper = newper;
            }

            /* Write the latest chunk of data, but don't   */
            /* write more than cmd->numout points.         */
            numtowrite = numread;
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            chkfwrite(outdata, sizeof(float), numtowrite, outfile);
            totwrote += numtowrite;

            /* Update the statistics */
            if (!padding) {
                for (ii = 0; ii < numtowrite; ii++)
                    update_stats(statnum + ii, outdata[ii], &min, &max, &avg, &var);
                statnum += numtowrite;
            }

            /* Stop if we have written out all the data we need to */
            if (totwrote == cmd->numout)
                break;

        } while (numread);

        datawrote = totwrote;

    } else {                    /* Main loop if we are barycentering... */

        double avgvoverc = 0.0, maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;
        double *bobsf = NULL, *btoa = NULL, *ttoa = NULL;

        /* What ephemeris will we use?  (Default is DE405) */
        strcpy(ephem, "DE405");

        /* Define the RA and DEC of the observation */
        ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
        ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

        /* Allocate some arrays */
        bobsf = gen_dvect(numchan);
        btoa = gen_dvect(numbarypts);
        ttoa = gen_dvect(numbarypts);
        voverc = gen_dvect(numbarypts);
        for (ii = 0; ii < numbarypts; ii++)
            ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

        /* Call TEMPO for the barycentering */
        printf("Generating barycentric corrections...\n");
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
        printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

        /* Determine the initial dispersion time delays for each channel */
        for (ii = 0; ii < numchan; ii++) {
            bobsf[ii] = doppler(tobsf[ii], avgvoverc);
            dispdt[ii] = delay_from_dm(cmd->dm, bobsf[ii]);
        }

        /* The highest frequency channel gets no delay                   */
        /* All other delays are positive fractions of bin length (dt)    */
        barydispdt = dispdt[numchan - 1];
        for (ii = 0; ii < numchan; ii++) {
            dispdt[ii] = (dispdt[ii] - barydispdt) / idata.dt;
            idispdt[ii] = (int) (dispdt[ii] + 0.5);
        }
        if (RAWDATA)
            worklen *= ((int) (dispdt[0]) / worklen) + 1;

        /* If the data is de-dispersed radio data... */
        if (!strcmp(idata.band, "Radio")) {
            printf("The DM of %.2f at the barycentric observing freq of %.3f MHz\n",
                   idata.dm, bobsf[numchan - 1]);
            printf("   causes a delay of %f seconds compared to infinite freq.\n",
                   barydispdt);
            printf("   This delay is removed from the barycented times.\n\n");
        }
        printf("Topocentric epoch (at data start) is:\n");
        printf("   %17.11f\n\n", tlotoa);
        printf("Barycentric epoch (infinite obs freq at data start) is:\n");
        printf("   %17.11f\n\n", blotoa - (barydispdt / SECPERDAY));

        /* Convert the bary TOAs to differences from the topo TOAs in  */
        /* units of bin length (dsdt) rounded to the nearest integer.  */
        dtmp = (btoa[0] - ttoa[0]);
        for (ii = 0; ii < numbarypts; ii++)
            btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

        {                       /* Find the points where we need to add or remove bins */

            int oldbin = 0, currentbin;
            double lobin, hibin, calcpt;

            numdiffbins = abs(NEAREST_LONG(btoa[numbarypts - 1])) + 1;
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

        printf("Massaging the data...\n\n");
        printf("Amount Complete = 0%%");

        /* Allocate our data array */

        outdata = gen_fvect(worklen);

        do {                    /* Loop to read and write the data */
            int numwritten = 0;
            double block_avg, block_var;

            if (RAWDATA)
                numread = read_psrdata(outdata, worklen, &s, idispdt, &padding,
                                       maskchans, &nummasked, &obsmask);
            else if (useshorts)
                numread = read_shorts(s.files[0], outdata, worklen, numchan);
            else
                numread = read_floats(s.files[0], outdata, worklen, numchan);
            if (numread == 0)
                break;

            /* Downsample if requested */
            if (cmd->downsamp > 1)
                numread = downsample(outdata, numread, cmd->downsamp);

            /* Determine the approximate local average */
            avg_var(outdata, numread, &block_avg, &block_var);

            /* Print percent complete */

            newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
            if (newper > oldper) {
                printf("\rAmount Complete = %3d%%", newper);
                fflush(stdout);
                oldper = newper;
            }

            /* Simply write the data if we don't have to add or */
            /* remove any bins from this batch.                 */
            /* OR write the amount of data up to cmd->numout or */
            /* the next bin that will be added or removed.      */

            numtowrite = abs(*diffbinptr) - datawrote;
            /* FIXME: numtowrite+totwrote can wrap! */
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            if (numtowrite > numread)
                numtowrite = numread;
            chkfwrite(outdata, sizeof(float), numtowrite, outfile);
            datawrote += numtowrite;
            totwrote += numtowrite;
            numwritten += numtowrite;

            /* Update the statistics */

            if (!padding) {
                for (ii = 0; ii < numtowrite; ii++)
                    update_stats(statnum + ii, outdata[ii], &min, &max, &avg, &var);
                statnum += numtowrite;
            }

            if ((datawrote == abs(*diffbinptr)) && (numwritten != numread) && (totwrote < cmd->numout)) {       /* Add/remove a bin */
                float favg;
                int skip, nextdiffbin;

                skip = numtowrite;

                do {            /* Write the rest of the data after adding/removing a bin  */

                    if (*diffbinptr > 0) {

                        /* Add a bin */

                        favg = (float) block_avg;
                        chkfwrite(&favg, sizeof(float), 1, outfile);
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
                    chkfwrite(outdata + skip, sizeof(float), numtowrite, outfile);
                    numwritten += numtowrite;
                    datawrote += numtowrite;
                    totwrote += numtowrite;

                    /* Update the statistics and counters */

                    if (!padding) {
                        for (ii = 0; ii < numtowrite; ii++)
                            update_stats(statnum + ii, outdata[skip + ii],
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

        } while (numread);

        /* Free the arrays used in barycentering */

        vect_free(bobsf);
        vect_free(btoa);
        vect_free(ttoa);
    }

    /* Calculate what the amount of padding we need  */

    if (cmd->numout > totwrote)
        padwrote = padtowrite = cmd->numout - totwrote;


    /* Write the new info file for the output data */

    if (!cmd->nobaryP) {
        idata.bary = 1;
        idata.mjd_i = (int) floor(blotoa - (barydispdt / SECPERDAY));
        idata.mjd_f = blotoa - (barydispdt / SECPERDAY) - idata.mjd_i;
    }
    if (cmd->downsamp > 1)
        idata.dt = dsdt;
    update_infodata(&idata, totwrote, padtowrite, diffbins, numdiffbins);
    writeinf(&idata);

    /* Set the padded points equal to the average data point */

    if (idata.numonoff >= 1) {
        int jj, index, startpad, endpad;

        for (ii = 0; ii < worklen; ii++)
            outdata[ii] = avg;
        fclose(outfile);
        outfile = chkfopen(datafilenm, "rb+");
        for (ii = 0; ii < idata.numonoff; ii++) {
            index = 2 * ii;
            startpad = idata.onoff[index + 1];
            if (ii == idata.numonoff - 1)
                endpad = idata.N - 1;
            else
                endpad = idata.onoff[index + 2];
            chkfseek(outfile, (startpad + 1) * sizeof(float), SEEK_SET);
            padtowrite = endpad - startpad;
            for (jj = 0; jj < padtowrite / worklen; jj++)
                chkfwrite(outdata, sizeof(float), worklen, outfile);
            chkfwrite(outdata, sizeof(float), padtowrite % worklen, outfile);
        }
    }
    vect_free(outdata);

    //  Close all the raw files and free their vectors
    close_rawfiles(&s);

    /* Print simple stats and results */

    var /= (datawrote - 1);

    /* Conver the '.dat' file to '.sdat' if requested */

    if (cmd->shortsP) {
        FILE *infile;
        int safe_convert = 1, bufflen = 65536;
        char *sdatafilenm;
        float *fbuffer;
        short *sbuffer;

        offset = (int) (floor(avg));
        if ((max - min) > (SHRT_MAX - SHRT_MIN)) {
            if ((max - min) < 1.5 * (SHRT_MAX - SHRT_MIN)) {
                printf("Warning:  There is more dynamic range in the data\n"
                       "          than can be handled perfectly:\n"
                       "               max - min = %.2f - %.2f = %.2f\n"
                       "          Clipping the low values...\n\n", max, min,
                       max - min);
                offset = max - SHRT_MAX;
            } else {
                printf("Error:  There is way too much dynamic range in the data:\n"
                       "               max - min = %.2f - %.2f = %.2f\n"
                       "        Not converting to shorts.\n\n", max, min, max - min);
                safe_convert = 0;
            }
        }

        if (safe_convert) {
            fbuffer = gen_fvect(bufflen);
            sbuffer = gen_svect(bufflen);
            sdatafilenm = (char *) calloc(slen, 1);
            sprintf(sdatafilenm, "%s.sdat", cmd->outfile);
            printf("\n\nConverting floats in '%s' to shorts in '%s'.",
                   datafilenm, sdatafilenm);
            fflush(NULL);

            infile = chkfopen(datafilenm, "rb");
            outfile = chkfopen(sdatafilenm, "wb");
            while ((numread = chkfread(fbuffer, sizeof(float), bufflen, infile))) {
                for (ii = 0; ii < numread; ii++)
                    sbuffer[ii] = (short) (fbuffer[ii] + 1e-20 - offset);
                chkfwrite(sbuffer, sizeof(short), numread, outfile);
            }
            fclose(infile);
            fclose(outfile);
            remove(datafilenm);
            vect_free(fbuffer);
            vect_free(sbuffer);
        }
    }

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
    if (cmd->shortsP && offset != 0)
        printf("          Offset applied to data:  %d\n", -offset);
    printf("\n");

    /* Cleanup */

    if (cmd->maskfileP) {
        free_mask(obsmask);
        vect_free(maskchans);
    }
    vect_free(tobsf);
    vect_free(dispdt);
    vect_free(idispdt);
    free(outinfonm);
    free(datafilenm);
    if (!cmd->nobaryP)
        vect_free(diffbins);
    return (0);
}


static int read_floats(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
    /* Read the raw data and return numbar read */

    return chkfread(data, sizeof(float), (numpts * numchan), file) / numchan;
}


static int read_shorts(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains short integer data.  */
/* The equivalent floats are placed in *data.               */
/* It returns the number of points read.                    */
{
    short *sdata;
    int ii, numread;

    sdata = (short *) malloc((size_t) (sizeof(short) * (numpts * numchan)));
    if (!sdata) {
        perror("\nError allocating short array in read_shorts()");
        printf("\n");
        exit(-1);
    }
    numread = chkfread(sdata, sizeof(short),
                       (unsigned long) (numpts * numchan), file) / numchan;
    for (ii = 0; ii < numread; ii++)
        data[ii] = (float) sdata[ii];
    free(sdata);
    return numread;
}

static int downsample(float outdata[], int numread, int downsampfact)
/* Downsample the floating point data by a factor downsampfact */
{
    float tmpout;
    int ii, jj, index;

    numread /= downsampfact;
    for (ii = 0, index = 0; ii < numread; ii++, index += downsampfact) {
        tmpout = 0.0;
        for (jj = 0; jj < downsampfact; jj++)
            tmpout += outdata[index + jj];
        outdata[ii] = tmpout;
    }
    return numread;
}


static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins)
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
