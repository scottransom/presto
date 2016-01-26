#include "presto.h"
#include "mask.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "fitsfile.h"
#include "fitshead.h"
#include "spigot2filterbank_cmd.h"
#include "slalib.h"

extern void get_calibrated_lags(void *rawlags, float *calibrated_lags);

/* extern void sla_obs_(int *N, char *scope, char *name,  */
/*                      double *lon, double *lat, double *hgt); */
/* extern void sla_oap_(char *type, double *azimuth, double *zendist,  */
/*                      double *MJD, double *dut, double *lon, double *lat,  */
/*                      double *hgt, double *xp, double *yp, double *temp,  */
/*                      double *atm, double *humid, double *microns,  */
/*                      double *tlr, double *rap, double *dap); */
/* extern void sla_amp_(double *rap, double *dap, double *MJD,  */
/*                      double *equinox, double *ramean, double *decmean); */

void spigot2sigprocfb(SPIGOT_INFO * spigot, sigprocfb * fb, char *filenmbase,
                      int lokill, int hikill, int downsamp,
                      int update_posn, double time_offset)
{
    int h_or_d, m;
    double s, dt;

    /* Set the time offset for the posn calc */
    if (update_posn) {
        dt = time_offset;
    } else {
        dt = 0.0;
    }
    strncpy(fb->inpfile, filenmbase, 40);
    strncpy(fb->source_name, spigot->object, 80);
    fb->nifs = 1;
    if (spigot->num_samplers == 1)
        strncpy(fb->ifstream, "YXXX", 8);
    else if (spigot->num_samplers == 2)
        strncpy(fb->ifstream, "YYXX", 8);
    fb->tstart = spigot->MJD_obs + spigot->elapsed_time / SECPERDAY;
    fb->tsamp = spigot->dt_us / 1e6 * downsamp;
    hours2hms(spigot->ra / 15.0, &h_or_d, &m, &s);
    fb->src_raj = h_or_d * 10000.0 + m * 100.0 + s;
    deg2dms(spigot->dec, &h_or_d, &m, &s);
    fb->src_dej = abs(h_or_d) * 10000.0 + abs(m) * 100.0 + fabs(s);
    if (spigot->dec < 0)
        fb->src_dej = -fb->src_dej;
    fb->az_start = spigot->az;
    fb->za_start = 90.0 - spigot->el;
    fb->nchans = spigot->lags_per_sample;
    fb->foff = spigot->bandwidth / fb->nchans;
    fb->fch1 = spigot->freq_ctr + (fb->nchans / 2 - 0.5) * fb->foff;
    fb->fch1 -= fb->foff * hikill;
    fb->nchans -= (hikill + lokill);
    fb->foff = -fb->foff;
    fb->machine_id = 7;
    fb->telescope_id = 6;
    fb->nbits = 8;
    fb->sumifs = spigot->summed_pols;
    if (fb->sumifs)
        fb->nifs = 1;
    else {
        if (spigot->num_samplers == 2)
            fb->nifs = 2;
        else
            fb->nifs = 1;
    }
    /* The following are not necessary for writing filterbank files */
    fb->headerlen = 0;
    fb->N = spigot->samples_per_file;

    /* Update the position if the GBT was not tracking */
    /* (i.e. for the driftscan surveys)                */
    if (update_posn && !spigot->tracking) {
        int N = 0;
        char scope[10] = { "GBT" }, type[] = {
        "A"};
        char name[40];
        double MJD, lon, lat, hgt, microns, az, zd, rap, dap, rmn, dmn;
        double dtmp = 0.0, atm = 1010.0, temp = 283.0;
        double humid = 0.5, tlr = 0.0065, eq = 2000.0;

        /* Compute the RA/DEC using SLALIB from the Az/El */

        slaObs(N, scope, name, &lon, &lat, &hgt);
        if (fabs(hgt - 880.0) > 0.01) {
            printf("Warning!:  SLALIB is not correctly identifying the GBT!\n\n");
        }
        //printf("slalib: %d '%s' '%s' %f  %f  %f\n", 
        //      N, scope, name, lon, lat, hgt);
        lon = -lon;
        az = fb->az_start * DEGTORAD;
        zd = fb->za_start * DEGTORAD;
        microns = 3e8 / (spigot->freq_ctr * 1e6) * 1e6;
        MJD = fb->tstart + dt / 86400.0;
        slaOap(type, az, zd, MJD, dtmp, lon, lat, hgt,
               dtmp, dtmp, temp, atm, humid, microns, tlr, &rap, &dap);
        //printf("slalib:  %.15f  %.15f\n", rap, dap);
        slaAmp(rap, dap, MJD, eq, &rmn, &dmn);
        //printf("slalib:  %.15f  %.15f\n", rmn, dmn);

        /* Now update the positions */
        hours2hms(rmn * RADTODEG / 15.0, &h_or_d, &m, &s);
        fb->src_raj = h_or_d * 10000.0 + m * 100.0 + s;
        deg2dms(dmn * RADTODEG, &h_or_d, &m, &s);
        fb->src_dej = abs(h_or_d) * 10000.0 + abs(m) * 100.0 + fabs(s);
        if (dmn < 0)
            fb->src_dej = -fb->src_dej;
    }
}


int main(int argc, char *argv[])
{
    FILE **infiles, *outfile = NULL, *scalingfile, *zerolagfile = NULL;
    int filenum, ptsperblock, numlags, numfiles, outnumlags;
    int bytes_per_read, scaling = 0, numscalings, firstfile, firstspec;
    int write_data = 1, update_posn = 0;
    long long N, numconverted = 0, numwritten = 0;
    char *path, *filenm, rawlags[4096];
    char outfilenm[200], filenmbase[200], zerolagnm[200], scalingnm[200];
    unsigned char output_samples[2048];
    float *scalings = NULL, lags[2048], tmp_floats[2048];
    double dt, T, time_offset;
    SPIGOT_INFO *spigots;
    sigprocfb fb;
    infodata idata;
    Cmdline *cmd;

    /* Call usage() if we have no command line arguments */
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }

    /* Parse the command line using the excellent program Clig */
    cmd = parseCmdline(argc, argv);
    // showOptionValues();
    numfiles = cmd->argc;

    /* If requesting that we skip values or only write          */
    /* a specific number of values, require an output file name */
    if ((cmd->skip > 0) || cmd->numoutP) {
        if (!cmd->outfileP && !cmd->stdoutP) {
            fprintf(stderr,
                    "\nspigot2filterbank ERROR:  You need to specify an output\n"
                    "     filename (or -stdout) when using the -skip and/or \n"
                    "     -numout options!\n\n");
            exit(1);
        }
    }

    /* Give an error if specifying stdout and an output file */
    if (cmd->stdoutP && cmd->outfileP) {
        fprintf(stderr,
                "\nspigot2filterbank ERROR:  You cannot specify both an output\n"
                "     file and that the data should go to STDOUT!\n\n");
        exit(1);
    }

    if (!cmd->stdoutP)
        printf
            ("\nConverting SPIGOT FITs lags into SIGPROC filterbank format...\n\n");

    /* Determine the filename base from the first spigot file */
    split_path_file(cmd->argv[0], &path, &filenm);
    strncpy(filenmbase, filenm, strlen(filenm) - 5);
    filenmbase[strlen(filenm) - 5] = '\0';
    free(path);
    free(filenm);

    /* Open a file to store the zerolags */
    if (cmd->zerolagsP) {
        if (cmd->outfileP) {
            sprintf(zerolagnm, "%s.zerolags", cmd->outfile);
        } else {
            sprintf(zerolagnm, "%s.zerolags", filenmbase);
        }
        zerolagfile = chkfopen(zerolagnm, "wb");
    }

    /* Attempt to read a file with lag scalings in it */
    sprintf(scalingnm, "%s.scaling", filenmbase);
    if ((scalingfile = fopen(scalingnm, "rb"))) {
        /* Determine the length of the file */
        numscalings = (int) chkfilelen(scalingfile, sizeof(float));
        /* Create the array and read 'em */
        scalings = gen_fvect(numscalings);
        chkfread(scalings, sizeof(float), numscalings, scalingfile);
        scaling = 1;
        /* close the scaling file */
        fclose(scalingfile);
        if (!cmd->stdoutP)
            printf("Scaling the lags with the %d values found in '%s'\n\n",
                   numscalings, scalingnm);
    }

    /* Read and convert the basic SPIGOT file information */
    if (!cmd->stdoutP)
        printf("Spigot card input file information:\n");
    spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
    infiles = (FILE **) malloc(sizeof(FILE *) * numfiles);
    for (filenum = 0; filenum < numfiles; filenum++) {
        if (!cmd->stdoutP)
            printf("  '%s'\n", cmd->argv[filenum]);
        infiles[filenum] = chkfopen(cmd->argv[filenum], "rb");
        read_SPIGOT_header(cmd->argv[filenum], spigots + filenum);
        rewind(infiles[filenum]);
    }
    if (!cmd->stdoutP)
        printf("\n");

    /* The following is necessary in order to initialize all the */
    /* static variables in spigot.c                              */
    get_SPIGOT_file_info(infiles, spigots, numfiles, 0, 0, &N,
                         &ptsperblock, &numlags, &dt, &T, &idata, !cmd->stdoutP);

    /* Compute the first file required */
    firstfile = cmd->skip / spigots[0].samples_per_file;
    firstspec = cmd->skip % spigots[0].samples_per_file;
    if (!cmd->numoutP) {
        cmd->numout = (N - cmd->skip) / cmd->downsamp;
    }
    bytes_per_read = numlags * spigots[0].bits_per_lag / 8;
    outnumlags = numlags - cmd->lokill - cmd->hikill;

    /* Step through the SPIGOT files */
    filenum = firstfile;
    while (numwritten < cmd->numout) {
        split_path_file(cmd->argv[filenum], &path, &filenm);
        strncpy(filenmbase, filenm, strlen(filenm) - 5);
        filenmbase[strlen(filenm) - 5] = '\0';
        if (cmd->outfileP) {
            if (filenum == firstfile) {
                sprintf(outfilenm, "%s", cmd->outfile);
                outfile = chkfopen(outfilenm, "wb");
            }
        } else {
            sprintf(outfilenm, "%s.fil", filenmbase);
            if (cmd->stdoutP) {
                outfile = stdout;
            } else {
                outfile = chkfopen(outfilenm, "wb");
            }
        }
        if (!cmd->stdoutP) {
            if (filenum == firstfile)
                printf("Reading Spigot lags from '%s' (starting at sample %d)\n",
                       cmd->argv[filenum], firstspec);
            else
                printf("Reading Spigot lags from '%s'\n", cmd->argv[filenum]);
            printf("Writing filterbank spectra to   '%s'\n\n", outfilenm);
        }

        /* Update the filterbank header information for each file */
        if (!spigots[filenum].tracking) {
            if (cmd->outfileP || cmd->stdoutP) {
                time_offset = 0.5 * cmd->numout *
                    spigots[filenum].dt_us * 1e-6 * cmd->downsamp;
            } else {            // Just normal files
                time_offset = 0.5 * spigots[filenum].file_duration;
            }
            update_posn = 1;
        } else {
            update_posn = 0;
            time_offset = 0.0;
        }
        /* Adjust the Spigot start time for the skip */
        spigots[filenum].elapsed_time += cmd->skip * spigots[filenum].dt_us * 1e-6;
        /* Determine the SIGPROC header */
        spigot2sigprocfb(&(spigots[filenum]), &fb, filenmbase,
                         cmd->lokill, cmd->hikill, cmd->downsamp,
                         update_posn, time_offset);
        /* Correct the structure if we are using floats */
        if (cmd->floatsP)
            fb.nbits = 32;

        /* Write a filterbank header if we have not been told not to.   */
        /* Don't write it, though, if using stdio or a specified output */
        /* file and the input file is not the first.                    */
        if (!cmd->nohdrP) {
            if ((!cmd->stdoutP && !cmd->outfileP) ||
                ((cmd->stdoutP || cmd->outfileP) && filenum == firstfile)) {
                write_filterbank_header(&fb, outfile);
            }
        }

        /* Go to the correct autocorrelation in the correct first FITs file */
        chkfseek(infiles[filenum],
                 spigots[filenum].header_len + bytes_per_read * firstspec, SEEK_SET);

        /* Loop over the samples in the file */
        while ((numwritten < cmd->numout) &&
               (chkfread(rawlags, bytes_per_read, 1, infiles[filenum]))) {
            if (cmd->zerolagsP) {
                /* Correct the lags so we can write the zerolag */
                get_calibrated_lags(rawlags, lags);
                chkfwrite(lags, sizeof(float), 1, zerolagfile);
            }
            if (scaling) {
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS,
                                     scalings[cmd->skip + numconverted]);
            } else {
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS, 1.0);
            }
            /* If downsampling, average the current spectra */
            if (cmd->downsamp > 1) {
                int ii;
                if (numconverted % cmd->downsamp == 0) {
                    write_data = 0;
                    /* Zero the array used for averaging */
                    for (ii = 0; ii < outnumlags; ii++)
                        tmp_floats[ii] = 0.0;
                } else {
                    /* Add the current data to the array used for averaging */
                    for (ii = 0; ii < outnumlags; ii++)
                        tmp_floats[ii] += (float) output_samples[ii + cmd->lokill];
                    /* If that was the last sample to be added, average them */
                    /* and put them back into output_samples */
                    if (numconverted % cmd->downsamp == (cmd->downsamp - 1)) {
                        write_data = 1;
                        for (ii = 0; ii < outnumlags; ii++) {
                            tmp_floats[ii] /= (float) cmd->downsamp;
                            output_samples[ii + cmd->lokill] =
                                (unsigned char) (tmp_floats[ii] + 0.5);
                        }
                    }
                }
            }
            numconverted++;

            if (write_data) {
                int ii;

                /* Invert the band so that the high freqs are first */
                /* This is how SIGPROC stores its data.             */
                if (cmd->floatsP) {
                    float tempzz = 0.0, *loptr, *hiptr;
                    loptr = tmp_floats; //  killed channels are already gone
                    hiptr = tmp_floats + outnumlags - 1;
                    for (ii = 0; ii < outnumlags / 2; ii++, loptr++, hiptr--) {
                        SWAP(*loptr, *hiptr);
                    }
                } else {
                    unsigned char tempzz = 0.0, *loptr, *hiptr;
                    loptr = output_samples + cmd->lokill;
                    hiptr = output_samples + cmd->lokill + outnumlags - 1;
                    for (ii = 0; ii < outnumlags / 2; ii++, loptr++, hiptr--) {
                        SWAP(*loptr, *hiptr);
                    }
                }

                /* Now actually write the data */
                if (cmd->floatsP) {
                    /* Copy the bytes to floats */
                    for (ii = 0; ii < outnumlags; ii++)
                        tmp_floats[ii] = (float) output_samples[ii + cmd->lokill];
                    chkfwrite(tmp_floats, sizeof(float), fb.nchans, outfile);
                } else {
                    chkfwrite(output_samples + cmd->lokill,
                              sizeof(unsigned char), fb.nchans, outfile);
                }
                numwritten++;
            }
        }
        if ((!cmd->stdoutP) && (!cmd->outfileP))
            fclose(outfile);
        fclose(infiles[filenum]);
        firstspec = 0;
        filenum++;
        free(path);
        free(filenm);
    }
    if ((!cmd->stdoutP) && (cmd->outfileP))
        fclose(outfile);
    if (cmd->zerolagsP)
        fclose(zerolagfile);
    if (!cmd->stdoutP)
        fprintf(stderr, "Converted %lld samples and wrote %lld.\n\n",
                numconverted, numwritten);
    if (scaling)
        vect_free(scalings);
    free(spigots);
    free(infiles);
    return 0;
}
