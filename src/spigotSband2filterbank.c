#include "presto.h"
#include "mask.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "fitsfile.h"
#include "fitshead.h"

void spigot2sigprocfb(SPIGOT_INFO * spigot, sigprocfb * fb, char *filenmbase)
{
    int h_or_d, m;
    double s;

    strncpy(fb->inpfile, filenmbase, 40);
    strncpy(fb->source_name, spigot->object, 80);
    fb->nifs = 1;
    if (spigot->num_samplers == 1)
        strncpy(fb->ifstream, "YXXX", 8);
    else if (spigot->num_samplers == 2)
        strncpy(fb->ifstream, "YYXX", 8);
    fb->tstart = spigot->MJD_obs + spigot->elapsed_time / SECPERDAY;
    fb->tsamp = spigot->dt_us / 1e6;
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
    fb->foff = -fb->foff;
    // Here is where we chop the bottom 1/4 (200MHz) of the band
    fb->nchans = fb->nchans - fb->nchans / 4;
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
}


int main(int argc, char *argv[])
{
    FILE **infiles, *outfile = NULL, *scalingfile;
    int filenum, argnum = 1, ii = 0, ptsperblock, numlags, numfiles;
    int bytes_per_read, scaling = 0, numscalings, output = 1;
    long long N;
    char *path, *filenm;
    char outfilenm[200], filenmbase[200], scalingnm[200], rawlags[4096];
    unsigned char output_samples[2048];
    float *scalings = NULL;
    double dt, T;
    SPIGOT_INFO *spigots;
    sigprocfb fb;
    infodata idata;

    if (argc == 1) {
        fprintf(stderr, "Usage: spigotSband2filterbank [-stdout] SPIGOT_files\n");
        exit(0);
    }

    if (!strcmp(argv[argnum], "-stdout")) {     /* Use STDOUT */
        argnum++;
        output = 0;
        outfile = stdout;
    } else {
        printf
            ("\nConverting raw SPIGOT S-band FITs data into SIGPROC\n"
             "filterbank format and throwing out the bottom 200MHz...\n\n");
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
        if (outfile != stdout)
            printf("Scaling the lags with the %d values found in '%s'\n\n",
                   numscalings, scalingnm);
    }

    /* Read and convert the basic SPIGOT file information */

    numfiles = argc - 1;
    if (outfile == stdout)
        numfiles--;
    else
        printf("Spigot card input file information:\n");
    spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
    infiles = (FILE **) malloc(sizeof(FILE *) * numfiles);
    for (filenum = 0; filenum < numfiles; filenum++, argnum++) {
        if (outfile != stdout)
            printf("  '%s'\n", argv[argnum]);
        infiles[filenum] = chkfopen(argv[argnum], "rb");
        read_SPIGOT_header(argv[argnum], spigots + filenum);
        rewind(infiles[filenum]);
    }
    if (outfile != stdout)
        printf("\n");

    /* The following is necessary in order to initialize all the */
    /* static variables in spigot.c                              */
    get_SPIGOT_file_info(infiles, spigots, numfiles, 0, 0, &N,
                         &ptsperblock, &numlags, &dt, &T, &idata, output);

    /* Step through the SPIGOT files */
    ii = 0;
    if (outfile == stdout)
        argnum = 2;
    else
        argnum = 1;
    for (filenum = 0; filenum < numfiles; filenum++, argnum++) {
        split_path_file(argv[argnum], &path, &filenm);
        strncpy(filenmbase, filenm, strlen(filenm) - 5);
        filenmbase[strlen(filenm) - 5] = '\0';
        sprintf(outfilenm, "%s.fil", filenmbase);
        if (outfile != stdout) {
            printf("Reading S-band Spigot lags from '%s'\n", argv[argnum]);
            printf("Writing filterbank spectra to   '%s'\n\n", outfilenm);
            outfile = chkfopen(outfilenm, "wb");
        }
        // Update the filterbank header information for each file
        spigot2sigprocfb(&(spigots[filenum]), &fb, filenmbase);
        write_filterbank_header(&fb, outfile);
        chkfseek(infiles[filenum], spigots[filenum].header_len, SEEK_SET);
        bytes_per_read = numlags * spigots[filenum].bits_per_lag / 8;

        /* Loop over the samples in the file */
        while (chkfread(rawlags, bytes_per_read, 1, infiles[filenum])) {
            if (scaling)
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS, scalings[ii]);
            else
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS, 1.0);
            ii++;
            /* Invert the band so that the high freqs are first */
            /* This is how SIGPROC stores its data.             */
            {
                int jj;
                unsigned char tempzz = 0.0, *loptr, *hiptr;
                loptr = output_samples + 0;
                hiptr = output_samples + numlags - 1;
                for (jj = 0; jj < numlags / 2; jj++, loptr++, hiptr--) {
                    SWAP(*loptr, *hiptr);
                }
            }
            chkfwrite(output_samples, sizeof(unsigned char), fb.nchans, outfile);
        }
        fclose(infiles[filenum]);
        if (outfile != stdout)
            fclose(outfile);
    }
    if (outfile != stdout)
        fprintf(stderr, "Converted and wrote %d samples.\n\n", ii);
    if (scaling)
        vect_free(scalings);
    free(spigots);
    free(path);
    free(filenm);
    free(infiles);
    return 0;
}
