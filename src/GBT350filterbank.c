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
    char *path, *filenm;

    split_path_file(filenmbase, &path, &filenm);
    strncpy(fb->inpfile, filenm, 40);
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
    if (h_or_d < 0) {
        h_or_d = abs(h_or_d);
        fb->src_dej = h_or_d * 10000.0 + m * 100.0 + s;
        fb->src_dej *= -1.0;
    } else {
        fb->src_dej = h_or_d * 10000.0 + m * 100.0 + s;
    }
    fb->az_start = 0.0;
    fb->za_start = 0.0;
    fb->nchans = spigot->lags_per_sample;
    fb->foff = spigot->bandwidth / fb->nchans;
    fb->fch1 = spigot->freq_ctr + (fb->nchans / 2 - 0.5) * fb->foff;
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
    free(path);
    free(filenm);
}


int main(int argc, char *argv[])
{
    FILE **infiles, *outfile = NULL, *offsetfile = NULL;
    int filenum, argnum = 1, ii = 0, ptsperblock, numlags, numfiles;
    int bytes_per_read, offset = 0, numoffsets, output = 1;
    long long N;
    char *path, *filenm;
    char outfilenm[200], filenmbase[200], offsetnm[200], rawlags[4096];
    unsigned char output_samples[2048];
    float *offsets = NULL;
    double dt, T;
    SPIGOT_INFO *spigots;
    sigprocfb fb;
    infodata idata;

    if (argc == 1) {
        fprintf(stderr, "Usage: GBT350filterbank SPIGOT_fits_files\n");
        exit(0);
    }

    /* Determine the filename base in order to read the offsets */
    strncpy(filenmbase, argv[argnum], strlen(argv[argnum]) - 10);
    filenmbase[strlen(argv[argnum]) - 10] = '\0';

    /* Attempt to read a file with lag offsets in it */

    /* First try the same path as the input file */
    sprintf(offsetnm, "%s.offset", filenmbase);
    offsetfile = fopen(offsetnm, "rb");
    if (!offsetfile) {
        /* Now try the current working directory */
        split_path_file(filenmbase, &path, &filenm);
        sprintf(offsetnm, "%s.offset", filenm);
        offsetfile = fopen(offsetnm, "rb");
        free(path);
        free(filenm);
    }
    if (offsetfile) {
        /* Determine the length of the file */
        numoffsets = (int) chkfilelen(offsetfile, sizeof(float));
        /* Create the array and read 'em */
        offsets = gen_fvect(numoffsets);
        chkfread(offsets, sizeof(float), numoffsets, offsetfile);
        offset = 1;
        /* close the offset file */
        fclose(offsetfile);
        printf("Offseting the spectra with the %d values found in '%s'\n\n",
               numoffsets, offsetnm);
    }

    /* Read and convert the basic SPIGOT file information */

    numfiles = argc - 1;
    printf("Spigot card input file information:\n");
    spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
    infiles = (FILE **) malloc(sizeof(FILE *) * numfiles);
    for (filenum = 0; filenum < numfiles; filenum++, argnum++) {
        printf("  '%s'\n", argv[argnum]);
        infiles[filenum] = chkfopen(argv[argnum], "rb");
        read_SPIGOT_header(argv[argnum], spigots + filenum);
        rewind(infiles[filenum]);
    }
    printf("\n");

    /* The following is necessary in order to initialize all the */
    /* static variables in spigot.c                              */
    get_SPIGOT_file_info(infiles, spigots, numfiles, 0, 0, &N,
                         &ptsperblock, &numlags, &dt, &T, &idata, output);
    spigot2sigprocfb(&(spigots[0]), &fb, filenmbase);
    fb.N = N;
    free(spigots);

    /* Determine the output file name */
    sprintf(outfilenm, "GBT350_%s_%.0f_%04d.fil",
            spigots[0].object, floor(spigots[0].MJD_obs), spigots[0].scan_number);
    printf("Writing data to file '%s'.\n\n", outfilenm);

    /* Write the header */
    outfile = chkfopen(outfilenm, "wb");
    write_filterbank_header(&fb, outfile);

    /* Step through the SPIGOT files */
    ii = 0;
    for (filenum = 0, argnum = 1; filenum < numfiles; filenum++, argnum++) {
        printf("Reading from file '%s'...\n", argv[argnum]);
        chkfseek(infiles[filenum], spigots[0].header_len, SEEK_SET);
        bytes_per_read = spigots[0].lags_per_sample * spigots[0].bits_per_lag / 8;

        /* Loop over the samples in the file */
        while (chkfread(rawlags, bytes_per_read, 1, infiles[filenum])) {
            if (offset)
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS, offsets[ii]);
            else
                convert_SPIGOT_point(rawlags, output_samples, SUMIFS, 1.0);
            ii++;
            /* Invert the band so that the high freqs are first */
            /* This is how SIGPROC stores its data.             */
            {
                int jj;
                unsigned char tempzz = 0.0, *loptr, *hiptr;
                loptr = output_samples + 0;
                hiptr = output_samples + fb.nchans - 1;
                for (jj = 0; jj < fb.nchans / 2; jj++, loptr++, hiptr--) {
                    SWAP(*loptr, *hiptr);
                }
            }
            chkfwrite(output_samples, sizeof(unsigned char), fb.nchans, outfile);
        }
        fclose(infiles[filenum]);
    }
    if (offset)
        vect_free(offsets);
    fclose(outfile);
    printf("Converted and wrote %d samples.\n\n", ii);
    free(infiles);
    return 0;
}
