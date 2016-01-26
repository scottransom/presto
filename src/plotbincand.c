#include "presto.h"
#include "plot2d.h"

#define ZOOMFACT       10
#define ZOOMNEIGHBORS  20

/* Function definitions */

void usage(void);

/* Main routine */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *fftfile, *candfile = NULL, *psfile = NULL;
    char filenm[100], candnm[100], psfilenm[120];
    float locpow, norm, powargr, powargi;
    float *powr, *spreadpow, *minizoompow, *freqs;
    fcomplex *data, *minifft, *minizoom, *spread;
    fcomplex *resp, *kernel;
    double T, dr, ftobinp;
    int ii, nbins, ncands, candnum, lofreq = 0, nzoom, numsumpow = 1;
    int numbetween, numkern, kern_half_width;
    binaryprops binprops;
    infodata idata;

    if (argc < 3 || argc > 6) {
        usage();
        exit(1);
    }
    printf("\n\n");
    printf("         Binary Candidate Display Routine\n");
    printf("              by Scott M. Ransom\n\n");

    /* Initialize the filenames: */

    sprintf(filenm, "%s.fft", argv[1]);
    sprintf(candnm, "%s_bin.cand", argv[1]);

    /* Read the info file */

    readinf(&idata, argv[1]);
    if (strlen(remove_whitespace(idata.object)) > 0) {
        printf("Plotting a %s candidate from '%s'.\n", idata.object, filenm);
    } else {
        printf("Plotting a candidate from '%s'.\n", filenm);
    }
    T = idata.N * idata.dt;

    /* Open the FFT file and get its length */

    fftfile = chkfopen(filenm, "rb");
    nbins = chkfilelen(fftfile, sizeof(fcomplex));

    /* Open the candidate file and get its length */

    candfile = chkfopen(candnm, "rb");
    ncands = chkfilelen(candfile, sizeof(binaryprops));

    /* The candidate number to examine */

    candnum = atoi(argv[2]);

    /* Check that candnum is in range */

    if ((candnum < 1) || (candnum > ncands)) {
        printf("\nThe candidate number is out of range.\n\n");
        exit(1);
    }
    /* The lowest freq present in the FFT file */

    if (argc >= 4) {
        lofreq = atoi(argv[3]);
        if ((lofreq < 0) || (lofreq > nbins - 1)) {
            printf("\n'lofreq' is out of range.\n\n");
            exit(1);
        }
    }
    /* Is the original FFT a sum of other FFTs with the amplitudes added */
    /* in quadrature?  (i.e. an incoherent sum)                          */

    if (argc >= 5) {
        numsumpow = atoi(argv[4]);
        if (numsumpow < 1) {
            printf("\nNumber of summed powers must be at least one.\n\n");
            exit(1);
        }
    }
    /* Initialize PGPLOT using Postscript if requested  */

    if ((argc == 6) && (!strcmp(argv[5], "ps"))) {
        sprintf(psfilenm, "%s_bin_cand_%d.ps", argv[1], candnum);
        cpgstart_ps(psfilenm, "landscape");
    } else {
        cpgstart_x("landscape");
    }

    /* Read the binary candidate */

    chkfileseek(candfile, (long) (candnum - 1), sizeof(binaryprops), SEEK_SET);
    chkfread(&binprops, sizeof(binaryprops), 1, candfile);
    fclose(candfile);

    /* Output the binary candidate */

    print_bin_candidate(&binprops, 2);

    /* Allocate some memory */

    powr = gen_fvect(binprops.nfftbins);
    minifft = gen_cvect(binprops.nfftbins / 2);
    spread = gen_cvect(binprops.nfftbins);
    spreadpow = gen_fvect(binprops.nfftbins);
    nzoom = 2 * ZOOMFACT * ZOOMNEIGHBORS;
    minizoom = gen_cvect(nzoom);
    minizoompow = gen_fvect(nzoom);

    /* Allocate and initialize our interpolation kernel */

    numbetween = 2;
    kern_half_width = r_resp_halfwidth(LOWACC);
    numkern = 2 * numbetween * kern_half_width;
    resp = gen_r_response(0.0, numbetween, numkern);
    kernel = gen_cvect(binprops.nfftbins);
    place_complex_kernel(resp, numkern, kernel, binprops.nfftbins);
    COMPLEXFFT(kernel, binprops.nfftbins, -1);

    /* Read the data from the FFT file */

    data = read_fcomplex_file(fftfile, binprops.lowbin - lofreq, binprops.nfftbins);

    /* Turn the Fourier amplitudes into powers */

    for (ii = 0; ii < binprops.nfftbins; ii++)
        powr[ii] = POWER(data[ii].r, data[ii].i);

    /* Chop the powers that are way above the median level */

    prune_powers(powr, binprops.nfftbins, numsumpow);

    /* Perform the minifft */

    memcpy((float *) minifft, powr, sizeof(float) * binprops.nfftbins);
    realfft((float *) minifft, binprops.nfftbins, -1);

    /* Calculate the normalization constant */

    norm = sqrt((double) binprops.nfftbins * (double) numsumpow) / minifft[0].r;
    locpow = minifft[0].r / binprops.nfftbins;

    /* Divide the original power spectrum by the local power level */

    for (ii = 0; ii < binprops.nfftbins; ii++)
        powr[ii] /= locpow;

    /* Now normalize the miniFFT */

    minifft[0].r = 1.0;
    minifft[0].i = 1.0;
    for (ii = 1; ii < binprops.nfftbins / 2; ii++) {
        minifft[ii].r *= norm;
        minifft[ii].i *= norm;
    }

    /* Interpolate the minifft and convert to power spectrum */

    corr_complex(minifft, binprops.nfftbins / 2, RAW,
                 kernel, binprops.nfftbins, FFT,
                 spread, binprops.nfftbins, kern_half_width,
                 numbetween, kern_half_width, CORR);
    for (ii = 0; ii < binprops.nfftbins; ii++)
        spreadpow[ii] = POWER(spread[ii].r, spread[ii].i);

    /* Plot the initial data set */

    freqs = gen_freqs(binprops.nfftbins, binprops.lowbin / T, 1.0 / T);
    xyline(binprops.nfftbins, freqs, powr, "Pulsar Frequency (hz)",
           "Power / Local Power", 1);
    vect_free(freqs);
    printf("The initial data set (with high power outliers removed):\n\n");

    /* Plot the miniFFT */

    freqs = gen_freqs(binprops.nfftbins, 0.0, T / (2 * binprops.nfftbins));
    xyline(binprops.nfftbins, freqs, spreadpow, "Binary Period (sec)",
           "Normalized Power", 1);
    vect_free(freqs);
    printf("The miniFFT:\n\n");

    /* Interpolate and plot the actual candidate peak */

    ftobinp = T / binprops.nfftbins;
    freqs = gen_freqs(nzoom, (binprops.rdetect - ZOOMNEIGHBORS) *
                      ftobinp, ftobinp / (double) ZOOMFACT);
    for (ii = 0; ii < nzoom; ii++) {
        dr = -ZOOMNEIGHBORS + (double) ii / ZOOMFACT;
        rz_interp(minifft, binprops.nfftbins / 2, binprops.rdetect + dr,
                  0.0, kern_half_width, &minizoom[ii]);
        minizoompow[ii] = POWER(minizoom[ii].r, minizoom[ii].i);
    }
    xyline(nzoom, freqs, minizoompow, "Binary Period (sec)", "Normalized Power", 1);
    vect_free(freqs);
    printf("The candidate itself:\n\n");
    printf("Done.\n\n");

    /* Cleanup */

    cpgend();
    vect_free(data);
    vect_free(powr);
    vect_free(resp);
    vect_free(kernel);
    vect_free(minifft);
    vect_free(spread);
    vect_free(spreadpow);
    vect_free(minizoom);
    vect_free(minizoompow);
    fclose(fftfile);
    if ((argc == 6) && (!strcmp(argv[5], "ps"))) {
        fclose(psfile);
    }
    return (0);
}

void usage(void)
{
    printf("\nUsage:  plotbincand file candnum [lofreq] [ps]\n\n");
    printf("  Mandatory arguments:\n");
    printf("       'file' = a string containing the FFT file's name.\n");
    printf("                You must have an '.inf' file of the same\n");
    printf("                name as well.  Do not add a suffix.\n");
    printf("    'candnum' = (int) The candidate number from search_bin\n");
    printf("                to examine.\n");
    printf("  Optional arguments:\n");
    printf("     'lofreq' = (int) Lowest Fourier bin in FFT file.  This is\n");
    printf("                useful if you chop a long FFT for space reasons.\n");
    printf("                If 'lofreq' is present and not equal to '0', \n");
    printf("                we will assume the 0th frequency bin = 1 for\n");
    printf("                normalization purposes.\n");
    printf("         'ps' = Send graphics output to a Postscript file\n");
    printf("                instead of the screen.\n");
    printf("        'sum' = Number of summed FFTs (added in quadrature)\n");
    printf("                that make up the original FFT to search.\n\n");
    printf("\n\n");
}
