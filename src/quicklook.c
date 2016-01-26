#include "presto.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *fftfile;
    double flook, dt, nph, maxz, pwr, hipow = 0.0;
    double zlo = -30.0, zhi = 30.0, dr, dz = 2.0;
    double hir = 0.0, hiz = 0.0, newhir, newhiz;
    fcomplex **ffdotplane, *data;
    float powargr, powargi;
    int startbin, numdata, nextbin, nr, nz, numkern;
    int i, j, kernel_half_width, numbetween = 4;
    int n, corrsize = 1024;
    char filenm[80], compare[200];
    rderivs derivs;
    fourierprops props;
    infodata idata;
    struct tms runtimes;
    double ttim, utim, stim, tott;


    tott = times(&runtimes) / (double) CLK_TCK;
    if (argc != 3) {
        printf("\nUsage:  'quicklook filename fftfreq dt'\n\n");
        printf("   'filename' = a string containing the FFT file's name.\n");
        printf("                (do not include the '.fft' suffix)\n");
        printf("    'fftfreq' = the central fourier frequency to examine.\n");
        printf("  Quicklook will search a region of the f-fdot plane\n");
        printf("  of a file containing a long, single precision FFT\n");
        printf("  using the Correlation method (i.e. Ransom and \n");
        printf("  Eikenberry, 1997, unpublished as of yet).\n");
        printf("  The search uses a spacing of 0.5 frequency bins in\n");
        printf("  the fourier frequency (r) direction, and 2 'bins' in\n");
        printf("  the fdot (z) direction.  The routine will output\n");
        printf("  statistics for the best candidate in the region.\n\n");
        printf("  The routine was written as a quick but useful hack\n");
        printf("  to show the power of the Correlation method, and the\n");
        printf("  forthcoming power of Scott Ransom's Pulsar Finder\n");
        printf("  Software.\n");
        printf("                                        2 March 2001\n\n");
        exit(0);
    }
    printf("\n\n");
    printf("  Quick-Look Pulsation Search\n");
    printf("     With database lookup.\n");
    printf("      by Scott M. Ransom\n");
    printf("         2 March, 2001\n\n");

    /*  Initialize our data: */

    sprintf(filenm, "%s.fft", argv[1]);
    readinf(&idata, argv[1]);
    flook = atof(argv[2]);
    dt = idata.dt;
    dr = 1.0 / (double) numbetween;
    fftfile = chkfopen(filenm, "r");
    nph = get_numphotons(fftfile);
    n = chkfilelen(fftfile, sizeof(float));
    nz = (int) ((zhi - zlo) / dz) + 1;

    /* Determine our starting frequency and get the data */

    maxz = (fabs(zlo) < fabs(zhi)) ? zhi : zlo;
    kernel_half_width = z_resp_halfwidth(maxz, HIGHACC);
    numkern = 2 * numbetween * kernel_half_width;
    while (numkern > 2 * corrsize)
        corrsize *= 2;
    startbin = (int) (flook) - corrsize / (2 * numbetween);
    numdata = corrsize / numbetween;
    data = read_fcomplex_file(fftfile, startbin, numdata);

    /*  Do the f-fdot plane correlations: */

    ffdotplane = corr_rz_plane(data, numdata, numbetween, kernel_half_width,
                               zlo, zhi, nz, corrsize, LOWACC, &nextbin);
    nr = corrsize - 2 * kernel_half_width * numbetween;

    /*  Search the resulting data set: */

    for (i = 0; i < nz; i++) {
        for (j = 0; j < nr; j++) {
            pwr = POWER(ffdotplane[i][j].r, ffdotplane[i][j].i);
            if (pwr > hipow) {
                hir = j * dr + kernel_half_width;
                hiz = i * dz + zlo;
                hipow = pwr;
            }
        }
    }

    /*  Maximize the best candidate: */

    hipow = max_rz_arr(data, numdata, hir, hiz, &newhir, &newhiz, &derivs);
    newhir += startbin;
    calc_props(derivs, newhir, newhiz, 0.0, &props);

    printf("Searched %d pts ", nz * nr);
    printf("(r: %.1f to %.1f, ", (double) startbin + kernel_half_width,
           (double) startbin + kernel_half_width + dr * (nr - 1));
    printf("z: %.1f to %.1f)\n\n", zlo, zhi);

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n",
           ttim, utim, stim);
    printf("  Total time: %.3f sec\n\n", tott);

    printf("The best candidate is:\n");

    print_candidate(&props, dt, n, nph, 2);
    comp_psr_to_cand(&props, &idata, compare, 1);
    printf("%s\n", compare);
    fclose(fftfile);

    /* Cleanup and exit */

    vect_free(ffdotplane[0]);
    vect_free(ffdotplane);
    vect_free(data);
    exit(0);
}
