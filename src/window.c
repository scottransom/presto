#include "presto.h"
#include "plot2d.h"

#define NPTS   262144

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *file;
    int ii = 0, numbetween, width, nbins, kern_half_width;
    unsigned long num;
    float powargr, powargi;
    float *data, *freqs, *tmp;
    fcomplex *window;
    double amp, freq, ont, offt, tb = 0.0, dt, dtmp;
    double onoffpairs[80], *onoffpt = onoffpairs;
    char filenm[100], psfilenm[100];
    infodata idata;

    if ((argc < 4) || (argc > 5)) {
        printf("\nUsage: window filename numbetween width [ps]\n");
        printf("     This routine displays and saves the window function\n");
        printf("     for an observation.  It takes it's information from\n");
        printf("     the '.inf' file of 'filename' (No suffixes).\n");
        printf("        Numbetween is the level of interpolation to use\n");
        printf("           between each fourier bin.  For example, if\n");
        printf("           numbetween=4, you will get responses at x.0, x.25, \n");
        printf("           x.5, and x.75, where x is some arbitrary frequency.\n");
        printf("        Width is the number of FFT bins to save on each\n");
        printf("           side of the center of the window.  Do not multiply\n");
        printf("           by numbetween (the # of interpolated points).\n");
        printf("        If 'ps' is present, you will get postscript plots.\n\n");
        exit(0);
    }
    printf("\n\n");
    printf("  Window Generation Program\n");
    printf("     by Scott M. Ransom\n");
    printf("        7 Nov, 1997\n\n");

    /* Determine basic parameters */

    data = gen_fvect(NPTS);
    numbetween = atoi(argv[2]);
    if ((numbetween < 1) || (numbetween > 128)) {
        printf("\nNumbetween is out of range.\n\n");
        exit(1);
    }
    width = atoi(argv[3]) * numbetween;
    if ((width < 1) || (width > (NPTS / 2 - 2 * width) / numbetween)) {
        printf("\nWidth is out of range.\n\n");
        exit(1);
    }
    sprintf(filenm, "%s.win%d", argv[1], numbetween);
    sprintf(psfilenm, "%s.win%d.ps", argv[1], numbetween);
    printf("Window data file will be stored as 2*'width'*'numbetween'\n");
    printf("   complex floats in %s, where the '%d'\n", filenm, numbetween);
    printf("   represents 'numbetween'.\n\n");
    readinf(&idata, argv[1]);
    file = chkfopen(filenm, "wb");
    freq = NPTS / 4.0;
    dt = 1.0 / NPTS;
    amp = dt * 2.0;
    num = (long) idata.N;

    /* Send plots to X-windows or a PS file */

    if ((argc == 5) && (!strcmp(argv[4], "ps"))) {
        cpgstart_ps(psfilenm, "landscape");
    } else {
        cpgstart_x("landscape");
    }

    /* Get and normalize the on-off pairs */

    if (idata.numonoff > 1) {
        ii = 0;
        do {
            onoffpairs[ii] = (double) idata.onoff[ii] / (double) (num - 1);
            onoffpairs[ii + 1] = (double) idata.onoff[ii + 1] / (double) (num - 1);
            ii += 2;
        } while (idata.onoff[ii - 1] != idata.N - 1);
        idata.numonoff = ii / 2;
    } else {
        onoffpairs[0] = 0.0;
        onoffpairs[1] = 1.0;
    }

    ont = (*onoffpt++);
    offt = (*onoffpt++);

    /* Main data loop */

    for (ii = 0; ii < NPTS; ii++) {

        tb = ii * dt;

        /*  Advance onoff pointers when signal turns off */

        if (tb >= offt)
            do {
                ont = (*onoffpt++);
                offt = (*onoffpt++);
            }
            while (tb >= offt);

        /*  Signal is on */

        if ((tb >= ont) && (tb < offt)) {
            data[ii] = amp * cos(TWOPI * freq * tb);
        }
        /*  Signal is off  */

        else
            data[ii] = 0.0;

    }

    /* FFT the data */

    realfft(data, NPTS, -1);
    nbins = 2 * width;

    /* Interpolate the area in question */

    window = gen_cvect(nbins);
    tmp = gen_fvect(nbins);
    kern_half_width = r_resp_halfwidth(HIGHACC);
    for (ii = 0; ii < nbins; ii++) {
        dtmp = freq - width / numbetween + (double) (ii) / numbetween;
        rz_interp((fcomplex *) data, NPTS / 2, dtmp, 0.0, kern_half_width,
                  &window[ii]);
        tmp[ii] = POWER(window[ii].r, window[ii].i);
    }

    /* Plot the window's power (normal) */

    freqs = gen_freqs(nbins, (double) (-width / numbetween), 1.0 / numbetween);
    xyline(nbins, freqs, tmp, "Fourier Frequency Offset", "Normalized Power", 1);
    printf("\nThe full resolution window function power.\n");

    /* Plot the window's power (Log_10) */

    for (ii = 0; ii < nbins; ii++) {
        tmp[ii] = log10(tmp[ii]);
    }
    xyline(nbins, freqs, tmp, "Fourier Frequency Offset",
           "Log_10 Normalized Power", 1);
    printf("\nThe full resolution window function power (Log base 10).\n");

    /* Write the window function to file */

    chkfwrite(window, sizeof(float), (unsigned long) (2 * nbins), file);
    fclose(file);

    cpgend();

    printf("\nDone.\n\n");

    vect_free(data);
    vect_free(window);
    vect_free(tmp);
    vect_free(freqs);
    exit(0);
}


#undef NPTS
