#include "plot2d.h"
#include <string.h>

int nice_output_2(char *output, double val, double err, int len);

void cpgstart_ps(const char *filenm, const char *orientation)
{
    char tmp[100];
    if (0 == strcmp(orientation, "portrait")) {
        sprintf(tmp, "%s/VCPS", filenm);
        if (cpgopen(tmp) <= 0)
            exit(EXIT_FAILURE);
    } else {
        sprintf(tmp, "%s/CPS", filenm);
        if (cpgopen(tmp) <= 0)
            exit(EXIT_FAILURE);
    }
}

void cpgstart_x(const char *orientation)
{
    if (cpgopen("/XWIN") <= 0)
        exit(EXIT_FAILURE);
    if (0 == strcmp(orientation, "portrait")) {
        cpgpap(8.5, 11.0 / 8.5);
    } else {
        cpgpap(11.0, 8.5 / 11.0);
    }
    cpgask(1);
}

void find_min_max_arr(int npts, float *arr, float *min, float *max)
{
    int i;

    *min = 1.1E38;
    *max = -1.1E38;
    for (i = 0; i < npts; i++) {
        if (arr[i] > *max)
            *max = arr[i];
        if (arr[i] < *min)
            *min = arr[i];
    }
}

void dfind_min_max_arr(int npts, double *arr, double *min, double *max)
{
    int i;

    *min = 1.1E100;
    *max = -1.1E100;
    for (i = 0; i < npts; i++) {
        if (arr[i] > *max)
            *max = arr[i];
        if (arr[i] < *min)
            *min = arr[i];
    }
}

void dxyline(int npts, double *x, double *y, const char *xlab,
             const char *ylab, int id)
/* Wrapper to plot double precision vectors */
{
    float *fx, *fy;
    long i;

    /* Copy our double vectors to float vectors */

    fx = (float *) malloc(sizeof(float) * npts);
    fy = (float *) malloc(sizeof(float) * npts);

    for (i = 0; i < npts; i++) {
        fx[i] = (float) x[i];
        fy[i] = (float) y[i];
    }

    /* Call xyline */

    xyline(npts, fx, fy, xlab, ylab, id);

    /* Free our memory */

    free(fx);
    free(fy);
}


void xyline(int npts, float *x, float *y, const char *xlab, const char *ylab, int id)
{
    float xmin, xmax, ymin, ymax;
    float overy, over = 0.1;

    /* Determine min and max values to plot and scaling: */
    find_min_max_arr(npts, x, &xmin, &xmax);
    find_min_max_arr(npts, y, &ymin, &ymax);
    overy = over * (ymax - ymin);
    ymax += overy;
    ymin -= overy;

    /* Setup the plot screen: */
    cpgenv(xmin, xmax, ymin, ymax, 0, 0);

    /* Choose the font: */
    cpgscf(2);

    /* Label the axes: */
    cpglab(xlab, ylab, "");

    /* Add ID line if required */
    if (id == 1)
        cpgiden();

    /* Plot the points: */
    cpgline(npts, x, y);

}


void dxybinned(int npts, double *x, double *y, const char *xlab,
               const char *ylab, int id)
/* Wrapper to plot double precision vectors */
{
    float *fx, *fy;
    long i;

    /* Copy our double vectors to float vectors */

    fx = (float *) malloc(sizeof(float) * npts);
    fy = (float *) malloc(sizeof(float) * npts);

    for (i = 0; i < npts; i++) {
        fx[i] = (float) x[i];
        fy[i] = (float) y[i];
    }

    /* Call xyline */

    xybinned(npts, fx, fy, xlab, ylab, id);

    /* Free our memory */

    free(fx);
    free(fy);
}


void xybinned(int npts, float *x, float *y, const char *xlab,
              const char *ylab, int id)
{
    float xmin, xmax, ymin, ymax;
    float overy, over = 0.1;

    /* Determine min and max values to plot and scaling: */
    find_min_max_arr(npts, x, &xmin, &xmax);
    find_min_max_arr(npts, y, &ymin, &ymax);
    overy = over * (ymax - ymin);
    ymax += overy;
    ymin -= overy;

    /* Setup the plot screen: */
    cpgenv(xmin, xmax, ymin, ymax, 0, 0);

    /* Choose the font: */
    cpgscf(2);

    /* Label the axes: */
    cpglab(xlab, ylab, "");

    /* Add ID line if required */
    if (id == 1)
        cpgiden();

    /* Plot the points: */
    cpgbin(npts, x, y, 0);
}


void xyline2lab(int npts, float *x, float *y, float *y2, const char *xlab,
                const char *ylab, const char *ylab2, int id)
{
    float xmin, xmax, ymin, ymax, ymin2, ymax2;
    float overy, over = 0.1;

    /* Determine min and max values to plot and scaling: */
    find_min_max_arr(npts, x, &xmin, &xmax);
    find_min_max_arr(npts, y, &ymin, &ymax);
    find_min_max_arr(npts, y2, &ymin2, &ymax2);
    overy = over * (ymax - ymin);
    ymax += overy;
    ymin -= overy;
    overy = over * (ymax2 - ymin2);
    ymax2 += overy;
    ymin2 -= overy;

    /* Choose the font: */
    cpgscf(2);

    /* Setup the plot screen for the first set of y's: */
    cpgpage();
    cpgvstd();
    cpgswin(xmin, xmax, ymin, ymax);
    cpgbox("BCNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 3.0, 0.5, 0.5, xlab);
    cpgmtxt("L", 2.6, 0.5, 0.5, ylab);

    /* Plot the points for the 1st y axis: */
    cpgline(npts, x, y);

    /* Setup the plot screen for the second set of y's: */
    cpgvstd();
    cpgswin(xmin, xmax, ymin2, ymax2);
    cpgbox("", 0.0, 0, "CMST", 0.0, 0);
    cpgmtxt("R", 3.0, 0.5, 0.5, ylab2);

    /* Plot the points for the 2nd y axis: */
    cpgline(npts, x, y2);

    /* Add ID line if required */
    if (id == 1)
        cpgiden();

}


void plot_spectrum(fcomplex * spect, int numspect,
                   double lor, double dr, double T, double average)
/* Plot a chunk of the Fourier power spectrum normalized by average  */
/* The viewing area is left defined with the xvals as _bins_.        */
{
    int ii;
    float *yvals, *xvals, maxy = 0.0;
    double offsetr, hir;
    char lab[100];

    if (lor > 10000) {
        offsetr = floor(lor / 1000.0) * 1000.0;
        sprintf(lab, "Fourier Frequency - %.0f (bins)", offsetr);
    } else {
        offsetr = 0.0;
        sprintf(lab, "Fourier Frequency (bins)");
    }
    lor = lor - offsetr;
    hir = lor + numspect * dr;
    xvals = (float *) malloc(sizeof(float) * numspect);
    yvals = (float *) malloc(sizeof(float) * numspect);
    for (ii = 0; ii < numspect; ii++) {
        xvals[ii] = lor + ii * dr;
        yvals[ii] =
            (spect[ii].r * spect[ii].r + spect[ii].i * spect[ii].i) / average;
        if (yvals[ii] > maxy)
            maxy = yvals[ii];
    }
    maxy *= 1.1;

    /* Setup the plot screen for the first set of y's: */
    cpgpage();
    cpgvstd();

    /* Draw the box for the frequency (Hz) axis */
    cpgswin((lor + offsetr) / T, (hir + offsetr) / T, 0.0, maxy);
    cpgbox("BNST", 0.0, 0, "BCNST", 0.0, 0);
    cpgmtxt("B", 2.5, 0.5, 0.5, "Frequency (Hz)");

    /* Draw the box for the Fourier frequency (bins) axis */
    cpgswin(lor, hir, 0.0, maxy);
    cpgbox("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 2.0, 0.5, 0.5, lab);
    cpgmtxt("L", 2.0, 0.5, 0.5, "Normalized Power");

    /* Plot the points */
    cpgline(numspect, xvals, yvals);

    free(xvals);
    free(yvals);
}


void plot_profile(int proflen, float *profile, const char *title,
                  const char *probtxt, const char *foldtxt,
                  int showerr, float *errors, int showid)
{
    int ii;
    float *x, overy, ymin, ymax;
    float errmin = 0.0, errmax = 0.0, offset, avg = 0.0, av[2];

    find_min_max_arr(proflen, profile, &ymin, &ymax);
    if (showerr)
        find_min_max_arr(proflen, errors, &errmin, &errmax);
    overy = 0.1 * (ymax + errmax - ymin - errmin);
    ymax = ymax + overy + errmax;
    ymin = ymin - overy - errmin;
    x = gen_fvect(proflen);
    for (ii = 0; ii < proflen; ii++)
        x[ii] = (float) ii / (float) proflen;
    cpgenv(0.0, 1.00001, ymin, ymax, 0, 0);
    cpgscf(2);
    cpglab("Pulse Phase", "Counts", "");
    if (showid)
        cpgiden();
    cpgslw(5);
    if (showerr) {
        cpgbin(proflen, x, profile, 0);
    } else {
        cpgline(proflen, x, profile);
    }
    cpgslw(1);
    if (showerr) {
        offset = 0.5 / (float) proflen;
        for (ii = 0; ii < proflen; ii++)
            x[ii] += offset;
        cpgerrb(6, proflen, x, profile, errors, 2);
        cpgpt(proflen, x, profile, 5);
    }
    for (ii = 0; ii < proflen; ii++)
        avg += profile[ii];
    avg /= proflen;
    cpgsls(4);
    x[0] = 0.0;
    x[1] = 1.0;
    av[0] = avg;
    av[1] = avg;
    cpgline(2, x, av);
    cpgsls(1);
    cpgsch(1.3);
    cpgmtxt("T", +2.0, 0.5, 0.5, title);
    cpgsch(1.0);
    cpgmtxt("T", +0.8, 0.5, 0.5, foldtxt);
    cpgmtxt("T", -1.5, 0.5, 0.5, probtxt);
    vect_free(x);
}


int cpgnice_output_2(char *out, double val, double err, int len)
/* Return a string that has the same formatting as       */
/* nice_output_2(), but for PGPLOT.  This way, exponents */
/* are actually in superscript!  Woo-hoo!                */
{
    int inlen, goodlen;
    char tempout[100], part[20], expon[20], *expptr;

    inlen = nice_output_2(tempout, val, err, len);

    /* See if the formatted output */

    expptr = strrchr(tempout, '^');

    /* Not in scientific notation */

    if (expptr == NULL) {
        strcpy(out, tempout);
        return inlen;
    } else {
        goodlen = expptr - tempout;
        strcpy(expon, tempout + goodlen + 1);
        strncpy(part, tempout, goodlen);
        part[goodlen] = 0;
        sprintf(out, "%s\\u%s\\d", part, expon);
        return strlen(out);
    }
}
