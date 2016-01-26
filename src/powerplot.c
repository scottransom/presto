#include "plot2d.h"

void powerplot(int npts, float *freqs, float *amp, float norm, int id)
{
    float *pwr, xmin, xmax, ymin, ymax;
    float overy, over = 0.1;
    int i, ptr;

    pwr = (float *) malloc((size_t) npts * sizeof(float));
    if (!pwr) {
        printf("Error allocating 'pwr' in powerplot.  Exiting.\n\n");
        exit(EXIT_FAILURE);
    }
    /* Turn the complex amps into a power series: */
    for (i = 0; i < npts; i++) {
        ptr = i * 2;
        pwr[i] = plot_power(amp[ptr], amp[ptr + 1]) / norm;
    }

    /* Determine min and max values to plot and scaling: */
    find_min_max_arr(npts, freqs, &xmin, &xmax);
    find_min_max_arr(npts, pwr, &ymin, &ymax);
    overy = over * (ymax - ymin);
    ymax += overy;

    /* Setup the plot screen: */
    cpgenv(xmin, xmax, ymin, ymax, 0, 1);

    /* Choose the font: */
    cpgscf(2);

    /* Label the axes: */
    cpglab("Frequency", "Power", "");

    /* Add ID line if required */
    if (id == 1)
        cpgiden();

    /* Plot the points: */
    cpgline(npts, freqs, pwr);

    free(pwr);
}

double plot_power(double rl, double im)
{
    return rl * rl + im * im;
}
