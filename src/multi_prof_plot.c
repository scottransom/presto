#include "plot2d.h"

void multi_prof_plot(int proflen, int numprofs, double *profiles,
                     double *sumprof, const char *xlbl,
                     double loly, double ldy, const char *lylbl,
                     double lory, double rdy, const char *rylbl)
{
    float *x, *y, yoffset, ynorm;
    float lox = 0.0, hix = 1.0, ddx = 0.01;
    double ymin = 0.0, ymax = 0.0;
    int i, j, index;

    x = gen_fvect(proflen + 1);
    y = gen_fvect(proflen + 1);
    for (i = 0; i <= proflen; i++)
        x[i] = (float) i / (float) proflen;

    /* The multiplots... */

    /* Define the Viewport */
    cpgsvp(0.1, 0.85, 0.2, 0.9);
    /* Define the Window */
    cpgswin(lox - ddx, hix + ddx, (float) (loly - ldy),
            (float) (loly + numprofs * ldy));
    /* Define the left border */
    cpgbox("CST", 0.2, 2, "BNST", 0.0, 0);
    /* Write the left-hand label */
    cpgmtxt("L", 2.6, 0.5, 0.5, lylbl);
    /* Re-Define the Window */
    cpgswin(lox - ddx, hix + ddx, (float) (lory - rdy),
            (float) (lory + numprofs * rdy));
    /* Define the right border */
    cpgbox("", 0.2, 2, "CMST", 0.0, 0);
    /* Write the right-hand label */
    cpgmtxt("R", 3.0, 0.5, 0.5, rylbl);

    /* Plot the individual channel profiles */

    for (i = 0; i < numprofs; i++) {
        /* Determine min and max values to plot */
        index = i * proflen;
        yoffset = lory + i * rdy;
        dfind_min_max_arr(proflen, profiles + index, &ymin, &ymax);
        ynorm = 0.9 * rdy / (ymax - ymin);
        for (j = 0; j < proflen; j++)
            y[j] = (profiles[index + j] - ymin) * ynorm + yoffset;
        y[proflen] = y[0];
        cpgbin(proflen + 1, x, y, 1);
    }

    /* The summed plot... */

    /* Define the Viewport */
    cpgsvp(0.1, 0.85, 0.1, 0.2);
    /* Define the Window */
    cpgswin(lox - ddx, hix + ddx, -0.1, 1.0);
    /* Define the border */
    cpgbox("BNST", 0.2, 2, "BC", 0.0, 0);
    /* Write the bottom label */
    cpgmtxt("B", 3.0, 0.5, 0.5, xlbl);
    /* Determine min and max values to plot */
    dfind_min_max_arr(proflen, sumprof, &ymin, &ymax);
    ynorm = 0.9 / (ymax - ymin);
    for (j = 0; j < proflen; j++)
        y[j] = (sumprof[j] - ymin) * ynorm;
    y[proflen] = y[0];
    /* Plot the summed profile */
    cpgbin(proflen + 1, x, y, 1);

    /* Cleanup */
    vect_free(x);
    vect_free(y);

}
