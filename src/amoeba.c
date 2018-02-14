/* Note:  This version of the Downhill Simplex Algorithm is     */
/*        a modified version on the NR in C 2ed version         */
/*        by Press et el.  Their copyright statement is below.  */
/*        Modification by Scott M. Ransom, 2 March 2001         */
/* (C) Copr. 1986-92 Numerical Recipes Software 3#1y-i.31-.     */

#include <math.h>
#include <stdio.h>
#include "ransomfft.h"

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

static double amotry(double p[3][2], double *y, double *psum,
                     double (*funk) (double[], fcomplex[], long[], float[], int[], int[]),
                     int ihi, double fac, fcomplex data[], long *numdata,
                     float *locpows, int *numharm, int *kernhw);

void amoeba(double p[3][2], double *y, double ftol,
            double (*funk) (double[], fcomplex[], long[], float[], int[], int[]),
            int *nfunk, fcomplex data[], long *numdata, float *locpows, int *numharm, int *kernhw)
{
    int ii, ihi, ilo, inhi;
    double rtol, ysave, ytry, psum[2], tempzz;

    *nfunk = 0;
    psum[0] = p[0][0] + p[1][0] + p[2][0];
    psum[1] = p[0][1] + p[1][1] + p[2][1];
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
        for (ii = 0; ii <= 2; ii++) {
            if (y[ii] <= y[ilo])
                ilo = ii;
            if (y[ii] > y[ihi]) {
                inhi = ihi;
                ihi = ii;
            } else if (y[ii] > y[inhi] && ii != ihi)
                inhi = ii;
        }
        rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + 1.0e-15);
        if (rtol < ftol) {
            SWAP(y[0], y[ilo])
                SWAP(p[0][0], p[ilo][0])
                SWAP(p[0][1], p[ilo][1])
                break;
        }
        if (*nfunk >= 5000) {
            /*
               printf("\nWarning:  amoeba() exceeded %d iterations for r=%f  z=%f.\n",
               *nfunk, p[0][0], p[0][1]);
             */
            return;
        }
        *nfunk += 2;
        ytry = amotry(p, y, psum, funk, ihi, -1.0, data,
                      numdata, locpows, numharm, kernhw);
        if (ytry <= y[ilo])
            ytry = amotry(p, y, psum, funk, ihi, 2.0, data,
                          numdata, locpows, numharm, kernhw);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry = amotry(p, y, psum, funk, ihi, 0.5, data,
                          numdata, locpows, numharm, kernhw);
            if (ytry >= ysave) {
                for (ii = 0; ii <= 2; ii++) {
                    if (ii != ilo) {
                        p[ii][0] = psum[0] = 0.5 * (p[ii][0] + p[ilo][0]);
                        p[ii][1] = psum[1] = 0.5 * (p[ii][1] + p[ilo][1]);
                        y[ii] = (*funk) (psum, data, numdata, locpows,
                                         numharm, kernhw);
                    }
                }
                *nfunk += 2;
                psum[0] = p[0][0] + p[1][0] + p[2][0];
                psum[1] = p[0][1] + p[1][1] + p[2][1];
            }
        } else
            --(*nfunk);
    }
}


static double amotry(double p[3][2], double *y, double *psum,
                     double (*funk) (double[], fcomplex[], long[], float[], int[], int[]),
                     int ihi, double fac, fcomplex data[], long *numdata, float *locpows,
                     int *numharm, int *kernhw)
{
    double fac1, fac2, ytry, ptry[2];

    fac1 = 0.5 * (1.0 - fac);
    fac2 = fac1 - fac;
    ptry[0] = psum[0] * fac1 - p[ihi][0] * fac2;
    ptry[1] = psum[1] * fac1 - p[ihi][1] * fac2;
    ytry = (*funk) (ptry, data, numdata, locpows, numharm, kernhw);
    if (ytry < y[ihi]) {
        y[ihi] = ytry;
        psum[0] += ptry[0] - p[ihi][0];
        p[ihi][0] = ptry[0];
        psum[1] += ptry[1] - p[ihi][1];
        p[ihi][1] = ptry[1];
    }
    return ytry;
}
