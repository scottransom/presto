#include "presto.h"

#define ZSCALE 4.0
#define WSCALE 20.0

static fcomplex *maxdata;
static int nummaxdata, max_kern_half_width, num_funct_calls;

extern void amoeba(double p[3][2], double y[], double ftol,
                   double (*funk) (double[]), int *nfunk);

static double power_call_rzw(double rzw[])
/*  Maximization function used with an array */
{
    double powargr, powargi;
    fcomplex ans;

    num_funct_calls++;
    rzw_interp(maxdata, nummaxdata, rzw[0], rzw[1] * ZSCALE,
               rzw[2] * WSCALE, max_kern_half_width, &ans);
    return -POWER(ans.r, ans.i);
}


double max_rzw_arr(fcomplex * data, int numdata, double rin, double zin,
                   double win, double *rout, double *zout,
                   double *wout, rderivs * derivs)
/* Return the Fourier frequency, f-dot, and fdotdot that    */
/* maximizes the power.                                     */
{
    double maxpower, x[3], locpow;

    maxdata = data;
    nummaxdata = numdata;

    /* Use a slightly larger working value for 'z' just incase */
    /* the true value of z is a little larger than z.  This    */
    /* keeps a little more accuracy.                           */

    max_kern_half_width = w_resp_halfwidth(fabs(zin) + 4.0, win, HIGHACC);
    x[0] = rin;
    x[1] = zin / ZSCALE;
    x[2] = win / WSCALE;

    /* Call the solver: */

    solvopt_options[0] = -0.1;
    maxpower = solvopt(3, x, power_call_rzw, NULL, solvopt_options, NULL, NULL);

    /* Re-run solvopt from the value obtained if needed. */

    if (solvopt_options[8] == -11.0) {
        solvopt_options[0] = -0.01;
        maxpower = solvopt(3, x, power_call_rzw, NULL, solvopt_options, NULL, NULL);
    }

    printf("\nCalled rzw_interp() %d times.\n", num_funct_calls);

    /* The following calculates derivatives at the peak           */

    x[1] *= ZSCALE;
    x[2] *= WSCALE;
    *rout = x[0];
    *zout = x[1];
    *wout = x[2];
    locpow = get_localpower3d(data, numdata, x[0], x[1], x[2]);
    get_derivs3d(data, numdata, x[0], x[1], x[2], locpow, derivs);
    return -maxpower;
}

double max_rzw_file(FILE * fftfile, double rin, double zin, double win,
                    double *rout, double *zout, double *wout, rderivs * derivs)
/* Return the Fourier frequency, f-dot, and fdotdot that    */
/* maximizes the power of the candidate in 'fftfile'.       */
{
    double maxpow, rin_int, rin_frac;
    int kern_half_width, filedatalen, startbin, extra = 10;
    fcomplex *filedata;

    rin_frac = modf(rin, &rin_int);
    kern_half_width = w_resp_halfwidth(fabs(zin) + 4.0, win, HIGHACC);
    filedatalen = 2 * kern_half_width + extra;
    startbin = (int) rin_int - filedatalen / 2;

    filedata = read_fcomplex_file(fftfile, startbin, filedatalen);
    maxpow = max_rzw_arr(filedata, filedatalen, rin_frac + filedatalen / 2,
                         zin, win, rout, zout, wout, derivs);
    *rout += startbin;
    vect_free(filedata);
    return maxpow;
}
