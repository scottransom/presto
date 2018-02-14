#include "presto.h"
#define FTOL 3.0e-8

/* Function declarations */

static double power_call_r(double r);
double fminbr(double a, double b, double (*f) (double x), double tol);


/* Some Static-Global Variables */
static fcomplex *maxdata;
static int max_kern_half_width;
static long nummaxdata;

static double power_call_r(double r)
/*  Maximization function used with an array */
{
    double powargr, powargi;
    fcomplex ans;

    rz_interp(maxdata, nummaxdata, r, 0.0, max_kern_half_width, &ans);
    powargr = (double) ans.r;
    powargi = (double) ans.i;
    return -POWER(powargr, powargi);
}


double max_r_arr(fcomplex * data, long numdata, double rin,
                 double *rout, rderivs * derivs)
/* Return the Fourier frequency that maximizes the power.  */
{
    double ax, bx, xmin, locpow;

    maxdata = data;
    nummaxdata = numdata;

    /*  Now prep and do the maximization at LOWACC for speed */

    max_kern_half_width = r_resp_halfwidth(LOWACC);
    ax = rin - 0.55;
    bx = rin + 0.55;

    xmin = fminbr(ax, bx, power_call_r, FTOL);

    /*  Restart at minimum using HIGHACC to get a better result */
    /*  Note:  Don't know if this is really necessary...        */

    max_kern_half_width = r_resp_halfwidth(HIGHACC);
    ax = xmin - 0.05;
    bx = xmin + 0.05;

    xmin = fminbr(ax, bx, power_call_r, FTOL);

    /* The following calculate derivatives at the peak        */
    /*    See characteristics.c for their definitions         */

    *rout = xmin;
    locpow = get_localpower(data, numdata, xmin);
    get_derivs3d(data, numdata, xmin, 0.0, 0.0, locpow, derivs);
    return derivs->pow;
}

#undef FTOL
