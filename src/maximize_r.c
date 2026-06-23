#include "presto.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#define FTOL 3.0e-8

/* Parameters passed to the function being minimized.  Using a   */
/* params struct (instead of static globals) keeps max_r_arr()   */
/* reentrant and matches GSL's gsl_function interface.           */
typedef struct {
    fcomplex *data;
    long numdata;
    int kern_half_width;
} power_call_r_params;


static double power_call_r(double r, void *params)
/*  Maximization function used with an array */
{
    double powargr, powargi;
    fcomplex ans;
    power_call_r_params *p = (power_call_r_params *) params;

    rz_interp(p->data, p->numdata, r, 0.0, p->kern_half_width, &ans);
    return -POWER(ans.r, ans.i);
}


static double brent_min_r(double ax, double bx,
                          power_call_r_params * params, double tol)
/* Find the location of the minimum of power_call_r() in [ax, bx]  */
/* using GSL's Brent minimization routine.                         */
{
    gsl_min_fminimizer *s;
    gsl_function F;
    double a = ax, b = bx, m = 0.5 * (ax + bx);
    int iter = 0, max_iter = 50, status;

    F.function = &power_call_r;
    F.params = params;

    s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

    /* gsl_min_fminimizer_set() requires that the initial guess 'm'   */
    /* have a lower function value than both endpoints.  Temporarily   */
    /* turn off the (aborting) GSL error handler so we can detect a    */
    /* bad bracket and fall back gracefully, as the old fminbr() did.  */
    {
        gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
        status = gsl_min_fminimizer_set(s, &F, m, a, b);
        gsl_set_error_handler(old_handler);
    }

    if (status != GSL_SUCCESS) {
        /* No interior minimum was bracketed:  return the better of    */
        /* the two endpoints (mirrors fminbr()'s fallback behavior).   */
        gsl_min_fminimizer_free(s);
        return (power_call_r(a, params) < power_call_r(b, params)) ? a : b;
    }

    do {
        iter++;
        gsl_min_fminimizer_iterate(s);
        m = gsl_min_fminimizer_x_minimum(s);
        a = gsl_min_fminimizer_x_lower(s);
        b = gsl_min_fminimizer_x_upper(s);
        // printf("%d %f %f %f\n", iter, a, b, m);
        status = gsl_min_test_interval(a, b, tol, 0.0);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free(s);
    return m;
}


double max_r_arr(fcomplex * data, long numdata, double rin,
                 double *rout, rderivs * derivs)
/* Return the Fourier frequency that maximizes the power.  */
{
    double ax, bx, xmin, locpow;
    power_call_r_params params;

    params.data = data;
    params.numdata = numdata;

    /*  Now prep and do the maximization at LOWACC for speed */

    params.kern_half_width = r_resp_halfwidth(LOWACC);
    ax = rin - 0.55;
    bx = rin + 0.55;

    xmin = brent_min_r(ax, bx, &params, FTOL);

    /*  Restart at minimum using HIGHACC to get a better result */
    /*  Note:  Don't know if this is really necessary...        */

    params.kern_half_width = r_resp_halfwidth(HIGHACC);
    ax = xmin - 0.05;
    bx = xmin + 0.05;

    xmin = brent_min_r(ax, bx, &params, FTOL);

    /* The following calculate derivatives at the peak        */
    /*    See characteristics.c for their definitions         */

    *rout = xmin;
    locpow = get_localpower(data, numdata, xmin);
    get_derivs3d(data, numdata, xmin, 0.0, 0.0, locpow, derivs);
    return derivs->pow;
}

#undef FTOL
