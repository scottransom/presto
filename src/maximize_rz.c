#include "presto.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#define ZSCALE 4.0

/* Parameters passed to the power functions being maximized.  Using a params  */
/* struct (instead of the old amoeba() function-pointer plumbing) keeps these  */
/* routines reentrant and matches GSL's gsl_multimin_function interface, as    */
/* already done for the 1-D Brent minimization in maximize_r.c.                */
typedef struct {
    fcomplex *data;
    long numdata;
    float *locpows;             /* per-harmonic local powers (harmonics only) */
    int numharm;
    int kernhw;
} rz_params;


static double power_call_rz(double r, double z_over_scale, rz_params * p)
/* Single-harmonic f-fdot plane power (negated, since we minimize).  Internal  */
/* coordinates are (r, z/ZSCALE), matching the old amoeba simplex.             */
{
    double powargr, powargi;
    fcomplex ans;

    rz_interp(p->data, p->numdata, r, z_over_scale * ZSCALE, p->kernhw, &ans);
    return -POWER(ans.r, ans.i);
}


static double power_call_rz_gsl(const gsl_vector * v, void *params)
{
    rz_params *p = (rz_params *) params;
    return power_call_rz(gsl_vector_get(v, 0), gsl_vector_get(v, 1), p);
}


static double power_call_rz_harmonics_gsl(const gsl_vector * v, void *params)
{
    rz_params *p = (rz_params *) params;
    double r = gsl_vector_get(v, 0), z_over_scale = gsl_vector_get(v, 1);
    double total_power = 0.0, powargr, powargi;
    fcomplex ans;
    int ii;

    for (ii = 0; ii < p->numharm; ii++) {
        int n = ii + 1;
        rz_interp(p->data, p->numdata, r * n, z_over_scale * ZSCALE * n,
                  p->kernhw, &ans);
        total_power += POWER(ans.r, ans.i) / p->locpows[ii];
    }
    /* Guard the simplex: GSL's Nelder-Mead aborts on a non-finite function   */
    /* value, whereas the old amoeba() just carried the bad value along.  A   */
    /* harmonic pushed out of the data band gives 0 power over a 0 local      */
    /* power (0/0 -> NaN); return a large finite penalty so the simplex is    */
    /* simply steered away from there.                                        */
    if (!isfinite(total_power))
        return 0.0;
    return -total_power;
}


static double simplex_max_rz(gsl_multimin_function * F, double *r,
                             double *z_over_scale, double step, double tol)
/* Maximize the power (i.e. minimize F) with GSL's Nelder-Mead simplex,        */
/* starting from (*r, *z_over_scale) with an initial simplex of size 'step'.   */
/* On return, (*r, *z_over_scale) hold the optimum and the minimum function    */
/* value is returned.  Replaces the old Numerical Recipes amoeba().            */
{
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector *ss = gsl_vector_alloc(2);
    int iter = 0, status;
    double size, fval;

    gsl_vector_set(x, 0, *r);
    gsl_vector_set(x, 1, *z_over_scale);
    gsl_vector_set_all(ss, step);
    gsl_multimin_fminimizer_set(s, F, x, ss);

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
            break;
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tol);
    } while (status == GSL_CONTINUE && iter < 5000);

    *r = gsl_vector_get(s->x, 0);
    *z_over_scale = gsl_vector_get(s->x, 1);
    fval = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    return fval;
}


double max_rz_arr(fcomplex * data, long numdata, double rin, double zin,
                  double *rout, double *zout, rderivs * derivs)
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */
{
    double r, z, fval, locpow;
    rz_params par;
    gsl_multimin_function F;

    par.data = data;
    par.numdata = numdata;
    par.locpows = NULL;
    par.numharm = 1;
    F.n = 2;
    F.f = &power_call_rz_gsl;
    F.params = &par;

    /*  Now prep the maximization at LOWACC for speed */

    /* Use a slightly larger working value for 'z' just incase */
    /* the true value of z is a little larger than z.  This    */
    /* keeps a little more accuracy.                           */

    par.kernhw = z_resp_halfwidth(fabs(zin) + 4.0, LOWACC);
    r = rin;
    z = zin / ZSCALE;
    fval = simplex_max_rz(&F, &r, &z, 0.4, 1.0e-4);

    /*  Restart at minimum using HIGHACC to get a better result */

    /* NOTE:  'z' here is the *internal* (z/ZSCALE) coordinate, so it must be */
    /* scaled back to a physical Fourier f-dot before sizing the kernel.      */
    par.kernhw = z_resp_halfwidth(fabs(z * ZSCALE) + 4.0, HIGHACC);
    fval = simplex_max_rz(&F, &r, &z, 0.01, 1.0e-7);

    /* The following calculates derivatives at the peak           */

    *rout = r;
    *zout = z * ZSCALE;
    locpow = get_localpower3d(data, numdata, *rout, *zout, 0.0);
    get_derivs3d(data, numdata, *rout, *zout, 0.0, locpow, derivs);
    return -fval;
}


double max_rz_file(FILE * fftfile, double rin, double zin,
                   double *rout, double *zout, rderivs * derivs)
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power of the candidate in 'fftfile'.       */
{
    double maxz, maxpow, rin_int, rin_frac;
    int kern_half_width, filedatalen, extra = 10;
    long startbin;
    fcomplex *filedata;

    maxz = fabs(zin) + 4.0;
    rin_frac = modf(rin, &rin_int);
    kern_half_width = z_resp_halfwidth(maxz, HIGHACC);
    filedatalen = 2 * kern_half_width + extra;
    startbin = (long) rin_int - filedatalen / 2;

    filedata = read_fcomplex_file(fftfile, startbin, filedatalen);
    maxpow = max_rz_arr(filedata, filedatalen, rin_frac + filedatalen / 2,
                        zin, rout, zout, derivs);
    *rout += startbin;
    vect_free(filedata);
    return maxpow;
}


void max_rz_arr_harmonics(fcomplex data[], long numdata,
                          int num_harmonics,
                          double rin, double zin,
                          double *rout, double *zout,
                          rderivs derivs[], double powers[])
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */
{
    double r, z;
    float *locpow;
    rz_params par;
    gsl_multimin_function F;
    int ii;

    locpow = gen_fvect(num_harmonics);

    for (ii = 0; ii < num_harmonics; ii++) {
        int n = ii + 1;
        locpow[ii] = get_localpower3d(data, numdata, rin * n, zin * n, 0.0);
    }

    par.data = data;
    par.numdata = numdata;
    par.locpows = locpow;
    par.numharm = num_harmonics;
    F.n = 2;
    F.f = &power_call_rz_harmonics_gsl;
    F.params = &par;

    /*  Now prep the maximization at LOWACC for speed */

    /* Use a slightly larger working value for 'z' just incase */
    /* the true value of z is a little larger than z.  This    */
    /* keeps a little more accuracy.                           */

    par.kernhw = z_resp_halfwidth(fabs(zin * num_harmonics) + 4.0, LOWACC);
    r = rin;
    z = zin / ZSCALE;
    simplex_max_rz(&F, &r, &z, 0.4, 1.0e-4);

    /*  Restart at minimum using HIGHACC to get a better result */

    /* As above, 'z' is the internal (z/ZSCALE) coordinate at this point. */
    par.kernhw = z_resp_halfwidth(fabs(z * ZSCALE * num_harmonics) + 4.0, HIGHACC);
    simplex_max_rz(&F, &r, &z, 0.01, 1.0e-7);

    /* The following calculates derivatives at the peak           */

    *rout = r;
    *zout = z * ZSCALE;
    for (ii = 0; ii < num_harmonics; ii++) {
        int n = ii + 1;
        locpow[ii] = get_localpower3d(data, numdata, *rout * n, *zout * n, 0.0);
        powers[ii] = -power_call_rz(*rout * n, *zout / ZSCALE * n, &par);
        get_derivs3d(data, numdata, *rout * n, *zout * n, 0.0, locpow[ii],
                     &(derivs[ii]));
    }
    vect_free(locpow);
}
