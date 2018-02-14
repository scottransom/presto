#include "presto.h"

#define ZSCALE 4.0

#define UNUSED(x) (void)(x)

void amoeba(double p[3][2], double *y, double ftol,
            double (*funk) (double[], fcomplex[], long[], float[], int[], int[]),
            int *nfunk, fcomplex data[], long *numdata, float *locpows, int *numharm, int *kernhw);
    
static double power_call_rz(double rz[], fcomplex data[], long *numdata,
                            float *locpows, int *numharm, int *kernhw)
/* f-fdot plane power function */
{
    double powargr, powargi;
    fcomplex ans;

    UNUSED(locpows);
    UNUSED(numharm);
    rz_interp(data, *numdata, rz[0], rz[1] * ZSCALE, *kernhw, &ans);
    return -POWER(ans.r, ans.i);
}

double max_rz_arr(fcomplex * data, long numdata, double rin, double zin,
                  double *rout, double *zout, rderivs * derivs)
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */
{
    double y[3], x[3][2], step = 0.4;
    float locpow = 0.0;
    int numeval = 0, numharm = 1, max_kernhw;

    /*  Now prep the maximization at LOWACC for speed */

    /* Use a slightly larger working value for 'z' just incase */
    /* the true value of z is a little larger than z.  This    */
    /* keeps a little more accuracy.                           */

    max_kernhw = z_resp_halfwidth(fabs(zin) + 4.0, LOWACC);

    /* Initialize the starting simplex */

    x[0][0] = rin - step;
    x[0][1] = zin / ZSCALE - step;
    x[1][0] = rin - step;
    x[1][1] = zin / ZSCALE + step;
    x[2][0] = rin + step;
    x[2][1] = zin / ZSCALE;

    /* Initialize the starting function values */

    y[0] = power_call_rz(x[0], data, &numdata, &locpow, &numharm, &max_kernhw);
    y[1] = power_call_rz(x[1], data, &numdata, &locpow, &numharm, &max_kernhw);
    y[2] = power_call_rz(x[2], data, &numdata, &locpow, &numharm, &max_kernhw);

    /* Call the solver: */

    numeval = 0;
    amoeba(x, y, 1.0e-7, power_call_rz, &numeval,
           data, &numdata, &locpow, &numharm, &max_kernhw);

    /*  Restart at minimum using HIGHACC to get a better result */

    max_kernhw = z_resp_halfwidth(fabs(x[0][1]) + 4.0, HIGHACC);

    /* Re-Initialize some of the starting simplex */

    x[1][0] = x[0][0] + 0.01;
    x[1][1] = x[0][1];
    x[2][0] = x[0][0];
    x[2][1] = x[0][1] + 0.01;

    /* Re-Initialize the starting function values */

    y[0] = power_call_rz(x[0], data, &numdata, &locpow, &numharm, &max_kernhw);
    y[1] = power_call_rz(x[1], data, &numdata, &locpow, &numharm, &max_kernhw);
    y[2] = power_call_rz(x[2], data, &numdata, &locpow, &numharm, &max_kernhw);

    /* Call the solver: */

    numeval = 0;
    amoeba(x, y, 1.0e-10, power_call_rz, &numeval,
           data, &numdata, &locpow, &numharm, &max_kernhw);

    /* The following calculates derivatives at the peak           */

    *rout = x[0][0];
    *zout = x[0][1] * ZSCALE;
    locpow = get_localpower3d(data, numdata, *rout, *zout, 0.0);
    get_derivs3d(data, numdata, *rout, *zout, 0.0, locpow, derivs);
    return -y[0];
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

static double power_call_rz_harmonics(double rz[], fcomplex data[], long *numdata,
                                      float *locpows, int *numharm, int *kernhw)
{
    int ii;
    double total_power = 0.;
    double powargr, powargi;
    fcomplex ans;

    for (ii = 0; ii < *numharm; ii++) {
        int n = ii + 1;
        rz_interp(data, *numdata, rz[0] * n, rz[1] * ZSCALE * n, *kernhw, &ans);
        total_power += POWER(ans.r, ans.i) / locpows[ii];
    }
    return -total_power;
}

void max_rz_arr_harmonics(fcomplex data[], long numdata,
                          int num_harmonics,
                          double rin, double zin,
                          double *rout, double *zout,
                          rderivs derivs[], double powers[])
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */
{
    double y[3], x[3][2], step = 0.4;
    float *locpow;
    int numeval, ii, max_kernhw;

    locpow = gen_fvect(num_harmonics);

    for (ii = 0; ii < num_harmonics; ii++) {
        int n = ii + 1;
        locpow[ii] = get_localpower3d(data, numdata, rin * n, zin * n, 0.0);
    }

    /*  Now prep the maximization at LOWACC for speed */

    /* Use a slightly larger working value for 'z' just incase */
    /* the true value of z is a little larger than z.  This    */
    /* keeps a little more accuracy.                           */

    max_kernhw = z_resp_halfwidth(fabs(zin * num_harmonics) + 4.0, LOWACC);

    /* Initialize the starting simplex */

    x[0][0] = rin - step;
    x[0][1] = zin / ZSCALE - step;
    x[1][0] = rin - step;
    x[1][1] = zin / ZSCALE + step;
    x[2][0] = rin + step;
    x[2][1] = zin / ZSCALE;

    /* Initialize the starting function values */

    y[0] = power_call_rz_harmonics(x[0], data, &numdata, locpow, &num_harmonics, &max_kernhw);
    y[1] = power_call_rz_harmonics(x[1], data, &numdata, locpow, &num_harmonics, &max_kernhw);
    y[2] = power_call_rz_harmonics(x[2], data, &numdata, locpow, &num_harmonics, &max_kernhw);

    /* Call the solver: */

    numeval = 0;
    amoeba(x, y, 1.0e-7, power_call_rz_harmonics, &numeval,
           data, &numdata, locpow, &num_harmonics, &max_kernhw);

    /*  Restart at minimum using HIGHACC to get a better result */

    max_kernhw = z_resp_halfwidth(fabs(x[0][1] * num_harmonics)
                                  + 4.0, HIGHACC);

    /* Re-Initialize some of the starting simplex */

    x[1][0] = x[0][0] + 0.01;
    x[1][1] = x[0][1];
    x[2][0] = x[0][0];
    x[2][1] = x[0][1] + 0.01;

    /* Re-Initialize the starting function values */

    y[0] = power_call_rz_harmonics(x[0], data, &numdata, locpow, &num_harmonics, &max_kernhw);
    y[1] = power_call_rz_harmonics(x[1], data, &numdata, locpow, &num_harmonics, &max_kernhw);
    y[2] = power_call_rz_harmonics(x[2], data, &numdata, locpow, &num_harmonics, &max_kernhw);

    /* Call the solver: */

    numeval = 0;
    amoeba(x, y, 1.0e-10, power_call_rz_harmonics, &numeval,
           data, &numdata, locpow, &num_harmonics, &max_kernhw);

    /* The following calculates derivatives at the peak           */

    *rout = x[0][0];
    *zout = x[0][1] * ZSCALE;
    for (ii = 0; ii < num_harmonics; ii++) {
        int n = ii + 1;
        locpow[ii] = get_localpower3d(data, numdata, *rout * n, *zout * n, 0.0);
        x[0][0] = *rout * n;
        x[0][1] = *zout / ZSCALE * n;
        powers[ii] = -power_call_rz(x[0], data, &numdata, locpow, &num_harmonics, &max_kernhw);
        get_derivs3d(data, numdata, *rout * n, *zout * n, 0.0, locpow[ii], &(derivs[ii]));
    }
    vect_free(locpow);
}
