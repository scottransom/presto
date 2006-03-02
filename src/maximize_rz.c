#include "presto.h"

#define ZSCALE 4.0

static fcomplex *maxdata;
static int nummaxdata, max_kern_half_width;

extern void amoeba(double p[3][2], double y[], double ftol,
                   double (*funk) (double[]), int *nfunk);

static double power_call_rz(double rz[])
/* f-fdot plane power function */
{
   double powargr, powargi;
   fcomplex ans;

   /* num_funct_calls++; */
   rz_interp(maxdata, nummaxdata, rz[0], rz[1] * ZSCALE, max_kern_half_width, &ans);
   return -POWER(ans.r, ans.i);
}


double max_rz_arr(fcomplex * data, int numdata, double rin, double zin,
                  double *rout, double *zout, rderivs * derivs)
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */
{
   double y[3], x[3][2], step = 0.4;
   float locpow;
   int numeval;

   maxdata = data;
   nummaxdata = numdata;
   locpow = get_localpower3d(data, numdata, rin, zin, 0.0);

   /*  Now prep the maximization at LOWACC for speed */

   /* Use a slightly larger working value for 'z' just incase */
   /* the true value of z is a little larger than z.  This    */
   /* keeps a little more accuracy.                           */

   max_kern_half_width = z_resp_halfwidth(fabs(zin) + 4.0, LOWACC);

   /* Initialize the starting simplex */

   x[0][0] = rin - step;
   x[0][1] = zin / ZSCALE - step;
   x[1][0] = rin - step;
   x[1][1] = zin / ZSCALE + step;
   x[2][0] = rin + step;
   x[2][1] = zin / ZSCALE;

   /* Initialize the starting function values */

   y[0] = power_call_rz(x[0]);
   y[1] = power_call_rz(x[1]);
   y[2] = power_call_rz(x[2]);

   /* Call the solver: */

   numeval = 0;
   amoeba(x, y, 1.0e-7, power_call_rz, &numeval);

   /*  Restart at minimum using HIGHACC to get a better result */

   max_kern_half_width = z_resp_halfwidth(fabs(x[0][1]) + 4.0, HIGHACC);

   /* Re-Initialize some of the starting simplex */

   x[1][0] = x[0][0] + 0.01;
   x[1][1] = x[0][1];
   x[2][0] = x[0][0];
   x[2][1] = x[0][1] + 0.01;

   /* Re-Initialize the starting function values */

   y[0] = power_call_rz(x[0]);
   y[1] = power_call_rz(x[1]);
   y[2] = power_call_rz(x[2]);

   /* Call the solver: */

   numeval = 0;
   amoeba(x, y, 1.0e-10, power_call_rz, &numeval);

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
   int kern_half_width, filedatalen, startbin, extra = 10;
   fcomplex *filedata;

   maxz = fabs(zin) + 4.0;
   rin_frac = modf(rin, &rin_int);
   kern_half_width = z_resp_halfwidth(maxz, HIGHACC);
   filedatalen = 2 * kern_half_width + extra;
   startbin = (int) rin_int - filedatalen / 2;

   filedata = read_fcomplex_file(fftfile, startbin, filedatalen);
   maxpow = max_rz_arr(filedata, filedatalen, rin_frac + filedatalen / 2,
                       zin, rout, zout, derivs);
   *rout += startbin;
   free(filedata);
   return maxpow;
}
