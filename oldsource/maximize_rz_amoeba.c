#include "presto.h"

/* Some Static-Global Variables */

static fcomplex *maxdata;
static int nummaxdata, max_kern_half_width;

/* Function definition for the minimization routine */

extern void amoeba(double **p, double y[], int ndim, double ftol,
		   double (*funk)(double []), int *nfunk);

static double power_call_rz(double rz[])
/*  Maximization function used with an array */
{
  double powargr, powargi;
  fcomplex ans;

  rz_interp(maxdata, nummaxdata, rz[0], rz[1], \
	    max_kern_half_width, &ans);
  powargr = (double) ans.r;
  powargi = (double) ans.i;
  return -POWER(powargr, powargi);
}


double max_rz_arr(fcomplex *data, int numdata, double rin, double zin, \
		  double *rout, double *zout, rderivs * derivs)
/* Return the Fourier frequency and Fourier f-dot that      */ 
/* maximizes the power.                                     */
{
  double y[3], **x;
  float locpow;
  int numeval;

  maxdata = data;
  nummaxdata = numdata;
  x = gen_dmatrix(3, 2);

  /*  Now prep the maximization at LOWACC for speed */

  /* Use a slightly larger working value for 'z' just incase */
  /* the true value of z is a little larger than z.  This    */
  /* keeps a little more accuracy.                           */

  max_kern_half_width = z_resp_halfwidth(fabs(zin) + 4.0, LOWACC);

  /* Initialize the starting simplex */

  x[0][0] = rin;
  x[0][1] = zin;
  x[1][0] = rin + 0.1;
  x[1][1] = zin;
  x[2][0] = rin;
  x[2][1] = zin + 0.4;

printf("r = %f  z = %f numeval = %d\n", x[0][0], x[0][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[1][0], x[1][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[2][0], x[2][1], numeval);

  /* Initialize the starting function values */

  y[0] = power_call_rz(x[0]);
  y[1] = power_call_rz(x[1]);
  y[2] = power_call_rz(x[2]);

  /* Call the solver: */

  amoeba(x, y, 2, 1.0e-10, power_call_rz, &numeval);

printf("r = %f  z = %f numeval = %d\n", x[0][0], x[0][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[1][0], x[1][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[2][0], x[2][1], numeval);

  /*  Restart at minimum using HIGHACC to get a better result */

  max_kern_half_width = z_resp_halfwidth(fabs(x[0][1]) + 4.0, HIGHACC);

  /* Re-Initialize some of the starting simplex */

  x[1][0] = x[0][0] + 0.01;
  x[1][1] = x[0][1];
  x[2][0] = x[0][0];
  x[2][1] = x[0][1] + 0.04;

  /* Re-Initialize the starting function values */

  y[0] = power_call_rz(x[0]);
  y[1] = power_call_rz(x[1]);
  y[2] = power_call_rz(x[2]);

  /* Call the solver: */

  amoeba(x, y, 2, 1.0e-10, power_call_rz, &numeval);

printf("r = %f  z = %f numeval = %d\n", x[0][0], x[0][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[1][0], x[1][1], numeval);
printf("r = %f  z = %f numeval = %d\n", x[2][0], x[2][1], numeval);

  /* The following calculates derivatives at the peak           */

  *rout = x[0][0];
  *zout = x[0][1];
  locpow = get_localpower2d(data, numdata, x[0][0], x[0][1]);
  get_derivs2d(data, numdata, x[0][0], x[0][1], locpow, derivs);
  free(x[0]);
  free(x);
  return -y[0];
}


double max_rz_file(FILE *fftfile, double rin, double zin, \
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
  maxpow = max_rz_arr(filedata, filedatalen, rin_frac + filedatalen / 2, \
		      zin, rout, zout, derivs);
  *rout += startbin;
  free(filedata);
  return maxpow;
}
