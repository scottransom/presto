#include "presto.h"
#define FTOL 3.0e-8

/* Function declarations */

double power_call_r(double r);
double power_call_r_file(double r);
double brent(double ax, double bx, double cx,
	     double (*f) (double), double tol, double *xmin);

/*  Global variables */
FILE *global_maxfile;
float *global_maxarr;
long global_numarrbins, global_m;
int global_ncom;
double *global_pcom, *global_xicom, (*global_minfunc) (double[]);

/*  Maximization function used with an array */

double power_call_r(double r)
{
  double drl, dim;
  float rl, im;

  rz_interp(global_maxarr, global_numarrbins, r, 0.0, \
	    LOWACC, &global_m, &rl, &im);
  drl = (double) rl;
  dim = (double) im;
  return 0.0 - (drl * drl + dim * dim);
}


/*  Maximization function used with a file */

double power_call_r_file(double r)
{
  double drl, dim;
  float rl, im;

  rz_interp_file(global_maxfile, r, 0.0, \
		 LOWACC, &global_m, &rl, &im);
  drl = (double) rl;
  dim = (double) im;
  return 0.0 - (drl * drl + dim * dim);
}


double max_r_file(FILE * file, double rin, double *rout, \
		  rderivs * derivs)
{
  double maxpower, ax, bx, cx, xmin = 0.0;
  float locpow;

  global_maxfile = file;

  /*  Now prep and do the maximization */

  global_m = get_r_resp_width(LOWACC);
  bx = rin;
  ax = bx - 0.55;
  cx = bx + 0.55;

  brent(ax, bx, cx, power_call_r_file, FTOL, &xmin);

  /*  Restart at minimum to ensure no "bad" steps */

  global_m = get_r_resp_width(HIGHACC);
  ax = xmin - 0.05;
  bx = xmin;
  cx = xmin + 0.05;

  brent(ax, bx, cx, power_call_r_file, FTOL, &xmin);
  maxpower = power_call_r_file(xmin);

  *rout = xmin;

  /* The following calculate derivatives and local power at the peak */
  /*    See characteristics.c for their difinitions                  */

  locpow = get_localpower2d_file(file, xmin, 0.0);
  get_derivs2d_file(derivs, file, xmin, 0.0, locpow);
  return -maxpower;
}


double max_r_arr(float *data, long numfreqs, double rin, 
		 double *rout, rderivs * derivs, double locpow)
{
  double maxpower, ax, bx, cx, xmin = 0.0;

  global_maxarr = data;
  global_numarrbins = numfreqs;

  /*  Now prep and do the maximization */

  global_m = get_r_resp_width(LOWACC);
  bx = rin;
  ax = bx - 0.55;
  cx = bx + 0.55;

  brent(ax, bx, cx, power_call_r, FTOL, &xmin);

  /*  Restart at minimum to ensure no "bad" steps */

  global_m = get_r_resp_width(HIGHACC);
  ax = xmin - 0.05;
  bx = xmin;
  cx = xmin + 0.05;

  brent(ax, bx, cx, power_call_r, FTOL, &xmin);
  maxpower = power_call_r(xmin);

  *rout = xmin;

  /* The following calculate derivatives at the peak        */
  /*    See characteristics.c for their difinitions         */

  get_derivs2d(derivs, data, global_numarrbins, \
	       xmin, 0.0, locpow);
  return -maxpower;
}

#undef FTOL
