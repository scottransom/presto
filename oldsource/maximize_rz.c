#include "presto.h"
#include "nrutil.h"
#include "nr.h"

#define NDIM 2
#define FTOL 3.0e-8
#define ITMAX 200
#define TOL 3.0e-8
#define BRENTITMAX 100
#define CGOLD 0.3819660112501051518
#define ZEPS 1.0e-14
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define GOLD  1.6180339887498948482
#define GLIMIT 100.0
#define TINY 1.0e-20

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Function declarations */

double power_call_rz(double rz[]);
double power_call_rz_file(double rz[]);
double brent(double ax, double bx, double cx,
	     double (*f) (double), double tol, double *xmin);
double f1dim(double x);
void linmin(double p[], double xi[], int n, double *fret,
	    double (*func) (double[]));
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	    double *fc, double (*func) (double));
void powell(double p[], double **xi, int n, double ftol, int *iter,
	    double *fret, double (*func) (double[]));


/*  Global variables */

extern FILE *global_maxfile;
extern float *global_maxarr;
extern long global_numarrbins, global_m;
extern int global_ncom;
extern double *global_pcom, *global_xicom, (*global_minfunc) (double[]);


/*  Maximization function used with an array */

double power_call_rz(double rz[])
{
  double drl, dim;
  float rl, im;

  rz_interp(global_maxarr, global_numarrbins, rz[1], rz[2], \
	    LOWACC, &global_m, &rl, &im);
  drl = (double) rl;
  dim = (double) im;
  return 0.0 - (drl * drl + dim * dim);
}


/*  Maximization function used with a file */

double power_call_rz_file(double rz[])
{
  double drl, dim;
  float rl, im;

  rz_interp_file(global_maxfile, rz[1], rz[2], \
		 LOWACC, &global_m, &rl, &im);
  drl = (double) rl;
  dim = (double) im;
  return 0.0 - (drl * drl + dim * dim);
}


double max_rz_file(FILE * file, double rin, double zin, double
		   *rout, double *zout, rderivs * derivs)
{
  int iter;
  double maxpower;
  double **xi, *p;
  float locpow;

  global_maxfile = file;

  /*  Now prep the maximization */

  global_m = get_z_resp_width(zin, LOWACC);
  p = dvector(1, NDIM);
  p[1] = rin;
  p[2] = zin;
  xi = dmatrix(1, NDIM, 1, NDIM);
  xi[1][1] = 0.1;
  xi[1][2] = 0.0;		/* Initial direction vector for freq */

  xi[2][1] = 0.0;
  xi[2][2] = 0.4;		/* Initial direction vector for fdot */

  powell(p, xi, NDIM, FTOL, &iter, &maxpower, power_call_rz_file);

  /*  Restart at minimum to ensure no "bad" steps */

  global_m = get_z_resp_width(zin, HIGHACC);
  xi[1][1] = 0.01;
  xi[1][2] = 0.0;		/* Initial direction vector for freq */

  xi[2][1] = 0.0;
  xi[2][2] = 0.04;		/* Initial direction vector for fdot */

  powell(p, xi, NDIM, FTOL, &iter, &maxpower, power_call_rz_file);
  *rout = p[1];
  *zout = p[2];

  /* The following calculate derivatives and local power at the peak */
  /*    See characteristics.c for their difinitions                  */

  locpow = get_localpower2d_file(file, p[1], p[2]);
  get_derivs2d_file(derivs, file, p[1], p[2], locpow);
  free_dvector(p, 1, NDIM);
  free_dmatrix(xi, 1, NDIM, 1, NDIM);
  return -maxpower;
}


double max_rz_arr(float *data, long numfreqs, double rin, double zin, \
		  double *rout, double *zout, rderivs * derivs)
{
  int iter;
  double maxpower;
  double **xi, *p;
  float locpow;

  global_maxarr = data;
  global_numarrbins = numfreqs;

  /*  Now prep the maximization */

  global_m = get_z_resp_width(zin, LOWACC);
  p = dvector(1, NDIM);
  p[1] = rin;
  p[2] = zin;
  xi = dmatrix(1, NDIM, 1, NDIM);
  xi[1][1] = 0.1;
  xi[1][2] = 0.0;		/* Initial direction vector for freq */

  xi[2][1] = 0.0;
  xi[2][2] = 0.4;		/* Initial direction vector for fdot */

  powell(p, xi, NDIM, FTOL, &iter, &maxpower, power_call_rz);

  /*  Restart at minimum to ensure no "bad" steps */

  global_m = get_z_resp_width(zin, HIGHACC);
  xi[1][1] = 0.01;
  xi[1][2] = 0.0;		/* Initial direction vector for freq */

  xi[2][1] = 0.0;
  xi[2][2] = 0.04;		/* Initial direction vector for fdot */

  powell(p, xi, NDIM, FTOL, &iter, &maxpower, power_call_rz);
  *rout = p[1];
  *zout = p[2];

  /* The following calculate derivatives and local power at the peak */
  /*    See characteristics.c for their difinitions                  */

  locpow = get_localpower2d(data, numfreqs, p[1], p[2]);
  get_derivs2d(derivs, data, numfreqs, p[1], p[2], locpow);
  free_dvector(p, 1, NDIM);
  free_dmatrix(xi, 1, NDIM, 1, NDIM);
  free_vector(global_maxarr, 0, 2 * global_numarrbins - 1);
  return -maxpower;
}


void powell(double p[], double **xi, int n, double ftol, int *iter,
	    double *fret, double (*func) (double[]))
{
  int i, ibig, j;
  double del, fp, fptt, t, *pt, *ptt, *xit;

  pt = dvector(1, n);
  ptt = dvector(1, n);
  xit = dvector(1, n);
  *fret = (*func) (p);
  for (j = 1; j <= n; j++)
    pt[j] = p[j];
  for (*iter = 1;; ++(*iter)) {
    fp = (*fret);
    ibig = 0;
    del = 0.0;
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++)
	xit[j] = xi[j][i];
      fptt = (*fret);
      linmin(p, xit, n, fret, func);
      if (fabs(fptt - (*fret)) > del) {
	del = fabs(fptt - (*fret));
	ibig = i;
      }
    }
    if (2.0 * fabs(fp - (*fret)) <= ftol * (fabs(fp) + fabs(*fret))) {
      free_dvector(xit, 1, n);
      free_dvector(ptt, 1, n);
      free_dvector(pt, 1, n);
      return;
    }
    if (*iter == ITMAX){
      printf("\n powell() exceeding maximum iterations.");
    }
    for (j = 1; j <= n; j++) {
      ptt[j] = 2.0 * p[j] - pt[j];
      xit[j] = p[j] - pt[j];
      pt[j] = p[j];
    }
    fptt = (*func) (ptt);
    if (fptt < fp) {
      t = 2.0 * (fp - 2.0 * (*fret) + fptt) * \
	  DSQR(fp - (*fret) - del) - del * DSQR(fp - fptt);
      if (t < 0.0) {
	linmin(p, xit, n, fret, func);
	for (j = 1; j <= n; j++) {
	  xi[j][ibig] = xi[j][n];
	  xi[j][n] = xit[j];
	}
      }
    }
  }
}


void linmin(double p[], double xi[], int n, double *fret, \
	    double (*func) (double[]))
{
  int j;
  double xx, xmin, fx, fb, fa, bx, ax;

  global_ncom = n;
  global_pcom = dvector(1, n);
  global_xicom = dvector(1, n);
  global_minfunc = func;
  for (j = 1; j <= n; j++) {
    global_pcom[j] = p[j];
    global_xicom[j] = xi[j];
  }
  ax = 0.0;
  xx = 1.0;
  mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
  *fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
  for (j = 1; j <= n; j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(global_xicom, 1, n);
  free_dvector(global_pcom, 1, n);
}


double f1dim(double x)
{
  int j;
  double f, *xt;

  xt = dvector(1, global_ncom);
  for (j = 1; j <= global_ncom; j++)
    xt[j] = global_pcom[j] + x * global_xicom[j];
  f = (*global_minfunc) (xt);
  free_dvector(xt, 1, global_ncom);
  return f;
}


double brent(double ax, double bx, double cx, double (*f) (double), double tol,
	     double *xmin)
{
  int iter;
  double a, b, d = 0.0, etemp, fu, fv, fw, fx, p, q;
  double e = 0.0, r, tol1, tol2, u, v, w, x, xm;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = (*f) (x);
  for (iter = 1; iter <= BRENTITMAX; iter++) {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
	p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
	d = CGOLD * (e = (x >= xm ? a - x : b - x));
      else {
	d = p / q;
	u = x + d;
	if (u - a < tol2 || b - u < tol2)
	  d = SIGN(tol1, xm - x);
      }
    } else {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }
    u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
    fu = (*f) (u);
    if (fu <= fx) {
      if (u >= x)
	a = x;
      else
	b = x;
      SHFT(v, w, x, u)
	  SHFT(fv, fw, fx, fu)
    } else {
      if (u < x)
	a = u;
      else
	b = u;
      if (fu <= fw || w == x) {
	v = w;
	w = u;
	fv = fw;
	fw = fu;
      } else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
  }
  printf("\nToo many iterations in brent().\n");
  return 0.0;
}


void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	    double (*func) (double))
{
  double ulim, u, r, q, fu, dum;

  *fa = (*func) (*ax);
  *fb = (*func) (*bx);
  if (*fb > *fa) {
    SHFT(dum, *ax, *bx, dum)
	SHFT(dum, *fb, *fa, dum)
  }
  *cx = (*bx) + GOLD * (*bx - *ax);
  *fc = (*func) (*cx);
  while (*fb > *fc) {
    r = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
	(2.0 * SIGN(DMAX(fabs(q - r), TINY), q - r));
    ulim = (*bx) + GLIMIT * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0) {
      fu = (*func) (u);
      if (fu < *fc) {
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      } else if (fu > *fb) {
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx) + GOLD * (*cx - *bx);
      fu = (*func) (u);
    } else if ((*cx - u) * (u - ulim) > 0.0) {
      fu = (*func) (u);
      if (fu < *fc) {
	SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
	    SHFT(*fb, *fc, fu, (*func) (u))
      }
    } else if ((u - ulim) * (ulim - *cx) >= 0.0) {
      u = ulim;
      fu = (*func) (u);
    } else {
      u = (*cx) + GOLD * (*cx - *bx);
      fu = (*func) (u);
    }
    SHFT(*ax, *bx, *cx, u)
	SHFT(*fa, *fb, *fc, fu)
  }
}


#undef ITMAX
#undef TOL
#undef FTOL
#undef BRENTITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef GOLD
#undef GLIMIT
#undef TINY
