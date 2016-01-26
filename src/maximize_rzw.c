#include "presto.h"

#define ZSCALE 4.0
#define WSCALE 24.0

/* Some Static-Global Variables */

static fcomplex *maxdata;
static int nummaxdata, max_kern_half_width, num_funct_calls = 0;
static double solvopt_options[13] = {
    -1.0, 1.e-4, 1.e-6, 15000.0, -1.0, 1.e-8, 2.5, 0.003,
    0.0, 0.0, 0.0, 0.0, 0.0
};

/* Explanation:
   solvopt_options[0]:  If < 0 minimize, if > 0 maximize the function.
                        The magnitude is related to the initial stepsize.
                        Default = -1.0
   solvopt_options[1]:  Relative error of the argument in terms of the 
                        max gradient norm.  The smoother your function
                        is, the smaller this number can be.
                        Default = 1.e-4
   solvopt_options[2]:  Relative error of the function value.
                        Default = 1.e-6
   solvopt_options[3]:  Limit for the number of iterations.
                        Default = 15000
   solvopt_options[4]:  Controls the display of intermediate results and
                        error/warning messages.
                        Default = 0 (No intermediate messages)
   solvopt_options[5]:  Admissible maximal residual for a set of constraints
                        Default = 1.e-8
   solvopt_options[6]:  The coefficient of space dilation.
                        Default = 2.5
   solvopt_options[7]:  The lower bound for the stepsize used for the finite
                        difference approximation of gradients.
                        Default = 1.e-11
   Returned optional values:
      solvopt_options[8]:  The number of iterations. (Returns < 0 if error.)
         -1 = allocation error,
         -2 = improper space dimension,
         -3 = <fun> returns an improper value,
         -4 = <grad> returns a zero or improper vector at the starting point,
         -5 = <func> returns an improper value,
         -6 = <gradc> returns an improper vector,
         -7 = function is unbounded,
         -8 = gradient is zero, but stopping criteria are not fulfilled,
         -9 = iterations limit exceeded,
         -11 = Premature stop is possible,
         -12 = Result may not provide the true optimum,
         -13 = Result may be inaccurate in view of a point.
         -14 = Result may be inaccurate in view of a function value,
      solvopt_options[9]:  The number of objective function evaluations.
      solvopt_options[10]: The number of gradient evaluations.
      solvopt_options[11]: The number of constraint function evaluations.
      solvopt_options[12]: The number of constraint gradient evaluations.
*/

/* Function definition for the minimization routine */

extern double solvopt(unsigned short n, double x[], double fun(),
                      void grad(), double options[], double func(), void gradc());

static double power_call_rzw(double rzw[])
/*  Maximization function used with an array */
{
    double powargr, powargi;
    fcomplex ans;

    num_funct_calls++;
    rzw_interp(maxdata, nummaxdata, rzw[0], rzw[1] * ZSCALE,
               rzw[2] * WSCALE, max_kern_half_width, &ans);
    powargr = (double) ans.r;
    powargi = (double) ans.i;
    return -POWER(powargr, powargi);
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
