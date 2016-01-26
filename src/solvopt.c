/* SOLVOPT version 1.1 (June, 1997)
   by Alexei Kuntsevich and Franz Kappel
   University of Graz, Austria

   The function SOLVOPT performs a modified version of Shor's r-algorithm in
   order to find a local minimum resp. maximum of a nonlinear function
   defined on the n-dimensional Euclidean space 
   or 
   a solution of a nonlinear constrained problem: 
   min { f(x): g(x) (<)= 0, g(x) in R(m), x in R(n) }
*/
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>


/* max(A,B) - larger (most +ve) of two numbers (generic) */
#define max(A,B) ((A)>(B)?(A):(B))
/* min(A,B) - smaller (least +ve) of two numbers (generic) */
#define min(A,B) ((A)<(B)?(A):(B))

#define errmes "\nSolvOpt error:"
#define wrnmes "\nSolvOpt warning:"
#define error2 "\nArgument X has to be a vector of dimension > 1."
#define error32 "\nFunction equals infinity at the point."
#define error42 "\nGradient equals infinity at the starting point."
#define error43 "\nGradient equals zero at the starting point."
#define error52 "\n<func> returns infinite value at the point."
#define error62 "\n<gradc> returns infinite vector at the point."
#define error63 "\n<gradc> returns zero vector at an infeasible point."
#define error5 "\nFunction is unbounded."
#define error6 "\nChoose another starting point."
#define warn1  "\nGradient is zero at the point, but stopping criteria are not fulfilled."
#define warn20 "\nNormal re-setting of a transformation matrix."
#define warn21 "\nRe-setting due to the use of a new penalty coefficient."
#define warn4  "\nIterations limit exceeded."
#define warn31 "\nThe function is flat in certain directions."
#define warn32 "\nTrying to recover by shifting insensitive variables."
#define warn09 "\nRe-run from recorded point."
#define warn08 "\nRavine with a flat bottom is detected."
#define termwarn0 "\nSolvOpt: Normal termination."
#define termwarn1 "\nSolvOpt: Termination warning:"
#define appwarn "\nThe above warning may be reasoned by inaccurate gradient approximation"
#define endwarn1 "\nPremature stop is possible. Try to re-run the routine from the obtained point."
#define endwarn2 "\nResult may not provide the optimum. The function apparently has many extremum points."
#define endwarn3 "\nResult may be inaccurate in the coordinates. The function is flat at the optimum."
#define endwarn4 "\nResult may be inaccurate in a function value. The function is extremely steep at the optimum."
#define allocerrstr "\nAllocation Error = "

double solvopt(unsigned short n,
               double x[],
               double fun(),
               void grad(), double options[], double func(), void gradc()
    )
{

/*solvopt  returns the optimum function value.

  Arguments to the function:
  n       is the space dimension,
  x       is an n-vector, the coordinates of the starting point
          at a call to the function and the optimizer at a regular return,
  fun     is the entry name of an external function which computes the value
          of the objective function 'fun' at a point x.
          synopsis: double fun(double x[])
  grad    is the entry name of an external function which computes the gradient
          vector of the objective function 'fun' at a point x.
          synopsis: void grad(double x[],double g[])
  options is a vector of optional parameters (see the description in SOLVOPT.H).
        Returned optional values:
        options[8], the number of iterations, if positive,
            or an abnormal stop code, if negative (see manual for more),
            -1: allocation error,
            -2: improper space dimension,
            -3: <fun> returns an improper value,
            -4: <grad> returns a zero or improper vector at the starting point,
            -5: <func> returns an improper value,
            -6: <gradc> returns an improper vector,
            -7: function is unbounded,
            -8: gradient is zero, but stopping criteria are not fulfilled,
            -9: iterations limit exceeded,
           -11: Premature stop is possible,
           -12: Result may not provide the true optimum,
           -13: Result may be inaccurate in view of a point.
           -14: Result may be inaccurate in view of a function value,
        options[9] , the number of objective function evaluations,    
        options[10], the number of gradient evaluations,
        options[11], the number of constraint function evaluations, and
        options[12], the number of constraint gradient evaluations.
____________________________________________________________________________*/

    double default_options[13] =
        { -1.0, 1.e-4, 1.e-6, 15000., 0.0, 1.e-8, 2.5, 1.e-11, 0.0, 0.0, 0.0,
        0.0, 0.0
    };
    void apprgrdn();
    unsigned short constr, app, appconstr = 0;
    unsigned short FsbPnt = 0, FsbPnt1 = 0, termflag, stopf;
    unsigned short stopping, dispwarn, Reset = 0, ksm, knan, obj;
    unsigned short kstore, knorms, k, kcheck, numelem;
    long ajp, ajpp;
    unsigned short ld, mxtc, termx, limxterm, nzero, krerun;
    unsigned short kflat, stepvanish, i, j, ni, ii, kd = 0, kj, kc, ip;
    unsigned short iterlimit, kg, k1, k2, kless = 0;
    short dispdata, warnno;
    double nsteps[3] = { 0.0, 0.0, 0.0 }, kk, nx;
    double gnorms[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ajb, ajs, des, dq, du20, du10, du03;
    double n_float, cnteps = 0.0;
    double low_bound, ZeroGrad, ddx, y;
    double lowxbound, lowfbound, detfr, detxr, grbnd;
    double f, fp = 0.0, fp1 = 0.0, fc = 0.0, f1, f2, fm, fopt, frec, fst, fp_rate;
    double PenCoef = 0.0, PenCoefNew;
    double gamma, w, wdef, h1, h, hp;
    double dx, ng, ngc, nng, ngt, nrmz, ng1, d, dd, laststep;
    double zero = 0., one = 1., two = 2., three = 3., four = 4., five =
        5., six = 6., seven = 7., eight = 8., nine = 9., ten = 10., hundr =
        100., infty = 1.e100, epsnorm = 1.e-15, epsnorm2 = 1.e-30, powerm12 = 1.e-12;
    double *B;                  /* space transformation matrix (allocatable)      */
    /* allocatable working arrays: */
    double *g, *g0, *g1, *gt, *gc, *z, *x1, *xopt, *xrec, *grec, *xx, *deltax;
    unsigned short *idx;
    char *endwarn = NULL;

/* Check the dimension: */
    if (n < 2) {
        printf(errmes);
        printf(error2);
        options[8] = -one;
        return (zero);
    }
    n_float = n;

/* Allocate the memory for working arrays: */

    B = (double *) calloc(n * n, sizeof(double));
    g = (double *) calloc(n, sizeof(double));
    g0 = (double *) calloc(n, sizeof(double));
    g1 = (double *) calloc(n, sizeof(double));
    gt = (double *) calloc(n, sizeof(double));
    gc = (double *) calloc(n, sizeof(double));
    z = (double *) calloc(n, sizeof(double));
    x1 = (double *) calloc(n, sizeof(double));
    xopt = (double *) calloc(n, sizeof(double));
    xrec = (double *) calloc(n, sizeof(double));
    grec = (double *) calloc(n, sizeof(double));
    xx = (double *) calloc(n, sizeof(double));
    deltax = (double *) calloc(n, sizeof(double));
    idx = (unsigned short *) calloc(n, sizeof(unsigned short));
    if (B == NULL || g == NULL || g0 == NULL || g1 == NULL || gt == NULL ||
        gc == NULL || z == NULL || x1 == NULL || xopt == NULL || xrec ==
        NULL || grec == NULL || xx == NULL || deltax == NULL || idx == NULL) {
        printf(allocerrstr);
        options[8] = -one;
        return (zero);
    }

/* ANALIZE THE ARGUMENTS PASSED
   User-supplied gradients: */
    if (grad == NULL)
        app = 1;
    else
        app = 0;
    if (func == NULL)
        constr = 0;
    else {
        constr = 1;
        if (gradc == NULL)
            appconstr = 1;
        else
            appconstr = 0;
    }
/* options: */
    for (i = 0; i <= 7; i++) {
        if (options[i] == zero)
            options[i] = default_options[i];
        else if (i == 1 || i == 2 || i == 5) {
            options[i] = max(options[i], powerm12);
            options[i] = min(options[i], one);
            if (i == 1)
                options[i] = max(options[i], options[8] * hundr);
        } else if (i == 6)
            options[6] = max(options[i], 1.5e0);
    }
    for (i = 8; i <= 12; i++)
        options[i] = zero;

    iterlimit = options[3];
/* Minimize resp. maximize the objective function:*/
    if (constr) {
        h1 = -one;              /* NLP: restricted to minimization */
        cnteps = options[5];
    } else if (options[0] < zero)
        h1 = -one;
    else
        h1 = one;
/* Multiplier for the matrix for the inverse of the space dilation: */
    wdef = one / options[6] - one;
/* Iterations counter: */
    k = 0;
/* Gamma control : */
    ajb = one + 0.1 / (n_float * n_float);
    ajp = 20;
    ajpp = ajp;
    ajs = 1.15e0;
    knorms = 0;
/* Display control : */
    if (options[4] <= zero) {
        dispdata = 0;
        if (options[4] == -one)
            dispwarn = 0;
        else
            dispwarn = 1;
    } else {
        dispdata = floor(options[4] + 0.1);
        dispwarn = 1;
    }
    ld = dispdata;
/* Stepsize control : */
    dq = 5.1;                   /* Step divider (at f_{i+1}>gamma*f_{i})  */
    du20 = two;
    du10 = 1.5;
    du03 = 1.05;                /* Step multipliers */
    kstore = 3;
    if (app)
        des = 6.3;              /* Desired number of steps per 1-D search */
    else
        des = 3.3;              /* Same for the case of analytical grads. */
    mxtc = 3;                   /* Number of trial cycles (wall detect)   */
    termx = 0;
    limxterm = 50;              /* Counter and limit for x-criterion      */

    ddx = max(1.e-11, options[7]);      /* stepsize for gradient approximation   */
    low_bound = -one + 1.e-4;   /* Lower bound cosine to detect a ravine */
    ZeroGrad = n_float * 1.e-16;        /* Lower bound for a gradient norm       */
    nzero = 0;                  /* Zero-gradient events counter          */

    lowxbound = max(options[1], 1.e-3); /* Low bound for the variables      */
    lowfbound = options[2] * options[2];        /* Lower bound for function values  */

    krerun = 0;                 /* Re-run events counter */
    detfr = options[2] * hundr; /* Relative error for f/f_{record} */
    detxr = options[1] * ten;   /* Relative error for norm(x)/norm(x_{record}) */
    warnno = 0;                 /* the number of a warn.mess. to end with */
    kflat = 0;                  /* counter for points of flatness */
    stepvanish = 0;             /* counter for vanished steps */
    stopf = 0;                  /* last-check flag */
/*  End of setting constants */
/*  End of the preamble */

/* COMPUTE THE OBJECTIVE FUNCTION (first time): */

    f = fun(x);
    options[9] += one;
    if (fabs(f) >= infty) {
        if (dispwarn) {
            printf(errmes);
            printf(error32);
            printf(error6);
        }
        options[8] = -three;
        goto endrun;
    }
    for (i = 0; i < n; i++)
        xrec[i] = x[i];
    frec = f;                   /* record the point */
    if (constr) {
        kless = 0;
        fp = f;
        fc = func(x);
        options[11] += one;
        if (fabs(fc) >= infty) {
            if (dispwarn) {
                printf(errmes);
                printf(error52);
                printf(error6);
            }
            options[8] = -five;
            goto endrun;
        }
        PenCoef = one;          /* first rough approximation */
        if (fc <= cnteps) {
            FsbPnt = 1;
            fc = zero;
        } else
            FsbPnt = 0;         /* infeasible point */
        f = f + PenCoef * fc;
    }
/* COMPUTE THE GRADIENT (first time): */
    if (app) {
        for (i = 0; i < n; i++)
            deltax[i] = h1 * ddx;
        obj = 1;
        if (constr)
            apprgrdn(n, g, x, fp, fun, deltax, obj);
        else
            apprgrdn(n, g, x, f, fun, deltax, obj);
        options[9] += n_float;
    } else {
        grad(x, g);
        options[10] += one;
    }
    ng = zero;
    for (i = 0; i < n; i++)
        ng += g[i] * g[i];
    ng = sqrt(ng);
    if (ng >= infty) {
        if (dispwarn) {
            printf(errmes);
            printf(error42);
            printf(error6);
        }
        options[8] = -four;
        goto endrun;
    } else if (ng < ZeroGrad) {
        if (dispwarn) {
            printf(errmes);
            printf(error43);
            printf(error6);
        }
        options[8] = -four;
        goto endrun;
    }

    if (constr) {
        if (!FsbPnt) {
            if (appconstr) {
                for (j = 0; j < n; j++) {
                    if (x[j] >= zero)
                        deltax[j] = ddx;
                    else
                        deltax[j] = -ddx;
                }
                obj = 0;
                apprgrdn(n, gc, x, fc, func, deltax, obj);
            } else
                gradc(x, gc);
            ngc = zero;
            for (i = 0; i < n; i++)
                ngc += gc[i] * gc[i];
            ngc = sqrt(ngc);
            if (ngc >= infty) {
                if (dispwarn) {
                    printf(errmes);
                    printf(error62);
                    printf(error6);
                }
                options[8] = -six;
                goto endrun;
            } else if (ngc < ZeroGrad) {
                if (dispwarn) {
                    printf(errmes);
                    printf(error63);
                }
                options[8] = -six;
                goto endrun;
            }
            ng = zero;
            for (i = 0; i < n; i++) {
                g[i] += PenCoef * gc[i];
                ng += g[i] * g[i];
            }
            ng = sqrt(ng);
        }
    }
    for (i = 0; i < n; i++)
        grec[i] = g[i];
    nng = ng;
/* INITIAL STEPSIZE : */
    d = zero;
    for (i = 0; i < n; i++) {
        if (d < fabs(x[i]))
            d = fabs(x[i]);
    }
    h = h1 * sqrt(options[1]) * d;      /* smallest possible stepsize */
    if (fabs(options[0]) != one)
        h = h1 * max(fabs(options[0]), fabs(h));        /* user-supplied stepsize */
    else
        h = h1 * max(one / log(ng + 1.1), fabs(h));     /* calculated stepsize */
/*--------------------------------------------------------------------
RESETTING LOOP */

    while (1) {
        kcheck = 0;             /* checkpoint counter */
        kg = 0;                 /* stepsizes stored */
        kj = 0;                 /* ravine jump counter */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                B[i * n + j] = zero;
            B[i * n + i] = one;
            g1[i] = g[i];
        }
        fst = f;
        dx = zero;
/*-----------------------------------------------------------------
  MAIN ITERATIONS  */

        while (1) {
            k += 1;
            kcheck += 1;
            laststep = dx;
            /* ADJUST GAMMA : */
            gamma = one + max(pow(ajb, (ajp - kcheck) * n), two * options[2]);
            gamma = min(gamma, pow(ajs, max(one, log10(nng + one))));

            /* Gradient in the transformed space (gt) : */
            ngt = zero;
            ng1 = zero;
            dd = zero;
            for (i = 0; i < n; i++) {
                d = zero;
                for (j = 0; j < n; j++)
                    d += B[j + i * n] * g[j];
                gt[i] = d;
                dd += d * g1[i];
                ngt += d * d;
                ng1 += g1[i] * g1[i];
            }
            ngt = sqrt(ngt);
            ng1 = sqrt(ng1);
            dd /= ngt * ng1;
            w = wdef;
            /* JUMPING OVER A RAVINE */
            if (dd < low_bound) {
                if (kj == 2)
                    for (i = 0; i < n; i++)
                        xx[i] = x[i];
                if (kj == 0)
                    kd = 4;
                kj += 1;
                w = -0.9;
                h *= two;
                if (kj > 2 * kd) {
                    kd += 1;
                    warnno = 1;
                    endwarn = endwarn1;
                    for (i = 0; i < n; i++) {
                        if (fabs(x[i] - xx[i]) < epsnorm * fabs(x[i])) {
                            if (dispwarn) {
                                printf(wrnmes);
                                printf(warn08);
                            }
                        }
                    }
                }
            } else
                kj = 0;

            /* DILATION : */
            nrmz = zero;
            for (i = 0; i < n; i++) {
                z[i] = gt[i] - g1[i];
                nrmz += z[i] * z[i];
            }
            nrmz = sqrt(nrmz);
            if (nrmz > epsnorm * ngt) {
                for (i = 0; i < n; i++)
                    z[i] /= nrmz;
                d = zero;
                for (i = 0; i < n; i++)
                    d += z[i] * gt[i];
                ng1 = zero;
                d *= w;
                for (i = 0; i < n; i++)
                    /* Make a space transformation:  g1=gt+w*(z*gt')*z: */
                {
                    dd = zero;
                    g1[i] = gt[i] + d * z[i];
                    ng1 += g1[i] * g1[i];
                    for (j = 0; j < n; j++)
                        dd += B[j * n + i] * z[j];
                    dd *= w;
                    /* New inverse matrix: B = B ( I + (1/alpha -1)zz' ) */
                    for (j = 0; j < n; j++)
                        B[j * n + i] += dd * z[j];
                }
                ng1 = sqrt(ng1);
            } else {
                for (i = 0; i < n; i++)
                    z[i] = zero;
                nrmz = zero;
            }
            for (i = 0; i < n; i++)
                gt[i] = g1[i] / ng1;
            /* Gradient in the non-transformed space: g0 = B' * gt   */
            for (i = 0; i < n; i++) {
                d = zero;
                for (j = 0; j < n; j++)
                    d += B[j * n + i] * gt[j];
                g0[i] = d;
            }

            /*  CHECK FOR THE NEED OF RESETTING   */
            if (kcheck > 1) {
                numelem = 0;
                for (i = 0; i < n; i++) {
                    if (fabs(g[i]) > ZeroGrad) {
                        idx[numelem] = i;
                        numelem += 1;
                    }
                }
                if (numelem > 0) {
                    grbnd = epsnorm * (numelem * numelem);
                    ii = 0;
                    for (i = 0; i < numelem; i++) {
                        j = idx[i];
                        if (fabs(g1[j]) <= fabs(g[j]) * grbnd)
                            ii += 1;
                    }
                    if (ii == n || nrmz == zero) {
                        if (dispwarn) {
                            printf(wrnmes);
                            printf(warn20);
                        }
                        if (fabs(fst - f) < fabs(f) * .01)
                            ajp -= 10 * n;
                        else
                            ajp = ajpp;
                        h = h1 * dx / three;
                        k = k - 1;
                        break;
                    }
                }
            }

            /* STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH  */

            for (i = 0; i < n; i++)
                xopt[i] = x[i];
            hp = h;
            fopt = f;
            k1 = 0;
            k2 = 0;
            ksm = 0;
            kc = 0;
            knan = 0;
            if (constr)
                Reset = 0;

            /* 1-D SEARCH   */

            while (1) {
                for (i = 0; i < n; i++)
                    x1[i] = x[i];
                f1 = f;
                if (f1 < zero)
                    dd = -one;
                else
                    dd = one;
                if (constr) {
                    FsbPnt1 = FsbPnt;
                    fp1 = fp;
                }

                /* Next point:   */
                for (i = 0; i < n; i++)
                    x[i] += hp * g0[i];
                ii = 0;
                for (i = 0; i < n; i++) {
                    if (fabs(x[i] - x1[i]) < fabs(x[i]) * epsnorm)
                        ii += 1;
                }

                /* COMPUTE THE FUNCTION VALUE AT A POINT:  */

                f = fun(x);
                options[9] += one;
                if (h1 * f >= infty) {
                    if (dispwarn) {
                        printf(errmes);
                        printf(error5);
                    }
                    options[8] = -seven;
                    goto endrun;
                }
                if (constr) {
                    fp = f;
                    fc = func(x);
                    options[11] += one;
                    if (fabs(fc) >= infty) {
                        if (dispwarn) {
                            printf(errmes);
                            printf(error52);
                            printf(error6);
                        }
                        options[8] = -five;
                        goto endrun;
                    }
                    if (fc <= cnteps) {
                        FsbPnt = 1;
                        fc = zero;
                    } else {
                        FsbPnt = 0;
                        fp_rate = fp - fp1;
                        if (fp_rate < -epsnorm) {
                            if (!FsbPnt1) {
                                d = zero;
                                for (i = 0; i < n; i++)
                                    d += (x[i] - x1[i]) * (x[i] - x1[i]);
                                d = sqrt(d);
                                PenCoefNew = -15. * fp_rate / d;
                                if (PenCoefNew > 1.2 * PenCoef) {
                                    PenCoef = PenCoefNew;
                                    Reset = 1;
                                    kless = 0;
                                    f += PenCoef * fc;
                                    break;
                                }
                            }
                        }
                        f += PenCoef * fc;
                    }
                }
                /* No function value at a point : */
                if (fabs(f) >= infty) {
                    if (dispwarn) {
                        printf(wrnmes);
                        printf(error32);
                    }
                    if (ksm || kc >= mxtc) {
                        options[8] = -three;
                        goto endrun;
                    } else {
                        k2 += 1;
                        k1 = 0;
                        hp /= dq;
                        for (i = 0; i < n; i++)
                            x[i] = x1[i];
                        f = f1;
                        knan = 1;
                        if (constr) {
                            FsbPnt = FsbPnt1;
                            fp = fp1;
                        }
                    }
                }
                /* STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM */
                else if (ii == n) {
                    stepvanish += 1;
                    if (stepvanish >= 5) {
                        if (dispwarn) {
                            printf(termwarn1);
                            printf(endwarn4);
                        }
                        options[8] = -14.;
                        goto endrun;
                    } else {
                        for (i = 0; i < n; i++)
                            x[i] = x1[i];
                        f = f1;
                        hp *= ten;
                        ksm = 1;
                        if (constr) {
                            FsbPnt = FsbPnt1;
                            fp = fp1;
                        }
                    }
                }
                /*  USE A SMALLER STEP:   */
                else if (h1 * f < h1 * pow(gamma, dd) * f1) {
                    if (ksm)
                        break;
                    k2 += 1;
                    k1 = 0;
                    hp /= dq;
                    for (i = 0; i < n; i++)
                        x[i] = x1[i];
                    f = f1;
                    if (constr) {
                        FsbPnt = FsbPnt1;
                        fp = fp1;
                    }
                    if (kc >= mxtc)
                        break;
                }
                /* 1-D OPTIMIZER IS LEFT BEHIND */
                else {
                    if (h1 * f <= h1 * f1)
                        break;
                    /* USE LARGER STEP */
                    k1 += 1;
                    if (k2 > 0)
                        kc += 1;
                    k2 = 0;
                    if (k1 >= 20)
                        hp *= du20;
                    else if (k1 >= 10)
                        hp *= du10;
                    else if (k1 >= 3)
                        hp *= du03;
                }
            }
            /* ------------------------  End of 1-D search  ------------------  */

            /* ADJUST THE TRIAL STEP SIZE : */
            dx = zero;
            for (i = 0; i < n; i++)
                dx += (xopt[i] - x[i]) * (xopt[i] - x[i]);
            dx = sqrt(dx);
            if (kg < kstore)
                kg += 1;
            if (kg >= 2)
                for (i = kg - 1; i > 0; i--)
                    nsteps[i] = nsteps[i - 1];
            d = zero;
            for (i = 0; i < n; i++)
                d += g0[i] * g0[i];
            d = sqrt(d);
            nsteps[0] = dx / (fabs(h) * d);
            kk = zero;
            d = zero;
            for (i = 1; i <= kg; i++) {
                dd = kg - i + 1;
                d += dd;
                kk += nsteps[i - 1] * dd;
            }
            kk /= d;
            if (kk > des) {
                if (kg == 1)
                    h *= kk - des + one;
                else
                    h *= sqrt(kk - des + one);
            } else if (kk < des)
                h *= sqrt(kk / des);

            if (ksm)
                stepvanish += 1;

            /* COMPUTE THE GRADIENT : */
            if (app) {
                for (j = 0; j < n; j++) {
                    if (g0[j] >= zero)
                        deltax[j] = h1 * ddx;
                    else
                        deltax[j] = -h1 * ddx;
                }
                obj = 1;
                if (constr)
                    apprgrdn(n, g, x, fp, fun, deltax, obj);
                else
                    apprgrdn(n, g, x, f, fun, deltax, obj);
                options[9] += n_float;
            } else {
                grad(x, g);
                options[10] += one;
            }
            ng = zero;
            for (i = 0; i < n; i++)
                ng += g[i] * g[i];
            ng = sqrt(ng);
            if (ng >= infty) {
                if (dispwarn) {
                    printf(errmes);
                    printf(error42);
                }
                options[8] = -four;
                goto endrun;
            } else if (ng < ZeroGrad) {
                if (dispwarn) {
                    printf(wrnmes);
                    printf(warn1);
                }
                ng = ZeroGrad;
            }
            /* Constraints: */
            if (constr) {
                if (!FsbPnt) {
                    if (ng < 0.01 * PenCoef) {
                        kless += 1;
                        if (kless >= 20) {
                            PenCoef /= ten;
                            Reset = 1;
                            kless = 0;
                        }
                    } else
                        kless = 0;
                    if (appconstr) {
                        for (j = 0; j < n; j++) {
                            if (x[j] >= zero)
                                deltax[j] = ddx;
                            else
                                deltax[j] = -ddx;
                        }
                        obj = 0;
                        apprgrdn(n, gc, x, fc, func, deltax, obj);
                        options[11] += n_float;
                    } else {
                        gradc(x, gc);
                        options[12] += one;
                    }
                    ngc = zero;
                    for (i = 0; i < n; i++)
                        ngc += gc[i] * gc[i];
                    ngc = sqrt(ngc);
                    if (ngc >= infty) {
                        if (dispwarn) {
                            printf(errmes);
                            printf(error62);
                        }
                        options[8] = -six;
                        goto endrun;
                    } else if (ngc < ZeroGrad && !appconstr) {
                        if (dispwarn) {
                            printf(errmes);
                            printf(error63);
                        }
                        options[8] = -six;
                        goto endrun;
                    }
                    ng = zero;
                    for (i = 0; i < n; i++) {
                        g[i] += PenCoef * gc[i];
                        ng += g[i] * g[i];
                    }
                    ng = sqrt(ng);
                    if (Reset) {
                        if (dispwarn) {
                            printf(wrnmes);
                            printf(warn21);
                        }
                        h = h1 * dx / three;
                        k -= 1;
                        nng = ng;
                        break;
                    }
                }
            }
            /* new record */
            if (h1 * f > h1 * frec) {
                frec = f;
                for (i = 0; i < n; i++) {
                    xrec[i] = x[i];
                    grec[i] = g[i];
                }
            }
            /* average gradient norm */
            if (ng > ZeroGrad) {
                if (knorms < 10)
                    knorms += 1;
                if (knorms >= 2) {
                    for (i = knorms - 1; i > 0; i--)
                        gnorms[i] = gnorms[i - 1];
                }
                gnorms[0] = ng;
                nng = one;
                for (i = 0; i < knorms; i++)
                    nng *= gnorms[i];
                nng = pow(nng, one / knorms);
            }
            /* Norm of X: */
            nx = zero;
            for (i = 0; i < n; i++)
                nx += x[i] * x[i];
            nx = sqrt(nx);

   /*-----------------------------------------------------------------
   DISPLAY THE CURRENT VALUES: */
            if (k == ld) {
                printf("\nIteration # ..... Function Value ..... "
                       "Step Value ..... Gradient Norm"
                       "\n     %5i     %13.5g      %13.5g       %13.5g", k, f, dx,
                       ng);
                ld += dispdata;
            }

   /*-----------------------------------------------------------------
   CHECK THE STOPPING CRITERIA: */

            termflag = 1;
            if (constr) {
                if (!FsbPnt)
                    termflag = 0;
            }
            if (kcheck <= 5 || (kcheck <= 12 && ng > one))
                termflag = 0;
            if ((kc >= mxtc) || knan)
                termflag = 0;
            /* ARGUMENT : */
            if (termflag) {
                ii = 0;
                stopping = 1;
                for (i = 0; i < n; i++) {
                    if (fabs(x[i]) >= lowxbound) {
                        idx[ii] = i;
                        ii += 1;
                        if (fabs(xopt[i] - x[i]) > options[1] * fabs(x[i]))
                            stopping = 0;
                    }
                }
                if (ii == 0 || stopping) {
                    stopping = 1;
                    termx += 1;
                    d = zero;
                    for (i = 0; i < n; i++)
                        d += (x[i] - xrec[i]) * (x[i] - xrec[i]);
                    d = sqrt(d);
                    /* FUNCTION : */
                    if (fabs(f - frec) > detfr * fabs(f) &&
                        fabs(f - fopt) >= options[2] * fabs(f) &&
                        krerun <= 3 && !constr) {
                        stopping = 0;
                        if (ii > 0) {
                            for (i = 0; i < ii; i++) {
                                j = idx[i];
                                if (fabs(xrec[j] - x[j]) > detxr * fabs(x[j])) {
                                    stopping = 1;
                                    break;
                                }
                            }
                        }
                        if (stopping) {
                            if (dispwarn) {
                                printf(wrnmes);
                                printf(warn09);
                            }
                            ng = zero;
                            for (i = 0; i < n; i++) {
                                x[i] = xrec[i];
                                g[i] = grec[i];
                                ng += g[i] * g[i];
                            }
                            ng = sqrt(ng);
                            f = frec;
                            krerun += 1;
                            h = h1 * max(dx, detxr * nx) / krerun;
                            warnno = 2;
                            endwarn = endwarn2;
                            break;
                        } else
                            h *= ten;
                    } else if (fabs(f - frec) > options[2] * fabs(f) &&
                               d < options[1] * nx && constr) {
                    } else if (fabs(f - fopt) <= options[2] * fabs(f) ||
                               fabs(f) <= lowfbound ||
                               (fabs(f - fopt) <= options[2] && termx >= limxterm)) {
                        if (stopf) {
                            if (dx <= laststep) {
                                if (warnno == 1 && ng < sqrt(options[2]))
                                    warnno = 0;
                                if (!app) {
                                    for (i = 0; i < n; i++) {
                                        if (fabs(g[i]) <= epsnorm2) {
                                            warnno = 3;
                                            endwarn = endwarn3;
                                            break;
                                        }
                                    }
                                }
                                if (warnno != 0) {
                                    options[8] = -warnno - ten;
                                    if (dispwarn) {
                                        printf(termwarn1);
                                        printf(endwarn);
                                        if (app)
                                            printf(appwarn);
                                    }
                                } else {
                                    options[8] = k;
                                    if (dispwarn)
                                        printf(termwarn0);
                                }
                                goto endrun;
                            }
                        } else
                            stopf = 1;
                    }

                    else if (dx < powerm12 * max(nx, one) && termx >= limxterm) {
                        options[8] = -14.;
                        if (dispwarn) {
                            printf(termwarn1);
                            printf(endwarn4);
                            if (app)
                                printf(appwarn);
                        }
                        f = frec;
                        for (i = 0; i < n; i++)
                            x[i] = xrec[i];
                        goto endrun;
                    }
                }               /* stopping */
            }
            /* termflag */
            /* ITERATIONS LIMIT */
            if (k == iterlimit) {
                options[8] = -nine;
                if (dispwarn) {
                    printf(wrnmes);
                    printf(warn4);
                }
                goto endrun;
            }
            /* ------------ end of the check ---------------- */
            /* ZERO GRADIENT : */
            if (constr) {
                if (ng <= ZeroGrad) {
                    if (dispwarn) {
                        printf(termwarn1);
                        printf(warn1);
                    }
                    options[8] = -eight;
                    goto endrun;
                }
            } else {
                if (ng <= ZeroGrad) {
                    nzero += 1;
                    if (dispwarn) {
                        printf(wrnmes);
                        printf(warn1);
                    }
                    if (nzero >= 3) {
                        options[8] = -eight;
                        goto endrun;
                    }
                    for (i = 0; i < n; i++)
                        g0[i] *= -h / two;
                    for (i = 1; i <= 10; i++) {
                        for (j = 0; j < n; j++)
                            x[j] += g0[j];
                        f = fun(x);
                        options[9] += one;
                        if (fabs(f) >= infty) {
                            if (dispwarn) {
                                printf(errmes);
                                printf(error32);
                            }
                            options[8] = -three;
                            goto endrun;
                        }
                        if (app) {
                            for (j = 0; j < n; j++) {
                                if (g0[j] >= zero)
                                    deltax[j] = h1 * ddx;
                                else
                                    deltax[j] = -h1 * ddx;
                            }
                            obj = 1;
                            apprgrdn(n, g, x, f, fun, deltax, obj);
                            options[9] += n_float;
                        } else {
                            grad(x, g);
                            options[10] += one;
                        }
                        ng = zero;
                        for (j = 0; j < n; j++)
                            ng += g[j] * g[j];
                        ng = sqrt(ng);
                        if (ng >= infty) {
                            if (dispwarn) {
                                printf(errmes);
                                printf(error42);
                            }
                            options[8] = -four;
                            goto endrun;
                        }
                        if (ng > ZeroGrad)
                            break;
                    }
                    if (ng <= ZeroGrad) {
                        if (dispwarn) {
                            printf(termwarn1);
                            printf(warn1);
                        }
                        options[8] = -eight;
                        goto endrun;
                    }
                    h = h1 * dx;
                    break;
                }
            }
            /* FUNCTION IS FLAT AT THE POINT : */
            if (!constr &&
                fabs(f - fopt) < fabs(fopt) * options[2] && kcheck > 5 && ng < one) {
                ni = 0;
                for (i = 0; i < n; i++) {
                    if (fabs(g[i]) <= epsnorm2) {
                        idx[ni] = i;
                        ni += 1;
                    }
                }
                if (ni >= 1 && ni <= n / 2 && kflat <= 3) {
                    kflat += 1;
                    if (dispwarn) {
                        printf(wrnmes);
                        printf(warn31);
                    }
                    warnno = 1;
                    endwarn = endwarn1;
                    for (i = 0; i < n; i++)
                        x1[i] = x[i];
                    fm = f;
                    for (i = 0; i < ni; i++) {
                        j = idx[i];
                        f2 = fm;
                        y = x[j];
                        if (y == zero)
                            x1[j] = one;
                        else if (fabs(y) < one) {
                            if (y < 0)
                                x1[j] = -one;
                            else
                                x1[j] = one;
                        } else
                            x1[j] = y;
                        for (ip = 1; ip <= 20; i++) {
                            x1[j] /= 1.15;
                            f1 = fun(x1);
                            options[9] += one;
                            if (fabs(f1) < infty) {
                                if (h1 * f1 > h1 * fm) {
                                    y = x1[j];
                                    fm = f1;
                                } else if (h1 * f2 > h1 * f1)
                                    break;
                                else if (f2 == f1)
                                    x1[j] /= 1.5;
                                f2 = f1;
                            }
                        }
                        x1[j] = y;
                    }
                    if (h1 * fm > h1 * f) {
                        if (app) {
                            for (j = 0; j < n; j++)
                                deltax[j] = h1 * ddx;
                            obj = 1;
                            apprgrdn(n, gt, x1, fm, fun, deltax, obj);
                            options[9] += n_float;
                        } else {
                            grad(x1, gt);
                            options[10] += one;
                        }
                        ngt = zero;
                        for (i = 0; i < n; i++)
                            ngt += gt[i] * gt[i];
                        if (ngt > epsnorm2 || ngt < infty) {
                            if (dispwarn)
                                printf(warn32);
                            for (i = 0; i < n; i++) {
                                x[i] = x1[i];
                                g[i] = gt[i];
                            }
                            ng = ngt;
                            f = fm;
                            h = h1 * dx / three;
                            options[2] /= five;
                            break;
                        }       /* regular gradient */
                    }           /* a better value has been found */
                }               /* function is flat */
            }
            /* pre-conditions are fulfilled */
        }                       /* end of the iteration cycle */
    }                           /*  end of the resetting cycle */
  endrun:
    /* deallocate working arrays: */
    free(idx);
    free(deltax);
    free(xx);
    free(grec);
    free(xrec);
    free(xopt);
    free(x1);
    free(z);
    free(gc);
    free(gt);
    free(g1);
    free(g0);
    free(g);
    free(B);
    return (f);
}
