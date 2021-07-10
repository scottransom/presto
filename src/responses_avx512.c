#include "presto.h"
#include <immintrin.h>
#include <math.h>

#define MIN_NUMDATA 131072
#define MIN_NUMORBPTS 2049      /* This should be a power-of-two + 1 */

/* Function declarations */
// int fresnl(double xxa, double *ssa, double *cca);

int r_resp_halfwidth(presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a standard Fourier signal.                  */
  /*  Arguments:                                                       */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
    if (accuracy == HIGHACC) {
        return ((NUMFINTBINS * 3) + (NUMLOCPOWAVG >> 1) + DELTAAVGBINS);
    } else {
        return NUMFINTBINS;
    }
}


int z_resp_halfwidth(double z, presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
  /*  a constant frequency derivative)                                 */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
    int m = (int) (0.5 * 1.1 * fabs(z));
    if (accuracy == HIGHACC) {
        m += NUMFINTBINS * 3;
    } else {
        m += NUMFINTBINS;
    }
    return m;
}

int w_resp_halfwidth(double z, double w, presto_interp_acc accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with an f-dot that (i.e    */
  /*  varies linearly in time -- a constant f-dotdot)                  */
  /*  Arguments:                                                       */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
    if (fabs(w) < 1.0e-7)
        return z_resp_halfwidth(z, accuracy);
    double r0 = 0.5 * (w / 6.0 - z); // Starting deviation from r_avg
    double r1 = 0.5 * (w / 6.0 + z); // Ending deviation from r_avg
    // We need to know the maximum deviation from r_avg
    double maxdev = fabs(r0) > fabs(r1) ? fabs(r0) : fabs(r1);
    // If the extrema of the parabola is within 0 < u < 1, then
    // it will be a new freq minimum or maximum
    double u_ext = 0.5 - z / w;
    if (u_ext > 0.0 && u_ext < 1.0) {
        double z0 = z - w / 2.0; // Starting z
        // Value of r at the extremum
        double r_ext =  0.5 * w * u_ext * u_ext + z0 * u_ext + r0;
        maxdev = fabs(r_ext) > maxdev ? fabs(r_ext) : maxdev;
    }
    if (accuracy == HIGHACC) {
        return (int) (1.1 * maxdev) + NUMFINTBINS * 3;
    } else {
        return (int) (1.1 * maxdev) + NUMFINTBINS;
    }
}


void binary_velocity(double T, orbitparams * orbit, double *minv, double *maxv)
  /*  Return the minimum and maximum orbital velocities of a pulsar    */
  /*  during an observation as a fraction of the speed of light.       */
  /*  Arguments:                                                       */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */
{
    *minv = 1.0;
    *maxv = -1.0;
    if (T >= orbit->p) {

        double c1, c2;

        c1 = TWOPI * orbit->x / (orbit->p * sqrt(1.0 - orbit->e * orbit->e));
        c2 = orbit->e * cos(orbit->w * DEGTORAD);
        *maxv = c1 * (c2 + 1.0);
        *minv = c1 * (c2 - 1.0);

    } else {

        double dtb, startE, *E;
        int ii, numpoints = 1025;
        orbitparams orb;

        dtb = T / (double) (numpoints - 1);
        orb.p = orbit->p;
        orb.x = orbit->x;
        orb.e = orbit->e;
        orb.w = orbit->w;
        orb.t = orbit->t;
        startE = keplers_eqn(orb.t, orb.p, orb.e, 1.0E-15);
        E = dorbint(startE, numpoints, dtb, &orb);
        E_to_v(E, numpoints, &orb);
        for (ii = 0; ii < numpoints; ii++) {
            E[ii] *= 1000.0 / SOL;
            if (E[ii] < *minv)
                *minv = E[ii];
            if (E[ii] > *maxv)
                *maxv = E[ii];
        }
        vect_free(E);
    }
    if (*maxv < *minv) {
        printf("Something is wrong in binary_velocity()\n");
    }
}

int bin_resp_halfwidth(double ppsr, double T, orbitparams * orbit)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a pulsar in a binary orbit.                 */
  /*  Arguments:                                                       */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2 * 'numbetween' to get the   */
  /*    length of the array required to hold such a kernel.            */
{
    double maxv, minvel, maxvel, maxdevbins;
    int retval;

    binary_velocity(T, orbit, &minvel, &maxvel);
    maxv = (fabs(minvel) > fabs(maxvel)) ? minvel : maxvel;
    maxdevbins = fabs(T * maxv / (ppsr * (1.0 + maxv)));
    retval = (int) floor(1.1 * maxdevbins + 0.5);
    return (retval < NUMFINTBINS) ? NUMFINTBINS : retval;
}


fcomplex *gen_r_response(double roffset, int numbetween, int numkern)
  /*  Generate a complex response function for Fourier interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{
    int ii;
    double tmp, sinc, s, c, alpha, beta, delta, startr, r;
    fcomplex *response;

    /* Check that the arguments are OK */

    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_r_response().\n\n", roffset);
        exit(-1);
    }
    if (numbetween < 1 || numbetween >= 20000) {
        printf("\n  numbetween = %d (out of bounds) in gen_r_response().\n\n",
               numbetween);
        exit(-1);
    }
    if (numkern < numbetween) {
        printf("\n  numkern = %d (out of bounds) in gen_r_response().\n\n", numkern);
        exit(-1);
    }
    if ((numkern % (2 * numbetween)) != 0) {
        printf("\n  numkern %% (2 * numbetween) != 0 in gen_r_response().\n\n");
        exit(-1);
    }

    /* Prep the recursion */

    response = gen_cvect(numkern);
    startr = PI * (numkern / (double) (2 * numbetween) + roffset);
    delta = -PI / numbetween;
    tmp = sin(0.5 * delta);
    alpha = -2.0 * tmp * tmp;
    beta = sin(delta);
    c = cos(startr);
    s = sin(startr);

    /* Generate the points */

    for (ii = 0, r = startr; ii < numkern; ii++, r += delta) {
        if (r == 0.0)
            sinc = 1.0;
        else
            sinc = s / r;
        response[ii].r = c * sinc;
        response[ii].i = s * sinc;
        c = alpha * (tmp = c) - beta * s + c;
        s = alpha * s + beta * tmp + s;
    }

    /* Correct for divide by zero when the roffset is close to zero */

    // printf("[lqq] numkern:%d\n", numkern);
    if (roffset < 1E-3) {
        response[numkern / 2].r = 1 - 6.579736267392905746 * (tmp =
                                                              roffset * roffset);
        response[numkern / 2].i = roffset * (PI - 10.335425560099940058 * tmp);
    }
    return response;
}


#define SQRT2         1.4142135623730950488016887242096980785696718753769
#define PI            3.1415926535897932384626433832795028841971693993751
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876

static double sn[6] = {
    -2.99181919401019853726E3,
    7.08840045257738576863E5,
    -6.29741486205862506537E7,
    2.54890880573376359104E9,
    -4.42979518059697779103E10,
    3.18016297876567817986E11,
};

static double sd[6] = {
    /* 1.00000000000000000000E0, */
        2.81376268889994315696E2,
        4.55847810806532581675E4,
        5.17343888770096400730E6,
        4.19320245898111231129E8,
        2.24411795645340920940E10,
        6.07366389490084639049E11,
};

/* C(x) for small x */
static double cn[6] = {
    -4.98843114573573548651E-8,
    9.50428062829859605134E-6,
    -6.45191435683965050962E-4,
    1.88843319396703850064E-2,
    -2.05525900955013891793E-1,
    9.99999999999999998822E-1,
};

static double cd[7] = {
    3.99982968972495980367E-12,
    9.15439215774657478799E-10,
    1.25001862479598821474E-7,
    1.22262789024179030997E-5,
    8.68029542941784300606E-4,
    4.12142090722199792936E-2,
    1.00000000000000000118E0,
};

/* Auxiliary function f(x) */
static double fn[10] = {
    4.21543555043677546506E-1,
    1.43407919780758885261E-1,
    1.15220955073585758835E-2,
    3.45017939782574027900E-4,
    4.63613749287867322088E-6,
    3.05568983790257605827E-8,
    1.02304514164907233465E-10,
    1.72010743268161828879E-13,
    1.34283276233062758925E-16,
    3.76329711269987889006E-20,
};

static double fd[10] = {
    /*  1.00000000000000000000E0, */
        7.51586398353378947175E-1,
        1.16888925859191382142E-1,
        6.44051526508858611005E-3,
        1.55934409164153020873E-4,
        1.84627567348930545870E-6,
        1.12699224763999035261E-8,
        3.60140029589371370404E-11,
        5.88754533621578410010E-14,
        4.52001434074129701496E-17,
        1.25443237090011264384E-20,
};

/* Auxiliary function g(x) */
static double gn[11] = {
    5.04442073643383265887E-1,
    1.97102833525523411709E-1,
    1.87648584092575249293E-2,
    6.84079380915393090172E-4,
    1.15138826111884280931E-5,
    9.82852443688422223854E-8,
    4.45344415861750144738E-10,
    1.08268041139020870318E-12,
    1.37555460633261799868E-15,
    8.36354435630677421531E-19,
    1.86958710162783235106E-22,
};

static double gd[11] = {
    /*  1.00000000000000000000E0, */
        1.47495759925128324529E0,
        3.37748989120019970451E-1,
        2.53603741420338795122E-2,
        8.14679107184306179049E-4,
        1.27545075667729118702E-5,
        1.04314589657571990585E-7,
        4.60680728146520428211E-10,
        1.10273215066240270757E-12,
        1.38796531259578871258E-15,
        8.39158816283118707363E-19,
        1.86958710162783236342E-22,
};


void fresnl_avx512(__m512d xxa_avx, __m512d* ssa, __m512d* cca, __mmask8 mask0)
{
    __m512d cc_avx, ss_avx, cc1_avx, ss1_avx, cc2_avx, ss2_avx, x_avx, x2_avx;
    cc_avx = _mm512_set1_pd(0.5);
    ss_avx = _mm512_set1_pd(0.5);
    x_avx = _mm512_abs_pd(xxa_avx);
    x2_avx = _mm512_mul_pd(x_avx, x_avx);

    // fresnl_avx512_f1(cc, ss, x, x2, max_i);
    __m512d t_avx, cn_avx, sn_avx, cd_avx, sd_avx;
	__mmask8 mask1 = _mm512_cmplt_pd_mask(x2_avx, _mm512_set1_pd(2.5625));
	mask1 &= mask0;
    if (mask1 != 0)
    {
        //t[i] = x2[i] * x2[i];
        t_avx = _mm512_mul_pd(x2_avx, x2_avx);
        //ss[i] = x[i] * x2[i] * polevl(t[i], sn, 5) / p1evl(t[i], sd, 6);
        //cc[i] = x[i] * polevl(t[i], cn, 5) / polevl(t[i], cd, 6);
        sn_avx = _mm512_set1_pd(sn[0]);
        sd_avx = _mm512_add_pd(t_avx, _mm512_set1_pd(sd[0]));
        cn_avx = _mm512_set1_pd(cn[0]);
        cd_avx = _mm512_set1_pd(cd[0]);
        sn_avx = _mm512_fmadd_pd(sn_avx, t_avx, _mm512_set1_pd(sn[1]));
        sd_avx = _mm512_fmadd_pd(sd_avx, t_avx, _mm512_set1_pd(sd[1]));
        cn_avx = _mm512_fmadd_pd(cn_avx, t_avx, _mm512_set1_pd(cn[1]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[1]));
        sn_avx = _mm512_fmadd_pd(sn_avx, t_avx, _mm512_set1_pd(sn[2]));
        sd_avx = _mm512_fmadd_pd(sd_avx, t_avx, _mm512_set1_pd(sd[2]));
        cn_avx = _mm512_fmadd_pd(cn_avx, t_avx, _mm512_set1_pd(cn[2]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[2]));
        sn_avx = _mm512_fmadd_pd(sn_avx, t_avx, _mm512_set1_pd(sn[3]));
        sd_avx = _mm512_fmadd_pd(sd_avx, t_avx, _mm512_set1_pd(sd[3]));
        cn_avx = _mm512_fmadd_pd(cn_avx, t_avx, _mm512_set1_pd(cn[3]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[3]));
        sn_avx = _mm512_fmadd_pd(sn_avx, t_avx, _mm512_set1_pd(sn[4]));
        sd_avx = _mm512_fmadd_pd(sd_avx, t_avx, _mm512_set1_pd(sd[4]));
        cn_avx = _mm512_fmadd_pd(cn_avx, t_avx, _mm512_set1_pd(cn[4]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[4]));
        sn_avx = _mm512_fmadd_pd(sn_avx, t_avx, _mm512_set1_pd(sn[5]));
        sd_avx = _mm512_fmadd_pd(sd_avx, t_avx, _mm512_set1_pd(sd[5]));
        cn_avx = _mm512_fmadd_pd(cn_avx, t_avx, _mm512_set1_pd(cn[5]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[5]));
        cd_avx = _mm512_fmadd_pd(cd_avx, t_avx, _mm512_set1_pd(cd[6]));
        cc1_avx = _mm512_div_pd(_mm512_mul_pd(x_avx, cn_avx), cd_avx);
        ss1_avx = _mm512_div_pd(_mm512_mul_pd(_mm512_mul_pd(x_avx, x2_avx), sn_avx), sd_avx);       
        ss_avx = _mm512_mask_mov_pd(ss_avx, mask1, ss1_avx);
        cc_avx = _mm512_mask_mov_pd(cc_avx, mask1, cc1_avx);
    }

    // fresnl_avx512_f2(xxa, cc, ss, x, x2, max_i);
    double sin_xx[8], cos_xx[8];
    __m512d f_avx, g_avx, c_avx, s_avx, u_avx, fn_avx, gn_avx, fd_avx, gd_avx;
	__mmask8 mask2 = _mm512_cmple_pd_mask(x_avx, _mm512_set1_pd(36974.0));
	mask2 &= mask0;
    mask2 &= ~mask1;
    if (mask2 != 0)
    {
        //t[i] = PI * x2[i];
        t_avx = _mm512_mul_pd(_mm512_set1_pd(PI), x2_avx);
        //u[i] = 1.0 / (t[i] * t[i]);
        //t[i] = 1.0 / t[i];
        t_avx = _mm512_div_pd(_mm512_set1_pd(1.0), t_avx);
		u_avx = _mm512_mul_pd(t_avx, t_avx);
        //f[i] = 1.0 - u[i] * polevl(u[i], fn, 9) / p1evl(u[i], fd, 10);
        //g[i] = t[i] * polevl(u[i], gn, 10) / p1evl(u[i], gd, 11);
        fn_avx = _mm512_set1_pd(fn[0]);
        fd_avx = _mm512_add_pd(u_avx, _mm512_set1_pd(fd[0]));
        gn_avx = _mm512_set1_pd(gn[0]);
        gd_avx = _mm512_add_pd(u_avx, _mm512_set1_pd(gd[0]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[1]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[1]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[1]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[1]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[2]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[2]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[2]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[2]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[3]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[3]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[3]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[3]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[4]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[4]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[4]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[4]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[5]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[5]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[5]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[5]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[6]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[6]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[6]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[6]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[7]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[7]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[7]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[7]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[8]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[8]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[8]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[8]));
        fn_avx = _mm512_fmadd_pd(fn_avx, u_avx, _mm512_set1_pd(fn[9]));
        fd_avx = _mm512_fmadd_pd(fd_avx, u_avx, _mm512_set1_pd(fd[9]));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[9]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[9]));
        f_avx = _mm512_sub_pd(_mm512_set1_pd(1.0), _mm512_div_pd(_mm512_mul_pd(u_avx, fn_avx), fd_avx));
        gn_avx = _mm512_fmadd_pd(gn_avx, u_avx, _mm512_set1_pd(gn[10]));
        gd_avx = _mm512_fmadd_pd(gd_avx, u_avx, _mm512_set1_pd(gd[10]));
        g_avx = _mm512_div_pd(_mm512_mul_pd(t_avx, gn_avx), gd_avx);
        //c[i] = cos_kx2[i];
        //s[i] = sin_kx2[i];
		__m512d kx2_avx = _mm512_mul_pd(_mm512_set1_pd(PIBYTWO), x2_avx);
		s_avx = _mm512_sincos_pd(&c_avx, kx2_avx);
        //cc[i] = 0.5 + (f[i] * s[i] - g[i] * c[i]) / (PI * x[i]);
        //ss[i] = 0.5 - (f[i] * c[i] + g[i] * s[i]) / (PI * x[i]);
		t_avx = _mm512_div_pd(_mm512_set1_pd(1.0 / PI), x_avx);
        cc2_avx = _mm512_add_pd(_mm512_set1_pd(0.5), _mm512_mul_pd(_mm512_fmsub_pd(f_avx, s_avx, _mm512_mul_pd(g_avx, c_avx)), t_avx));
        ss2_avx = _mm512_sub_pd(_mm512_set1_pd(0.5), _mm512_mul_pd(_mm512_fmadd_pd(f_avx, c_avx, _mm512_mul_pd(g_avx, s_avx)), t_avx));
        cc_avx = _mm512_mask_mov_pd(cc_avx, mask2, cc2_avx);
        ss_avx = _mm512_mask_mov_pd(ss_avx, mask2, ss2_avx);
    }
    __mmask8 sign_mask = _mm512_cmplt_pd_mask(xxa_avx, _mm512_set1_pd(0.0));
    const unsigned long long NEG_MASK = 0x8000000000000000;
    __m512d neg_mask_avx = _mm512_set1_pd(*(double*)&NEG_MASK);
    cc_avx = _mm512_mask_mov_pd(cc_avx, sign_mask, _mm512_xor_pd(cc_avx, neg_mask_avx));
    ss_avx = _mm512_mask_mov_pd(ss_avx, sign_mask, _mm512_xor_pd(ss_avx, neg_mask_avx));
    *cca = cc_avx;
    *ssa = ss_avx;
}

fcomplex *gen_z_response(double roffset, int numbetween, double z, int numkern)
/*  Generate the response function for Fourier f-dot interpolation.  */
/*  Arguments:                                                       */
/*    'roffset' is the offset in Fourier bins for the full response  */
/*       (i.e. At this point, the response would equal 1.0)          */
/*    'numbetween' is the number of points to interpolate between    */
/*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
/*    'z' is the Fourier Frequency derivative (# of bins the signal  */
/*       smears over during the observation).                        */
/*    'numkern' is the number of complex points that the kernel will */
/*       contain.                                                    */
{
	int ii, i, numkernby2;
	double useless;
    fcomplex *response;

	/* Check that the arguments are OK */

	if (roffset < 0.0 || roffset >= 1.0) {
		printf("\n  roffset = %f (out of bounds) in gen_z_response().\n\n", roffset);
		exit(-1);
	}
	if (numbetween < 1 || numbetween >= 20000) {
		printf("\n  numbetween = %d (out of bounds) in gen_z_response().\n\n",
			numbetween);
		exit(-1);
	}
	if (numkern < numbetween) {
		printf("\n  numkern = %d (out of bounds) in gen_z_response().\n\n", numkern);
		exit(-1);
	}
	if ((numkern % (2 * numbetween)) != 0) {
		printf("\n  numkern %% (2 * numbetween) != 0 in gen_z_response().\n\n");
		exit(-1);
	}

	/* If z~=0 use the normal Fourier interpolation kernel */
	double absz = fabs(z);
	if (absz < 1E-4) {
		response = gen_r_response(roffset, numbetween, numkern);
		return response;
	}

    response = gen_cvect(numkern);

	/* Begin the calculations */

	double startr = roffset - (0.5 * z);
	startr += numkern / (double)(2 * numbetween);

    const unsigned long long NEG_MASK = 0x8000000000000000;
    __m512d neg_mask_avx = _mm512_set1_pd(*(double*)&NEG_MASK);
    __m512d signz_mask_avx = (z < 0.0) ? neg_mask_avx : _mm512_set1_pd(0.0);
    double zd = ((z < 0.0) ? -SQRT2 : SQRT2) / sqrt(absz);
    double cons = zd / 2.0;
    double pibyz = PI / z;

    __m128d r_avx = _mm_set_sd(startr);
    __m128d delta_avx = _mm_set_sd(-1.0 / numbetween);
    for (ii = 0; ii < numkern; ii += 8)
    {
        int max_i = numkern - ii;
        if (max_i > 8)
            max_i = 8;
        __mmask8 mask0 = (__mmask8)((1 << max_i) - 1);

        double rr[8];
        for (i = 0; i < max_i; i++)
        {
            _mm_store_sd(&rr[i], r_avx);
            r_avx = _mm_add_sd(r_avx, delta_avx);
        }
        __m512d xx_avx = _mm512_maskz_loadu_pd(mask0, rr);
        __m512d yy_avx = _mm512_mul_pd(xx_avx, _mm512_set1_pd(zd));
        __m512d zz_avx = _mm512_fmadd_pd(_mm512_set1_pd(z), _mm512_set1_pd(zd), yy_avx);
		__m512d x2_avx = _mm512_mul_pd(xx_avx, xx_avx);
		__m512d kx2_avx = _mm512_mul_pd(_mm512_set1_pd(pibyz), x2_avx);
		__m512d sin_xx_avx, cos_xx_avx;
		sin_xx_avx = _mm512_sincos_pd(&cos_xx_avx, kx2_avx);

        __m512d fressy_avx, frescy_avx, fressz_avx, frescz_avx;
        fresnl_avx512(yy_avx, &fressy_avx, &frescy_avx, mask0);
        fresnl_avx512(zz_avx, &fressz_avx, &frescz_avx, mask0);
        __m512d tmprl_avx = _mm512_xor_pd(_mm512_sub_pd(frescz_avx, frescy_avx), signz_mask_avx);
        __m512d tmpim_avx = _mm512_sub_pd(fressy_avx, fressz_avx);
        __m512d response_r_avx = _mm512_mul_pd(_mm512_fmsub_pd(tmprl_avx, cos_xx_avx, _mm512_mul_pd(tmpim_avx, sin_xx_avx)), _mm512_set1_pd(cons));
        __m512d response_i_avx = _mm512_mul_pd(_mm512_fmadd_pd(tmprl_avx, sin_xx_avx, _mm512_mul_pd(tmpim_avx, cos_xx_avx)), _mm512_set1_pd(-cons));
        double response_r[8], response_i[8];
        _mm512_mask_storeu_pd(response_r, mask0, response_r_avx);
        _mm512_mask_storeu_pd(response_i, mask0, response_i_avx);
        for (i = 0; i < max_i; i++)
        {
            response[ii + i].r = response_r[i];
            response[ii + i].i = response_i[i];
        }
    }

	/* Correct for divide by zero when the roffset and z is close to zero */
	double startroffset = (startr < 0) ? 1.0 + modf(startr, &useless) : modf(startr, &useless);

	if (startroffset < 1E-3 && absz < 1E-3) {
		double zz = z * z;
		double xx = startroffset * startroffset;
		numkernby2 = numkern / 2;
		response[numkernby2].r = 1.0 - 0.16449340668482264365 * zz;
		response[numkernby2].i = -0.5235987755982988731 * z;
		response[numkernby2].r += startroffset * 1.6449340668482264365 * z;
		response[numkernby2].i += startroffset * (PI - 0.5167712780049970029 * zz);
		response[numkernby2].r += xx * (-6.579736267392905746
			+ 0.9277056288952613070 * zz);
		response[numkernby2].i += xx * (3.1006276680299820175 * z);
	}
    return response;
}



// Tests conducted checking the fractional deviation of the amplitudes
// of the w-response calculation using different num_pts_wdat,
// compared to 262144.  roffset=[0,1], z=[-200,200], w=[-1000,1000]
//
// NUM_PTS_WDAT  MinFracDev   MedFracDev  MaxFracDev
//   131072      1.5983e-05   6.4267e-05   0.002060
//    65536      5.1875e-05   0.00021747   0.005147
//    32768      0.00012699   0.00051079   0.012568
//    16384      0.00027375   0.00112215   0.026279
//     8192      0.00054102   0.00221496   0.053507
//     4096      0.00104040   0.00410371   0.101785
//     2048      0.00244757   0.00875644   0.224530
//     1024      0.00427585   0.01669957   0.497524

fcomplex *gen_w_response(double roffset, int numbetween, double z,
                         double w, int numkern)
  /*  Generate the response for Fourier f, f-dot, f-dotdot interp.     */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
  /*  This version uses zero-padding to get the "numbetween"           */
{
    int ii, fbar, num_pts_wdat;
    float *data;
    double amp, f, fd, fdd, dt, t, phase, dfbar;
    fcomplex *response;

    /* Check that the arguments are OK */
    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_w_response().\n\n", roffset);
        exit(-1);
    }
    if (numbetween < 1 || numbetween >= 20000) {
        printf("\n  numbetween = %d (out of bounds) in gen_w_response().\n\n",
               numbetween);
        exit(-1);
    }
    if (numkern < numbetween) {
        printf("\n  numkern = %d (out of bounds) in gen_w_response().\n\n", numkern);
        exit(-1);
    }
    if ((numkern % (2 * numbetween)) != 0) {
        printf("\n  numkern %% (2 * numbetween) != 0 in gen_w_response().\n\n");
        exit(-1);
    }

    /* If w~=0 use the normal F-dot Fourier interpolation kernel */
    if (fabs(w) < 1E-4) {
        response = gen_z_response(roffset, numbetween, z, numkern);
        return response;
    }

    /* Cheeose num_pts_wdat so that there is plenty of Freq range */
    /* outside of the RZW response. */
    num_pts_wdat = next2_to_n(6 * w_resp_halfwidth(z, w, LOWACC) +
                              200 + numkern / numbetween);

    /* Otherwise initialize some data */
    dt = 1.0 / (double) num_pts_wdat;
    amp = 2.0 * dt;
    fbar = num_pts_wdat / 4;  // num_pts_wdat / 4 is average freq
    dfbar = (double) fbar + roffset;
    // r_o = rbar - zbar/2 + w/12  where _o is initial and bar is average
    // z_o = zbar - w/2
    f = dfbar - 0.5 * z + w / 12.0;     //  This shifts the initial f appropriately
    fd = (z - 0.5 * w) / 2.0;   // z - w/2 is the initial z value
    fdd = w / 6.0;

    /* Generate the data set.  Use zero-padding to do the interpolation. */
    data = gen_fvect(num_pts_wdat * numbetween);
    for (ii = 0; ii < num_pts_wdat * numbetween; ii++)
        data[ii] = 0.0;
    for (ii = 0; ii < num_pts_wdat; ii++) {
        t = ii * dt;
        phase = TWOPI * (t * (t * (t * fdd + fd) + f));
        data[ii] = amp * cos(phase);
    }

    /* FFT the data */
    realfft(data, num_pts_wdat * numbetween, -1);

    /* Generate the final response */
    response = gen_cvect(numkern);

    /* Chop off the contaminated ends and/or the extra data */
    memcpy(response, data + 2 * (fbar * numbetween - numkern / 2),
           sizeof(fcomplex) * numkern);

    /* cleanup */
    vect_free(data);
    return response;
}


fcomplex *gen_w_response2(double roffset, int numbetween, double z,
                          double w, int numkern)
  /*  Generate the response for Fourier f, f-dot, f-dotdot interp.     */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
  /*  This version uses Fourier interpolation to get the "numbetween"  */
{
    int ii, fbar, num_pts_wdat;
    float *data;
    double amp, f, fd, fdd, dt, t, phase, dfbar;
    fcomplex *response = NULL;
    static int firsttime = 1, old_numkern = 0, old_numbetween = 1, old_fftlen = 0;
    static fcomplex *kernelarray = NULL;

    /* Check that the arguments are OK */
    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_w_response().\n\n", roffset);
        exit(-1);
    }
    if (numbetween < 1 || numbetween >= 20000) {
        printf("\n  numbetween = %d (out of bounds) in gen_w_response().\n\n",
               numbetween);
        exit(-1);
    }
    if (numkern < numbetween) {
        printf("\n  numkern = %d (out of bounds) in gen_w_response().\n\n", numkern);
        exit(-1);
    }
    if ((numkern % (2 * numbetween)) != 0) {
        printf("\n  numkern %% (2 * numbetween) != 0 in gen_w_response().\n\n");
        exit(-1);
    }

    /* If w~=0 use the normal F-dot Fourier interpolation kernel */
    if (fabs(w) < 1E-4) {
        response = gen_z_response(roffset, numbetween, z, numkern);
        return response;
    }

    /* Cheeose num_pts_wdat so that there is plenty of Freq range */
    /* outside of the RZW response. */
    num_pts_wdat = next2_to_n(6 * w_resp_halfwidth(z, w, LOWACC) +
                              200 + numkern / numbetween);

    /* Otherwise initialize some data */
    dt = 1.0 / (double) num_pts_wdat;
    amp = 2.0 * dt;
    fbar = num_pts_wdat / 4; // num_pts_wdat / 4 is average freq
    dfbar = (double) fbar + roffset;
    // r_o = rbar - zbar/2 + w/12  where _o is initial and bar is average
    // z_o = zbar - w/2
    f = dfbar - 0.5 * z + w / 12.0; //  This shifts the initial f appropriately
    fd = 0.5 * (z - 0.5 * w); // z - w/2 is the initial z value
    fdd = w / 6.0;

    /* Generate the data set. */
    data = gen_fvect(num_pts_wdat);
    for (ii = 0; ii < num_pts_wdat; ii++) {
        t = ii * dt;
        phase = TWOPI * (t * (t * (t * fdd + fd) + f));
        data[ii] = amp * cos(phase);
    }

    /* FFT the data */
    realfft(data, num_pts_wdat, -1);

    /* The following block saves us from having to re-compute */
    /* the Fourier interpolation kernels if 'numkern' is the  */
    /* same length as on prior calls.                         */

    if (numbetween > 1) {
        int numintkern, fftlen, beginbin;
        fcomplex *tmpresponse, *rresp, *dataarray;

        beginbin = (int)(fbar) - numkern / (2 * numbetween);
        numintkern = 2 * numbetween * r_resp_halfwidth(HIGHACC);
        fftlen = next2_to_n(numkern + numintkern);
        //if (old_fftlen > fftlen)
        //    fftlen = old_fftlen;
        if (firsttime ||
            old_numkern != numkern ||
            old_numbetween != numbetween ||
            old_fftlen != fftlen) {

            /* Free the old kernelarray if one exists */
            if (!firsttime)
                vect_free(kernelarray);

            /* Generate the interpolating kernel array */
            kernelarray = gen_cvect(fftlen);
            rresp = gen_r_response(0.0, numbetween, numintkern);
            place_complex_kernel(rresp, numintkern, kernelarray, fftlen);
            vect_free(rresp);

            /* FFT the kernel array */
            COMPLEXFFT(kernelarray, fftlen, -1);

            /* Set our new static variables */
            old_numkern = numkern;
            old_numbetween = numbetween;
            old_fftlen = fftlen;
            firsttime = 0;
        }

        /* Generate the data array */
        dataarray = gen_cvect(fftlen);
        spread_no_pad(((fcomplex *) (data + 2 * beginbin)),
                      numkern / numbetween, dataarray, fftlen, numbetween);

        /* Generate the final response */
        response = gen_cvect(numkern);
        tmpresponse = complex_corr_conv(dataarray, kernelarray,
                                        fftlen, FFTD, CORR);

        /* Chop off the extra data */
        memcpy(response, tmpresponse, sizeof(fcomplex) * numkern);
        vect_free(tmpresponse);
        vect_free(dataarray);

    } else {

        /* Generate the final response */
        response = gen_cvect(numkern);

        /* Chop off the contaminated ends and/or the extra data */
        memcpy(response, data + 2 * (fbar * numbetween - numkern / 2),
               sizeof(fcomplex) * numkern);
    }

    /* cleanup */
    vect_free(data);
    return response;
}


fcomplex *gen_bin_response(double roffset, int numbetween, double ppsr,
                           double T, orbitparams * orbit, int numkern)
  /*  Generate the Fourier response function for a sinusoidal PSR      */
  /*  signal from a binary orbit.                                      */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{

    int fftlen, ii, beginbin, numintkern, index, datar;
    int numdata = MIN_NUMDATA;
    int numorbpts = MIN_NUMORBPTS;      /* This should be a power-of-two + 1 */
    float *data;
    double *phi = NULL, startE;
    double amp, f, dt, dtb, t, tp, fpart, ipart, dtemp;
    static int old_numbetween = 0, old_numkern = 0, old_fftlen = 0, firsttime = 1;
    static fcomplex *kernelarray = NULL;
    fcomplex *response, *tmpresponse, *rresp, *dataarray;
    orbitparams orb;

    /* Check that the arguments are OK */

    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_bin_response().\n\n",
               roffset);
        exit(-1);
    }
    if (numbetween < 1 || numbetween >= 20000) {
        printf("\n  numbetween = %d (out of bounds) in gen_bin_response().\n\n",
               numbetween);
        exit(-1);
    }
    if (numkern < numbetween) {
        printf("\n  numkern = %d (out of bounds) in gen_bin_response().\n\n",
               numkern);
        exit(-1);
    }
    if ((numkern % (2 * numbetween)) != 0) {
        printf("\n  numkern %% (2 * numbetween) != 0 in gen_bin_response().\n\n");
        exit(-1);
    }

    /* Initialize some data */

    datar = numdata / 4;
    if (numkern > datar) {
        numdata = next2_to_n(numkern * 4);
        datar = numdata / 4;
    }
    dt = 1.0 / numdata;
    orb.p = orbit->p / T;
    if (orb.p < 1.0)
        numorbpts = next2_to_n(MIN_NUMORBPTS / orb.p) + 1;
    dtb = 1.0 / (numorbpts - 1);
    amp = 2.0 * dt;
    f = TWOPI * datar;
    orb.x = orbit->x / (ppsr * datar);
    orb.e = orbit->e;
    orb.w = orbit->w;
    orb.t = orbit->t / T;

    /* Generate the orbit */

    startE = keplers_eqn(orb.t, orb.p, orb.e, 1.0E-15);
    phi = dorbint(startE, numorbpts, dtb, &orb);
    E_to_phib(phi, numorbpts, &orb);

    /* Generate the data set */

    data = gen_fvect(numdata);
    for (ii = 0; ii < numdata; ii++) {
        t = ii * dt;
        /* The following 4 lines simply linearly interpolate   */
        /* the orbital solution and add the delay to the time. */
        fpart = modf(t / dtb, &ipart);
        index = (int) (ipart + DBLCORRECT);
        dtemp = phi[index];
        tp = t - ((fpart * (phi[index + 1] - dtemp)) + dtemp);
        data[ii] = amp * cos(f * tp);
    }
    vect_free(phi);

    /* FFT the data */

    realfft(data, numdata, -1);
    beginbin = numdata / 4 - numkern / (2 * numbetween);

    /* The following block saves us from having to re-compute */
    /* the Fourier interpolation kernels if 'numkern' is the  */
    /* same length as on prior calls.                         */

    if (numbetween > 1) {
        numintkern = 2 * numbetween * r_resp_halfwidth(HIGHACC);
        fftlen = next2_to_n(numkern + numintkern);
        if (fftlen > numdata) {
            printf("WARNING:  fftlen > numdata in gen_bin_response().\n");
        }
        if (firsttime ||
            old_numkern != numkern ||
            old_numbetween != numbetween || old_fftlen != fftlen) {

            /* Generate an interpolation kernel for the data */

            rresp = gen_r_response(0.0, numbetween, numintkern);

            /* Free the old kernelarray if one exists */

            if (!firsttime)
                vect_free(kernelarray);

            /* Generate the interpolating kernel array */

            kernelarray = gen_cvect(fftlen);
            place_complex_kernel(rresp, numintkern, kernelarray, fftlen);
            vect_free(rresp);

            /* FFT the kernel array */

            COMPLEXFFT(kernelarray, fftlen, -1);

            /* Set our new static variables */

            old_numkern = numkern;
            old_numbetween = numbetween;
            old_fftlen = fftlen;
            firsttime = 0;
        }

        /* Generate the data array */

        dataarray = gen_cvect(fftlen);
        if (fftlen / numbetween >= numdata - beginbin) {
            printf("WARNING:  fftlen too large in gen_bin_response().\n");
        }
        spread_no_pad(((fcomplex *) data) + beginbin, numkern / numbetween,
                      dataarray, fftlen, numbetween);

        /* Generate the final response */

        response = gen_cvect(numkern);
        tmpresponse = complex_corr_conv(dataarray, kernelarray, fftlen, FFTD, CORR);

        /* Chop off the extra data */

        memcpy(response, tmpresponse, sizeof(fcomplex) * numkern);
        vect_free(tmpresponse);
        vect_free(dataarray);

    } else {

        /* Generate the final response */

        response = gen_cvect(numkern);
        memcpy(response, ((fcomplex *) data) + beginbin, sizeof(fcomplex) * numkern);
    }
    vect_free(data);
    return response;
}
