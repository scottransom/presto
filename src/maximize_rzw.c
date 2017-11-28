#include "accel.h"

double max_rzw_arr(fcomplex * data, long numdata, double rin, double zin,
                   double win, double *rout, double *zout,
                   double *wout, rderivs * derivs)
/* Return the Fourier frequency, f-dot, and fdotdot that    */
/* maximizes the power.                                     */
/* This isn't a true optimization, but just maximizing the  */
/* oversampled F/Fdot/Fdotdot volume.                       */
{
    float locpow = get_localpower3d(data, numdata, rin, zin, win);
    float maxpow = 0, pow, powargr, powargi;
    int kern_half_width, extra = 10, startbin = (int)(rin) - 1;
    int numz, numw, nextbin, fftlen, numbetween;
    // The factor beyond ACCEL_DR, ACCEL_DZ, ACCEL_DW will interpolate
    int interpfac = 4, wind = 0, zind = 0, rind = 0;
    double ifrac = 1.0/interpfac;
    fcomplex ***vol, amp;

    kern_half_width = w_resp_halfwidth(zin, win, LOWACC);
    numz = numw = 2 * interpfac + 1;
    numbetween = interpfac*ACCEL_RDR;
    fftlen = next2_to_n(numbetween * (2 * kern_half_width + extra));
    // printf("fftlen = %d\n", fftlen);
    vol = corr_rzw_vol(data, numdata, numbetween, startbin, 
                       zin-ACCEL_DZ, zin+ACCEL_DZ, numz,
                       win-ACCEL_DW, win+ACCEL_DW, numw, 
                       fftlen, LOWACC, &nextbin);
    {
        int ii, jj, kk;
        for (ii = 0; ii < numw; ii++) {
            for (jj = 0; jj < numz; jj++) {
                for (kk = 0; kk < 3*interpfac*ACCEL_RDR; kk++) {
                    amp = vol[ii][jj][kk];
                    pow = POWER(amp.r, amp.i);
                    if (pow > maxpow) {
                        maxpow = pow;
                        wind = ii;
                        zind = jj;
                        rind = kk;
                    }
                }
            }
        }
    }
    *rout = startbin + rind * ACCEL_DR * ifrac;
    *zout = zin-ACCEL_DZ + zind * ACCEL_DZ * ifrac;
    *wout = win-ACCEL_DW + wind * ACCEL_DW * ifrac;
    vect_free(vol[0][0]);
    vect_free(vol[0]);
    vect_free(vol);
    get_derivs3d(data, numdata, *rout, *zout, *wout, locpow, derivs);
    return maxpow;
}

double max_rzw_file(FILE * fftfile, double rin, double zin, double win,
                    double *rout, double *zout, double *wout, rderivs * derivs)
/* Return the Fourier frequency, f-dot, and fdotdot that    */
/* maximizes the power of the candidate in 'fftfile'.       */
{
    double maxpow, rin_int, rin_frac;
    int kern_half_width, filedatalen, extra = 10;
    long startbin;
    fcomplex *filedata;

    rin_frac = modf(rin, &rin_int);
    kern_half_width = w_resp_halfwidth(zin, win, HIGHACC);
    filedatalen = 2 * kern_half_width + extra;
    startbin = (long) rin_int - filedatalen / 2;

    filedata = read_fcomplex_file(fftfile, startbin, filedatalen);
    maxpow = max_rzw_arr(filedata, filedatalen, rin_frac + filedatalen / 2,
                         zin, win, rout, zout, wout, derivs);
    *rout += startbin;
    vect_free(filedata);
    return maxpow;
}
