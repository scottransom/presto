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
    int extra = 10, startbin = (int)(rin) - 1;
    int numr, numz, numw, nextbin, fftlen, numbetween;
    // The factor beyond ACCEL_DR, ACCEL_DZ, ACCEL_DW will interpolate
    int interpfac = 4, wind = 0, zind = 0, rind = 0;
    fcomplex ***vol, amp;
    // Use a more conservative length kernel
    int kern_half_width = w_resp_halfwidth(fabs(zin)+2, fabs(win)+20, LOWACC);

    numr = 3 * interpfac * ACCEL_RDR;
    numz = numw = 2 * interpfac + 1;

    numbetween = interpfac * ACCEL_RDR;
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
                for (kk = 0; kk < numr; kk++) {
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
    *rout = startbin + rind * ACCEL_DR / (double) interpfac;
    *zout = zin-ACCEL_DZ + zind * ACCEL_DZ / (double) interpfac;
    *wout = win-ACCEL_DW + wind * ACCEL_DW  / (double) interpfac;
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
    int filedatalen, extra = 10;
    long startbin;
    fcomplex *filedata;
    // Use a more conservative length kernel
    int kern_half_width = w_resp_halfwidth(fabs(zin)+2, fabs(win)+20, LOWACC);

    rin_frac = modf(rin, &rin_int);
    filedatalen = 2 * kern_half_width + extra;
    startbin = (long) rin_int - filedatalen / 2;

    filedata = read_fcomplex_file(fftfile, startbin, filedatalen);
    maxpow = max_rzw_arr(filedata, filedatalen, rin_frac + filedatalen / 2,
                         zin, win, rout, zout, wout, derivs);
    *rout += startbin;
    vect_free(filedata);
    return maxpow;
}


void max_rzw_arr_harmonics(fcomplex data[], long numdata,
                           int num_harmonics,
                           double rin, double zin, double win,
                           double *rout, double *zout, double *wout,
                           rderivs derivs[], double powers[])
/* Return the Fourier frequency, f-dot, and f-dotdot that       */
/* maximizes the *summed* power of the multi-harmonic candidate */
{
    float maxpow = 0, pow, powargr, powargi, locpow;
    long lobin, hhlobin = 0;
    int hind, ii, jj, kk, extra = 10, nextbin, fftlen;
    int interpfac = 4; // Factor beyond ACCEL_D[RZW] we will interpolate
    int numwid = 2; // The number of full peak widths we will search
    int numbetween = interpfac * ACCEL_RDR;
    int numr = 3 * numbetween; // Search 3 full Fourier bins around high harm
    int numz = 2 * numwid * interpfac + 1;
    int numw = numz;
    int rind = 0, zind = 0, wind = 0;
    double dr = 1.0 / (double) numbetween;

    // The summed power spectrum, initialized to zeros
    float ***powsum = gen_f3Darr(numw, numz, numr);
    for (ii = 0; ii < numw * numz * numr; ii++)
        powsum[0][0][ii] = 0.0;

    for (hind = 0; hind < num_harmonics; hind++) {
        int n = num_harmonics - hind; // harmonic number, starting from highest
        double rh = rin * n, zh = zin * n, wh = win * n;
        // Use a more conservative length kernel
        int kern_half_width = w_resp_halfwidth(fabs(zh)+2, fabs(wh)+20, LOWACC);
        fcomplex ***vol, amp;
        double rh_int, rh_frac, hfrac = n / (double) num_harmonics;

        locpow = get_localpower3d(data, numdata, rh, zh, wh);
        rh_frac = modf(rh, &rh_int);
        // Will do 1+ bins below and 1+ bins above rin
        lobin = (long) rh_int - 1;
        if (hind==0) hhlobin = lobin;
        fftlen = next2_to_n(numbetween * (2 * kern_half_width + extra));
        // Create the RZW volume for the harmonic.
        // Note that we are computing the z and w values in exact harmonic
        // ratios.  But the r values are on a power-of-two grid.
        vol = corr_rzw_vol(data, numdata, numbetween, lobin,
                           zh-numwid*hfrac*ACCEL_DZ,
                           zh+numwid*hfrac*ACCEL_DZ, numz,
                           wh-numwid*hfrac*ACCEL_DW,
                           wh+numwid*hfrac*ACCEL_DW, numw,
                           fftlen, LOWACC, &nextbin);
        // Calculate and sum the powers
        for (ii = 0; ii < numw; ii++) {
            for (jj = 0; jj < numz; jj++) {
                for (kk = 0; kk < numr; kk++) {
                    rind = (long) round(((hhlobin + kk * dr) * hfrac
                                        - lobin) * numbetween);
                    amp = vol[ii][jj][rind];
                    powsum[ii][jj][kk] += POWER(amp.r, amp.i) / locpow;
                }
            }
        }
        vect_free(vol[0][0]);
        vect_free(vol[0]);
        vect_free(vol);
    }

    // Now search the power sums for the highest value
    rind = zind = wind = 0;
    for (ii = 0; ii < numw; ii++) {
        for (jj = 0; jj < numz; jj++) {
            for (kk = 0; kk < numr; kk++) {
                pow = powsum[ii][jj][kk];
                if (pow > maxpow) {
                    maxpow = pow;
                    wind = ii;
                    zind = jj;
                    rind = kk;
                }
            }
        }
    }

    // Calculate the proper r, z, and w peak values
    {
        double nh = num_harmonics;
        *rout = (hhlobin + rind * dr) / nh;
        *zout = ((zin * nh) - numwid * ACCEL_DZ +
                 zind * ACCEL_DZ / (double) interpfac) / nh;
        *wout = ((win * nh) - numwid * ACCEL_DW +
                 wind * ACCEL_DW / (double) interpfac) / nh;
    }

    // Now calculate the derivatives at the peak
    for (ii = 0; ii < num_harmonics; ii++) {
        int hh = ii + 1;
        double rr = *rout * hh, zz = *zout * hh, ww = *wout * hh;
        locpow = get_localpower3d(data, numdata, rr, zz, ww);
        get_derivs3d(data, numdata, rr, zz, ww, locpow, &(derivs[ii]));
        powers[ii] = derivs[ii].pow;
    }
    /*
    printf("numr = %d  numz = %d  numw = %d\n", numr, numz, numw);
    printf("rind = %d  zind = %d  wind = %d\n", rind, zind, wind);
    printf("rin  = %f  zin  = %f  win  = %f\n", rin , zin , win);
    printf("rout = %f  zout = %f  wout = %f\n", *rout, *zout, *wout);
    */
    vect_free(powsum[0][0]);
    vect_free(powsum[0]);
    vect_free(powsum);
}
