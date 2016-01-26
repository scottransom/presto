#define COPYSIGN(a,b) ((b) < 0.0 ? (-fabs(a)) : (fabs(a)))

/*  NOTE:  See presto.h for function definitions. */
#include "presto.h"

int nice_output_1(char *output, double val, double err, int len)
/* Generates a string in "output" of length len with "val" rounded  */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 1 decimal place.     */
/* Note:  len should be ~ 20 to show full double precision          */
/*   if the base 10 exponent of the error needs to be shown.        */
/*   If len == 0, left-justified minimum length string is returned. */
/*   If len > 0, the string returned has is right justified.        */
{
    int nint, nfrac, totprec;
    int errexp, errval, outexp;
    double rndval, outmant;
    char temp[50];

    sprintf(temp, "There is a problem with 'nice_output()'.\n");
    if (fabs(err) == 0.0) {
        errexp = 0;
    } else {
        errexp = (int) floor(log10(fabs(err)));
    }

    /* 1 digit error value:  */

    errval = (int) floor(fabs(err) * pow(10.0, (double) (-errexp)) +
                         DBLCORRECT + 0.5);
    if (errval == 10) {
        errval = 1;
        errexp++;
    }
    /* val rounded to the appropriate decimal place due to err: */

    rndval = pow(10.0, (double) errexp) *
        floor(val * pow(10.0, (double) (-errexp)) + 0.5);

    /* Space needed for integer part: */

    if (fabs(val) == 0.0) {
        nint = 1;
    } else {
        nint = (int) (ceil(log10(fabs(val))));
        if (nint == 0)
            nint++;
    }

    /* Space needed for fractional part: */

    nfrac = -errexp;

    /* Total number of digits of precision in output value: */

    totprec = nint + nfrac;

    /* Base 10 exponent of output value: */

    if (fabs(rndval) == 0.0) {
        outexp = 0;
    } else {
        outexp = (int) floor(log10(fabs(rndval)));
    }

    /* Unsigned base 10 mantissa of output value: */

    outmant = rndval * pow(10.0, (double) (-outexp));
    if (fabs(1.0 - outmant) < DBLCORRECT || fabs(-1.0 - outmant) < DBLCORRECT)
        totprec++;

    /* Use scientific notation:  */

    if ((outexp >= 0 && errexp > 0) && outexp > errexp)
        sprintf(temp, "% .*f(%d)x10^%d", totprec - 1,
                COPYSIGN(outmant, rndval), errval, outexp);

    /* Use scientific notation but with integer mantissa */

    else if ((outexp >= 0 && errexp > 0) && outexp == errexp)
        sprintf(temp, "% d(%d)x10^%d",
                (int) (COPYSIGN(outmant, rndval)), errval, outexp);

    /* Use scientific notation for real small numbers: */

    else if (outexp < -4 && outexp >= errexp)
        sprintf(temp, "% .*f(%d)x10^%d", totprec - 1,
                COPYSIGN(outmant, rndval), errval, outexp);

    /* Use scientific notation but with integer mantissa */

    else if (outexp < errexp && errexp != 0)
        sprintf(temp, "% d(%d)x10^%d",
                (int) (COPYSIGN(outmant, rndval) + DBLCORRECT), errval, errexp);

    /* Use regular notation: */

    else if (nfrac == 0 && fabs(rndval) < 1.0e-15)
        sprintf(temp, "% d(%d)", (int) fabs(rndval), errval);
    else if (fabs(rndval) <= DBLCORRECT && errexp < -5)
        sprintf(temp, "0.0(%d)x10^%d", errval, errexp + 1);
    else
        sprintf(temp, "% .*f(%d)", nfrac, rndval, errval);

    if (len == 0) {             /* Left-justify  */
        sprintf(output, "%s", temp);
    } else {                    /* Right-justify with a length of len */
        sprintf(output, "%*s", len, temp);
    }
    return strlen(output);
}


int nice_output_2(char *output, double val, double err, int len)
/* Generates a string in "output" of length len with "val" rounded  */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 2 decimal places.    */
/* Note:  len should be ~ 20 to show full double precision          */
/*   if the base 10 exponent of the error needs to be shown.        */
/*   If len == 0, left-justified minimum length string is returned. */
/*   If len > 0, the string returned has is right justified.        */
{
    int errexp, errval, outexp, valexp, tmp;
    double rndval, outmant;
    char temp[50];

    sprintf(temp, "There is a problem with 'nice_output2()'.\n");
    if (fabs(err) == 0.0) {
        errexp = 0;
    } else {
        errexp = (int) floor(log10(fabs(err)));
    }

    /* 2 digit error value: */

    errval = (int) floor(fabs(err) * pow(10.0, (double) (-errexp + 1))
                         + DBLCORRECT + 0.5);
    if (errval == 100) {
        errval = 10;
        errexp++;
    }

    /* val rounded to the appropriate decimal place due to err: */

    rndval = pow(10.0, (double) (errexp - 1)) *
        floor(val * pow(10.0, (double) (-errexp + 1)) + 0.5);

    /* Base 10 exponent of output value: */

    if (rndval == 0.0)
        outexp = 0;
    else
        outexp = (int) floor(log10(fabs(rndval)));

    /* Base 10 exponent of original value: */

    if (val == 0.0)
        valexp = 0;
    else
        valexp = (int) floor(log10(fabs(val)));

    /* Signed base 10 mantissa of output value: */

    outmant = rndval * pow(10.0, (double) (-outexp));
    if (outexp < -4) {
        tmp = outexp - errexp + 1;
        if (tmp < 2) {
            if (tmp == 0) {
                sprintf(temp, "% .1f(%d.%d)x10^%d", outmant / 10,
                        errval / 10, errval % 10, outexp + 1);
            } else
                sprintf(temp, "% .1f(%d.%d)x10^%d", outmant,
                        errval / 10, errval % 10, outexp);
        } else
            sprintf(temp, "% .*f(%.2d)x10^%d", tmp, outmant, errval, outexp);
    } else if (errexp - valexp > 1 && (errexp != 0 && errexp != 1)) {
        sprintf(temp, "0.0(%d.%d)x10^%d", errval / 10, errval % 10, errexp);
    } else if (errexp == 0) {
        if (fabs(rndval) < 0.1)
            rndval = fabs(rndval);
        sprintf(temp, "% .1f(%d.%d)", rndval, errval / 10, errval % 10);
    } else if (errexp < 0) {
        if (fabs(rndval) * pow(10.0, (double) (-errexp + 1)) < 0.1)
            rndval = fabs(rndval);
        if (fabs(rndval) <= DBLCORRECT && errexp < -4) {
            sprintf(temp, "0.0(%d.%d)x10^%d", errval / 10, errval % 10, errexp);
        } else {
            sprintf(temp, "% .*f(%.2d)", -errexp + 1, rndval, errval);
        }
    } else if (errexp - outexp == 0) {
        if (errexp == 1) {
            sprintf(temp, "% d(%.2d)", (int) rndval, errval);
        } else
            sprintf(temp, "% .1f(%d.%d)x10^%d", outmant, errval / 10,
                    errval % 10, errexp);
    } else if (errexp - outexp == 1) {
        if (fabs(outmant) < 0.99999999999999)
            outmant = fabs(outmant);
        sprintf(temp, "% .1f(%d.%d)x10^%d", outmant / 10.0, errval / 10,
                errval % 10, errexp);
    } else if ((errexp > 1 && outexp > 1) && errexp - outexp <= -1) {
        sprintf(temp, "% .*f(%.2d)x10^%d", outexp - errexp + 1, outmant,
                errval, outexp);
    } else if (errexp == 1 && outexp > 0) {
        sprintf(temp, "% d(%.2d)", (int) rndval, errval);
    } else {
        printf("This is an undefined condition in nice_output2().\n");
        exit(1);
    }

    if (len == 0) {             /* Left-justify  */
        sprintf(output, "%s", temp);
    } else {                    /* Right-justify with a length of len */
        sprintf(output, "%*s", len, temp);
    }
    return strlen(output);
}


/*  Need to modify so that it looks in the ".inf" file for dt, N, and nph */

void print_candidate(fourierprops * cand, double dt, long N,
                     double nph, int numerrdigits)
{
    double T, T2, T3, f2, fd2, fderr2, p, pd, pdd, f, fd, fdd;
    double perr, pderr, pdderr, ferr, fderr, fdderr;
    double pownph, signph, pownpherr, rawpowerr, temp;
    int width = 19;
    char output[40], output2[40];
    int (*nice_output) (char *, double, double, int);

    if (numerrdigits == 1)
        nice_output = nice_output_1;
    else if (numerrdigits == 2)
        nice_output = nice_output_2;
    else {
        printf("print_candidate() can't handle that many digits of\n");
        printf("   precision. Use either 1 or 2.");
        exit(1);
    }
    printf("\n");
    T = dt * N;
    T2 = T * T;
    T3 = T2 * T;
    pownph = cand->rawpow / nph;
    pownpherr = sqrt(2.0 * pownph);
    signph = sqrt(2.0 * pownph - log(PI * pownph));
    rawpowerr = sqrt(2.0 * cand->rawpow);

    /*  Calculate values for "frequency" format: */

    f = cand->r / T;
    fd = cand->z / T2;
    fdd = cand->w / T3;
    ferr = cand->rerr / T;
    fderr = cand->zerr / T2;
    fdderr = cand->werr / T3;

    /*  Calculate values for "period" format: */

    p = 1.0 / f;
    f2 = f * f;
    pd = -fd / f2;
    fd2 = fd * fd;
    if (fdd == 0.0)
        pdd = 0.0;
    else
        pdd = (2.0 * fd2 / f - fdd) / f2;
    perr = ferr / f2;
    pderr = sqrt(4.0 * fd2 * ferr * ferr / f2 + (fderr2 = fderr * fderr)) / f2;
    temp = (-6.0 * fd2 / f2 + 2.0 * fd / f) * ferr;
    pdderr = sqrt(temp * temp + 16.0 * fd2 * fderr2 / f2 + fdderr * fdderr);

    /*  Now output it... */

    (*nice_output) (output, p, perr, width);
    (*nice_output) (output2, cand->pow, cand->powerr, width);
    printf("Period (s)       %s  Pow/Local Pow %s\n", output, output2);
    (*nice_output) (output, pd, pderr, width);
    printf("Pdot (s/s)       %s  Sigma (Local) %*.2f\n", output, width, cand->sig);
    (*nice_output) (output, pdd, pdderr, width);
    if (pownph > 1.0e7) {
        sprintf(output, "%12.5g", pownph);
        sprintf(output2, "%*s", width, output);
    } else {
        (*nice_output) (output2, pownph, pownpherr, width);
    }
    (*nice_output) (output, pdd, pdderr, width);
    printf("Pdotdot (s/s^2)  %s  Pow/Freq 0 Pow%s\n", output, output2);
    (*nice_output) (output, f, ferr, width);
    printf("Frequency (hz)   %s  Sigma (Freq 0)%*.2f\n", output, width, signph);
    if (cand->rawpow > 1.0e7) {
        sprintf(output, "%12.5g", cand->rawpow);
        sprintf(output2, "%*s", width, output);
    } else {
        (*nice_output) (output2, cand->rawpow, rawpowerr, width);
    }
    (*nice_output) (output, fd, fderr, width);
    printf("Fdot (hz/s)      %s  Raw Power     %s\n", output, output2);
    (*nice_output) (output, fdd, fdderr, width);
    (*nice_output) (output2, cand->phs, cand->phserr, width);
    printf("Fdotdot (hz/s^2) %s  Phase (rad)   %s\n", output, output2);
    (*nice_output) (output, cand->r, cand->rerr, width);
    (*nice_output) (output2, cand->phs * RADTODEG, cand->phserr * RADTODEG, width);
    printf("DFT Frequency    %s  Phase (deg)   %s\n", output, output2);
    (*nice_output) (output, cand->z, cand->zerr, width);
    (*nice_output) (output2, cand->cen, cand->cenerr, width);
    printf("DFT Fdot         %s  Centroid      %s\n", output, output2);
    (*nice_output) (output, cand->w, cand->werr, width);
    (*nice_output) (output2, cand->pur, cand->purerr, width);
    printf("DFT Fdotdot      %s  Purity        %s\n", output, output2);
    printf("\n");
}

void print_bin_candidate(binaryprops * cand, int numerrdigits)
/* Outputs a 2 column summary of all the properties of a fourier peak  */
{
    int width = 19;
    char output[40], output2[40];
    int (*nice_output) (char *, double, double, int);

    if (numerrdigits == 1)
        nice_output = nice_output_1;
    else if (numerrdigits == 2)
        nice_output = nice_output_2;
    else {
        printf("print_bin_candidate() can't handle that many digits of\n");
        printf("   precision. Use either 1 or 2.");
        exit(1);
    }
    printf("\n");

    /*  Output it... */

    (*nice_output) (output, cand->ppsr, cand->ppsrerr, width);
    (*nice_output) (output2, cand->pow, cand->powerr, width);
    printf("PSR Period (s)   %s  Power (Norm)  %s\n", output, output2);
    (*nice_output) (output, cand->fpsr, cand->fpsrerr, width);
    printf("PSR Freq (hz)    %s  Sigma         %*.2f\n", output, width, cand->sig);
    (*nice_output) (output, cand->rpsr, cand->rpsrerr, width);
    printf("PSR DFT Freq     %s  MiniFFT Length%*ld\n",
           output, width, cand->nfftbins);
    (*nice_output) (output, cand->pbin / SECPERDAY, cand->pbinerr / SECPERDAY,
                    width);
    printf("Bin Period (days)%s  MiniFFT Lowbin%*ld\n", output, width, cand->lowbin);
    (*nice_output) (output, cand->pbin, cand->pbinerr, width);
    (*nice_output) (output2, cand->rdetect, cand->rdetecterr, width);
    printf("Bin Period (sec) %s  MiniFFT Freq  %s\n", output, output2);
    (*nice_output) (output, cand->rbin, cand->rbinerr, width);
    (*nice_output) (output2, cand->cen, cand->cenerr, width);
    printf("Bin DFT Freq     %s  Centroid      %s\n", output, output2);
    (*nice_output) (output, cand->asinic, cand->asinicerr, width);
    (*nice_output) (output2, cand->pur, cand->purerr, width);
    printf("Bin a_sin_i/c (s)%s  Purity        %s\n", output, output2);
    (*nice_output) (output, cand->z, cand->zerr, width);
    (*nice_output) (output2, cand->phs, cand->phserr, width);
    printf("Freq Mod \"z\"     %s  Phase (rad)   %s\n", output, output2);
    printf("\n");
}


#undef COPYSIGN
