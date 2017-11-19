#include "makedata.h"
#include "makeinf.h"
#include "randlib.h"

extern char bands[NUMBANDS][40];
extern char scopes[NUMSCOPES][40];
double fwhm;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    /* Basic variables */
    /* long numpad; */
    double tb = 0.0, T, phase, amp, signal, orbmaxt = 0.0;
    infodata idata;
    makedata mdata;

    /* Initialization strings */
    char datafilenm[200], infofilenm[200];
    FILE *datfile;

    /* For binary orbit integration */
    double *phib = NULL, startE, orbdt = 1.0;
    long numpoints;

    /* Signal on and off */
    double ont, offt;
    double *onoffpt;
    double (*ptype_ptr) (double val) = modsin;

    /* Looping and temp variables */
    long ct, i;

    /* Writing to files and status */
    int buffloc, oldper = 0, newper;
    float tempsig[BUFFSIZE];


    printf("\n\n");
    printf("  Pulsation Data Generation Program\n");
    printf("        by Scott M. Ransom\n");
    printf("           23 Jun, 1998\n\n");

    /* Reads input data from a .mak file if a     */
    /*   filename was given as a parameter in     */
    /*   the command line.                        */

    if (argc > 1) {

        read_mak_file(argv[1], &mdata);
        printf("\nDescription:  %s\n\n", mdata.description);

    } else {

        read_mak_input(&mdata);
        write_mak_file(&mdata);

    }

    sprintf(datafilenm, "%s.dat", mdata.basefilenm);
    sprintf(infofilenm, "%s.inf", mdata.basefilenm);

    printf("Generating data...\n");

    /* Open the data file */

    datfile = chkfopen(datafilenm, "w+b");

    /* Determine the orbital phase adjustments if needed */

    if (mdata.binary) {
        startE = keplers_eqn(mdata.orb.t, mdata.orb.p, mdata.orb.e, 1.0E-15);
        if (mdata.T > 2048)
            orbdt = 0.5;
        else
            orbdt = mdata.T / 4096.0;
        numpoints = (long) floor(mdata.T / orbdt + 0.5) + 2;
        phib = dorbint(startE, numpoints, orbdt, &mdata.orb);
        orbmaxt = orbdt * numpoints;
        E_to_phib(phib, numpoints, &mdata.orb);
    }

    /* Initialize Data Generation Loop */

    switch (mdata.pnum) {
    case 1:
        ptype_ptr = modsin;
        break;
    case 2:
        ptype_ptr = crab;
        break;
    case 3:
        ptype_ptr = spike;
        fwhm = mdata.fwhm;
        break;
    case 4:
        ptype_ptr = gauss;
        fwhm = mdata.fwhm;
        break;
    }

    printf("Calculating data points...\n\n");
    printf("Amount Complete = 0%%");
    fflush(stdout);
    amp = mdata.amp;
    mdata.phs /= 360.0;
    mdata.fd /= 2.0;
    mdata.fdd /= 6.0;
    mdata.ampmodf *= TWOPI;
    mdata.ampmodp *= DEGTORAD;
    onoffpt = mdata.onoff;
    ont = (*onoffpt++) * mdata.T;
    offt = (*onoffpt++) * mdata.T;

    /* Main data loop */

    for (ct = 0, buffloc = 1; ct < mdata.N; ct++, buffloc++) {

        /* Calculate percentage complete */

        tb = ct * mdata.dt;

        newper = (int) (tb / mdata.T * 100.0) + 1;
        if (newper > oldper) {
            printf("\rAmount Complete = %3d%%", newper);
            fflush(stdout);
            oldper = newper;
        }
        /*  Advance onoff pointers when signal turns off */

        if (tb >= offt)
            do {
                ont = (*onoffpt++) * mdata.T;
                offt = (*onoffpt++) * mdata.T;
            }
            while (tb >= offt);

        /*  Signal is on */

        if ((tb >= ont) && (tb < offt)) {

            T = tb;

            /*  Get additional time due to binary orbit */

            if (mdata.binary)
                T -= lin_interp_E(phib, tb, 0.0, orbdt, orbmaxt);

            /*  Get current pulsar phase (cyclic) */

            phase = T * (T * (T * mdata.fdd + mdata.fd) + mdata.f) + mdata.phs;

            /*  Get current pulsation amplitude */
            /*     Question:  Should this use time T or tb? */

            if (mdata.ampmod)
                amp = mdata.ampmoda * cos(mdata.ampmodf * tb + mdata.ampmodp);

            /*  Calculate the signal */

            signal = amp * (*ptype_ptr) (phase) + mdata.dc;
        }
        /*  Signal is off  */

        else
            signal = mdata.dc;

        /*  Add Poissonian noise  */

        if (mdata.noise == 1)
            signal = (float) ignpoi(signal);

        /*  Add Gaussian noise or no noise */

        else if (mdata.noisesig != 0.0) {
            signal = gennor(signal, mdata.noisesig);

            /*  Rounds if needed */

            if (mdata.roundnum)
                signal = floor(signal + 0.5);
        }
        /*  Save data in buffer */

        tempsig[buffloc - 1] = (float) signal;

        /*  Write the data when needed */

        if ((buffloc == BUFFSIZE) || (ct == mdata.N - 1)) {
            chkfwrite(&tempsig, sizeof(float), (unsigned long) buffloc, datfile);
            buffloc = 0;
        }
    }

    /*  Add dc padding if necessary  */

/*   if (mdata.next2_to_n != mdata.N) { */
/*     numpad = mdata.next2_to_n - mdata.N; */
/*     for (ct = 1, buffloc = 1; ct <= numpad; ct++, buffloc++) { */
/*       tempsig[buffloc - 1] = mdata.dc; */
/*       if ((buffloc == BUFFSIZE) || (ct == numpad)) { */
/* 	chkfwrite(&tempsig, sizeof(float), \ */
/* 		  (unsigned long) buffloc, datfile); */
/* 	buffloc = 0; */
/*       } */
/*     } */
/*   } */

    printf("\n\nData saved in binary floating point format ");
    printf("in \"%s\".\n\n", datafilenm);

    /* Prep infofile data */

    strcpy(idata.name, mdata.basefilenm);
    strcpy(idata.telescope, scopes[0]);
    strcpy(idata.analyzer, "Unknown");
    idata.N = mdata.N;
    idata.dt = mdata.dt;
    idata.numonoff = mdata.numonoff;
    i = 0;
    do {
        idata.onoff[i] = mdata.onoff[i] * mdata.N;
        idata.onoff[i + 1] = mdata.onoff[i + 1] * mdata.N;
        i += 2;
    }
    while (mdata.onoff[i - 1] != 1.0);
    idata.onoff[i - 1]--;
    if (idata.onoff[i - 2] == mdata.N)
        idata.onoff[i - 2]--;

    /* Close files and deallocate array memory */

    if (mdata.binary) {
        sprintf(idata.notes,
                "   T (s) = %-12.4f, r = %-14.4f, z_bin= %-10.4f, "
                "r_bin = %-8.4f, w_bin (deg) = %-7.3f\n", mdata.T,
                mdata.r, TWOPI * mdata.orb.x * mdata.f,
                mdata.T / mdata.orb.p, mdata.orb.w);
        vect_free(phib);
    } else {
        sprintf(idata.notes, "   T = %f and r = %f\n", mdata.T, mdata.r);
    }

    if (fclose(datfile) != 0)
        printf("Cannot close datafile.\n");

    /* Write infofile data */

    writeinf(&idata);
    exit(0);

}

double modsin(double val)
{
    /*  Normalized to area per pulse = 1 */

    return (cos(TWOPI * val) + 1.0);
}


double crab(double val)
{
    /*  Normalized to area per pulse = 1 */

    static double phsval, normalize = 1.0 / 0.135, integral;
    static double phistart, m, b;

    phsval = modf(val, &integral);
    if ((phsval <= 0.15) || (phsval > 0.8))
        return (0.0);
    else if ((phsval > 0.35) && (phsval <= 0.55))
        return (0.74074074);
    else if (phsval <= 0.2) {
        phistart = 0.15;
        m = 2.0;
        b = 0.0;
    } else if (phsval <= 0.25) {
        phistart = 0.2;
        m = 18.0;
        b = 0.1;
    } else if (phsval <= 0.3) {
        phistart = 0.25;
        m = -16.0;
        b = 1.0;
    } else if (phsval <= 0.35) {
        phistart = 0.3;
        m = -2.0;
        b = 0.2;
    } else if (phsval <= 0.6) {
        phistart = 0.55;
        m = 2.0;
        b = 0.1;
    } else if (phsval <= 0.65) {
        phistart = 0.6;
        m = 4.0;
        b = 0.2;
    } else if (phsval <= 0.7) {
        phistart = 0.65;
        m = -4.0;
        b = 0.4;
    } else if (phsval <= 0.8) {
        phistart = 0.7;
        m = -2.0;
        b = 0.2;
    }
    return (normalize * ((phsval - phistart) * m + b));
}


double spike(double val)
{
    /*  Normalized to area per pulse = 1 */

    double phsval, integral;
    static double peakst = 0.0, peakend = 0.0, normalize = 0.0;
    static int firsttime = 1;

    if (firsttime) {
        peakst = 0.5 - fwhm;
        peakend = 0.5 + fwhm;
        normalize = 1.0 / fwhm;
        firsttime = 0;
    }

    phsval = modf(val, &integral);
    if ((phsval <= peakst) || (phsval > peakend))
        return (0.0);
    else if (phsval <= 0.5)
        return ((phsval - peakst) * normalize * normalize);
    else
        return ((1.0 - (phsval - 0.5) * normalize) * normalize);
}


double gauss(double val)
{
    /*  Normalized to area per pulse = 1 */

    double phsval, dtmp, integral;
    static double sigma, sigfact;
    static int firsttime = 1;

    if (firsttime) {
        sigma = fwhm / 2.35482;
        sigfact = 1.0 / (sigma * sqrt(TWOPI));
        firsttime = 0;
    }

    phsval = modf(val, &integral);
    dtmp = (phsval - 0.5) / sigma;
    return exp(-0.5 * dtmp * dtmp) * sigfact;
}
