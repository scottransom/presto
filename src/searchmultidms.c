#include "presto.h"
#include "plot2d.h"
#include "randlib.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *infile;
    int numchan, numtrials = 1000;
    long ii, jj, kk, proflen;
    float *chisq, *dms;
    double *profs, *sumprof, *tmpprof, nc, pl, tt, bt, p, f, df;
    double dm = 0.0, lodm = 0.0, ddm = 0.2, bindelay = 0.0;
    double chixmeas = 0.0, chitmp = 0.0, maxchi = 0.0, maxdm = 0.0;
    double varph = 0.0, avgph = 0.0, skew = 0.0, kurt = 0.0;
    char device[200], output[200], ylab[200];

    if (argc < 2 || argc > 6) {
        printf("\nusage:  searchmultidms file [trials] [ddm] [lodm] [dev]\n");
        printf("  'file'    (required):  The multi-profile save file.\n");
        printf("  'trials'  (optional):  Number of trials.  Default = 2000\n");
        printf("  'ddm'     (optional):  The DM step.  Default = 0.2\n");
        printf("  'lodm'    (optional):  The low DM to check.  Default = 0.0\n");
        printf("  'dev'     (optional):  Pgplot device to use ('x' or 'ps').\n");
        printf("                            Default = 'ps'\n\n");
        exit(1);
    }

    infile = chkfopen(argv[1], "rb");
    sprintf(output, "%s_dmsearch.ps", argv[1]);
    chkfread(&nc, sizeof(double), 1, infile);
    chkfread(&pl, sizeof(double), 1, infile);
    chkfread(&p, sizeof(double), 1, infile);
    chkfread(&tt, sizeof(double), 1, infile);
    chkfread(&bt, sizeof(double), 1, infile);
    chkfread(&f, sizeof(double), 1, infile);
    chkfread(&df, sizeof(double), 1, infile);
    numchan = nc;
    proflen = pl;

    if (argc == 2) {
        strcpy(device, "ps");
    } else if (argc == 3) {
        strcpy(device, "ps");
        numtrials = strtol(argv[2], (char **) NULL, 10);
    } else if (argc == 4) {
        strcpy(device, "ps");
        numtrials = strtol(argv[2], (char **) NULL, 10);
        ddm = strtod(argv[3], (char **) NULL);
    } else if (argc == 5) {
        strcpy(device, "ps");
        numtrials = strtol(argv[2], (char **) NULL, 10);
        ddm = strtod(argv[3], (char **) NULL);
        lodm = strtod(argv[4], (char **) NULL);
    } else if (argc == 6) {
        strcpy(device, "ps");
        numtrials = strtol(argv[2], (char **) NULL, 10);
        ddm = strtod(argv[3], (char **) NULL);
        lodm = strtod(argv[4], (char **) NULL);
        strcpy(device, argv[5]);
    }

    printf("\n   Multi-Profile DM Optimization Program\n");
    printf("              Scott M. Ransom\n");
    printf("               7 March 1999\n");
    printf("\nProfile properties:\n");
    printf("Initial folding period   (s)  =  %-15.13f\n", p);
    printf("Topocentric time     (start)  =  %-15.10f\n", tt);
    printf("Barycentric time     (start)  =  %-15.10f\n", bt);
    printf("Profile length        (bins)  =  %-ld\n", proflen);
    printf("Number of channels            =  %-d\n", numchan);
    printf("Channel 1 frequency    (MHz)  =  %-10.5f\n", f);
    printf("Channel freq width     (MHz)  =  %-10.5f\n\n", df);

    /* Read the profiles. */

    profs = gen_dvect(proflen * numchan);
    chkfread(profs, sizeof(double), (unsigned long) (numchan * proflen), infile);
    fclose(infile);

    /* Create a Summed-Profile vector */

    sumprof = gen_dvect(proflen);
    for (ii = 0; ii < proflen; ii++)
        sumprof[ii] = 0.0;
    tmpprof = gen_dvect(proflen);
    for (ii = 0; ii < proflen; ii++)
        tmpprof[ii] = 0.0;

    /* Create a Chi-Squared Storage Vector */

    chisq = gen_fvect(numtrials);
    dms = gen_fvect(numtrials);

    /* Randomly Rotate the vectors and sum the profiles  */
    /* in order to estimate the mean and variance of the */
    /* summed profile.                                   */

    for (ii = 0; ii < numchan; ii++) {
        bindelay = delay_from_dm(dm, f + ii * df) * (double) proflen / p;
        memcpy(tmpprof, &profs[ii * proflen], sizeof(double) * proflen);
        drotate(tmpprof, proflen, -genunf(0.0, (float) proflen - 1));
        for (jj = 0; jj < proflen; jj++)
            sumprof[jj] += tmpprof[jj];
    }

    /* Determine some stats for this trial */

    dstats(sumprof, proflen, &avgph, &varph, &skew, &kurt);

    for (ii = 0; ii < numtrials; ii++) {

        /* Set the DM */

        dm = lodm + ii * ddm;
        dms[ii] = dm;

        /* Reset the summed-profile to all zeros */

        for (jj = 0; jj < proflen; jj++)
            sumprof[jj] = 0.0;

        /* Rotate the vectors and sum the profiles */

        for (jj = 0; jj < numchan; jj++) {
            bindelay = delay_from_dm(dm, f + jj * df) * (double) proflen / p;
            memcpy(tmpprof, &profs[jj * proflen], sizeof(double) * proflen);
            drotate(tmpprof, proflen, bindelay);
            for (kk = 0; kk < proflen; kk++)
                sumprof[kk] += tmpprof[kk];
        }

        /* Compute the Chi-Squared probability that there is a signal */
        /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */

        chixmeas = 0.0;
        for (jj = 0; jj < proflen; jj++) {
            chitmp = sumprof[jj] - avgph;
            chixmeas += chitmp * chitmp;
        }
        chixmeas /= (varph * proflen);
        chisq[ii] = (float) chixmeas;
        printf("DM = %f  Reduced Chi-Square = %f\n", dm, chixmeas);
        if (chixmeas > maxchi) {
            maxchi = chixmeas;
            maxdm = dm;
        }
    }

    printf("\n\nThe maximum Chi-Squared of %g occured at DM = %f\n\n", maxchi,
           maxdm);

    /* Plot the Chi-Squared Results */

    if (0 == strcmp("x", device))
        cpgstart_x("landscape");
    else
        cpgstart_ps(output, "landscape");
    sprintf(ylab, "Reduced Chi-Square (NPTS = %ld)", proflen);
    xyline(numtrials, dms, chisq, "Dispersion Measure", ylab, 1);
    cpgend();

    /* Cleanup */

    vect_free(profs);
    vect_free(sumprof);
    vect_free(tmpprof);
    vect_free(chisq);
    vect_free(dms);

    return 0;
}
