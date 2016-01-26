#include "presto.h"

void hunt(double *xx, int n, double x, int *jlo);
int compare_birds(const void *ca, const void *cb);

int get_birdies(char *zapfilenm, double T, double avg_vel,
                double **lobins, double **hibins)
/* Open, read, and close a text file containing frequencies (Hz)   */
/* and widths (Hz) to ignore in a pulsar search.  The text file    */
/* should have one frequency and width per line.  Lines beginning  */
/* with '#' are ignored, and so may be used as comments.           */
/* 'T' is the total length in seconds of the observation that was  */
/* FFTd.  'avg_vel' is the avg topocentric velocity (in units      */
/* of c) towards the target during the obs.  The returned arrays   */
/* are sorted in order of increasing 'lobins' and contain the low  */
/* and high Fourier freqs that mark the boundaries of the birdies  */
/* (based on 'T'and 'avg_vel').                                    */
{
    FILE *zapfile;
    double freq, width;
    char line[200];
    int ii, numzap;
    bird *birds;

    zapfile = chkfopen(zapfilenm, "r");

    /* Read the input file once to count the birdies */

    numzap = 0;
    while (!feof(zapfile)) {
        fgets(line, 200, zapfile);
        if (line[0] == '#')
            continue;
        else
            numzap++;
    }
    numzap--;

    /* Allocate the birdie arrays */

    birds = (bird *) malloc(numzap * sizeof(bird));

    /* Rewind and read the birdies for real */

    rewind(zapfile);
    ii = 0;
    while (ii < numzap) {
        fgets(line, 200, zapfile);
        if (line[0] == '#')
            continue;
        else {
            if (line[0] == 'B') {
                sscanf(line, "B%lf %lf\n", &freq, &width);
                birds[ii].lobin = (freq - 0.5 * width) * T;
                birds[ii].hibin = (freq + 0.5 * width) * T;
            } else {
                sscanf(line, "%lf %lf\n", &freq, &width);
                birds[ii].lobin = (freq - 0.5 * width) * T * (1.0 + avg_vel);
                birds[ii].hibin = (freq + 0.5 * width) * T * (1.0 + avg_vel);
            }
            /* Insure that all birds are at least 1 bin wide */
            if ((birds[ii].hibin - birds[ii].lobin) < 1.0) {
                double avgbin;

                avgbin = 0.5 * (birds[ii].hibin + birds[ii].lobin);
                birds[ii].lobin = avgbin - 0.5;
                birds[ii].hibin = avgbin + 0.5;
            }

            ii++;
        }
    }
    fclose(zapfile);

    /* Sort the birds and then transfer them to the individual arrays */

    qsort(birds, numzap, sizeof(bird), compare_birds);
    *lobins = gen_dvect(numzap);
    *hibins = gen_dvect(numzap);
    for (ii = 0; ii < numzap; ii++) {
        (*lobins)[ii] = birds[ii].lobin;
        (*hibins)[ii] = birds[ii].hibin;
        /* printf("%15g %15g\n", birds[ii].lobin, birds[ii].hibin); */
    }
    free(birds);

    printf("Read %d birdies from '%s'.\n", numzap, zapfilenm);
    return numzap;
}


int get_std_birds(char *zapfilenm, double T, double avg_vel,
                  double **basebin, int **numharm)
/* Open, read, and close a text file containing frequencies (Hz)   */
/* and the number of harmonics to zap from a FFT.  The text file   */
/* should have one frequency and number of harmonics per line.     */
/* Lines beginning with '#' are ignored (i.e. used as comments).   */
/* 'T' is the total length in seconds of the observation that was  */
/* FFTd.  'avg_vel' is the avg topocentric velocity (in units      */
/* of c) towards the target during the obs.  The returned arrays   */
/* are sorted in order of increasing 'basebins' and contain the    */
/* base Fourier freq and the number of harmonics to check.  The    */
/* base freqs are adjusted based on avg_vel.                       */
{
    FILE *zapfile;
    double freq, harm;
    char line[200];
    int ii, numzap;
    bird *birds;

    zapfile = chkfopen(zapfilenm, "r");

    /* Read the input file once to count the birdies */

    numzap = 0;
    while (!feof(zapfile)) {
        fgets(line, 200, zapfile);
        if (line[0] == '#')
            continue;
        else
            numzap++;
    }
    numzap--;

    /* Allocate the birdie arrays */

    birds = (bird *) malloc(numzap * sizeof(bird));

    /* Rewind and read the birdies for real */

    rewind(zapfile);
    ii = 0;
    while (ii < numzap) {
        fgets(line, 200, zapfile);
        if (line[0] == '#')
            continue;
        else {
            sscanf(line, "%lf %lf\n", &freq, &harm);
            birds[ii].lobin = freq * T * (1.0 + avg_vel);
            birds[ii].hibin = harm;
            ii++;
        }
    }
    fclose(zapfile);

    /* Sort the birds and then transfer them to the individual arrays */

    printf("Read %d birdies from '%s':\n", numzap, zapfilenm);
    qsort(birds, numzap, sizeof(bird), compare_birds);
    *basebin = gen_dvect(numzap);
    *numharm = gen_ivect(numzap);
    for (ii = 0; ii < numzap; ii++) {
        (*basebin)[ii] = birds[ii].lobin;
        (*numharm)[ii] = (int) birds[ii].hibin;
        /*
           printf("  %12.7g Hz for %2d harmonics\n", (*basebin)[ii]/T, 
           (*numharm)[ii]);
         */
    }
    free(birds);

    return numzap;
}


int check_to_zap(double candbin, double *lobins, double *hibins, int numzap)
/* Look at the closest birdies from the zapfile to see if our  */
/* candidate matches one of them.  If it does, return '1' for  */
/* TRUE.  If it doesn't match, return a '0' for FALSE.  Note   */
/* that the zapfreqs _must_ be in increasing order of 'lobins' */
/* since this routine keeps track of its place in the file.    */
/* Also, numzap _must be >= 2.                                 */
{
    static int index = 0;

    if (numzap < 2) {
        printf("\n\n'numzap' = %d must be >= 2 in check_to_zap().", numzap);
        printf("  Exiting.\n\n");
        exit(1);
    }

    /* If we are beyond the end of the list, return a '0' */

    if (candbin > hibins[numzap - 1])
        return 0;

    /* Find the proper index for the highest freq birdie with a */
    /* freq below that of our candidate.                        */

    index++;
    hunt(lobins - 1, numzap, candbin, &index);
    printf("%15g  %15g  %15g\n", lobins[index], candbin, lobins[index + 1]);

    /* Check the birdie freqs to see if they match. */

    if (index == numzap)
        return (candbin < hibins[index - 1]);
    else
        return (candbin < hibins[index]);
}


void hunt(double *xx, int n, double x, int *jlo)
{
    int jm, jhi, inc;
    int ascnd;

    ascnd = (xx[n] >= xx[1]);
    if (*jlo <= 0 || *jlo > n) {
        *jlo = 0;
        jhi = n + 1;
    } else {
        inc = 1;
        if ((x >= xx[*jlo]) == ascnd) {
            if (*jlo == n)
                return;
            jhi = (*jlo) + 1;
            while ((x >= xx[jhi]) == ascnd) {
                *jlo = jhi;
                inc += inc;
                jhi = (*jlo) + inc;
                if (jhi > n) {
                    jhi = n + 1;
                    break;
                }
            }
        } else {
            if (*jlo == 1) {
                *jlo = 0;
                return;
            }
            jhi = (*jlo)--;
            while ((x < xx[*jlo]) == ascnd) {
                jhi = (*jlo);
                inc <<= 1;
                if (inc >= jhi) {
                    *jlo = 0;
                    break;
                } else
                    *jlo = jhi - inc;
            }
        }
    }
    while (jhi - (*jlo) != 1) {
        jm = (jhi + (*jlo)) >> 1;
        if ((x >= xx[jm]) == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }
    if (x == xx[n])
        *jlo = n - 1;
    if (x == xx[1])
        *jlo = 1;
}
