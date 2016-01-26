/**********************************************************************
 * 
 *  dedisperse SUBROUTINE using the tree algorithm
 *
 *  this version is simply a rename of dedisperse. (mam 6 Aug 1996)
 **********************************************************************/

#include <math.h>

void fastfold_stage(float *indata, float *outdata, int maxshift,
                    int numpts, int numchan, int stage);
void tree_stage(float **, int, int, int, float **, int);


void dedisperse_tree(loop, shift, n, m, maxlag, output, ddata)
int n, m, maxlag, loop, shift;
float **ddata, **output;
{
/*
 *   REMOVES DISPERSION USING TAYLOR TREE ALGORITHM
 *
 *   ORIGINAL VERSION: JIM CORDES, IN FORTRAN
 *   ddata:   INPUT ARRAY
 *   n:       NUMBER OF TIME SAMPLES
 *   m:       NUMBER OF FREQUENCY CHANNELS
 *   shift:   NUMBER OF SAMPLES SHIFTED, MAXIMUM
 *            (e.g. shift = 1 IF LARGEST DM USED CORRESPONDS TO
 *            TIME SHIFT OF ONE TIME SAMPLE BETWEEN ADJACENT 
 *            FREQUENCY CHANNELS; LARGER DMs CAN BE SEARCH FOR SHIFT > 1)
 *   output:  OUTPUT ARRAY (TEMPORARILY USED)
 */

    int nstages, ns, nn, mm, j, i;
    float jsum[50];

/* # TRIAL DISPERSIONS = # FREQUENCY CHANNELS = m */
/* # STAGES IN ALGORITHM IS log_2(maxlag)         */
    nstages = (int) (log((float) maxlag) / log(2.0) + 0.5);


/* GO THROUGH LOOP nstage TIMES, ALTERNATING ddata AND output AS
   INPUT AND OUTPUT TO ROUTINE (TO SAVE SPACE)                 */

    for (ns = 1; ns < nstages + 1; ns++) {
        if (ns % 2 == 1)
            tree_stage(ddata, n, m, ns, output, shift);
        else
            tree_stage(output, n, m, ns, ddata, shift);
    }


/* SINCE ddata SHOULD CONTAIN DEDISPERSED DATA,
   WE NEED TO WRITE output INTO ddata IF nstages IS ODD */

    if (nstages % 2 == 1) {
        for (nn = 0; nn < n; nn++) {
            for (mm = 0; mm < m; mm++)
                ddata[nn][mm] = output[nn][mm];
        }
    }

}


/**************************************************************************/


void tree_stage(ddata, n, m, ns, output, shift)
int n, m, shift;
float **ddata, **output;
{

/*
 *   REMOVES DISPERSION USING TAYLOR TREE ALGORITHM
 *
 *   ORIGINAL VERSION: JIM CORDES, IN FORTRAN
 *   ddata:   INPUT ARRAY
 *   n:       NUMBER OF TIME SAMPLES
 *   m:       NUMBER OF FREQUENCY CHANNELS
 *   shift:   NUMBER OF SAMPLES SHIFTED, MAXIMUM
 *            (e.g. shift = 1 IF LARGEST DM USED CORRESPONDS TO
 *            TIME SHIFT OF ONE TIME SAMPLE BETWEEN ADJACENT 
 *            FREQUENCY CHANNELS; LARGER DMs CAN BE SEARCH FOR SHIFT > 1)
 *   output:  OUTPUT ARRAY (TEMPORARILY USED)
 */


    int k, j, l, l2, lj1, ns;
    int j_index1, j_index2;
    int k_index1, k_index2;
    int ngroupsize, ngroups, ndiff;
    float x;

/* FREQUENCY CHANNELS ARE PROCESSED IN GROUPS WHOSE SIZE 
   IS DETERMINED BY WHICH STAGE IS BEING PROCESSED    */

    ngroupsize = (int) pow(2.0, (float) ns);
    ngroups = m / ngroupsize;
    ndiff = ngroupsize / 2;


/* LOOP OVER TIME SAMPLES */

    for (k = 0; k < n; k++) {

        /* LOOP OVER GROUP */

        for (l = 1; l < ngroupsize + 1; l++) {

            k_index1 = k;
            k_index2 = k - l * shift / 2;

            l2 = l / 2;
            lj1 = l - 1 - l2;

            /* LOOP OVER NUMBER OF GROUPS */

            for (j = 0; j < m; j = j + ngroupsize) {

                j_index1 = j + lj1;
                j_index2 = j_index1 + ndiff;

                /* NEED TO GUARD AGAINST GOING OUT OF ARRAY BOUNDARIES
                   ON THE TIME SAMPLES (FREQUENCY CHANNELS SHOULD BE OK) */

                x = ddata[k_index1][j_index1];

                if (k_index2 >= 0 && k_index2 < n)
                    x = x + ddata[k_index2][j_index2];

                output[k][j + l - 1] = x;

            }
        }
    }
}


/**************************************************************************/


void fast_fold(float *indata, float *outdata, int maxshift, int numpts, int numchan)
/* Folds or dedisperses using the Taylor Tree algorithm,               */
/* (AA Supp, 15, p 367, 1974), which is based on the Fast-folding      */
/* algorithm (see eg. Hankins and Rickett, 1975).                      */
/* Based on code written by Jim Cordes and Maura McLaughlin.           */
/* Arguments:                                                          */
/*    'indata':  float array with 'numchan' vectors of 'numpts' each.  */
/*       (These are the raw data to be folded/de-dispersed)            */
/*    'outdata':  float array with 'numchan' vectors of 'numpts' each. */
/*       (This are the folded/de-dispersed time series)                */
/*    'maxshift':  integer describing number of bins to use for the    */
/*       unit shift required by the algorithm.                         */
/*    'numpts':  the number of points per input/output vector          */
/*    'numchan':  is the number of input/output vectors                */
{
    int nstages = 0, ii;

    /* Number of stages in the algorithm (log_2(numchan) */

    ii = numchan;
    while (ii > 1) {
        ii >>= 1;
        nstages++;
    }

    /* Call the fold loop nstages times. */

    for (ii = 1; ii <= nstages; ii++) {
        if (ii % 2)
            fastfold_stage(indata, outdata, maxshift, numpts, numchan, ii);
        else
            fastfold_stage(outdata, indata, maxshift, numpts, numchan, ii);
    }

    /* Since 'outdata' should contain the folded/dedispersed data, */
    /* we must copy the data into 'outdata' if 'nstages' is even.  */

    if (!(nstages % 2))
        memcpy(outdata, indata, sizeof(float) * numchan * numpts);
}


/**************************************************************************/


void fastfold_stage(float *indata, float *outdata, int maxshift,
                    int numpts, int numchan, int stage)
/* Performs one stage of the Taylor Tree de-dispersion algorithm,      */
/* (AA Supp, 15, p 367, 1974), which is based on the Fast-folding      */
/* algorithm (see eg. Hankins and Rickett, 1975).                      */
/* Based on code written by Jim Cordes and Maura McLaughlin.           */
/* Arguments:                                                          */
/*    'indata':  float array with 'numchan' vectors of 'numpts' each.  */
/*       (These are the raw data to be folded/de-dispersed)            */
/*    'outdata':  float array with 'numchan' vectors of 'numpts' each. */
/*       (This are the folded/de-dispersed time series)                */
/*    'maxshift':  integer describing number of bins to use for the    */
/*       unit shift required by the algorithm.                         */
/*    'numpts':  the number of points per input/output vector          */
/*    'numchan':  is the number of input/output vectors                */
/*    'stage':  is the current stage we are on in the process          */
{
    int ii, jj, jj2, kk, maxshift2, ii_index, jj_index;
    int kk_index1, kk_index2, ii_offset, ii_offset2;
    int ngroupsize, ngroups, ndiff;
    float x;

    /* Channels are processed in groups whose size depends on 'stage' */

    ngroupsize = 1 << stage;
    ngroups = numchan / ngroupsize;
    ndiff = ngroupsize / 2;
    maxshift2 = maxshift / 2;

    /* Loop over the time samples */

    for (ii = 0; ii < numpts; ii++) {

        ii_offset = ii * numpts;

        /* Loop over the group members */

        for (jj = 0; jj < ngroupsize; jj++) {

            ii_index = ii - (jj - 1) * maxshift2;
            ii_offset2 = ii_index * numpts;
            jj2 = (jj - 1) / 2;
            jj_index = jj - jj2;

            /* Loop over the number of groups */

            for (kk = 0; kk < numchan; kk += ngroupsize) {

                kk_index1 = kk + jj_index;
                kk_index2 = kk_index1 + ndiff;
                x = indata[ii_offset + kk_index1];
                if (ii_index >= 0 && ii_index < numpts)
                    x += indata[ii_offset2 + kk_index2];
                outdata[ii_offset + kk + jj] = x;
            }
        }
    }
}
