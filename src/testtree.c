#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/vectors.h"

void dedisperse_tree(int, int, int, int, int, float **, float **);
void fast_fold(float *indata, float *outdata, int maxshift, int numpts, int numchan);

#define NUMCHAN 4
#define NUMPTS 4

int main(void)
{
    float **indata, **outdata;
    int ii, jj;

    indata = gen_fmatrix(NUMCHAN, NUMPTS);
    outdata = gen_fmatrix(NUMCHAN, NUMPTS);

    printf("\nindata:\n");
    for (ii = 0; ii < NUMCHAN; ii++) {
        for (jj = 0; jj < NUMPTS; jj++) {
            indata[ii][jj] = ii * NUMPTS + jj;
            printf(" %3.0f", indata[ii][jj]);
        }
        printf("\n");
    }

/*   dedisperse_tree(0, 1, NUMPTS, NUMCHAN, NUMCHAN, outdata, indata); */

/*   printf("\noutdata (Cornell):\n"); */
/*   for (ii=0; ii<NUMCHAN; ii++){ */
/*     for (jj=0; jj<NUMPTS; jj++){ */
/*         printf(" %3.0f", indata[ii][jj]); */
/* 	indata[ii][jj] = ii * NUMPTS + jj; */
/*     } */
/*     printf("\n"); */
/*   } */

    fast_fold(indata[0], outdata[0], 1, NUMPTS, NUMCHAN);

    printf("\noutdata (Cornell):\n");
    for (ii = 0; ii < NUMCHAN; ii++) {
        for (jj = 0; jj < NUMPTS; jj++) {
            printf(" %3.0f", outdata[ii][jj]);
        }
        printf("\n");
    }
    exit(0);
}
