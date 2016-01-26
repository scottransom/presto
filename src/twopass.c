#include "ransomfft.h"
#define Maxblocksize          67108864

/* Optimized "two-pass" mass storage FFT function for complex data  */
/* This version utilizes a scratch file the same size as the        */
/*     original data set.                                           */

long long find_blocksize(long long n1, long long n2)
{
    long long ii, minb = 8, maxb, b = 0, b1, b2;

    b1 = Maxblocksize / n1;
    b2 = Maxblocksize / n2;
    maxb = (b1 < b2) ? b1 : b2;
    for (ii = minb; ii <= maxb; ii++) {
        if (!(n1 % ii) && !(n2 % ii))
            b = ii;
    }
    return b;
}

void twopassfft_scratch(multifile * infile, multifile * scratch,
                        long long nn, int isign)
{
    long long n1, n2, bb, fp, ii, jj, kk, kind, df;
    int move_size;
    unsigned char *move;
    rawtype *data, *dp;
    double tmp1, tmp2, wtemp, wpi, wpr, wi, wr, delta;

    if (nn < 2)
        return;

    /* Treat the input data as a n1 (rows) x n2 (cols) */
    /* matrix.  Make sure that n2 >= n1.               */

    n1 = good_factor(nn);
    if (n1 == 0) {
        printf("\nLength of FFT in twopassfft_scratch() must be factorable\n\n");
        exit(1);
    }
    n2 = nn / n1;
    bb = find_blocksize(n1, n2);
    if (bb == 0) {
        printf("\nCan't factor the FFT length in twopassfft_scratch()\n");
        printf("   into useful sizes.\n\n");
        exit(1);
    }

    data = gen_rawvect(bb * n2);
    move_size = (bb + n2) / 2;
    move = (unsigned char *) malloc(move_size);

    /* First do n2 transforms of length n1 by  */
    /* fetching size bb x n1 blocks in memory. */

    for (ii = 0; ii < n2; ii += bb) {

        /* Read a n1 (rows) x bb (cols) block of data */

        dp = data;
        fp = sizeof(rawtype) * ii;
        df = sizeof(rawtype) * n2;
        for (jj = 0; jj < n1; jj++) {
            fseek_multifile(infile, fp, SEEK_SET);
            fread_multifile(dp, sizeof(rawtype), bb, infile);
            dp += bb;           /* Data ptr */
            fp += df;           /* File ptr */
        }

        /* Transpose the n1 (rows) x bb (cols) block of data */

        transpose_fcomplex(data, n1, bb, move, move_size);

        /* Do bb transforms of length n1 */

        for (jj = 0; jj < bb; jj++)
            COMPLEXFFT(data + jj * n1, n1, isign);

        /* Multiply the matrix A(ii,jj) by exp(isign 2 pi i jj ii / nn). */
        /* Use recursion formulas from Numerical Recipes.                */

        for (jj = 0; jj < bb; jj++) {
            delta = isign * TWOPI * (ii + jj) / nn;
            wr = cos(delta);
            wi = sin(delta);
            wtemp = sin(0.5 * delta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = wi;
            kind = jj * n1 + 1;
            for (kk = 1; kk < n1; kk++, kind++) {
                tmp1 = data[kind].r;
                tmp2 = data[kind].i;
                data[kind].r = tmp1 * wr - tmp2 * wi;
                data[kind].i = tmp2 * wr + tmp1 * wi;
                wtemp = wr;
                wr = wtemp * wpr - wi * wpi + wr;
                wi = wi * wpr + wtemp * wpi + wi;
            }
        }
        fwrite_multifile(data, sizeof(rawtype), bb * n1, scratch);
    }

    /* Now do n1 transforms of length n2 by fetching  */
    /* groups of size n2 (rows) x bb (cols) blocks.   */

    for (ii = 0; ii < n1; ii += bb) {

        /* Read two n2 (rows) x bb (cols) blocks from the file */

        dp = data;
        fp = sizeof(rawtype) * ii;
        df = sizeof(rawtype) * n1;
        for (jj = 0; jj < n2; jj++) {
            fseek_multifile(scratch, fp, SEEK_SET);
            fread_multifile(dp, sizeof(rawtype), bb, scratch);
            dp += bb;           /* Data ptr */
            fp += df;           /* File ptr */
        }

        /* Transpose the n2 (rows) x bb (cols) block of data */

        transpose_fcomplex(data, n2, bb, move, move_size);

        /* Do bb transforms of length n2 */

        for (jj = 0; jj < bb; jj++)
            COMPLEXFFT(data + jj * n2, n2, isign);

        /* Transpose the bb (rows) x n2 (cols) block of data */

        transpose_fcomplex(data, bb, n2, move, move_size);

        /* Scale the data if needed */

        if (isign == 1) {
            tmp1 = 1.0 / (double) nn;
            for (jj = 0; jj < n2 * bb; jj++) {
                data[jj].r *= tmp1;
                data[jj].i *= tmp1;
            }
        }


        /* Write n2 (rows) x bb (cols) blocks to the file  */

        dp = data;
        fp = sizeof(rawtype) * ii;
        df = sizeof(rawtype) * n1;
        for (jj = 0; jj < n2; jj++) {
            fseek_multifile(infile, fp, SEEK_SET);
            fwrite_multifile(dp, sizeof(rawtype), bb, infile);
            dp += bb;           /* Data ptr */
            fp += df;           /* File ptr */
        }
    }
    free(move);
    vect_free(data);
}
