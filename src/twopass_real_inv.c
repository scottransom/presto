#include "ransomfft.h"

/* Optimized "two-pass" mass storage FFT function for real data   */
/* This version utilizes a scratch file the same size as the      */
/*     original data set.                                         */
/* This is the inverse FFT.                                       */

extern long long find_blocksize(long long n1, long long n2);

void realfft_scratch_inv(multifile * infile, multifile * scratch, long long nn)
{
    long long n1, n2, bb, bb2, fp1, fp2, df, ii, jj, kk, kind;
    int i1, i2, move_size;
    unsigned char *move;
    rawtype *data, *dp;
    double tmp1, tmp2, h1r, h1i, h2r, h2i;
    double h2rwr, h2rwi, h2iwr, h2iwi;
    double wtemp, wpi, wpr, wpi2, wpr2, wi, wr, theta = 0.0, delta;

    if (nn < 2)
        return;

    /* Treat the input data as a n1 x n2 matrix */
    /* n2 >= n1 */

    if (nn % 4 != 0) {
        printf("\nLength of FFT in twopassfft_real_inv()\n");
        printf("   must be divisible by 4.\n\n");
        exit(1);
    }
    n2 = good_factor(nn / 4) * 2;
    if (n2 == 0) {
        printf("\nLength of FFT in twopassfft_real_inv()\n");
        printf("   must be factorable\n\n");
        exit(1);
    }
    n1 = nn / (2 * n2);
    bb = find_blocksize(n1, n2);
    if (bb == 0 || bb % 2 || n1 % 2 || n2 % 2) {
        printf("\nCan't factor the FFT length in twopassfft_real_inv()\n");
        printf("   into useful sizes.\n\n");
        exit(1);
    }

    /* First 'split' the data as per Numerical Recipes. */

    data = gen_rawvect(bb * n2);

    /* Some values for the trig recursion below: */

    delta = TWOPI / nn;
    wtemp = sin(0.5 * delta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(delta);

    /* Correct the first n1 values that will not be corrected */
    /*   later due to the asymetry in the recombination.      */

    /* Read the n1 data points */

    for (jj = 0; jj < n1; jj++) {
        fseek_multifile(infile, sizeof(rawtype) * jj * n2, SEEK_SET);
        fread_multifile(data + jj, sizeof(rawtype), 1, infile);
    }

    /* Do the special cases of freq=0 and Nyquist freq */

    tmp1 = data[0].r;
    data[0].r = 0.5 * (tmp1 + data[0].i);
    data[0].i = 0.5 * (tmp1 - data[0].i);

    for (jj = 1, kk = n1 - 1; jj < n1 / 2; jj++, kk--) {
        theta = delta * n2 * jj;
        wr = cos(theta);
        wi = sin(theta);
        h1r = 0.5 * (data[jj].r + data[kk].r);
        h1i = 0.5 * (data[jj].i - data[kk].i);
        h2r = -0.5 * (data[jj].i + data[kk].i);
        h2i = 0.5 * (data[jj].r - data[kk].r);
        h2rwr = h2r * wr;
        h2rwi = h2r * wi;
        h2iwr = h2i * wr;
        h2iwi = h2i * wi;
        data[jj].r = h1r + h2rwr - h2iwi;
        data[jj].i = h1i + h2iwr + h2rwi;
        data[kk].r = h1r - h2rwr + h2iwi;
        data[kk].i = -h1i + h2iwr + h2rwi;
    }

    /* FFT and write the array: */

    COMPLEXFFT(data, n1, 1);
    fwrite_multifile(data, sizeof(rawtype), n1, scratch);

    /* Now do n2 transforms of length n1 by fetching      */
    /* groups of 2 size n1 (rows) x bb2 (cols) blocks     */
    /* after splitting the data from the real transform.  */

    bb2 = bb >> 1;

    /* transpose scratch space */

    move_size = (bb + n2) / 2;
    move = (unsigned char *) malloc(move_size);
    move_size = (bb + n1) / 2;
    df = sizeof(rawtype) * n2;

    for (ii = 0; ii < (n2 / 2); ii += bb2) {

        /* Read two n1 (rows) x bb2 (cols) blocks from the file */
        /* The first block comes from the start of the file and */
        /* the second block comes from the end.                 */
        /* Note:  The first block is shifted by one complex     */
        /*        value in order to make the complex->real      */
        /*        assembly quite a bit more elegant...          */

        dp = data;
        fp1 = sizeof(rawtype) * (ii + 1);       /* File ptr */
        fp2 = sizeof(rawtype) * (n2 - ii - bb2);        /* File ptr */
        for (jj = 0; jj < n1; jj++) {
            fseek_multifile(infile, fp1, SEEK_SET);
            fread_multifile(dp, sizeof(rawtype), bb2, infile);
            dp += bb2;          /* Data ptr */
            fp1 += df;          /* File ptr */
            fseek_multifile(infile, fp2, SEEK_SET);
            fread_multifile(dp, sizeof(rawtype), bb2, infile);
            dp += bb2;          /* Data ptr */
            fp2 += df;          /* File ptr */
        }

        /* Begin the re-assembly of the realFFT */

        for (jj = 0; jj < n1; jj++) {

            /* Start the trig recursion: */

            theta = (jj * n2 + ii + 1) * delta;
            wr = cos(theta);
            wi = sin(theta);

            /* Combine n and N/2-n terms as per Numerical Recipes. */

            i1 = jj * bb;       /* n     */
            i2 = bb * n1 - i1 - 1;      /* N/2-n */
            for (kk = 0; kk < bb2; kk++, i1++, i2--) {
                h1r = 0.5 * (data[i1].r + data[i2].r);
                h1i = 0.5 * (data[i1].i - data[i2].i);
                h2r = -0.5 * (data[i1].i + data[i2].i);
                h2i = 0.5 * (data[i1].r - data[i2].r);
                h2rwr = h2r * wr;
                h2rwi = h2r * wi;
                h2iwr = h2i * wr;
                h2iwi = h2i * wi;
                data[i1].r = h1r + h2rwr - h2iwi;
                data[i1].i = h1i + h2iwr + h2rwi;
                data[i2].r = h1r - h2rwr + h2iwi;
                data[i2].i = -h1i + h2iwr + h2rwi;
                wtemp = wr;
                wr = wtemp * wpr - wi * wpi + wr;
                wi = wi * wpr + wtemp * wpi + wi;
            }
        }

        /* Transpose the n1 (rows) x bb (cols) block of data */

        transpose_fcomplex(data, n1, bb, move, move_size);

        /* Do bb transforms of length n1 */

        for (jj = 0; jj < bb; jj++)
            COMPLEXFFT(data + jj * n1, n1, 1);

        /* Multiply the matrix A(ii,jj) by exp(2 pi i jj ii / nn).   */
        /* Use recursion formulas from Numerical Recipes.            */

        for (jj = 0; jj < bb; jj++) {
            theta = (jj < bb2) ?
                TWOPI * (ii + jj + 1) / (double) (nn >> 1) :
                TWOPI * (n2 - ii - bb + jj) / (double) (nn >> 1);
            wr = cos(theta);
            wi = sin(theta);
            wtemp = sin(0.5 * theta);
            wpr2 = -2.0 * wtemp * wtemp;
            wpi2 = wi;
            kind = jj * n1 + 1;
            for (kk = 1; kk < n1; kk++, kind++) {
                tmp1 = data[kind].r;
                tmp2 = data[kind].i;
                data[kind].r = tmp1 * wr - tmp2 * wi;
                data[kind].i = tmp2 * wr + tmp1 * wi;
                wtemp = wr;
                wr = wtemp * wpr2 - wi * wpi2 + wr;
                wi = wi * wpr2 + wtemp * wpi2 + wi;
            }
        }
        fseek_multifile(scratch, sizeof(rawtype) * (ii + 1) * n1, SEEK_SET);
        fwrite_multifile(data, sizeof(rawtype), bb2 * n1, scratch);
        fseek_multifile(scratch, sizeof(rawtype) * (n2 - (ii + bb2)) * n1, SEEK_SET);
        fwrite_multifile(data + bb2 * n1, sizeof(rawtype), bb2 * n1, scratch);
    }

    /* Now do n1 transforms of length n2 by fetching  */
    /* groups of size n2 (rows) x bb (cols) blocks.   */

    move_size = (bb + n2) / 2;

    for (ii = 0; ii < n1; ii += bb) {

        /* Read two n2 (rows) x bb (cols) blocks from the file */

        dp = data;
        fp1 = sizeof(rawtype) * ii;
        df = sizeof(rawtype) * n1;
        for (jj = 0; jj < n2; jj++) {
            fseek_multifile(scratch, fp1, SEEK_SET);
            fread_multifile(dp, sizeof(rawtype), bb, scratch);
            dp += bb;           /* Data ptr */
            fp1 += df;          /* File ptr */
        }

        /* Transpose the n2 (rows) x bb (cols) block of data */

        transpose_fcomplex(data, n2, bb, move, move_size);

        /* Do bb transforms of length n2 */

        for (jj = 0; jj < bb; jj++)
            COMPLEXFFT(data + jj * n2, n2, 1);

        /* Transpose the bb (rows) x n2 (cols) block of data */

        transpose_fcomplex(data, bb, n2, move, move_size);

        /* Scale and write the data */

        tmp1 = 2.0 / (double) nn;
        for (jj = 0; jj < n2 * bb; jj++) {
            data[jj].r *= tmp1;
            data[jj].i *= tmp1;
        }
        dp = data;
        fp1 = sizeof(rawtype) * ii;
        df = sizeof(rawtype) * n1;
        for (jj = 0; jj < n2; jj++) {
            fseek_multifile(infile, fp1, SEEK_SET);
            fwrite_multifile(dp, sizeof(rawtype), bb, infile);
            dp += bb;           /* Data ptr */
            fp1 += df;          /* File ptr */
        }
    }
    free(move);
    vect_free(data);
}
