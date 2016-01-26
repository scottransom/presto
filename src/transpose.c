#include "presto.h"
#include "fftw3.h"

fftwf_plan plan_transpose(int rows, int cols, float *in, float *out)
{
// FFTW can be tricked into doing *very* fast transposes
// (initial testing showed 6-7x faster than TOMs!)
// http://agentzlerich.blogspot.com/2010/01/using-fftw-for-in-place-matrix.html
// http://www.fftw.org/faq/section3.html#transpose
    const unsigned flags = FFTW_MEASURE;        /* other flags are possible */
    fftwf_iodim howmany_dims[2];
    howmany_dims[0].n = rows;
    howmany_dims[0].is = cols;
    howmany_dims[0].os = 1;
    howmany_dims[1].n = cols;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = rows;
    return fftwf_plan_guru_r2r( /*rank= */ 0, /*dims= */ NULL,
                               /*howmany_rank= */ 2, howmany_dims,
                               in, out, /*kind= */ NULL, flags);
}

static int TOMS_gcd(int a, int b)
/* Return the greatest common denominator of 'a' and 'b' */
{
    int r;
    do {
        r = a % b;
        a = b;
        b = r;
    } while (r != 0);

    return a;
}

short transpose_bytes(unsigned char *a, int nx, int ny, unsigned char *move,
                      int move_size)
/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 *
 * "a" is a 1D array of length ny*nx which contains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */
{
    int i, j, im, mn;
    unsigned char b, c, d;
    int ncount;
    int k;

    /* check arguments and initialize: */
    if (ny < 0 || nx < 0)
        return -1;
    if (ny < 2 || nx < 2)
        return 0;
    if (move_size < 1)
        return -2;

    if (ny == nx) {
        /*
         * if matrix is square, exchange elements a(i,j) and a(j,i):
         */
        for (i = 0; i < nx; ++i)
            for (j = i + 1; j < nx; ++j) {
                b = a[i + j * nx];
                a[i + j * nx] = a[j + i * nx];
                a[j + i * nx] = b;
            }
        return 0;
    }
    ncount = 2;                 /* always at least 2 fixed points */
    k = (mn = ny * nx) - 1;

    for (i = 0; i < move_size; ++i)
        move[i] = 0;

    if (ny >= 3 && nx >= 3)
        ncount += TOMS_gcd(ny - 1, nx - 1) - 1; /* # fixed points */

    i = 1;
    im = ny;

    while (1) {
        int i1, i2, i1c, i2c;
        int kmi;

    /** Rearrange the elements of a loop
	and its companion loop: **/

        i1 = i;
        kmi = k - i;
        b = a[i1];
        i1c = kmi;
        c = a[i1c];

        while (1) {
            i2 = ny * i1 - k * (i1 / nx);
            i2c = k - i2;
            if (i1 < move_size)
                move[i1] = 1;
            if (i1c < move_size)
                move[i1c] = 1;
            ncount += 2;
            if (i2 == i)
                break;
            if (i2 == kmi) {
                d = b;
                b = c;
                c = d;
                break;
            }
            a[i1] = a[i2];
            a[i1c] = a[i2c];
            i1 = i2;
            i1c = i2c;
        }
        a[i1] = b;
        a[i1c] = c;

        if (ncount >= mn)
            break;              /* we've moved all elements */

    /** Search for loops to rearrange: **/

        while (1) {
            int max;

            max = k - i;
            ++i;
            if (i > max)
                return i;
            im += ny;
            if (im > k)
                im -= k;
            i2 = im;
            if (i == i2)
                continue;
            if (i >= move_size) {
                while (i2 > i && i2 < max) {
                    i1 = i2;
                    i2 = ny * i1 - k * (i1 / nx);
                }
                if (i2 == i)
                    break;
            } else if (!move[i])
                break;
        }
    }

    return 0;
}



short transpose_float(float *a, int nx, int ny, unsigned char *move, int move_size)
/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 *
 * "a" is a 1D array of length ny*nx which contains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */
{
    int i, j, im, mn;
    float b, c, d;
    int ncount;
    int k;

    /* check arguments and initialize: */
    if (ny < 0 || nx < 0)
        return -1;
    if (ny < 2 || nx < 2)
        return 0;
    if (move_size < 1)
        return -2;

    if (ny == nx) {
        /*
         * if matrix is square, exchange elements a(i,j) and a(j,i):
         */
        for (i = 0; i < nx; ++i)
            for (j = i + 1; j < nx; ++j) {
                b = a[i + j * nx];
                a[i + j * nx] = a[j + i * nx];
                a[j + i * nx] = b;
            }
        return 0;
    }
    ncount = 2;                 /* always at least 2 fixed points */
    k = (mn = ny * nx) - 1;

    for (i = 0; i < move_size; ++i)
        move[i] = 0;

    if (ny >= 3 && nx >= 3)
        ncount += TOMS_gcd(ny - 1, nx - 1) - 1; /* # fixed points */

    i = 1;
    im = ny;

    while (1) {
        int i1, i2, i1c, i2c;
        int kmi;

    /** Rearrange the elements of a loop
	and its companion loop: **/

        i1 = i;
        kmi = k - i;
        b = a[i1];
        i1c = kmi;
        c = a[i1c];

        while (1) {
            i2 = ny * i1 - k * (i1 / nx);
            i2c = k - i2;
            if (i1 < move_size)
                move[i1] = 1;
            if (i1c < move_size)
                move[i1c] = 1;
            ncount += 2;
            if (i2 == i)
                break;
            if (i2 == kmi) {
                d = b;
                b = c;
                c = d;
                break;
            }
            a[i1] = a[i2];
            a[i1c] = a[i2c];
            i1 = i2;
            i1c = i2c;
        }
        a[i1] = b;
        a[i1c] = c;

        if (ncount >= mn)
            break;              /* we've moved all elements */

    /** Search for loops to rearrange: **/

        while (1) {
            int max;

            max = k - i;
            ++i;
            if (i > max)
                return i;
            im += ny;
            if (im > k)
                im -= k;
            i2 = im;
            if (i == i2)
                continue;
            if (i >= move_size) {
                while (i2 > i && i2 < max) {
                    i1 = i2;
                    i2 = ny * i1 - k * (i1 / nx);
                }
                if (i2 == i)
                    break;
            } else if (!move[i])
                break;
        }
    }

    return 0;
}



short transpose_fcomplex(fcomplex * a, int nx, int ny, unsigned char *move,
                         int move_size)
/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 *
 * "a" is a 1D array of length ny*nx which contains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */
{
    int i, j, im, mn;
    fcomplex b, c, d;
    int ncount;
    int k;

    /* check arguments and initialize: */
    if (ny < 0 || nx < 0)
        return -1;
    if (ny < 2 || nx < 2)
        return 0;
    if (move_size < 1)
        return -2;

    if (ny == nx) {
        /*
         * if matrix is square, exchange elements a(i,j) and a(j,i):
         */
        for (i = 0; i < nx; ++i)
            for (j = i + 1; j < nx; ++j) {
                b = a[i + j * nx];
                a[i + j * nx] = a[j + i * nx];
                a[j + i * nx] = b;
            }
        return 0;
    }
    ncount = 2;                 /* always at least 2 fixed points */
    k = (mn = ny * nx) - 1;

    for (i = 0; i < move_size; ++i)
        move[i] = 0;

    if (ny >= 3 && nx >= 3)
        ncount += TOMS_gcd(ny - 1, nx - 1) - 1; /* # fixed points */

    i = 1;
    im = ny;

    while (1) {
        int i1, i2, i1c, i2c;
        int kmi;

    /** Rearrange the elements of a loop
	and its companion loop: **/

        i1 = i;
        kmi = k - i;
        b = a[i1];
        i1c = kmi;
        c = a[i1c];

        while (1) {
            i2 = ny * i1 - k * (i1 / nx);
            i2c = k - i2;
            if (i1 < move_size)
                move[i1] = 1;
            if (i1c < move_size)
                move[i1c] = 1;
            ncount += 2;
            if (i2 == i)
                break;
            if (i2 == kmi) {
                d = b;
                b = c;
                c = d;
                break;
            }
            a[i1] = a[i2];
            a[i1c] = a[i2c];
            i1 = i2;
            i1c = i2c;
        }
        a[i1] = b;
        a[i1c] = c;

        if (ncount >= mn)
            break;              /* we've moved all elements */

    /** Search for loops to rearrange: **/

        while (1) {
            int max;

            max = k - i;
            ++i;
            if (i > max)
                return i;
            im += ny;
            if (im > k)
                im -= k;
            i2 = im;
            if (i == i2)
                continue;
            if (i >= move_size) {
                while (i2 > i && i2 < max) {
                    i1 = i2;
                    i2 = ny * i1 - k * (i1 / nx);
                }
                if (i2 == i)
                    break;
            } else if (!move[i])
                break;
        }
    }

    return 0;
}
