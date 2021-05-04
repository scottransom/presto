#include "vectors.h"

// following is in bytes and must be
// 1) a power of 2 and
// 2) a multiple of the size of (void *), 8 bytes on 64 bit Linux 
#define ALIGNSIZE 64

float *gen_fvect(long length)
{
    float *v;

#ifdef USE_FFTW_MALLOC
    v = (float *) fftwf_malloc((size_t) (sizeof(float) * length));
#else
    v = (float *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(float) * length));
#endif
    if (!v) {
        perror("\nError in gen_fvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


double *gen_dvect(long length)
{
    double *v;

#ifdef USE_FFTW_MALLOC
    v = (double *) fftwf_malloc((size_t) (sizeof(double) * length));
#else
    v = (double *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(double) * length));
#endif
    if (!v) {
        perror("\nError in gen_dvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


fcomplex *gen_cvect(long length)
{
    fcomplex *v;

#ifdef USE_FFTW_MALLOC
    v = (fcomplex *) fftwf_malloc((size_t) (sizeof(fcomplex) * length));
#else
    v = (fcomplex *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(fcomplex) * length));
#endif
    if (!v) {
        perror("\nError in gen_cvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


short *gen_svect(long length)
{
    short *v;

#ifdef USE_FFTW_MALLOC
    v = (short *) fftwf_malloc((size_t) (sizeof(short) * length));
#else
    v = (short *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(short) * length));
#endif
    if (!v) {
        perror("\nError in gen_svect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


int *gen_ivect(long length)
{
    int *v;

#ifdef USE_FFTW_MALLOC
    v = (int *) fftwf_malloc((size_t) (sizeof(int) * length));
#else
    v = (int *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(int) * length));
#endif
    if (!v) {
        perror("\nError in gen_ivect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


long *gen_lvect(long length)
{
    long *v;

#ifdef USE_FFTW_MALLOC
    v = (long *) fftwf_malloc((size_t) (sizeof(long) * length));
#else
    v = (long *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(long) * length));
#endif
    if (!v) {
        perror("\nError in gen_lvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


unsigned char *gen_bvect(long length)
{
    unsigned char *v;

#ifdef USE_FFTW_MALLOC
    v = (unsigned char *) fftwf_malloc((size_t) (sizeof(unsigned char) * length));
#else
    v = (unsigned char *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(unsigned char) * length));
#endif
    if (!v) {
        perror("\nError in gen_bvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


rawtype *gen_rawvect(long length)
/* Allocate a vector of length 'length' rawtype points */
{
    rawtype *v;

#ifdef USE_FFTW_MALLOC
    v = (rawtype *) fftwf_malloc((size_t) (sizeof(rawtype) * length));
#else
    v = (rawtype *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(rawtype) * length));
#endif
    if (!v) {
        perror("\nError in gen_rawvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}


unsigned char **gen_bmatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_bmatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    unsigned char **m;

#ifdef USE_FFTW_MALLOC
    m = (unsigned char **) fftwf_malloc((size_t) (nrows * sizeof(unsigned char *)));
#else
    m = (unsigned char **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(unsigned char *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_bmatrix()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    m[0] = (unsigned char *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(unsigned char)));
#else
    m[0] = (unsigned char *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(unsigned char)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_bmatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


short **gen_smatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_smatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    short **m;

#ifdef USE_FFTW_MALLOC
    m = (short **) fftwf_malloc((size_t) (nrows * sizeof(short *)));
#else
    m = (short **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(short *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_smatrix()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    m[0] = (short *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(short)));
#else
    m[0] = (short *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(short)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_smatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


int **gen_imatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_imatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    int **m;

#ifdef USE_FFTW_MALLOC
    m = (int **) fftwf_malloc((size_t) (nrows * sizeof(int *)));
#else
    m = (int **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(int *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_imatrix()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    m[0] = (int *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(int)));
#else
    m[0] = (int *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(int)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_imatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


float **gen_fmatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_fmatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    float **m;

#ifdef USE_FFTW_MALLOC
    m = (float **) fftwf_malloc((size_t) (nrows * sizeof(float *)));
#else
    m = (float **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(float *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_fmatrix()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    m[0] = (float *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(float)));
#else
    m[0] = (float *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(float)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_fmatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


double **gen_dmatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_dmatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    double **m;

#ifdef USE_FFTW_MALLOC
    m = (double **) fftwf_malloc((size_t) (nrows * sizeof(double *)));
#else
    m = (double **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(double *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_dmatrix()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    m[0] = (double *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(double)));
#else
    m[0] = (double *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(double)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_dmatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;
    return m;
}


fcomplex **gen_cmatrix(long nrows, long ncols)
{
    /* Note:  To free this matrix, assuming you called it with:    */
    /*             x = gen_cmatrix(10,10);                         */
    /*        all you need to do is the following:                 */
    /*             free(x[0]) ; free(x) ;                          */
    /*        The order is important!                              */

    long i;
    fcomplex **m;

    /* allocate pointers to rows */

#ifdef USE_FFTW_MALLOC
    m = (fcomplex **) fftwf_malloc((size_t) (nrows * sizeof(fcomplex *)));
#else
    m = (fcomplex **) aligned_alloc(ALIGNSIZE, (size_t) (nrows * sizeof(fcomplex *)));
#endif
    if (!m) {
        perror("\nError in 1st malloc() in gen_cmatrix()");
        printf("\n");
        exit(-1);
    }
    /* allocate rows and set pointers to them */

#ifdef USE_FFTW_MALLOC
    m[0] = (fcomplex *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(fcomplex)));
#else
    m[0] = (fcomplex *) aligned_alloc(ALIGNSIZE, (size_t) ((nrows * ncols) * sizeof(fcomplex)));
#endif
    if (!m[0]) {
        perror("\nError in 2nd malloc() in gen_cmatrix()");
        printf("\n");
        exit(-1);
    }
    for (i = 1; i < nrows; i++)
        m[i] = m[i - 1] + ncols;

    /* return pointer to array of pointers to rows */

    return m;
}

float ***gen_f3Darr(long nhgts, long nrows, long ncols)
{
    /* Note:  To free this 3D array, assuming you called it with:  */
    /*             x = gen_f3Darr(10,10,10);                       */
    /*        all you need to do is the following:                 */
    /*             free(x[0][0]) ; free(x[0]) ; free(x) ;          */
    /*        The order is important!                              */

    long i, j;
    float ***c;

#ifdef USE_FFTW_MALLOC
    c = (float ***) fftwf_malloc((size_t) (nhgts * sizeof(float **)));
#else
    c = (float ***) aligned_alloc(ALIGNSIZE, (size_t) (nhgts * sizeof(float **)));
#endif
    if (!c) {
        perror("\nError in 1st malloc() in gen_f3Darr()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    c[0] = (float **) fftwf_malloc((size_t) ((nhgts * nrows) * sizeof(float *)));
#else
    c[0] = (float **) aligned_alloc(ALIGNSIZE, (size_t) ((nhgts * nrows) * sizeof(float *)));
#endif
    if (!c[0]) {
        perror("\nError in 2nd malloc() in gen_f3Darr()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    c[0][0] = (float *) fftwf_malloc((size_t) ((nhgts * nrows * ncols) * sizeof(float)));
#else
    c[0][0] = (float *) aligned_alloc(ALIGNSIZE, (size_t) ((nhgts * nrows * ncols) * sizeof(float)));
#endif
    if (!c[0][0]) {
        perror("\nError in 3rd malloc() in gen_f3Darr()");
        printf("\n");
        exit(-1);
    }

    for (j = 1; j < nrows; j++)
        c[0][j] = c[0][j - 1] + ncols;

    for (i = 1; i < nhgts; i++) {
        c[i] = c[i - 1] + nrows;
        c[i][0] = c[i - 1][0] + nrows * ncols;
        for (j = 1; j < nrows; j++)
            c[i][j] = c[i][j - 1] + ncols;
    }

    return c;
}

fcomplex ***gen_c3Darr(long nhgts, long nrows, long ncols)
{
    /* Note:  To free this 3D array, assuming you called it with:  */
    /*             x = gen_f3Darr(10,10,10);                       */
    /*        all you need to do is the following:                 */
    /*             free(x[0][0]) ; free(x[0]) ; free(x) ;          */
    /*        The order is important!                              */

    long i, j;
    fcomplex ***c;

#ifdef USE_FFTW_MALLOC
    c = (fcomplex ***) fftwf_malloc((size_t) (nhgts * sizeof(fcomplex **)));
#else
    c = (fcomplex ***) aligned_alloc(ALIGNSIZE, (size_t) (nhgts * sizeof(fcomplex **)));
#endif
    if (!c) {
        perror("\nError in 1st malloc() in gen_c3Darr()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    c[0] = (fcomplex **) fftwf_malloc((size_t) ((nhgts * nrows) * sizeof(fcomplex *)));
#else
    c[0] = (fcomplex **) aligned_alloc(ALIGNSIZE, (size_t) ((nhgts * nrows) * sizeof(fcomplex *)));
#endif
    if (!c[0]) {
        perror("\nError in 2nd malloc() in gen_c3Darr()");
        printf("\n");
        exit(-1);
    }
#ifdef USE_FFTW_MALLOC
    c[0][0] = (fcomplex *) fftwf_malloc((size_t) ((nhgts * nrows * ncols) * sizeof(fcomplex)));
#else
    c[0][0] = (fcomplex *) aligned_alloc(ALIGNSIZE, (size_t) ((nhgts * nrows * ncols) * sizeof(fcomplex)));
#endif
    if (!c[0][0]) {
        perror("\nError in 3rd malloc() in gen_c3Darr()");
        printf("\n");
        exit(-1);
    }

    for (j = 1; j < nrows; j++)
        c[0][j] = c[0][j - 1] + ncols;

    for (i = 1; i < nhgts; i++) {
        c[i] = c[i - 1] + nrows;
        c[i][0] = c[i - 1][0] + nrows * ncols;
        for (j = 1; j < nrows; j++)
            c[i][j] = c[i][j - 1] + ncols;
    }

    return c;
}

void vect_free(void *vect)
{
#ifdef USE_FFTW_MALLOC
    fftwf_free(vect);
#else
    free(vect);
#endif
}
