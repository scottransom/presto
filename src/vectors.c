#include "vectors.h"

float *gen_fvect(long length)
{
    float *v;

    v = (float *) fftwf_malloc((size_t) (sizeof(float) * length));
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

    v = (double *) fftwf_malloc((size_t) (sizeof(double) * length));
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

    v = (fcomplex *) fftwf_malloc((size_t) (sizeof(fcomplex) * length));
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

    v = (short *) fftwf_malloc((size_t) (sizeof(short) * length));
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

    v = (int *) fftwf_malloc((size_t) (sizeof(int) * length));
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

    v = (long *) fftwf_malloc((size_t) (sizeof(long) * length));
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

    v = (unsigned char *) fftwf_malloc((size_t) (sizeof(unsigned char) * length));
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

    v = (rawtype *) fftwf_malloc((size_t) (sizeof(rawtype) * length));
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

    m = (unsigned char **) fftwf_malloc((size_t) (nrows * sizeof(unsigned char *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_bmatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (unsigned char *)
        fftwf_malloc((size_t) ((nrows * ncols) * sizeof(unsigned char)));
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

    m = (short **) fftwf_malloc((size_t) (nrows * sizeof(short *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_smatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (short *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(short)));
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

    m = (int **) fftwf_malloc((size_t) (nrows * sizeof(int *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_imatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (int *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(int)));
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

    m = (float **) fftwf_malloc((size_t) (nrows * sizeof(float *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_fmatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (float *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(float)));
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

    m = (double **) fftwf_malloc((size_t) (nrows * sizeof(double *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_dmatrix()");
        printf("\n");
        exit(-1);
    }
    m[0] = (double *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(double)));
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

    m = (fcomplex **) fftwf_malloc((size_t) (nrows * sizeof(fcomplex *)));
    if (!m) {
        perror("\nError in 1st malloc() in gen_cmatrix()");
        printf("\n");
        exit(-1);
    }
    /* allocate rows and set pointers to them */

    m[0] = (fcomplex *) fftwf_malloc((size_t) ((nrows * ncols) * sizeof(fcomplex)));
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

    c = (float ***) fftwf_malloc((size_t) (nhgts * sizeof(float **)));
    if (!c) {
        perror("\nError in 1st malloc() in gen_f3Darr()");
        printf("\n");
        exit(-1);
    }
    c[0] = (float **) fftwf_malloc((size_t) ((nhgts * nrows) * sizeof(float *)));
    if (!c[0]) {
        perror("\nError in 2nd malloc() in gen_f3Darr()");
        printf("\n");
        exit(-1);
    }
    c[0][0] = (float *) fftwf_malloc((size_t) ((nhgts * nrows * ncols) * sizeof(float)));
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

    c = (fcomplex ***) fftwf_malloc((size_t) (nhgts * sizeof(fcomplex **)));
    if (!c) {
        perror("\nError in 1st malloc() in gen_c3Darr()");
        printf("\n");
        exit(-1);
    }
    c[0] = (fcomplex **) fftwf_malloc((size_t) ((nhgts * nrows) * sizeof(fcomplex *)));
    if (!c[0]) {
        perror("\nError in 2nd malloc() in gen_c3Darr()");
        printf("\n");
        exit(-1);
    }
    c[0][0] = (fcomplex *) fftwf_malloc((size_t) ((nhgts * nrows * ncols) * sizeof(fcomplex)));
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
    fftwf_free(vect);
}
