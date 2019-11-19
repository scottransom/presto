#include "ransomfft.h"

/* Local functions */

/* The following are misc FFTs and other required routines  */
void tablefft(fcomplex * data, long nn, int isign);
void tablefftraw(float data[], double table[], long n);
void tablesplitfft(fcomplex * data, long nn, int isign);
void tablesplitfftraw(float data[], double table[], long n, int isign);
double *maketable(long nn, int isign);
void fft_scramble(float data[], long nn);

long long good_factor(long long nn)
/* Return the factor of a number that is closest to its sqrt. */
/* If the number is prime, return 0.                          */
{
    long long pp;

    /* Optimal factoring is one factor twice the size of the other */
    /* Try this brute force...                                     */

    pp = (long long) sqrt((double) (nn / 2));
    if (2 * pp * pp == nn)
        return pp;

    /* Calculate the best (closest to each other) factors */
    /* This is certainly not the best way to do this...   */

    pp = (long long) sqrt((double) nn);
    while (pp > 1) {
        if (nn % pp == 0)
            return pp;
        pp--;
    }
    return 0;
}

void tablesixstepfft(fcomplex * indata, long nn, int isign)
/*  This is a modified version of a six-step-FFT routine from the    */
/*  apfloat() package.  It is a standard complex FFT.                */
/*  It uses a split-radix, table-look-up, six-step FFT.              */
/*  It is very fast for huge transforms due to high memory locality. */
/*  The forward transform (i.e. normal FFT) is isign=-1              */
{
    long n1, n2, jj, kk, kind;
    double wpr, wpi, wr, wi, theta, tmp1, tmp2;
    int move_size;
    unsigned char *move;
    fftwf_plan plan;
    static fftwf_plan last_plan = NULL;
    static int firsttime = 1, lastn = 0, lastisign = 0;

    /* If calling for the first time, read the wisdom file */
    if (firsttime)
        read_wisdom();

    if (nn < 2)
        return;

    /* Treat the input data as a n1 (rows) x n2 (cols) */
    /* matrix.  Make sure that n2 >= n1.               */

    n1 = good_factor(nn);
    if (n1 == 0) {
        printf("\nLength of FFT in tablesixstepfft() must be factorable\n\n");
        exit(0);
    }
    n2 = nn / n1;

    /* transpose scratch space */

    move_size = (n1 + n2) / 2;
    move = (unsigned char *) malloc(move_size);

    /* first transpose the matrix */

    transpose_fcomplex(indata, n1, n2, move, move_size);

    /* then do n2 transforms of length n1 */

    /* Use FFTW for the small transforms if available. */

    if (n1 == lastn && isign == lastisign) {
        plan = last_plan;
    } else {
        const int N1[1] = { n1 };
        if (firsttime)
            firsttime = 0;
        else
            fftwf_destroy_plan(last_plan);
        plan = fftwf_plan_many_dft(1, N1, n2,
                                   (fftwf_complex *) indata, NULL, 1, n1,
                                   (fftwf_complex *) indata, NULL, 1, n1,
                                   isign, FFTW_ESTIMATE);
        last_plan = plan;
        lastn = n1;
        lastisign = isign;
    }
    fftwf_execute(plan);

    /* transpose the matrix */

    transpose_fcomplex(indata, n2, n1, move, move_size);

    /* then multiply the matrix A_jk by exp(isign * 2 pi i j k / nn) */
    /* Use recursion formulas from NR */

    for (jj = 1; jj < n1; jj++) {
        theta = isign * jj * TWOPI / nn;
        wr = cos(theta);
        wi = sin(theta);
        tmp1 = sin(0.5 * theta);
        wpr = -2.0 * tmp1 * tmp1;
        wpi = wi;
        kind = jj * n2 + 1;
        for (kk = 1; kk < n2; kk++, kind++) {
            tmp1 = indata[kind].r;
            tmp2 = indata[kind].i;
            indata[kind].r = tmp1 * wr - tmp2 * wi;
            indata[kind].i = tmp2 * wr + tmp1 * wi;
            tmp1 = wr;
            wr = tmp1 * wpr - wi * wpi + wr;
            wi = wi * wpr + tmp1 * wpi + wi;
        }
    }

    /* then do n1 transforms of length n2 */

    /* Use FFTW for the small transforms if available. */

    if (n2 == lastn && isign == lastisign) {
        plan = last_plan;
    } else {
        const int N2[1] = { n2 };
        fftwf_destroy_plan(last_plan);
        plan = fftwf_plan_many_dft(1, N2, n1,
                                   (fftwf_complex *) indata, NULL, 1, n2,
                                   (fftwf_complex *) indata, NULL, 1, n2,
                                   isign, FFTW_ESTIMATE);
        last_plan = plan;
        lastn = n2;
        lastisign = isign;
    }
    fftwf_execute(plan);

    // Comment this out so it matches FFTW
    // Scale the FFT if it is an inverse FFT
    //if (isign == 1) {
    //   tmp1 = 1.0 / (double) nn;
    //   for (jj = 0; jj < n1 * n2; jj++) {
    //      indata[jj].r *= tmp1;
    //      indata[jj].i *= tmp1;
    //   }
    //}

    /* last transpose the matrix */

    transpose_fcomplex(indata, n1, n2, move, move_size);
    free(move);
}


void realfft(float idata[], long n, int isign)
/*  This is a modified version of the NR routine with correct (-)  */
/*  exponent.  It uses the above tablesixstepfft making it very    */
/*  fast.  The forward transform (i.e. normal FFT) is isign=-1     */
{
    long nby2, il, ih;
    double cc, h1r, h1i, h2r, h2i, h2rwr, h2iwr, h2rwi, h2iwi;
    double wr, wi, wpr, wpi, tmp1, theta;
    fcomplex *data;

    if (n % 2) {
        printf("\nrealfft() arrays lengths must be evenly divisible by 2.\n\n");
        exit(-1);
    }
    nby2 = n >> 1;
    data = (fcomplex *) idata;
    if (isign == -1) {
        cc = -0.5;
        theta = -TWOPI / (double) n;
        COMPLEXFFT(data, nby2, -1);
    } else {
        cc = 0.5;
        theta = TWOPI / (double) n;
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
    }
    /* Prep the trig recursion */
    wr = cos(theta);
    wi = sin(theta);
    tmp1 = sin(0.5 * theta);
    wpr = -2.0 * tmp1 * tmp1;
    wpi = wi;
    il = 1;                     /* n     */
    ih = nby2 - il;             /* N/2-n */
    for (; il <= (n >> 2); il++, ih--) {
        h1r = 0.5 * (data[il].r + data[ih].r);
        h1i = 0.5 * (data[il].i - data[ih].i);
        h2r = -cc * (data[il].i + data[ih].i);
        h2i = cc * (data[il].r - data[ih].r);
        h2rwr = h2r * wr;
        h2rwi = h2r * wi;
        h2iwr = h2i * wr;
        h2iwi = h2i * wi;
        data[il].r = h1r + h2rwr - h2iwi;
        data[il].i = h1i + h2iwr + h2rwi;
        data[ih].r = h1r - h2rwr + h2iwi;
        data[ih].i = -h1i + h2iwr + h2rwi;
        tmp1 = wr;
        wr = tmp1 * wpr - wi * wpi + wr;
        wi = wi * wpr + tmp1 * wpi + wi;
    }
    if (isign == -1) {
        /* Set data[0].r to Freq 0 value  */
        /* Set data[0].i to Nyquist value */
        tmp1 = data[0].r;
        data[0].r = tmp1 + data[0].i;
        data[0].i = tmp1 - data[0].i;
    } else {
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
        tmp1 = data[0].r;
        data[0].r = 0.5 * (tmp1 + data[0].i);
        data[0].i = 0.5 * (tmp1 - data[0].i);
        COMPLEXFFT(data, nby2, 1);
        tmp1 = 2.0 / (double) n;
        for (il = 0; il < n; il++)
            idata[il] *= tmp1;
    }
}


/* Various FFT routines and aux. routines */


void tablesplitfft(fcomplex * data, long nn, int isign)
{

/*  This is a split-radix Decimation in Frequency FFT */

    double *table;

    table = maketable(nn, 1);
    tablesplitfftraw((float *) data, table, nn, isign);
    fft_scramble((float *) data, nn);
    vect_free(table);
}


void tablefft(fcomplex * data, long nn, int isign)
{

/*  This is a radix-2 Gentleman-Sande or Decimation in Frequency FFT */

    double *table;

    table = maketable(nn, isign);
    tablefftraw((float *) data, table, nn);
    fft_scramble((float *) data, nn);
    vect_free(table);
}


double *maketable(long nn, int isign)
{
    long i, n;
    double tmp1, wr, wpr, wpi, wi, theta;
    double *table;

    n = (nn << 1);
    table = gen_dvect(n);
    table[0] = 1.0;
    table[1] = 0.0;
    theta = isign * (TWOPI / nn);
    wr = cos(theta);
    wi = sin(theta);
    tmp1 = sin(0.5 * theta);
    wpr = -2.0 * tmp1 * tmp1;
    wpi = wi;
    for (i = 2; i < n; i += 2) {
        table[i] = wr;
        table[i + 1] = wi;
        wr = (tmp1 = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + tmp1 * wpi + wi;
    }
    /* To check trig recursion above...
       for (i = 0; i < n; i += 2) {
       theta = isign*i*(PI/nn);
       table[i] = cos(theta);
       table[i + 1] = sin(theta);
       }
     */
    return table;
}


void fft_scramble(float data[], long nn)
{
    long i, j, m, n;
    float tempzz;

    data--;
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}


void tablefftraw(float data[], double table[], long n)
{

/*  This is a radix-2 Gentleman-Sande or Decimation in Frequency FFT */

    register long j, m = 2, p, q, k, n2 = n, n1, nn;
    register double c, s, rtmp, itmp;

    nn = (n << 1);
    while (m < nn) {
        n1 = n2;
        n2 >>= 1;
        for (j = 0, q = 0; j < n1; j += 2) {
            c = table[q];
            s = table[q + 1];
            q += m;
            for (k = j; k < nn; k += n1 * 2) {
                p = (k + n1);
                rtmp = data[k] - data[p];
                itmp = data[k + 1] - data[p + 1];
                data[k] += data[p];
                data[k + 1] += data[p + 1];
                data[p] = c * rtmp - s * itmp;
                data[p + 1] = c * itmp + s * rtmp;
            }
        }
        m <<= 1;
    }
}


void tablesplitfftraw(float data[], double table[], long n, int isign)
{

/*  This is a split-radix Decimation in Frequency FFT */

    int m, n2, j, is, id;
    register int i0, n4, n3;
    register int i0i, i1i, i2i, i3i;
    double r1, r2, s1, s2, s3, cc1, ss1, cc3, ss3;
    int a, a3, ai, a3i, ndec = n - 1;
    float *x;

    /* The following is a total HACK.  See below also. */
    if (isign == 1)
        for (j = 1; j < n * 2; j += 2)
            data[j] = -data[j];
    x = data - 2;
    n2 = n << 1;
    m = 1;
    while (m < n / 2) {
        n2 >>= 1;
        n4 = n2 >> 2;
        n3 = n2 >> 1;
        a = 0;
        for (j = 1; j <= n4; j++) {
            ai = a << 1;
            a3 = (a + (a << 1)) & ndec;
            a3i = a3 << 1;
            cc1 = table[ai];
            ss1 = table[ai + 1];
            cc3 = table[a3i];
            ss3 = table[a3i + 1];
            a = (a + m) & ndec;
            is = j;
            id = n2 << 1;
            do {
                for (i0 = is; i0 <= n - 1; i0 += id) {
                    i0i = i0 << 1;
                    i1i = i0i + n3;
                    i2i = i1i + n3;
                    i3i = i2i + n3;
                    r1 = x[i0i] - x[i2i];
                    x[i0i] += x[i2i];
                    r2 = x[i1i] - x[i3i];
                    x[i1i] += x[i3i];
                    s1 = x[i0i + 1] - x[i2i + 1];
                    x[i0i + 1] += x[i2i + 1];
                    s2 = x[i1i + 1] - x[i3i + 1];
                    x[i1i + 1] += x[i3i + 1];
                    s3 = r1 - s2;
                    r1 += s2;
                    s2 = r2 - s1;
                    r2 += s1;
                    x[i2i] = r1 * cc1 - s2 * ss1;
                    x[i2i + 1] = -s2 * cc1 - r1 * ss1;
                    x[i3i] = s3 * cc3 + r2 * ss3;
                    x[i3i + 1] = r2 * cc3 - s3 * ss3;
                }
                is = (id << 1) - n2 + j;
                id <<= 2;
            }
            while (is < n);
        }
        m <<= 1;
    }
    is = 1;
    id = 4;
    do {
        for (i0 = is; i0 <= n; i0 += id) {
            i0i = i0 << 1;
            i1i = i0i + 2;
            r1 = x[i0i];
            x[i0i] = r1 + x[i1i];
            x[i1i] = r1 - x[i1i];
            r1 = x[i0i + 1];
            x[i0i + 1] = r1 + x[i1i + 1];
            x[i1i + 1] = r1 - x[i1i + 1];
        }
        is = (id << 1) - 1;
        id <<= 2;
    }
    while (is < n);
    /* The following is a total HACK. */
    if (isign == 1)
        for (j = 1; j < n * 2; j += 2)
            data[j] = -data[j];
}
