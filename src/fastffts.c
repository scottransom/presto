#include "ransomfft.h"

/* Local functions */

/* The following are misc FFTs and other required routines  */
void tablefft(float data[], long nn, int isign);
void tablefftraw(float data[], double table[], long n);
void tablesplitfft(float data[], long nn, int isign);
void tablesplitfftraw(float data[], double table[], long n, int isign);
double *maketable(long nn, int isign);
void fft_scramble(float data[], long nn);

unsigned long good_factor(unsigned long nn)
/* Return the factor of a number that is closest to its sqrt. */
/* If the number is prime, return 0.                          */
{
  unsigned long pp;
  
  /* Optimal factoring is one factor twice the size of the other */
  /* Try this brute force...                                     */
  
  pp = (unsigned long) sqrt(nn/2);
  if (2 * pp * pp == nn)
    return pp;
  
  /* Calculate the best (closest to each other) factors */
  /* This is certainly not the best way to do this...   */
  
  pp = (unsigned long) sqrt(nn);
  while (pp > 1){
    if (nn % pp == 0)
      return pp;
    pp--;
  }
  return 0;
}


void tablesixstepfft(fcomplex *indata, long nn, int isign)
/*  This is a modified version of a six-step-FFT routine from the    */
/*  apfloat() package.  It is a standard complex FFT.                */
/*  It uses a split-radix, table-look-up, six-step FFT.              */
/*  It is very fast for huge transforms due to high memory locality. */
/*  The forward transform (i.e. normal FFT) is isign=-1              */

{
  long n1, n2, twon2, j, k, kind;
  double wpr, wpi, wr, wi, wtemp, theta, tmp = 0.0;
  float *p1;
  rawtype *data;
  int move_size;
  unsigned char *move;

#if defined USEFFTW

  FILE *wisdomfile;
  fftw_plan plan_forward, plan_inverse;
  static fftw_plan last_plan_forward = NULL, last_plan_inverse = NULL;
  static int firsttime = 1, lastn = 0;
  static char wisdomfilenm[120];

#else

  double *table;

#endif

#if defined USEFFTW

  /* If calling for the first time, read the wisdom file */

  if (firsttime) {
    sprintf(wisdomfilenm, "%s/fftw_wisdom.txt", DATABASE);
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
      printf("Error opening '%s'.  Run makewisdom again.\n", \
	     wisdomfilenm);
      printf("Exiting.\n");
      exit(1);
    }
    if (FFTW_FAILURE == fftw_import_wisdom_from_file(wisdomfile)) {
      printf("Error importing FFTW wisdom.\n");
      printf("Exiting.\n");
      exit(1);
    }
    fclose(wisdomfile);
  }

#endif

  if (nn < 2)
    return;

  data = (rawtype *) indata;

  /* treat the input data as a n1 x n2 matrix */
  /* with n2 >= n1                            */

  n1 = good_factor(nn);
  if (n1 == 0){
    printf("\nLength of FFT in tablesixstepfft() must be factorable\n\n");
    exit(0);
  }
  n2 = nn / n1;

  /* transpose scratch space */

  move_size = (n1 + n2) / 2;
  move = (unsigned char *)malloc(move_size);

  /* first transpose the matrix */

  transpose_fcomplex((fcomplex *) data, n1, n2, move, move_size); 

  /* then do n2 transforms of length n1 */

#if defined USEFFTW

  /* Use FFTW for the small transforms if available. */

  if (n1 == lastn){
    plan_forward = last_plan_forward;
    plan_inverse = last_plan_inverse;
  } else {
    if (!firsttime){
      fftw_destroy_plan(last_plan_forward);
      fftw_destroy_plan(last_plan_inverse);
    }
    plan_forward = fftw_create_plan(n1, -1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    plan_inverse = fftw_create_plan(n1, +1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    last_plan_forward = plan_forward;
    last_plan_inverse = plan_inverse;
    lastn = n1;
  }
  firsttime = 0;

  if (isign == -1) {
    fftw(plan_forward, n2, (FFTW_COMPLEX *) indata, 1, n1, \
	 NULL, 1, n1);
  } else {
    fftw(plan_inverse, n2, (FFTW_COMPLEX *) indata, 1, n1, \
	 NULL, 1, n1);
  }

#else

  table = maketable(n1, 1);
  for (k = 0, p1 = (float *)indata; k < n2; k++, p1 += 2 * n1) {
    tablesplitfftraw(p1, table, n1, isign);
    fft_scramble(p1, n1);
  }
  free(table);

#endif

  /* transpose the matrix */

  transpose_fcomplex((fcomplex *) data, n2, n1, move, move_size); 

  /* then multiply the matrix A_jk by exp(isign * 2 pi i j k / nn) */
  /* Use recursion formulas from NR */

  twon2 = 2 * n2;

  for (j = 1, p1 = (float *)indata + twon2; j < n1; j++, p1 += twon2) {
    theta = isign * j * TWOPI / nn;
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    for (k = 1, kind = 2; k < n2; k++, kind += 2) {
      p1[kind] = (tmp = p1[kind]) * wr - p1[kind + 1] * wi;
      p1[kind + 1] = p1[kind + 1] * wr + tmp * wi;
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
  }

  /* then do n1 transforms of length n2 */

#if defined USEFFTW

  /* Use FFTW for the small transforms if available. */

  if (n2 == lastn){
    plan_forward = last_plan_forward;
    plan_inverse = last_plan_inverse;
  } else {
    if (!firsttime){
      fftw_destroy_plan(last_plan_forward);
      fftw_destroy_plan(last_plan_inverse);
    }
    plan_forward = fftw_create_plan(n2, -1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    plan_inverse = fftw_create_plan(n2, +1, FFTW_MEASURE | \
				    FFTW_USE_WISDOM | \
				    FFTW_IN_PLACE);
    last_plan_forward = plan_forward;
    last_plan_inverse = plan_inverse;
    lastn = n2;
  }

  if (isign == -1) {
    fftw(plan_forward, n1, (FFTW_COMPLEX *) indata, 1, n2, \
	 NULL, 1, n2);
  } else {
    fftw(plan_inverse, n1, (FFTW_COMPLEX *) indata, 1, n2, \
	 NULL, 1, n2);
  }

#else

  table = maketable(n2, 1);
  for (k = 0, p1 = (float *)indata; k < n1; k++, p1 += 2 * n2) {
    tablesplitfftraw(p1, table, n2, isign);
    fft_scramble(p1, n2);
  }
  free(table);

#endif

  /* last transpose the matrix */

  transpose_fcomplex((fcomplex *) data, n1, n2, move, move_size); 
  free(move);
}



void realfft(float idata[], long n, int isign)
/*  This is a modified version of the NR routine with correct (-)  */
/*  exponent.  It uses the above tablesixstepfft making it very    */
/*  fast.  The forward transform (i.e. normal FFT) is isign=-1     */

{
  long i, i1, i2, i3, i4, np3;
  float c1=0.5, c2, h1r, h1i, h2r, h2i;
  double wr, wi, wpr, wpi, wtemp, theta;
  float *data;

  if (n % 2){
    printf("\nrealfft() can only handle arrays of even length.\n\n");
    exit(-1);
  }
  data = idata - 1;
  theta = TWOPI / (double) (n);
  if (isign == -1) {
    c2 = -0.5;
    COMPLEXFFT((fcomplex *)idata, n >> 1, -1);
    theta = -theta;
  } else {
    c2 = 0.5;
    /* Numerical Recipes gives a sign error for */
    /* the imaginary part of frequency n/2.     */
    if ((n+2)%4)
      data[(n>>1)+2] = -data[(n>>1)+2];
  }
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n + 3;
  for (i = 2; i <= ((n + 2) >> 2); i++) {
    i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
    h1r = c1 * (data[i1] + data[i3]);
    h1i = c1 * (data[i2] - data[i4]);
    h2r = -c2 * (data[i2] + data[i4]);
    h2i = c2 * (data[i1] - data[i3]);
    data[i1] = h1r + wr * h2r - wi * h2i;
    data[i2] = h1i + wr * h2i + wi * h2r;
    data[i3] = h1r - wr * h2r + wi * h2i;
    data[i4] = -h1i + wr * h2i + wi * h2r;
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
  }
  if (isign == -1) {
    /* Numerical Recipes gives a sign error for */
    /* the imaginary part of frequency n/2.     */
    if ((n+2)%4)
      data[(n>>1)+2] = -data[(n>>1)+2];
    data[1] = (h1r = data[1]) + data[2];
    /* This sets data[2]=Nyquist Freq value */
    data[2] = h1r - data[2];
  } else {
    double norm;
    data[1] = c1 * ((h1r = data[1]) + data[2]);
    data[2] = c1 * (h1r - data[2]);
    COMPLEXFFT((fcomplex *)idata, n >> 1, 1);
    norm = 2.0 / (double) n;
    for (i = 0; i < n; i++)
      idata[i] *= norm;
  }
}


/* Various FFT routines and aux. routines */


void tablesplitfft(float data[], long nn, int isign)
{

/*  This is a split-radix Decimation in Frequency FFT */

  double *table;

  table = maketable(nn, 1);
  tablesplitfftraw(data, table, nn, isign);
  fft_scramble(data, nn);
  free(table);
}


void tablefft(float data[], long nn, int isign)
{

/*  This is a radix-2 Gentleman-Sande or Decimation in Frequency FFT */

  double *table;

  table = maketable(nn, isign);
  tablefftraw(data, table, nn);
  fft_scramble(data, nn);
  free(table);
}


double *maketable(long nn, int isign)
{
  long i, n;
  double wtemp, wr, wpr, wpi, wi, theta;
  double *table;

  n = (nn << 1);
  table = gen_dvect(n);
  table[0] = 1.0;
  table[1] = 0.0;
  theta = isign * (TWOPI / nn);
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1 + wpr;
  wi = wpi;
  for (i = 2; i < n; i += 2) {
    table[i] = wr;
    table[i + 1] = wi;
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
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




