#include "ransomfft.h"
#include "chkio.h"

/* Optimized "two-pass" mass storage FFT function for real data   */
/* This version utilizes a scratch file the same size as the      */
/*     original data set.                                         */

void realfft_scratch_fwd(FILE * infile, FILE * scratch, long nn)
{
  long n1, n2, j, k, m, b, b2, s, s1, s2, ct;
  long ind, lorowind, hirowind, n2by2;
  int i, i1, i2, i3, i4;
  rawtype *data, *p1, *p2, *p3, *p4;
  float *f1, *datap;
  double tmp = 0.0, h1r, h1i, h2r, h2i;
  double h2rwr, h2rwi, h2iwr, h2iwi;
  double wtemp, wpi, wpr, wi, wr, theta, delta;

  if (nn < 2)
    return;

  for (n1 = 1, n2 = 0; n1 < nn / 2; n1 <<= 1, n2++);
  n1 = n2 >> 1;
  n2 -= n1;

  n1 = 1 << n1;
  n2 = 1 << n2;

  /* n2 >= n1 */

  /* treat the input data as a n1 x n2 matrix */

  /* first do n2 transforms of length n1 */
  /* by fetching n1 x b blocks in memory */

  b = Maxblocksize / n1;
  if (b > n1)
    b = n1;

  data = gen_rawvect(nn / 2 < Maxblocksize ? nn / 2 : Maxblocksize);
  datap = (float *) data;

  for (k = 0; k < n2; k += b) {
    /* read the data from the input file in b x b blocks */
    for (j = 0, p1 = data, s = k; j < n1; j += b, p1 += b) {
      for (m = 0, p2 = p1; m < b; m++, p2 += n1, s += n2) {
	chkfseek(infile, (long) (sizeof(rawtype) * s), SEEK_SET);
	chkfread(p2, sizeof(rawtype), (unsigned long) b, infile);
      }

      /* transpose the b x b block */

      transposesquare(p1, b, n1);
    }

    /* do b transforms of size n1 */

    for (j = 0, p1 = data; j < b; j++, p1 += n1)
      tablesixstepfft((fcomplex *) p1, n1, -1);

    /* then multiply the matrix A_jk by exp(- 2 pi i j k / nn) */
    /* Use recursion formulas from NR                          */

    for (j = 0, f1 = (float *) data; j < b; j++, f1 += 2 * n1) {
      theta = -(j + k) * 6.2831853071795864769 / (nn >> 1);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0 + wpr;
      wi = wpi;
      for (ct = 1, ind = 2; ct < n1; ct++, ind += 2) {
	f1[ind] = (tmp = f1[ind]) * wr - f1[ind + 1] * wi;
	f1[ind + 1] = f1[ind + 1] * wr + tmp * wi;
	wr = (wtemp = wr) * wpr - wi * wpi + wr;
	wi = wi * wpr + wtemp * wpi + wi;
      }
    }

    /* write the data to the scratch file */
    chkfwrite(data, sizeof(rawtype), (unsigned long) (b * n1), scratch);
  }

  /* then do n1 transforms of length n2          */
  /* by fetching 2 x (n2 x b) blocks in memory   */
  /* and recombine blocks to make real transform */

  b2 = Maxblocksize / n2;
  b = b2 >> 1;
  if (b > n1 / 2)
    b = n1 / 2;

  /* Some values for the trig recursion below: */

  delta = -6.2831853071795864769 / nn;
  wtemp = sin(0.5 * delta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(delta);

  for (k = 0; k < (n1 / 2); k += b) {

    /* read the data from the input file in b x b blocks */

    for (j = 0, s1 = k + 1, s2 = n1 - k - b, p1 = data; \
	 j < n2; j += b2, p1 += b2) {

      for (m = 0, p2 = p1; m < b2; m++, p2 += n2) {
	chkfseek(scratch, (long) (sizeof(rawtype) * s1), SEEK_SET);
	chkfread(p2, sizeof(rawtype), (unsigned long) b, scratch);
	chkfseek(scratch, (long) (sizeof(rawtype) * s2), SEEK_SET);
	chkfread(p2 + b, sizeof(rawtype), (unsigned long) b, scratch);
	s1 += n1;
	s2 += n1;
      }

      /* transpose the b2 x b2 blocks */

      transposesquare(p1, b2, n2);
    }

    /* do 2*b transforms of size n2y */

    for (j = 0, p1 = data; j < b2; j++, p1 += n2)
      tablesixstepfft((fcomplex *) p1, n2, -1);

    /* transpose the b2 x b2 blocks */

    for (j = 0, p1 = data; j < n2; j += b2, p1 += b2)
      transposesquare(p1, b2, n2);

    /* File pointers: */
    s1 = k + 1;
    s2 = nn / 2 - k - b;

    /* Data pointers: */
    p1 = data;
    p3 = data + n2 * b2 - b;

    for (j = 0; j < n2; j += b2, p1 += b2, p3 -= b2) {

      p2 = p1;
      p4 = p3;

      for (m = 0; m < b2; m++, p2 += n2, p4 -= n2) {

	/* Start the trig recursion: */

	theta = ((j + m) * n1 + k + 1) * delta;
	wr = cos(theta);
	wi = sin(theta);

	/* Combine n and N/2-n terms as per Numerical Recipes. */

	lorowind = 2 * (m * n2 + j);
	hirowind = 2 * (n2 * (b2 - m) - j - b);
	i1 = lorowind;
	i2 = i1 + 1;
	i3 = hirowind + b2 - 2;
	i4 = i3 + 1;

	for (i = 0; i < b; i++) {
	  h1r = 0.5 * (datap[i1] + datap[i3]);
	  h1i = 0.5 * (datap[i2] - datap[i4]);
	  h2r = 0.5 * (datap[i2] + datap[i4]);
	  h2i = -0.5 * (datap[i1] - datap[i3]);
	  h2rwr = h2r * wr;
	  h2rwi = h2r * wi;
	  h2iwr = h2i * wr;
	  h2iwi = h2i * wi;
	  datap[i1] = h1r + h2rwr - h2iwi;
	  datap[i2] = h1i + h2iwr + h2rwi;
	  datap[i3] = h1r - h2rwr + h2iwi;
	  datap[i4] = -h1i + h2iwr + h2rwi;
	  wr = (wtemp = wr) * wpr - wi * wpi + wr;
	  wi = wi * wpr + wtemp * wpi + wi;
	  i1 += 2;
	  i2 += 2;
	  i3 -= 2;
	  i4 -= 2;
	}

	/* Write the "low" freqs: */
	chkfseek(infile, (long) (sizeof(rawtype) * s1), SEEK_SET);
	chkfwrite(p2, sizeof(rawtype), (unsigned long) b, infile);

	/* Write the "high" freqs: */
	chkfseek(infile, (long) (sizeof(rawtype) * s2), SEEK_SET);
	chkfwrite(p4, sizeof(rawtype), (unsigned long) b, infile);

	s1 += n1;
	s2 -= n1;
      }
    }
  }
  free(data);

  /* Now correct the final n2 values that have not been corrected */
  /*   due to the asymetry in the recombination.                  */

  datap = gen_fvect(2 * n2);

  /* Read the n2 data points */

  for (j = 0, s1 = 0, f1 = datap; j < n2; j++, s1 += n1, f1 += 2) {
    chkfseek(scratch, (long) (sizeof(rawtype) * s1), SEEK_SET);
    chkfread(f1, sizeof(float), (unsigned long) 2, scratch);
  }

  /* FFT the array: */

  tablesixstepfft((fcomplex *) datap, n2, -1);

  /* Do the special cases of freq=0 and Nyquist freq */

  datap[0] = (tmp = datap[0]) + datap[1];
  datap[1] = tmp - datap[1];

  /* Some values for the trig recursion below: */

  theta = -6.2831853071795864769 * n1 / nn;
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = cos(theta);
  wi = wpi;

  n2by2 = n2 >> 1;
  i1 = 2;
  i2 = 3;
  i3 = n2 * 2 - 2;
  i4 = i3 + 1;

  for (j = 1; j < n2by2; j++) {
    h1r = 0.5 * (datap[i1] + datap[i3]);
    h1i = 0.5 * (datap[i2] - datap[i4]);
    h2r = 0.5 * (datap[i2] + datap[i4]);
    h2i = -0.5 * (datap[i1] - datap[i3]);
    h2rwr = h2r * wr;
    h2rwi = h2r * wi;
    h2iwr = h2i * wr;
    h2iwi = h2i * wi;
    datap[i1] = h1r + h2rwr - h2iwi;
    datap[i2] = h1i + h2iwr + h2rwi;
    datap[i3] = h1r - h2rwr + h2iwi;
    datap[i4] = -h1i + h2iwr + h2rwi;
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
    i1 += 2;
    i2 += 2;
    i3 -= 2;
    i4 -= 2;
  }

  /* Write the n2 data points and clean up */

  for (j = 0, s1 = 0, f1 = datap; j < n2; j++, s1 += n1, f1 += 2) {
    chkfseek(infile, (long) (sizeof(rawtype) * s1), SEEK_SET);
    chkfwrite(f1, sizeof(float), (unsigned long) 2, infile);
  }

  free(datap);
}
