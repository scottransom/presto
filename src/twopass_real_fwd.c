#include "ransomfft.h"

/* Optimized "two-pass" mass storage FFT function for real data   */
/* This version utilizes a scratch file the same size as the      */
/*     original data set.                                         */

extern long long find_blocksize(long long n1, long long n2);

void realfft_scratch_fwd(multifile* infile, multifile* scratch, 
			 long long nn)
{
  long long n1, n2, b, b2, s, s1, s2, j, k;
  int i, i1, i2, i3, i4, move_size;
  int m, ct, ind, lorowind, hirowind, n2by2;
  unsigned char *move;
  rawtype *data, *p1, *p2, *p3, *p4;
  float *f1, *datap;
  double tmp = 0.0, h1r, h1i, h2r, h2i;
  double h2rwr, h2rwi, h2iwr, h2iwi;
  double wtemp, wpi, wpr, wi, wr, theta, delta;

  if (nn < 2)
    return;

  /* treat the input data as a n1 x n2 matrix */
  /* n2 >= n1 */

   for (n1 = 1, n2 = 0; n1 < nn / 2; n1 <<= 1, n2++);
   n1 = n2 >> 1;
   n2 -= n1;
   
   n1 = 1 << n1;
   n2 = 1 << n2;
   b = Maxblocksize / n1;
   if (b > n1)
     b = n1;
   printf("b = %lld, n1 = %lld, n2 = %lld, nn = %lld (%lld)\n", 
	 b, n1, n2, nn, n1*n2);
   
  if (nn % 4 != 0){
    printf("\nLength of FFT in twopassfft_real_fwd()\n");
    printf("   must be divisible by 4.\n\n");
    exit(1);
  }
  n2 = good_factor(nn / 4) * 2;
  if (n2 == 0){
    printf("\nLength of FFT in twopassfft_real_fwd()\n");
    printf("   must be factorable\n\n");
    exit(1);
  }
  n1 = nn / (2 * n2);
  b = find_blocksize(n1, n2);
  printf("b = %lld, n1 = %lld, n2 = %lld, nn = %lld (%lld)\n", 
	 b, n1, n2, nn, n1*n2);
  if (b==0 || b % 2){
    printf("\nCan't factor the FFT length in twopassfft_real_fwd()\n");
    printf("   into useful sizes.\n\n");
    exit(1);
  }

  /* first do n2 transforms of length n1 */
  /* by fetching n1 x b blocks in memory */

  data = gen_rawvect(nn / 2 < b * n2 ? nn / 2 : b * n2);
  datap = (float *) data;

  /* transpose scratch space */

  move_size = (b + n1) / 2;
  move = (unsigned char *)malloc(move_size);

  for (k = 0; k < n2; k += b) {
    /* read the data from the input file in b x b blocks */
    for (j = 0, p1 = data, s = k; j < n1; j += b, p1 += b) {
      for (m = 0, p2 = p1; m < b; m++, p2 += n1, s += n2) {
	fseek_multifile(infile, sizeof(rawtype) * s, SEEK_SET);
	fread_multifile(p2, sizeof(rawtype), b, infile);
      }

      /* transpose the b x b block */

      transpose_fcomplex((fcomplex *) p1, b, n1, move, move_size); 
    }

    /* do b transforms of size n1 */

    for (j = 0, p1 = data; j < b; j++, p1 += n1)
      tablesixstepfft((fcomplex *) p1, n1, -1);

    /* then multiply the matrix A_jk by exp(- 2 pi i j k / nn) */
    /* Use recursion formulas from NR                          */

    for (j = 0, f1 = (float *) data; j < b; j++, f1 += 2 * n1) {
      theta = -(j + k) * TWOPI / (nn >> 1);
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
    fwrite_multifile(data, sizeof(rawtype), b * n1, scratch);
  }
  free(move);

  /* then do n1 transforms of length n2          */
  /* by fetching 2 x (n2 x b) blocks in memory   */
  /* and recombine blocks to make real transform */

  b2 = b;
  b = b2 >> 1;

  /* Some values for the trig recursion below: */

  delta = -TWOPI / nn;
  wtemp = sin(0.5 * delta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(delta);

  /* transpose scratch space */

  move_size = (b + n2) / 2;
  move = (unsigned char *)malloc(move_size);
  
  for (k = 0; k < (n1 / 2); k += b) {

    /* read the data from the input file in b x b blocks */

    for (j = 0, s1 = k + 1, s2 = n1 - k - b, p1 = data; \
	 j < n2; j += b2, p1 += b2) {

      for (m = 0, p2 = p1; m < b2; m++, p2 += n2) {
	fseek_multifile(scratch, sizeof(rawtype) * s1, SEEK_SET);
	fread_multifile(p2, sizeof(rawtype), b, scratch);
	fseek_multifile(scratch, sizeof(rawtype) * s2, SEEK_SET);
	fread_multifile(p2 + b, sizeof(rawtype), b, scratch);
	s1 += n1;
	s2 += n1;
      }

      /* transpose the b2 x b2 blocks */

      transpose_fcomplex((fcomplex *) p1, b, n2, move, move_size); 
    }

    /* do 2*b transforms of size n2y */

    for (j = 0, p1 = data; j < b2; j++, p1 += n2)
      tablesixstepfft((fcomplex *) p1, n2, -1);

    /* transpose the b2 x b2 blocks */

    for (j = 0, p1 = data; j < n2; j += b2, p1 += b2)
      transpose_fcomplex((fcomplex *) p1, b, n2, move, move_size); 

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
	fseek_multifile(infile, sizeof(rawtype) * s1, SEEK_SET);
	fwrite_multifile(p2, sizeof(rawtype), b, infile);

	/* Write the "high" freqs: */
	fseek_multifile(infile, sizeof(rawtype) * s2, SEEK_SET);
	fwrite_multifile(p4, sizeof(rawtype), b, infile);

	s1 += n1;
	s2 -= n1;
      }
    }
  }
  free(move);
  free(data);

  /* Now correct the final n2 values that have not been corrected */
  /*   due to the asymetry in the recombination.                  */

  datap = gen_fvect(2 * n2);

  /* Read the n2 data points */

  for (j = 0, s1 = 0, f1 = datap; j < n2; j++, s1 += n1, f1 += 2) {
    fseek_multifile(scratch, sizeof(rawtype) * s1, SEEK_SET);
    fread_multifile(f1, sizeof(float), 2, scratch);
  }

  /* FFT the array: */

  tablesixstepfft((fcomplex *) datap, n2, -1);

  /* Do the special cases of freq=0 and Nyquist freq */

  datap[0] = (tmp = datap[0]) + datap[1];
  datap[1] = tmp - datap[1];

  /* Some values for the trig recursion below: */

  theta = -TWOPI * n1 / nn;
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
    fseek_multifile(infile, sizeof(rawtype) * s1, SEEK_SET);
    fwrite_multifile(f1, sizeof(float), 2, infile);
  }

  free(datap);
}
