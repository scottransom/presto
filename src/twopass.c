#include "ransomfft.h"

/* Optimized "two-pass" mass storage FFT function for complex data  */
/* This version utilizes a scratch file the same size as the        */
/*     original data set.                                           */

long long find_blocksize(long long n1, long long n2)
{
  long long ii, minb=8, maxb, b=0, b1, b2;

  b1 = Maxblocksize / n1;
  b2 = Maxblocksize / n2;
  maxb = (b1 < b2) ? b1 : b2;
  for (ii=minb; ii<=maxb; ii++){
    if (!(n1 % ii) && !(n2 % ii))
      b = ii;
  }
  return b;
}


void twopassfft_scratch(multifile* infile, multifile* scratch, \
			long long nn, int isign)
{
  long long n1, n2, b, s, j, k;
  int m, ct, ind, move_size;
  unsigned char *move;
  rawtype *data, *p1, *p2;
  float *f1;
  double tmp = 0.0, wtemp, wpi, wpr, wi, wr, theta;

  if (nn < 2)
    return;

  /* treat the input data as a n1 x n2 matrix */
  /* n2 >= n1 */

  if (nn % 4 != 0){
    printf("\nLength of FFT in twopassfft_scratch() must be divisible by 4.\n\n");
    exit(1);
  }
  n1 = good_factor(nn / 4) * 2;
  if (n1 == 0){
    printf("\nLength of FFT in twopassfft_scratch() must be factorable\n\n");
    exit(1);
  }
  n2 = nn / n1;
  b = find_blocksize(n1, n2);
  if (b==0){
    printf("\nCan't factor the FFT length in twopassfft_scratch()\n");
    printf("   into useful sizes.\n\n");
    exit(1);
  }

  /* first do n2 transforms of length n1 */
  /* by fetching n1 x b blocks in memory */

  data = gen_rawvect(nn < Maxblocksize ? nn : Maxblocksize);

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
      tablesixstepfft((fcomplex *) p1, n1, isign);

    /* then multiply the matrix A_jk by exp(isign * 2 pi i j k / nn) */
    /* Use recursion formulas from NR                                */

    for (j = 0, f1 = (float *) data; j < b; j++, f1 += 2 * n1) {
      theta = isign * (j + k) * TWOPI / nn;
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

  /* then do n1 transforms of length n2  */
  /* by fetching n2 x b blocks in memory */

  /* transpose scratch space */

  move_size = (b + n2) / 2;
  move = (unsigned char *)malloc(move_size);

  for (k = 0; k < n1; k += b) {
    /* read the data from the input file in b x b blocks */

    for (j = 0, p1 = data, s = k; j < n2; j += b, p1 += b) {
      for (m = 0, p2 = p1; m < b; m++, p2 += n2, s += n1) {
	fseek_multifile(scratch, sizeof(rawtype) * s, SEEK_SET);
	fread_multifile(p2, sizeof(rawtype), b, scratch);
      }

      /* transpose the b x b block */

      transpose_fcomplex((fcomplex *) p1, b, n2, move, move_size); 
    }

    /* do b transforms of size n2 */

    for (j = 0, p1 = data; j < b; j++, p1 += n2)
      tablesixstepfft((fcomplex *) p1, n2, isign);

    /* write the data to the original file */

    for (j = 0, p1 = data, s = k; j < n2; j += b, p1 += b) {
      /* transpose the b x b block */

      transpose_fcomplex((fcomplex *) p1, b, n2, move, move_size); 

      for (m = 0, p2 = p1; m < b; m++, p2 += n2, s += n1) {
	fseek_multifile(infile, sizeof(rawtype) * s, SEEK_SET);
	fwrite_multifile(p2, sizeof(rawtype), b, infile);
      }
    }
  }
  free(move);
  free(data);
}
