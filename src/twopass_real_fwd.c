#include "ransomfft.h"

/* Optimized "two-pass" mass storage FFT function for real data   */
/* This version utilizes a scratch file the same size as the      */
/*     original data set.                                         */

extern long long find_blocksize(long long n1, long long n2);

void realfft_scratch_fwd(multifile* infile, multifile* scratch, 
			 long long nn)
{
  long long n1, n2, bb, bb2, fp1, fp2, ii, jj, kk, ll, s1, s2;
  int i1, i2, i3, i4, move_size, n2by2;
  unsigned char *move;
  rawtype *data, *dp, *p1, *p2, *p3, *p4;
  float *f1, *datap;
  double tmp1, tmp2, h1r, h1i, h2r, h2i;
  double h2rwr, h2rwi, h2iwr, h2iwi;
  double wtemp, wpi, wpr, wi, wr, theta, delta;

  if (nn < 2)
    return;

  /* Treat the input data as a n1 (rows) x n2 (cols) */
  /* matrix.  Make sure that n2 >= n1.               */

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
  bb = find_blocksize(n1, n2);
  printf("bb = %lld, n1 = %lld, n2 = %lld, nn = %lld (%lld)\n", 
	 bb, n1, n2, nn, n1*n2);
  if (bb==0 || bb % 2){
    printf("\nCan't factor the FFT length in twopassfft_real_fwd()\n");
    printf("   into useful sizes.\n\n");
    exit(1);
  }

  /* First do n2 transforms of length n1 by  */
  /* fetching size bb x n1 blocks in memory. */

  data = gen_rawvect(bb * n2);

  /* Transpose scratch space */

  move_size = (bb + n1) / 2;
  move = (unsigned char *)malloc(move_size);

  for (ii=0; ii<n2; ii+=bb){

    /* Read a n1 (rows) x bb (cols) block of data */

    dp = data;
    for (jj=0; jj<n1; jj++){
      fp1 = sizeof(rawtype) * (jj * n2 + ii); /* File ptr */
      fseek_multifile(infile, fp1, SEEK_SET);
      fread_multifile(dp, sizeof(rawtype), bb, infile);
      dp += bb; /* Data ptr */
    }

    /* Transpose the n1 (rows) x bb (cols) block of data */

    transpose_fcomplex(data, n1, bb, move, move_size); 

    /* Do bb transforms of length n1 */

    for (jj=0; jj<bb; jj++)
      fftwcall(data + jj * n1, n1, -1);

    /* Multiply the matrix A(ii,jj) by exp(- 2 pi i jj ii / nn). */
    /* Use recursion formulas from Numerical Recipes.            */

    f1 = (float *)data;
    for (jj=0; jj<bb; jj++){
      theta = -(jj + ii) * TWOPI / (nn >> 1);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0 + wpr;
      wi = wpi;
      i1 = 2;
      i2 = 3;
      for (kk=1; kk<n1; kk++){
	tmp1 = f1[i1];
	tmp2 = f1[i2];
	f1[i1] = tmp1 * wr - tmp2 * wi;
	f1[i2] = tmp2 * wr + tmp1 * wi;
	wtemp = wr;
	wr = wtemp * wpr - wi * wpi + wr;
	wi = wi * wpr + wtemp * wpi + wi;
	i1 += 2;
	i2 += 2;
      }
      f1 += 2 * n1;
    }

    /* Write the data to the scratch file */
    
    fwrite_multifile(data, sizeof(rawtype), bb * n1, scratch);

  }

  free(move);

  /* Now do n1 transforms of length n2 by fetching      */
  /* groups of 2 size n2 (rows) x bb2 (cols) blocks and */
  /* then recombining blocks to make a real transform.  */

  bb2 = bb >> 1;

  /* Some values for the trig recursion below: */

  delta = -TWOPI / nn;
  wtemp = sin(0.5 * delta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(delta);

  /* transpose scratch space */

  move_size = (bb + n2) / 2;
  move = (unsigned char *)malloc(move_size);
  tmp1 = sizeof(rawtype) * n1;

  for (ii=0; ii<(n1/2); ii+=bb2) {

    /* Read two n2 (rows) x bb2 (cols) blocks from the file */
    /* The first block comes from the start of the file and */
    /* the second block comes from the end.                 */

    dp = data;
    fp1 = sizeof(rawtype) * ii;              /* File ptr */
    fp2 = sizeof(rawtype) * (n1 - ii - bb2); /* File ptr */
    for (jj=0; jj<n2; jj++){
      fseek_multifile(scratch, fp1, SEEK_SET);
      fread_multifile(dp, sizeof(rawtype), bb2, scratch);
      dp += bb2;   /* Data ptr */
      fp1 += tmp1; /* File ptr */
      fseek_multifile(scratch, fp2, SEEK_SET);
      fread_multifile(dp, sizeof(rawtype), bb2, scratch);
      dp += bb2;   /* Data ptr */
      fp2 += tmp1; /* File ptr */
    }

    /* Transpose the n2 (rows) x bb (cols) block of data */

    transpose_fcomplex(data, n2, bb, move, move_size); 

    /* Do bb transforms of length n2 */

    for (jj=0; jj<bb; jj++)
      fftwcall(data + jj * n2, n2, -1);

    /* Transpose the bb (rows) x n2 (cols) block of data */

    transpose_fcomplex(data, bb, n2, move, move_size); 

    /* File pointers: */

    s1 = ii + 1;
    s2 = nn / 2 - ii - bb2;

    /* Data pointers: */

    datap = (float *) data;
    p1 = data;
    p3 = data + n2 * bb - bb2;

    for (jj=0; jj<n2; jj+=bb, p1+=bb, p3-=bb){

      p2 = p1;
      p4 = p3;

      for (kk=0; kk<bb; kk++, p2+=n2, p4-=n2){

	/* Start the trig recursion: */

	theta = ((jj + kk) * n1 + ii + 1) * delta;
	wr = cos(theta);
	wi = sin(theta);

	/* Combine n and N/2-n terms as per Numerical Recipes. */

	i1 = 2 * (kk * n2 + jj);
	i2 = i1 + 1;
	i3 = 2 * (n2 * (bb - kk) - jj - bb2) + bb - 2;
	i4 = i3 + 1;

	for (ll=0; ll<bb2; ll++){
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
	  wtemp = wr;
	  wr = wtemp * wpr - wi * wpi + wr;
	  wi = wi * wpr + wtemp * wpi + wi;
	  i1 += 2;
	  i2 += 2;
	  i3 -= 2;
	  i4 -= 2;
	}

	/* Write the "low" freqs: */
	fseek_multifile(infile, sizeof(rawtype) * s1, SEEK_SET);
	fwrite_multifile(p2, sizeof(rawtype), bb2, infile);

	/* Write the "high" freqs: */
	fseek_multifile(infile, sizeof(rawtype) * s2, SEEK_SET);
	fwrite_multifile(p4, sizeof(rawtype), bb2, infile);

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

  for (jj=0, s1=0, f1=datap; jj<n2; jj++, s1+=n1, f1+=2){
    fseek_multifile(scratch, sizeof(rawtype) * s1, SEEK_SET);
    fread_multifile(f1, sizeof(float), 2, scratch);
  }

  /* FFT the array: */

  tablesixstepfft((fcomplex *) datap, n2, -1);

  /* Do the special cases of freq=0 and Nyquist freq */

  tmp1 = datap[0];
  datap[0] = tmp1 + datap[1];
  datap[1] = tmp1 - datap[1];

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

  for (jj = 1; jj < n2by2; jj++){
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

  for (jj=0, s1=0, f1=datap; jj<n2; jj++, s1+=n1, f1+=2){
    fseek_multifile(infile, sizeof(rawtype) * s1, SEEK_SET);
    fwrite_multifile(f1, sizeof(float), 2, infile);
  }

  free(datap);
}
