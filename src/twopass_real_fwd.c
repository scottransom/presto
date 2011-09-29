#include "ransomfft.h"

/* Optimized "two-pass" mass storage FFT function for real data   */
/* This version utilizes a scratch file the same size as the      */
/*     original data set.                                         */
/* This is the forward FFT.                                       */

extern long long find_blocksize(long long n1, long long n2);

void realfft_scratch_fwd(multifile * infile, multifile * scratch, long long nn)
{
   long long n1, n2, bb, bb2, fp1, fp2, df, ii, jj, kk, kind;
   int i1, i2, move_size;
   unsigned char *move;
   rawtype *data, *dp;
   double tmp1, tmp2, h1r, h1i, h2r, h2i;
   double h2rwr, h2rwi, h2iwr, h2iwi;
   double wtemp, wpi, wpr, wi, wr, theta, delta;

   if (nn < 2)
      return;

   /* Treat the input data as a n1 (rows) x n2 (cols) */
   /* matrix.  Make sure that n2 >= n1.               */

   if (nn % 4 != 0) {
      printf("\nLength of FFT in twopassfft_real_fwd()\n");
      printf("   must be divisible by 4.\n\n");
      exit(1);
   }
   n2 = good_factor(nn / 4) * 2;
   if (n2 == 0) {
      printf("\nLength of FFT in twopassfft_real_fwd()\n");
      printf("   must be factorable\n\n");
      exit(1);
   }
   n1 = nn / (2 * n2);
   bb = find_blocksize(n1, n2);
   if (bb == 0 || bb % 2 || n1 % 2 || n2 % 2) {
      printf("\nCan't factor the FFT length in twopassfft_real_fwd()\n");
      printf("   into useful sizes.\n\n");
      exit(1);
   }

   /* First do n2 transforms of length n1 by  */
   /* fetching size bb x n1 blocks in memory. */

   data = gen_rawvect(bb * n2);

   /* Transpose scratch space */

   move_size = (bb + n2) / 2;
   move = (unsigned char *) malloc(move_size);
   move_size = (bb + n1) / 2;

   for (ii = 0; ii < n2; ii += bb) {

      /* Read a n1 (rows) x bb (cols) block of data */

      dp = data;
      fp1 = sizeof(rawtype) * ii;
      df = sizeof(rawtype) * n2;
      for (jj = 0; jj < n1; jj++) {
         fseek_multifile(infile, fp1, SEEK_SET);
         fread_multifile(dp, sizeof(rawtype), bb, infile);
         dp += bb;              /* Data ptr */
         fp1 += df;             /* File ptr */
      }

      /* Transpose the n1 (rows) x bb (cols) block of data */

      transpose_fcomplex(data, n1, bb, move, move_size);

      /* Do bb transforms of length n1 */

      for (jj = 0; jj < bb; jj++)
         COMPLEXFFT(data + jj * n1, n1, -1);

      /* Multiply the matrix A(ii,jj) by exp(- 2 pi i jj ii / nn). */
      /* Use recursion formulas from Numerical Recipes.            */

      for (jj = 0; jj < bb; jj++) {
         delta = -TWOPI * (ii + jj) / (double) (nn >> 1);
         wr = cos(delta);
         wi = sin(delta);
         wtemp = sin(0.5 * delta);
         wpr = -2.0 * wtemp * wtemp;
         wpi = wi;
         kind = jj * n1 + 1;
         for (kk = 1; kk < n1; kk++, kind++) {
            tmp1 = data[kind].r;
            tmp2 = data[kind].i;
            data[kind].r = tmp1 * wr - tmp2 * wi;
            data[kind].i = tmp2 * wr + tmp1 * wi;
            wtemp = wr;
            wr = wtemp * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
         }
      }
      fwrite_multifile(data, sizeof(rawtype), bb * n1, scratch);
   }

   /* Now do n1 transforms of length n2 by fetching      */
   /* groups of 2 size n2 (rows) x bb2 (cols) blocks and */
   /* then recombining blocks to make a real transform.  */

   bb2 = bb >> 1;

   /* Some values for the trig recursion below: */

   delta = -TWOPI / nn;
   wtemp = sin(0.5 * delta);
   wpr = -2.0 * wtemp * wtemp;
   wpi = sin(delta);

   /* Some Values for later */

   move_size = (bb + n2) / 2;
   df = sizeof(rawtype) * n1;

   for (ii = 0; ii < (n1 / 2); ii += bb2) {

      /* Read two n2 (rows) x bb2 (cols) blocks from the file */
      /* The first block comes from the start of the file and */
      /* the second block comes from the end.                 */
      /* Note:  The first block is shifted by one complex     */
      /*        value in order to make the complex->real      */
      /*        assembly quite a bit more elegant...          */

      dp = data;
      fp1 = sizeof(rawtype) * (ii + 1); /* File ptr */
      fp2 = sizeof(rawtype) * (n1 - ii - bb2);  /* File ptr */
      for (jj = 0; jj < n2; jj++) {
         fseek_multifile(scratch, fp1, SEEK_SET);
         fread_multifile(dp, sizeof(rawtype), bb2, scratch);
         dp += bb2;             /* Data ptr */
         fp1 += df;             /* File ptr */
         fseek_multifile(scratch, fp2, SEEK_SET);
         fread_multifile(dp, sizeof(rawtype), bb2, scratch);
         dp += bb2;             /* Data ptr */
         fp2 += df;             /* File ptr */
      }

      /* Transpose the n2 (rows) x bb (cols) block of data */

      transpose_fcomplex(data, n2, bb, move, move_size);

      /* Do bb transforms of length n2 */

      for (jj = 0; jj < bb; jj++)
         COMPLEXFFT(data + jj * n2, n2, -1);

      /* Transpose the bb (rows) x n2 (cols) block of data */

      transpose_fcomplex(data, bb, n2, move, move_size);

      /* Begin the re-assembly of the realFFT */

      for (jj = 0; jj < n2; jj++) {

         /* Start the trig recursion: */

         theta = (jj * n1 + ii + 1) * delta;
         wr = cos(theta);
         wi = sin(theta);

         /* Combine n and N/2-n terms as per Numerical Recipes. */

         i1 = jj * bb;          /* n     */
         i2 = bb * n2 - i1 - 1; /* N/2-n */
         for (kk = 0; kk < bb2; kk++, i1++, i2--) {
            h1r = 0.5 * (data[i1].r + data[i2].r);
            h1i = 0.5 * (data[i1].i - data[i2].i);
            h2r = 0.5 * (data[i1].i + data[i2].i);
            h2i = -0.5 * (data[i1].r - data[i2].r);
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

      /* Write two n2 (rows) x bb2 (cols) blocks to the file  */
      /* The first block goes to the start of the file and    */
      /* the second block goes to the end.                    */
      /* Note:  The first block is shifted by one complex     */
      /*        value in order to make the complex->real      */
      /*        assembly quite a bit more elegant...          */

      dp = data;
      fp1 = sizeof(rawtype) * (ii + 1); /* File ptr */
      fp2 = sizeof(rawtype) * (n1 - ii - bb2);  /* File ptr */
      for (jj = 0; jj < n2; jj++) {
         fseek_multifile(infile, fp1, SEEK_SET);
         fwrite_multifile(dp, sizeof(rawtype), bb2, infile);
         dp += bb2;             /* Data ptr */
         fp1 += df;             /* File ptr */
         fseek_multifile(infile, fp2, SEEK_SET);
         fwrite_multifile(dp, sizeof(rawtype), bb2, infile);
         dp += bb2;             /* Data ptr */
         fp2 += df;             /* File ptr */
      }
   }

   /* Now correct the final n2 values that have not been corrected */
   /*   due to the asymetry in the recombination.                  */

   /* Read the n2 data points */

   for (jj = 0; jj < n2; jj++) {
      fseek_multifile(scratch, sizeof(rawtype) * jj * n1, SEEK_SET);
      fread_multifile(data + jj, sizeof(rawtype), 1, scratch);
   }

   /* FFT the array: */

   COMPLEXFFT(data, n2, -1);

   /* Do the special cases of freq=0 and Nyquist freq */

   tmp1 = data[0].r;
   data[0].r = tmp1 + data[0].i;
   data[0].i = tmp1 - data[0].i;

   for (jj = 1, kk = n2 - 1; jj < n2 / 2; jj++, kk--) {
      theta = delta * n1 * jj;
      wr = cos(theta);
      wi = sin(theta);
      h1r = 0.5 * (data[jj].r + data[kk].r);
      h1i = 0.5 * (data[jj].i - data[kk].i);
      h2r = 0.5 * (data[jj].i + data[kk].i);
      h2i = -0.5 * (data[jj].r - data[kk].r);
      h2rwr = h2r * wr;
      h2rwi = h2r * wi;
      h2iwr = h2i * wr;
      h2iwi = h2i * wi;
      data[jj].r = h1r + h2rwr - h2iwi;
      data[jj].i = h1i + h2iwr + h2rwi;
      data[kk].r = h1r - h2rwr + h2iwi;
      data[kk].i = -h1i + h2iwr + h2rwi;
   }

   /* Write the n2 data points and clean up */

   for (jj = 0; jj < n2; jj++) {
      fseek_multifile(infile, sizeof(rawtype) * jj * n1, SEEK_SET);
      fwrite_multifile(data + jj, sizeof(rawtype), 1, infile);
   }
   free(move);
   vect_free(data);
}
