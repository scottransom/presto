/*    FFT Timing and Testing Program       */
/*          by Scott Ransom                */
/*            Version 3.0                  */
/*              11Sep98                    */

#include <time.h>
#include <sys/times.h>
#include <stdio.h>
#include <math.h>
#include "vectors.h"
#include "ransomfft.h"
#include "meminfo.h"
#include "randlib.h"
#include "clk_tck.h"

#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif

/* #define USERAWFFTW */

int main(int argc, char *argv[])
{
  float *data1, *data2;
  fcomplex *ptr1, *ptr2;
  long n, npts, tmp = 0, ct, plimit, prn = 0;
  long i, isign = -1;
  double err = 0.0;
#if defined USERAWFFTW
  FILE *wisdomfile;
  fftw_plan plan_forward, plan_inverse;
  static char wisdomfilenm[120];
#endif
  struct tms runtimes;
  double ttim, stim, utim, tott;
  
  if (argc <= 1 || argc > 4) {
    printf("\nUsage:  testffts [sign (1/-1)] [print (0/1)] [frac err tol]\n\n");
    exit(0);
  } else if (argc == 2) {
    isign = atoi(argv[1]);
    prn = 0;
    err = 0.02;
  } else if (argc == 3) {
    isign = atoi(argv[1]);
    prn = atoi(argv[2]);
    err = 0.02;
  }
  if (argc == 4) {
    isign = atoi(argv[1]);
    prn = atoi(argv[2]);
    err = atof(argv[3]);
  }

  /* import the wisdom for FFTW */
#if defined USERAWFFTW
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
#endif

  for (i = 0; i <= 8; i++) {
    
    /* npts = 1 << (i + 14);        # of points in FFT */
    /*      npts = 1 << 16;	 # of points in FFT */
    /*      npts = 4096;  	 # of points in FFT */
    /*      npts = 524288;   	 # of points in FFT */
    
    npts = 300000 * (i + 1);

    n = npts << 1;	       	/* # of float vals */
    
    data1 = gen_fvect(n);
    data2 = gen_fvect(n);
    ptr1 = (fcomplex *)data1;
    ptr2 = (fcomplex *)data2;
    
    /*      make the data = {1,1,1,1,-1,-1,-1,-1} (all real) */
    /*
      for (ct = 0; ct < npts/2; ct++) {
      tmp = 2 * ct;
      data1[tmp] = 1.0;
      data1[tmp + 1] = 0.0;
      data1[tmp + npts] = -1.0;
      data1[tmp + npts + 1] = 0.0;
      data2[tmp] = 1.0;
      data2[tmp + 1] = 0.0;
      data2[tmp + npts] = -1.0;
      data2[tmp + npts + 1] = 0.0;
      }
    */
    
    /*      make the data a sin wave of fourier freq 12.12345... */
    /*
      for (ct = 0; ct < npts; ct++) {
      tmp = 2 * ct;
      data1[tmp] = sin(2.0*3.14159265358979*ct*12.12345/npts)+1.0;
      data2[tmp] = data1[tmp];
      data1[tmp+1] = 0.0;
      data2[tmp+1] = data1[tmp+1];
      }
    */
    
    /*      make the data a sin wave of fourier freq 12.12345... with noise */
    
    for (ct = 0; ct < npts; ct++) {
      tmp = 2 * ct;
      data1[tmp] = 10.0 * sin(TWOPI * ct * 12.12345 / npts) + 100.0;
      data1[tmp] = gennor(data1[tmp], 10.0);
      data2[tmp] = data1[tmp];
      data1[tmp + 1] = gennor(100.0, 10.0);
      data2[tmp + 1] = data1[tmp + 1];
    }
    
    printf("\nCalculating...\n");
    
    /*  The challenger... */
    
    tott = times(&runtimes) / (double) CLK_TCK;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;

    tablesixstepfft(ptr1, npts, isign);
    /* tablesixstepfft(plan1, plan2, ptr1, npts, isign); */
    /*  sixstepfft(ptr1, npts, isign);       */
    /*  four1(ptr1 - 1, npts, isign);        */
    /*  tablefft(ptr1, npts, isign);         */
    /*  tablesplitfft(ptr1, npts, isign);    */
    /*  realfft(ptr1, n, isign);             */
    /*  fftw(plan, 1, in, 1, 0, out, 1, 0);  */
    
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    printf("Timing summary (Ransom)  npts = %ld:\n", npts);
    utim = runtimes.tms_utime / (double) CLK_TCK - utim;
    stim = runtimes.tms_stime / (double) CLK_TCK - stim;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);
    
    /*  The "Standard" FFT... */
    
    /* The following is for the fftw FFT */

    /* Create new plans */
#if defined USERAWFFTW
    plan_forward = fftw_create_plan(npts, -1, FFTW_MEASURE | \
                                           FFTW_USE_WISDOM | \
                                           FFTW_IN_PLACE);
    plan_inverse = fftw_create_plan(npts, +1, FFTW_MEASURE | \
                                           FFTW_USE_WISDOM | \
                                           FFTW_IN_PLACE);
#endif

    tott = times(&runtimes) / (double) CLK_TCK;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;

    /*  four1(ptr2 - 1, npts, isign);        */
    /*  tablefft(ptr2, npts, isign);         */
    /*  tablesplitfft(ptr1, npts, isign);    */
    /*  tablesixstepfft(ptr2, npts, isign);  */
    /*  realft(ptr2 - 1, n, isign);          */
    fftwcall(ptr2, npts, -1);

#if defined USERAWFFTW
    if (isign == -1) {
      fftw(plan_forward, 1, (FFTW_COMPLEX *) ptr2, 1, 1, NULL, 1, 1);
    } else {
      fftw(plan_inverse, 1, (FFTW_COMPLEX *) ptr2, 1, 1, NULL, 1, 1);
    }
#endif

    tott = times(&runtimes) / (double) CLK_TCK - tott;
    printf("Timing summary (FFTW)  npts = %ld:\n", npts);
    utim = runtimes.tms_utime / (double) CLK_TCK - utim;
    stim = runtimes.tms_stime / (double) CLK_TCK - stim;
    ttim = utim + stim;
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);
    
    /* The following is for the fftw FFT */

#if defined USERAWFFTW
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
#endif
        
    /* Check if correct with fractional errors... */
    
    for (ct = 0; ct < n; ct++) {
      if (data2[ct] != 0.0) {
	if (fabs((1.0 - (data1[ct] / data2[ct]))) > err) {
	  if ((ct % 2) == 1) {
	    printf("Values at freq %ld do not match to %4.2f%% fractional error:\n", (ct - 1) / 2, err * 100);
	    printf("  rl1 = %f  im1 = %f   rl2 = %f  im2 = %f\n",
		   data1[ct - 1], data1[ct], data2[ct - 1], data2[ct]);
	  } else {
	    printf("Values at freq %ld do not match to %4.2f%% fractional error:\n", ct / 2, err * 100);
	    printf("  rl1 = %f  im1 = %f   rl2 = %f  im2 = %f\n", data1[ct],
		   data1[ct + 1], data2[ct], data2[ct + 1]);
	  }
	}
      }
    }
    
    if (npts >= 64)
      plimit = 64;
    else
      plimit = npts;
    
    /* Print the output... */
    
    if (prn) {
      printf("\n   #1:  Challenger FFT...                      ");
      printf("#2:  Standard...\n");
      for (ct = 0; ct < plimit; ct++) {
	printf(" %3ld  rl = %12.3f   ", ct, data1[2 * ct]);
	printf("im = %12.3f    rl = %12.3f   im = %12.3f\n", \
	       data1[2 * ct + 1], data2[2 * ct], data2[2 * ct + 1]);
      }
    }

    free(data1);
    free(data2);
  }
  
  return 0;
  
}

