#include "fftapps.h"
#include "plotting.h"

/*    Example main program to test FFT applications.       */
/*    Requires files test_finterp.fft, test_fdot.fft,      */
/*               and test_plot.fft                         */


void main(int argc, char *argv[])
{
    /* Theoretical values from test_finterp.mak and test_fdot.mak: */

    unsigned long npts = 131072;
    double dt = 7.62939453125e-06;
    double theo_freq = 2345.75;
    double theo_fdot = 23.456789;
    double theo_fint_pow = 4.294967295e+09;

    /* Adjustable parameters: */
    int fftlen = 16384;
    int numbetween = 8;
    int numf = 71;
    int numfdot = 71;
    double df = 0.25;
    double dfdot = 1.5;
    char notes[20];

    /* Don't adjust these: */
    FILE *file;
    char palette[40], go[20];
    double time1, time2, time3;
    double nph, f, w, sincfac, startf, timefact;
    double lof, hif, lofdot, hifdot, temp, fdot;
    double maxedpow, maxedr, maxedz;
    clock_t start, finish;
    rderivs derivs;
    fourierprops props, cands[2];
    float *theo_vec, *reg_vec, *corr_vec, **reg_mat, **corr_mat;
    float rl, im, *resp;
    int ct, pt, numintpts, numcorr;
    int numdisplay = 10, offset = 45.0;		/* params for interp */

    unsigned long nextbin, lobin;
    long numbet, r;

/*
   //This is a test for the nice_output_*() functions

   double val,err;
   int i,j;
   char output1[40], output2[40];


   for (i=0;i!=15;i++){
   val = -1.23456789000000000000*pow(10,7-i);
   for (j=0;j!=15;j++){
   err = 9.9567000000000000000*pow(10,7-j);
   printf("val=%20.10f  err=%20.10f ",val,err);
   nice_output_1(output1,val,err,20);
   printf("%s  ",output1);
   nice_output_2(output2,val,err,20);
   printf("%s  \n",output2);
   }
   printf("\n");
   }

   exit(0);

 */
    
    if (argc > 1) {
	strcpy(palette, argv[1]);
    } else {			/*  If no palette, set to gray. */
	strcpy(palette, "gray");
    }

    printf("\n\nFirst check the fourier interpolation:\n");
    printf("    The following test data is from the file\n");
    printf("       test_finterp.fft.\n\n");

    file = chkfopen("test_finterp.fft", "r+b");
    nph = get_numphotons(file);
    numintpts = fftlen - 2 * numbetween * NUMFINTBINS;
    theo_vec = vector(0, numintpts * 2 - 1);
    startf = 2300.0;
    for (ct = 0, f = startf, pt = 0; ct < numintpts;
	 ct++, f += 1.0 / numbetween, pt += 2) {
	w = PI * (theo_freq - f);
	if (w == 0.0)
	    sincfac = 0.5 * npts;
	else
	    sincfac = 0.5 * npts * sin(w) / w;
	theo_vec[pt] = sincfac * cos(w - PIBYTWO);
	theo_vec[pt + 1] = sincfac * sin(w - PIBYTWO);
    }

    start = clock();
    reg_vec = finterp_file(file, startf,
	    startf + ((double) (numintpts - 1)) / numbetween, numintpts);
    finish = clock();
    time1 = (double) (finish - start) / CLOCKS_PER_SEC;
    printf("    Standard brute force finterp took %f secs for %d points.\n",
	   time1, numintpts);

    start = clock();
    lobin = (unsigned long) (startf + DBLCORRECT);
    corr_vec = corr_file_finterp(file, lobin, numbetween, fftlen, &nextbin);
    finish = clock();
    time2 = (double) (finish - start) / CLOCKS_PER_SEC;
    numcorr = (nextbin - lobin) * numbetween;
    printf("  Correlation finterp (w/ setup) took %f secs for %d points.\n",
	   time2, numcorr);

    start = clock();
    lobin = (unsigned long) (startf + DBLCORRECT);
    corr_vec = corr_file_finterp(file, lobin, numbetween, fftlen, &nextbin);
    finish = clock();
    time3 = (double) (finish - start) / CLOCKS_PER_SEC;
    printf("  Correlation finterp (no setup) took %f secs for %d points.\n",
	   time3, numcorr);

    timefact = time1 * numcorr / (double) (time3 * numintpts);
    printf("\n\tResult:  Correlation method is over %6.2f times faster.\n\n",
	   timefact);

    printf("    Freq        Theo: Pow   Phs     Slow:");
    printf(" Pow   Phs     Fast: Pow   Phs\n");
    printf("-----------------------------------------");
    printf("-------------------------------\n");
    for (ct = 0; ct < numdisplay; ct++) {
	pt = 2 * (offset * numbetween + ct);
	f = startf + offset + ct / (double) numbetween;
	printf("%12.6f    %8.2f  %6.2f    ", f,
	       power(theo_vec[pt], theo_vec[pt + 1]) / nph,
	       phase(theo_vec[pt], theo_vec[pt + 1]));
	printf("%8.2f  %6.2f    ",
	       power(reg_vec[pt], reg_vec[pt + 1]) / nph,
	       phase(reg_vec[pt], reg_vec[pt + 1]));
	printf("%8.2f  %6.2f\n",
	       power(corr_vec[pt], corr_vec[pt + 1]) / nph,
	       phase(corr_vec[pt], corr_vec[pt + 1]));
    }
    free_vector(theo_vec, 0, numintpts * 2 - 1);
    free_vector(reg_vec, 0, numintpts * 2 - 1);

    nph = get_numphotons(file);
    printf("\n");
    printf("Check out the maximization and candidate property routines:\n\n");
    maxedpow = max2d(file, 2346.0, 0.0, &maxedr, &maxedz, &derivs);
    props = calc_props(derivs, maxedr, maxedz, 0.0);
    print_candidate(&props, dt, npts, nph, 2);

    cands[0] = props;

    fclose(file);

    printf("\n\nNow check the f/fdot plane generator:\n");
    printf("    The following test data is from the file\n");
    printf("       test_plot.fft.\n\n");

    file = chkfopen("test_plot.fft", "r+b");
    nph = get_numphotons(file);
    temp = (numf - 1) * df * 0.5;
    lof = floor(theo_freq - temp + DBLCORRECT);
    hif = lof + 2 * temp;
    temp = (numfdot - 1) * dfdot * 0.5;
    lofdot = theo_fdot - temp;
    hifdot = lofdot + 2 * temp;

    /*  Check maximization procedures: */
    printf("Check out the maximization and candidate property routines:\n\n");
    file_readffdot(file, theo_freq, theo_fdot, 1.0, &rl, &im);
    printf("Power at the theoretical f/fdot = %f\n", power(rl, im));
    printf("              Theoretical power = %f\n", theo_fint_pow);
    printf("        Total number of photons = %f\n", nph);
    printf("                         Amp rl = %f\n", rl);
    printf("                         Amp im = %f\n", im);
    printf("                  Theo max freq = %f\n", theo_freq);
    printf("                  Theo max fdot = %f\n\n", theo_fdot);
    maxedpow = max2d(file, (long) theo_freq, (long) theo_fdot, &maxedr, \
		     &maxedz, &derivs);
    props = calc_props(derivs, maxedr, maxedz, 0.0);
    print_candidate(&props, dt, npts, nph, 2);

    cands[1] = props;
    file_reg_candidates(cands, notes, 2, dt, npts, nph, "testfft");

    reg_mat = file_ffdot_plane(file, theo_freq - 2, theo_freq + 2, 41,
			       theo_fdot - 16, theo_fdot + 16, 41, 1.0);
    plot3d_complex(reg_mat, 41, 41, palette);

    /*    exit(0);  */

    start = clock();
    reg_mat = file_ffdot_plane(file, lof, hif, numf,
			       lofdot, hifdot, numfdot, 1.0);
    finish = clock();
    time1 = (double) (finish - start) / CLOCKS_PER_SEC;
    printf("    Standard brute force f/fdot took %f secs for %d points.\n",
	   time1, numf * numfdot);

    fftlen /= 4;
    lobin = (unsigned long) (lof + DBLCORRECT);
    numbet = (long) (1 / df + DBLCORRECT);
    start = clock();
    corr_mat = corr_file_ffdot(file, lobin, numbet, lofdot, hifdot,
			       numfdot, 1.0, fftlen, &nextbin);
    finish = clock();
    time2 = (double) (finish - start) / CLOCKS_PER_SEC;
    numcorr = numfdot * (nextbin - lobin) * numbet;
    printf("  Correlation f/fdot (w/ setup) took %f secs for %d points.\n",
	   time2, numcorr);

    start = clock();
    corr_mat = corr_file_ffdot(file, lobin, numbet, lofdot, hifdot,
			       numfdot, 1.0, fftlen, &nextbin);
    finish = clock();
    time3 = (double) (finish - start) / CLOCKS_PER_SEC;
    printf("  Correlation f/fdot (no setup) took %f secs for %d points.\n",
	   time3, numcorr);

    timefact = time1 * numcorr / (double) (time3 * numf * numfdot);
    printf("\n\tResult:  Correlation method is over %6.2f times faster.\n\n",
	   timefact);


    file_readffdot(file, theo_freq, theo_fdot, 1.0, &rl, &im);
    printf("Power at the theoretical f/fdot = %f\n", power(rl, im) / nph);
    printf("        Total number of photons = %f\n", nph);
    printf("                  Theo max freq = %f\n", theo_freq);
    printf("                  Theo max fdot = %f\n\n", theo_fdot);

    printf("    Freq        Fdot       Slow:");
    printf(" Pow   Phs     Fast: Pow   Phs\n");
    printf("-----------------------------------------");
    printf("-------------------------------\n");
    r = numf - 1;
    f = lof + df * 0.5 * (numf - 1);
    for (ct = 0; ct < numdisplay; ct++) {
	pt = (numfdot - 1) / 2 - numdisplay / 2 + ct;
	fdot = theo_fdot + (ct - 0.5 * numdisplay) * dfdot;
	printf("%12.6f  %8.3f     ", f, fdot);
	printf("%8.2f  %6.2f    ",
	       power(reg_mat[pt][r], reg_mat[pt][r + 1]) / nph,
	       phase(reg_mat[pt][r], reg_mat[pt][r + 1]));
	printf("%8.2f  %6.2f\n",
	       power(corr_mat[pt][r], corr_mat[pt][r + 1]) / nph,
	       phase(corr_mat[pt][r], corr_mat[pt][r + 1]));
    }
    fclose(file);

    printf("\n\nTesting 3-D plotting routines:  ");
    printf("(Press 'c' and enter to change plots.)\n\n");
    printf("     Theoretical response data from slow method.\n");
    plot3d_complex(reg_mat, numfdot, numf, palette);

    file = chkfopen("test_plot.fft", "r+b");
    scanf("%c", go);
    printf("     Noisy \"pulsar\" data using fast method \"quick plot\".\n");
    quick3d_complex(file, theo_freq, palette);
    fclose(file);

    printf("\n\nDone.\n");
    cleanup();
}
