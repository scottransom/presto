#include "fftapps.h"
#include "plot2d.h"

#define NUMHARM 4

void main(int argc, char *argv[])
{
    FILE *file;
    float *data, *pow, *tmppow, *freqs, *locpow, *totstretch, nph;
    float *stretch[NUMHARM];
    long i, j, num=1, start, ct, pt, N, itmp=0;
    double maxpow=0.0, maxfreq=0.0, maxf=0.0, tmp, dt, T;
    char filenm[100], psfilenm[100];
    rderivs derivs;
    fourierprops props;
    infodata idata;

    if (argc==1) {
      printf("Usage: test_binary filename fourier_bin_# log2_fftsize [ps]\n\n");
      exit(1);
    }

    sprintf(filenm, "%s.fft",argv[1]);
    sprintf(psfilenm, "%s.ps",argv[1]);
    idata = readinf(argv[1]);
    file = chkfopen(filenm, "r+b");
    num <<= atoi(argv[3]);
    N = idata.numbins;
    dt = idata.bindt;
    T = N*dt;
    start = (int) (atof(argv[2]) - num/2);
    printf("Number of FFT bins %ld\n",num);
    if ((argc == 5) && (!strcmp(argv[4],"ps"))){
      cpgstart_ps(psfilenm);
    } else {
      cpgstart_x();
    }

    nph = get_numphotons(file);

    /* Initialize data arrays */

    data = vector(0, 2*num-1);
    pow = vector(0, num-1);
    tmppow = vector(0, num-1);
    freqs = vector(0, num-1);

    /* Get the raw FFT data */

    readfftarray(file, start, start+num-1, data);

    /* Turn the complex data into raw powers */

    for(ct=0, pt=0; ct<num; ct++, pt+=2) {
      freqs[ct] = (float) ct;
      pow[ct] = power(data[pt],data[pt+1]);
    }

    /* Do the Mini-FFT */

    realfft(pow,num,-1);

    /* Determine the local power values of the Mini-FFT */

    locpow = corr_loc_pow(pow,num/2);

    /* Plot the raw power */
    tmppow[0] = pow[0];
    tmppow[1] = pow[1];
    for(ct=2, pt=2; ct<num/2; ct++, pt+=2) {
      tmppow[ct] = power(pow[pt],pow[pt+1]);
    }
    xyline(num/2-2, freqs+2, tmppow+2, "Fourier Freq", "Power",1); 
    printf("Raw power.\n");

    /* Plot the raw power divided by local power */
    tmppow[0] = locpow[0];
    tmppow[1] = locpow[1];
    for(ct=0; ct<num/2; ct++) {
      tmppow[ct] /= locpow[ct];
    }
    xyline(num/2-2, freqs+2, tmppow+2, "Fourier Freq", "Pow/Local Pow",1); 
    printf("Raw Power/Local Power.\n");

    /* Correct the raw Mini-FFT to local-pow levels */

    pow[0] /= sqrt(locpow[0]);
    pow[1] /= sqrt(locpow[1]);
    for(ct=2, pt=2; ct<num/2; ct++, pt+=2) {
      tmp = 1.0/sqrt(locpow[ct]);
      pow[pt] *= tmp;
      pow[pt+1] *= tmp;
    }

    /* Generate interpolated harmonic arrays */

    totstretch = vector(0,num/2-1);
    for(ct=0; ct<num/2 ; ct++) totstretch[i]=tmppow[i];
    for(i=0 ; i<NUMHARM ; i++){
      stretch[i] = vector(0,num-1);
      stretch[i] = stretch_fft(pow,num/2,i+2, stretch[i]);
      for(ct=0, pt=0; ct<num/2; ct++, pt+=2) {
	tmppow[ct] = power(stretch[i][pt],stretch[i][pt+1]);
	totstretch[ct] += tmppow[ct];
	freqs[ct] = (float) (ct)/(i+2);
      }
      xyline(num/2, freqs, tmppow, "Fourier Freq", "Pow/Local Pow",1);
      printf("Stretch interpolation (stretch=%ld)\n",i+2);
    }


    /* Show the sum of the stretched values */
    for(ct=0 ; ct<num/2; ct++) {
       freqs[ct] = (float) (ct);
    }
    xyline(num/2, freqs, totstretch, "Fourier Freq", "Sum Loc Pow",1);
    printf("Stretch totals\n");

    cpgend();

    /*  resp = vector(0, 399);
	fmresp_dq(-100.8915, 1.0, 200, 1082.480541, 4.5,\
	DEGTORAD*(34.56789+90.0), 5, resp);
	for (i = 0, w = -100.8915; i < 200; i++, w += 1.0) {
	printf("%d resp[%7.3f] = %10.6f%+10.6fi\n",i, \
	w, resp[2 * i]*2317040.9, resp[2 * i + 1]*2317040.9);
	}
	free_vector(resp, 0, 399);
	*/

    for(i=0 ; i<NUMHARM ; i++){
      free_vector(stretch[i],0,num-1);
    }
    free_vector(totstretch, 0, num/2-1);
    free_vector(data, 0, 2*num-1);
    free_vector(freqs, 0, num-1);
    free_vector(pow, 0, num-1);
    free_vector(tmppow, 0, num-1);
    free_vector(locpow, 0, num-1);
    fclose(file);
    printf("\n\nDone.\n");
}
