#include "fftapps.h"
#include "plot2d.h"

void main(int argc, char *argv[])
{
    FILE *file;
    /*    FILE *file2;*/
    float *data, *pow, *pow2, *freqs, *interp, nph;
    int num=1, start, ct, pt, lobin, N, itmp=0;
    double maxpow, maxfreq=0.0, maxf=0.0, tmp, dt, T;
    char filenm[100], psfilenm[100];
    rderivs derivs;
    fourierprops props;
    infodata idata;

    if (argc==1) {
      printf("Usage: tester filename fourier_bin_# log2_fftsize [ps]\n\n");
      exit(1);
    }

    sprintf(filenm, "%s.fft",argv[1]);
    sprintf(psfilenm, "%s.ps",argv[1]);
    idata = readinf(argv[1]);
    file = chkfopen(filenm, "r+b");
    num <<= atoi(argv[3]);;
    N = idata.numbins;
    dt = idata.bindt;
    T = N*dt;
    lobin = 0;
    start = (int) (atof(argv[2]) - num/2);
    printf("Number of FFT bins %d\n",num);
    if ((argc == 5) && (!strcmp(argv[4],"ps"))){
      cpgstart_ps(psfilenm);
    } else {
      cpgstart_x();
    }

    nph = get_numphotons(file);

    data = vector(0, 2*num-1);

    readfftarray(file, start, start+num-1, data);

    pow = vector(0, num-1);
    freqs = vector(0, num-1);

    for(ct=0, pt=0; ct<num; ct++, pt+=2) {
      freqs[ct] = (float) ct;
      /* Convert complex amplitudes into normalized power spectrum */
      pow[ct] = power(data[pt],data[pt+1])/nph;
    }

    /* Plot original power spectrum */
    powerplot(num, freqs, data, nph, 1);

    /* Use realfft on power spectrum */
    realfft(pow,num,-1);

    for(ct=1000, pt=2000, maxpow = 0.0; ct < num/2; ct++, pt+=2) {
      tmp = power(pow[pt],pow[pt+1]);
      if (tmp > maxpow) {
	maxpow = tmp;
	maxf = ct;
      }
    }

    /* Test maximize1d functions */
    /*    file2 = chkfopen("tempfile.fft", "wb+");
    chkfwrite(pow,sizeof(float),num,file2);
    */

    maxpow = max1d_arr(pow,num,maxf,&maxfreq,&derivs);
    props = calc_props(derivs, maxfreq, 0.0, 0.0);
    print_bin_candidate(&props, idata.bindt, N, pow[0],\
			num, start+lobin, 2);

    /*    maxpow = max1d_file(file2,maxf,&maxfreq,&derivs);
    props = calc_props(derivs, maxfreq, 0.0, 0.0);
    print_bin_candidate(&props, idata.bindt, N, pow[0],\
			num, start+lobin, 2);
    */

    /* Plot the resulting power spectrum */

    pow2 = vector(0, num/2 - 1);
    for (ct = 0 ; ct < num/2 ; ct++){
      freqs[ct] = (float) ct*T/num;
      pow2[ct] = power(pow[(itmp = 2*ct)],pow[itmp+1])/pow[0];
    }
    
    xyline(num/2-2, freqs+2, pow2+2, "Binary Period (sec)", "Power/Freq 0 Power",1); 

    /* Interpolate and plot the peak */

    start = maxf-50.0; 
    interp = finterp_array(pow, start, start+100, 3000); 
    for(ct = 0; ct < 3000; ct++) { 
      freqs[ct] = (ct*(100.0/3000.0)+start)*T/num;
      pow2[ct] = power(interp[2*ct],interp[2*ct+1])/pow[0];
    }
    xyline(3000, freqs, pow2, "Binary Period (sec)", "Power/Freq 0 Power",1); 
    
    cpgend();
    
    free_vector(data, 0, 2*num-1);
    free_vector(freqs, 0, num-1);
    free_vector(pow, 0, num-1);
    free_vector(pow2, 0, num/2-1);
    free_vector(interp, 0, 5999);
    fclose(file);
    /*fclose(file2);*/
    cleanup();
    printf("\n\nDone.\n");
}







