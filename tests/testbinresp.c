#include "presto.h"
#include "plot2d.h"

int main(int argc, char** argv)
{
  FILE *datafile;
  char datafilenm[100];
  double N, T, dt, ppsr, apsr;
  double rint, z, roffset, qlo, dq, ro;
  int nq, tol=5;
  long i, numbetween = 3, m = 0, numout = 200, filelen;
  float *resp=NULL, *data=NULL, *freqs=NULL, norm, nph;
  struct tms runtimes;
  double ttim, stim, utim, tott;
  orbitparams orbit;
  makedata mdata;
  /* double ti, tr, q; */

  if (argc != 2){
    printf("usage:  testbinresp filename\n\n");
    exit(0);
  }

  printf("\n     Binary Response Test Program\n");
  printf("           by Scott M. Ransom\n");
  printf("              7 March 1998\n\n");

  /* Determine the filenames to use */

  sprintf(datafilenm, "%s.fft", argv[1]);
  sprintf(makefilenm, "%s.mak", argv[1]);

  /* read the makefile */

  read_mak_file(argv[1], &mdata);

/*
orbit.t = 0.0;
orbit.p = 36075.5;
orbit.x = 0.55;
orbit.w = 230.0;
*/


  /* Calculate some of the parameters */

  T = N * dt;
  norm = apsr * apsr * N * N / 4.0;
  z = TWOPI * orbit.x / ppsr;
  ro = T / ppsr;
  roffset = modf(ro, &rint);

  printf("Integration time  = %f sec\n",T);
  printf("\nPulsar period     = %f sec\n", ppsr);
  printf("Pulsar FFT freq   = %f cycles\n", ro);
  printf("Pulsar amplitude  = %f \n", apsr);
  printf("Response roffset  = %f cycles\n", roffset);
  printf("\nBinary period     = %f sec\n", orbit.p);
  printf("Binary FFT freq   = %f cycles\n", T / orbit.p);
  printf("Binary a*sin(i)/c = %f lt-sec\n", orbit.x);
  printf("Binary e          = %f\n", orbit.e);
  printf("Binary w          = %f deg\n", orbit.w);
  printf("Binary To         = %f sec\n", orbit.t);
  printf("\nModulation z      = %f\n", z);


  /* Generate the theoretical response */

  printf("\n\nGenerating the theoretical response...\n");

  tott = times(&runtimes) / (double) CLOCKS_PER_SEC;

  roffset = 0.0;

  resp = gen_bin_response(roffset, numbetween, ppsr, T, \
			  &orbit, tol, &m);
/*   resp = gen_r_response(roffset, numbetween, LOWACC, &m);  */

  qlo = m / (double) numbetween + roffset;
  dq = -1.0 / (double) numbetween;
  nq = 2 * m;
  printf("m = %ld\n",m);
  
  printf("Timing summary for theoretical response creation:\n");
  tott = times(&runtimes) / (double) CLOCKS_PER_SEC - tott;
  utim = runtimes.tms_utime / (double) CLOCKS_PER_SEC;
  stim = runtimes.tms_stime / (double) CLOCKS_PER_SEC;
  ttim = utim + stim;
  printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	 ttim, utim, stim);
  printf("Total time elapsed:  %.3f sec.\n\n", tott);

  /* Start PGPLOT */

  cpgstart_x("landscape");

  /* Plot the theoretical response */

  freqs = gen_freqs(nq, qlo, dq);
  powerplot(nq, freqs, resp, 1.0, 1);
  printf("The theoretical response powers.\n");
  free(freqs);

  /* Open and prep the FFT file */

  datafile = chkfopen(datafilenm, "rb");
  filelen = get_filelen(datafile, sizeof(float));
  nph = get_numphotons(datafile);
/*   data = gen_cvect(nq); */
/*   read_fcomplex_file(datafile, ro - qlo, ro - qlo - nq * dq - 1, data); */

  /* Correlate... */

  data = gen_cvect(numout*numbetween);
  corr_complex_file(datafile, filelen, (long) floor(ro-numout/2), \
		    numout, numbetween, resp, m, data, 0, 0);

  /* Plot the data */


  freqs = gen_freqs(numout*numbetween, (float) -numout/2, \
		    1.0 / (float) numbetween);
  powerplot(numout*numbetween, freqs, data, norm, 1);
  free(freqs);

/*
  freqs = gen_freqs(nq, qlo, dq);
  powerplot(nq, freqs, data, norm, 1);
  free(freqs);
*/

  printf("The measured response powers.\n");

  /* Print the complex values */

  for (i=0 ; i<2*nq ; i+=2){
    /*    printf("fft = %10.7f + %10.7fi,   theo = %10.7f + %10.7fi\n", \
	  data[i]/norm, data[i+1]/norm, resp[i], resp[i+1]);*/
/*
    printf("fft: pow =  %13.11f  phs = %9.4f,  ", \
 	   POWER(data[i], data[i+1]) / norm, \
 	   PHASE(data[i], data[i+1]));
    printf("theo: pow = %13.11f  phs = %9.4f\n", \
	   POWER(resp[i], resp[i+1]), \
 	   PHASE(resp[i], resp[i+1]));
*/
    
  }

  /* Cleanup */

  fclose(datafile);
  free(resp);
  free(data);
  if (mdata.onoff) free(mdata.onoff);
  while (getchar() == '\n');
  return(0);
}

