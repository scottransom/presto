#include "fftapps.h"
#include "plot2d.h"

void readmake(FILE * makefile, char *makefilenm, double *N, \
              double *dt, double *ppsr, double *apsr, double *p_orb, \
	      double *x_orb, double *e_orb, double *t_orb, double *w_orb);

int main(int argc, char** argv)
{
  FILE *makefile, *datafile;
  char makefilenm[100], datafilenm[100];
  double N, T, dt, qlo, dq, ppsr, apsr;
  double pb, z, e, x, w, to, ro, rint;
  int nq, tol=5;
  long i;
  float powargr, powargi;
  float phsargr, phsargi, phstmp;
  float *resp=NULL, *data=NULL, *freqs=NULL, norm, nph;
  float dpow1, dphs1, dpow2, dphs2;
  struct tms runtimes;
  double ttim, stim, utim, tott;
  orbitparams orbit;
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

  makefile = chkfopen(makefilenm, "r");
  readmake(makefile, makefilenm, &N, &dt, &ppsr, &apsr, &pb, &x, &e, &to, &w);
  fclose(makefile);

  /* Calculate some of the parameters */

  T = N * dt;
  z = TWOPI * x / ppsr;
  nq = (int) floor(T / pb * z * 2.2 + 0.5);

  if (nq & 1) nq++;
  ro = T / ppsr;
  qlo = (modf(ro, &rint) + nq / 2.0);
  dq = -1.0;

  orbit.p = pb;
  orbit.x = x;
  orbit.w = w / DEGTORAD;
  orbit.e = e;
  orbit.t = to;

  norm = apsr * apsr * N * N / 4.0;

  printf("Integration time  = %f sec\n",T);
  printf("\nPulsar period     = %f sec\n", ppsr);
  printf("Pulsar FFT freq   = %f cycles\n", T / ppsr);
  printf("Pulsar amplitude  = %f \n", apsr);
  printf("\nBinary period     = %f sec\n", pb);
  printf("Binary FFT freq   = %f cycles\n", T / pb);
  printf("Binary a*sin(i)/c = %f lt-sec\n", x);
  printf("Binary e          = %f\n", e);
  printf("Binary w          = %f radians\n", w);
  printf("Binary To         = %f sec\n", to);
  printf("\nModulation z      = %f\n", z);
  printf("Modulation phi    = %f\n", w);

  /* Generate the theoretical response */

  printf("\n\nGenerating the theoretical response...\n");
  resp = gen_fvect(2*nq);

  tott = times(&runtimes) / (double) CLK_TCK;

  printf("qlo = %f, dq = %f, lobin(r) = %f\n", qlo, dq, ro-qlo);
  
  fmresp_dq(qlo, dq, nq, z, T / pb, w, tol, resp);
  /*
  for (i = 0, q = qlo ; i < nq ; i++, q += dq){
    fmresp(q, z, T / pb, w, tol, &tr, &ti);
    resp[2*i] = (float) tr;
    resp[2*i+1] = (float) ti;
  }
  */
  printf("Timing summary for theoretical response creation:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	 ttim, utim, stim);
  printf("Total time elapsed:  %.3f sec.\n\n", tott);

  /* Start PGPLOT */

  cpgstart_x();

  /* Plot the theoretical response */

  dpow2 = POWER(resp[0], resp[1]);
  dphs2 = PHASE(resp[0], resp[1]);
  freqs = gen_freqs(nq, qlo, dq);
  powerplot(nq, freqs, resp, 1.0, 1);
  printf("The theoretical response powers.\n");

  /* Read the data from the FFT file */

  datafile = chkfopen(datafilenm, "rb");
  nph = get_numphotons(datafile);
  data = gen_fvect(2*nq);
  readfftarray(datafile, ro - qlo, ro - qlo - nq * dq - 1, data);
  fclose(datafile);
  dpow1 = POWER(data[0], data[1]);
  dphs1 = PHASE(data[0], data[1]);

  /* Plot the data */

  powerplot(nq, freqs, data, norm, 1);
  free(freqs);
  printf("The measured response powers.\n");

  /* Print the complex values */

  for (i=0 ; i<2*nq ; i+=2){
    /*    printf("fft = %10.7f + %10.7fi,   theo = %10.7f + %10.7fi\n", \
	   data[i]/norm, data[i+1]/norm, resp[i], resp[i+1]);*/
    printf("fft: pow =  %13.11f  phs = %9.4f,  theo: pow = %13.11f  phs = %9.4f\n", \
	   POWER(data[i], data[i+1]) / norm,
	   PHASE(data[i], data[i+1]),
	   POWER(resp[i], resp[i+1]),
	   PHASE(resp[i], resp[i+1]));
  }

  free(resp);
  free(data);
  while (getchar() == '\n');
  return(0);
}


void readmake(FILE * makefile, char *makefilenm, double *N, \
              double *dt, double *ppsr, double *apsr, double *p_orb, \
	      double *x_orb, double *e_orb, double *t_orb, double *w_orb)
{
  double f_psr;
  char temp[200];

  fscanf(makefile, "%[^\n]", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *N = atof(temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *dt = atof(temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  f_psr = atof(temp);
  *ppsr = 1.0 / f_psr;
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *apsr = atof(temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *p_orb = atof(temp);
  if (*p_orb == 0.0) {
    printf("The data defined in %s is not a binary.  Exiting\n\n", makefilenm);
    exit(1);
  }
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *x_orb = atof(temp);
  if (*x_orb == 0.0) {
    printf("The data defined in %s is not a binary.  Exiting\n\n", makefilenm);
    exit(1);
  }
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *e_orb = atof(temp);
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *w_orb = atof(temp) * DEGTORAD;
  fscanf(makefile, "%*[^=] %*s %s", temp);
  *t_orb = atof(temp);
  *t_orb = fmod(*t_orb, *p_orb);
}
