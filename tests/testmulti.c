#include "presto.h"
#include "plot2d.h"

double crab(double val);

int main(int argc, char *argv[])
{
  int proflen = 5, numchan = 60;
  long i, j;
  double *profs, *sumprof;
  float f = 1400, df = 0.25;

  profs = gen_dvect(proflen * numchan);
  sumprof = gen_dvect(proflen);

  for (i = 0 ; i < proflen ; i++){
/*     profs[i] = crab((double) i / (double) proflen); */
    profs[i] = 0.0;
    sumprof[i] = 0.0;
  }
  profs[2] = 1.0;
  profs[3] = 1.0;
  for (i = 1 ; i < numchan ; i++){
    printf("\n%f\n", 0.25 * i - 10.0);
    memcpy(&profs[i*proflen], profs, sizeof(double) * proflen);
    drotate(&profs[i*proflen], proflen, 0.25 * i - 10.0);
    for (j = 0 ; j < proflen ; j++){
      printf("%f ", profs[i*proflen+j]);
      sumprof[j] += profs[i*proflen+j];
    }
  }
  /*
  multi_prof_plot(proflen, numchan, profs, sumprof, f, df, \
		  "Pulse Phase (Periods)", "Channel Number", \
		  "Frequency (MHz)", "/XWIN");
  multi_prof_plot(proflen, numchan, profs, sumprof, f, df, \
		  "Pulse Phase (Periods)", "Channel Number", \
		  "Frequency (MHz)", "/VCPS");
  */
  free(profs);
  free(sumprof);

  return 0;
}



double crab(double val)
{
  /*  Normalized to area per pulse = 1 */

  static double phsval, normalize = 1.0 / 0.135, integral;
  static double phistart, m, b;

  phsval = modf(val, &integral);
  if ((phsval <= 0.15) || (phsval > 0.8))
    return (0.0);
  else if ((phsval > 0.35) && (phsval <= 0.55))
    return (0.74074074);
  else if (phsval <= 0.2) {
    phistart = 0.15;
    m = 2.0;
    b = 0.0;
  } else if (phsval <= 0.25) {
    phistart = 0.2;
    m = 18.0;
    b = 0.1;
  } else if (phsval <= 0.3) {
    phistart = 0.25;
    m = -16.0;
    b = 1.0;
  } else if (phsval <= 0.35) {
    phistart = 0.3;
    m = -2.0;
    b = 0.2;
  } else if (phsval <= 0.6) {
    phistart = 0.55;
    m = 2.0;
    b = 0.1;
  } else if (phsval <= 0.65) {
    phistart = 0.6;
    m = 4.0;
    b = 0.2;
  } else if (phsval <= 0.7) {
    phistart = 0.65;
    m = -4.0;
    b = 0.4;
  } else if (phsval <= 0.8) {
    phistart = 0.7;
    m = -2.0;
    b = 0.2;
  }
  return (normalize * ((phsval - phistart) * m + b));
}
