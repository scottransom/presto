#include "presto.h"
#include "plot2d.h"

int main (int argc, char *argv[])
{
  long i, N = 1048576, numfold = 0;
  double dt = 0.125, t = 0.0, f = 0.31234, fdot = -5.1e-10, phs = 0.234; 
  double prof[100], *freqs, onoff[2];
  float data[1048576];

  for (i = 0 ; i < N ; i++){
    t = i*dt;
    data[i] = 0.1234 * cos(TWOPI * (f + 0.5 * fdot * t) * t + phs);
  }
  for (i = 0 ; i < 100 ; i++){
    prof[i] = 0.0;
  }

  onoff[0] = 0;
  onoff[1] = 1023;
  numfold = 0;

  printf("Data is complete...\n");

  for (i = 0 ; i < 1024 ; i++){
    fold(data+i*1024, 1024, dt, i*1024*dt, prof, \
	 100, f, fdot, 0.0, 0, NULL, 0, 0, 0, onoff, &numfold);
  }

  cpgstart_x("landscape");
  
  freqs = gen_dfreqs(100, 0, 0.01);

  dxyline(100, freqs, prof, "Phase", "Relative Intensity", 1);

  cpgend();

  free(freqs);
  return 0;

}
