
#include "fftapps.h"
#include "plot3d.h"
#include "plot2d.h"

int main(void)
{
  FILE *file;
  long i, j, m = 0, fftlen = 1024, rlo = 10000, numbetween = 5, npts = 32768;
  long nextbin;
  double ro = 10000.5, ravg, z = 0.6, r;
  double dt = 0.000030517578125, dz;
  double maxedpow, maxedr, maxedz;
  float *interp_new, *interp_new2, *data, *freqs, *response, *kernel,
  **plane;
  float rl, im, nph;
  rderivs derivs;
  fourierprops props;

  file = chkfopen("testresp.fft", "rb");
  nph = get_numphotons(file);

  ravg = ro + 0.5 * z;
  maxedz = 0.0;
  interp_new = corr_rz_interp_file(file, numbetween, rlo, z, fftlen, \
				   LOWACC, &m, &nextbin);

  kernel = gen_cvect(fftlen);
  response = gen_z_response(0.0, numbetween, z, LOWACC, &m);
  place_complex_kernel(response, m, fftlen, kernel);
  free(response);
  COMPLEXFFT(kernel, fftlen, -1);
  interp_new2 = corr_complex_raw_file(file, numbetween, rlo, \
				      kernel, fftlen, m, &nextbin);

  maxedpow = max_r_file(file, ravg, &maxedr, &derivs);
  printf("power = %f  r = %f\n", maxedpow, maxedr);
  props = calc_props(derivs, maxedr, maxedz, 0.0);
  print_candidate(&props, dt, npts, nph, 2);

  maxedpow = max_rz_file(file, ravg, z, &maxedr, &maxedz, &derivs);
  printf("power = %f  r = %f  z = %f\n", maxedpow, maxedr, maxedz);
  props = calc_props(derivs, maxedr, maxedz, 0.0);
  print_candidate(&props, dt, npts, nph, 2);

  /*plane = corr_rz_plane_file(file, numbetween, rlo, -2.0, 2.0, 3, fftlen, \
     LOWACC, &m, &nextbin); */

  for (i = 0; i < 7; i++) {
    r = ravg - 0.0015 + i * 0.0005;
    m = 0;
    rz_interp_file(file, r, z, LOWACC, &m, &rl, &im);
    printf("freq = %f   rl = %f  im = %f\n", r, rl, im);
  }

  /*  cpgstart_x();
     cpgask(1);
     for (i = 0 ; i < 5 ; i++){
     m=10000;
     dz=0.5;
     response = gen_z_response(0.0, m, i*dz, LOWACC, &m);
     data = gen_fvect(2*m);
     for (j = 0 ; j < 2*m ; j++){
     data[j] = power(response[2*j],response[2*j+1]);
     }
     free(response);
     freqs = gen_freqs(2*m, -1.0, 1.0/m);
     xyline(2*m, freqs, data, "Fourier Offset", "Phase (Rads)", 1);
     powerplot(2*m, freqs, data, 1.0, 1);
     }
     free(data);

     cpgend(); */

  /*  quick3d_complex(file, maxedr, "hue"); */

  freqs = gen_freqs(fftlen, rlo, 1.0 / numbetween);

  for (i = 0; i < 20; i++) {
    printf("%5d r = %6.3f   rl = %12.5f  im = %12.5f", i, freqs[i], \
	   interp_new[2 * i], interp_new[2 * i + 1]);
    printf("  rl = %12.5f  im = %12.5f\n", \
	   interp_new2[2 * i], interp_new2[2 * i + 1]);
    /*    printf("   pow = %12.5f  phs = %12.5f\n",\
       power(interp_new[2*i], interp_new[2*i+1]),\
       phase(interp_new[2*i], interp_new[2*i+1])); */
  }

  free(kernel);
  /* free(response);
     free(plane[0]);
     free(plane); */
  free(freqs);
  free(interp_new);
  free(interp_new2);
  fclose(file);
  return (1);
}
