#include "presto.h"

int read_database(void)
{
  FILE *psrdatafile;
  int i, np, ct = 0, fsize;
  double per, ra, dec, width = 0.0;
  psrdatabase data;
  char jname[12];

  /* Open the binary data file */
  psrdatafile = chkfopen("bincat.dat", "rb");

  /* This is the length of the binary file in bytes */
  chkfread(&fsize, sizeof(int), 1, psrdatafile);
  printf("Data length = %d bytes\n", fsize);

  /* Read the pulsar data and the number of pulsars */
  chkfread(&data, sizeof(data), 1, psrdatafile);
  chkfread(&np, sizeof(int), 1, psrdatafile);
  printf("Number of pulsars in the data set = %d\n\n", np);

  for (i = 0; i < np; i++) {
    if (\
//        data.ntype[i]&1    &&   /* Globular Cluster */
    //      data.ntype[i]&2    &&   /* SNR */
    //      data.ntype[i]&4    &&   /* Glitches */
    //      data.ntype[i]&8    &&   /* Binary or multiple system */
	data.ntype[i] & 16 &&	/* ms Pulsar */
//      data.ntype[i]&32   &&   /* Recycled Pulsar */
    //      data.ntype[i]&64   &&   /* Interpulse */
    //      data.ntype[i]&128  &&   /* Optical, X-ray, or Gamma-Ray */
	1) {
      printf("%4d  PSR %.12s (PSR %.8s)  p_psr = %13.12f s\n", \
	     i + 1, data.jname + i * 12, data.bname + i * 8, data.w50[i]);
      if (fabs(data.w50[i]) > 0.001) {
	width += data.w50[i] / (1000.0 * data.p[i]);
	ct++;
      }
    }
  }
  printf("\nAvg pulse width at peak half-height is %5.4f units of phase.\n",
	 width / ct);
  printf("For %d pulsars.\n\n", ct);

  fclose(psrdatafile);
}
