#include "presto.h"

float *gen_comb(double roffset, double numteeth, double teethspace, \
		int numbetween)
{
  /* Generate a binary-like comb function to be used to correlate         */
  /* against a power spectrum.  (Note:  _power_, not complex amplitudes)  */
  /* A total of numteeth are generated, with the teeth centered around    */
  /* the data point (*m)/2 + roffset and separated by teethspace bins     */
  /* (Note:  bin != point -- see numbetween).  numbetween is the          */
  /* number of points that together equal one fourier bin (this is used   */
  /* when correlating with a fourier interpolated power spectrum)         */
  /* m is the # of pts returned on each side of the center point (as per  */
  /* the other gen_*_response() functions.  Note:  numteeth must be odd.  */

  long i, nptssinc, nptsresp, sincm;
  double x, dx;
  float *sincresp, *teethresp, *teetharray, *sincarray;


  /* Check that roffset and numbetween are in bounds */

  if (roffset < 0.0 || roffset >= 1.0) {
    printf("   roffset out of bounds in gen_comb().\n");
    exit(1);
  }
  if (numbetween < 1 || numbetween >= 20000) {
    printf("   numbetween out of bounds in gen_comb().\n");
    exit(1);
   }
  if (numteeth & 1) {
    printf("   numteeth must be odd in gen_comb().\n");
    exit(1);
  }

  /* Generate a sinc function kernel */

  sincm = 20;
  nptssinc = 2 * sincm * numbetween;
  x = - PI * sincm;
  dx = PI / numbetween;
  sincresp = gen_fvect(nptssinc);
  for (i = 0 ; i < nptssinc ; i++, x += dx){
    sincresp[i] = (fabs(x) > 0.00001) ? \
      sin(x)/x : 1.0 - x * x / 6.0;
  }

  /* Generate the teetharray */

  nptsresp = (long) floor((2 * sincm + (numteeth - 1) * teethspace
  teethresp = gen_fvect(nptsresp);
  for (i = 0, r = startr; i < 2 * (*m); i++, r += delta) {
    sinc = s / r;
    if (r != 0.0) sinc = s / r;
    else sinc = 1.0;
    s = alpha * s + beta * tmp + s;
  }

  /* Correct for divide by zero when the roffset is close to zero */

  if (roffset < 1E-3) {
    itmp = 2 * (*m);
    response[itmp] = 1 - 6.579736267392905746 * (tmp = roffset * roffset);
    response[itmp + 1] = roffset * (PI - 10.335425560099940058 * tmp);
  }
  return response;
}


