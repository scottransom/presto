#include "presto.h"

void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);

int read_zapfile(char *zapfilenm, double **zapfreqs, double **zapwidths)
/* Open, read, and close a text file containing frequencies (Hz)  */
/* and widths (bins) to ignore in a pulsar search.  The text file */
/* should have one frequency and width per line.  Lines beginning */
/* with '#' are ignored, and so may be used as comments.          */
{
  FILE *zapfile;
  char line[200];
  int ii, numzap;

  zapfile = chkfopen(zapfilenm, "r");

  /* Read the input file once to count TOAs */
  
  numzap = 0;
  while (!feof(zapfile)){
    fgets(line, 200, zapfile);
    if (line[0]=='#') continue;
    else numzap++;
  }
  numzap--;

  /* Allocate the birdie arrays */

  *zapfreqs = gen_dvect(numzap);
  *zapwidths = gen_dvect(numzap);

  /* Rewind and read the TOAs for real */

  rewind(zapfile);
  ii = 0;
  while(ii < numzap){
    fgets(line, 200, zapfile);
    if (line[0]=='#') continue;
    else {
      sscanf(line, "%lf %lf\n", &(*zapfreqs)[ii], &(*zapwidths)[ii]);
      ii++;
    }
  }

  /* Close the file and exit */

  fclose(zapfile);
  printf("Read %d 'birdie' pairs from '%s'.\n\n", 
	 numzap, zapfilenm);
  return numzap;
}


int check_to_zap(double candfreq, double *zapfreqs, double *zapwidths, 
		 int numzap)
/* Look at the closest birdies in the zapfile to see if our candidate  */
/* matches one of them.  If it does, return '1' for TRUE.  If it       */
/* doesn't match, return a '0' for FALSE.  Note that the zapfreqs      */
/* _must_ be in increasing order since this routine keeps track of its */
/* place in the file.  Also, numzap _must be >= 2.                     */
{
  static unsigned long index=0;

  if (numzap < 2){
    printf("\n\n'numzap' = %d must be >= 2 in check_to_zap().", 
	   numzap);
    printf("  Exiting.\n\n");
    exit(1);
  }

  /* If we are beyond the end of the list, return a '0' */

  if (candfreq > zapfreqs[numzap-1] + 0.5 * zapwidths[numzap-1]) 
    return 0;

  /* Find the indices for the birdies above and below our candidate */

  index++;
  hunt(zapfreqs-1, numzap, candfreq, &index);
  index--;

  /* Check the lower and high birdie freqs to see if they match. */
  
  if ((fabs(zapfreqs[index] - candfreq) < 0.5 * zapwidths[index]) ||
      (fabs(zapfreqs[index+1] - candfreq) < 0.5 * zapwidths[index+1]))
    return 1;
  else return 0;
}


void hunt(double *xx, unsigned long n, double x, unsigned long *jlo)
{
  unsigned long jm,jhi,inc;
  int ascnd;
  
  ascnd=(xx[n] >= xx[1]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if ((x >= xx[*jlo]) == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while ((x >= xx[jhi]) == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo)--;
      while ((x < xx[*jlo]) == ascnd) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if ((x >= xx[jm]) == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
  if (x == xx[n]) *jlo=n-1;
  if (x == xx[1]) *jlo=1;
}


