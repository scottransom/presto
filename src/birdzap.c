#include "presto.h"

int read_zapfile(char *zapfilenm, double **zapfreqs, double **zapwidths)
/* Open, read, and close a text file containing frequencies (Hz)  */
/* and widths (Hz) to ignore in a pulsar search.  The text file   */
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
  printf("Read %d freq/width pairs from the zapfile '%s'.\n\n", 
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
  static double max_freq = 0.0, current_freq = 0.0,  last_freq = 0.0;
  static int index = 1, firsttime = 1;

  if (firsttime){
    if (numzap < 2){
      printf("\n\n'numzap' = %d must be >= 2 in check_to_zap().", 
	     numzap);
      printf("  Exiting.\n\n");
      exit(1);
    }
    current_freq = zapfreqs[index] + 0.5 * zapwidths[index];
    last_freq = zapfreqs[index-1] + 0.5 * zapwidths[index-1];
    max_freq = zapfreqs[numzap-1] + 0.5 * zapwidths[numzap-1];
    firsttime = 0;
  }

  /* If we are beyond the end of the list, return a '0' */

  if (candfreq > max_freq) return 0;

  /* Shift our index so that we are pointing to the closest */
  /* 'birdie' above candfreq.                               */

  while (candfreq > current_freq){
    index++;
    last_freq = current_freq;
    current_freq = zapfreqs[index] + 0.5 * zapwidths[index];
  }

  /* Check the previous frequency as well as the current */
  /* frequency to see if they match.                     */

  if (((zapfreqs[index-1] - candfreq) < zapwidths[index-1]) ||
      ((zapfreqs[index] - candfreq) < zapwidths[index]))
    return 1;
  else return 0;
}
