#include "presto.h"

double delay_from_dm(double dm, double freq_emitted)
/* Return the delay in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, and the emitted    */
/* frequency (freq_emitted) of the pulsar in MHz.           */
{
  if (freq_emitted==0.0)
    return 0.0;
  else
    return dm / (0.000241 * freq_emitted * freq_emitted);
}


double dm_from_delay(double delay, double freq_emitted)
/* Return the Dispersion Measure in cm-3 pc, that would     */
/* cause a pulse emitted at frequency 'freq_emitted' to be  */
/* delayed by 'delay' seconds.                              */
{
  if (freq_emitted==0.0)
    return 0.0;
  else
    return delay * 0.000241 * freq_emitted * freq_emitted;
}


void dedisp(unsigned char *data, unsigned char *lastdata, long numpts,
	    double *dispdelays, long numchan, float *result)
/* De-disperse a stretch of data with numpts * numchan points. */
/* The delays (in bins) are in dispdelays for each channel.    */
/* The result is returned in result.  The input data and       */
/* dispdelays are always in ascending frequency order.         */
/* Data are ordered in time, with each channel in a row for    */
/* each time point.                                            */
{
  static int approx_mean, firsttime = 1, *offset;
  long ii, jj, kk;

  if (firsttime){
    offset = gen_ivect(numchan);
    for (ii = 0; ii < numchan; ii++)
      offset[ii] = (int) dispdelays[ii];
    approx_mean = -(numchan / 2 - 1);
    firsttime = 0;
  }

  /* Set the result array to negative of numchan / 2. */
  /* This will result in data with approx zero mean.  */

  for (ii = 0; ii < numpts; ii++)
    result[ii] = approx_mean;

  /* De-disperse */

  for (ii = 0; ii < numchan; ii++){
    jj = ii + offset[ii] * numchan;
    for (kk = 0; kk < numpts - offset[ii]; kk++, jj += numchan)
      result[kk] += lastdata[jj];
    jj = ii;
    for (; kk < numpts; kk++, jj += numchan)
      result[kk] += data[jj];
  }
}








