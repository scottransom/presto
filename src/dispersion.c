#include "presto.h"

double delay_from_dm(double dm, double freq_emitted)
/* Return the delay in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, and the emitted    */
/* frequency (freq_emitted) of the pulsar in MHz.           */
{
  return dm / (0.000241 * freq_emitted * freq_emitted);
}


double dm_from_delay(double delay, double freq_emitted)
/* Return the Dispersion Measure in cm-3 pc, that would     */
/* cause a pulse emitted at frequency 'freq_emitted' to be  */
/* delayed by 'delay' seconds.                              */
{
  return delay * 0.000241 * freq_emitted * freq_emitted;
}


void dedisp(float *data, float *lastdata, long numpts, \
	    double *dispdelays, long numchan, float *result)
/* De-disperse a stretch of data with numpts * numchan points. */
/* The delays (in bins) are in dispdelays for each channel.    */
/* The result is returned in result.  The input data and       */
/* dispdelays are always in ascending frequency order.         */
/* Data are ordered in time, with each channel in a row for    */
/* each time point.                                            */
{
  static double *lofrac, *hifrac;
  static float *tmpdat;
  static int firsttime = 1, *offset;
  long ii, jj, ptr, tmpptr;

  if (firsttime) {
    lofrac = gen_dvect(numchan);
    hifrac = gen_dvect(numchan);
    offset = gen_ivect(numchan);
    tmpdat = gen_fvect(numpts + 1);
    for (ii = 0; ii < numchan; ii++) {
      offset[ii] = (int) floor(dispdelays[ii]);
      lofrac[ii] = dispdelays[ii] - offset[ii];
      hifrac[ii] = 1.0 - lofrac[ii];
    }
    firsttime = 0;
  }

  /* Reset the result array to 0's */

  for (ii = 0; ii < numpts; ii++)
    result[ii] = 0.0;

  /* De-disperse */

  for (ii = 0; ii < numchan; ii++) {

    /* Organize the input data from *lastdata */
    
    ptr = ii + offset[ii] * numchan;
    for (jj = 0; jj < numpts - offset[ii]; jj++, ptr += numchan)
      tmpdat[jj] = lastdata[ptr];

    /* Organize the input data from *data */
    
    tmpptr = numpts - offset[ii];
    ptr = ii;
    for (jj = 0; jj <= offset[ii]; jj++, ptr += numchan)
      tmpdat[jj + tmpptr] = data[ptr];

    /* Now combine the input data with the proper weights  */
    
    for (jj = 0; jj < numpts; jj++)
      result[jj] += hifrac[ii] * tmpdat[jj] + lofrac[ii] * tmpdat[jj + 1];
  }
}








