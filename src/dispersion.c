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
	    double *dispdelays, long numchan, float *result, int order)
/* De-disperse a stretch of data with numpts * numchan points. */
/* The delays (in bins) are in dispdelays for each channel.    */
/* The result is returned in result.  If order is positive,    */
/* the channels are in order from lowest freq to highest.  If  */
/* order is negative, the reverse is true.  Note that the      */
/* dispdelays are always in ascending freq order though.       */
/* Data are ordered in time, with each channel in a row for    */
/* each time point.                                            */
{
  static double *lofrac, *hifrac;
  static float *tmpdat;
  static int firsttime = 1, *offset, chan = 0, dchan = 0;
  long i, j, ptr, tmpptr;

  if (firsttime) {
    lofrac = gen_dvect(numchan);
    hifrac = gen_dvect(numchan);
    offset = gen_ivect(numchan);
    tmpdat = gen_fvect(numpts + 1);

    for (i = 0; i < numchan; i++) {
      offset[i] = (int) floor(dispdelays[i]);
      lofrac[i] = dispdelays[i] - offset[i];
      hifrac[i] = 1.0 - lofrac[i];
    }

    if (order > 0) {
      chan = 0;
      dchan = 1;
    } else {
      chan = numchan - 1;
      dchan = -1;
    }
    
    firsttime = 0;
  }
  /* Reset the result array to 0's */

  for (i = 0; i < numpts; i++) {
    result[i] = 0.0;
  }

  /* De-disperse */

  for (i = 0; i < numchan; i++) {

    /* Organize the input data from *lastdata */
    
    ptr = chan + i * dchan + offset[i] * numchan;
    for (j = 0; j < numpts - offset[i]; j++, ptr += numchan)
      tmpdat[j] = lastdata[ptr];

    /* Organize the input data from *data */
    
    tmpptr = numpts - offset[i];
    ptr = chan + i * dchan;
    for (j = 0; j <= offset[i]; j++, ptr += numchan)
      tmpdat[j + tmpptr] = data[ptr];

    /* Now combine the input data with the proper weights  */
    
    for (j = 0; j < numpts; j++)
      result[j] += hifrac[i] * tmpdat[j] + lofrac[i] * tmpdat[j + 1];
  }
}








