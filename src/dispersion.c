#include "presto.h"

double tree_max_dm(int numchan, double dt, double lofreq, double hifreq)
/* Return the maximum Dispersion Measure (dm) in cm-3 pc, the  */
/* tree de-dispersion technique can correct for given a sample */
/* interval 'dt', the number of channels 'numchan', and the    */
/* low and high observation frequencies in MHz.                */
{
  if (lofreq==0.0 || hifreq==0.0)
    return 0.0;
  else
    return 0.000241 * (numchan - 1) * dt / 
      ((1.0 / (lofreq * lofreq)) - (1.0 / (hifreq * hifreq)));
}


double smearing_from_bw(double dm, double center_freq, double bandwidth)
/* Return the smearing in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, the central frequency */
/* and the bandwith of the observation in MHz.                 */
{
  if (center_freq==0.0)
    return 0.0;
  else
    return dm * bandwidth / 
      (0.0001205 * center_freq * center_freq * center_freq);
}


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


double *dedisp_delays(int numchan, double dm, double lofreq, 
		      double chanwidth)
/* Return an array of delays (sec) for dedispersing 'numchan'  */
/* channels at a DM of 'dm'.  'lofreq' is the center frequency */
/* in MHz of the lowest frequency channel.  'chanwidth' is the */
/* width in MHz of each channel.  The returned array is        */
/* allocated by this routine.                                  */
{
  int ii;
  double *delays;

  delays = gen_dvect(numchan);
  for (ii = 0; ii < numchan; ii++)
    delays[ii] = delay_from_dm(dm, lofreq + ii * chanwidth);
  return delays;
}


void dedisp(unsigned char *data, unsigned char *lastdata, int numpts,
	    int numchan, double *dispdelays, float *result)
/* De-disperse a stretch of data with numpts * numchan points. */
/* The delays (in bins) are in dispdelays for each channel.    */
/* The result is returned in result.  The input data and       */
/* dispdelays are always in ascending frequency order.         */
/* Input data are ordered in time, with the channels stored    */
/* together at each time point.                                */ 
{
  static int approx_mean, firsttime = 1, *offset;
  int ii, jj, kk;

  if (firsttime){
    offset = gen_ivect(numchan);
    for (ii = 0; ii < numchan; ii++){
      if (dispdelays[ii] < 0.0){
	printf("\ndispdelays[%d] = %f is < 0.0 in dedisp().\n\n",
	       ii, dispdelays[ii]);
	exit(-1);
      }
      offset[ii] = (int) dispdelays[ii];
    }
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


double *subband_search_delays(int numchan, int numsubbands, double dm, 
			      double lofreq, double chanwidth)
/* Return an array of delays (sec) for a subband DM search.  The      */
/* delays are calculated normally for each of the 'numchan' channels  */
/* using the appropriate frequencies at the 'dm'.  Then the delay     */
/* from the highest frequency channel of each of the 'numsubbands'    */
/* subbands is subtracted from each subband.  This gives the subbands */
/* the correct delays for each freq in the subband, but the subbands  */
/* themselves are offset as if they had not been de-dispersed.  This  */
/* way, we can call dedisp() on the group of subbands if needed.      */
/* 'lofreq' is the center frequency in MHz of the lowest frequency    */
/* channel.  'chanwidth' is the width in MHz of each channel.  The    */
/* returned array is allocated by this routine.                       */
/* Note:  When performing a subband search, the delays for each       */
/*   subband must be calculated with the frequency of the highest     */
/*   channel in each subband, _not_ the center subband frequency.     */
{
  int ii, jj, chan_per_subband;
  double *delays, *subbanddelays, subbandwidth, losub_hifreq;

  chan_per_subband = numchan / numsubbands;
  subbandwidth = chanwidth * chan_per_subband;
  losub_hifreq = lofreq + subbandwidth - chanwidth;

  /* Calculate the appropriate delays to subtract from each subband */

  subbanddelays = dedisp_delays(numsubbands, dm, losub_hifreq, 
				subbandwidth);

  /* Calculate the appropriate delays for each channel */

  delays = dedisp_delays(numchan, dm, lofreq, chanwidth);
  for (ii = 0; ii < numsubbands; ii++)
    for (jj = 0; jj < chan_per_subband; jj++)
      delays[ii * chan_per_subband + jj] -= subbanddelays[ii];
  free(subbanddelays);

  return delays;
}


void dedisp_subbands(unsigned char *data, unsigned char *lastdata,
		     int numpts, int numchan, double *dispdelays,
		     int numsubbands, float *result)
/* De-disperse a stretch of data with numpts * numchan points    */
/* into numsubbands subbands.  Each time point for each subband  */
/* is a float in the result array.  The result array order is    */
/* subbands of increasing frequency together at each time pt.    */
/* The delays (in bins) are in dispdelays for each channel.      */
/* The input data and dispdelays are always in ascending         */
/* frequency order.  Input data are ordered in time, with the    */
/* channels stored together at each time point.                  */ 
{
  static int approx_mean, firsttime = 1, *offset, chan_per_subband;
  int ii, jj, kk, ll, chan;

  if (firsttime){
    offset = gen_ivect(numchan);
    for (ii = 0; ii < numchan; ii++){
      if (dispdelays[ii] < 0.0){
	printf("\ndispdelays[%d] = %f is < 0.0 in dedisp_subbands().\n\n",
	       ii, dispdelays[ii]);
	exit(-1);
      }
      offset[ii] = (int) dispdelays[ii];
    }
    chan_per_subband = numchan / numsubbands;
    approx_mean = -(numchan / 2 - 1);
/*     approx_mean = -(chan_per_subband / 2); */
    firsttime = 0;
  }

  /* Set the result array to negative of numchan / 2. */
  /* This will result in data with approx zero mean.  */

  for (ii = 0; ii < numpts * numsubbands; ii++)
    result[ii] = approx_mean;

  /* De-disperse into the subbands */

  for (ii = 0; ii < numsubbands; ii++){
    chan = ii * chan_per_subband;
    for (jj = 0; jj < chan_per_subband; jj++, chan++){
      kk = chan + offset[chan] * numchan;
      for (ll = 0; ll < (numpts - offset[chan]) * numsubbands; 
	   ll += numsubbands, kk += numchan)
	result[ll + ii] += lastdata[kk];
      kk = chan;
      for (; ll < (numpts * numsubbands); 
	   ll += numsubbands, kk += numchan)
	result[ll + ii] += data[kk];
    }
  }
}







