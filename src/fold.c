#include "presto.h"

/* The number of points to work with at a time from the input file */
#define WORKLEN  8192

void foldfile(FILE *datafile, double dt, double *prof, long proflen, \
	      double fo, double fdot, double fdotdot, int binary, \
	      double *delays, double orbto, double orbdt, long numdelays, \
	      int poisson, double *avg, double *var, float *chiarr, \
	      double *onoffpairs, long *totnumfolded)
/*
 * This routine is a general pulsar folding algorithm.  It will fold
 * pulsar data for a pulsar with frequency derivatives and double
 * derivatives as well as pulsars in binary orbits.  Additional
 * applications can be added by changing the makeup of the arrays
 * *delays and *delaytimes.  The array *delays descirbes the sum
 * of all time-of-arrival delays that act on the pulses.  For a binary
 * system, the dominant component is the light-propagation delay caused
 * by the pulsar orbit.  The array *delaytimes simply describes the
 * times where the delays were sampled.  The delays are linearly
 * interpolated using both *delays and *delaytimes.  The function
 * arguments are as follows:
 *
 * FILE *datafile:         The input data file pointer containing the
 *                         floating point data.
 * double dt:              The integration length of each bin in
 *                         the input file.
 * double *prof:           A double precision array of length proflen
 *                         to contain the profile.
 * long proflen:           The length of the profile array.
 * double fo:              The starting frequency (hz) to fold.
 * double fdot:            The starting frequency derivative (hz/s).
 * double fdotdot:         The frequency double-derivative (hz/s^2).
 * int binary:             1 = Use the delays, 0 = do not use the delays.
 * double *delays:         An array of pulse TOA delays (s).
 * double orbto:           The starting time of the evenly spaced *delays
 * double orbdt:           The time interval used in sampling the orbit (s).
 * long numdelays:         The number of delays in *delays.
 * int poisson:            True if we think the data is Poisson distributed.
 * double *avg:            (Return val) The average value per time series bin.
 * double *var:            (Return val) The variance per time series bin.
 * float *chiarr:          An array containing the instant chi-square
 *                         during the folding (1 point each WORKLEN data).
 * double *onoffpairs:     An array containing pairs of numbers that represent
 *                         the bin numbers when to actively add to the profile.
 *                         for the whole file, the array should be [0,N-1].
 * long *totnumfolded:     (Return val) The total number of bins folded.
 */
{
  int loprofbin = 0, hiprofbin, numbins, modhiprofbin, stopit = 0;
  float data[WORKLEN];
  unsigned long N;
  long i, j, numread, numreads, pt = 1;
  double phasenext, loprofphase = 0.0, hiprofphase, deltaphase, T;
  double phase = 0.0, Tnext, tbnext = 0.0, tbkeep, integer;
  double *profbinphases, profbinwidth, dtmp, dev = 0.0;
  double lopart = 0.0, midpart = 0.0, hipart = 0.0;
  double avgph = 0.0, varph = 0.0, chitmp = 0.0;
  double *onptr, *offptr, tbon, tboff, tbtmp, orbmaxt;

  /* Get the data file length */

  N = chkfilelen(datafile, sizeof(float));
  T = N * dt;
  numreads = N / WORKLEN;
  orbmaxt = orbto + numdelays * orbdt;
  *totnumfolded = 0;
  *avg = 0.0;
  *var = 0.0;
  
  /* Initiate the on-off pointers and variables */

  onptr = onoffpairs;
  offptr = onoffpairs + 1;
  tbon = *onptr * dt;
  tboff = *offptr * dt;
  tbkeep = tbon;

  /* Save some floating point ops later... */

  fdot /= 2.0;
  fdotdot /= 6.0;
  profbinwidth = 1.0 / proflen;

  /* Set up the initial phase values and the file */

  tbtmp = tbon;
  if (binary){
    tbtmp -= lin_interp_E(delays, tbon, orbto, orbdt, orbmaxt);
  }
  phase = tbtmp * (tbtmp * (tbtmp * fdotdot + fdot) + fo);
  loprofphase = (phase >= 0.0) ? \
    modf(phase, &integer) : 1.0 + modf(phase, &integer);

  /* Move to the first data point to fold */

  chkfileseek(datafile, (long) (*onptr),  sizeof(float), SEEK_SET);

  /* Generate an array that contains the starting cyclic (0.0-1.0) */
  /* pulse phase that each profile array bin represents.           */

  profbinphases = gen_dvect(proflen+1);
  for (i = 0; i <= proflen; i++)
    profbinphases[i] = (double) (i) / proflen;

  /* Generate the profile */

  while ((numread = chkfread(data, sizeof(float),			
			     WORKLEN, datafile)) == WORKLEN) {		

    for (i = 0; i < WORKLEN; i++) {

      /* Check if we are within the current onoff boundary.      */
      /* Advance the pointers and other variables if we are not. */

      if (tbnext >= tboff){
	onptr += 2;
	offptr += 2;
	if (*onptr == 0.0){
	  stopit = 1;
	  break;
	}
	tbon = *onptr * dt;
	tboff = *offptr * dt;
	tbnext = tbon;
	chkfileseek(datafile, (long) (*onptr), sizeof(float), SEEK_SET);
	tbtmp = tbon;
	if (binary){
	  tbtmp -= lin_interp_E(delays, tbon, orbto, orbdt, orbmaxt);
	}
	phase = tbtmp * (tbtmp * (tbtmp * fdotdot + fdot) + fo);
	loprofphase = (phase >= 0.0) ? \
	  modf(phase, &integer) : 1.0 + modf(phase, &integer);
	loprofbin = (int) floor(proflen * loprofphase + DBLCORRECT);
	break;
      }

      /* Calculate the barycentric time for the next point. */
      
      tbnext = tbkeep + (i + 1) * dt;
      
      /* Correct the next barycentric time to pulsar proper time. */
      
      if (binary)
	Tnext = tbnext - \
	  lin_interp_E(delays, tbnext, orbto, orbdt, orbmaxt);
      else 
	Tnext = tbnext;

      /* Get the pulsar phase (cyclic) for the next point. */
      
      phasenext = Tnext * (Tnext * (Tnext * fdotdot + fdot) + fo);

      /* How much total phase does the data point cover? */

      deltaphase = phasenext - phase;

      /* Find the highest numbered bin we will add data to.   */
      /* Note:  This number will be used modulo proflen so it */
      /*        could be greater than proflen.                */
      
      hiprofphase = loprofphase + deltaphase;
      hiprofbin = (int) floor(proflen * hiprofphase + DBLCORRECT);
      modhiprofbin = hiprofbin % proflen;
      
      /* How many profile bins we will spread the data over? */

      numbins = hiprofbin - loprofbin + 1;

      /* Spread the data into the proper bins. */

      if (numbins >= 3){

	/* Data point will be spread over 3 or more profile bins */

	hipart = data[i] * \
	  (profbinphases[loprofbin + 1] - loprofphase) / deltaphase;
	dtmp = modf(hiprofphase, &integer);
	if (dtmp == 0.0) dtmp = 1.0;
	lopart = data[i] * \
	  (dtmp - profbinphases[modhiprofbin]) / deltaphase;
	midpart = data[i] * profbinwidth / deltaphase;
	prof[loprofbin] += hipart;
	prof[modhiprofbin] += lopart;
	for (j = loprofbin + 1; j < hiprofbin; j++){
	  prof[j % proflen] += midpart;
	}

      } else if (numbins == 2) {

	/* Data point will be spread over 2 profile bins */

	dtmp = profbinphases[modhiprofbin];
	if (dtmp == 0.0) dtmp = 1.0;
	hipart = data[i] * (dtmp - loprofphase) / deltaphase;
	lopart = data[i] - hipart;
	prof[loprofbin] += hipart;
	prof[modhiprofbin] += lopart;

      } else {

	/* Data point will go into only 1 profile bin */

	prof[loprofbin] += data[i];

      }
	
      /* Update variables */
      
      loprofphase = modf(hiprofphase, &integer);
      loprofbin = (int) floor(proflen * loprofphase + DBLCORRECT);
      phase = phasenext;

      /* Use clever single pass mean and variance calculation */
      
      dev = data[i] - *avg;
      *avg += dev / (*totnumfolded + i + 1);
      *var += dev * (data[i] - *avg);
      
    }
    
    *totnumfolded += i;
    tbkeep = tbnext;

    /* Average value of a profile bin... */
    
    avgph = *avg * *totnumfolded * profbinwidth;

    /* Variance of a profile bin... */

    varph = *var * profbinwidth;

    /* Compute the Chi-Squared probability that there is a signal */
    /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */

    pt = (long) floor((tbkeep / T) * (numreads+1));

    if (poisson){
      for (i = 0 ; i < proflen ; i++){
	dtmp = prof[i];
	chitmp = dtmp - avgph;
	dtmp = (dtmp == 0.0) ? 1.0 : dtmp;
	chiarr[pt] += ((chitmp * chitmp) / dtmp);
      }
    } else {
      for (i = 0 ; i < proflen ; i++){
	chitmp = prof[i] - avgph;
	chiarr[pt] += chitmp * chitmp;
      }
      chiarr[pt] /= varph;
    }      
    chiarr[pt] /= (proflen - 1.0);

    /* Break out of the loop if we are done due to onoff */

    if (stopit) break;

  }
  
  *var /= (*totnumfolded - 1.0);

  /* Cleanup */

  free(profbinphases);

}



void fold(float *data, long N, double dt, double tb, double *prof, \
	  long proflen, double fo, double fdot, double fdotdot, \
	  int binary, double *delays, double orbto, double orbdt, \
	  long numdelays, double *onoffpairs, long *totnumfolded)
/*
 * This routine is a general pulsar folding algorithm.  It will fold
 * pulsar data for a pulsar with frequency derivatives and double
 * derivatives as well as pulsars in binary orbits.  Additional
 * applications can be added by changing the makeup of the arrays
 * *delays and *delaytimes.  The array *delays descirbes the sum
 * of all time-of-arrival delays that act on the pulses.  For a binary
 * system, the dominant component is the light-propagation delay caused
 * by the pulsar orbit.  The array *delaytimes simply describes the
 * times where the delays were sampled.  The delays are linearly
 * interpolated using both *delays and *delaytimes.  The profile
 * will always have tb = 0.0 at profile or fold phase 0.0.  The function
 * arguments are as follows:
 *
 * float *data:            A float pointer containing the data to fold.
 * long N:                 The number of points in *data.  (or number to fold)
 * double dt:              The integration length of each data bin.
 * double tb:              The barycentric time at the start of data (s).
 *                         This time must correspond to correct values 
 *                         of fo, fdot, fdotdot, and *delaytimes.
 * double *prof:           A double precision array of length proflen
 *                         to contain the profile.
 * long proflen:           The length of the profile array.
 * double fo:              The starting frequency (hz) to fold.
 * double fdot:            The starting frequency derivative (hz/s).
 * double fdotdot:         The frequency double-derivative (hz/s^2).
 * int binary:             1 = Use the delays, 0 = do not use the delays.
 * double *delays:         An array of pulse TOA delays (s).
 * double orbto:           The starting time of the evenly spaced *delays
 * double orbdt:           The time interval used in sampling the orbit (s).
 * long numdelays:         The number of delays in *delays.
 * double *onoffpairs:     An array containing pairs of numbers that represent
 *                         the bins when we will actively add to the profile.
 *                         For the whole array, onoffpairs should be [0,N-1].
 * long *totnumfolded:     (Return val) The total number of bins folded.
 */
{
  int loprofbin = 0, hiprofbin, numbins, modhiprofbin, stopit = 0;
  long i, j;
  double phase = 0.0, phasenext = 0.0, deltaphase = 0.0, Tnext = 0.0, tbnext;
  double *profbinphases, profbinwidth, loprofphase = 0.0, hiprofphase = 0.0;
  double lopart = 0.0, midpart = 0.0, hipart = 0.0, integer, dtmp;
  double *onptr, *offptr, tbtmp = 0.0, orbmaxt;

  /* Save some floating point ops later... */

  orbmaxt = orbto + numdelays * orbdt;
  fdot /= 2.0;
  fdotdot /= 6.0;
  profbinwidth = 1.0 / proflen;

  /* Initiate the on-off pointers and variables */

  onptr = onoffpairs;
  offptr = onoffpairs + 1;
  tbnext = tb;

  /* Correct the starting barycentric time to pulsar proper time. */
  
  if (binary)
    tbtmp = tb - lin_interp_E(delays, tb, orbto, orbdt, orbmaxt);
  else 
    tbtmp = tb;

  /* Get the starting pulsar phase (cyclic). */
  
  phase = tbtmp * (tbtmp * (tbtmp * fdotdot + fdot) + fo);
  loprofphase = (phase >= 0.0) ? \
    modf(phase, &integer) : 1.0 + modf(phase, &integer);
  loprofbin = (int) floor(proflen * loprofphase + DBLCORRECT);
  
  /* Generate an array that contains the starting cyclic (0.0-1.0) */
  /* pulse phase that each profile array bin represents.           */

  profbinphases = gen_dvect(proflen+1);
  for (i = 0; i <= proflen; i++)
    profbinphases[i] = (double) (i) / proflen;

  /* Generate the profile */

  for (i = 0; i < N; i++) {
    
    /* Check if we are within the current onoff boundary.      */
    /* Advance the pointers and other variables if we are not. */
    
    if (i > *offptr){
      onptr += 2;
      offptr += 2;
      if (*onptr == 0.0){
	stopit = 1;
	break;
      }
      tbnext = *onptr * dt + tb;
      i = *onptr;
      tbtmp = tbnext;
      if (binary){
	tbtmp -= lin_interp_E(delays, tbnext, orbto, orbdt, orbmaxt);
      }
      phase = tbtmp * (tbtmp * (tbtmp * fdotdot + fdot) + fo);
      loprofphase = (phase >= 0.0) ? \
	modf(phase, &integer) : 1.0 + modf(phase, &integer);
      loprofbin = (int) floor(proflen * loprofphase + DBLCORRECT);
      break;
    }
    
    /* Calculate the barycentric time for the next point. */
    
    tbnext = tb + (i + 1) * dt;
    
    /* Correct the next barycentric time to pulsar proper time. */
    
    if (binary)
      Tnext = tbnext - lin_interp_E(delays, tbnext, orbto, orbdt, orbmaxt);
    else 
      Tnext = tbnext;
    
    /* Get the pulsar phase (cyclic) for the next point. */
    
    phasenext = Tnext * (Tnext * (Tnext * fdotdot + fdot) + fo);
    
    /* How much total phase does the data point cover? */
    
    deltaphase = phasenext - phase;

    /* Find the highest numbered bin we will add data to.   */
    /* Note:  This number will be used modulo proflen so it */
    /*        could be greater than proflen.                */
    
    hiprofphase = loprofphase + deltaphase;
    hiprofbin = (int) floor(proflen * hiprofphase + DBLCORRECT);
    modhiprofbin = hiprofbin % proflen;
    
    /* How many profile bins we will spread the data over? */
    
    numbins = hiprofbin - loprofbin + 1;
    
    /* Spread the data into the proper bins. */
    
    if (numbins >= 3){
      
      /* Data point will be spread over 3 or more profile bins */
      
      hipart = data[i] * \
	(profbinphases[loprofbin + 1] - loprofphase) / deltaphase;
      dtmp = modf(hiprofphase, &integer);
      if (dtmp == 0.0) dtmp = 1.0;
      lopart = data[i] * \
	(dtmp - profbinphases[modhiprofbin]) / deltaphase;
      midpart = data[i] * profbinwidth / deltaphase;
      prof[loprofbin] += hipart;
      prof[modhiprofbin] += lopart;
      for (j = loprofbin + 1; j < hiprofbin; j++){
	prof[j % proflen] += midpart;
      }
      
    } else if (numbins == 2) {
      
      /* Data point will be spread over 2 profile bins */
      
      dtmp = profbinphases[modhiprofbin];
      if (dtmp == 0.0) dtmp = 1.0;
      hipart = data[i] * (dtmp - loprofphase) / deltaphase;
      lopart = data[i] - hipart;
      prof[loprofbin] += hipart;
      prof[modhiprofbin] += lopart;
      
    } else {
      
      /* Data point will go into only 1 profile bin */
      
      prof[loprofbin] += data[i];
      
    }
    
    /* Update variables */
    
    loprofphase = modf(hiprofphase, &integer);
    loprofbin = (int) floor(proflen * loprofphase + DBLCORRECT);
    phase = phasenext;
    (*totnumfolded)++;
    
  }
  
  /* Cleanup */

  free(profbinphases);
}



