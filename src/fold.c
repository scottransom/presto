#include "presto.h"

/* The number of points to work with at a time         */
/* This must be the same as the WORKLEN in profile.c!  */
#define WORKLEN   16384

/* Some macros to make the flag checking easier */
#define DELAYS (flags % 2 == 1)
#define ONOFF (flags > 1)

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);

double foldfile(FILE *datafile, double dt, double tlo, 
		double *prof, int numprof, double startphs, 
		double fo, double fdot, double fdotdot, int flags, 
		double *delays, double *delaytimes, int numdelays, 
		double *onoffpairs, foldstats *stats, float *chiarr)
/* This routine is a general pulsar folding algorithm.  It will fold  */
/* data for a pulsar with single and double frequency derivatives and */
/* with arbitrary pulse delays (for example: variable time delays     */
/* due to light travel time in a binary).  These delays are described */
/* in the arrays '*delays' and '*delaytimes'. The folding may also be */
/* turned on and off throughout the data by using 'onoffpairs'. The   */
/* profile will have the data corresponding to time 'tlo' placed at   */
/* the phase corresponding to time 'tlo' using 'fo', 'fdot', and      */
/* 'fdotdot' plus 'startphs' and the appropriate delay.               */
/* Arguments:                                                         */
/*    'datafile' is FILE ptr for the input floating point data file.  */
/*    'dt' is the time duration of each data bin.                     */
/*    'tlo' is the time of the start of the 1st data pt.              */
/*    'prof' is a double prec array to contain the profile.           */
/*    'numprof' is the length of the profile array.                   */
/*    'startphs'is the phase offset [0-1] for the first point.        */
/*    'fo' the starting frequency to fold.                            */
/*    'fdot' the starting frequency derivative.                       */
/*    'fdotdot' the frequency second derivative.                      */
/*    'flags' is an integer containing flags of how to fold:          */
/*            0 = No *delays and no *onoffpairs                       */
/*            1 = Use *delays but no *onoffpairs                      */
/*            2 = No *delays but use *onoffpairs                      */
/*            3 = Use *delays and use *onoffpairs                     */
/*    'delays' is an array of time delays.                            */
/*    'delaytimes' are the times where 'delays' were calculated.      */
/*    'numdelays' is how many points are in 'delays' and 'delaytimes' */
/*    'onoffpairs' is array containing pairs of normalized times      */
/*            that represent the bins when we will actively add       */
/*            to the profile.  To fold the full file,                 */
/*            onoffpairs should be [0.0, T_obs].                      */
/*    'stats' are statistics of the data and the profile.             */
/*    'chiarr' is an array containing the instant reduced chi-square  */
/*            during the folding (1 point each WORKLEN data).  This   */
/*            array must have been allocated and set to 0.            */
/* Notes:  fo, fdot, and fdotdot correspon to 'tlo' = 0.0             */
/*    (i.e. to the beginning of the first data point)                 */
{
  float data[WORKLEN];
  double *onoffptr=NULL, phase=0.0;
  int ourflags;
  unsigned long ii, N, onbin, offbin, numbins;
  unsigned long remainbins, binstoread, numreads;

  /* Get the data file length and initialize some variables */

  N = chkfilelen(datafile, sizeof(float));
  if (ONOFF) onoffptr = onoffpairs;
  stats->numdata = stats->data_avg = stats->data_var = 0.0;
  if (DELAYS) ourflags = 1;
  else ourflags = 0;

  do { /* Loop over the on-off pairs */

    /* Set the on-off variables */
    
    if (ONOFF){
      onbin = (unsigned long) (*onoffptr * N + DBLCORRECT);
      offbin = (unsigned long) (*(onoffptr + 1) * N  + DBLCORRECT);
      if (offbin) offbin--;
      onoffptr += 2;
    } else {
      onbin = 0;
      offbin = N - 1;
    }
    numbins = (offbin == onbin) ? 0 : offbin - onbin + 1;
    numreads = numbins / WORKLEN;
    remainbins = numbins % WORKLEN;
    if (remainbins) numreads++;
    binstoread = WORKLEN;

    /* Skip to the correct file location */

    chkfileseek(datafile, onbin, sizeof(float), SEEK_SET);

    /* Loop over the number of reads we have to perform for */
    /* the current on-off pair.                             */

    for (ii = 0; ii < numreads; ii++, onbin += binstoread){
      
      /* Correct for the fact that our last read might be short */
      
      if (remainbins && (ii == numreads - 1))
	binstoread = remainbins;

      /* Read the current chunk of data */

      chkfread(data, sizeof(float), binstoread, datafile);

      /* Fold the current chunk of data */

      phase = fold(data, binstoread, dt, tlo + onbin * dt, prof, 
		   numprof, startphs, fo, fdot, fdotdot, ourflags, 
		   delays, delaytimes, numdelays, NULL, stats);

      /* Set the current chiarr value */

      chiarr[onbin / WORKLEN] = (float) stats->redchi;
    }

  } while (offbin < N - 1 && offbin != 0);

  /* Return the ending phase from folding */

  return phase;
}


double simplefold(float *data, int numdata, double dt, double tlo,
		  double *prof, int numprof, double startphs, 
		  double fo, double fdot, double fdotdot)
/* This routine is a simplified pulsar folding algorithm.  It    */
/* folds data for a pulsar with single and double frequency      */
/* derivatives.  The profile will have the data corresponding    */
/* to time 'tlo' placed at the phase corresponding to time 'tlo' */
/* using 'fo', 'fdot', and 'fdotdot' plus 'startphs'.            */
/* Arguments:                                                    */
/*    'data' is a float array containing the data to fold.       */
/*    'numdata' is the number of points in *data.                */
/*    'dt' is the time duration of each data bin.                */
/*    'tlo' is the time of the start of the 1st data pt.         */
/*    'prof' is a double prec array to contain the profile.      */
/*    'numprof' is the length of the profile array.              */
/*    'startphs'is the phase offset [0-1] for the first point.   */
/*    'fo' the starting frequency to fold.                       */
/*    'fdot' the starting frequency derivative.                  */
/*    'fdotdot' the frequency second derivative.                 */
/* Notes:  fo, fdot, and fdotdot correspon to 'tlo' = 0.0        */
/*    (i.e. to the beginning of the first data point)            */
{
  int ii, jj, loprofbin, hiprofbin, numbins, modhiprofbin;
  double phase, phasenext = 0.0, deltaphase, rdeltaphase, Tnext;
  double profbinwidth, loprofphase, hiprofphase, tmpphase;
  double lopart, midpart, hipart, dtmp;

  /* Save some floating point ops later... */

  fdot /= 2.0;
  fdotdot /= 6.0;
  profbinwidth = 1.0 / numprof;

  /* Get the starting pulsar phase (cyclic). */
  
  phase = tlo * (tlo * (tlo * fdotdot + fdot) + fo) + startphs;
  loprofphase = (phase < 0.0) ? 
    1.0 + phase - (int) phase : phase - (int) phase;
  loprofbin = (int) (loprofphase * numprof + DBLCORRECT);
  
  /* Generate the profile */

  for (ii = 0; ii < numdata; ii++) {
    
    /* Get the pulsar phase (cyclic) for the next point. */
    
    Tnext = tlo + (ii + 1) * dt;
    phasenext = Tnext * (Tnext * (Tnext * fdotdot + fdot) 
			 + fo) + startphs;
    
    /* How much total phase does the data point cover? */
    
    deltaphase = phasenext - phase;
    rdeltaphase = 1.0 / deltaphase;

    /* Find the highest numbered bin we will add data to.   */
    /* Note:  This number will be used modulo numprof so it */
    /*        could be greater than numprof.                */
    
    hiprofphase = loprofphase + deltaphase;
    hiprofbin = (int) (hiprofphase * numprof + DBLCORRECT);
    modhiprofbin = hiprofbin % numprof;
    
    /* How many profile bins we will spread the data over? */
    
    numbins = hiprofbin - loprofbin + 1;
    
    /* Spread the data into the proper bins. */
    
    if (numbins >= 3){
      
      /* Data point will be spread over 3 or more profile bins */
      
      dtmp = data[ii] * rdeltaphase;
      hipart = dtmp * ((loprofbin + 1) * profbinwidth - loprofphase);
      tmpphase = hiprofphase - (int) hiprofphase;
      tmpphase = (tmpphase == 0.0) ? 1.0 : tmpphase;
      lopart = dtmp * (tmpphase - modhiprofbin * profbinwidth);
      midpart = dtmp * profbinwidth;
      prof[loprofbin] += hipart;
      prof[modhiprofbin] += lopart;
      for (jj = loprofbin + 1; jj < hiprofbin; jj++)
	prof[jj % numprof] += midpart;
      
    } else if (numbins == 2) {
      
      /* Data point will be spread over 2 profile bins */
      
      tmpphase = modhiprofbin * profbinwidth;
      tmpphase = (tmpphase == 0.0) ? 1.0 : tmpphase;
      hipart = data[ii] * (tmpphase - loprofphase) * rdeltaphase;
      lopart = data[ii] - hipart;
      prof[loprofbin] += hipart;
      prof[modhiprofbin] += lopart;
      
    } else {
      
      /* Data point will go into only 1 profile bin */
      
      prof[loprofbin] += data[ii];
    }
    
    /* Update variables */
    
    loprofphase = hiprofphase - (int) hiprofphase;
    loprofbin = (int) (loprofphase * numprof + DBLCORRECT);
    phase = phasenext;

  }

  phasenext = (phasenext < 0.0) ? 
    1.0 + phasenext - (int) phasenext : phasenext - (int) phasenext;
  return phasenext;
}


double fold(float *data, int numdata, double dt, double tlo, 
	    double *prof, int numprof, double startphs, 
	    double fo, double fdot, double fdotdot, int flags, 
	    double *delays, double *delaytimes, int numdelays, 
	    int *onoffpairs, foldstats *stats)
/* This routine is a general pulsar folding algorithm.  It will fold  */
/* data for a pulsar with single and double frequency derivatives and */
/* with arbitrary pulse delays (for example: variable time delays     */
/* due to light travel time in a binary).  These delays are described */
/* in the arrays '*delays' and '*delaytimes'. The folding may also be */
/* turned on and off throughout the data by using 'onoffpairs'. The   */
/* profile will have the data corresponding to time 'tlo' placed at   */
/* the phase corresponding to time 'tlo' using 'fo', 'fdot', and      */
/* 'fdotdot' plus 'startphs' and the appropriate delay.               */
/* Arguments:                                                         */
/*    'data' is a float array containing the data to fold.            */
/*    'numdata' is the number of points in *data.                     */
/*    'dt' is the time duration of each data bin.                     */
/*    'tlo' is the time of the start of the 1st data pt.              */
/*    'prof' is a double prec array to contain the profile.           */
/*    'numprof' is the length of the profile array.                   */
/*    'startphs'is the phase offset [0-1] for the first point.        */
/*    'fo' the starting frequency to fold.                            */
/*    'fdot' the starting frequency derivative.                       */
/*    'fdotdot' the frequency second derivative.                      */
/*    'flags' is an integer containing flags of how to fold:          */
/*            0 = No *delays and no *onoffpairs                       */
/*            1 = Use *delays but no *onoffpairs                      */
/*            2 = No *delays but use *onoffpairs                      */
/*            3 = Use *delays and use *onoffpairs                     */
/*    'delays' is an array of time delays.                            */
/*    'delaytimes' are the times where 'delays' were calculated.      */
/*    'numdelays' is how many points are in 'delays' and 'delaytimes' */
/*    'onoffpairs' is array containing pairs of numbers that          */
/*            represent the bins when we will actively add            */
/*            to the profile.  To fold the whole array,               */
/*            onoffpairs should be [0, numdata-1].                    */
/*    'stats' are statistics of the data that were folded as well     */
/*            as the folded profile itself.  If this                  */
/*            routine is used on consecutive pieces of the            */
/*            same data, fold() will use the current values           */
/*            and update them at the end of each call.                */
/*            So each parameter must be set to 0.0 before             */
/*            fold() is called for the first time.                    */
/* Notes:  fo, fdot, and fdotdot correspon to 'tlo' = 0.0             */
/*    (i.e. to the beginning of the first data point)                 */
{
  int loprofbin, hiprofbin, numbins, modhiprofbin;
  int ii, jj, onbin, offbin, *onoffptr=NULL;
  unsigned long arrayoffset=0;
  double phase, phasenext=0.0, deltaphase, T, Tnext, TD, TDnext;
  double profbinwidth, loprofphase, hiprofphase;
  double lopart, midpart, hipart, tmpphase, dtmp, rdeltaphase;
  double dev, delaytlo=0.0, delaythi=0.0, delaylo=0.0, delayhi=0.0;
  double *delayptr=NULL, *delaytimeptr=NULL;

  /* Initialize some variables and save some FLOPs later... */

  fdot /= 2.0;
  fdotdot /= 6.0;
  profbinwidth = 1.0 / numprof;
  if (ONOFF) onoffptr = onoffpairs;
  stats->numprof = (double) numprof;
  stats->data_var *= (stats->numdata - 1.0);

  do { /* Loop over the on-off pairs */

    /* Set the on-off pointers and variables */
    
    if (ONOFF){
      onbin = *onoffptr;
      offbin = *(onoffptr + 1);
      onoffptr += 2;
    } else {
      onbin = 0;
      offbin = numdata - 1;
    }
    
    /* Initiate the folding start time */
    
    T = tlo + onbin * dt;
    TD = T;
    
    /* Set the delay pointers and variables */
    
    if (DELAYS){

      /* Guess that the next delay we want is the next available */

      arrayoffset += 2;  /* Beware nasty NR zero-offset kludges! */
      hunt(delaytimes-1, numdelays, T, &arrayoffset);
      arrayoffset--;
      delaytimeptr = delaytimes + arrayoffset;
      delayptr = delays + arrayoffset;
      delaytlo = *delaytimeptr;
      delaythi = *(delaytimeptr + 1);
      delaylo = *delayptr;
      delayhi = *(delayptr + 1);
      
      /* Adjust the folding start time for the delays */
      
      TD -= LININTERP(TD, delaytlo, delaythi, delaylo, delayhi);
    }

    /* Get the starting pulsar phase (cyclic). */
    
    phase = TD * (TD * (TD * fdotdot + fdot) + fo) + startphs;
    loprofphase = (phase < 0.0) ?
      1.0 + phase - (int) phase : phase - (int) phase;
    loprofbin = (int) (loprofphase * numprof + DBLCORRECT);

    /* Generate the profile for this onoff pair */
    
    for (ii = onbin; ii <= offbin; ii++) {

      /* Calculate the barycentric time for the next point. */
      
      Tnext = tlo + (ii + 1) * dt;
      TDnext = Tnext;

      /* Set the delay pointers and variables */
    
      if (DELAYS){
	if (Tnext > delaythi){

	  /* Guess that the next delay we want is the next available */

	  arrayoffset += 2;  /* Beware nasty NR zero-offset kludges! */
	  hunt(delaytimes-1, numdelays, Tnext, &arrayoffset);
	  arrayoffset--;
	  delaytimeptr = delaytimes + arrayoffset;
	  delayptr = delays + arrayoffset;
	  delaytlo = *delaytimeptr;
	  delaythi = *(delaytimeptr + 1);
	  delaylo = *delayptr;
	  delayhi = *(delayptr + 1);
	}

	/* Adjust the folding start time for the delays */

	TDnext -= LININTERP(Tnext, delaytlo, delaythi, delaylo, delayhi);
      }

      /* Get the pulsar phase (cyclic) for the next point. */
      
      phasenext = TDnext * (TDnext * (TDnext * fdotdot + fdot) 
			    + fo) + startphs;
      
      /* How much total phase does the data point cover? */
      
      deltaphase = phasenext - phase;
      rdeltaphase = 1.0 / deltaphase;
      
      /* Find the highest numbered bin we will add data to.   */
      /* Note:  This number will be used modulo numprof so it */
      /*        could be greater than numprof.                */
      
      hiprofphase = loprofphase + deltaphase;
      hiprofbin = (int) (hiprofphase * numprof + DBLCORRECT);
      modhiprofbin = hiprofbin % numprof;
      
      /* How many profile bins we will spread the data over? */
      
      numbins = hiprofbin - loprofbin + 1;
      
      /* Spread the data into the proper bins. */
      
      if (numbins >= 3){
	
	/* Data point will be spread over 3 or more profile bins */
	
	dtmp = data[ii] * rdeltaphase;
	hipart = dtmp * ((loprofbin + 1) * profbinwidth - loprofphase);
	tmpphase = hiprofphase - (int) hiprofphase;
	tmpphase = (tmpphase == 0.0) ? 1.0 : tmpphase;
	lopart = dtmp * (tmpphase - modhiprofbin * profbinwidth);
	midpart = dtmp * profbinwidth;
	prof[loprofbin] += hipart;
	prof[modhiprofbin] += lopart;
	for (jj = loprofbin + 1; jj < hiprofbin; jj++)
	  prof[jj % numprof] += midpart;
	
      } else if (numbins == 2) {
	
	/* Data point will be spread over 2 profile bins */
	
	tmpphase = modhiprofbin * profbinwidth;
	tmpphase = (tmpphase == 0.0) ? 1.0 : tmpphase;
	hipart = data[ii] * (tmpphase - loprofphase) * rdeltaphase;
	lopart = data[ii] - hipart;
	prof[loprofbin] += hipart;
	prof[modhiprofbin] += lopart;
	
      } else {
	
	/* Data point will go into only 1 profile bin */
	
	prof[loprofbin] += data[ii];
      }
    
      /* Update variables */
      
      loprofphase = hiprofphase - (int) hiprofphase;
      loprofbin = (int) (loprofphase * numprof + DBLCORRECT);
      phase = phasenext;

      /* Use clever single pass mean and variance calculation */
      
      stats->numdata += 1.0;
      dev = data[ii] - stats->data_avg;
      stats->data_avg += dev / (stats->numdata + 1.0);
      stats->data_var += dev * (data[ii] - stats->data_avg);
    }

  } while (offbin < numdata - 1 && offbin != 0);
  
  /* Update and correct the statistics */

  stats->prof_avg = 0.0;
  for (ii = 0; ii < numprof; ii++)
    stats->prof_avg += prof[ii];
  stats->prof_avg /= numprof;

  /* Compute the Chi-Squared probability that there is a signal */
  /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */
  
  stats->redchi = 0.0;
  for (ii = 0 ; ii < numprof ; ii++){
    dtmp = prof[ii] - stats->prof_avg;
    stats->redchi += dtmp * dtmp;
  }
  stats->data_var /= (stats->numdata - 1.0);
  stats->prof_var = stats->data_var * dt * stats->numdata * fo;
  stats->redchi /= (stats->prof_var * (numprof - 1));

  phasenext = (phasenext < 0.0) ? 
    1.0 + phasenext - (int) phasenext : phasenext - (int) phasenext;
  return(phasenext);
}

#undef WORKLEN
#undef DELAYS
#undef ONOFF
#undef LININTERP


void shift_prof(double *prof, int proflen, int shift, double *outprof)
/* Rotates a profile 'prof' by an integer 'shift' places.    */
/* If 'shift' < 0 then shift left, 'shift' > 0, shift right. */ 
/* Place the shifted  profile in 'outprof'.                  */
{
  int wrap=0, nowrap=0;

  wrap = shift % proflen;

  /* no shift */

  if (wrap==0){
    memcpy(outprof, prof, proflen * sizeof(double));
    return;

  /* Convert a left shift into the equivalent right shift */

  } else if (wrap < 0)
    wrap += proflen;

  /* Perform a right shift */

  nowrap = proflen - wrap;
  memcpy(outprof, prof + nowrap, wrap * sizeof(double));
  memcpy(outprof + wrap, prof, nowrap * sizeof(double));
}


void combine_profs(double *profs, foldstats *instats, int numprofs, 
		   int proflen, int shift, double *outprof,
		   foldstats *outstats)
/* Combine a series of 'numprofs' profiles, each of length 'proflen', */
/* into a single profile of length 'proflen'.  The profiles are       */
/* summed after being shifted (+:right, -:left) by an an appropriate  */
/* amount such that the phase would drift 'shift' bins over the time  */
/* represented by all of the profiles.  Returns the summed profile in */
/* 'outprof'.  Note that 'profs' must contain all of the profiles     */
/* arranged end-to-end.  Also, 'outprof' must already be allocated.   */
/* The input profile stats in 'instats' are combined and placed in    */
/* the 'outstats' structure.                                          */
{
  int ii, jj, kk, wrap, nowrap, profoffset;

  /* Note:  This routine uses only integer arithmetic in order to   */
  /*        speed up the computations.                              */

  /* Initiate the output statistics */

  outstats->numdata = instats[0].numdata;
  outstats->data_avg = instats[0].data_avg;
  outstats->data_var = instats[0].data_var;
  outstats->numprof = proflen;
  outstats->prof_avg = instats[0].prof_avg;
  outstats->prof_var = instats[0].prof_var;

  /* Set the output array to the first profile */

  memcpy(outprof, profs, sizeof(double) * proflen);

  /* Loop over the profiles */
  
  for (ii = 1; ii < numprofs; ii++){
    
    /* Loop over the low index elements in each profile */
    
    profoffset = ii * proflen;
    wrap = (int)((double) (ii * shift) / 
		 ((double) numprofs) + 0.5) % proflen;
    wrap = (wrap < 0) ? wrap + proflen : wrap;
    nowrap = proflen - wrap;
    for (jj = 0, kk = profoffset + nowrap; jj < wrap; jj++, kk++)
      outprof[jj] += profs[kk];
    
    /* Loop over the high index elements in each profile */
    
    for (kk = profoffset; jj < proflen; jj++, kk++)
      outprof[jj] += profs[kk];
    
    /* Update the output statistics structure */
    
    outstats->numdata += instats[ii].numdata;
    outstats->data_avg += instats[ii].data_avg;
    outstats->data_var += instats[ii].data_var;
    outstats->prof_avg += instats[ii].prof_avg;
    outstats->prof_var += instats[ii].prof_var;
  }

  /* Profile information gets added together, but */
  /* data set info gets averaged together.        */

  outstats->data_avg /= numprofs;
  outstats->data_var /= numprofs;

  /* Calculate the reduced chi-squared */

  outstats->redchi = chisqr(outprof, proflen, outstats->prof_avg, 
			    outstats->prof_var) / (proflen - 1.0);
}


void initialize_foldstats(foldstats *stats)
/* Zeroize all of the components of stats */
{
  stats->numdata = 0.0;   /* Number of data bins folded         */
  stats->data_avg = 0.0;  /* Average level of the data bins     */
  stats->data_var = 0.0;  /* Variance of the data bins          */
  stats->numprof = 0.0;   /* Number of bins in the profile      */
  stats->prof_avg = 0.0;  /* Average level of the profile bins  */
  stats->prof_var = 0.0;  /* Variance of the profile bins       */
  stats->redchi = 0.0;    /* Reduced chi-squared of the profile */
}
