#include "presto.h"
#include "float.h"

/* The number of points to work with at a time         */
/* This must be the same as the WORKLEN in profile.c!  */
#define WORKLEN 16384

/* Some macros to make the flag checking easier */
#define DELAYS (flags % 2 == 1)
#define ONOFF (flags > 1)

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) \
    ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Macro to test if a float is equal to 1 */
#define TEST_ONE(a) (fabs((a) - 1.0) <= 2 * DBL_EPSILON ? 1 : 0)

/* Macro to add data to the buffer */
#define ADDTOBUFFER(part, bin, phase)           \
    newphaseadded = *phaseadded + partphase;    \
    if (newphaseadded > 1.0){                   \
      dump_buffer(prof, buffer, numprof, phase, \
		  partphase, part, phaseadded); \
    } else {                                    \
      *phaseadded = newphaseadded;              \
      buffer[bin] += part;                      \
    }

/* Macro to add data directly to the profile */
/* (i.e. this is the old method)             */
/* #define ADDTOBUFFER(part, bin, phase) \ */
/*     prof[bin] += part; */

/* Some helper functions */

void hunt(double *xx, int n, double x, int *jlo);
static void add_to_prof(double *prof, double *buffer, int numprof,
                        double lophase, double deltaphase,
                        double dataval, double *phaseadded);
static void dump_buffer(double *prof, double *buffer, int numprof,
                        double lophase, double partphase,
                        double dataval, double *phaseadded);


static void dump_buffer(double *prof, double *buffer, int numprof,
                        double lophase, double partphase,
                        double dataval, double *phaseadded)
/* Fill the profile buffer and transfer its contents */
/* to the real profile.  Reset the buffer.           */
{
   int ii;
   double addpart, deltaphase;

   /* Check for bad input values */

   if (*phaseadded < 0.0 || *phaseadded > 1.0)
      printf("Ack!  phaseadded is %17.15g in dump_buffer()!\n", *phaseadded);
   if (lophase >= 1.0)
      lophase -= (int) lophase;

   /* Squeeze in the last remaining bit of data we can */

   deltaphase = 1.0 - *phaseadded;
   addpart = dataval * deltaphase / partphase;
   add_to_prof(prof, buffer, numprof, lophase, deltaphase, addpart, phaseadded);

   /* Dump the buffer into the profile array */
   /* and reset the buffer array to zeros.   */

   for (ii = 0; ii < numprof; ii++) {
      prof[ii] += buffer[ii];
      buffer[ii] = 0.0;
   }

   /* Add the remaining bit of data to the fresh buffer */

   lophase += deltaphase;
   if (lophase >= 1.0)
      lophase -= (int) lophase;
   deltaphase = partphase - deltaphase;
   addpart = dataval - addpart;
   add_to_prof(prof, buffer, numprof, lophase, deltaphase, addpart, phaseadded);

   /* Return the new value of phaseadded */

   *phaseadded = deltaphase;
}


static void add_to_prof(double *prof, double *buffer, int numprof,
                        double lophase, double deltaphase,
                        double dataval, double *phaseadded)
/* Note:  phaseadded should be set to 0.0 the first time this */
/*        function is called.                                 */
{
   int ii, numbins, loprofbin, hiprofbin;
   double dtmp, part, partphase, profbinwidth;
   double newphaseadded, phaseleft, hiphase, phase;

   /* Reset the phase counter if needed */

   if (TEST_ONE(*phaseadded))
      *phaseadded = 0.0;

   /* Dump (and reset) the buffer if add_to_prof() gets called with lophase==1 */

   if (TEST_ONE(lophase)) {
      for (ii = 0; ii < numprof; ii++) {
         prof[ii] += buffer[ii];
         buffer[ii] = 0.0;
      }
      lophase = 0.0;
   }

   /* Find the bins we will add data to.                 */
   /* Note:  hiprofbin will be used modulo numprof so    */
   /*        it could be greater than numprof.           */

   hiphase = lophase + deltaphase;
   loprofbin = (int) (lophase * numprof);
   hiprofbin = (int) (hiphase * numprof);

   /* How many profile bins we will spread the data over? */

   numbins = hiprofbin - loprofbin + 1;

   /* Spread the data into the proper bins. */

   if (numbins >= 3) {

      /* Data point will be spread over 3 or more profile bins */

      profbinwidth = 1.0 / numprof;
      dtmp = dataval / deltaphase;

      /* Fill the current (low) buffer bin */
      partphase = (loprofbin + 1) * profbinwidth - lophase;
      part = dtmp * partphase;
      ADDTOBUFFER(part, loprofbin, lophase)
          phase = lophase + partphase;
      phaseleft = deltaphase - partphase;

      /* Fill the middle buffer bins */
      partphase = profbinwidth;
      part = dtmp * partphase;
      for (ii = loprofbin + 1; ii < hiprofbin; ii++) {
         ADDTOBUFFER(part, ii % numprof, phase)
             phase += partphase;
      }

      /* Add what's left into the high beffer bin */
      partphase = hiphase - phase;
      part = dtmp * partphase;
      ADDTOBUFFER(part, hiprofbin % numprof, phase)

   } else if (numbins == 2) {

      /* Data point will be spread over 2 profile bins */

      /* Low bin */
      partphase = (double) hiprofbin / numprof - lophase;
      part = dataval * partphase / deltaphase;
      ADDTOBUFFER(part, loprofbin, lophase)
          phase = lophase + partphase;

      /* High bin */
      partphase = deltaphase - partphase;
      part = dataval - part;
      ADDTOBUFFER(part, hiprofbin % numprof, phase)

   } else {

      /* Data point will go into only 1 profile bin */

      partphase = deltaphase;
      part = dataval;
      ADDTOBUFFER(part, loprofbin, lophase)
   }
}


void fold_errors(double *prof, int proflen, double dt, double N,
                 double datavar, double p, double pd, double pdd,
                 double *perr, double *pderr, double *pdderr)
/* Calculate estimates for the errors in period p-dot and    */
/* p-dotdot using Middleditch's error formula.  The routine  */
/* calculates the errors for each Fourier harmonic present   */
/* in the profile that is significant.  Then it combines     */
/* the errors for the harmonics into an error for the        */
/* fundamental.                                              */
/*   Arguments:                                              */
/*      'prof' is and array pointing to the profile          */
/*      'proflen' is the number of bins in 'prof'            */
/*      'dt' is the sample interval of the original data     */
/*      'N' is the total number of points folded             */
/*      'datavar' is the variance of the original data       */
/*      'p' is the folding period                            */
/*      'pd' is the folding period derivative                */
/*      'pdd' is the folding period 2nd dervivative          */
/*      'perr' is the returned period standard deviation     */
/*      'pderr' is the returned p-dot standard deviation     */
/*      'pdderr' is the returned p-dotdot standard deviation */
{
   int ii, gotone = 0;
   double T, pwr, norm, sigpow = 2.7, r2, r4, z2, sr2, sz2, errvar = 0.0;
   double err, r, z, w, pwrfact = 0.0, pwrerr = 0.0, rerr, zerr, werr, dtmp;
   double rerrn = 0.0, zerrn = 0.0, werrn = 0.0, rerrd = 0.0, zerrd = 0.0, werrd =
       0.0;
   float powargr, powargi;
   fcomplex *fftprof;

   /* Total length in time of data set */

   T = N * dt;

   /* Convert p, pd, and pdd into r, z, and w */

   r = T / p;
   z = -pd * r * r;
   if (pdd == 0.0)
      w = 0.0;
   else
      w = 2.0 * z * z / r - pdd * r * r * T;

   /* Calculate the normalization constant which converts the raw */
   /* powers into normalized powers -- just as if we had FFTd the */
   /* full data set.                                              */

   norm = 1.0 / (N * datavar);

   /* Place the profile into a complex array */

   fftprof = gen_cvect(proflen);
   for (ii = 0; ii < proflen; ii++) {
      fftprof[ii].r = (float) prof[ii];
      fftprof[ii].i = 0.0;
   }

   /* FFT the profile */

   COMPLEXFFT(fftprof, proflen, -1);

   /* Step through the powers and find the significant ones.  */
   /* Estimate the error of the fundamental using each one.   */
   /* Combine these errors into a unified error of the freq.  */
   /* Note:  In our case the errors are the data points and   */
   /*        we are combining them using a weighted mean.     */
   /*        The weights come from the fact that the powers   */
   /*        have a measurements error = sqrt(2 * P).  This   */
   /*        causes an error in our estimates of rerr.        */

   for (ii = 1; ii < proflen / 2; ii++) {
      pwr = POWER(fftprof[ii].r, fftprof[ii].i) * norm;
      pwrerr = sqrt(2.0 * pwr);
      pwrfact = 1.0 / (sqrt(pwr) * ii);
      if (pwr > sigpow) {
         gotone = 1;
         /* r error */
         err = 0.38984840062 * pwrfact;
         errvar = err / pwrerr;
         errvar = errvar * errvar;
         rerrn += err / errvar;
         rerrd += 1.0 / errvar;
         /* z error */
         err = 3.01975272627 * pwrfact;
         errvar = err / pwrerr;
         errvar = errvar * errvar;
         zerrn += err / errvar;
         zerrd += 1.0 / errvar;
         /* w error */
         err = 19.5702343923 * pwrfact;
         errvar = err / pwrerr;
         errvar = errvar * errvar;
         werrn += err / errvar;
         werrd += 1.0 / errvar;
      }
   }

   if (gotone) {

      /* Calculate the standard deviations */

      rerr = rerrn / rerrd;
      zerr = zerrn / zerrd;
      werr = werrn / werrd;

      /* Help protect against really low significance profiles.  */
      /* And note that this is probably _underestimating_ the    */
      /* errors in this case...                                  */

   } else {
      rerr = 0.5;
      zerr = 7.8;
      werr = 50.2;
   }

   /* Some useful values */

   r2 = r * r;
   sr2 = rerr * rerr;
   r4 = r2 * r2;
   z2 = z * z;
   sz2 = zerr * zerr;
   dtmp = r * w - 3.0 * z2;

   /* Convert the standard deviations to periods */

   *perr = T * rerr / r2;
   *pderr = sqrt(4.0 * z2 * sr2 / (r4 * r2) + sz2 / r4);
   *pdderr = sqrt((werr * werr * r4 + 16 * z2 * sz2 * r2 +
                   4.0 * dtmp * dtmp * sr2) / (r4 * r4 * T * T));

   /* Free our FFT array */

   vect_free(fftprof);
}


double foldfile(FILE * datafile, double dt, double tlo,
                double *prof, int numprof, double startphs,
                double fo, double fdot, double fdotdot, int flags,
                double *delays, double *delaytimes, int numdelays,
                double *onoffpairs, foldstats * stats, float *chiarr)
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
   double *onoffptr = NULL, *buffer, phase = 0.0, phaseadded = 0.0;
   int ourflags;
   unsigned long ii, N, onbin, offbin, numbins;
   unsigned long remainbins, binstoread, numreads;

   /* Get the data file length and initialize some variables */

   N = chkfilelen(datafile, sizeof(float));
   if (ONOFF)
      onoffptr = onoffpairs;
   stats->numdata = stats->data_avg = stats->data_var = 0.0;
   if (DELAYS)
      ourflags = 1;
   else
      ourflags = 0;

   /* Create and initialize the buffer needed by fold() */

   buffer = gen_dvect(numprof);
   for (ii = 0; ii < (unsigned long) numprof; ii++)
      buffer[ii] = 0.0;

   do {                         /* Loop over the on-off pairs */

      /* Set the on-off variables */

      if (ONOFF) {
         onbin = (unsigned long) (*onoffptr * N + DBLCORRECT);
         offbin = (unsigned long) (*(onoffptr + 1) * N + DBLCORRECT);
         if (offbin)
            offbin--;
         onoffptr += 2;
      } else {
         onbin = 0;
         offbin = N - 1;
      }
      numbins = (offbin == onbin) ? 0 : offbin - onbin + 1;
      numreads = numbins / WORKLEN;
      remainbins = numbins % WORKLEN;
      if (remainbins)
         numreads++;
      binstoread = WORKLEN;

      /* Skip to the correct file location */

      chkfileseek(datafile, onbin, sizeof(float), SEEK_SET);

      /* Loop over the number of reads we have to perform for */
      /* the current on-off pair.                             */

      for (ii = 0; ii < numreads; ii++, onbin += binstoread) {

         /* Correct for the fact that our last read might be short */

         if (remainbins && (ii == numreads - 1))
            binstoread = remainbins;

         /* Read the current chunk of data */

         chkfread(data, sizeof(float), binstoread, datafile);

         /* Fold the current chunk of data */

         phase = fold(data, binstoread, dt, tlo + onbin * dt, prof,
                      numprof, startphs, buffer, &phaseadded,
                      fo, fdot, fdotdot, ourflags,
                      delays, delaytimes, numdelays, NULL, stats);

         /* Set the current chiarr value */

         chiarr[onbin / WORKLEN] = (float) stats->redchi;
      }

   } while (offbin < N - 1 && offbin != 0);

   /* Free the buffer and return */

   vect_free(buffer);
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
   loprofphase = (phase < 0.0) ? 1.0 + phase - (int) phase : phase - (int) phase;
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

      if (numbins >= 3) {

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
            double *buffer, double *phaseadded,
            double fo, double fdot, double fdotdot, int flags,
            double *delays, double *delaytimes, int numdelays,
            int *onoffpairs, foldstats * stats)
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
/*    'buffer' is a double prec array of numprof values containing    */
/*            data that hasn't made it into the prof yet.             */
/*    'phaseadded' is the address to a variable showing how much      */
/*            has been added to the buffer [0-1] (must start as 0.0)  */
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
   int ii, onbin, offbin, *onoffptr = NULL;
   int arrayoffset = 0;
   double phase, phasenext = 0.0, deltaphase, T, Tnext, TD, TDnext;
   double profbinwidth, lophase, hiphase, dtmp;
   double dev, delaytlo = 0.0, delaythi = 0.0, delaylo = 0.0, delayhi = 0.0;
   double *delayptr = NULL, *delaytimeptr = NULL;

   /* Initialize some variables and save some FLOPs later... */

   fdot /= 2.0;
   fdotdot /= 6.0;
   profbinwidth = 1.0 / numprof;
   if (ONOFF)
      onoffptr = onoffpairs;
   stats->numprof = (double) numprof;
   stats->data_var *= (stats->numdata - 1.0);

   do {                         /* Loop over the on-off pairs */

      /* Set the on-off pointers and variables */

      if (ONOFF) {
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

      if (DELAYS) {

         /* Guess that the next delay we want is the next available */

         arrayoffset += 2;      /* Beware nasty NR zero-offset kludges! */
         hunt(delaytimes - 1, numdelays, T, &arrayoffset);
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
      lophase = (phase < 0.0) ? 1.0 + phase - (int) phase : phase - (int) phase;

      /* Generate the profile for this onoff pair */

      for (ii = onbin; ii <= offbin; ii++) {

         /* Calculate the barycentric time for the next point. */

         Tnext = tlo + (ii + 1) * dt;
         TDnext = Tnext;

         /* Set the delay pointers and variables */

         if (DELAYS) {
            if (Tnext > delaythi) {

               /* Guess that the next delay we want is the next available */

               arrayoffset += 2;        /* Beware nasty NR zero-offset kludges! */
               hunt(delaytimes - 1, numdelays, Tnext, &arrayoffset);
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

         /* Add the current point to the buffer or the profile */

         add_to_prof(prof, buffer, numprof, lophase,
                     deltaphase, data[ii], phaseadded);

         /* Update variables */

         hiphase = lophase + deltaphase;
         lophase = hiphase - (int) hiphase;
         phase = phasenext;

         /* Use clever single pass mean and variance calculation */

         stats->numdata += 1.0;
         dev = data[ii] - stats->data_avg;
         stats->data_avg += dev / stats->numdata;
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
   for (ii = 0; ii < numprof; ii++) {
      dtmp = prof[ii] - stats->prof_avg;
      stats->redchi += dtmp * dtmp;
   }
   stats->data_var /= (stats->numdata - 1.0);
   stats->prof_var = stats->data_var * stats->numdata * profbinwidth;
   stats->redchi /= (stats->prof_var * (numprof - 1));

   phasenext = (phasenext < 0.0) ?
       1.0 + phasenext - (int) phasenext : phasenext - (int) phasenext;

   return (phasenext);
}

#undef WORKLEN
#undef DELAYS
#undef ONOFF
#undef LININTERP
#undef TEST_ONE


void shift_prof(double *prof, int proflen, int shift, double *outprof)
/* Rotates a profile 'prof' by an integer 'shift' places.    */
/* If 'shift' < 0 then shift left, 'shift' > 0, shift right. */
/* Place the shifted  profile in 'outprof'.                  */
{
   int wrap = 0, nowrap = 0;

   wrap = shift % proflen;

   if (prof == outprof) {
      double *tmpprof;
      /* no shift */
      if (wrap == 0) {
         return;
         /* Convert a left shift into the equivalent right shift */
      } else if (wrap < 0)
         wrap += proflen;
      /* Perform a right shift */
      nowrap = proflen - wrap;
      tmpprof = gen_dvect(proflen);
      memcpy(tmpprof, prof + nowrap, wrap * sizeof(double));
      memcpy(tmpprof + wrap, prof, nowrap * sizeof(double));
      memcpy(outprof, tmpprof, proflen * sizeof(double));
      vect_free(tmpprof);
   } else {
      /* no shift */
      if (wrap == 0) {
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
}


void combine_profs(double *profs, foldstats * instats, int numprofs,
                   int proflen, double *delays, double *outprof,
                   foldstats * outstats)
/* Combine a series of 'numprofs' profiles, each of length 'proflen',   */
/* into a single profile of length 'proflen'.  The profiles are         */
/* summed after the appropriate 'delays' are added to each profile.     */
/* The result is a profile in 'outprof' (which must be pre-allocated)   */
/* The input stats in 'instats' are combined and placed in 'outstats'   */
{
   int ii, jj, kk, index = 0, offset;
   double *local_delays;

   /* Initiate the output statistics */
   initialize_foldstats(outstats);
   outstats->numprof = proflen;
   local_delays = gen_dvect(numprofs);

   /* Convert all the delays to positive offsets from   */
   /* the phase=0 profile bin, in units of profile bins */
   /* Note:  The negative sign refers to the fact that  */
   /*        we want positiev numbers to represent      */
   /*        shifts _to_ the right not _from_ the right */

   for (ii = 0; ii < numprofs; ii++) {
      local_delays[ii] = fmod(-delays[ii], proflen);
      if (local_delays[ii] < 0.0)
         local_delays[ii] += proflen;
   }

   /* Set the output array to zeros */
   for (ii = 0; ii < proflen; ii++)
      outprof[ii] = 0.0;

   /* Loop over the profiles */
   for (ii = 0; ii < numprofs; ii++) {

      /* Calculate the appropriate offset into the profile array */
      offset = (int) (local_delays[ii] + 0.5);

      /* Sum the profiles */
      for (jj = 0, kk = proflen - offset; jj < offset; jj++, kk++, index++)
         outprof[kk] += profs[index];
      for (kk = 0; jj < proflen; jj++, kk++, index++)
         outprof[kk] += profs[index];

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
   vect_free(local_delays);
}


void initialize_foldstats(foldstats * stats)
/* Zeroize all of the components of stats */
{
   stats->numdata = 0.0;        /* Number of data bins folded         */
   stats->data_avg = 0.0;       /* Average level of the data bins     */
   stats->data_var = 0.0;       /* Variance of the data bins          */
   stats->numprof = 0.0;        /* Number of bins in the profile      */
   stats->prof_avg = 0.0;       /* Average level of the profile bins  */
   stats->prof_var = 0.0;       /* Variance of the profile bins       */
   stats->redchi = 0.0;         /* Reduced chi-squared of the profile */
}
