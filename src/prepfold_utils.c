#include "prepfold.h"

int read_floats(FILE *file, float *data, int numpts,
		double *dispdelays, int numsubbands, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
  /* The following 2 lines just get rid of some compiler warnings */

  *dispdelays = *dispdelays;
  numsubbands = numsubbands;

  /* Read the raw data and return numbar read */

  return chkfread(data, sizeof(float),
		  (unsigned long) (numpts * numchan), file) / numchan;
}


void fold_errors(double *prof, int proflen, double dt, double N, 
		 double datavar, double p, double pd, double pdd, 
		 double *perr, double *pderr, double *pdderr)
/* Calculate estimates for the errors in period p-dot and   */
/* p-dotdot using Middleditch's error formula.  The routine */
/* calculates the errors for each Fourier harmonic present  */
/* in the profile that is significant.  Then it combines    */
/* the errors for the harmonics into an error for the       */
/* fundamental.                                             */
/*   Arguments:                                             */
/*      'prof' is and array pointing to the profile         */
/*      'proflen' is the number of bins in 'prof'           */
/*      'dt' is the sample interval of the original data    */
/*      'N' is the total number of points folded            */
/*      'datavar' is the variance of the original data      */
/*      'p' is the folding period                           */
/*      'pd' is the folding period derivative               */
/*      'pdd' is the folding period 2nd dervivative         */
/*      'perr' is the returned period error                 */
/*      'pderr' is the returned p-dot error                 */
/*      'pdderr' is the returned p-dotdot error             */
{
  int ii, gotone=0;
  double T, T2, pwr, norm, sigpow=6.6, r2, r4, z2, sr2, sz2;
  double dtmp, r, z, w, pwrfact=0.0, rerr, zerr, werr;
  double rerrn=0.0, zerrn=0.0, werrn=0.0, rerrd=0.0, zerrd=0.0, werrd=0.0;
  float powargr, powargi;
  fcomplex *fftprof;

  /* Total length in time of data set */

  T = N * dt;

  /* Convert p, pd, and pdd into r, z, and w */

  dtmp = p * p;
  T2 = T * T;
  r = T / p;
  z = T2 * -pd / dtmp;
  if (pdd==0.0)
    w = 0.0;
  else 
    w = T2 * T * (2.0 * pd * pd / (dtmp * p) - pdd / dtmp);

  /* Calculate the normalization constant which converts the raw */
  /* powers into normalized powers -- just as if we had FFTd the */
  /* full data set.                                              */

  norm = 1.0 / (N * datavar);

  /* Place the profile into a complex array */

  fftprof = gen_cvect(proflen);
  for (ii = 0; ii < proflen; ii++){
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

  ii = 1;
  for (ii = 1; ii < proflen / 2; ii++){
    pwr = POWER(fftprof[ii].r, fftprof[ii].i) * norm;
    pwrfact = 2.0 * pwr;
    if (pwr > sigpow){
      gotone = 1;
      dtmp = 3.0 / ((PI * sqrt(6.0 * pwr)) * ii);
      rerrn += pwrfact / dtmp;
      rerrd += pwrfact / (dtmp * dtmp);
      dtmp = 3.0 * sqrt(10.0) / ((PI * sqrt(pwr)) * ii);
      zerrn += pwrfact / dtmp;
      zerrd += pwrfact / (dtmp * dtmp);
      dtmp = 6.0 * sqrt(105.0) / ((PI * sqrt(pwr) * ii));
      werrn += pwrfact / dtmp;
      werrd += pwrfact / (dtmp * dtmp);
    }
  }

  /* Help protect against really low significance profiles...*/

  if (!gotone){
    rerr = zerr = werr = 0.5;
  }

  /* Calculate the standard deviations */

  rerr = rerrn / rerrd;
  zerr = zerrn / zerrd;
  werr = werrn / werrd;

  /* Some useful values */

  r2 = r * r;
  sr2 = rerr * rerr;
  r4 = r2 * r2;
  z2 = z * z;
  sz2 = zerr * zerr;
  dtmp = r * w - 3 * z2;

  /* Convert the standard deviations to periods */
  
  *perr = T * rerr / r2;
  *pderr = sqrt(4 * z2 * sr2 / (r4 * r2) + sz2 / r4);
  *pdderr = sqrt((werr * werr * r4 + 16 * sz2 * r2 + \
		  4 * dtmp * dtmp * sr2) / (r4 * r4 * T2));

  /* Free our FFT array */

  free(fftprof);
}


int bary2topo(double *topotimes, double *barytimes, int numtimes, 
	      double fb, double fbd, double fbdd, 
	      double *ft, double *ftd, double *ftdd)
/* Convert a set of barycentric pulsar spin parameters (fb, fbd, fbdd) */
/* into topocentric spin parameters (ft, ftd, ftdd) by performing      */
/* a linear least-squares fit (using LAPACK routine DGELS).  The       */
/* routine equates the pulse phase using topcentric parameters and     */
/* times to the pulse phase using barycentric parameters and times.    */
{
  double *work, *aa, *bb, dtmp;
  int ii, mm=3, nn, nrhs=1, lwork, info, index;
  char trans='T';

  if (numtimes < 4){
    printf("\n'numtimes' < 4 in bary2topo():  Cannot solve.\n\n");
    exit(0);
  }
  nn = numtimes; 
  lwork = mm + nn * 9;
  aa = gen_dvect(mm * nn);
  bb = gen_dvect(nn);
  work = gen_dvect(lwork);
  for (ii = 0; ii < nn; ii++){
    index = ii * 3;
    dtmp = (topotimes[ii] - topotimes[0]) * SECPERDAY;
    aa[index] = dtmp;
    aa[index+1] = 0.5 * dtmp * dtmp;
    aa[index+2] = dtmp * dtmp * dtmp / 6.0;
    dtmp = (barytimes[ii] - barytimes[0]) * SECPERDAY;
    bb[ii] = dtmp * (fb + dtmp * (0.5 * fbd + fbdd * dtmp / 6.0));
  }
  dgels_(&trans, &mm, &nn, &nrhs, aa, &mm, bb, &nn, work, &lwork, &info);
  *ft = bb[0];
  *ftd = bb[1];
  *ftdd = bb[2];
  free(aa);
  free(bb);
  free(work);
  return info;
}


void init_prepfoldinfo(prepfoldinfo *in)
/* Set all values to 0 or NULL */
{
  in->rawfolds = NULL;
  in->dms = NULL;
  in->periods = NULL;
  in->pdots = NULL;
  in->stats = NULL;
  in->numdms = 0;
  in->numperiods = 0;
  in->numpdots = 0;
  in->nsub = 0;
  in->npart = 0;
  in->proflen = 0;
  in->numchan = 1;
  in->ndmfact = 2;
  in->npfact = 1;
  in->step = 1;
  in->filenm = NULL;
  in->candnm = NULL;
  in->telescope = NULL;
  in->pgdev = NULL;
  in->dt = 0.0;
  in->tepoch = 0.0;
  in->bepoch = 0.0;
  in->avgvoverc = 0.0;
  in->lofreq = 0.0;
  in->chan_wid = 0.0;
  in->bestdm = 0.0;
  in->topo.pow = 0.0;
  in->topo.p1 = 0.0;
  in->topo.p2 = 0.0;
  in->topo.p3 = 0.0;
  in->bary.pow = 0.0;
  in->bary.p1 = 0.0;
  in->bary.p2 = 0.0;
  in->bary.p3 = 0.0;
  in->fold.pow = 0.0;
  in->fold.p1 = 0.0;
  in->fold.p2 = 0.0;
  in->fold.p3 = 0.0;
  in->orb.p = 0.0;
  in->orb.e = 0.0;
  in->orb.x = 0.0;
  in->orb.w = 0.0;
  in->orb.t = 0.0;
  in->orb.pd = 0.0;
  in->orb.wd = 0.0;
}

void delete_prepfoldinfo(prepfoldinfo *in)
/* Free all dynamic arrays in the prepfold array */
{
  free(in->rawfolds);
  if (in->nsub > 1) free(in->dms);
  free(in->periods);
  free(in->pdots);
  free(in->stats);
  free(in->filenm);
  free(in->candnm);
  free(in->telescope);
  free(in->pgdev);
}

void double2float(double *in, float *out, int numpts)
/* Copy a double vector into a float vector */
{
  int ii;

  for (ii = 0; ii < numpts; ii++)
    out[ii] = (float) in[ii];
}


int cpgnice_output_2(char *out, double val, double err, int len)
/* Return a string that has the same formatting as       */
/* nice_output_2(), but for PGPLOT.  This way, exponents */
/* are actually in superscript!  Woo-hoo!                */
{
  int inlen, goodlen;
  char tempout[100], part[20], expon[20], *expptr;

  inlen = nice_output_2(tempout, val, err, len);

  /* See if the formatted output*/

  expptr = strrchr(tempout, '^');

  /* Not in scientific notation */

  if (expptr==NULL){
    strcpy(out, tempout);
    return inlen;
  } else {
    goodlen = expptr - tempout;
    strcpy(expon, tempout + goodlen + 1);
    strncpy(part, tempout, goodlen);
    part[goodlen] = 0;
    sprintf(out, "%s\\u%s\\d", part, expon);
    return strlen(out);
  }
}
