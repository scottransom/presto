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
		 double *perr, double *pderr, double *pdderr){
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
/*      'pderr' is the returned p-dotdot error              */
  int ii;
  double T, pwr, norm, hferr, hfderr, hfdderr;
  double r, z, w;
  fcomplex *fftprof;

  /* Total length in time of data set */

  T = datanum * dt;

  /* Convert p, pd, and pdd into r, z, and w */

  r = T / p;
  z = -pd * r * r;

  /* calculate the normalization constant which converts the raw */
  /* powers into normalized powers -- just as if we had FFTd the */
  /* full data set.                                              */

  norm = 1.0 / (datanum * datavar);

  /* Place the profile into a complex array */

  fftprof = gen_cvect(proflen);
  for (ii = 0; ii < proflen; ii++){
    fftprof[ii].r = (float) prof[ii];
    fftprof[ii].r = 0.0;
  }

  /* FFT the profile */

  COMPLEXFFT(fftprof, proflen, -1);

  /* Step through the powers and find the significant ones */

  for (ii = 0; ii < proflen; ii++){
    pwr = POWER(fftprof[ii].r, fftprof[ii].i) * norm;
    if (pwr > 3.0){
      hferr
    }
  }
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
  in->filenm = NULL;
  in->candnm = NULL;
  in->telescope = NULL;
  in->pgdev = NULL;
  in->tepoch = 0.0;
  in->bepoch = 0.0;
  in->dt = 0.0;
  in->bestdm = 0.0;
  in->topo.pow = 0.0;
  in->topo.p1 = 0.0;
  in->topo.p2 = 0.0;
  in->topo.p3 = 0.0;
  in->bary.pow = 0.0;
  in->bary.p1 = 0.0;
  in->bary.p2 = 0.0;
  in->bary.p3 = 0.0;
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
  free(in->dms);
  free(in->periods);
  free(in->pdots);
  free(in->stats);
  free(in->filenm);
  free(in->candnm);
  free(in->telescope);
  free(in->pgdev);
}
