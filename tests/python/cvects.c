#include <cvects.h>

float *complex_arr(long N)
{
  float *arr;
  long i;

  arr = gen_cvect(N);
  for (i=0; i<2*N; i++)
    arr[i] = (float) i;
  return arr;
}

void mult_arr(float *arr, float val, long N)
{
  long i;

  for (i=0; i<2*N; i+=2)
    arr[i] *= val;
}


void dgenrotate_1d(double *data, long numbins, double bins_to_left)
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */
{
    double *tmp, lopart, hipart, intpart;
    long i, index;

    bins_to_left = fmod(bins_to_left, numbins);
    if (bins_to_left < 0.0) bins_to_left += numbins;
    tmp = gen_dvect(numbins);
    lopart = modf(bins_to_left, &intpart);
    hipart = 1.0 - lopart;
    index = (long) floor(intpart + 1.0e-20);
    for (i = 0 ; i < numbins ; i++){
      tmp[i] = hipart * data[(index + i) % numbins] + \
        lopart * data[(index + i + 1) % numbins];
    }
    memcpy(data, tmp, sizeof(double) * numbins);
    free(tmp);
}


void dstats(double *x, int n, double *mean, double *var, 
            double *skew, double *kurt)
/* For a double precision vector, *x, of length n, this routine  */
/* returns the mean, variance, skewness, and kurtosis of *x.     */
{
  long i;
  double an=0.0, an1=0.0, dx, t1, t2, stdevmle;
  
  /*  Modified (29 June 98) C version of the following:        */
  /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
  /*  Returned values were checked with Mathematica 3.01       */

  if (n < 1) {
    printf("\vVector length must be > 0 in moments().  Exiting\n");
    exit(1);
  } else {
    *mean = (double) x[0];
    *var = 0.0;
    *skew = 0.0;
    *kurt = 0.0;
  }
  
  for (i = 1 ; i < n ; i++){
    an = (double) (i + 1);
    an1 = (double) (i);
    dx = (x[i] - *mean) / an;
    t1 = an1 * dx * dx;
    t2 = an * t1;
    *kurt -= dx * (4.0 * *skew - dx * \
                   (6.0 * *var + t1 * (1.0 + an1 * an1 * an1)));
    *skew -= dx * (3.0 * *var - t2 * (an - 2.0));
    *var += t2;
    *mean += dx;
  }

  if (n > 1){
    stdevmle = sqrt(*var/an);
    t1 = stdevmle * stdevmle * stdevmle;
    *skew /= t1 * an;
    *kurt = *kurt / (t1 * stdevmle * an) - 3.0;
    *var /= an1;
  }

  return;
}


long get_filelen(FILE *file)
/* Return the length of the file *file in bytes without */
/* losing your place in the file.                       */
{
  int errorval;
  long tmploc, returnval;

  tmploc = ftell(file);
  if ((errorval = fseek(file, 0L, SEEK_END))){
    printf("\nError %d in fseek().  Exiting.\n\n",errorval);
  }
  returnval = ftell(file);
  if ((errorval = fseek(file, tmploc, SEEK_SET))){
    printf("\nError %d in fseek().  Exiting.\n\n",errorval);
  }
  return returnval;
}
