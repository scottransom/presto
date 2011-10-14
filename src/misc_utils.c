#include <ctype.h>
#include "misc_utils.h"

#ifndef ARCSEC2RAD
/* arcseconds to radians */
#define ARCSEC2RAD 4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
/* seconds of time to radians */
#define SEC2RAD    7.2722052166430399038487115353692196393452995355905e-5
#endif

/* return the max of two double values */
#define DMAX(a,b) ((a)>(b)?(a):(b))

/* like strcpy, except guaranteed to work with overlapping strings      */
#define strMove(d,s) memmove(d,s,strlen(s)+1)


char *rmtrail(char *str)
/* Removes trailing space from a string */
{
   int i;

   if (str && 0 != (i = strlen(str))) {
      while (--i >= 0) {
         if (!isspace(str[i]))
            break;
      }
      str[++i] = '\0';
   }
   return str;
}

char *rmlead(char *str)
/* Removes leading space from a string */
{
   char *obuf;

   if (str) {
      for (obuf = str; *obuf && isspace(*obuf); ++obuf);
      if (str != obuf)
         strMove(str, obuf);
   }
   return str;
}

char *remove_whitespace(char *str)
/* Remove leading and trailing space from a string */
{
   return rmlead(rmtrail(str));
}

char *strlower(char *str)
/* Convert a string to lower case */
{
   char *ss;

   if (str) {
      for (ss = str; *ss; ++ss)
         *ss = tolower(*ss);
   }
   return str;
}


void split_path_file(char *input, char **path, char **file)
/* This routine splits an input string into a path and */
/* a filename.  Since it allocates the memory for the  */
/* path and filename dynamically, the calling program  */
/* must free both "path" and "file".                   */
{
   char *sptr = NULL, stmp[200];
   unsigned int len, pathlen = 0, filelen = 0;

   len = strlen(input);
   sptr = strrchr(input, '/');
   if (sptr == NULL) {
      getcwd(stmp, 200);
      if (stmp == NULL) {
         printf("\nCurrent directory name is too long.\n");
         printf("Exiting\n\n");
         exit(1);
      }
      pathlen = strlen(stmp);
      *path = (char *) calloc(pathlen + 1, sizeof(char));
      *file = (char *) calloc(len + 1, sizeof(char));
      strcpy(*path, stmp);
      strncpy(*file, input, len);
   } else {
      pathlen = sptr - input;
      filelen = len - pathlen - 1;
      *path = (char *) calloc(pathlen + 1, sizeof(char));
      *file = (char *) calloc(filelen + 1, sizeof(char));
      strncpy(*path, input, pathlen);
      strncpy(*file, sptr + 1, filelen);
   }
}

int split_root_suffix(char *input, char **root, char **suffix)
/* This routine splits an input string into a root name */
/* + suffix.  Since it allocates the memory for the     */
/* root and suffix dynamically, the calling program     */
/* must free both "root" and "suffix".                  */
/* If the routine finds a suffix, it returns 1, else 0. */
{
   char *sptr = NULL;
   unsigned int len, rootlen = 0, suffixlen = 0;

   len = strlen(input);
   sptr = strrchr(input, '.');
   if (sptr == NULL) {
      *root = (char *) calloc(len + 1, sizeof(char));
      strncpy(*root, input, len);
      return 0;
   } else {
      rootlen = sptr - input;
      *root = (char *) calloc(rootlen + 1, sizeof(char));
      strncpy(*root, input, rootlen);
      suffixlen = len - rootlen - 1;
      *suffix = (char *) calloc(suffixlen + 1, sizeof(char));
      strncpy(*suffix, sptr + 1, suffixlen);
      return 1;
   }
}

void strtofilename(char *string)
/* Trim spaces off the end of *input and convert */
/* all other spaces into underscores.            */
{
   int ii;

   ii = strlen(string) - 1;
   do {
      if (string[ii] == ' ')
         string[ii] = '\0';
      else
         break;
   } while (ii--);
   do {
      if (string[ii] == ' ')
         string[ii] = '_';
   } while (ii--);
}


void telescope_to_tempocode(char *inname, char *outname, char*obscode)
// Return the 2 character TEMPO string for an observatory
// whose name is in the string "inname".  Return a nice
// name in "outname".
{
    char scope[40];
    
    strncpy(scope, inname, 40);
    strlower(scope);
    if (strcmp(scope, "gbt") == 0) {
        strcpy(obscode, "GB");
        strcpy(outname, "GBT");
    } else if (strcmp(scope, "arecibo") == 0) {
        strcpy(obscode, "AO");
        strcpy(outname, "Arecibo");
    } else if (strcmp(scope, "vla") == 0) {
        strcpy(obscode, "VL");
        strcpy(outname, "VLA");
    } else if (strcmp(scope, "parkes") == 0) {
        strcpy(obscode, "PK");
        strcpy(outname, "Parkes");
    } else if (strcmp(scope, "jodrell") == 0) {
        strcpy(obscode, "JB");
        strcpy(outname, "Jodrell Bank");
    } else if ((strcmp(scope, "gb43m") == 0) ||
               (strcmp(scope, "gb 140ft") == 0)){
        strcpy(obscode, "G1");
        strcpy(outname, "GB43m");
    } else if (strcmp(scope, "nancay") == 0) {
        strcpy(obscode, "NC");
        strcpy(outname, "Nancay");
    } else if (strcmp(scope, "effelsberg") == 0) {
        strcpy(obscode, "EF");
        strcpy(outname, "Effelsberg");
    } else if (strcmp(scope, "wsrt") == 0) {
        strcpy(obscode, "WT");
        strcpy(outname, "WSRT");
    } else if (strcmp(scope, "gmrt") == 0) {
        strcpy(obscode, "GM");
        strcpy(outname, "GMRT");
    } else if (strcmp(scope, "lofar") == 0 ) {
      strcpy(obscode, "LF");
      strcpy(outname, "LOFAR");
    } else if (strcmp(scope, "geocenter") == 0) {
        strcpy(obscode, "EC");
        strcpy(outname, "Geocenter");
    } else {
        printf("\nWARNING!!!:  I don't recognize the observatory (%s)!\n",
               inname);
        printf("                 Defaulting to the Geocenter for TEMPO.\n");
        strcpy(obscode, "EC");
        strcpy(outname, "Unknown");
    }
}


float invsqrt(float x)
// See http://en.wikipedia.org/wiki/Fast_inverse_square_root
{
    union {
        float f;
        int i;
    } tmp;
    tmp.f = x;
    tmp.i = 0x5f3759df - (tmp.i >> 1);
    float y = tmp.f;
    return y * (1.5f - 0.5f * x * y * y);
}


long long next2_to_n(long long x)
/* Return the first value of 2^n >= x */
{
   long long i = 1;

   while (i < x)
      i <<= 1;
   return i;
}

#if 0
int gcd(int a, int b)
/* Return the greatest common divisor of a and b */
{
   int aa, bb, tmpa;

   aa = abs(a);
   bb = abs(b);
   while (a) {
      tmpa = a;
      a = b % a;
      b = tmpa;
   }
   return b;
}
#endif

float *gen_freqs(long numfreqs, double lof, double df)
/* This routine generates a float vector of length numfreqs */
/* with values set to lof, lof+df, lof+2df, ...             */
/* It is normally used when generating a list of freqs      */
/* for an x-y plot of spectral data.                        */
{
   long i;
   float *freqs;

   freqs = gen_fvect(numfreqs);
   for (i = 0; i < numfreqs; i++) {
      freqs[i] = lof + i * df;
   }
   return freqs;
}


double *gen_dfreqs(long numfreqs, double lof, double df)
/* This routine generates a double vector of length numfreqs */
/* with values set to lof, lof+df, lof+2df, ...              */
/* It is normally used when generating a list of freqs       */
/* for an x-y plot of spectral data.                         */
{
   long i;
   double *freqs;

   freqs = gen_dvect(numfreqs);
   for (i = 0; i < numfreqs; i++) {
      freqs[i] = lof + i * df;
   }
   return freqs;
}


void i_to_n(int n, double *rl, double *im)
/* Return the real and imaginary portions of i^n */
{
   n %= 4;
   if (n < 0)
      n += 4;
   if (n == 0) {
      *rl = 1.0;
      *im = 0.0;
   } else if (n == 1) {
      *rl = 0.0;
      *im = 1.0;
   } else if (n == 2) {
      *rl = -1.0;
      *im = 0.0;
   } else {
      *rl = 0.0;
      *im = -1.0;
   }
}


void rotate_1d(float *data, long numbins, long bins_to_left)
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of FLOATING points to move.       */
/* frotate is better.  Use it.                        */
{
   float *tmp;

   if (bins_to_left < 0 || bins_to_left >= numbins) {
      printf("\nNumber of bins to rotate array in rotate_1d is\n");
      printf("\nout of bounds.  Tried to rotate %ld bins.  Exiting.\n",
             bins_to_left);
      exit(1);
   }
   tmp = gen_fvect(bins_to_left);
   memcpy(tmp, data, sizeof(float) * bins_to_left);
   memmove(data, data + bins_to_left, sizeof(float) * (numbins - bins_to_left));
   memcpy(data + bins_to_left, tmp, sizeof(float) * bins_to_left);
   vect_free(tmp);
}


void drotate_1d(double *data, long numbins, long bins_to_left)
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */
/* drotate is better.  Use it.                       */
{
   double *tmp;

   if (bins_to_left < 0 || bins_to_left >= numbins) {
      printf("\nNumber of bins to rotate array in rotate_1d is\n");
      printf("\nout of bounds.  Tried to rotate %ld bins.  Exiting.\n",
             bins_to_left);
      exit(1);
   }
   tmp = gen_dvect(bins_to_left);
   memcpy(tmp, data, sizeof(double) * bins_to_left);
   memmove(data, data + bins_to_left, sizeof(double) * (numbins - bins_to_left));
   memcpy(data + bins_to_left, tmp, sizeof(double) * bins_to_left);
   vect_free(tmp);
}


void frotate(float *data, long numbins, float bins_to_left)
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of FLOATING points to move.       */
{
   float *tmp;
   double lopart, hipart, intpart;
   long i, index;

   bins_to_left = fmod(bins_to_left, (double) numbins);
   if (bins_to_left < 0.0)
      bins_to_left += numbins;
   tmp = gen_fvect(numbins);
   lopart = modf(bins_to_left, &intpart);
   hipart = 1.0 - lopart;
   index = (long) floor(intpart + 1.0E-20);
   for (i = 0; i < numbins; i++)
      tmp[i] = hipart * data[(index + i) % numbins] +
          lopart * data[(index + i + 1) % numbins];
   memcpy(data, tmp, sizeof(float) * numbins);
   vect_free(tmp);
}


void drotate(double *data, long numbins, double bins_to_left)
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */
{
   double *tmp, lopart, hipart, intpart;
   long i, index;

   bins_to_left = fmod(bins_to_left, (double) numbins);
   if (bins_to_left < 0.0)
      bins_to_left += numbins;
   tmp = gen_dvect(numbins);
   lopart = modf(bins_to_left, &intpart);
   hipart = 1.0 - lopart;
   index = (long) floor(intpart + 1.0E-20);
   for (i = 0; i < numbins; i++)
      tmp[i] = hipart * data[(index + i) % numbins] +
          lopart * data[(index + i + 1) % numbins];
   memcpy(data, tmp, sizeof(double) * numbins);
   vect_free(tmp);
}


void stats(float *x, int n, double *mean, double *var, double *skew, double *kurt)
/* For a floating point vector, *x, of length n, this routine  */
/* returns the mean, variance, skewness, and kurtosis of *x.   */
{
   long i;
   double an = 0.0, an1 = 0.0, dx, t1, t2, stdevmle;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in stats().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
      *skew = 0.0;
      *kurt = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = ((double) x[i] - *mean) / an;
      t1 = an1 * dx * dx;
      t2 = an * t1;
      *kurt -= dx * (4.0 * *skew - dx * (6.0 * *var + t1 * (1.0 + an1 * an1 * an1)));
      *skew -= dx * (3.0 * *var - t2 * (an - 2.0));
      *var += t2;
      *mean += dx;
   }

   if (n > 1) {
      stdevmle = sqrt(*var / an);
      t1 = stdevmle * stdevmle * stdevmle;
      *skew /= t1 * an;
      *kurt = *kurt / (t1 * stdevmle * an) - 3.0;
      *var /= an1;
   }

   return;
}


void dstats(double *x, int n, double *mean, double *var, double *skew, double *kurt)
/* For a double precision vector, *x, of length n, this routine  */
/* returns the mean, variance, skewness, and kurtosis of *x.     */
{
   long i;
   double an = 0.0, an1 = 0.0, dx, t1, t2, stdevmle;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in dstats().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
      *skew = 0.0;
      *kurt = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      t1 = an1 * dx * dx;
      t2 = an * t1;
      *kurt -= dx * (4.0 * *skew - dx * (6.0 * *var + t1 * (1.0 + an1 * an1 * an1)));
      *skew -= dx * (3.0 * *var - t2 * (an - 2.0));
      *var += t2;
      *mean += dx;
   }

   if (n > 1) {
      stdevmle = sqrt(*var / an);
      t1 = stdevmle * stdevmle * stdevmle;
      *skew /= t1 * an;
      *kurt = *kurt / (t1 * stdevmle * an) - 3.0;
      *var /= an1;
   }

   return;
}


void avg_var(float *x, int n, double *mean, double *var)
/* For a float vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.               */
{
   long i;
   double an = 0.0, an1 = 0.0, dx;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in avg_var().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      *var += an * an1 * dx * dx;
      *mean += dx;
   }

   if (n > 1)
      *var /= an1;

   return;
}


void davg_dvar(double *x, int n, double *mean, double *var)
/* For a double vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.                */
{
   long i;
   double an = 0.0, an1 = 0.0, dx;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in avg_var().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      *var += an * an1 * dx * dx;
      *mean += dx;
   }

   if (n > 1)
      *var /= an1;

   return;
}


inline void update_stats(int N, double x, double *min, double *max,
                         double *avg, double *var)
/* Update time series statistics using one-pass technique */
{
   double dev;

   /* Check the max and min values */

   if (x > *max)
      *max = x;
   if (x < *min)
      *min = x;

   /* Use clever single pass mean and variance calculation */

   dev = x - *avg;
   *avg += dev / (N + 1.0);
   *var += dev * (x - *avg);
}


void ra_dec_to_string(char *radec, int h_or_d, int m, double s)
/* Return a properly formatted string containing RA or DEC values   */
/*   radec is a string with J2000 RA  in the format 'hh:mm:ss.ssss' */
/*   or a string with J2000 DEC in the format 'dd:mm:ss.ssss'       */
{
   int offset=0;

   if (h_or_d==0 && (m < 0 || s < 0.0)) {
       radec[0] = '-';
       offset = 1;
   }
   if (fabs(s) >= 10.0) {
       sprintf(radec+offset, "%.2d:%.2d:%.4f", h_or_d, abs(m), fabs(s));
   } else {
       sprintf(radec+offset, "%.2d:%.2d:0%.4f", h_or_d, abs(m), fabs(s));
   }
}


void ra_dec_from_string(char *radec, int *h_or_d, int *m, double *s)
/* Return a values for hours or degrees, minutes and seconds        */
/* given a properly formatted RA or DEC string.                     */
/*   radec is a string with J2000 RA  in the format 'hh:mm:ss.ssss' */
/*   or a string with J2000 DEC in the format 'dd:mm:ss.ssss'       */
{
   radec = remove_whitespace(radec);
   sscanf(radec, "%d:%d:%lf\n", h_or_d, m, s);
   if (radec[0]=='-' && *h_or_d==0) {
      *m = -*m;
      *s = -*s;
   }
}


double hms2hours(int hours, int min, double sec)
/* Convert hours, minutes, and seconds of time to hours */
{
   return (double) hours + (((double) min + sec / 60.0) / 60.0);
}

void hours2hms(double hours, int *h, int *m, double *s)
/* Convert decimal hours to hours, minutes, and seconds */
{
   double tmp;

   *h = (int) floor(hours);
   tmp = (hours - *h) * 60.0;
   *m = (int) floor(tmp);
   *s = (tmp - *m) * 60.0;
}

void deg2dms(double degrees, int *d, int *m, double *s)
/* Convert decimal degrees to degrees, minutes, and seconds */
{
   int sign = 1;
   double tmp;

   if (degrees < 0.0)
      sign = -1;
   *d = (int) floor(fabs(degrees));
   tmp = (fabs(degrees) - *d) * 60.0;
   *m = (int) floor(tmp);
   *s = (tmp - *m) * 60.0;
   *d *= sign;
   if (*d==0) {
       *m *= sign;
       *s *= sign;
   }
}

double dms2rad(int deg, int min, double sec)
/* Convert degrees, minutes, and seconds of arc to radians */
{
   double sign = 1.0;

   if (deg < 0) sign = -1.0;
   if (deg==0 && (min < 0 || sec < 0.0)) sign = -1.0;
   return sign * ARCSEC2RAD * (60.0 * (60.0 * (double) abs(deg)
                                       + (double) abs(min)) + fabs(sec));
}

double hms2rad(int hour, int min, double sec)
/* Convert hours, minutes, and seconds of arc to radians */
{
   return SEC2RAD * (60.0 * (60.0 * (double) hour + (double) min) + sec);
}

double mjd_sec_diff(int int1, double frac1, int int2, double frac2)
/* Return the difference in seconds between two MJDs (1 - 2) */
{
   int idiff;
   double fdiff;

   idiff = (int1 - int2) * 86400;
   fdiff = (frac1 - frac2) * 86400.0;
   return (double) (idiff + fdiff);
}

double sphere_ang_diff(double ra1, double dec1, double ra2, double dec2)
/* Returns the angular difference in radians between two sets */
/* of RA and DEC (in radians).                                */
{

   int i;
   double d, vec1[3], vec2[3], s2, c2, cosb;

   /* Convert coordinates from spherical to Cartesian */

   cosb = cos(dec1);
   vec1[0] = cos(ra1) * cosb;
   vec1[1] = sin(ra1) * cosb;
   vec1[2] = sin(dec1);
   cosb = cos(dec2);
   vec2[0] = cos(ra2) * cosb;
   vec2[1] = sin(ra2) * cosb;
   vec2[2] = sin(dec2);

   /* Modulus squared of half the difference vector */

   s2 = 0.0;
   for (i = 0; i < 3; i++) {
      d = vec1[i] - vec2[i];
      s2 += d * d;
   }
   s2 /= 4.0;

   /* Angle between the vectors */

   c2 = 1.0 - s2;
   return 2.0 * atan2(sqrt(s2), sqrt(DMAX(0.0, c2)));
}

#undef DMAX
