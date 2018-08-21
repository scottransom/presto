#include <ctype.h>
#include "misc_utils.h"
#include "slamac.h"

#ifndef ARCSEC2RAD
/* arcseconds to radians */
#define ARCSEC2RAD 4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
/* seconds of time to radians */
#define SEC2RAD    7.2722052166430399038487115353692196393452995355905e-5
#endif
/* Speef of light in m/s */
#ifndef SOL
#define SOL           299792458.0
#endif
/* Radians to degrees */
#ifndef RADTODEG
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#endif


/* return the max of two double values */
#define DMAX(a,b) ((a)>(b)?(a):(b))

/* like strcpy, except guaranteed to work with overlapping strings      */
#define strMove(d,s) memmove(d,s,strlen(s)+1)

void slaDjcl(double djm, int *iy, int *im, int *id, double *fd, int *j);

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
    char *sptr = NULL;
    unsigned int len, pathlen = 0, filelen = 0;

    len = strlen(input);
    sptr = strrchr(input, '/');
    // The 2nd part of the following handles relative paths,
    // but in a strange way, in that the file will have the
    // relative part in it.  The path will be the absolute
    // path of the current working directory.  That's usually
    // euivalent to what we want, though.
    if ((sptr == NULL) || (input[0] == '.')) {
        sptr = getcwd(NULL, 0);
        if (sptr == NULL) {
            perror
                ("Error:  could not get current directory name (too long?) in split_path_file()\n");
            exit(-1);
        }
        pathlen = strlen(sptr);
        *path = (char *) calloc(pathlen + 1, sizeof(char));
        *file = (char *) calloc(len + 1, sizeof(char));
        strcpy(*path, sptr);
        strncpy(*file, input, len);
        free(sptr);
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


void mjd_to_datestr(double mjd, char *datestr)
// Convert an MJD to a PSRFITS-style DATEOBS
{
    int year, month, day, hour, min, err;
    double fracday, dh, dm, sec;

    slaDjcl(mjd, &year, &month, &day, &fracday, &err);
    if (err == -1) {
        printf("Error in mjd_to_datestr:  Bad MJD '%.12f'.\n", mjd);
        exit(1);
    }
    dh = fracday * 24.0;
    hour = (int) dh;
    dm = (dh - hour) * 60.0;
    min = (int) dm;
    sec = (dm - min) * 60.0;
    if (sec < 10.0) {
        sprintf(datestr, "%4d-%02d-%02dT%02d:%02d:0%.6g",
                year, month, day, hour, min, sec);
    } else {
        sprintf(datestr, "%4d-%02d-%02dT%02d:%02d:%.6g",
                year, month, day, hour, min, sec);
    }
}


void telescope_to_tempocode(char *inname, char *outname, char *obscode)
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
               (strcmp(scope, "gb 140ft") == 0) || (strcmp(scope, "nrao20") == 0)) {
        strcpy(obscode, "G1");
        strcpy(outname, "GB43m");
    } else if (strcmp(scope, "nancay") == 0) {
        strcpy(obscode, "NC");
        strcpy(outname, "Nancay");
    } else if (strcmp(scope, "effelsberg") == 0) {
        strcpy(obscode, "EF");
        strcpy(outname, "Effelsberg");
    } else if (strcmp(scope, "srt") == 0) {
        strcpy(obscode, "SR");
        strcpy(outname, "Sardinia Radio Telescope");
    } else if (strcmp(scope, "fast") == 0) {
        strcpy(obscode, "FA");
        strcpy(outname, "FAST");
    } else if (strcmp(scope, "wsrt") == 0) {
        strcpy(obscode, "WT");
        strcpy(outname, "WSRT");
    } else if (strcmp(scope, "gmrt") == 0) {
        strcpy(obscode, "GM");
        strcpy(outname, "GMRT");
    } else if (strcmp(scope, "chime") == 0) {
        strcpy(obscode, "CH");
        strcpy(outname, "CHIME");
    } else if (strcmp(scope, "lofar") == 0) {
        strcpy(obscode, "LF");
        strcpy(outname, "LOFAR");
    } else if (strcmp(scope, "lwa") == 0) {
        strcpy(obscode, "LW");
        strcpy(outname, "LWA1");
    } else if (strcmp(scope, "mwa") == 0 ) {
        strcpy(obscode, "MW");
        strcpy(outname, "MWA128T");
    } else if (strcmp(scope, "meerkat") == 0 ) {
        strcpy(obscode, "MK");
        strcpy(outname, "MeerKAT");
    } else if (strcmp(scope, "ata") == 0) {
        strcpy(obscode, "AT");
        strcpy(outname, "ATA");
    } else if (strcmp(scope, "k7") == 0 ) {
        strcpy(obscode, "K7");
        strcpy(outname, "KAT-7");
    } else if (strcmp(scope, "geocenter") == 0) {
        strcpy(obscode, "0 ");
        strcpy(outname, "Geocenter");
    } else {
        printf("\nWARNING!!!:  I don't recognize the observatory (%s)!\n", inname);
        printf("                 Defaulting to the Geocenter for TEMPO.\n");
        strcpy(obscode, "0 ");
        strcpy(outname, "Unknown");
    }
}


float invsqrtf(float x)
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


int is_power_of_10(long long n)
/* Check whether n is a power of 10 or not.  Return 0 or 1 */
{
    while (n > 9L && n % 10L == 0)
        n /= 10L;
    return n == 1;
}


long long choose_good_N(long long orig_N)
// Choose a time series length that is larger than the input value but
// that is highly factorable.
{
    char ctmp20[20], ctmp5[5];
    long long two_N = 2, small_N;
    int ii, first4;
    int goodfactors[114] = {1000, 1008, 1024, 1056, 1120, 1152, 1200, 1232,
                            1280, 1296, 1344, 1408, 1440, 1536, 1568, 1584,
                            1600, 1680, 1728, 1760, 1792, 1920, 1936, 2000,
                            2016, 2048, 2112, 2160, 2240, 2304, 2352, 2400,
                            2464, 2560, 2592, 2640, 2688, 2800, 2816, 2880,
                            3024, 3072, 3136, 3168, 3200, 3360, 3456, 3520,
                            3584, 3600, 3696, 3840, 3872, 3888, 3920, 4000,
                            4032, 4096, 4224, 4320, 4400, 4480, 4608, 4704,
                            4752, 4800, 4928, 5040, 5120, 5184, 5280, 5376,
                            5488, 5600, 5632, 5760, 5808, 6000, 6048, 6144,
                            6160, 6272, 6336, 6400, 6480, 6720, 6912, 7040,
                            7056, 7168, 7200, 7392, 7680, 7744, 7776, 7840,
                            7920, 8000, 8064, 8192, 8400, 8448, 8624, 8640,
                            8800, 8960, 9072, 9216, 9408, 9504, 9600, 9680,
                            9856, 10000};
    if (orig_N <= 0) return 0;
    // Get the number represented by the first 4 digits of orig_N
    sprintf(ctmp20, "%lld", orig_N);
    first4 = atoi(strncpy(ctmp5, ctmp20, 4));
    // Now get the number that is just bigger than orig_N
    // that has its first 4 digits equal to "factor"
    for (ii = 0; ii < 114; ii++) {
        small_N = goodfactors[ii];
        //  check to see if orig_N is a goodfactor times a power of 10
        if (small_N == first4 &&
            orig_N % small_N == 0 &&
            is_power_of_10(orig_N/small_N)) break;
        if (small_N > first4) break;
    }
    while (small_N < orig_N) small_N *= 10;
    // Finally, compare new_N to the closest power_of_two
    // greater than orig_N.  Take the closest.
    while (two_N < orig_N) two_N *= 2;
    if (two_N < small_N) return two_N;
    else return small_N;
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


float beam_halfwidth(float freq, float dish_diam)
// Return the beam halfwidth in arcsec when freq
// is in MHz and dish_diam is in meters
{
    return 0.5 * 1.2 * SOL / (freq * 1e6) / dish_diam * RADTODEG * 3600.0;
}


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
        *kurt -=
            dx * (4.0 * *skew - dx * (6.0 * *var + t1 * (1.0 + an1 * an1 * an1)));
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
        *kurt -=
            dx * (4.0 * *skew - dx * (6.0 * *var + t1 * (1.0 + an1 * an1 * an1)));
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


void ra_dec_to_string(char *radec, int h_or_d, int m, double s)
/* Return a properly formatted string containing RA or DEC values   */
/*   radec is a string with J2000 RA  in the format 'hh:mm:ss.ssss' */
/*   or a string with J2000 DEC in the format 'dd:mm:ss.ssss'       */
{
    int offset = 0;

    if (h_or_d == 0 && (m < 0 || s < 0.0)) {
        radec[0] = '-';
        offset = 1;
    }
    sprintf(radec + offset, "%.2d:%.2d:%07.4f", h_or_d, abs(m), fabs(s));
}


void ra_dec_from_string(char *radec, int *h_or_d, int *m, double *s)
/* Return a values for hours or degrees, minutes and seconds        */
/* given a properly formatted RA or DEC string.                     */
/*   radec is a string with J2000 RA  in the format 'hh:mm:ss.ssss' */
/*   or a string with J2000 DEC in the format 'dd:mm:ss.ssss'       */
{
    int retval;

    radec = remove_whitespace(radec);
    retval = sscanf(radec, "%d:%d:%lf\n", h_or_d, m, s);
    if (retval != 3) {
        char tmp[100];
        sprintf(tmp,
                "Error:  can not convert '%s' to RA or DEC in ra_dec_from_string()\n",
                radec);
        perror(tmp);
        exit(1);
    }
    if (radec[0] == '-' && *h_or_d == 0) {
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
    if (*d == 0) {
        *m *= sign;
        *s *= sign;
    }
}

double dms2rad(int deg, int min, double sec)
/* Convert degrees, minutes, and seconds of arc to radians */
{
    double sign = 1.0;

    if (deg < 0)
        sign = -1.0;
    if (deg == 0 && (min < 0 || sec < 0.0))
        sign = -1.0;
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
