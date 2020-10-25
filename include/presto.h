#include <sys/times.h>
#include "chkio.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include "clk_tck.h"
#include "ransomfft.h"
#include "meminfo.h"
#include "vectors.h"
#include "makeinf.h"
#include "orbint.h"
#include "misc_utils.h"
#include "cdflib.h"
#include "database.h"
#include "makedata.h"

#ifndef SQRT2
#define SQRT2         1.4142135623730950488016887242096980785696718753769
#endif
#ifndef DBLCORRECT
#define DBLCORRECT    1e-14
#endif
#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#endif
#ifndef PIBYTWO
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#endif
#ifndef SOL
#define SOL           299792458.0
#endif
#ifndef SECPERJULYR
#define SECPERJULYR   31557600.0
#endif
#ifndef SECPERDAY
#define SECPERDAY     86400.0
#endif
#ifndef ARCSEC2RAD
#define ARCSEC2RAD    4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
#define SEC2RAD       7.2722052166430399038487115353692196393452995355905e-5
#endif
#ifndef __GNUC__
#define __inline__
#endif
/* Maximum number of input files to try and patch together */
#define MAXPATCHFILES 100
/* Blocksize to use when reading datafiles or subbands */
#define SUBSBLOCKLEN 1024

/* various function-like macros */

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

#ifndef POWER
/* Returns unnormalized Fourier power  */
/*   Requires the following variables in calling function */
/*   double powargr, powargi; */
#define POWER(r,i) (powargr=(r),powargi=(i),\
		    powargr*powargr+powargi*powargi)
#endif

#ifndef PHASE
/* Returns Fourier phase (degrees)  */
/*   Requires the following variables in calling function */
/*   double phsargr, phsargi, phstmp; */
#define PHASE(r,i) (phsargr=(r),phsargi=(i),\
		    ((phstmp=RADTODEG*atan2(phsargi,phsargr)) > 0.0) ? \
		    phstmp : phstmp+360.0)
#endif
  
#ifndef RADIAN_PHASE
/* Returns Fourier phase (radians)  */
/*   Requires the following variables in calling function */
/*   double radargr, radargi, radtmp; */
#define RADIAN_PHASE(r,i) (radargr=(r),radargi=(i),\
		    ((radtmp=atan2(radargi,radargr)) > 0.0) ? \
		    radtmp : radtmp+TWOPI)
#endif

#define GET_BIT(c, n) (*(c+(n>>3)) >> (7-(n&7)) & 1)
#define SET_BIT(c, n) (*(c+(n>>3)) |= 1 << (7-(n&7)))
#define UNSET_BIT(c, n) (*(c+(n>>3)) &= ~(1 << (7-(n&7))))

/*   Number of bins (total) to average for local power:   */
/*     Must be an even number (1/2 on each side).         */
#define NUMLOCPOWAVG  20

/*  Number of bins next to freq in question to            */
/*    ignore (on each side) when determining local power. */
#define DELTAAVGBINS  5

/*  Number of bins on each side of the central frequency  */
/*    to sum for Fourier interpolation (low accuracy)     */
#define NUMFINTBINS   16

/* Used for raw-data handling */
typedef enum {
  IF0, IF1, SUMIFS
} IFs;

/*  Constants used in the interpolation routines */
typedef enum {
  LOWACC, HIGHACC
} presto_interp_acc;

/*  Constants used in the binary search routines */
typedef enum {
  INTERBIN, INTERPOLATE
} presto_interptype;

typedef enum {
  NO_CHECK_ALIASED, CHECK_ALIASED
} presto_checkaliased;

/*  Constants used in the correlation/convolution routines */
typedef enum {
  CONV, CORR, INPLACE_CONV, INPLACE_CORR
} presto_optype;

typedef enum {
  FFTDK, FFTD, FFTK, NOFFTS
} presto_ffts;

typedef enum {
  RAW, PREPPED, FFT, SAME
} presto_datainf;

typedef struct FFTCAND {
  int nsum;   /* Number of harmonics summed */
  double p;   /* Summed normalized power */
  double r;   /* Frequency of lowest harmonic */
  double sig; /* Significance of the candidate taking N into account */
} fftcand;

typedef struct POSITION {
    float pow;	    /* Power normalized with local power             */
    double p1;	    /* r if rzw search, r_startfft if bin search     */
    double p2;	    /* z if rzw search, r_freqmod if bin search      */
    double p3;	    /* w if rzw search, numfftbins if bin search     */
} position;

typedef struct POWINDEX {
    float pow;	    /* Power of Fourier bin                          */
    long ind;	    /* Array index of that Fourier bin               */
} powindex;

typedef struct RDERIVS {
    double pow;	     /* Power normalized with local power   */
    double phs;	     /* Signal phase                        */
    double dpow;     /* 1st deriv of power wrt fourier freq */
    double dphs;     /* 1st deriv of phase wrt fourier freq */
    double d2pow;    /* 2nd deriv of power wrt fourier freq */
    double d2phs;    /* 2nd deriv of power wrt fourier freq */
    double locpow;   /* Local mean power level              */
} rderivs;

typedef struct FOURIERPROPS {
    double r;	    /* Fourier freq at start of run                  */
    float rerr;	    /* Error in fourier freq                         */
    double z;	    /* Fourier freq deriv at start of run = T^2*fdot */
    float zerr;	    /* Error in fourier freq deriv                   */
    double w;	    /* Fourier 2nd freq deriv = T^3*fdotdot          */
    float werr;	    /* Error in 2nd fourier freq deriv               */
    float pow;	    /* Power normalized with local power             */
    float powerr;   /* Error in power                                */
    float sig;	    /* Approx. Sigma level                           */
    float rawpow;   /* Raw signal power (unnormalized)               */
    float phs;	    /* Signal phase                                  */
    float phserr;   /* Error in phase                                */
    float cen;	    /* Centroid of signal pulsations over run        */
    float cenerr;   /* Error in centroid                             */
    float pur;	    /* Purity:  rms duration of signal over run      */
    float purerr;   /* Error in purity                               */
    float locpow;   /* Local mean power level                        */
} fourierprops;

typedef struct BINARYPROPS {
    double ppsr;       /* Pulsar period (s)                            */
    double fpsr;       /* Pulsar freq (hz)                             */
    double rpsr;       /* Pulsar Fourier freq (bins)                   */
    double pbin;       /* Binary period (s)                            */
    double rbin;       /* Binary Fourier freq (bins)                   */
    double z;	       /* Amplitude of frequency modulation            */
    double asinic;     /* Projected semi-major axis of orbit (lt-sec)  */
    double rdetect;    /* Mini-FFT candidate detected freq (bins)      */
    long   nfftbins;   /* Number of FFT bins in mini-fft               */
    long   lowbin;     /* Start bin of mini-FFT in original FFT        */
    float  ppsrerr;    /* Error in pulsar period (s)                   */
    float  fpsrerr;    /* Error in pulsar freq (hz)                    */
    float  rpsrerr;    /* Error in pulsar Fourier freq (bins)          */
    float  pbinerr;    /* Error in binary period (s)                   */
    float  rbinerr;    /* Error in pulsar Fourier freq (bins)          */
    float  zerr;       /* Error in freq modulation ampltude            */
    float  asinicerr;  /* Error in a*sin(i)/c in (lt-sec)              */
    float  rdetecterr; /* Error in Mini-FFT detected freq (bins)       */
    float  sig;	       /* Approx. Sigma level                          */
    float  phs;	       /* Signal phase                                 */
    float  phserr;     /* Error in phase                               */
    float  cen;	       /* Centroid of signal pulsations over run       */
    float  cenerr;     /* Error in centroid                            */
    float  pur;	       /* Purity:  rms duration of signal over run     */
    float  purerr;     /* Error in purity                              */
    float  pow;        /* Power normalized with Nph equivalent         */
    float  powerr;     /* Error in signal power                        */
} binaryprops;

typedef struct RAWBINCAND {
  double full_N;       /* Number of points in original time series  */
  double full_T;       /* Length (s) of original time series        */
  double full_lo_r;    /* Lowest Fourier bin that was miniFFTd      */
  double mini_N;       /* Number of points in short (mini) FFT      */
  double mini_r;       /* Candidate Fourier bin in miniFFT          */
  double mini_power;   /* Candidate normalized power in miniFFT     */
  double mini_numsum;  /* Number of powers summed to get candidate  */
  double mini_sigma;   /* Equivalent candidate sigma (for sum pow)  */
  double psr_p;        /* Approx PSR period (miniFFT center bin)    */
  double orb_p;        /* Approx orbital period (s)                 */
} rawbincand;

typedef struct RZWERRS {
    double p;	    /* Pulsation period at start of run (s)          */
    double pd;	    /* Period derivative at start of run (s/s)       */
    double pdd;	    /* Period 2nd derivative (s/s^2)                 */
    double f;	    /* Pulsation frequency at start of run (hz)      */
    double fd;	    /* Frequency derivative at start of run (hz/s)   */
    double fdd;	    /* Frequency 2nd derivative (hz/s^2)             */
    float perr;	    /* Error in pulsation period (s)                 */
    float pderr;    /* Error in period derivative (s/s)              */
    float pdderr;   /* Error in 2nd period derivative (s/s^2)        */
    float ferr;	    /* Error in frequency (hz)                       */
    float fderr;    /* Error in frequency derivative (hz/s)          */
    float fdderr;   /* Error in 2nd frequency derivative (hh/s^2)    */
} rzwerrs;

#ifndef ORBITPARAMS_TYPE
typedef struct orbitparams {
    double p;	    /* Orbital period (s)                            */
    double e;	    /* Orbital eccentricity                          */
    double x;	    /* Projected semi-major axis (lt-sec)            */
    double w;	    /* Longitude of periapsis (deg)                  */
    double t;	    /* Time since last periastron passage (s)        */
    double pd;	    /* Orbital period derivative (s/yr)              */
    double wd;	    /* Advance of longitude of periapsis (deg/yr)    */
} orbitparams;
#define ORBITPARAMS_TYPE 1
#endif

typedef struct foldstats {
  double numdata;     /* Number of data bins folded         */
  double data_avg;    /* Average level of the data bins     */
  double data_var;    /* Variance of the data bins          */
  double numprof;     /* Number of bins in the profile      */
  double prof_avg;    /* Average level of the profile bins  */
  double prof_var;    /* Variance of the profile bins       */
  double redchi;      /* Reduced chi-squared of the profile */
} foldstats;

typedef struct bird{
  double lobin;
  double hibin;
} bird;

/*****  Function Prototypes    *****/

/* From swapendian.c: */

float swap_float(float var);
double swap_double(double var);
long double swap_longdouble(long double var);
long long swap_longlong(long long var);
int swap_int(int var);
unsigned int swap_uint(unsigned int var);
short swap_short(short var);
unsigned short swap_ushort(unsigned short var);

/* From responses.c: */

int r_resp_halfwidth(presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a standard Fourier signal.                  */
  /*  Arguments:                                                       */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

int z_resp_halfwidth(double z, presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
  /*  a constant frequency derivative)                                 */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

int w_resp_halfwidth(double z, double w, presto_interp_acc accuracy);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with an f-dot that (i.e    */
  /*  varies linearly in time -- a constant f-dotdot)                  */
  /*  Arguments:                                                       */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */

void binary_velocity(double T, orbitparams * orbit, 
		     double *minv, double *maxv);
  /*  Return the minimum and maximum orbital velocities of a pulsar    */
  /*  during an observation as a fraction of the speed of light.       */
  /*  Arguments:                                                       */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */

int bin_resp_halfwidth(double ppsr, double T, orbitparams * orbit);
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a pulsar in a binary orbit.                 */
  /*  Arguments:                                                       */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2 * 'numbetween' to get the   */
  /*    length of the array required to hold such a kernel.            */

fcomplex *gen_r_response(double roffset, int numbetween, int numkern);
  /*  Generate a complex response function for Fourier interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

fcomplex *gen_z_response(double roffset, int numbetween, double z, \
			 int numkern);
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

fcomplex *gen_w_response(double roffset, int numbetween, double z, \
                         double w, int numkern);
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

fcomplex *gen_w_response2(double roffset, int numbetween, double z, \
                          double w, int numkern);
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */

fcomplex *gen_bin_response(double roffset, int numbetween, double ppsr, \
			   double T, orbitparams * orbit, int numkern);
  /*  Generate the Fourier response function for a sinusoidal PSR      */
  /*  signal from a binary orbit.                                      */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 1 = interbins) */
  /*    'ppsr' is the period of the pusar in seconds.                  */
  /*    'T' is the length of the observation in seconds.               */
  /*    'orbit' is a ptr to a orbitparams structure containing the     */
  /*       Keplerian orbital parameters of the binary system.          */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */


/* From characteristics.c:    */
/* Routines to calculate characteristics of a candidate signal */


float get_numphotons(FILE * file);
  /* Return the total number of photons in the FFT file      */
  /* i.e.  it returns the value of the 0th frequency bin.    */
  /* Arguments:                                              */
  /*   'file' is a pointer to the file you want to access.   */

double get_localpower(fcomplex *data, long numdata, double r);
  /* Return the local power level at specific FFT frequency.  */
  /* Arguments:                                               */
  /*   'data' is a pointer to a complex FFT.                  */
  /*   'numdata' is the number of complex points in 'data'.   */
  /*   'r' is the Fourier frequency in data that we want to   */
  /*      interpolate.                                        */

double get_localpower3d(fcomplex *data, long numdata, double r, \
			double z, double w);
  /* Return the local power level around a specific FFT           */
  /* frequency, f-dot, and f-dotdot.                              */
  /* Arguments:                                                   */
  /*   'data' is a pointer to a complex FFT.                      */
  /*   'numdata' is the number of complex points in 'data'.       */
  /*   'r' is the Fourier frequency in data that we want to       */
  /*      interpolate.                                            */
  /*   'z' is the Fourier Frequency derivative (# of bins the     */
  /*       signal smears over during the observation).            */
  /*   'w' is the Fourier Frequency 2nd derivative (change in the */
  /*       Fourier f-dot during the observation).                 */

void get_derivs3d(fcomplex *data, long numdata, double r, \
		  double z, double w, double localpower, \
		  rderivs *result);
  /* Return an rderives structure that contains the power,      */
  /* phase, and their first and second derivatives at a point   */
  /* in the F/F-dot/F-dortdot volume.                           */  
  /* Arguments:                                                 */
  /*   'data' is a pointer to a complex FFT.                    */
  /*   'numdata' is the number of complex points in 'data'.     */
  /*   'r' is the Fourier frequency in data that we want to     */
  /*      interpolate.                                          */
  /*   'z' is the Fourier Frequency derivative (# of bins the   */
  /*       signal smears over during the observation).          */
  /*   'w' is the Fourier Frequency 2nd derivative (change in   */
  /*       the Fourier f-dot during the observation).           */
  /*   'localpower' is the local power level around the signal. */
  /*   'result' is a pointer to an rderivs structure that will  */
  /*       contain the results.                                 */

void calc_props(rderivs data, double r, double z, double w, \
		fourierprops * result);
  /* Return a fourierprops structure that contains the various  */
  /* properties of a signal described by Middleditch, Deich,    */ 
  /* and Kulkarni in _Isolated_Pulsars_, 1993, p372.            */  
  /* Arguments:                                                 */
  /*   'data' is a pointer to an rderivs structure containing   */
  /*       derivative information about the peak in question.   */
  /*   'r' is the Fourier frequency in data that we want to     */
  /*      interpolate.                                          */
  /*   'z' is the Fourier Frequency derivative (# of bins the   */
  /*       signal smears over during the observation).          */
  /*   'w' is the Fourier Frequency second derivative.          */
  /*   'result' is a pointer to an fourierprops structure that  */
  /*       will contain the results.                            */

void calc_binprops(fourierprops * props, double T, int lowbin, \
		   int nfftbins, binaryprops * result);
  /* Return a binaryprops structure that contains the various     */
  /* estimates of the binary pulsar system from a mini-FFT.       */
  /* Arguments:                                                   */
  /*   'props' is a pointer to the candidate's fourierprops.      */
  /*   'T' is the total length (sec) of the original time series. */
  /*   'lowbin' is the Fourier bin number from the original FFT   */
  /*      the lowest bin in the mini-FFT.                         */
  /*   'nfftbins' is the number of bins in the mini-FFT.          */
  /*   'absnorm' is the value of the power normalization          */
  /*      constant for this mini-FFT.                             */
  /*   'result' is the returned binaryprops structure.            */

void calc_rzwerrs(fourierprops *props, double T, rzwerrs *result);
  /* Calculate periods, frequencies, their derivatives        */
  /* and their errors.                                        */
  /* Arguments:                                               */
  /*   'props' is a pointer to a fourierprops structure.      */
  /*   'T' is the length of the data set in sec (i.e. N*dt).  */
  /*   'result' is a pointer to the returned rzwerrs struct.  */

double equivalent_gaussian_sigma(double logp);
/* Return the approximate significance in Gaussian sigmas */
/* corresponding to a natural log probability logp        */

double chi2_logp(double chi2, int dof);
/* Return the natural log probability corresponding to a chi^2 value */
/* of chi2 given dof degrees of freedom. */

double chi2_sigma(double chi2, int dof);
/* Return the approximate significance in Gaussian sigmas        */
/* sigmas of a chi^2 value of chi2 given dof degrees of freedom. */

double candidate_sigma(double power, int numsum, double numtrials);
/* Return the approximate significance in Gaussian       */
/* sigmas of a candidate of numsum summed powers,        */
/* taking into account the number of independent trials. */

double power_for_sigma(double sigma, int numsum, double numtrials);
/* Return the approximate summed power level required */
/* to get a Gaussian significance of 'sigma', taking  */
/* into account the number of independent trials.     */

double chisqr(double *data, int numdata, double avg, double var);
/* Calculates the chi-square of the 'data' which has average */
/* 'avg', and variance 'var'.                                */

void switch_f_and_p(double in, double ind, double indd,
		    double *out, double *outd, double *outdd);
/* Convert p, p-dot, and p-dotdot into f, f-dot, */
/* and f-dotdot or vise-versa.                   */


/*  dispersion.c:  */
/*  Functions to de-disperse data */


double tree_max_dm(int numchan, double dt, double lofreq, double hifreq);
/* Return the maximum Dispersion Measure (dm) in cm-3 pc, the  */
/* tree de-dispersion technique can correct for given a sample */
/* interval 'dt', the number of channels 'numchan', and the    */
/* low and high observation frequencies in MHz.                */

double smearing_from_bw(double dm, double center_freq, double bandwidth);
/* Return the smearing in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, the central frequency */
/* and the bandwith of the observation in MHz.                 */

double delay_from_dm(double dm, double freq_emitted);
/* Return the delay in seconds caused by dispersion, given  */
/* a Dispersion Measure (dm) in cm-3 pc, and the emitted    */
/* frequency (freq_emitted) of the pulsar in MHz.           */

double dm_from_delay(double delay, double freq_emitted);
/* Return the Dispersion Measure in cm-3 pc, that would     */
/* cause a pulse emitted at frequency 'freq_emitted' to be  */
/* delayed by 'delay' seconds.                              */

double *dedisp_delays(int numchan, double dm, double lofreq, 
		      double chanwidth, double voverc);
/* Return an array of delays (sec) for dedispersing 'numchan'    */
/* channels at a DM of 'dm'.  'lofreq' is the center frequency   */
/* in MHz of the lowest frequency channel.  'chanwidth' is the   */
/* width in MHz of each channel.  'voverc' is the observatory's  */
/* velocity towards or away from the source.  This is to adjust  */
/* the frequencies for doppler effects (for no correction use    */
/* voverc=0).  The returned array is allocated by this routine.  */

void dedisp(unsigned char *data, unsigned char *lastdata, int numpts,
	    int numchan, double *dispdelays, float *result);
/* De-disperse a stretch of data with numpts * numchan points. */
/* The delays (in bins) are in dispdelays for each channel.    */
/* The result is returned in result.  The input data and       */
/* dispdelays are always in ascending frequency order.         */
/* Input data are ordered in time, with the channels stored    */
/* together at each time point.                                */ 

double *subband_delays(int numchan, int numsubbands, double dm, 
		       double lofreq, double chanwidth, 
		       double voverc);
/* Return an array of delays (sec) for the highest frequency  */
/* channels of each subband used in a subband de-dispersion.  */
/* These are the delays described in the 'Note:' in the       */
/* description of subband_search_delays().  See the comments  */
/* for dedisp_delays() for more info.                         */

double *subband_search_delays(int numchan, int numsubbands, double dm, 
			      double lofreq, double chanwidth, 
			      double voverc);
/* Return an array of delays (sec) for a subband DM search.  The      */
/* delays are calculated normally for each of the 'numchan' channels  */
/* using the appropriate frequencies at the 'dm'.  Then the delay     */
/* from the highest frequency channel of each of the 'numsubbands'    */
/* subbands is subtracted from each subband.  This gives the subbands */
/* the correct delays for each freq in the subband, but the subbands  */
/* themselves are offset as if they had not been de-dispersed.  This  */
/* way, we can call float_dedisp() on the subbands if needed.         */
/* 'lofreq' is the center frequency in MHz of the lowest frequency    */
/* channel.  'chanwidth' is the width in MHz of each channel.  The    */
/* returned array is allocated by this routine.  'voverc' is used to  */
/* correct the input frequencies for doppler effects.  See the        */
/* comments in dedisp_delays() for more info.                         */
/* Note:  When performing a subband search, the delays for each       */
/*   subband must be calculated with the frequency of the highest     */
/*   channel in each subband, _not_ the center subband frequency.     */

void dedisp_subbands(float *data, float *lastdata,
                     int numpts, int numchan, 
                     int *delays, int numsubbands, float *result);
// De-disperse a stretch of data with numpts * numchan points into
// numsubbands subbands.  Each time point for each subband is a float
// in the result array.  The result array order is subbands of
// increasing frequency together at each time pt.  The delays (in
// bins) are in delays for each channel.  The input data and
// dispdelays are always in ascending frequency order.  Input data are
// ordered in time, with the channels stored together at each time
// point.

void float_dedisp(float *data, float *lastdata,
                  int numpts, int numchan,
                  int *delays, float approx_mean, float *result);
// De-disperse a stretch of data with numpts * numchan points. The
// delays (in bins) are in delays for each channel.  The result is
// returned in result.  The input data and delays are always in
// ascending frequency order.  Input data are ordered in time, with
// the channels stored together at each time point.

void combine_subbands(double *inprofs, foldstats *stats, 
		      int numparts, int numsubbands, int proflen, 
		      int *delays, double *outprofs, 
		      foldstats *outprofstats);
/* Combine 'nparts' sets of 'numsubbands' profiles, each of length     */
/* 'proflen' into a 'nparts' de-dispersed profiles.  The de-dispersion */
/* uses the 'delays' (of which there are 'numsubbands' many) to        */
/* show how many bins to shift each profile to the right.  Only        */
/* positive numbers may be used (left shifts may be accomplished using */
/* the shift modulo 'proflen').  The 'stats' about the profiles are    */
/* combined as well and the combined stats are returned in             */
/* 'outprofstats'. All arrays must be pre-allocated.                   */


/*  output.c:  */
/*  Functions for text-based ouput of information  */


int nice_output_1(char *output, double val, double err, int len);
/* Generates a string in "output" of length len with "val" rounded  */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 1 decimal place.     */
/* Note:  len should be ~ 20 to show full double precision          */
/*   if the base 10 exponent of the error needs to be shown.        */
/*   If len == 0, left-justified minimum length string is returned. */
/*   If len > 0, the string returned has is right justified.        */

int nice_output_2(char *output, double val, double err, int len);
/* Generates a string in "output" of length len with "val" rounded  */
/*   to the appropriate decimal place and the error in parenthesis  */
/*   as in scientific journals.  The error has 2 decimal places.    */
/* Note:  len should be ~ 20 to show full double precision          */
/*   if the base 10 exponent of the error needs to be shown.        */
/*   If len == 0, left-justified minimum length string is returned. */
/*   If len > 0, the string returned has is right justified.        */

void print_candidate(fourierprops * cand, double dt, long N, \
		     double nph, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

void print_bin_candidate(binaryprops * cand, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

void file_reg_candidates(fourierprops cand[], char *notes, int numcands, \
			 double dt, long N, double nph, 
			 char name[], char longname[]);
/* Outputs a .ps file describing all the candidates from a search.  */

void file_bin_candidates(binaryprops cand[], char *notes, \
			 int numcands, char name[]);
/* Outputs a .ps file describing all the binary candidates from a    */
/*   binary search.*/


/*  get_candidates.c:  */
/*  Functions for manipulating candidate files */

int read_rzw_cand(FILE *file, fourierprops *cands);
/* Read the next rzw candidate from the file */
/* If successful, return 1, else 0           */

int read_bin_cand(FILE *file, binaryprops *cands);
/* Read the next binary candidate from the file */
/* If successful, return 1, else 0              */

int read_rawbin_cand(FILE *file, rawbincand *cands);
/* Read the next rawbin candidate from the file */
/* If successful, return 1, else 0              */

void get_rzw_cand(char *filenm, int candnum, fourierprops * cand);
/*  Read the rzw candidate file 'filenm' and return a        */
/*  pointer to the fourierprops that describes it.           */

void get_bin_cand(char *filenm, int candnum, binaryprops * cand);
/*  Read the bin candidate file 'filenm' and return a        */
/*  pointer to the binaryprops that describes it.            */

void get_rawbin_cand(char *filenm, int candnum, rawbincand * cand);
/*  Read the rawbin candidate file 'filenm' and return a     */
/*  pointer to the rawbincand that describe it.              */

/*  read_fft.c:  */
/*  Functions for getting information from an FFT file  */

fcomplex *read_fcomplex_file(FILE *file, long firstpt, long numpts);
/* Return an fcomplex vector with complex data taken from a file. */
/* Argumants:                                                     */
/*   'file' is a pointer to the file you want to access.          */
/*   'firstpt' is the number of the first point to get. (0 = 1st  */
/*       point in the file).  If < 0, the resulting array will    */
/*       be zero padded.                                          */
/*   'numpts' is the number of points to get from the file.       */
/*       If the number of bins to read takes us past the end of   */
/*       file, the returned vector will be zero padded.           */

float *read_float_file(FILE *file, long firstpt, long numpts);
/* Return a float vector with complex data taken from a file.     */
/* Argumants:                                                     */
/*   'file' is a pointer to the file you want to access.          */
/*   'firstpt' is the number of the first point to get. (0 = 1st  */
/*       point in the file).  If < 0, the resulting array will    */
/*       be zero padded.                                          */
/*   'numpts' is the number of points to get from the file.       */
/*       If the number of bins to read takes us past the end of   */
/*       file, the returned vector will be zero padded.           */

/* The following routines are used by the routines above to do         */
/*   their thing..                                                     */


/* In select.c */

int prune_powers(float *arr, int n, int numsumpows);
/* Sets powers that are more than approx PRUNELEV standard */
/* devs above the median value to NEWLEV times the median  */
/* value.  This helps 'clean' the spectrum of high power   */
/* signals that probably have nothing to do with a phase   */
/* modulation spectrum (i.e. they are RF noise or strong   */
/* solitary pulsars.                                       */

void hpselect(unsigned long m, unsigned long n, \
	      float arr[], powindex heap[]);
/* Selects the m largest values from the array arr          */
/* and stores them and their indices in heap and heapindex. */

/* In median.c */

float median(float arr[], int n);
/* Finds the median (but messes up the array order) */

int comp_psr_to_cand(fourierprops * cand, infodata * idata, char *output, \
		     int full);
/* Compares a pulsar candidate defined by its properties found in   */
/*   *cand, and *idata with all of the pulsars in the pulsar        */
/*   database.  It returns a string (verbose if full==1) describing */
/*   the results of the search in *output.                          */

int comp_bin_to_cand(binaryprops * cand, infodata * idata, \
		     char *output, int full);
/* Compares a binary PSR candidate defined by its props found in    */
/*   *cand, and *idata with all of the pulsars in the pulsar        */
/*   database.  It returns a string (verbose if full==1) describing */
/*   the results of the search in *output.                          */

double dms2rad(int deg, int min, double sec);
/* Convert degrees, minutes and seconds into radians */

double hms2rad(int hour, int min, double sec);
/* Convert hours, minutes and seconds into radians */

double sphere_ang_diff(double ra1, double dec1, double ra2, double dec2);
/* Return the spherical angle (radians) between two RA and DECS */

/* In sorter.c */

float percolate(position * list, int nlist, int spot);
  /*  Pushes a position structure as far up a sorted list of positions */
  /*  as it needs to go to keep the list sorted.  Returns the new low  */
  /*  power in the list.                                               */

float percolate_bin(binaryprops * list, int nlist);
  /*  Pushes a binaryprops structure as far up a sorted list of structs*/
  /*  as it needs to go to keep the list sorted.  Returns the new low  */
  /*  power in the list.                                               */

/* From prep_corr.c  */

int next_good_fftlen(int N);
/* Return one of the shortest, yet best performing, FFT lengths larger
 * than N.  This assumes FFTW. */

int fftlen_from_kernwidth(int kernwidth);
  /* return the length of the optimal FFT to use for correlations with
   * some kernel width kernwidth.  This assumes FFTW. */

void spread_with_pad(fcomplex *data, int numdata,   \
		     fcomplex *result, int numresult, \
		     int numbetween, int numpad);
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array and padding it.                   */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */
  /*   'numpad' is the number of bins to use as zero padding     */

void spread_no_pad(fcomplex *data, int numdata, \
		   fcomplex *result, int numresult, \
		   int numbetween);
  /* Prepare the data array for correlation by spreading         */
  /*      the input data array.                                  */
  /* Arguments:                                                  */
  /*   'data' is the FFT array to be prepared                    */
  /*   'numdata' is the number of complex points in 'data'       */
  /*   'result' is the prepped data array                        */
  /*   'numresult' is the number of complex points in 'result'   */
  /*   'numbetween' is the number of interpolated pts per bin    */


void paddata(fcomplex *data, int numdata, int numpad);
  /* Pad the last 'numpad' bins of 'data' with zeros.         */
  /* Arguments:                                               */
  /*   'data' is the FFT array to be padded                   */
  /*   'numdata' is the number of complex points in 'data'    */
  /*   'numpad' is the number of bins to use as zero padding  */

void place_complex_kernel(fcomplex *kernel, int numkernel, \
			  fcomplex *result, int numresult);
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a complex response function.  Bin zero  */
  /*      response is in bin numkernel/2.                  */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */

void place_real_kernel(float *kernel, int numkernel, \
		       float *result, int numresult);
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a real-valued response function.  Bin   */
  /*      zero response is in bin numkernel/2.             */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */

void chop_complex_ends(fcomplex *data, int numdata, \
		       fcomplex *result, int numresult, \
		       int chopbins);
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */

void chop_real_ends(float *data, int numdata, \
		    float *result, int numresult, \
		    int chopbins);
  /* Chop the contaminated ends off of an array that has  */
  /* been correlated/convolved.                           */
  /* Arguments:                                           */
  /*   'data' is the array to chop.                       */
  /*   'numdata' is the number of points in data.         */
  /*   'result' is the resultant array.                   */
  /*   'numresult' is the number of points in the result. */
  /*   'chopbins' is the number of bins to chop on each   */
  /*      end of the data array.                          */



/* In correlations.c */

fcomplex *complex_corr_conv(fcomplex *data, fcomplex *kernel, \
			    int numdata, \
			    presto_ffts ffts, presto_optype type);
  /* Perform and return a complex correlation or convolution.       */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */

float *real_corr_conv(float *data, float *kernel, int numdata, \
		      presto_ffts ffts, presto_optype type);
  /* Perform and return a real-valued correlation or convolution.   */
  /* Arguments:                                                     */
  /*   'data' is the complex array to correlate/convolve.           */
  /*   'kernel' is the correlation/convolution kernel.              */
  /*   'numdata' is the length of 'data', 'kernel' and the result.  */
  /*   'ffts' describes how to perform the convolution/correlation. */
  /*      'ffts' = FFTDK:  FFT both the 'data' and the 'kernel'.    */
  /*      'ffts' = FFTD:  FFT only the 'data' not the 'kernel'.     */
  /*      'ffts' = FFTK:  FFT only the 'kernel' not the 'data'.     */
  /*      'ffts' = NOFFTS:  Don't FFT the 'data' or the 'kernel'.   */
  /*   'type' is the type of operation to perform.                  */
  /*      'type' = CONV:  return a convolution in a new vector.     */
  /*      'type' = CORR:  return a correlation in a new vector.     */
  /*      'type' = INPLACE_CONV:  convolution over-writes 'data'.   */
  /*      'type' = INPLACE_CORR:  correlation over-writes 'data'.   */


/* In corr_routines.c */

int corr_complex(fcomplex *data, int numdata, presto_datainf datainf, \
		 fcomplex *kern, int numkern, presto_datainf kerninf, \
		 fcomplex *result, int numresult, int lobin, \
		 int numbetween, int kern_half_width, presto_optype optype);
  /* This routine is a general correlation or convolution routine    */
  /* for complex data.  It can perform convolutions or correlations  */
  /* on raw complex data, data that is prepared for a convolution/   */
  /* correlation but not FFTd, or already FFTd data.  The kernel     */
  /* that it uses can also be raw, prepped, or FFTd.  If you call    */
  /* the routine multiple times with either the same kernel or data  */
  /* array, it uses a saved version of the array from the previous   */
  /* call to cut down on many processing steps. The return value     */
  /* tells how many usable (i.e.  non-contaminated) points were      */
  /* returned in the result array (the first value will be that of   */
  /* 'lobin').  This routine will _not_ perform in-place             */
  /* correlations or convolutions (i.e. it ignores those choices     */
  /* for 'optype').                                                  */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points in 'data'.          */
  /*   'datainf' is one of the following that describes the data:    */
  /*              RAW = Normal un-altered complex data.              */
  /*              PREPPED = Data has been padded and spread based    */
  /*                        on 'kern_half_width' and 'numbetween'    */
  /*                        and is ready to be FFTd.                 */
  /*              FFT = Data has already been prepared and FFTd.     */
  /*              SAME = Data is the same as the previous call.      */
  /*                        The routine uses its saved data.         */
  /*   'kern' is the correlation kernel.                             */
  /*   'numkern' is the number of complex points in 'kern'.          */
  /*   'kerninf' is one of the same choices as 'datainf' above.      */
  /*   'result' is the resulting complex array (must already exist). */
  /*   'numresult' is the number of complex points in 'result'.      */
  /*   'lobin' is the lowest fourier bin to convolve/correlate.      */
  /*   'numbetween' is the number of bins to spread the data points. */
  /*   'kern_half_width' is half the width (bins) of the raw kernel. */
  /*   'optype' is either CORR or CONV (correlation or convolution). */
  /* Notes:                                                          */
  /*   If either 'datainf' or 'kerninf' are of type PREPPED or FFT,  */
  /*   then the length of the FFTs used in the correlation/          */
  /*   convolution calculations will be of length 'numdata' or       */
  /*   'numkern'.  If both 'datainf' and 'kerninf' are of type       */
  /*   PREPPED or FFT then 'numdata' and 'numkern' must have the     */
  /*   same value.  In order for SAME values of 'datainf' and        */
  /*   'kerninf' to help out, the routine must be called with the    */
  /*   same values for 'kern_half_width' and 'numbetween' as well.   */

void stretch_fft(fcomplex *data, int numdata, \
		 fcomplex *result, int numresult);
  /* This routine stretches and/or interpolates an FFT of length    */
  /* numdata.  It zeros 'result' where end-effects have occurred.   */
  /* This routine is usually used to co-add stretched FFTs to       */
  /* increase the signal-to-noise ratios of a detection.            */
  /* Arguments:                                                     */
  /*   'data' is a pointer to a complex FFT.                        */
  /*   'numdata' is the number of complex points in 'data'.         */
  /*   'result' is a pointer to the complex stretched array.        */
  /*   'numresult' is the number of complex points in 'result'.     */
  /* Notes:                                                         */
  /*   The ratio of 'numresult' to 'numdata' determines the amount  */
  /*   of stretching that will take place.  For example, if         */
  /*   'numresult' is twice 'numdata', then the data will be        */
  /*   stretched by a factor of two (i.e. interbinned).             */

float *corr_loc_pow(float *powers, int numpowers);
  /* This routine determines the local power levels for every         */
  /* frequency in an FFT containing 'numpowers' complex frequencies.  */
  /* It sets the areas where end-effects are a problem to the         */
  /* local power level of the closest bin without end effect          */
  /* problems.  It returns a vector with the local power levels.      */
  /* Arguments:                                                       */
  /*   'powers' is a pointer to a fcomplex vector containing the FFT. */
  /*   'numpowers' is the number of complex points in 'powers'.       */


/* In rzinterp.c */

fcomplex **corr_rz_plane(fcomplex *data, int numdata, int numbetween, \
			 int startbin, double zlo, double zhi, \
			 int numz, int fftlen, \
			 presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do Fourier          */
  /* complex interpolations of the f-fdot plane.                     */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'zlo' is the lowest fdot to use (z=f-dot/T^2)                 */
  /*   'zhi' is the highest fdot to use (z=f-dot/T^2)                */
  /*   'numz' is the number of z values to use to make the plane     */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

fcomplex *corr_rz_interp(fcomplex *data, int numdata, int numbetween, \
			 int startbin, double z, int fftlen, \
			 presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation of a slice of the f-fdot plane.           */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'z' is the fdot to use (z=f-dot/T^2).                         */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

fcomplex ***corr_rzw_vol(fcomplex * data, int numdata, int numbetween, \
                         int startbin, double zlo, double zhi, int numz, \
                         double wlo, double whi, int numw, int fftlen, \
                         presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do Fourier          */
  /* complex interpolations of the f-fdot-fdotdot volume.            */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'zlo' is the lowest fdot to use (z=f-dot*T^2)                 */
  /*   'zhi' is the highest fdot to use (z=f-dot*T^2)                */
  /*   'numz' is the number of z values to use to make the volume    */
  /*   'wlo' is the lowest fdotdot to use (w=f-dotdot*T^3)           */
  /*   'whi' is the highest fdotdot to use (w=f-dotdot*T^3)          */
  /*   'numw' is the number of w values to use to make the volume    */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

    fcomplex *corr_rzw_interp(fcomplex *data, int numdata, int numbetween, \
			  int startbin, double z, double w, int fftlen, \
			  presto_interp_acc accuracy, int *nextbin);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation of part of the f-fdot-fdotdot volume.     */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'numbetween' is the number of points to interpolate per bin.  */
  /*   'startbin' is the first bin to use in data for interpolation. */
  /*   'z' is the fdot to use (z=f-dot*T^2).                         */
  /*   'w' is the fdotdot to use (z=f-dotdot*T^3).                   */
  /*   'fftlen' is the # of complex pts in kernel and result.        */
  /*   'accuracy' is either HIGHACC or LOWACC.                       */
  /*   'nextbin' will contain the bin number of the first bin not    */
  /*      interpolated in data.                                      */

void rz_interp(fcomplex *data, int numdata, double r, double z, \
	       int kern_half_width, fcomplex *ans);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation at a single point in the f-fdot plane.    */
  /* It does the correlations manually. (i.e. no FFTs)               */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'r' is the Fourier frequency in data that we want to          */
  /*      interpolate.  This can (and should) be fractional.         */
  /*   'z' is the fdot to use (z=f-dot/T^2 (T is integration time)). */
  /*   'kern_half_width' is the half-width of the kernel in bins.    */
  /*   'ans' is the complex answer.                                  */

void rzw_interp(fcomplex *data, int numdata, double r, double z, \
	       double w, int kern_half_width, fcomplex *ans);
  /* This routine uses the correlation method to do a Fourier        */
  /* complex interpolation at a single point in the f-fdot plane.    */
  /* It does the correlations manually. (i.e. no FFTs)               */
  /* Arguments:                                                      */
  /*   'data' is a complex array of the data to be interpolated.     */
  /*   'numdata' is the number of complex points (bins) in data.     */
  /*   'r' is the Fourier frequency in data that we want to          */
  /*      interpolate.  This can (and should) be fractional.         */
  /*   'z' is the fdot to use (z=f-dot*T^2 (T is integration time)). */
  /*   'w' is the fdotdot to use (z=f-dotdot*T^3).                   */
  /*   'kern_half_width' is the half-width of the kernel in bins.    */
  /*   'ans' is the complex answer.                                  */

/* In maximize_r.c and maximize_rw.c */

double max_r_arr(fcomplex *data, long numdata, double rin, 
		 double *rout, rderivs *derivs);
  /* Return the Fourier frequency that maximizes the power.  */


double max_rz_arr(fcomplex *data, long numdata, double rin, double zin, \
		  double *rout, double *zout, rderivs * derivs);
  /* Return the Fourier frequency and Fourier f-dot that      */ 
  /* maximizes the power.                                     */

double max_rz_file(FILE *fftfile, double rin, double zin, \
		   double *rout, double *zout, rderivs * derivs);
  /* Return the Fourier frequency and Fourier f-dot that      */ 
  /* maximizes the power of the candidate in 'fftfile'.       */


void max_rz_arr_harmonics(fcomplex data[], long numdata,
                          int num_harmonics,
                          double rin, double zin,
                          double *rout, double *zout,
                          rderivs derivs[], double powers[]);
/* Return the Fourier frequency and Fourier f-dot that      */
/* maximizes the power.                                     */

void max_rzw_arr_harmonics(fcomplex data[], long numdata,
                           int num_harmonics,
                           double rin, double zin, double win,
                           double *rout, double *zout, double *wout,
                           rderivs derivs[], double powers[]);
/* Return the Fourier frequency, f-dot, and f-dotdot that       */
/* maximizes the *summed* power of the multi-harmonic candidate */

double max_rzw_arr(fcomplex *data, long numdata, double rin, double zin, \
                   double win, double *rout, double *zout,              \
                   double *wout, rderivs * derivs);
/* Return the Fourier frequency, f-dot, and fdotdot that    */ 
/* maximizes the power.                                     */

double max_rzw_file(FILE * fftfile, double rin, double zin, double win, \
                    double *rout, double *zout, double *wout, rderivs * derivs);
/* Return the Fourier frequency, f-dot, and fdotdot that    */
/* maximizes the power of the candidate in 'fftfile'.       */


/* In fold.c */

void fold_errors(double *prof, int proflen, double dt, double N, 
		 double datavar, double p, double pd, double pdd, 
		 double *perr, double *pderr, double *pdderr);
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

double foldfile(FILE *datafile, double dt, double tlo, 
		double *prof, int numprof, double startphs, 
		double fo, double fdot, double fdotdot, int flags, 
		double *delays, double *delaytimes, int numdelays, 
		double *onoffpairs, foldstats *stats, float *chiarr);
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


double simplefold(float *data, int numdata, double dt, double tlo,
		  double *prof, int numprof, double startphase, 
		  double fo, double fdot, double fdotdot, int standard);
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
/*    'standard' If true, uses classic prepfold 'drizzling'      */
/*       Otherwise, adds full sample to nearest bin.             */
/* Notes:  fo, fdot, and fdotdot correspond to 'tlo' = 0.0       */
/*    (i.e. to the beginning of the first data point)            */

double fold(float *data, int numdata, double dt, double tlo, 
	    double *prof, int numprof, double startphs, 
	    double *buffer, double *phaseadded, 
	    double fo, double fdot, double fdotdot, int flags, 
	    double *delays, double *delaytimes, int numdelays, 
	    int *onoffpairs, foldstats *stats, int standard);
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
/*    'standard' If true, uses classic prepfold 'drizzling'           */
/*            Otherwise, adds full sample to nearest bin.             */
/* Notes:  fo, fdot, and fdotdot correspond to 'tlo' = 0.0            */
/*    (i.e. to the beginning of the first data point)                 */

void shift_prof(double *prof, int proflen, int shift, double *outprof);
/* Rotates a profile 'prof' by an integer 'shift' places.    */
/* If 'shift' < 0 then shift left, 'shift' > 0, shift right. */ 
/* Place the shifted  profile in 'outprof'.                  */

void combine_profs(double *profs, foldstats *instats, int numprofs, 
		   int proflen, double *delays, double *outprof,
		   foldstats *outstats);
/* Combine a series of 'numprofs' profiles, each of length 'proflen',   */
/* into a single profile of length 'proflen'.  The profiles are         */
/* summed after the appropriate 'delays' are added to each profile.     */
/* The result is a profile in 'outprof' (which must be pre-allocated)   */
/* The input stats in 'instats' are combined and placed in 'outstats'   */

void initialize_foldstats(foldstats *stats);
/* Zeroize all of the components of stats */


double doppler(double freq_observed, double voverc);
  /* This routine returns the frequency emitted by a pulsar */
  /* (in MHz) given that we observe the pulsar at frequency */
  /* freq_observed (MHz) while moving with radial velocity  */
  /* (in units of v/c) of voverc wrt the pulsar.            */


void barycenter(double *topotimes, double *barytimes, \
		double *voverc, long N, char *ra, char *dec, \
		char *obs, char *ephem);
  /* This routine uses TEMPO to correct a vector of           */
  /* topocentric times (in *topotimes) to barycentric times   */
  /* (in *barytimes) assuming an infinite observation         */
  /* frequency.  The routine also returns values for the      */
  /* radial velocity of the observation site (in units of     */
  /* v/c) at the barycentric times.  All three vectors must   */
  /* be initialized prior to calling.  The vector length for  */
  /* all the vectors is 'N' points.  The RA and DEC (J2000)   */
  /* of the observed object are passed as strings in the      */
  /* following format: "hh:mm:ss.ssss" for RA and             */
  /* "dd:mm:s.ssss" for DEC.  The observatory site is passed  */
  /* as a 2 letter ITOA code.  This observatory code must be  */
  /* found in obsys.dat (in the TEMPO paths).  The ephemeris  */
  /* found in obsys.dat (in the TEMPO paths).  The ephemeris  */
  /* is the full name of an ephemeris supported by TEMPO,     */
  /* examples include DE200, DE421, or DE436.                 */

fftcand *search_fft(fcomplex *fft, int numfft, int lobin, int hibin, 
		    int numharmsum, int numbetween, 
		    presto_interptype interptype,
		    float norm, float sigmacutoff, int *numcands, 
		    float *powavg, float *powvar, float *powmax);
/* This routine searches a short FFT of 'numfft' complex freqs      */
/* and returns a candidate vector of fftcand structures containing  */
/* information about the best candidates found.                     */
/* The routine uses either interbinning or interpolation as well    */
/* as harmonic summing during the search.                           */
/* The number of candidates returned is either 'numcands' if != 0,  */
/* or is determined automatically by 'sigmacutoff' -- which         */
/* takes into account the number of bins searched.                  */
/* The returned vector is sorted in order of decreasing power.      */
/* Arguments:                                                       */
/*   'fft' is the FFT to search (complex valued)                    */
/*   'numfft' is the number of complex points in 'fft'              */
/*   'lobin' is the lowest Fourier freq to search                   */
/*   'hibin' is the highest Fourier freq to search                  */
/*   'numharmsum' the number of harmonics to sum during the search  */
/*   'numbetween' the points to interpolate per bin                 */
/*   'interptype' is either INTERBIN or INTERPOLATE.                */
/*      INTERBIN = (interbinning) is fast but less sensitive.       */
/*      INTERPOLATE = (Fourier interpolation) is slower but more    */
/*        sensitive.                                                */
/*   'norm' is the normalization constant to multiply each power by */
/*   'sigmacutoff' if the number of candidates will be determined   */
/*      automatically, is the minimum Gaussian significance of      */
/*      candidates to keep -- taking into account the number of     */
/*      bins searched                                               */
/*   'numcands' if !0, is the number of candates to return.         */
/*      if 0, is a return value giving the number of candidates.    */
/*   'powavg' is a return value giving the average power level      */
/*   'powvar' is a return value giving the power level variance     */
/*   'powmax' is a return value giving the maximum power            */

void search_minifft(fcomplex *minifft, int numminifft, 
		    double min_orb_p, double max_orb_p,
		    rawbincand *cands, int numcands, int numharmsum,
		    int numbetween, double numfullfft, double timefullfft,
		    double lorfullfft, presto_interptype interptype,
		    presto_checkaliased checkaliased);
  /* This routine searches a short FFT (usually produced using the   */
  /* MiniFFT binary search method) and returns a candidte vector     */
  /* containing information about the best binary candidates found.  */
  /* The routine uses either interbinning or interpolation as well   */
  /* as harmonic summing during the search.                          */
  /* Arguments:                                                      */
  /*   'minifft' is the FFT to search (complex valued)               */
  /*   'numminifft' is the number of complex points in 'minifft'     */
  /*   'min_orb_p' is the minimum orbital period (s) to search       */
  /*   'max_orb_p' is the maximum orbital period (s) to search       */
  /*   'cands' is a pre-allocated vector of rawbincand type in which */
  /*      the sorted (in decreasing sigma) candidates are returned   */
  /*   'numcands' is the length of the 'cands' vector                */
  /*   'numharmsum' the number of harmonics to sum during the search */
  /*   'numbetween' the points to interpolate per bin                */
  /*   'numfullfft' the number of points in the original long FFT    */
  /*   'timefullfft' the duration of the original time series (s)    */
  /*   'lorfullfft' the 1st bin of the long FFT that was miniFFT'd   */
  /*   'interptype' is either INTERBIN or INTERPOLATE.               */
  /*      INTERBIN = (interbinning) is fast but less sensitive.      */
  /*      INTERPOLATE = (Fourier interpolation) is slower but more   */
  /*        sensitive.                                               */
  /*   'checkaliased' is either CHECK_ALIASED or NO_CHECK_ALIASED.   */
  /*      NO_CHECK_ALIASED = harmonic summing does not include       */
  /*        aliased freqs making it faster but less sensitive.       */
  /*      CHECK_ALIASED = harmonic summing includes aliased freqs    */
  /*        making it slower but more sensitive.                     */

void print_rawbincand(rawbincand cand);

int get_birdies(char *zapfilenm, double T, double avg_vel,
		double **lobins, double **hibins);
/* Open, read, and close a text file containing frequencies (Hz)   */
/* and widths (Hz) to ignore in a pulsar search.  The text file    */
/* should have one frequency and width per line.  Lines beginning  */
/* with '#' are ignored, and so may be used as comments.           */
/* 'T' is the total length in seconds of the observation that was  */
/* FFTd.  'avg_vel' is the avg topocentric velocity (in units      */
/* of c) towards the target during the obs.  The returned arrays   */
/* are sorted in order of increasing 'lobins' and contain the low  */
/* and high Fourier freqs that mark the boundaries of the birdies  */
/* (based on 'T'and 'avg_vel').                                    */

int get_std_birds(char *zapfilenm, double T, double avg_vel,
		  double **basebin, int **numharm);
/* Open, read, and close a text file containing frequencies (Hz)   */
/* and the number of harmonics to zap from a FFT.  The text file   */
/* should have one frequency and number of harmonics per line.     */
/* Lines beginning with '#' are ignored (i.e. used as comments).   */
/* 'T' is the total length in seconds of the observation that was  */
/* FFTd.  'avg_vel' is the avg topocentric velocity (in units      */
/* of c) towards the target during the obs.  The returned arrays   */
/* are sorted in order of increasing 'basebins' and contain the    */
/* base Fourier freq and the number of harmonics to check.  The    */
/* base freqs are adjusted based on avg_vel.                       */

int check_to_zap(double candbin, double *lobins, double *hibins, 
		 int numzap);
/* Look at the closest birdies from the zapfile to see if our  */
/* candidate matches one of them.  If it does, return '1' for  */
/* TRUE.  If it doesn't match, return a '0' for FALSE.  Note   */
/* that the zapfreqs _must_ be in increasing order of 'lobins' */
/* since this routine keeps track of its place in the file.    */
/* Also, numzap _must be >= 2.                                 */

short transpose_float(float *a, int nx, int ny, unsigned char *move, 
		      int move_size);
/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 *
 * "a" is a 1D array of length ny*nx which contains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */

/* NEW Clipping Routine (uses channel running averages) */
int new_clip_times(float *rawdata, int ptsperblk, int numchan, 
                   float clip_sigma, float *good_chan_levels);
// Perform time-domain clipping of rawdata.  This is a 2D array with
// ptsperblk*numchan points, each of which is a float.  The clipping
// is done at clip_sigma sigma above/below the running mean.  The
// up-to-date running averages of the channels are returned in
// good_chan_levels (which must be pre-allocated).

/* Old Clipping Routine (uses channel medians) */
int clip_times(float *rawdata, int ptsperblk, int numchan, 
               float clip_sigma, float *good_chan_levels);
// Perform time-domain clipping of rawdata.  This is a 2D array with
// ptsperblk*numchan points, each of which is a float.  The clipping
// is done at clip_sigma sigma above/below the running mean.  The
// up-to-date running averages of the channels are returned in
// good_chan_levels (which must be pre-allocated).

double *events_fdot_correct(double *events, int Nevents, 
                            double freq, double fdot);
/* Correct a set of sorted events (in sec) for a specific */
/* 'fdot' at the frequency 'freq' as per Chandler et al., */
/* 2001.  tnew_i = t_i + 0.5*fdot/freq*t_i^2.  Return a   */
/* new array of events.                                   */

fcomplex *atwood_search(double *events, double *weights, 
                        int Nevents, int Nwin, double dt);
/* Perform the time-differencing, incoherent, autocorrelation-like */
/* search for sparse event data described in                       */
/* Atwood et al. 2006, ApJL, 652, 49                               */
/*    events:  a sorted, double prec, array of event times in sec  */
/*    weights:  a weight factor (0-1) for each of the events       */
/*    Nevents:  the number of events                               */
/*    Nwin:  number of bins that make up a "window" (the FFT len)  */
/*    dt:  the time duration to use for the binning                */



