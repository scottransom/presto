%module presto

%{
#define SWIG_FILE_WITH_INIT
%}

%{
#include "presto.h"
#include "errno.h"

// A few function declarations from some functions not in headers
int fresnl(double xxa, double *ssa, double *cca);
int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);
%}

%include "numpy.i"
%include "typemaps.i"

%init %{
    import_array();
%}

typedef struct FCOMPLEX {
    float r, i;
} fcomplex;

#define SQRT2         1.4142135623730950488016887242096980785696718753769
#define PI            3.1415926535897932384626433832795028841971693993751
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#define SOL           299792458.0
#define SECPERJULYR   31557600.0
#define SECPERDAY     86400.0
#define ARCSEC2RAD    4.8481368110953599358991410235794797595635330237270e-6
#define SEC2RAD       7.2722052166430399038487115353692196393452995355905e-5

typedef enum {
  LOWACC, HIGHACC
} presto_interp_acc;

typedef enum {
  INTERBIN, INTERPOLATE
} presto_interptype;

typedef enum {
  NO_CHECK_ALIASED, CHECK_ALIASED
} presto_checkaliased;

typedef enum {
  CONV, CORR, INPLACE_CONV, INPLACE_CORR
} presto_optype;

typedef enum {
  FFTDK, FFTD, FFTK, NOFFTS
} presto_ffts;

typedef enum {
  RAW, PREPPED, FFT, SAME
} presto_datainf;

// Some functions from ransomfft.h

void read_wisdom(void);
long long good_factor(long long nn);

// setup a typemap for fcomplex to numpy complex
%numpy_typemaps(fcomplex, NPY_CFLOAT, long)
%apply (fcomplex* INPLACE_ARRAY1, long DIM1) {(fcomplex *indata, long nn)};
void fftwcall(fcomplex *indata, long nn, int isign);
void tablesixstepfft(fcomplex *indata, long nn, int isign);
%clear (fcomplex *indata, long nn);

%numpy_typemaps(float, NPY_FLOAT, long)
%apply (float* INPLACE_ARRAY1, long DIM1) {(float *data, long n)};
void realfft(float *data, long n, int isign);
%clear (float *data, long n);

typedef struct INFODATA {
    double ra_s;          /* Right ascension seconds (J2000)       */
    double dec_s;         /* Declination seconds (J2000)           */ 
    double N;             /* Number of bins in the time series     */
    double dt;            /* Width of each time series bin (sec)   */
    double fov;           /* Diameter of Beam or FOV in arcsec     */
    double mjd_f;         /* Epoch of observation (MJD) frac part  */
    double dm;            /* Radio -- Dispersion Measure (cm-3 pc) */
    double freq;          /* Radio -- Low chan central freq (Mhz)  */
    double freqband;      /* Radio -- Total Bandwidth (Mhz)        */
    double chan_wid;      /* Radio -- Channel Bandwidth (Mhz)      */
    double wavelen;       /* IR,Opt,UV -- central wavelength (nm)  */
    double waveband;      /* IR,Opt,UV -- bandpass (nm)            */
    double energy;        /* x-ray,gamma -- central energy (kev)   */
    double energyband;    /* x-ray,gamma -- energy bandpass (kev)  */
    // Don't handle onoff for right now
    // double *onoff;        /* Bin number pairs where obs is "on"    */
    int num_chan;         /* Radio -- Number Channels              */
    int mjd_i;            /* Epoch of observation (MJD) int part   */
    int ra_h;             /* Right ascension hours (J2000)         */
    int ra_m;             /* Right ascension minutes (J2000)       */
    int dec_d;            /* Declination degrees (J2000)           */
    int dec_m;            /* Declination minutes (J2000)           */  
    int bary;             /* Barycentered?  1=yes, 0=no            */
    int numonoff;         /* The number of onoff pairs in the data */ 
    %extend {
        char *notes;      /* Any additional notes                  */
        char *name;       /* Data file name without suffix         */
        char *object;     /* Object being observed                 */ 
        char *instrument; /* Instrument used                       */
        char *observer;	  /* Observer[s] for the data set          */
        char *analyzer;	  /* Who analyzed the data                 */
        char *telescope;  /* Telescope used                        */
        char *band;       /* Type of observation (EM band)         */
        char *filt;       /* IR,Opt,UV -- Photometric Filter       */
    }
} infodata;
%{
    //double __getitem__(int index) {
    //    return self->onoff[index];
    //}
    //void __setitem__(int index, double val) {
    //    self->onoff[index] = val;
    //}
    char *INFODATA_notes_get(infodata *p) {
        return p->notes;
    }
    void INFODATA_notes_set(infodata *p, char *val) {
        strncpy(p->notes,val,500);
    }
    char *INFODATA_name_get(infodata *p) {
        return p->name;
    }
    void INFODATA_name_set(infodata *p, char *val) {
        strncpy(p->name,val,200);
    }
    char *INFODATA_object_get(infodata *p) {
        return p->object;
    }
    void INFODATA_object_set(infodata *p, char *val) {
        strncpy(p->object,val,100);
    }
    char *INFODATA_instrument_get(infodata *p) {
        return p->instrument;
    }
    void INFODATA_instrument_set(infodata *p, char *val) {
        strncpy(p->instrument,val,100);
    }
    char *INFODATA_observer_get(infodata *p) {
        return p->observer;
    }
    void INFODATA_observer_set(infodata *p, char *val) {
        strncpy(p->observer,val,100);
    }
    char *INFODATA_analyzer_get(infodata *p) {
        return p->analyzer;
    }
    void INFODATA_analyzer_set(infodata *p, char *val) {
        strncpy(p->analyzer,val,100);
    }
    char *INFODATA_telescope_get(infodata *p) {
        return p->telescope;
    }
    void INFODATA_telescope_set(infodata *p, char *val) {
        strncpy(p->telescope,val,40);
    }
    char *INFODATA_band_get(infodata *p) {
        return p->band;
    }
    void INFODATA_band_set(infodata *p, char *val) {
        strncpy(p->band,val,40);
    }
    char *INFODATA_filt_get(infodata *p) {
        return p->filt;
    }
    void INFODATA_filt_set(infodata *p, char *val) {
        strncpy(p->filt,val,7);
    }
%}

void readinf(infodata *data, char *filenm);
/* Read a ".inf" file into memory. */

void writeinf(infodata *data);
/* Write a ".inf" file to disk. */

typedef struct orbitparams {
  double p;         /* Orbital period (s)                            */
  double e;         /* Orbital eccentricity                          */
  double x;         /* Projected semi-major axis (lt-sec)            */
  double w;         /* Longitude of periapsis (deg)                  */
  double t;         /* Time since last periastron passage (s)        */
  double pd;        /* Orbital period derivative (s/yr)              */
  double wd;        /* Advance of longitude of periapsis (deg/yr)    */
} orbitparams;

/* This is a structure that contains the 'key' pulsar info  */
typedef struct PSRPARAMS {
    %extend {
        char *jname;           /* The PSRs J2000 name         */
        char *bname;           /* The PSRs B1950 name         */
        char *alias;           /* An alias for the pulsar     */
    }
    double ra2000;             /* J2000 RA                    */
    double dec2000;            /* J2000 DEC                   */
    double dm;                 /* Dispersion Measure          */
    double timepoch;           /* MJD epoch for timing        */
    double p;                  /* PSR period (s)              */
    double pd;                 /* PSR period deriv (s/s)      */
    double pdd;                /* Period 2nd deriv (s/s^2)    */
    double f;                  /* PSR frequency (hz)          */
    double fd;                 /* PSR frequency deriv (s^-2)  */
    double fdd;                /* Frequency 2nd deriv (s^-3)  */
    orbitparams orb;            /* Orbital parameters          */
} psrparams;
%{
    char *PSRPARAMS_jname_get(psrparams *p) {
        return p->jname;
    }
    void PSRPARAMS_jname_set(psrparams *p, char *val) {
        strncpy(p->jname,val,13);
    }
    char *PSRPARAMS_bname_get(psrparams *p) {
        return p->bname;
    }
    void PSRPARAMS_bname_set(psrparams *p, char *val) {
        strncpy(p->bname,val,9);
    }
    char *PSRPARAMS_alias_get(psrparams *p) {
        return p->jname;
    }
    void PSRPARAMS_alias_set(psrparams *p, char *val) {
        strncpy(p->jname,val,10);
    }
%}

int get_psr_at_epoch(char *psrname, double epoch, psrparams *psr);
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.                                      */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */

int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);
/* Converts info from a "par" file to the "current" epoch.  */
/* Returned values go in *psr.  The name of the parfile is  */
/* in *parfilenm. epoch is the time in question in MJD.     */
/* The int returned is 1 if successful, 0 otherwise.        */

void mjd_to_datestr(double mjd, char *datestr);
// Convert an MJD to a PSRFITS-style DATEOBS

%apply double *OUTPUT { double *ssa,  double *cca };
int fresnl(double xxa, double *ssa, double *cca);
// Return the Fresnel inegrals

typedef struct RDERIVS {
  double pow;       /* Power normalized with local power             */
  double phs;       /* Signal phase                                  */
  double dpow;      /* 1st deriv of power wrt fourier freq           */
  double dphs;      /* 1st deriv of phase wrt fourier freq           */
  double d2pow;     /* 2nd deriv of power wrt fourier freq           */
  double d2phs;     /* 2nd deriv of power wrt fourier freq           */
  double locpow;    /* Local mean power level                        */
} rderivs;

typedef struct FOURIERPROPS {
  double r;        /* Fourier freq at start of run                  */
  float rerr;      /* Error in fourier freq                         */
  double z;        /* Fourier freq deriv at start of run = T^2*fdot */
  float zerr;      /* Error in fourier freq deriv                   */
  double w;        /* Fourier 2nd freq deriv = T^3*fdotdot          */
  float werr;      /* Error in 2nd fourier freq deriv               */
  float pow;       /* Power normalized with local power             */
  float powerr;    /* Error in power                                */
  float sig;       /* Approx. Sigma level                           */
  float rawpow;    /* Raw signal power (unnormalized)               */
  float phs;       /* Signal phase                                  */
  float phserr;    /* Error in phase                                */
  float cen;       /* Centroid of signal pulsations over run        */
  float cenerr;    /* Error in centroid                             */
  float pur;       /* Purity:  rms duration of signal over run      */
  float purerr;    /* Error in purity                               */
  float locpow;    /* Local mean power level                        */
} fourierprops;

typedef struct foldstats {
  double numdata;     /* Number of data bins folded         */
  double data_avg;    /* Average level of the data bins     */
  double data_var;    /* Variance of the data bins          */
  double numprof;     /* Number of bins in the profile      */
  double prof_avg;    /* Average level of the profile bins  */
  double prof_var;    /* Variance of the profile bins       */
  double redchi;      /* Reduced chi-squared of the profile */
} foldstats;

%apply (float** ARGOUTVIEWM_ARRAY1, long* DIM1) {(float** vect, long *nn)}
%apply (fcomplex** ARGOUTVIEWM_ARRAY1, long* DIM1) {(fcomplex** vect, long *nn)}

%exception
{
    errno = 0;
    $action

    if (errno != 0)
    {
        switch(errno)
        {
            case ENOMEM:
                PyErr_Format(PyExc_MemoryError, "Failed malloc()");
                break;
            default:
                PyErr_Format(PyExc_Exception, "Unknown exception");
        }
        SWIG_fail;
    }
}

%rename (gen_fvect) wrap_gen_fvect;
%rename (gen_cvect) wrap_gen_cvect;
%inline %{
void wrap_gen_fvect(long nl, float** vect, long *nn)
{
    *vect = gen_fvect(nl);
    *nn = nl;
}
void wrap_gen_cvect(long nl, fcomplex** vect, long *nn)
{
    *vect = gen_cvect(nl);
    *nn = nl;
}
%}
%clear (float **vect, long *nn);
%clear (fcomplex **vect, long *nn);


%apply (float** ARGOUTVIEWM_ARRAY1, long* DIM1) {(float** vect, long *nn)}
%apply (fcomplex* INPLACE_ARRAY1, long DIM1) {(fcomplex *dft, long N)};
%rename (power_arr) wrap_power_arr;
%rename (phase_arr) wrap_phase_arr;
%inline %{
void wrap_power_arr(fcomplex *dft, long N, float **vect, long *nn){
    /* Determine the spectral powers of the Fourier amplitudes 'dft'*/
    float powargr, powargi, *powers;
    long ii;
    
    powers = gen_fvect(N);
    for (ii=0; ii<N; ii++)
        powers[ii] = POWER(dft[ii].r, dft[ii].i);
    *vect = powers;
    *nn = N;
}

void wrap_phase_arr(fcomplex *dft, long N, float **vect, long *nn){
    /* Determine the spectral phases of the Fourier amplitudes 'dft'*/
    float phsargr, phsargi, phstmp, *phases;
    long ii;
    
    phases = gen_fvect(N);
    for (ii=0; ii<N; ii++)
        phases[ii] = PHASE(dft[ii].r, dft[ii].i);
    *vect = phases;
    *nn = N;
}
%}
%clear (float **vect, long *nn);
%clear (fcomplex *dft, long N);

%numpy_typemaps(float, NPY_FLOAT, long)
%apply (float* INPLACE_ARRAY1, long DIM1) {(float *data, long numbins)};
void frotate(float *data, long numbins, float bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of FLOATING points to move.       */
%clear (float *data, long numbins);

%numpy_typemaps(double, NPY_DOUBLE, long)
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *data, long numbins)};
void drotate(double *data, long numbins, double bins_to_left);
/* Rotates a vector by bins_to_left places to the left.    */
/* numbins is the number of DOUBLE points to move.         */
%clear (double *data, long numbins);

/*
 *   For all the following:
 *   E           = Eccentric anomaly (radians)
 *   Eo          = initial condition of E
 *   Eacc        = accuracy (1e-14 should be plenty good)
 *   t           = Time (sec)
 *   to          = orbital integration start time
 *   tend        = orbital integration ending time
 *   dt          = sampling interval for the integrator
 *   f_orb       = orbital frequency in hertz
 *   p_orb       = orbital period in seconds
 *   p_psr       = pulsar period in units of choice
 *   e or e_orb  = orbital eccentricity
 *   x or x_orb  = projected semi-major axis of orbit in lt-sec
 *   w or w_orb  = longitude of periastron (degrees)
 */

double keplers_eqn(double t, double p_orb, double e, double Eacc);
/* This routine solves Kepler's Equation at a single time t (sec) and  */
/* returns the value of the eccentric anomaly.  The orbital period (s) */
/* is in p_orb and the orbital eccentricity is in e.  Eacc is the      */
/* absolute accuracy in E that we want to achieve.  t is the time in   */
/* seconds since the last periapsis.  Uses Newton-Raphson.             */

%numpy_typemaps(double, NPY_DOUBLE, long)
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *E, long numpoints)};
void E_to_phib(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Phase delays. */

void E_to_v(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Pulsar line-of-sight velocity (km/s). */

void E_to_p(double *E, long numpoints, double p_psr, orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Pulse period. */

void E_to_z(double *E, long numpoints, double p_psr, double T,
            orbitparams *orb);
/* Convert eccentric anomalies returned by dorbint into */
/* Fourier f-dot. */

void E_to_phib_BT(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalys (*E) to time delays */
/* using Blanford and Teukolsky Equations         */
%clear (double *E, long numpoints);

%apply (double** ARGOUTVIEWM_ARRAY1, long* DIM1) {(double** vect, long *nn)}
%rename (dorbint) wrap_dorbint;
/* This routine integrates Kepler's Equation and returns a double      */
/* vector of the eccentric anomalys (E) for each point.  The initial   */
/* value for eccentric anomaly (usually determined by using            */
/* keplers_equation()) goes in Eo.  The time increment to use is dt,   */
/* total number of pts goes in 'numpts' and all of the various orbital */
/* parameters are found in *orb.  The routine uses 4th order Runge-    */
/* Kutta in a dumb mode (no adaptive step-size) since all we want is   */
/* tabulated results with even intervals.                              */
%inline %{
void wrap_dorbint(double **vect, long *nn, double Eo, long numpts, 
                 double dt, orbitparams *orb)
{
    *vect = dorbint(Eo, numpts, dt, orb);
    *nn = numpts;
}
%}
%clear (double **vect, long *nn);

%apply double *OUTPUT { double *minv,  double *maxv };
void binary_velocity(double T, orbitparams * orbit, 
                     double *minv, double *maxv);
/*  Return the minimum and maximum orbital velocities of a pulsar    */
/*  during an observation as a fraction of the speed of light.       */
/*  Arguments:                                                       */
/*    'ppsr' is the period of the pusar in seconds.                  */
/*    'T' is the length of the observation in seconds.               */
/*    'orbit' is a ptr to a orbitparams structure containing the     */
/*       Keplerian orbital parameters of the binary system.          */

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

%apply (fcomplex** ARGOUTVIEWM_ARRAY1, long* DIM1) {(fcomplex **vect, long *nn)};
%rename (gen_r_response) wrap_gen_r_response;
%rename (gen_z_response) wrap_gen_z_response;
%rename (gen_w_response) wrap_gen_w_response;
%rename (gen_w_response2) wrap_gen_w_response2;
%rename (gen_bin_response) wrap_gen_bin_response;
%inline %{
void wrap_gen_r_response(double roffset, int numbetween, int numkern,
                         fcomplex **vect, long *nn){
    *vect = gen_r_response(roffset, numbetween, numkern);
    *nn = numkern;
}

void wrap_gen_z_response(double roffset, int numbetween, int numkern,
                         double z,
                         fcomplex **vect, long *nn){
    *vect = gen_z_response(roffset, numbetween, z, numkern);
    *nn = numkern;
}

void wrap_gen_w_response(double roffset, int numbetween, int numkern,
                         double z, double w,
                         fcomplex **vect, long *nn){
    *vect = gen_w_response(roffset, numbetween, z, w, numkern);
    *nn = numkern;
}

void wrap_gen_w_response2(double roffset, int numbetween, int numkern,
                          double z, double w,
                         fcomplex **vect, long *nn){
    *vect = gen_w_response2(roffset, numbetween, z, w, numkern);
    *nn = numkern;
}

void wrap_gen_bin_response(double roffset, int numbetween, int numkern,
                           double ppsr, double T, orbitparams *orbit,
                           fcomplex **vect, long *nn){
    *vect = gen_bin_response(roffset, numbetween, ppsr, T, 
                             orbit, numkern);
    *nn = numkern;
}
%}
%clear (fcomplex **vect, long *nn);

%apply (fcomplex* IN_ARRAY1, long DIM1) {(fcomplex *data, long numdata)};
float get_localpower(fcomplex *data, long numdata, double r);
  /* Return the local power level at specific FFT frequency.  */
  /* Arguments:                                               */
  /*   'data' is a pointer to a complex FFT.                  */
  /*   'numdata' is the number of complex points in 'data'.   */
  /*   'r' is the Fourier frequency in data that we want to   */
  /*      interpolate.                                        */

float get_localpower3d(fcomplex *data, long numdata, double r,
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
%clear (fcomplex *data, long numdata);

void get_derivs3d(fcomplex *data, long numdata, double r,
                  double z, double w, float localpower,
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

void calc_props(rderivs data, double r, double z, double w, 
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

void calc_binprops(fourierprops * props, double T, int lowbin, 
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

%{
    double extended_equiv_gaussian_sigma(double logp);
    double log_asymtotic_incomplete_gamma(double a, double z);
    double log_asymtotic_gamma(double z);
%}
double extended_equiv_gaussian_sigma(double logp);
double log_asymtotic_incomplete_gamma(double a, double z);
double log_asymtotic_gamma(double z);

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

%apply double *OUTPUT { double *out, double *outd, double *outdd };
void switch_f_and_p(double in, double ind, double indd,
                    double *out, double *outd, double *outdd);
/* Convert p, p-dot, and p-dotdot into f, f-dot, */
/* and f-dotdot or vise-versa.                   */

%apply (double* IN_ARRAY1, int DIM1) {(double *data, int numdata)};
double chisqr(double *data, int numdata, double avg, double var);
  /* Calculates the chi-square of the 'data' which has average */
  /* 'avg', and variance 'var'.                                */
%clear (double *data, int numdata);

void print_candidate(fourierprops * cand, double dt, unsigned long N, 
		     float nph, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

void print_bin_candidate(binaryprops * cand, int numerrdigits);
/* Outputs a 2 column summary of all the properties or a fourier peak  */

int read_rzw_cand(FILE *file, fourierprops *cands);
/* Read the next rzw candidate from the file */
/* If successful, return 1, else 0           */

void get_rzw_cand(char *filenm, int candnum, fourierprops * cand);
/*  Read the rzw candidate file 'filenm' and return a        */
/*  pointer to the fourierprops that describes it.           */

int read_bin_cand(FILE *file, binaryprops *cands);
/* Read the next binary candidate from the file */
/* If successful, return 1, else 0              */

void get_bin_cand(char *filenm, int candnum, binaryprops * cand);
/*  Read the bin candidate file 'filenm' and return a        */
/*  pointer to the binaryprops that describes it.            */

long long next2_to_n(long long x);
/* Return the first value of 2^n >= x */

int is_power_of_10(long long n);
/* Check whether n is a power of 10 or not.  Return 0 or 1 */

long long choose_good_N(long long orig_N);
/* Choose a time series length that is larger than the input */
/* value but that is highly factorable.                      */

double dms2rad(int deg, int min, double sec);
/* Convert degrees, minutes and seconds into radians */

double hms2rad(int hour, int min, double sec);
/* Convert hours, minutes and seconds into radians */

%apply int *OUTPUT { int *xx, int *mm };
%apply double *OUTPUT { double *ss };
void hours2hms(double hours, int *xx, int *mm, double *ss);
/* Convert decimal hours to hours, minutes, and seconds */
void deg2dms(double degrees, int *xx, int *mm, double *ss);
/* Convert decimal degrees to degrees, minutes, and seconds */

double sphere_ang_diff(double ra1, double dec1, double ra2, double dec2);
/* Return the spherical angle (radians) between two RA and DECS */

%apply (fcomplex** ARGOUTVIEWM_ARRAY2, long* DIM1, long* DIM2) \
    {(fcomplex **arr, long *nr, long *nc)};
%apply (fcomplex* INPLACE_ARRAY1, long DIM1) {(fcomplex *data, long N)};
%rename (corr_rz_plane) wrap_corr_rz_plane;
%inline %{
    void wrap_corr_rz_plane(fcomplex *data, long N, int numbetween,
                            int startbin, double zlo, double zhi,
                            int numz, int fftlen,
                            presto_interp_acc accuracy,
                            fcomplex **arr, long *nr, long *nc){
        fcomplex **outarr;
        int nextbin;
        outarr = corr_rz_plane(data, N, numbetween, startbin, zlo, zhi,
                               numz, fftlen, accuracy, &nextbin);
        *arr = outarr[0];
        *nr = numz;
        *nc = (nextbin - startbin) * numbetween;
        vect_free(outarr);
    }
%}
%clear (fcomplex **arr, long *nr, long *nc);
%clear (fcomplex *data, long N);

%apply (fcomplex** ARGOUTVIEWM_ARRAY3, long* DIM1, long* DIM2, long* DIM3) \
    {(fcomplex **arr, long *nh, long *nr, long *nc)};
%apply (fcomplex* INPLACE_ARRAY1, long DIM1) {(fcomplex *data, long N)};
%rename (corr_rzw_vol) wrap_corr_rzw_vol;
%inline %{
    void wrap_corr_rzw_vol(fcomplex *data, long N, int numbetween, \
                           int startbin, double zlo, double zhi, int numz, \
                           double wlo, double whi, int numw, int fftlen, \
                           presto_interp_acc accuracy, \
                           fcomplex **arr, long *nh, long *nr, long *nc){
        fcomplex ***outarr;
        int nextbin;
        outarr = corr_rzw_vol(data, N, numbetween, startbin, zlo, zhi,
                              numz, wlo, whi, numw, fftlen, accuracy, &nextbin);
        *arr = outarr[0][0];
        *nh = numw;
        *nr = numz;
        *nc = (nextbin - startbin) * numbetween;
        vect_free(outarr[0]);
        vect_free(outarr);
    }
%}
%clear (fcomplex **arr, long *nh, long *nr, long *nc);
%clear (fcomplex *data, long N);

%apply double *OUTPUT { double *rout, double *powout };
%apply (fcomplex* INPLACE_ARRAY1, long DIM1) {(fcomplex *data, long numdata)};
%rename (max_r_arr) wrap_max_r_arr;
%inline %{
    void wrap_max_r_arr(fcomplex * data, long numdata, double rin,
                         rderivs * derivs, double *rout, double *powout){
        double pow;
        pow = max_r_arr(data, numdata, rin, rout, derivs);
        *powout = pow;
    }
%}

%apply double *OUTPUT { double *rout, double *zout, double *powout };
%rename (max_rz_arr) wrap_max_rz_arr;
%inline %{
    void wrap_max_rz_arr(fcomplex * data, long numdata, double rin, double zin,
                         rderivs * derivs, double *rout, double *zout, double *powout){
        double pow;
        pow = max_rz_arr(data, numdata, rin, zin, rout, zout, derivs);
        *powout = pow;
    }
%}

%apply double *OUTPUT { double *rout, double *zout };
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *derivdata, int len)};
%rename (max_rz_arr_harmonics) wrap_max_rz_arr_harmonics;
%inline %{
    void wrap_max_rz_arr_harmonics(fcomplex *data, long numdata,
                                   double rin, double zin,
                                   double *derivdata, int len,
                                   double *rout, double *zout){
        int ii, numharm = len / 7;
        double *powers = gen_dvect(numharm);
        rderivs *derivs = (rderivs *)malloc(sizeof(rderivs) * numharm);
        
        max_rz_arr_harmonics(data, numdata, numharm, rin, zin, rout, zout, derivs, powers);
        vect_free(powers);
        // Hack to effectively return a array of rderivs
        for (ii = 0 ; ii < numharm ; ii++) {
            derivdata[ii*7+0] = derivs[ii].pow;
            derivdata[ii*7+1] = derivs[ii].phs;
            derivdata[ii*7+2] = derivs[ii].dpow;
            derivdata[ii*7+3] = derivs[ii].dphs;
            derivdata[ii*7+4] = derivs[ii].d2pow;
            derivdata[ii*7+5] = derivs[ii].d2phs;
            derivdata[ii*7+6] = derivs[ii].locpow;
        }
        free(derivs);
    }
%}

%apply double *OUTPUT { double *rout, double *zout, double *wout };
%rename (max_rzw_arr_harmonics) wrap_max_rzw_arr_harmonics;
%inline %{
    void wrap_max_rzw_arr_harmonics(fcomplex *data, long numdata,
                                    double rin, double zin, double win,
                                    double *derivdata, int len,
                                    double *rout, double *zout, double *wout){
        int ii, numharm = len / 7;
        double *powers = gen_dvect(numharm);
        rderivs *derivs = (rderivs *)malloc(sizeof(rderivs) * numharm);
        
        max_rzw_arr_harmonics(data, numdata, numharm, rin, zin, win,
                              rout, zout, wout, derivs, powers);
        vect_free(powers);
        // Hack to effectively return a array of rderivs
        for (ii = 0 ; ii < numharm ; ii++) {
            derivdata[ii*7+0] = derivs[ii].pow;
            derivdata[ii*7+1] = derivs[ii].phs;
            derivdata[ii*7+2] = derivs[ii].dpow;
            derivdata[ii*7+3] = derivs[ii].dphs;
            derivdata[ii*7+4] = derivs[ii].d2pow;
            derivdata[ii*7+5] = derivs[ii].d2phs;
            derivdata[ii*7+6] = derivs[ii].locpow;
        }
        free(derivs);
    }
%}
%clear (double *derivdata, int len);

%apply double *OUTPUT { double *rout, double *zout, double *wout, double *powout };
%rename (max_rzw_arr) wrap_max_rzw_arr;
%inline %{
    void wrap_max_rzw_arr(fcomplex * data, long numdata, double rin, double zin, double win,
                          rderivs * derivs, double *rout, double *zout, double *wout, double *powout){
        double pow;
        pow = max_rzw_arr(data, numdata, rin, zin, win, rout, zout, wout, derivs);
        *powout = pow;
    }
%}
%clear (fcomplex *data, long numdata);

%apply (double* INPLACE_ARRAY1, long DIM1) {(double *topotimes, long N1)};
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *barytimes, long N2)};
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *voverc, long N3)};
%rename (barycenter) wrap_barycenter;
%inline %{
    void wrap_barycenter(double *topotimes, long N1,
                         double *barytimes, long N2,
                         double *voverc, long N3,
                         char *ra, char *dec, char *obs, char *ephem){
        barycenter(topotimes, barytimes, voverc, N1, ra, dec, obs, ephem);
    }
%}
%clear (double *topotimes, long N1);
%clear (double *barytimes, long N2);
%clear (double *voverc, long N3);

%apply (float* INPLACE_ARRAY1, long DIM1) {(float *data, long N1)};
%apply (double* INPLACE_ARRAY1, long DIM1) {(double *prof, long N2)};
%rename (simplefold) wrap_simplefold;
%inline %{
    double wrap_simplefold(float *data, long N1, double dt, double tlo,
                           double *prof, long N2, double startphs,
                           double f0, double fdot, double fdotdot,
                           int standard){
        return simplefold(data, N1, dt, tlo,
                          prof, N2, startphs,
                          f0, fdot, fdotdot, standard);
    }
%}
%clear (float *data, long N1);
%clear (double *prof, long N2);

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
