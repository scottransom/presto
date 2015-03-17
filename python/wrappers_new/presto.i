%module newpresto

%{
#define SWIG_FILE_WITH_INIT
%}

%{
#include "presto.h"
%}

%include "numpy.i"

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

%apply (long* DIM1, float** ARGOUTVIEWM_ARRAY1) {(long *nn, float** vect)}
%apply (long* DIM1, fcomplex** ARGOUTVIEWM_ARRAY1) {(long *nn, fcomplex** vect)}
%rename (gen_fvect) wrap_gen_fvect;
%rename (gen_cvect) wrap_gen_cvect;
%inline %{
#include <errno.h>
void wrap_gen_fvect(long nl, long *nn, float** vect)
{
    float *temp;
    temp = gen_fvect(nl);
    if (temp == NULL)
        errno = ENOMEM;
    *vect = temp;
    *nn = nl;
}
void wrap_gen_cvect(long nl, long *nn, fcomplex** vect)
{
    fcomplex *temp;
    temp = gen_cvect(nl);
    if (temp == NULL)
        errno = ENOMEM;
    *vect = temp;
    *nn = nl;
}
%}
%clear (long *nn, float** vect);
%clear (long *nn, fcomplex** vect);
