
%inline %{
typedef struct {
     double *dptr;
} DoubleArray;
%}

%typemap(memberin) char name[200] {
        strncpy($target,$source,199);
}       
%typemap(memberin) char object[100] {
        strncpy($target,$source,99);
}       
%typemap(memberin) char telescope[40] {
        strncpy($target,$source,39);
}       
%typemap(memberin) char instrument[100] {
        strncpy($target,$source,99);
}       
%typemap(memberin) char observer[100] {
        strncpy($target,$source,99);
}
%typemap(memberin) char band[40] {
        strncpy($target,$source,39);
}       
%typemap(memberin) char filt[7] {
        strncpy($target,$source,6);
}       
%typemap(memberin) char analyzer[100] {
        strncpy($target,$source,99);
}       
%typemap(memberin) char notes[500] {
        strncpy($target,$source,499);
}       

%addmethods DoubleArray {
  DoubleArray(double *dptr) {
    DoubleArray *d = (DoubleArray *) malloc(sizeof(DoubleArray));
    d->dptr = dptr;
    return d;
  }
  ~DoubleArray() {
    free(self);
  }
  double __getitem__(int index) {
    return self->dptr[index];
  }
  void __setitem__(int index, double val) {
    self->dptr[index] = val;
  }
};

// Now a "clever" typemap

// Return a DoubleArray instead of a double *
%typemap(memberout) DoubleArray *onoff {
  static DoubleArray d;
  d.dptr = $source;
  $target = &d;
}

// Set a DoubleArray
%typemap(memberin) DoubleArray *onoff {
  $target = $source->dptr;
}

typedef struct INFODATA {
  char name[200];	/* Data file name without suffix         */
  char object[100];	/* Object being observed                 */
  int ra_h;		/* Right ascension hours (J2000)         */
  int ra_m;		/* Right ascension minutes (J2000)       */
  double ra_s;		/* Right ascension seconds (J2000)       */
  int dec_d;		/* Declination degrees (J2000)           */
  int dec_m;		/* Declination minutes (J2000)           */
  double dec_s;		/* Declination seconds (J2000)           */
  char telescope[40];	/* Telescope used                        */
  char instrument[100];	/* Instrument used                       */
  char observer[100];	/* Observer[s] for the data set          */
  double N;		/* Number of bins in the time series     */
  double dt;	 	/* Width of each time series bin (sec)   */
  int numonoff;		/* Any breaks in the data?  1=yes, 0=no  */
  DoubleArray *onoff;	/* Bin number pairs where obs is "on"    */
  double fov;		/* Diameter of Beam or FOV in arcsec     */
  int mjd_i;		/* Epoch of observation (MJD) int part   */
  double mjd_f;		/* Epoch of observation (MJD) frac part  */
  int bary;		/* Barycentered?  1=yes, 0=no            */
  char band[40];	/* Type of observation (EM band)         */
  double dm;		/* Radio -- Dispersion Measure (cm-3 pc) */
  double freq;		/* Radio -- Low chan central freq (Mhz)  */
  double freqband;	/* Radio -- Total Bandwidth (Mhz)        */
  int num_chan;		/* Radio -- Number Channels              */
  double chan_wid;	/* Radio -- Channel Bandwidth (Mhz)      */
  char filt[7];		/* IR,Opt,UV -- Photometric Filter       */
  double wavelen;	/* IR,Opt,UV -- central wavelength (nm)  */
  double waveband;	/* IR,Opt,UV -- bandpass (nm)            */
  double energy;	/* x-ray,gamma -- central energy (kev)   */
  double energyband;	/* x-ray,gamma -- energy bandpass (kev)  */
  char analyzer[100];	/* Who analyzed the data                 */
  char notes[500];	/* Any additional notes                  */
  %addmethods {
    infodata(){
      infodata *p = (infodata *)malloc(sizeof(infodata));
      return p;
    }
    ~infodata(){
      free(self);
    }
  }
} infodata;

void readinf(infodata * data, char *filenm);
/* Read a ".inf" file into memory. */

void writeinf(infodata *data);
/* Write a ".inf" file to disk. */


%typemap(memberin) char basefilenm[200] {
        strncpy($target,$source,199);
}       
%typemap(memberin) char description[200] {
        strncpy($target,$source,199);
}       
%typemap(memberin) char ptype[20] {
        strncpy($target,$source,19);
}       
%typemap(memberin) char round[20] {
        strncpy($target,$source,19);
}       
%typemap(memberin) char noisetype[20] {
        strncpy($target,$source,19);
}       

typedef struct MAKEDATA {
  char basefilenm[200];	        /* Data file name without suffix         */
  char description[200];	/* Data description                      */
  long N;		        /* Number of bins in the time series     */
  long next2_to_n;              /* The next power-of-2 >= N              */
  double dt;	 	        /* Width of each time series bin (sec)   */
  double T;                     /* Total length of time series (sec)     */
  char ptype[20];		/* Pulsetype (Sine, Crab-like, Spike)    */
  int pnum;                     /* 1=Sine, 2=Crab-like, 3=Spike          */
  double fwhm;                  /* FWHM Phase (0-1) if ptype is Spike    */
  char round[20];		/* Rounding Format (Whole, Fractional)   */
  int roundnum;                 /* 1=Whole Numbers, 0=Fractional         */
  double f;                     /* Pulsar frequency (hz)                 */
  double fd;                    /* Pulsar frequency deriv (hz/s)         */
  double fdd;                   /* Pulsar 2nd frequency deriv (hz/s^2)   */
  double p;                     /* Pulsar period (s)                     */
  double pd;                    /* Pulsar period deriv (s/s)             */
  double pdd;                   /* Pulsar 2nd period deriv (s/s^2)       */
  double r;                     /* Pulsar FFT freq (bin)                 */
  double z;                     /* Pulsar FFT freq deriv (bins/T)        */
  double w;                     /* Pulsar FFT 2nd freq deriv (bins/T^2)  */
  double amp;                   /* Pulse amplitude                       */
  double phs;                   /* Pulse phase (deg)                     */
  double dc;                    /* Average DC background level           */
  int binary;			/* Binary pulsar?  (1=Yes, 0=No)         */
  orbitparams orb;              /* Binary pulsar data                    */
  int ampmod;                   /* Amplitude modulation? (1=Yes, 0=No)   */
  double ampmoda;               /* Amplitude modulation amplitude        */
  double ampmodf;               /* Amplitude modulation frequency (hz)   */
  double ampmodp;               /* Amplitude modulation phase (deg)      */
  char noisetype[20];		/* Noise type (Standard, Other)          */
  int noise;                    /* 1=Standard, 0=Other                   */
  double noisesig;              /* Noise standard deviation              */
  int numonoff;                 /* The number of onoff pairs in the data */
  DoubleArray *onoff;		/* Bin number pairs where obs is "on"    */
  %addmethods {
    makedata(){
      makedata *p = (makedata *)malloc(sizeof(makedata));
      return p;
    }
    ~makedata(){
      free(self);
    }
  }
} makedata;


void read_mak_input(makedata *mdata);
/* Get information for makedata from stdin.      */

void read_mak_file(char *basefilenm, makedata *mdata);
/* Read the data for makedata from the makefile. */

void write_mak_file(makedata *mdata);
/* Write the data for makedata to the makefile.  */

