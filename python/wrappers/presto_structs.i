
%inline %{
typedef struct {
     double *dptr;
} DoubleArray;
%}

%typemap(memberin) char *name {
        strncpy($1,$input,200);
}       
%typemap(memberin) char *object {
        strncpy($1,$input,100);
}       
%typemap(memberin) char *telescope {
        strncpy($1,$input,40);
}       
%typemap(memberin) char *instrument {
        strncpy($1,$input,100);
}       
%typemap(memberin) char *observer {
        strncpy($1,$input,100);
}
%typemap(memberin) char *band {
        strncpy($1,$input,40);
}       
%typemap(memberin) char *filt {
        strncpy($1,$input,7);
}       
%typemap(memberin) char *analyzer {
        strncpy($1,$input,100);
}       
%typemap(memberin) char *notes {
        strncpy($1,$input,500);
}       

%extend DoubleArray {
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
  d.dptr = $1;
  $result = &d;
}

// Set a DoubleArray
%typemap(memberin) DoubleArray *onoff {
  $1 = $input->dptr;
}

typedef struct infodata {
  double ra_s;		/* Right ascension seconds (J2000)       */
  double dec_s;		/* Declination seconds (J2000)           */ 
  double N;		/* Number of bins in the time series     */
  double dt;	 	/* Width of each time series bin (sec)   */
  double fov;		/* Diameter of Beam or FOV in arcsec     */
  double mjd_f;		/* Epoch of observation (MJD) frac part  */
  double dm;		/* Radio -- Dispersion Measure (cm-3 pc) */
  double freq;		/* Radio -- Low chan central freq (Mhz)  */
  double freqband;	/* Radio -- Total Bandwidth (Mhz)        */
  double chan_wid;	/* Radio -- Channel Bandwidth (Mhz)      */
  double wavelen;	/* IR,Opt,UV -- central wavelength (nm)  */
  double waveband;	/* IR,Opt,UV -- bandpass (nm)            */
  double energy;	/* x-ray,gamma -- central energy (kev)   */
  double energyband;	/* x-ray,gamma -- energy bandpass (kev)  */
  double onoff[40];     /* Bin number pairs where obs is "on"    */
  int num_chan;		/* Radio -- Number Channels              */
  int mjd_i;		/* Epoch of observation (MJD) int part   */
  int ra_h;		/* Right ascension hours (J2000)         */
  int ra_m;		/* Right ascension minutes (J2000)       */
  int dec_d;		/* Declination degrees (J2000)           */
  int dec_m;		/* Declination minutes (J2000)           */  
  int bary;		/* Barycentered?  1=yes, 0=no            */
  int numonoff;		/* The number of onoff pairs in the data */ 
  char *notes;	        /* Any additional notes                  */
  char *name;	        /* Data file name without suffix         */
  char *object;	        /* Object being observed                 */ 
  char *instrument;	/* Instrument used                       */
  char *observer;	/* Observer[s] for the data set          */
  char *analyzer;	/* Who analyzed the data                 */
  char *telescope;	/* Telescope used                        */
  char *band;	        /* Type of observation (EM band)         */
  char *filt;		/* IR,Opt,UV -- Photometric Filter       */
} infodata;

void readinf(infodata * data, char *filenm);
/* Read a ".inf" file into memory. */

void writeinf(infodata *data);
/* Write a ".inf" file to disk. */


%typemap(memberin) char *basefilenm {
        strncpy($1,$input,200);
}       
%typemap(memberin) char *description {
        strncpy($1,$input,200);
}       
%typemap(memberin) char *ptype {
        strncpy($1,$input,20);
}       
%typemap(memberin) char *round {
        strncpy($1,$input,20);
}       
%typemap(memberin) char *noisetype {
        strncpy($1,$input,20);
}       

typedef struct makedata {
  char *basefilenm;	        /* Data file name without suffix         */
  char *description;	        /* Data description                      */
  long N;		        /* Number of bins in the time series     */
  long next2_to_n;              /* The next power-of-2 >= N              */
  double dt;	 	        /* Width of each time series bin (sec)   */
  double T;                     /* Total length of time series (sec)     */
  char *ptype;		        /* Pulsetype (Sine, Crab-like, Spike)    */
  int pnum;                     /* 1=Sine, 2=Crab-like, 3=Spike          */
  double fwhm;                  /* FWHM Phase (0-1) if ptype is Spike    */
  char *round;		        /* Rounding Format (Whole, Fractional)   */
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
  char *noisetype;		/* Noise type (Standard, Other)          */
  int noise;                    /* 1=Standard, 0=Other                   */
  double noisesig;              /* Noise standard deviation              */
  int numonoff;                 /* The number of onoff pairs in the data */
  DoubleArray *onoff;		/* Bin number pairs where obs is "on"    */
} makedata;


%typemap(memberin) char *jname {
        strncpy($1,$input,13);
}       

%typemap(memberin) char *bname {
        strncpy($1,$input,9);
}       

%typemap(memberin) char *alias {
        strncpy($1,$input,10);
}       

/* This is a structure that contains the 'key' pulsar info  */
typedef struct psrparams {
  char *jname;               /* The PSRs J2000 name         */
  char *bname;               /* The PSRs B1950 name         */
  char *alias;               /* An alias for the pulsar     */
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
  orbitparams orb;           /* Orbital parameters          */
} psrparams;

/*
%extend psrparams {
  psrparams *Psrparams(void) {
    psrparams *p;
    p = (psrparams *p) calloc(sizeof(psrparams), 1);
    return p;
  }
  ~Psrparams() {
    free(self);
  }
}
*/

void read_mak_input(makedata *mdata);
/* Get information for makedata from stdin.      */

void read_mak_file(char *basefilenm, makedata *mdata);
/* Read the data for makedata from the makefile. */

void write_mak_file(makedata *mdata);
/* Write the data for makedata to the makefile.  */

