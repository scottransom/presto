#include "wapp_key.h"

/* Length of the header in bytes */
#define MAX_WAPP_HEADER_SIZE 4096
/* Maximum number of samples to process at a time */
#define WAPP_MAXPTSPERBLOCK 64
/* Maximum number of WAPPs we can handle at once */
#define WAPP_MAXNUMWAPPS 8
/* time between correlator dumps in us */
#define WAPP_DEADTIME 0.34
/* Maximum number of lags we can have per WAPP */
#define WAPP_MAXLAGS 1024
/* Maximum data block length in bytes (the 2 is for 2IFs) */
#define WAPP_MAXDATLEN WAPP_MAXPTSPERBLOCK*WAPP_MAXLAGS*WAPP_MAXNUMWAPPS*2
/* Maximum block length in bytes for the raw lags */
#define WAPP_MAXLAGLEN WAPP_MAXDATLEN*4

typedef struct WAPP_HEADERv1{
  int header_version;  /* some integer that increments with each revision */
  int header_size;     /* size (in bytes) of this header (nom = 2048) */
  /*
    The following are obtained from current telescope status display
    note that start AST/LST are for reference purposes only and should
    not be taken as accurate time stamps. The time stamp can be derived
    from the obs_date/start_time variables further down in the structure.
  */
  double src_ra;       /* requested ra J2000 (10000*hr+100*min+sec) */
  double src_dec;      /* requested dec J2000 (10000*deg+100*min+sec) */
  double start_az;     /* telescope azimuth at start of scan (deg) */
  double start_za;     /* telescope zenith angle at start of scan (deg) */
  double start_ast;    /* AST at start of scan (sec) */
  double start_lst;    /* local siderial time at start of scan (sec) */
  /*
    In the following, anything not supplied/requested by the user
    is assumed to be calculated by WAPP when it writes the header
  */
  double cent_freq;    /* user-supplied band center frequency (MHz) */
  double obs_time;     /* user-requested length of this integration (s) */
  double samp_time;    /* user-requested sample time (us) */
  double wapp_time;    /* actual sample time (us) i.e. requested+dead time */
  double bandwidth;    /* total bandwidth (MHz) for this observation */
  int num_lags;        /* user-requested number of lags per dump per spect */
  int scan_number;     /* built by WAPP from year+daynumber+3-digit-number */
  char src_name[24];   /* user-supplied source name (usually pulsar name) */
  char obs_date[24];   /* Date of OBS in format yyyymmdd */
  char start_time[24]; /* UT of OBS start in hh:mm:ss (start on 1-sec tick) */
  char project_id[24]; /* user-supplied AO proposal number (XYYYY) */
  char observers[24];  /* observer(s) name(s) */
  int nifs;            /* user-requested: number of IFs to be recorded */
  int level;           /* user-requested: 1 means 3-level; 2 mean 9-level */
  int sum;             /* user-requested: 1 means that data is sum of IFs */
  int freqinversion;   /* 1 band is inverted, else band is not inverted */
  long long timeoff;   /* number of reads between obs start and snap block */
  int lagformat;       /* 0=16 bit unsigned ints, 1=32 bit unsigned ints */
  int lagtrunc;        /* (reserved) number bits truncated (0 if no trunc) */
  double power_analog[2];   /* Power measured by Analog Detector */
} WAPP_HEADERv1;

typedef struct WAPP_HEADERv2plus{
  int header_version;  /* some integer that increments with each revision */
  int header_size;     /* size (in bytes) of this header (nom = 2048) */
  char obs_type[24];   /* what kind of observation is this: */
                       /* PULSAR_SEARCH, PULSAR_FOLDING, SPECTRA_TOTALPOWER */
  /*
    The following are obtained from current telescope status display
    note that start AST/LST are for reference purposes only and should
    not be taken as accurate time stamps. The time stamp can be derived
    from the obs_date/start_time variables further down in the structure.
  */
  double src_ra;       /* requested ra J2000 (10000*hr+100*min+sec) */
  double src_dec;      /* requested dec J2000 (10000*deg+100*min+sec) */
  double start_az;     /* telescope azimuth at start of scan (deg) */
  double start_za;     /* telescope zenith angle at start of scan (deg) */
  double start_ast;    /* AST at start of scan (sec) */
  double start_lst;    /* local siderial time at start of scan (sec) */
  /*
    In the following, anything not supplied/requested by the user
    is assumed to be calculated by WAPP when it writes the header
  */
  double cent_freq;    /* user-supplied band center frequency (MHz) */
  double obs_time;     /* user-requested length of this integration (s) */
  double samp_time;    /* user-requested sample time (us) */
  double wapp_time;    /* actual sample time (us) i.e. requested+dead time */
  double bandwidth;    /* total bandwidth (MHz) for this observation */
  int num_lags;        /* user-requested number of lags per dump per spect */
  int scan_number;     /* built by WAPP from year+daynumber+3-digit-number */
  char src_name[24];   /* user-supplied source name (usually pulsar name) */
  char obs_date[24];   /* Date of OBS in format yyyymmdd */
  char start_time[24]; /* UT of OBS start in hh:mm:ss (start on 1-sec tick) */
  char project_id[24]; /* user-supplied AO proposal number (XYYYY) */
  char observers[24];  /* observer(s) name(s) */
  int nifs;            /* user-requested: number of IFs to be recorded */
  int level;           /* user-requested: 1 means 3-level; 2 mean 9-level */
  int sum;             /* user-requested: 1 means that data is sum of IFs */
  int freqinversion;   /* 1/0 whether the WAPP inverted the band or not */
  long long timeoff;   /* number of reads between obs start and snap block */
  int lagformat;       /* 0=16 bit uint lags , 1=32 bit uint lags */
                       /* 2=32 bit float lags, 3=32 bit float spectra */
  int lagtrunc;        /* if we truncate data (0 no trunc) */
                       /* for 16 bit lagmux modes, selects which 16 bits */
                       /* of the 32 are included as data */
                       /* 0 is bits 15-0 1,16-1 2,17-2...7,22-7 */
  int firstchannel;    /* 0 when correlator channel a is first, 1 if b */
  int nbins;           /* number of time bins for pulsar folding mode */
                       /* doulbles as maxrecs for snap mode */
  int isfolding;       /* is folding selected                              */
  int isalfa;          /* is ALFA selected                                 */
  double dumptime;     /* folded integrations for this period of time */
  double power_analog[2];   /* Power measured by Analog Detector */
  char padding[2020];  /* Padding that ignores the timing mode params */
  int iflo_flip;       /* consider entire iflo and determine upper/lower sb */
} WAPP_HEADERv2plus;

typedef struct WAPP_HEADERv9{
  int header_version;  /* some integer that increments with each revision */
  int header_size;     /* size (in bytes) of this header (nom = 2048) */
  char obs_type[24];   /* what kind of observation is this: */
                       /* PULSAR_SEARCH, PULSAR_FOLDING, SPECTRA_TOTALPOWER */
  /*
    The following are obtained from current telescope status display
    note that start AST/LST are for reference purposes only and should
    not be taken as accurate time stamps. The time stamp can be derived
    from the obs_date/start_time variables further down in the structure.
  */
  double src_ra;       /* requested ra J2000 (10000*hr+100*min+sec) */
  double src_dec;      /* requested dec J2000 (10000*deg+100*min+sec) */
  double start_az;     /* telescope azimuth at start of scan (deg) */
  double start_za;     /* telescope zenith angle at start of scan (deg) */
  double start_ast;    /* AST at start of scan (sec) */
  double start_lst;    /* local siderial time at start of scan (sec) */
  /*
    In the following, anything not supplied/requested by the user
    is assumed to be calculated by WAPP when it writes the header
  */
  double cent_freq;    /* user-supplied band center frequency (MHz) */
  double obs_time;     /* user-requested length of this integration (s) */
  double samp_time;    /* user-requested sample time (us) */
  double wapp_time;    /* actual sample time (us) i.e. requested+dead time */
  double bandwidth;    /* total bandwidth (MHz) for this observation */
  int num_lags;        /* user-requested number of lags per dump per spect */
  int scan_number;     /* built by WAPP from year+daynumber+3-digit-number */
  char src_name[24];   /* user-supplied source name (usually pulsar name) */
  char obs_date[24];   /* Date of OBS in format yyyymmdd */
  char start_time[24]; /* UT of OBS start in hh:mm:ss (start on 1-sec tick) */
  char project_id[24]; /* user-supplied AO proposal number (XYYYY) */
  char observers[24];  /* observer(s) name(s) */
  char frontend[24];   /* receiver abbrev name */
  char backendmode[24];/* backend mode description */
  char caltype[8];     /* calibrator type */
  int nifs;            /* user-requested: number of IFs to be recorded */
  int level;           /* user-requested: 1 means 3-level; 2 mean 9-level */
  int sum;             /* user-requested: 1 means that data is sum of IFs */
  int freqinversion;   /* 1/0 whether the WAPP inverted the band or not */
  long long timeoff;   /* number of reads between obs start and snap block */
  int lagformat;       /* 0=16 bit uint lags , 1=32 bit uint lags */
                       /* 2=32 bit float lags, 3=32 bit float spectra */
  int lagtrunc;        /* if we truncate data (0 no trunc) */
                       /* for 16 bit lagmux modes, selects which 16 bits */
                       /* of the 32 are included as data */
                       /* 0 is bits 15-0 1,16-1 2,17-2...7,22-7 */
  int firstchannel;    /* 0 when correlator channel a is first, 1 if b */
  int nbins;           /* number of time bins for pulsar folding mode */
                       /* doulbles as maxrecs for snap mode */
  int isfolding;       /* is folding selected */
  int isalfa;          /* is ALFA selected */
  int isdual;          /* are WAPPs in dual board mode */
  int fold_bits;       /* 0 if 16 bits (old default) 1 if 32 bit folding */
  int iflo_flip[2];    /* is signal flipped in each board */
  int attena;          /* first board parallel port value */
  int attenb;          /* second board parallel port value */
  double dumptime;     /* folded integrations for this period of time */
  double power_analog[2];   /* power measured by analog detector */
/*
    In the following, pulsar-specific information is recorded for use
    by folding programs e.g. the quick-look software. This is passed to
    WAPP by psrcontrol at the start of the observation.

    The apparent pulse phase and frequency at time "dt" minutes with
    respect to the start of the observation are then calculated as:

    phase = rphase + dt*60*f0 + coeff[0] + dt*coeff[1] + dt*dt*coeff[2] + ...
    freq(Hz) = f0 + (1/60)*(coeff[1] + 2*dt*coeff[2] + 3*dt*dt*coeff[3] + ...)

    where the C notation has been used (i.e. coeff[0] is first coefficient etc)
    for details, see TEMPO notes (http://www.naic.edu/~pulsar/docs/tempo.txt)
*/
  double psr_dm;        /* pulsar's dispersion measure (cm-3 pc) */
  double rphase[16];    /* reference phase of pulse (0-1) */
  double psr_f0[16];    /* pulse frequency at reference epoch (Hz) */
  double poly_tmid[16]; /* mid point of polyco in (MJD) modified Julian date */
  double coeff[192];    /* polynomial coefs made by TEMPO, 16 sets of 12 */
  int num_coeffs[16];   /* number of coefficients */
  char hostname[24];    /* ascii hostname of machine that took this data */
  /* ALFA info */
  double rfeed_offaz[7];/* deg az rotated offset all alfa beam to center */
  double rfeed_offza[7];/* deg za rotated offset all alfa beam to center */
  double prfeed_offaz;  /* deg az offset to center of alfa beam */
  double prfeed_offza;  /* deg za offset to center of alfa beam */
  double alfa_raj[7];   /* hr starting actual ra position of alfa beams */
  double alfa_decj[7];  /* deg starting actual dec position of alfa beams */
  double alfa_az[7];    /* deg starting actual az position of alfa beams */
  double alfa_za[7];    /* deg starting actual za postion of alfa beams */
  double alfa_ang;      /* deg alfa rotation angle */
  double para_ang;      /* deg paralactic angle of center beam */
  /* add IF/LO data */
  double syn1;          /* upstairs synthesizer freq Hz */
  double synfrq[8];     /* downstairs synthesizers (Hz) */
  unsigned char prfeed;     /* centered alfa beam */
  unsigned char shcl;       /* true if receiver shutter closed */
  unsigned char sbshcl;     /* true if Sband receiver shutter closed */
  unsigned char rfnum;      /* position of the receiver selectror */
  unsigned char zmnormal;   /* transfer switch to reverse channels, true normal */
  unsigned char rfattn[2];  /* attenuator position */
  unsigned char ifnum;      /* if selector, 1/300 2/750, 3/1500, 4/10GHz1500, 5-thru */
  unsigned char ifattn[2];  /* IF attenuator positions */
  unsigned char fiber;      /* true if fiber is chosen (always the case) */
  unsigned char ac2sw;      /* ac power to various instruments and other stuff */
  unsigned char if750nb;    /* narrow band 750 filter selected */
  unsigned char phbsig;     /* converter combiner signal phase adjust */
  unsigned char hybrid;     /* converter combiner hybrid */
  unsigned char phblo;      /* convert combiner lo phase adjust */
  /* downstairs */
  unsigned char xfnormal;   /* transfer switch true if normal */
  unsigned char noise;      /* noise on */
  unsigned char ampgain[2]; /* gain of downstairs amps */
  unsigned char inpfrq;     /* input distributor position */
  unsigned char mixer[8];   /* mixer source switches */
  unsigned char vlbainp;    /* vlba input switch position */
  unsigned char syndest[8]; /* synthesizer destinations */
  unsigned char calsrc;     /* cal source bit */
  unsigned char vis30mhz;   /* greg 1 ch 0 */
  unsigned char pwrmet;     /* power meter input switch */
  unsigned char blank430;   /* 430 blanking on */
  unsigned char fill[6];    /* fill */
} WAPP_HEADERv9;

/* wapp.c */
void print_WAPP_hdr(struct HEADERP *hdr);
void get_WAPP_file_info(FILE *files[], int numwapps, int numfiles, int usewindow,
			float clipsig, long long *N, int *ptsperblock, 
			int *numchan, double *dt, double *T, 
			infodata *idata, int output);
void WAPP_update_infodata(int numfiles, infodata *idata);
