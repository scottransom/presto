/* Length of the header in bytes */
#define MAX_WAPP_HEADER_SIZE 4096
/* Maximum number of samples to process at a time */
#define WAPP_MAXPTSPERBLOCK 64
/* Maximum number of WAPPs we can handle at once */
#define WAPP_MAXNUMWAPPS 4
/* time between correlator dumps in us */
#define WAPP_DEADTIME 0.34
/* Maximum number of lags we can have per WAPP */
#define WAPP_MAXLAGS 1024
/* Maximum data block length in bytes (the 2 is for 2IFs) */
#define WAPP_MAXDATLEN WAPP_MAXPTSPERBLOCK*WAPP_MAXLAGS*WAPP_MAXNUMWAPPS*2
/* Maximum block length in bytes for the raw lags */
#define WAPP_MAXLAGLEN WAPP_MAXDATLEN*4

typedef struct WAPP_HEADERv1{
  long header_version; /* some integer that increments with each revision */
  long header_size;    /* size (in bytes) of this header (nom = 2048) */
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
  long num_lags;       /* user-requested number of lags per dump per spect */
  long scan_number;    /* built by WAPP from year+daynumber+3-digit-number */
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
  long header_version; /* some integer that increments with each revision */
  long header_size;    /* size (in bytes) of this header (nom = 2048) */
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
  long num_lags;       /* user-requested number of lags per dump per spect */
  long scan_number;    /* built by WAPP from year+daynumber+3-digit-number */
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

/* wapp.c */
void set_WAPP_padvals(float *fpadvals, int good_padvals);
void print_WAPP_hdr(char *hdr);
int read_WAPP_rawblock(FILE *infiles[], int numfiles, unsigned char *data, 
		       int *padding, IFs ifs);
int read_WAPP_rawblocks(FILE *infiles[], int numfiles, unsigned char rawdata[], 
			int numblocks, int *padding, IFs ifs);
int read_WAPP(FILE *infiles[], int numfiles, float *data, int numpts, 
	      double *dispdelays, int *padding, int *maskchans, 
	      int *nummasked, mask *obsmask, IFs ifs);
void get_WAPP_channel(int channum, float chandat[], 
		      unsigned char rawdata[], int numblocks);
void get_WAPP_file_info(FILE *files[], int numwapps, int numfiles, int usewindow,
			float clipsig, long long *N, int *ptsperblock, 
			int *numchan, double *dt, double *T, 
			infodata *idata, int output);
int prep_WAPP_subbands(unsigned char *rawdata, float *data, 
		       double *dispdelays, int numsubbands, 
		       int transpose, int *maskchans, 
		       int *nummasked, mask *obsmask);
int read_WAPP_subbands(FILE *infiles[], int numfiles, float *data, 
		       double *dispdelays, int numsubbands, int transpose, 
		       int *padding, int *maskchans, int *nummasked, 
		       mask *obsmask, IFs ifs);
void convert_WAPP_point(void *rawdata, unsigned char *bytes, IFs ifs);
void WAPP_update_infodata(int numfiles, infodata *idata);
int skip_to_WAPP_rec(FILE *infiles[], int numfiles, int rec);
int check_WAPP_byteswap(char *hdr);
