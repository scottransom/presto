/* Length of the header in bytes */
#define WAPP_HEADER_SIZE 2048
/* Number of samples to process at a time */
#define WAPP_PTSPERBLOCK 64
/* time between correlator dumps in us */
#define WAPP_DEADTIME 0.34
/* Maximum number of lags we can have */
#define WAPP_MAXLAGS 1024
/* Maximum data block length in bytes */
#define WAPP_MAXDATLEN WAPP_PTSPERBLOCK*WAPP_MAXLAGS
/* Maximum lag block length in bytes */
#define WAPP_MAXLAGLEN WAPP_MAXDATLEN*4

typedef struct WAPP_HEADER{
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
  /*
    In the following, pulsar-specific information is recorded for use
    by folding programs e.g. the quick-look software. This is passed to
    WAPP by psrcontrol at the start of the observation.
    
    The apparent pulse phase and frequency at time "dt" minutes with
    respect to the start of the observation are then calculated as:
    
    phase=rphase+dt*60*f0+coeff[0]+dt*coeff[1]+dt*dt*coeff[2]+...
    freq(Hz)=f0+(1/60)*(coeff[1]+2*dt*coeff[2]+3*dt*dt*coeff[3]+...)
 
    where C notation has been used (i.e. coeff[0] is first coeff etc)
    For details, see TEMPO notes (http://www.naic.edu/~pulsar/docs/tempo.txt)
  */
  double psr_dm;       /* pulsar's dispersion measure (cm-3 pc) */
  double rphase[9];    /* reference phase of pulse (0-1) */
  double psr_f0[9];    /* pulse frequency at reference epoch (Hz) */
  double poly_tmid[9]; /* mid point of polyco in (MJD) modified Julian date */
  double coeff[144];   /* polycos calculated by TEMPO effectively [9][16] */
  int num_coeffs[9];   /* number of coefficients */
  char filler[364];    /* this pads out the header to 2048 bytes */
} WAPP_HEADER;

/* wapp.c */
void WAPP_hdr_to_inf(WAPP_HEADER *hdr, infodata *idata);
void get_WAPP_file_info(FILE *files[], int numfiles, float clipsig, 
			long long *N, int *ptsperblock, int *numchan, 
			double *dt, double *T, infodata *idata, int output);
void WAPP_update_infodata(int numfiles, infodata *idata);
int skip_to_WAPP_rec(FILE *infiles[], int numfiles, int rec);
void print_WAPP_hdr(WAPP_HEADER *hdr);
int read_WAPP_rawblock(FILE *infiles[], int numfiles, 
		       unsigned char *data, int *padding);
int read_WAPP_rawblocks(FILE *infiles[], int numfiles, 
			unsigned char rawdata[], int numblocks, 
			int *padding);
int read_WAPP(FILE *infiles[], int numfiles, float *data, int numpts, 
	      double *dispdelays, int *padding, int *maskchans, 
	      int *nummasked, mask *obsmask);
void get_WAPP_channel(int channum, float chandat[], 
		      unsigned char rawdata[], int numblocks);
int read_WAPP_subbands(FILE *infiles[], int numfiles, float *data, 
		       double *dispdelays, int numsubbands, int transpose, 
		       int *padding, int *maskchans, int *nummasked, 
		       mask *obsmask);
void convert_WAPP_point(void *rawdata, unsigned char *bytes);
int check_WAPP_byteswap(WAPP_HEADER *hdr);
