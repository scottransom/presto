typedef struct SIGPROCFB {
  char inpfile[80];      /* Input filename */
  char source_name[80];  /* Source name */
  char ifstream[8];      /* Y=IF present, X=IF not present (4 possible values) */
  double tstart;         /* MJD start time */
  double tsamp;          /* Sampling time in sec */
  double src_raj;        /* Source RA  (J2000) in hhmmss.ss */
  double src_dej;        /* Source DEC (J2000) in ddmmss.ss */
  double az_start;       /* Starting azimuth in deg */
  double za_start;       /* Starting zenith angle in deg */
  double fch1;           /* Highest channel frequency (MHz) */
  double foff;           /* Channel stepsize (MHz) */
  int machine_id;        /* Instrument ID (see backend_name() */
  int telescope_id;      /* Telescope ID (see telescope_name() */
  int nchans;            /* Number of finterbank channels */
  int obits;             /* Number of bits in the filterbank samples */
  int nifs;              /* Number if IFs present */
  int sumifs;            /* Whether the IFs are summed or not */
  int headerlen;         /* Header length in bytes */
} sigprocfb;
