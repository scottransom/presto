#include "backend_common.h"

typedef struct SIGPROCFB {
  char inpfile[80];      /* Input filename */
  char source_name[80];  /* Source name */
  char ifstream[8];      /* Y=IF present, X=IF not present (4 possible values) */
  long long N;           /* Number of points (in time) in the file */
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
  int nbits;             /* Number of bits in the filterbank samples */
  int nifs;              /* Number of IFs present */
  int nbeams;            /* Number of beams in the observing system */
  int ibeam;             /* Beam number used for this data */
  int sumifs;            /* Whether the IFs are summed or not */
  int signedints;        /* Whether the integer data is signed or not */
  int headerlen;         /* Header length in bytes */
} sigprocfb;

/* sigproc_fb.c */
void get_telescope_name(int telescope_id, struct spectra_info *s);
void get_backend_name(int machine_id, struct spectra_info *s);
void write_filterbank_header(sigprocfb *fb, FILE *outfile);
int read_filterbank_header(sigprocfb *fb, FILE *inputfile);
void read_filterbank_files(struct spectra_info *s);
long long offset_to_filterbank_spectra(long long specnum, struct spectra_info *s);
int get_filterbank_rawblock(float *fdata, struct spectra_info *s, int *padding);
void convert_filterbank_block(float *outdata, unsigned char *indata, int numread, struct spectra_info *s);
