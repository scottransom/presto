#define RECLEN 49792
#define DATLEN 49152
#define HDRLEN 640

/* The following is from Manchester, 9 Feb 98 */

/* Structure defining the tape header 

 * $Log: multibeam.h,v $
 * Revision 1.6  2000/12/06 13:41:37  ransom
 * Added new function to multibeam.h
 * Added fftcand structure in presto.h
 * Modified calling sequence to power/sigma routines.
 *
 * Revision 1.5  1999/12/02 20:26:49  ransom
 * Made a multibeam subband de-dispersion routine.
 *
 * Revision 1.4  1999/11/16 19:43:51  ransom
 * Changed a few calling parameters.
 *
 * Revision 1.3  1999/11/16 03:51:57  ransom
 * Simplified interface and code in dedisp.c
 * Minor modifications to multibeam code (mostly asthetic)
 *
 * Revision 1.2  1999/07/12 20:42:20  ransom
 * Changed convert_multibeam_point() to accept a variable number of channels.
 *
 * Revision 1.1.1.1  1999/07/02 02:04:45  ransom
 * Imported PRESTO sources
 *
 * Revision 1.1.1.1  1997/08/20 05:07:11  rmanches
 * CVS import
 *
 *
 * Revision 1.1  1997/04/23 04:37:03  rmanches
 * Initial revision
 *
 */

typedef struct MULTIBEAM_TAPEHDR {

/*
 *  Length of full header is 640 bytes (all char)
 */

  char program[6];      /* Program name */
  char version[4];      /* Version number "vv.v" */
  char tape_lbl[6];     /* Tape label */
  char file_cntr[4];    /* File number on tape */
  char blk_cntr[8];     /* Block counter per beam in file */
  char date[8];         /* UT Date "yymmdd  " */
  char mjd[8];          /* MJD at 00h UT */
  char ut_start[16];    /* UT at file start  " hh:mm:ss.ssssss" */
  char lst_start[12];   /* LST at file start " hh:mm:ss.ss" */
  char coord_type[2];   /* Coord system " 4" = Galactic, " 5" = J2000 */
  char ra_start[16];    /* RA at file start  " hh:mm:ss.ssss" */
  char dec_start[16];   /* Dec at file start "-dd:mm:ss.sss " */
  char gl_start[8];     /* Galactic long at file start (deg) "ddd.dddd" */
  char gb_start[8];     /* Galactic lat. at file start (deg) "-dd.dddd" */
  char fangle[8];       /* FA, PA or GPA (deg) "-ddd.ddd" */
  char obs_time[8];     /* Obs length "ssss.sss"  */
  char comment[64];     /* Comment */
  char nfilter[2];      /* Number of filter systems, " 1" or " 2" */
  char chanbw[2][8];    /* Channel incrt (MHz, -ve if inverted) "-b.bbbbb" */
  char nchan[2][4];     /* Nr of channels in each filter system */
  char freq_ch1[2][12]; /* RF of first channel centre "ffffff.fffff" */
  char samp_int[2][12]; /* Samp intval in ms (1 = 2 always) "mm.mmmmmmmmm" */
  char samp_grp[2][4];  /* Samp per group (For backwards compat., now = 1) */
  char samp_blk[8];     /* Samples per block (Was groups per block) */
  char blk_sec[8];      /* Seconds per tape block per beam "ss.sssss" */
  char fdctrl[2];       /* 0 = none, 1 = fixed FA, 2 = fixed PA or GPA */
  char data[2][1];      /* Data type code ("2"=pol, "3"=normal, "5"=dedisp) */
  char ut_blk[16];      /* UT at blk start " hh:mm:ss.ssssss" */
  char lst_blk[12];     /* LST at blk start " hh:mm:ss.ss" */
  char ra_blk[16];      /* RA at blk start " hh:mm:ss.ssss" */
  char dec_blk[16];     /* Dec at blk start "-dd:mm:ss.sss" */
  char gl_blk[8];       /* Gal long. at blk start "ddd.dddd" */
  char gb_blk[8];       /* Gal lat at blk start "-dd.dddd" */
  char zen_blk[8];      /* Zenith angle at blk start (deg) "ddd.dddd" */
  char az_blk[8];       /* Azimuth angle at blk start (deg) "ddd.dddd" */
  char atten[4][4];     /* Attenuator settings (A,B) for Sys 1,2 "dd.d" */
  char tpower[20][4];   /* Total powers for 32-channel groups "iiii"  */
  char nblk_read[8];    /* Nr of tape blocks in disk file */
  char scan_ratel[8];   /* Scan rate in `long' (deg/min) "-r.rrrrr" */
  char scan_rateb[8];   /* Scan rate in `lat'  (deg/min) "-r.rrrrr" */
  char nbeam[4];        /* Total number of beams */
  char ibeam[4];        /* Beam number */
  char pname[16];       /* Psrname (B or J) or pointing name */
  char cfg_file[8];     /* Config file name */
  char schd_file[8];    /* Schedule file name */
  char nbits[2][2];     /* Nr of bits per sample */
  char ddm[8];          /* DM for dedispersion "ddddd.dd" */
  char nddch[2][4];     /* Nr of channels per dedispersed band */
  char move[2];         /* Tel status "k0": On posn k(0-9), "01": On grid posn,
                           "02": Off position, "03": GC scan, "04" SC scan */
  char pnterr[6];       /* Pointing error (arcmin) "mmm.mm" */
  char tree[2];         /* "  ": Normal, " T": Tree dedisp, " D": PDM dedisp */
  char ifilter[2];      /* Filter system for disk file */
  char telid[10];       /* Telescope ID (as in tel.def) */
  char pangle[8];       /* Parallactic angle (deg) "-ddd.ddd" */
  char bsw[8];          /* Nr of samples per beam switch */
  char cal_cycle[4];    /* Cal cycle period in samples */
  char spare[22];       /* Spare */
} multibeam_tapehdr;


/* Routines to read and convert Parkes Multibeam data files */

void convert_multibeam_point(unsigned char *rec, unsigned char *data,
			     int numchan, int decreasing_f);
/* This routine converts 1 bit digitized data with 'numchan' */
/* channels to an array of 'numchan' floats.                 */

void multibeam_hdr_to_inf(multibeam_tapehdr *hdr, infodata *idata);
/* Convert appropriate Multibeam header portions to make */
/* a basic '.inf' file                                   */

int skip_to_multibeam_rec(FILE * infile, int rec);
/* This routine skips to the record 'rec' in the input file */
/* *infile.  *infile contains 1 bit digitized data from the */
/* multibeam receiver at Parkes                             */
/* Returns the record that was skipped to.                  */

int read_rawmultibeam(FILE *infile, float *data, int numchan, 
		      int numblocks);
/* This routine reads numblocks PKMB records with numchan */
/* channels each from the input file *infile.  The number */
/* of blocks read is returned.                            */

int read_multibeam_recs(FILE * infile, multibeam_tapehdr * hdr,
			unsigned char *data, int blocks_to_read);
/* This routine reads blocks_to_read multibeam records from   */
/* the input file *infile which contains 1 bit digitized data */
/* from the multibeam correlator at Parkes.                   */
/* Length of a multibeam record is 640 bytes for the header   */
/* plus 48k of data = 49792 bytes.                            */
/* The header of the first record read is placed in hdr.      */
/* *data must be pre-allocated with size 48k * blocks_to_read */

int read_multibeam(FILE * file, float *data, int numpts,
		   double *dispdelays, int numchan);
/* This routine reads a numpts record with numchan each from */
/* the input file *file which contains 1 bit digitized data  */
/* from the multibeam correlator at Parkes.                  */
/* It returns the number of points read.                     */

int read_multibeam_subbands(FILE * infile, float *data, int numpts,
			    double *dispdelays, int numsubbands, 
			    int numchan);
/* This routine reads a numpts record with numchan each from     */
/* the input file *infile which contains 1 bit digitized data    */
/* from the multibeam correlator at Parkes.  The routine uses    */
/* dispersion delays in 'dispdelays' to de-disperse the data     */
/* into 'numsubbands' subbands.  It stores the resulting         */
/* data in vector 'data' of length numsubbands * numpts.  The    */
/* low frequency subband is stored first, then the next highest  */
/* subband etc, with 'numpts' floating points per subband.       */
/* It returns the number of points (pts_per_read) read if        */
/* succesful or 0 if unsuccessful.                               */

void print_multibeam_hdr(multibeam_tapehdr *hdr);
/* Output in human readable form a multibeam header. */
