#define RECLEN 49792
#define DATLEN 49152
#define HDRLEN 640
#define SAMPPERBLK DATLEN*8

/* The following is from Manchester, 9 Feb 98 */

/* Structure defining the tape header 

 * $Log: multibeam.h,v $
 * Revision 1.9  2000/12/17 06:00:53  ransom
 * Added mask.c and mask.h
 * Re-ran CLIG after modifying parameters in a fwe CLIG files.
 * "Finished" the reworking of the PKMB code to accept multiple
 *    files.  Compiles cleanly.  Now must test...
 * Added lots of stuff pertaining to mask determination.  Still
 *    need to actually compute the mask, though.  Also need to
 *    add the masking routines to the rawdata extracting files.
 *
 * Revision 1.8  2000/12/15 22:17:24  ransom
 * Added MJD difference routine to misc_utils.
 *
 * Revision 1.7  2000/12/15 04:47:39  ransom
 * Re-ran CLIG after modifying rfifind_cmd.cli for multiple files and
 *    small clean-ups.
 * Added a function that determines basic info for PKMB files that
 *    need to be combined (in multibeam.c).
 * Began modifying rfifind for multiple input files.
 *
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

typedef struct PKMB_TAPEHDR {

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
} PKMB_tapehdr;


/* Routines to read and convert Parkes Multibeam data files */

void get_PKMB_file_info(FILE *files[], int numfiles, long long *N, 
			int *ptsperblock, int *numchan, double *dt, 
			double *T, int output);
/* Read basic information into static variables and make padding      */
/* calculations for a set of PKMB rawfiles that you want to patch     */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */

int skip_to_PKMB_rec(FILE * infile, int rec);
/* This routine skips to the record 'rec' in the input file */
/* *infile.  *infile contains 1 bit digitized data from the */
/* PKMB backend at Parkes.  Returns the record skipped to.  */

int read_PKMB_rawblock(FILE *infiles[], int numfiles, 
		       PKMB_tapehdr *hdr, unsigned char *data);
/* This routine reads a single record from the         */
/* input files *infiles which contain 1 bit digitized  */
/* data from the PKMB pulsar backend at Parkes.        */
/* Length of a PKMB record is 640 bytes for the header */
/* plus 48k of data = 49792 bytes.                     */
/* The header of the record read is placed in hdr.     */
/* *data must be pre-allocated with a size of 48k.     */

int read_PKMB_rawblocks(FILE *infiles[], int numfiles, 
			unsigned char rawdata[], int numblocks);
/* This routine reads numblocks PKMB records from the input */
/* files *infiles.  The raw bit data is returned in rawdata */
/* which must have a size of numblocks*DATLEN.  The number  */
/* of blocks read is returned.                              */

int read_PKMB(FILE *infiles[], int numfiles, float *data, 
	      double *dispdelays);
/* This routine reads a PKMB record from the input     */
/* files *infiles.  These files contain 1 bit data     */
/* from the PKMB backend at Parkes.  Time delays and   */
/* and a mask are applied to each channel.  It returns */
/* the # of points read if succesful, 0 otherwise.     */

void get_PKMB_channel(int channum, float chandat[], 
		      unsigned char rawdata[], int numblocks);
/* Return the values for channel 'channum' of a block of       */
/* 'numblocks' raw PKMB data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                */
/* read_PKMB_rawblocks(), and 'chandat' must have at least     */
/* 'numblocks' * 'ptsperblk_st' spaces.                        */
/* Channel 0 is assumed to be the lowest freq channel.         */

int read_PKMB_subbands(FILE *infiles[], int numfiles, 
		       float *data, double *dispdelays, int numsubbands);
/* This routine reads a record from the input files *infiles[]   */
/* which contain data from the PKMB system.  The routine uses    */
/* dispersion delays in 'dispdelays' to de-disperse the data     */
/* into 'numsubbands' subbands.  It stores the resulting data    */
/* in vector 'data' of length 'numsubbands' * 'ptsperblk_st'.    */
/* The low freq subband is stored first, then the next highest   */
/* subband etc, with 'ptsperblk_st' floating points per subband. */
/* It returns the # of points read if succesful, 0 otherwise.    */

void PKMB_hdr_to_inf(PKMB_tapehdr * hdr, infodata * idata);
/* Convert PKMB header into an infodata structure */

void print_PKMB_hdr(PKMB_tapehdr * hdr);
/* Output a PKMB header in human readable form */

void convert_PKMB_point(unsigned char *bits, unsigned char *bytes);
/* This routine converts 1 bit digitized data */
/* into an array of 'numchan' bytes.          */
