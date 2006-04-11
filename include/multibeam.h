#define RECLEN 49792
#define DATLEN 49152
#define HDRLEN 640
#define SAMPPERBLK DATLEN*8

/* Structure defining the tape header */

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

void get_PKMB_file_info(FILE *files[], int numfiles, float clipsig, 
                        long long *N, int *ptsperblock, int *numchan, 
                        double *dt, double *T, int output);
void PKMB_update_infodata(int numfiles, infodata *idata);
int skip_to_PKMB_rec(FILE *infiles[], int numfiles, int rec);
int read_PKMB_rawblock(FILE *infiles[], int numfiles, 
		       PKMB_tapehdr *hdr, unsigned char *data,
		       int *padding);
int read_PKMB_rawblocks(FILE *infiles[], int numfiles, 
			unsigned char rawdata[], int numblocks,
			int *padding);
int read_PKMB(FILE *infiles[], int numfiles, float *data, 
	      int numpts, double *dispdelays, int *padding,
	      int *maskchans, int *nummasked, mask *obsmask);
void get_PKMB_channel(int channum, float chandat[], 
		      unsigned char rawdata[], int numblocks);
int prep_PKMB_subbands(unsigned char *rawdata, float *data, 
		       double *dispdelays, int numsubbands, 
		       int transpose, int *maskchans, 
		       int *nummasked, mask *obsmask);
int read_PKMB_subbands(FILE *infiles[], int numfiles, float *data, 
		       double *dispdelays, int numsubbands, 
		       int transpose, int *padding, 
		       int *maskchans, int *nummasked, mask *obsmask);
void PKMB_hdr_to_inf(PKMB_tapehdr * hdr, infodata * idata);
void print_PKMB_hdr(PKMB_tapehdr * hdr);
void convert_PKMB_point(unsigned char *bits, unsigned char *bytes);
