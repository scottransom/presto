#define GOODDATA 0x00
#define PADDING  0x01
#define OLDMASK  0x02
#define USERCHAN 0x04
#define USERINTS 0x08
#define BAD_POW  0x10
#define BAD_STD  0x20
#define BAD_AVG  0x40
#define BADDATA  (BAD_POW|BAD_STD|BAD_AVG)
#define USERZAP  (USERCHAN|USERINTS)

#ifndef MASK_DEFINED
typedef struct MASK { 
  double timesigma;        /* Cutoff time-domain sigma               */
  double freqsigma;        /* Cutoff freq-domain sigma               */
  double mjd;              /* MJD of time zero                       */
  double dtint;            /* Time in sec of each interval           */
  double lofreq;           /* Freq (MHz) of low channel              */
  double dfreq;            /* Channel width (MHz)                    */
  int numchan;             /* Number of channels                     */
  int numint;              /* Number of intervals                    */
  int ptsperint;           /* Points per interval                    */
  int num_zap_chans;       /* Number of full channels to zap         */
  int *zap_chans;          /* The full channels to zap               */
  int num_zap_ints;        /* Number of full intervals to zap        */
  int *zap_ints;           /* The full intervals to zap              */
  int *num_chans_per_int;  /* Number of channels zapped per interval */
  int **chans;             /* The channels zapped                    */
} mask;
#define MASK_DEFINED
#endif

void fill_mask(double timesigma, double freqsigma, double mjd, 
	       double dtint, double lofreq, double dfreq, 
	       int numchan, int numint, int ptsperint, 
	       int num_zap_chans, int *zap_chans, int num_zap_ints, 
	       int *zap_ints, unsigned char **bytemask, mask *obsmask);
/* Fill a mask structure with the appropriate values */

void set_oldmask_bits(mask *oldmask, unsigned char **bytemask);
/* Sets the oldmask bit in the appropriate bytes in bytemask */

void unset_oldmask_bits(mask *oldmask, unsigned char **bytemask);
/* Unsets the oldmask bits in bytemask */

void free_mask(mask obsmask);
/* Free the contents of an mask structure */

void read_mask(char *maskfilenm, mask *obsmask);
/* Read the contents of a mask structure from a file */

void write_mask(char *maskfilenm, mask *obsmask);
/* Write the contents of an mask structure to a file */

int check_mask(double starttime, double duration, mask *obsmask, 
	       int *maskchans);
/* Return value is the number of channels to mask.  The */
/* channel numbers are placed in maskchans (which must  */
/* have a length of numchan).  If -1 is returned, all   */
/* channels should be masked.                           */

void calc_avgmedstd(float *arr, int numarr, float fraction, 
		    int step, float *avg, float *med, float *std);
/* Calculates the median and middle-'fraction' std deviation  */
/* and average of the array 'arr'.  Values are returned in    */
/* 'avg', 'med' and 'std'.  The array is not modified.        */

int determine_padvals(char *maskfilenm, mask * obsmask, float *padvals);
// Determine reasonable padding values from the rfifind produced
// *.stats file if it is available.  The pre-allocated vector (of
// length numchan) is in padvals.  Return a '1' if the routine used
// the stats file, return 0 if the padding was set to aeros.
