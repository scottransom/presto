typedef struct MASK { 
  double sigma;            /* Cutoff sigma                           */
  double mjd;              /* MJD of time zero                       */
  double dtint;            /* Time in sec of each interval           */
  double lofreq;           /* Freq (MHz) of low channel              */
  double dfreq;            /* Channel width (MHz)                    */
  int numchan;             /* Number of channels                     */
  int numint;              /* Number of intervals                    */
  int ptsperint;           /* Points per interval                    */
  int num_user_chans;      /* Number of channels the user zapped     */
  int *user_chans;         /* The channels the user zapped           */
  int num_user_ints;       /* Number of intervals the user zapped    */
  int *user_ints;          /* The intervals the user zapped          */
  int *num_chans_per_int;  /* Number of channels zapped per interval */
  int **chans;             /* The channels zapped                    */
} mask;


void fill_mask(double sigma, double mjd, double dtint, double lofreq, 
	       double dfreq, int numchan, int numint, int ptsperint, 
	       int num_user_chans, int *user_chans, int num_user_ints, 
	       int *user_ints, unsigned char **bytemask, mask *obsmask);
/* Fill a mask structure with the appropriate values */

void free_mask(mask obsmask);
/* Free the contents of an mask structure */

void read_mask(char *maskfilenm, mask *obsmask);
/* Read the contents of a mask structure from a file */

void write_mask(char *maskfilenm, mask *obsmask);
/* Write the contents of an mask structure to a file */

