#define RFI_NUMHARMSUM 4
#define RFI_NUMBETWEEN 2
#define RFI_LOBIN      5
#define RFI_FRACTERROR 0.001
#define NUM_RFI_OBS    30

#define GET_BIT(c, n) (*(c+(n>>3)) >> (7-(n&7)) & 1)
#define SET_BIT(c, n) (*(c+(n>>3)) |= 1 << (7-(n&7)))
#define UNSET_BIT(c, n) (*(c+(n>>3)) &= ~(1 << (7-(n&7))))

typedef enum {
  GOOD, USER, POW, AVG, VAR
} mask_flags;

typedef struct RFI_INSTANCE {
  float freq;     /* Frequency (Hz) of the RFI           */
  float power;    /* Normalized power of measurement     */
  float fftbin;   /* FFT bin number of detection         */
  float fftbins;  /* Number of FFT bins                  */
  float inttime;  /* Time duration of integration (s)    */
  float sigma;    /* True significance                   */
  int channel;    /* Channel of detection (0--numchan-1) */
  int intnum;     /* Integration detected (0--numint-1)  */
  int numsum;     /* Number of bins that were summed     */
} rfi_instance;

typedef struct RFI_OBS {
  float freq_avg;    /* Average frequency of detection (Hz)  */
  float freq_var;    /* Variance of the measured frequencies */
  int number;        /* Number of times the RFI was detected */
  rfi_instance *rfi; /* The individual detections            */
} rfi_obs;

rfi_obs *create_rfi_obs(rfi_instance rfi);
/* Create and initialize a rfi_obs structure */

void free_rfi_obs(rfi_obs *rfi);
/* Free an rfi_obs structure and its contents */

void free_rfi_obs_vector(rfi_obs **rfi_vect, int num_rfi_vect);
/* Free a vector that holds rfi_obs structures */

void add_rfi_instance(rfi_obs *old, rfi_instance new);
/* Add an instance of RFI to a rfi_obs structure */

int find_rfi(rfi_obs **rfi_vect, int num_rfi_vect, 
	     double freq, double fract_error);
/* Try to find a birdie in an rfi_obs vector.  Compare     */
/* all currently known birdies with the new freq.  If it   */
/* finds one with a freq within fractional error, it       */
/* returns the number of the birdie -- otherwise, -1.      */

int compare_rfi_obs(void *ca, void *cb);
/*  Used as compare function for qsort() */

void percolate_rfi_obs(rfi_obs **list, int nlist);
/*  Pushes an rfi_obs structure as far up a sorted list of */
/*  structs as it needs to go to keep the list sorted.     */
/*  The new rfi_obs structure is in position nlist-1.      */
