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
  int channel;    /* Channel of detection (0--numchan-1) */
  int intnum;     /* Integration detected (0--numint-1)  */
} rfi_instance;

typedef struct RFI_OBS {
  float freq_avg;    /* Average frequency of detection (Hz)  */
  float freq_var;    /* Variance of the measured frequencies */
  int number;        /* Number of times the RFI was detected */
  rfi_instance *rfi; /* The individual detections            */
} rfi_obs;
