#include "presto.h"
#include "accelsearch_cmd.h"

/* #define ACCEL_USELEN 64000 */
#define ACCEL_USELEN 32000
#define ACCEL_NUMBETWEEN 2
/* Stepsize in Fourier Freq */
#define ACCEL_DR  0.5
/* Reciprocal of ACCEL_DR */
#define ACCEL_RDR 2
/* Stepsize in Fourier F-dot */
#define ACCEL_DZ  2
/* Reciprocal of ACCEL_DZ */
#define ACCEL_RDZ 0.5

typedef struct accelobs{
  long long N;         /* Number of data points in observation */
  long long numbins;   /* Number of spectral bins in the file */
  long long lobin;     /* Lowest spectral bin present in the file */
  long long highestbin;/* Highest spectral bin present in the file */
  int fftlen;          /* Length of short FFTs to us in search */
  int numharm;         /* Number of harmonics searched */
  int numz;            /* Number of f-dots searched */
  int numbetween;      /* Highest fourier freq resolution (2=interbin) */
  int numzap;          /* Number of birdies to zap */
  double dt;           /* Data sample length (s) */           
  double T;            /* Total observation length */
  double rlo;          /* Minimum fourier freq to search */
  double rhi;          /* Maximum fourier freq to search */
  double dr;           /* Stepsize in fourier freq (1/numbetween) */
  double zlo;          /* Minimum fourier fdot to search */
  double zhi;          /* Maximum fourier fdot to search */
  double dz;           /* Stepsize in fourier freq (1/numbetween) */
  double baryv;        /* Average barycentric velocity during observation */
  float nph;           /* Freq 0 level if present (1 otherwise) */
  float sigma;         /* Cutoff sigma to choose a candidate */
  float *powcut;       /* Cutoff powers to choose a cand (per harmsummed) */
  double *lobins;      /* The low Fourier freq boundaries to zap (RFI) */
  double *hibins;      /* The high Fourier freq boundaries to zap (RFI) */
  long long *numindep; /* Number of independent spectra (per harmsummed) */
} accelobs;

typedef struct accelcand{
  float power;         /* Power level (normalized) */
  float sigma;         /* Equivalent sigma based on numindep (above) */
  int numharm;         /* Number of harmonics summed */
  double r;            /* Fourier freq */
  double z;            /* Fourier f-dot */
} accelcand;

typedef struct kernel{
  int z;               /* The fourier f-dot of the kernel */
  int fftlen;          /* Number of complex points in the kernel */
  int numgoodbins;     /* The number of good points you can get back */
  int numbetween;      /* Fourier freq resolution (2=interbin) */
  int kern_half_width; /* Half width (bins) of the raw kernel. */
  fcomplex *data;      /* The FFTd kernel itself */
} kernel;

typedef struct subharminfo{
  int numharm;       /* The number of sub-harmonics */
  int harmnum;       /* The sub-harmonic number (fundamental = numharm) */
  int zmax;          /* The maximum Fourier f-dot for this harmonic */
  int numkern;       /* Number of kernels in the vector */
  kernel *kern;      /* The kernels themselves */
} subharminfo;


subharminfo *create_subharminfo_vect(int numharm, int zmax);
void free_subharminfo_vect(int numharm, subharminfo *shi);
void create_accelobs(FILE *infile, accelobs *obs, 
		     infodata *idata, Cmdline *cmd);
void free_accelobs(accelobs *obs);

