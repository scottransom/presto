#include "presto.h"

typedef struct accel_obs{
  long long N;         /* Number of data points in observation */
  long long numbins;   /* Number of spectral bins in the file */
  long long lobin;     /* Lowest spectral bin present in the file */
  int fftlen;          /* Length of short FFTs to us in search */
  int numharm;         /* Number of harmonics searched */
  int numz;            /* Number of f-dots searched */
  int numbetween;      /* Highest fourier freq resolution (2=interbin) */
  int numzap;
  double dt;           /* Data sample length (s) */           
  double T;            /* Total observation length */
  double rlo;          /* Minimum fourier freq to search */
  double rhi;          /* Maximum fourier freq to search */
  double dr;           /* Stepsize in fourier freq (1/numbetween) */
  double zlo;          /* Minimum fourier fdot to search */
  double zhi;          /* Maximum fourier fdot to search */
  double dz;           /* Stepsize in fourier freq (1/numbetween) */
  double baryv;        /* Average barycentric velocity during observation */
  float sigma;         /* Cutoff sigma to choose a candidate */
  float *powcut;       /* Cutoff powers to choose a cand (per harmsummed) */
  double *zapfreqs;    /* The central frequencies to zap (RFI) */
  double *zapwidths;   /* The frequency widths to zap (RFI) */
  long long *numindep; /* Number of independent spectra (per harmsummed) */
} accel_obs;

typedef struct accel_cand{
  float power;         /* Power level (normalized) */
  float sigma;         /* Equivalent sigma based on numindep (above) */
  double r;            /* Fourier freq */
  float z;             /* Fourier f-dot */
  float numharm;       /* Number of harmonics summed */
} accel_cand;

typedef struct accel_candvec{
  int numcands;        /* Number of candidates in the vector */
  accel_cand *cand;    /* The candidates themselves */
} accel_candvec;

typedef struct accel_kern{
  int fftlen;          /* Number of points in the work ffts */
  int z;               /* The fourier f-dot of the kernel */
  int harmnum;         /* The harmonic number of this kernel */
  int numgoodbins;     /* The number of good points you get back */
  int numbetween;      /* Fourier freq resolution (2=interbin) */
  int kern_half_width; /* Half width (bins) of the raw kernel. */
  fcomplex *kern;      /* The FFTd kernel itself */
} accel_kern;

typedef struct accel_kernvec{
  int numkerns;        /* Number of kernels in the vector */
  accel_kern *kern;    /* The kernels themselves */
} accel_kernvec;


