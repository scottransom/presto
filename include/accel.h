#include "presto.h"

typedef struct accelobs{
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
} accelobs;

typedef struct accelcand{
  float power;         /* Power level (normalized) */
  float sigma;         /* Equivalent sigma based on numindep (above) */
  int numharm;         /* Number of harmonics summed */
  double r;            /* Fourier freq */
  double z;            /* Fourier f-dot */
} accelcand;

typedef struct accelkern{
  double z;            /* The fourier f-dot of the kernel */
  int fftlen;          /* Number of complex points in the kernel */
  int harmnum;         /* The harmonic number of this kernel */
  int numgoodbins;     /* The number of good points you get back */
  int numbetween;      /* Fourier freq resolution (2=interbin) */
  int kern_half_width; /* Half width (bins) of the raw kernel. */
  fcomplex *kern;      /* The FFTd kernel itself */
} accelkern;

typedef struct accelkernvec{
  int numkerns;        /* Number of kernels in the vector */
  accelkern *kern;     /* The kernels themselves */
} accelkernvec;


