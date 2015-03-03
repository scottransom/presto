/* ACCEL_USELEN must be less than 65536 since we  */
/* use unsigned short ints to index our arrays... */
/* #define ACCEL_USELEN 32160 */
/* #define ACCEL_USELEN 16000 */
#define ACCEL_USELEN 7560 
#define ACCEL_NUMBETWEEN 2
/* Stepsize in Fourier Freq */
#define ACCEL_DR  0.5
/* Reciprocal of ACCEL_DR */
#define ACCEL_RDR 2
/* Stepsize in Fourier F-dot */
#define ACCEL_DZ  2
/* Reciprocal of ACCEL_DZ */
#define ACCEL_RDZ 0.5
/* Closest candidates we will accept as independent */
#define ACCEL_CLOSEST_R 15.0
/* Padding for .dat file reading so that we don't SEGFAULT */
#define ACCEL_PADDING 2000


#ifndef DBLCORRECT
#define DBLCORRECT    1e-14
#endif

#ifndef _ACCEL_UTILS_GPU_
#define _ACCEL_UTILS_GPU_

typedef struct accel_cand_gpu{
	float pow;					/*pow of selected candidate*/
	int		nof_cand;			/*number of candidates in sub_array/plane */
	int		z_ind;				/*z_index of the selected candidate*/
	int		r_ind;				/*r_index of the selected candidate*/
}	accel_cand_gpu;

#ifndef _FCOMPLEX_DECLARED_
typedef struct fcomplex {
    float r, i;
} fcomplex;
#define _FCOMPLEX_DECLARED_
#endif				/* _FCOMPLEX_DECLARED_ */

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
  unsigned short *rinds; /* Table of indices for Fourier Freqs */
} subharminfo;


/*  Constants used in the correlation/convolution routines */
typedef enum {
  CONV, CORR, INPLACE_CONV, INPLACE_CORR
} presto_optype;

typedef enum {
  FFTDK, FFTD, FFTK, NOFFTS
} presto_ffts;

#endif
