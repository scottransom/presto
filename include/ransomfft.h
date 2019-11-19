#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include "meminfo.h"
#include "vectors.h"
#include "multifiles.h"

#ifndef SWAP
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

#ifndef _FCOMPLEX_DECLARED_
typedef struct FCOMPLEX {
    float r, i;
} fcomplex;
#define _FCOMPLEX_DECLARED_
#endif				/* _FCOMPLEX_DECLARED_ */

#ifndef _DCOMPLEX_DECLARED_
typedef struct DCOMPLEX {
    double r, i;
} dcomplex;
#define _DCOMPLEX_DECLARED_
#endif				/* _DCOMPLEX_DECLARED_ */

#ifndef _RAWTYPE_DECLARED_
typedef fcomplex rawtype;
#define _RAWTYPE_DECLARED_
#endif				/* _RAWTYPE_DECLARED_ */

/* Figure out which FFTs we are using */

#include "fftw3.h"
void fftwcall(fcomplex *indata, long nn, int isign);
void read_wisdom(void);
#ifndef COMPLEXFFT
#define COMPLEXFFT    fftwcall
#endif

#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef PIBYTWO
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#endif

/* Function Prototypes   */

/* The FAST FFTs we will use                                     */
/* Note:  A standard forward transform (with negative exponent)  */
/*        is called by setting isign=-1                          */

long long good_factor(long long nn);
void tablesixstepfft(fcomplex *indata, long nn, int isign);
void realfft(float data[], long n, int isign);
void twopassfft_scratch(multifile* infile, multifile* scratch,
			long long nn, int isign);
void realfft_scratch_fwd(multifile* infile, multifile* scratch, 
			 long long nn);
void realfft_scratch_inv(multifile* infile, multifile* scratch, 
			 long long nn);
void realsingfft(FILE * bigfft[5], long numdata, int isign, \
		 char *inpath, char *outpath);
void fourew(FILE * file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE * file[5], unsigned long nn[], int ndim, int isign);
short transpose_float(float *a, int nx, int ny, unsigned char *move, 
		      int move_size);
short transpose_fcomplex(fcomplex *a, int nx, int ny, unsigned char *move, 
			 int move_size);
rawtype *gen_rawvect(long length);
