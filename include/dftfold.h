#include "rawtype.h"

typedef struct DFTVECTOR {
  int n;            /* Number of data points in each vector segment   */
  int numvect;      /* Number of vectors in the fold (vector length)  */
  double dt;        /* Duration (s) of each data point                */
  double r;         /* Fourier frequency folded                       */
  double norm;      /* Normalization used on vector                   */
  double T;         /* Total duration of original data set.           */
  fcomplex *vector; /* Calculated DFT vector                          */
} dftvector;

void read_dftvector(dftvector *data, char *filename);
/* Read a dftvector data structure from a binary file */

void write_dftvector(dftvector *data, char *filename);
/* Write a dftvector data structure to a binary file */

void init_dftvector(dftvector *data, int n, int numvect, 
		    double dt, double r, double norm,
		    double T);
/* Initialize a dftvector and allocate its vector part */

void free_dftvector(dftvector *data);
/* Free the dynamically allocated vector in data */
