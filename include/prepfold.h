#include "presto.h"
#include "prepfold_cmd.h"
#include "multibeam.h"

/* This causes the barycentric motion to be calculated once per second */

#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

/* structure used to pass information to plotting routine */

typedef struct PREPFOLDINFO {
  double *rawfolds;   /* Raw folds (nsub * npart * proflen points) */
  double *dms;        /* DMs used in the trials */
  double *periods;    /* Periods used in the trials */
  double *pdots;      /* P-dots used in the trials */
  foldstats *stats;   /* Statistics for the raw folds */
  int numdms;         /* Number of 'dms' */
  int numperiods;     /* Number of 'periods' */
  int numpdots;       /* Number of 'pdots' */
  int nsub;           /* Number of frequency subbands folded */
  int npart;          /* Number of folds in time over integration */
  int proflen;        /* Number of bins per profile */
  int numchan;        /* Number of channels for radio data */
  int pstep;          /* Minimum period stepsize in profile phase bins */
  int pdstep;         /* Minimum p-dot stepsize in profile phase bins */
  int dmstep;         /* Minimum DM stepsize in profile phase bins */
  int ndmfact;        /* 2*ndmfact*proflen+1 DMs to search */
  int npfact;         /* 2*npfact*proflen+1 periods and p-dots to search */
  char *filenm;       /* Filename of the folded data */
  char *candnm;       /* String describing the candidate */
  char *telescope;    /* Telescope where observation took place */
  char *pgdev;        /* PGPLOT device to use */
  double dt;          /* Sampling interval of the data */
  double startT;      /* Fraction of observation file to start folding */
  double endT;        /* Fraction of observation file to stop folding */
  double tepoch;      /* Topocentric eopch of data in MJD */
  double bepoch;      /* Barycentric eopch of data in MJD */
  double avgvoverc;   /* Average topocentric velocity */
  double lofreq;      /* Center of low frequency radio channel */
  double chan_wid;    /* Width of each radio channel in MHz */
  double bestdm;      /* Best DM */
  position topo;      /* Best topocentric p, pd, and pdd */
  position bary;      /* Best barycentric p, pd, and pdd */
  position fold;      /* f, fd, and fdd used to fold the initial data */
  orbitparams orb;    /* Barycentric orbital parameters used in folds */
} prepfoldinfo;

/* Some function definitions */

int (*readrec_ptr)(FILE * infile, float *data, int numpts,
		   double *dispdelays, int numsubbands, int numchan);

int read_resid_rec(FILE * file, double *toa, double *obsf);

int read_floats(FILE *file, float *data, int numpts,
		double *dispdelays, int numsubbands, int numchan);
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */

void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);

int dgels_(char *trans, int *mm, int *nn, int *nrhs, 
	   double *aa, int *lda, double *bb, int *ldb, 
	   double *work, int *lwork, int *info);

double switch_pfdot(double pf, double pfdot);

double fdot2phasedelay(double fdot, double time);

double phasedelay2fdot(double phasedelay, double time);

void double2float(double *in, float *out, int numpts);
/* Copy a double vector into a float vector */

void prepfold_plot(prepfoldinfo *in, int xwin);
/* Make the beautiful 1 page prepfold output */

int bary2topo(double *topotimes, double *barytimes, int numtimes, 
	      double fb, double fbd, double fbdd, 
	      double *ft, double *ftd, double *ftdd);
/* Convert a set of barycentric pulsar spin parameters (fb, fbd, fbdd) */
/* into topocentric spin parameters (ft, ftd, ftdd) by performing      */
/* a linear least-squares fit (using LAPACK routine DGELS).  The       */
/* routine equates the pulse phase using topcentric parameters and     */
/* times to the pulse phase using barycentric parameters and times.    */

void init_prepfoldinfo(prepfoldinfo *in);
/* Set all values to 0 or NULL */

void delete_prepfoldinfo(prepfoldinfo *in);
/* Free all dynamic arrays in the prepfold array */

void print_prepfoldinfo(prepfoldinfo *in);
/* Print a prepfoldinfo data structure to STDOUT */

void write_prepfoldinfo(prepfoldinfo *in, char *filename);
/* Write a prepfoldinfo data structure to a binary file */

void read_prepfoldinfo(prepfoldinfo *in, char *filename);
/* Read a prepfoldinfo data structure from a binary file */

int cpgnice_output_2(char *out, double val, double err, int len);
/* Return a string that has the same formatting as       */
/* nice_output_2(), but for PGPLOT.  This way, exponents */
/* are actually in superscript!  Woo-hoo!                */


		   
