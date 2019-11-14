/*    Pulsar Data Generation Program       */
/*          by Scott Ransom                */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "chkio.h"
#include "orbint.h"

#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif

#define BUFFSIZE 16384		/* Number of floats to write at a time */

typedef struct MAKEDATA {
  char basefilenm[200];	        /* Data file name without suffix         */
  char description[200];	/* Data description                      */
  long N;		        /* Number of bins in the time series     */
  long next2_to_n;              /* The next power-of-2 >= N              */
  double dt;	 	        /* Width of each time series bin (sec)   */
  double T;                     /* Total length of time series (sec)     */
  char ptype[20];		/* Pulsetype (Sine, Crab, Spike, Gauss)  */
  int pnum;                     /* 1=Sine, 2=Crab-like, 3=Spike, 4=Gauss */
  double fwhm;                  /* FWHM Phase (0-1) if ptype is Spike    */
  char round[20];		/* Rounding Format (Whole, Fractional)   */
  int roundnum;                 /* 1=Whole Numbers, 0=Fractional         */
  double f;                     /* Pulsar frequency (hz)                 */
  double fd;                    /* Pulsar frequency deriv (hz/s)         */
  double fdd;                   /* Pulsar 2nd frequency deriv (hz/s^2)   */
  double p;                     /* Pulsar period (s)                     */
  double pd;                    /* Pulsar period deriv (s/s)             */
  double pdd;                   /* Pulsar 2nd period deriv (s/s^2)       */
  double r;                     /* Pulsar FFT freq (bin)                 */
  double z;                     /* Pulsar FFT freq deriv (bins/T)        */
  double w;                     /* Pulsar FFT 2nd freq deriv (bins/T^2)  */
  double amp;                   /* Pulse amplitude                       */
  double phs;                   /* Pulse phase (deg)                     */
  double dc;                    /* Average DC background level           */
  int binary;			/* Binary pulsar?  (1=Yes, 0=No)         */
  orbitparams orb;              /* Binary pulsar data                    */
  int ampmod;                   /* Amplitude modulation? (1=Yes, 0=No)   */
  double ampmoda;               /* Amplitude modulation amplitude        */
  double ampmodf;               /* Amplitude modulation frequency (hz)   */
  double ampmodp;               /* Amplitude modulation phase (deg)      */
  char noisetype[20];		/* Noise type (Standard, Other)          */
  int noise;                    /* 1=Standard, 0=Other                   */
  double noisesig;              /* Noise standard deviation              */
  int numonoff;			/* The number of onoff pairs in the data */
  double *onoff;		/* Bin number pairs where obs is "on"    */
} makedata;

/*  Function Prototypes */

void read_mak_input(makedata *mdata);
/* Get information for makedata from stdin.      */

void read_mak_file(char basefilenm[], makedata *mdata);
/* Read the data for makedata from the makefile. */

void write_mak_file(makedata *mdata);
/* Write the data for makedata to the makefile.  */

/* The various functions to describe the fake data pulse shape. */
double modsin(double val);
double crab(double val);
double spike(double val);
double gauss(double val);
