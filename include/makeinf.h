/*    Pulsar Data Set Info File Creator     */
/*          by Scott Ransom                 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "chkio.h"

#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif
#define NUMSCOPES  8
#define NUMBANDS   6
#define NUMFILTERS 11

typedef struct INFODATA {
    char name[200];		/* Data file name without suffix         */
    char object[100];		/* Object being observed                 */
    int ra_h;			/* Right ascension hours (J2000)         */
    int ra_m;			/* Right ascension minutes (J2000)       */
    double ra_s;		/* Right ascension seconds (J2000)       */
    int dec_d;			/* Declination degrees (J2000)           */
    int dec_m;			/* Declination minutes (J2000)           */
    double dec_s;		/* Declination seconds (J2000)           */
    char telescope[40];		/* Telescope used                        */
    char instrument[100];	/* Instrument used                       */
    char observer[100];		/* Observer[s] for the data set          */
    double N;		        /* Number of bins in the time series     */
    double dt;	 	        /* Width of each time series bin (sec)   */
    int numonoff;		/* The number of onoff pairs in the data */
    double *onoff;		/* Bin number pairs where obs is "on"    */
    double fov;			/* Diameter of Beam or FOV in arcsec     */
    int mjd_i;			/* Epoch of observation (MJD) int part   */
    double mjd_f;		/* Epoch of observation (MJD) frac part  */
    int bary;			/* Barycentered?  1=yes, 0=no            */
    char band[40];		/* Type of observation (EM band)         */
    double dm;			/* Radio -- Dispersion Measure (cm-3 pc) */
    double freq;		/* Radio -- Low chan central freq (Mhz)  */
    double freqband;		/* Radio -- Total Bandwidth (Mhz)        */
    int num_chan;		/* Radio -- Number Channels              */
    double chan_wid;		/* Radio -- Channel Bandwidth (Mhz)      */
    char filt[7];		/* IR,Opt,UV -- Photometric Filter       */
    double wavelen;		/* IR,Opt,UV -- central wavelength (nm)  */
    double waveband;		/* IR,Opt,UV -- bandpass (nm)            */
    double energy;		/* x-ray,gamma -- central energy (kev)   */
    double energyband;		/* x-ray,gamma -- energy bandpass (kev)  */
    char analyzer[100];		/* Who analyzed the data                 */
    char notes[500];		/* Any additional notes                  */
} infodata;

void readinf(infodata * data, char *filenm);
void writeinf(infodata * data);
