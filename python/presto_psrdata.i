// %module database

// %include numpy.i

// %inline %{
// #include "presto.h"
// %}

%typemap(memberin) char jname[12] {
  sprintf($target, "%.12s", $source);
}       
%typemap(memberin) char bname[8] {
  sprintf($target, "%.8s", $source);
}       
%typemap(memberout) char jname[12] {
  sprintf($target, "%.12s", $source);
}       
%typemap(memberout) char bname[8] {
  sprintf($target, "%.8s", $source);
}       

/* A complete pulsar database entry for a single pulsar */
typedef struct psrdata{
  /* See the PSR database documentation for details */
  double ra2000, ra1950, rae, dec2000, dec1950, dece,
    dmin__, dmax__, dist, ldeg, bdeg, pmra, pmrae, pmdec,
    pmdece, posepoch, p, pe, pdot, pdote, f2, f2e, f3, f3e,
    epoch, dm, dme, rm, rme, we, w50, w10, s400, s600, s1400,
    tau, t408, distmod, lum, bsurf, age, edot, pb, pbe, a1,
    a1e, om, ome, omdot, omdote, e, ee, t0, t0e, gamma,
    gammae, pbdot, pbdote, si, sie, r__, re, pb2, pb2e, a12,
    a12e, om2, om2e, omdot2, omdot2e, e2, e2e, t02, t02e,
    gamma2, gamma2e, pbdot2, pbdot2e, si2, si2e, r2, r2e;
  int nscode, ndflag, ntauflag, ntype, modcode, limcode, ibin;
  char jname[12];
  char bname[8];
  char lcode, ucode;
} psrdata;


%typemap(memberin) char jname[20] {
  sprintf($target, "%.20s", $source);
}       
%typemap(memberin) char bname[20] {
  sprintf($target, "%.20s", $source);
}       
%typemap(memberout) char jname[20] {
  sprintf($target, "%.20s", $source);
}       
%typemap(memberout) char bname[20] {
  sprintf($target, "%.20s", $source);
}       

// typedef struct orbitparams {
//    double p;       /* Orbital period (s)                            */
//    double e;       /* Orbital eccentricity                          */
//    double x;       /* Projected semi-major axis (lt-sec)            */
//    double w;       /* Longitude of periapsis (deg)                  */
//    double t;       /* Time since last periastron passage (s)        */
//    double pd;      /* Orbital period derivative (s/yr)              */
//    double wd;      /* Advance of longitude of periapsis (deg/yr)    */
// } orbitparams;

/* This is a structure that contains the 'key' pulsar info  */
typedef struct psrparams {
  char jname[20];            /* The PSRs J2000 name         */
  char bname[20];            /* The PSRs B1950 name         */
  int ntype;                 /* Pulsar type (see below)     */
  double ra2000;             /* J2000 RA                    */
  double dec2000;            /* J2000 DEC                   */
  double dm;                 /* Dispersion Measure          */
  double dist;               /* Adopted distance (kpc)      */
  double fwhm;               /* FWHM pulse width in ms      */
  double timepoch;           /* MJD epoch for timing        */
  double p;                  /* PSR period (s)              */
  double pd;                 /* PSR period deriv (s/s)      */
  double pdd;                /* Period 2nd deriv (s/s^2)    */
  double f;                  /* PSR frequency (hz)          */
  double fd;                 /* PSR frequency deriv (s^-2)  */
  double fdd;                /* Frequency 2nd deriv (s^-3)  */
  orbitparams orb;           /* Orbital parameters          */
  /*                                                        */
  /* Key to ntype:  (ntype & 1)   =  Globular cluster       */
  /*                (ntype & 2)   =  In a SNR               */
  /*                (ntype & 4)   =  Glitches               */
  /*                (ntype & 8)   =  Binary (or more)       */
  /*                (ntype & 16)  =  Millisecond PSR        */
  /*                (ntype & 32)  =  Recycled PSR           */
  /*                (ntype & 64)  =  Has an interpulse      */
  /*                (ntype & 128) =  High-energy pulsations */
} psrparams;


int num_psrs_in_database(void);
/* Returns the number of entries in the database */

void get_psrdata(psrdata * psr, char * psrname);
/* Read a full pulsar database entry for pulsar psrname. */
/* Return the data in a psrdata structure.               */


void get_psrdata_by_num(psrdata * psr, int pnum);
/* Read a full pulsar database entry for the pulsar pnum */
/* in the database.  Return the data in psrdata.         */

int return_psrparams_at_epoch(psrparams * psr, char * psrname,
			      double epoch);
/* Reads info from the pulsar database and converts returned values */
/* to epoch 'epoch'.  Returned values go in psr (psrparams).        */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */

int get_psr_at_epoch(char * psrname, double epoch, psrdatabase * pdata,
		     psrparams * psr);
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.  The database data is in *pdata.     */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */






