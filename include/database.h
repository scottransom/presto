/* Number of "slots" in PSR database */
#define NP 800

/* The full pulsar database is stored in the following structure */
typedef struct PSRDATABASE {
  char jname[9600], bname[6400];
  double ra2000[NP], ra1950[NP], rae[NP], dec2000[NP], \
         dec1950[NP], dece[NP];
  int nscode[NP];
  double dmin__[NP], dmax__[NP], dist[NP];
  int ndflag[NP];
  char lcode[NP], ucode[NP];
  double ldeg[NP], bdeg[NP], pmra[NP], pmrae[NP], pmdec[NP], \
         pmdece[NP], posepoch[NP], p[NP], pe[NP], pdot[NP], \
         pdote[NP], f2[NP], f2e[NP], f3[NP], f3e[NP], epoch[NP], \
         dm[NP], dme[NP], rm[NP], rme[NP], we[NP], w50[NP], \
         w10[NP], s400[NP], s600[NP], s1400[NP], tau[NP];
  int ntauflag[NP];
  double t408[NP];
  int ntype[NP], modcode[NP], limcode[NP];
  double distmod[NP], lum[NP], bsurf[NP], age[NP], edot[NP];
  int ibin[NP];
  double pb[NP], pbe[NP], a1[NP], a1e[NP], om[NP], ome[NP], \
         omdot[NP], omdote[NP], e[NP], ee[NP], t0[NP], t0e[NP], \
         gamma[NP], gammae[NP], pbdot[NP], pbdote[NP], si[NP], \
         sie[NP], r__[NP], re[NP], pb2[NP], pb2e[NP], a12[NP], \
         a12e[NP], om2[NP], om2e[NP], omdot2[NP], omdot2e[NP], \
         e2[NP], e2e[NP], t02[NP], t02e[NP], gamma2[NP], \
         gamma2e[NP], pbdot2[NP], pbdot2e[NP], si2[NP], \
         si2e[NP], r2[NP], r2e[NP];
} psrdatabase;


/* A complete pulsar database entry for a single pulsar */
typedef struct PSRDATA {
  double ra2000, ra1950, rae, dec2000, dec1950, dece, \
    dmin__, dmax__, dist, ldeg, bdeg, pmra, pmrae, pmdec, \
    pmdece, posepoch, p, pe, pdot, pdote, f2, f2e, f3, f3e, \
    epoch, dm, dme, rm, rme, we, w50, w10, s400, s600, s1400, \
    tau, t408, distmod, lum, bsurf, age, edot, pb, pbe, a1, \
    a1e, om, ome, omdot, omdote, e, ee, t0, t0e, gamma, \
    gammae, pbdot, pbdote, si, sie, r__, re, pb2, pb2e, a12, \
    a12e, om2, om2e, omdot2, omdot2e, e2, e2e, t02, t02e, \
    gamma2, gamma2e, pbdot2, pbdot2e, si2, si2e, r2, r2e;
  int nscode, ndflag, ntauflag, ntype, modcode, limcode, ibin;
  char jname[12], bname[8], lcode, ucode;
} psrdata;


/* This is a structure that contains the 'key' pulsar info  */
typedef struct PSRPARAMS {
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

int read_database(psrdatabase * pdata);
/* Reads the full pulsar database into the structure *pdata */

void collect_psrdata(psrdata ** psrs, int *np);
/* Read the full pulsar database and place it in an array    */
/* of psrdata structures.  Return the number of psrs stored. */

void get_psrdata_by_num(psrdata * psr, int pnum);
/* Read a full pulsar database entry for the pulsar pnum */
/* in the database.  Return the data in psrdata.         */

void get_psrdata(psrdata * psr, char * psrname);
/* Read a full pulsar database entry for pulsar psrname. */
/* Return the data in a psrdata structure.               */

void get_psr(int psrnumber, psrdata * psrinfo, psrdatabase * pdata);
/* Returns a full database entry for the pulsar #psrnumber in */
/* the database variable pdata.  Returns *psrinfo completed.  */

int psr_number_from_name(char * psrname, psrdatabase * pdata);
/* Returns the pulsar number of psrname from the database */

int return_psrparams_at_epoch(psrparams * psr, char * psrname, \
			      double epoch);
/* Reads info from the pulsar database and converts returned values */
/* to epoch 'epoch'.  Returned values go in psr (psrparams).        */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */

int get_psr_at_epoch(char * psrname, double epoch, psrdatabase * pdata, \
		     psrparams * psr);
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.  The database data is in *pdata.     */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */

