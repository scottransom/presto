/* Number of entries in PSR database */
#define NP  3000
#define NBP 200

/* This is a structure that contains the "normal" information in the database */
typedef struct PSRDATA {
  char jname[13];            /* The PSRs J2000 name          */
  char bname[9];             /* The PSRs B1950 name          */
  char alias[10];            /* An alias for the pulsar      */
  double ra2000;             /* J2000 RA (radians)           */
  double ra2000err;          /* J2000 RA error (radians)     */
  double dec2000;            /* J2000 DEC (radians)          */
  double dec2000err;         /* J2000 DEC error (radians)    */
  double p;                  /* PSR period (s)               */
  double perr;               /* PSR period error (s)         */
  double pd;                 /* PSR period deriv (s/s)       */
  double pderr;              /* PSR period deriv error (s/s) */
  double dm;                 /* Dispersion Measure           */
  double dmerr;              /* Dispersion Measure error     */
  double timepoch;           /* MJD epoch for timing         */
  double binary;             /* Flag for binary pulsars      */
} psrdata;

typedef struct BINPSRDATA {
  double pb;                 /* orbital period (d)           */
  double pberr;              /* orbital period error (d)     */
  double x;                  /* asin(i)/c (s)                */
  double xerr;               /* asin(i)/c error (s)          */
  double e;                  /* eccentricity                 */
  double eerr;               /* eccentricity error           */
  double w;                  /* angle of peri (deg)          */
  double werr;               /* angle of peri error (deg)    */
  double To;                 /* Time of peri (MJD)           */
  double Toerr;              /* Time of peri error (d)       */
} binpsrdata;


/* This is a structure that contains the 'key' pulsar info  */
typedef struct PSRPARAMS {
  char jname[13];            /* The PSRs J2000 name         */
  char bname[9];             /* The PSRs B1950 name         */
  char alias[10];            /* An alias for the pulsar     */
  double ra2000;             /* J2000 RA                    */
  double dec2000;            /* J2000 DEC                   */
  double dm;                 /* Dispersion Measure          */
  double timepoch;           /* MJD epoch for timing        */
  double p;                  /* PSR period (s)              */
  double pd;                 /* PSR period deriv (s/s)      */
  double pdd;                /* Period 2nd deriv (s/s^2)    */
  double f;                  /* PSR frequency (hz)          */
  double fd;                 /* PSR frequency deriv (s^-2)  */
  double fdd;                /* Frequency 2nd deriv (s^-3)  */
  orbitparams orb;           /* Orbital parameters          */
} psrparams;

int read_database(void);
/* Reads the full pulsar database into the static array psrdata */

void get_psrparams(psrparams *psr, char *psrname);
/* Read a full pulsar database entry for pulsar psrname. */
/* Return the data in a psrparams structure.             */

void get_psr(int psrnumber, psrparams *psr);
/* Returns a full database entry for the pulsar #psrnumber in */
/* the database psrdata.  Returns *psr completed.             */

int psr_number_from_name(char *psrname);
/* Returns the pulsar number of psrname from the database */
/* This number can be from zero to the total number       */
/* of pulsars minus 1.  This way you can use this number  */
/* as an index from the result of collect_psrparams().    */
/* Return -1 if no pulsar is found.                       */

int get_psr_at_epoch(char *psrname, double epoch, psrparams *psr);
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.                                      */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */

