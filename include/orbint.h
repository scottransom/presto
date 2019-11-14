#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef SOL
#define SOL 299792458.0
#endif
#ifndef DBLCORRECT
#define DBLCORRECT    1e-14
#endif
#ifndef MAXIT
#define MAXIT 1000
#endif

#ifndef ORBITPARAMS_TYPE
typedef struct orbitparams {
    double p;	    /* Orbital period (s)                            */
    double e;	    /* Orbital eccentricity                          */
    double x;	    /* Projected semi-major axis (lt-sec)            */
    double w;	    /* Longitude of periapsis (deg)                  */
    double t;	    /* Time since last periastron passage (s)        */
    double pd;	    /* Orbital period derivative (s/yr)              */
    double wd;	    /* Advance of longitude of periapsis (deg/yr)    */
} orbitparams;
#define ORBITPARAMS_TYPE 1
#endif

/* Function declarations */

/*
 *   For all the following:
 *   E           = Eccentric anomaly (radians)
 *   Eo          = initial condition of E
 *   Eacc        = accuracy (1e-14 should be plenty good)
 *   t           = Time (sec)
 *   to          = orbital integration start time
 *   tend        = orbital integration ending time
 *   dt          = sampling interval for the integrator
 *   f_orb       = orbital frequency in hertz
 *   p_orb       = orbital period in seconds
 *   p_psr       = pulsar period in units of choice
 *   e or e_orb  = orbital eccentricity
 *   x or x_orb  = projected semi-major axis of orbit in lt-sec
 *   w or w_orb  = longitude of periastron (degrees)
 */

double *dorbint(double Eo, long numpts, double dt, orbitparams *orb); 
/* This routine integrates Kepler's Equation and returns a double       */
/* vector of the eccentric anomalys (E) for each point.  The initial    */
/* value for eccentric anomaly (usually determined by using             */
/* keplers_equation()) goes in Eo.  The time increment to use is dt,    */
/* total number of pts goes in 'numpts' and all of the various orbital  */
/* parameters are found in *orb.  The routine uses 4th order Runge-     */
/* Kutta in a dumb mode (no adaptive step-size) since all we want is    */
/* tabulated results with even intervals.                               */

double keplers_eqn(double t, double p_orb, double e, double Eacc);
/* This routine solves Kepler's Equation at a single time t (sec) and  */
/* returns the value of the eccentric anomaly.  The orbital period (s) */
/* is in p_orb and the orbital eccentricity is in e.  Eacc is the      */
/* absolute accuracy in E that we want to achieve.  t is the time in   */
/* seconds since the last periapsis.  Uses Newton-Raphson.             */


double lin_interp_E(double *E, double currenttime, double to, \
		    double dt, double maxt);
/* Return a linearly interpolated value of E at time 'time'.  */
/* to is the starting time and dt is the time interval of the */
/* evenly tabulated dat vector *E.                            */


/* Convert eccentric anomalies returned by dorbint into: */

void E_to_phib(double *E, long numpoints, orbitparams *orb);
/*   Phase delays: */

void E_to_v(double *E, long numpoints, orbitparams *orb);
/*   Pulsar line-of-sight velocity (km/s): */

void E_to_p(double *E, long numpoints, double p_psr, orbitparams *orb);
/*   Pulse period: */

void E_to_z(double *E, long numpoints, double p_psr, double T, \
	    orbitparams *orb);
/*   Fourier f-dot */

void E_to_phib_BT(double *E, long numpoints, orbitparams *orb);
/* Convert eccentric anomalys (*E) to time delays */
/* using Blanford and Teukolsky Equations         */
/* This model is NOT currently in use.            */


