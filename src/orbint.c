#include "orbint.h"
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif

/* The first time derivative of eccentric anomaly.             */
/* Z = eccentric anomaly.  twopif = TWOPI * orbital freq.      */
/* e = eccentricity.                                           */
#ifndef EDOT
#define EDOT(Z) ((twopif) / (1.0 - (e) * cos(Z)))
#endif


double *dorbint(double Eo, long numpts, double dt, orbitparams * orb)
/* This routine integrates Kepler's Equation and returns a double      */
/* vector of the eccentric anomalys (E) for each point.  The initial   */
/* value for eccentric anomaly (usually determined by using            */
/* keplers_equation()) goes in Eo.  The time increment to use is dt,   */
/* total number of pts goes in 'numpts' and all of the various orbital */
/* parameters are found in *orb.  The routine uses 4th order Runge-    */
/* Kutta in a dumb mode (no adaptive step-size) since all we want is   */
/* tabulated results with even intervals.                              */
{
    long ii;
    double k1, k2, k3, k4, dt2, twopif, *E, e, Etmp;

    E = gen_dvect(numpts);
    E[0] = Eo;
    e = orb->e;
    twopif = TWOPI / orb->p;
    dt2 = 0.5 * dt;
    for (ii = 0; ii < numpts - 1; ii++) {
        Etmp = E[ii];
        k1 = EDOT(Etmp);
        k2 = EDOT(Etmp + dt2 * k1);
        k3 = EDOT(Etmp + dt2 * k2);
        k4 = EDOT(Etmp + dt * k3);
        E[ii + 1] = Etmp + dt * (((k1 + k4) * 0.5 + k2 + k3) / 3.0);
    }
    return E;
}


double keplers_eqn(double t, double p_orb, double e, double Eacc)
{
/* This routine solves Kepler's Equation at a single time t (sec) and  */
/* returns the value of the eccentric anomaly.  The orbital period (s) */
/* is in p_orb and the orbital eccentricity is in e.  Eacc is the      */
/* absolute accuracy in E that we want to achieve.  t is the time in   */
/* seconds since the last periapsis.                                   */

    int ii;
    double E1 = 0.0, E2 = PI, z, twopif;
    double df, dE, dEold, f, fh, fl;
    double temp, Eh, El, rEs;

    twopif = TWOPI / p_orb;
    z = twopif * t;
    fl = E1 - e * sin(E1) - z;
    fh = E2 - e * sin(E2) - z;
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        E1 = PI;
        E2 = TWOPI;
        fl = E1 - e * sin(E1) - z;
        fh = E2 - e * sin(E2) - z;
        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
            printf("Problem with the initial conditions in keplers_eqn().\n");
            printf("  t = %g  p_orb = %g  e = %g  Eacc = %g\n", t, p_orb, e, Eacc);
            exit(1);
        }
    }
    if (fl == 0.0)
        return E1;
    if (fh == 0.0)
        return E2;
    if (fl < 0.0) {
        El = E1;
        Eh = E2;
    } else {
        Eh = E1;
        El = E2;
    }
    rEs = 0.5 * (E1 + E2);
    dEold = fabs(E2 - E1);
    dE = dEold;
    f = rEs - e * sin(rEs) - z;
    df = EDOT(rEs);
    for (ii = 0; ii < MAXIT; ii++) {
        if ((((rEs - Eh) * df - f) * ((rEs - El) * df - f) >= 0.0)
            || (fabs(2.0 * f) > fabs(dEold * df))) {
            dEold = dE;
            dE = 0.5 * (Eh - El);
            rEs = El + dE;
            if (El == rEs)
                return rEs;
        } else {
            dEold = dE;
            dE = f / df;
            temp = rEs;
            rEs -= dE;
            if (temp == rEs)
                return rEs;
        }
        if (fabs(dE) < Eacc)
            return rEs;
        f = rEs - e * sin(rEs) - z;
        df = EDOT(rEs);
        if (f < 0.0)
            El = rEs;
        else
            Eh = rEs;
    }
    printf("Maximum number of iterations exceeded in keplers_eqn.\n");
    return 0.0;
}


void E_to_phib_BT(double *E, long numpoints, orbitparams * orb)
/* Convert eccentric anomalys (*E) to time delays */
/* using Blanford and Teukolsky Equations         */
{
    long i;
    double c1, c2, cE, sE, wrad;

    wrad = orb->w * DEGTORAD;
    c1 = orb->x * sin(wrad);
    c2 = orb->x * cos(wrad) * sqrt(1 - orb->e * orb->e);
    for (i = 0; i < numpoints; i++) {
        cE = cos(E[i]);
        sE = sin(E[i]);
        E[i] = (c1 * (cE - orb->e) + c2 * sE) *
            (1.0 - TWOPI * (c2 * cE - c1 * sE) / orb->p) / (1.0 - orb->e * cE);
    }
}


void E_to_v(double *E, long numpoints, orbitparams * orb)
/* Convert eccentric anomalys (*E) to pulsar velocities (km/s) */
{
    long i;
    double c1, c2, c3, tmp = 0.0, wrad;

    wrad = orb->w * DEGTORAD;
    c1 = TWOPI * orb->x / orb->p;
    c2 = cos(wrad) * sqrt(1 - orb->e * orb->e);
    c3 = sin(wrad);
    for (i = 0; i < numpoints; i++) {
        tmp = cos(E[i]);
        E[i] = SOL / 1000.0 * c1 * (c2 * tmp - c3 * sin(E[i]))
            / (1.0 - orb->e * tmp);
    }
}


void E_to_p(double *E, long numpoints, double p_psr, orbitparams * orb)
/* Convert eccentric anomalys (*E) to new pulsar periods */
{
    long i;
    double c1, c2, c3, tmp = 0.0, wrad;

    wrad = orb->w * DEGTORAD;
    c1 = TWOPI * orb->x / orb->p;
    c2 = cos(wrad) * sqrt(1 - orb->e * orb->e);
    c3 = sin(wrad);
    for (i = 0; i < numpoints; i++) {
        tmp = cos(E[i]);
        E[i] = p_psr * (1 + c1 * (c2 * tmp - c3 * sin(E[i]))
                        / (1.0 - orb->e * tmp));
    }
}


void E_to_phib(double *E, long numpoints, orbitparams * orb)
/* Convert eccentric anomalys (*E) to time delays */
/* This is the pulsar orbit Roemer delay.         */
{
    long i;
    double c1, c2, wrad;

    wrad = orb->w * DEGTORAD;
    c1 = orb->x * sin(wrad);
    c2 = orb->x * cos(wrad) * sqrt(1 - orb->e * orb->e);
    for (i = 0; i < numpoints; i++)
        E[i] = c1 * (cos(E[i]) - orb->e) + c2 * sin(E[i]);
}


void E_to_z(double *E, long numpoints, double p_psr, double T, orbitparams * orb)
/* Convert eccentric anomalys (*E) to Fourier f-dots */
{
    long i;
    double c1, c2, c3, tmp = 0.0, wrad;

    wrad = orb->w * DEGTORAD;
    c1 = -TWOPI * TWOPI * T * T * orb->x / (orb->p * orb->p * p_psr);
    c2 = cos(wrad) * sqrt(1 - orb->e * orb->e);
    c3 = sin(wrad);
    for (i = 0; i < numpoints; i++) {
        tmp = cos(E[i]);
        E[i] = c1 * (c2 * sin(E[i]) + c3 * (tmp - orb->e))
            / pow(orb->e * tmp - 1.0, 3.0);
    }
}


double lin_interp_E(double *E, double currenttime, double to, double dt, double maxt)
/* Return a linearly interpolated value of E at time 'currenttime'.  */
/* to is the starting time and dt is the time interval of the        */
/* evenly tabulated dat vector *E.                                   */
{
    int ipart;
    double fpart, dtemp;

    if (currenttime < to || currenttime > maxt) {
        printf("\nInput currenttime out-of-bounds in lin_interp_E().\n");
        printf("Exiting.\n\n");
        exit(1);
    }
    dtemp = (currenttime - to) / dt;
    ipart = (int) (dtemp);
    fpart = dtemp - ipart;
    dtemp = E[ipart];
    return (fpart * (E[ipart + 1] - dtemp)) + dtemp;
}
