      SUBROUTINE sla_UE2PV ( DATE, U, PV, JSTAT )
*+
*     - - - - - -
*      U E 2 P V
*     - - - - - -
*
*  Heliocentric position and velocity of a planet, asteroid or comet,
*  starting from orbital elements in the "universal variables" form.
*
*  Given:
*     DATE     d       date, Modified Julian Date (JD-2400000.5)
*
*  Given and returned:
*     U        d(13)   universal orbital elements (updated; Note 1)
*
*       given    (1)   combined mass (M+m)
*         "      (2)   total energy of the orbit (alpha)
*         "      (3)   reference (osculating) epoch (t0)
*         "    (4-6)   position at reference epoch (r0)
*         "    (7-9)   velocity at reference epoch (v0)
*         "     (10)   heliocentric distance at reference epoch
*         "     (11)   r0.v0
*     returned  (12)   date (t)
*         "     (13)   universal eccentric anomaly (psi) of date
*
*  Returned:
*     PV       d(6)    position (AU) and velocity (AU/s)
*     JSTAT    i       status:  0 = OK
*                              -1 = radius vector zero
*                              -2 = failed to converge
*
*  Notes
*
*  1  The "universal" elements are those which define the orbit for the
*     purposes of the method of universal variables (see reference).
*     They consist of the combined mass of the two bodies, an epoch,
*     and the position and velocity vectors (arbitrary reference frame)
*     at that epoch.  The parameter set used here includes also various
*     quantities that can, in fact, be derived from the other
*     information.  This approach is taken to avoiding unnecessary
*     computation and loss of accuracy.  The supplementary quantities
*     are (i) alpha, which is proportional to the total energy of the
*     orbit, (ii) the heliocentric distance at epoch, (iii) the
*     outwards component of the velocity at the given epoch, (iv) an
*     estimate of psi, the "universal eccentric anomaly" at a given
*     date and (v) that date.
*
*  2  The companion routine is sla_EL2UE.  This takes the conventional
*     orbital elements and transforms them into the set of numbers
*     needed by the present routine.  A single prediction requires one
*     one call to sla_EL2UE followed by one call to the present routine;
*     for convenience, the two calls are packaged as the routine
*     sla_PLANEL.  Multiple predictions may be made by again
*     calling sla_EL2UE once, but then calling the present routine
*     multiple times, which is faster than multiple calls to sla_PLANEL.
*
*     It is not obligatory to use sla_EL2UE to obtain the parameters.
*     However, it should be noted that because sla_EL2UE performs its
*     own validation, no checks on the contents of the array U are made
*     by the present routine.
*
*  3  DATE is the instant for which the prediction is required.  It is
*     in the TT timescale (formerly Ephemeris Time, ET) and is a
*     Modified Julian Date (JD-2400000.5).
*
*  4  The universal elements supplied in the array U are in canonical
*     units (solar masses, AU and canonical days).  The position and
*     velocity are not sensitive to the choice of reference frame.  The
*     sla_EL2UE routine in fact produces coordinates with respect to the
*     J2000 equator and equinox.
*
*  5  The algorithm was originally adapted from the EPHSLA program of
*     D.H.P.Jones (private communication, 1996).  The method is based
*     on Stumpff's Universal Variables.
*
*  Reference:  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
*
*  P.T.Wallace   Starlink   22 October 2005
*
*  Copyright (C) 2005 Rutherford Appleton Laboratory
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      DOUBLE PRECISION DATE,U(13),PV(6)
      INTEGER JSTAT

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Canonical days to seconds
      DOUBLE PRECISION CD2S
      PARAMETER (CD2S=GCON/86400D0)

*  Test value for solution and maximum number of iterations
      DOUBLE PRECISION TEST
      INTEGER NITMAX
      PARAMETER (TEST=1D-13,NITMAX=25)

      INTEGER I,NIT,N

      DOUBLE PRECISION CM,ALPHA,T0,P0(3),V0(3),R0,SIGMA0,T,PSI,DT,W,
     :                 TOL,PSJ,PSJ2,BETA,S0,S1,S2,S3,
     :                 FF,R,FLAST,PLAST,F,G,FD,GD



*  Unpack the parameters.
      CM = U(1)
      ALPHA = U(2)
      T0 = U(3)
      DO I=1,3
         P0(I) = U(I+3)
         V0(I) = U(I+6)
      END DO
      R0 = U(10)
      SIGMA0 = U(11)
      T = U(12)
      PSI = U(13)

*  Approximately update the universal eccentric anomaly.
      PSI = PSI+(DATE-T)*GCON/R0

*  Time from reference epoch to date (in Canonical Days: a canonical
*  day is 58.1324409... days, defined as 1/GCON).
      DT = (DATE-T0)*GCON

*  Refine the universal eccentric anomaly, psi.
      NIT = 1
      W = 1D0
      TOL = 0D0
      DO WHILE (ABS(W).GE.TOL)

*     Form half angles until BETA small enough.
         N = 0
         PSJ = PSI
         PSJ2 = PSJ*PSJ
         BETA = ALPHA*PSJ2
         DO WHILE (ABS(BETA).GT.0.7D0)
            N = N+1
            BETA = BETA/4D0
            PSJ = PSJ/2D0
            PSJ2 = PSJ2/4D0
         END DO

*     Calculate Universal Variables S0,S1,S2,S3 by nested series.
         S3 = PSJ*PSJ2*((((((BETA/210D0+1D0)
     :                      *BETA/156D0+1D0)
     :                      *BETA/110D0+1D0)
     :                      *BETA/72D0+1D0)
     :                      *BETA/42D0+1D0)
     :                      *BETA/20D0+1D0)/6D0
         S2 = PSJ2*((((((BETA/182D0+1D0)
     :                  *BETA/132D0+1D0)
     :                  *BETA/90D0+1D0)
     :                  *BETA/56D0+1D0)
     :                  *BETA/30D0+1D0)
     :                  *BETA/12D0+1D0)/2D0
         S1 = PSJ+ALPHA*S3
         S0 = 1D0+ALPHA*S2

*     Undo the angle-halving.
         TOL = TEST
         DO WHILE (N.GT.0)
            S3 = 2D0*(S0*S3+PSJ*S2)
            S2 = 2D0*S1*S1
            S1 = 2D0*S0*S1
            S0 = 2D0*S0*S0-1D0
            PSJ = PSJ+PSJ
            TOL = TOL+TOL
            N = N-1
         END DO

*     Values of F and F' corresponding to the current value of psi.
         FF = R0*S1+SIGMA0*S2+CM*S3-DT
         R = R0*S0+SIGMA0*S1+CM*S2

*     If first iteration, create dummy "last F".
         IF ( NIT.EQ.1) FLAST = FF

*     Check for sign change.
         IF ( FF*FLAST.LT.0D0 ) THEN

*        Sign change:  get psi adjustment using secant method.
            W = FF*(PLAST-PSI)/(FLAST-FF)
         ELSE

*        No sign change:  use Newton-Raphson method instead.
            IF (R.EQ.0D0) GO TO 9010
            W = FF/R
         END IF

*     Save the last psi and F values.
         PLAST = PSI
         FLAST = FF

*     Apply the Newton-Raphson or secant adjustment to psi.
         PSI = PSI-W

*     Next iteration, unless too many already.
         IF (NIT.GT.NITMAX) GO TO 9020
         NIT = NIT+1
      END DO

*  Project the position and velocity vectors (scaling velocity to AU/s).
      W = CM*S2
      F = 1D0-W/R0
      G = DT-CM*S3
      FD = -CM*S1/(R0*R)
      GD = 1D0-W/R
      DO I=1,3
         PV(I) = P0(I)*F+V0(I)*G
         PV(I+3) = CD2S*(P0(I)*FD+V0(I)*GD)
      END DO

*  Update the parameters to allow speedy prediction of PSI next time.
      U(12) = DATE
      U(13) = PSI

*  OK exit.
      JSTAT = 0
      GO TO 9999

*  Null radius vector.
 9010 CONTINUE
      JSTAT = -1
      GO TO 9999

*  Failed to converge.
 9020 CONTINUE
      JSTAT = -2

 9999 CONTINUE
      END
