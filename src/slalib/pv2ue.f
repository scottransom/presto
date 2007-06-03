      SUBROUTINE sla_PV2UE (PV, DATE, PMASS, U, JSTAT)
*+
*     - - - - - -
*      P V 2 U E
*     - - - - - -
*
*  Construct a universal element set based on an instantaneous position
*  and velocity.
*
*  Given:
*     PV        d(6)   heliocentric x,y,z,xdot,ydot,zdot of date,
*                      (AU,AU/s; Note 1)
*     DATE      d      date (TT Modified Julian Date = JD-2400000.5)
*     PMASS     d      mass of the planet (Sun=1; Note 2)
*
*  Returned:
*     U         d(13)  universal orbital elements (Note 3)
*
*                 (1)  combined mass (M+m)
*                 (2)  total energy of the orbit (alpha)
*                 (3)  reference (osculating) epoch (t0)
*               (4-6)  position at reference epoch (r0)
*               (7-9)  velocity at reference epoch (v0)
*                (10)  heliocentric distance at reference epoch
*                (11)  r0.v0
*                (12)  date (t)
*                (13)  universal eccentric anomaly (psi) of date, approx
*
*     JSTAT     i      status:  0 = OK
*                              -1 = illegal PMASS
*                              -2 = too close to Sun
*                              -3 = too slow
*
*  Notes
*
*  1  The PV 6-vector can be with respect to any chosen inertial frame,
*     and the resulting universal-element set will be with respect to
*     the same frame.  A common choice will be mean equator and ecliptic
*     of epoch J2000.
*
*  2  The mass, PMASS, is important only for the larger planets.  For
*     most purposes (e.g. asteroids) use 0D0.  Values less than zero
*     are illegal.
*
*  3  The "universal" elements are those which define the orbit for the
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
*  Reference:  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
*
*  P.T.Wallace   Starlink   18 March 1999
*
*  Copyright (C) 1999 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION PV(6),DATE,PMASS,U(13)
      INTEGER JSTAT

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Canonical days to seconds
      DOUBLE PRECISION CD2S
      PARAMETER (CD2S=GCON/86400D0)

*  Minimum allowed distance (AU) and speed (AU per canonical day)
      DOUBLE PRECISION RMIN,VMIN
      PARAMETER (RMIN=1D-3,VMIN=1D-3)

      DOUBLE PRECISION T0,CM,X,Y,Z,XD,YD,ZD,R,V2,V,ALPHA,RDV


*  Reference epoch.
      T0 = DATE

*  Combined mass (mu=M+m).
      IF (PMASS.LT.0D0) GO TO 9010
      CM = 1D0+PMASS

*  Unpack the state vector, expressing velocity in AU per canonical day.
      X = PV(1)
      Y = PV(2)
      Z = PV(3)
      XD = PV(4)/CD2S
      YD = PV(5)/CD2S
      ZD = PV(6)/CD2S

*  Heliocentric distance, and speed.
      R = SQRT(X*X+Y*Y+Z*Z)
      V2 = XD*XD+YD*YD+ZD*ZD
      V = SQRT(V2)

*  Reject unreasonably small values.
      IF (R.LT.RMIN) GO TO 9020
      IF (V.LT.VMIN) GO TO 9030

*  Total energy of the orbit.
      ALPHA = V2-2D0*CM/R

*  Outward component of velocity.
      RDV = X*XD+Y*YD+Z*ZD

*  Construct the universal-element set.
      U(1) = CM
      U(2) = ALPHA
      U(3) = T0
      U(4) = X
      U(5) = Y
      U(6) = Z
      U(7) = XD
      U(8) = YD
      U(9) = ZD
      U(10) = R
      U(11) = RDV
      U(12) = T0
      U(13) = 0D0

*  Exit.
      JSTAT = 0
      GO TO 9999

*  Negative PMASS.
 9010 CONTINUE
      JSTAT = -1
      GO TO 9999

*  Too close.
 9020 CONTINUE
      JSTAT = -2
      GO TO 9999

*  Too slow.
 9030 CONTINUE
      JSTAT = -3

 9999 CONTINUE
      END
