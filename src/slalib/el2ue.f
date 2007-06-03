      SUBROUTINE sla_EL2UE (DATE, JFORM, EPOCH, ORBINC, ANODE,
     :                      PERIH, AORQ, E, AORL, DM,
     :                      U, JSTAT)
*+
*     - - - - - -
*      E L 2 U E
*     - - - - - -
*
*  Transform conventional osculating orbital elements into "universal"
*  form.
*
*  Given:
*     DATE    d      epoch (TT MJD) of osculation (Note 3)
*     JFORM   i      choice of element set (1-3, Note 6)
*     EPOCH   d      epoch (TT MJD) of the elements
*     ORBINC  d      inclination (radians)
*     ANODE   d      longitude of the ascending node (radians)
*     PERIH   d      longitude or argument of perihelion (radians)
*     AORQ    d      mean distance or perihelion distance (AU)
*     E       d      eccentricity
*     AORL    d      mean anomaly or longitude (radians, JFORM=1,2 only)
*     DM      d      daily motion (radians, JFORM=1 only)
*
*  Returned:
*     U       d(13)  universal orbital elements (Note 1)
*
*               (1)  combined mass (M+m)
*               (2)  total energy of the orbit (alpha)
*               (3)  reference (osculating) epoch (t0)
*             (4-6)  position at reference epoch (r0)
*             (7-9)  velocity at reference epoch (v0)
*              (10)  heliocentric distance at reference epoch
*              (11)  r0.v0
*              (12)  date (t)
*              (13)  universal eccentric anomaly (psi) of date, approx
*
*     JSTAT   i      status:  0 = OK
*                            -1 = illegal JFORM
*                            -2 = illegal E
*                            -3 = illegal AORQ
*                            -4 = illegal DM
*                            -5 = numerical error
*
*  Called:  sla_UE2PV, sla_PV2UE
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
*  2  The companion routine is sla_UE2PV.  This takes the set of numbers
*     that the present routine outputs and uses them to derive the
*     object's position and velocity.  A single prediction requires one
*     call to the present routine followed by one call to sla_UE2PV;
*     for convenience, the two calls are packaged as the routine
*     sla_PLANEL.  Multiple predictions may be made by again calling the
*     present routine once, but then calling sla_UE2PV multiple times,
*     which is faster than multiple calls to sla_PLANEL.
*
*  3  DATE is the epoch of osculation.  It is in the TT timescale
*     (formerly Ephemeris Time, ET) and is a Modified Julian Date
*     (JD-2400000.5).
*
*  4  The supplied orbital elements are with respect to the J2000
*     ecliptic and equinox.  The position and velocity parameters
*     returned in the array U are with respect to the mean equator and
*     equinox of epoch J2000, and are for the perihelion prior to the
*     specified epoch.
*
*  5  The universal elements returned in the array U are in canonical
*     units (solar masses, AU and canonical days).
*
*  6  Three different element-format options are available:
*
*     Option JFORM=1, suitable for the major planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = longitude of perihelion, curly pi (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e (range 0 to <1)
*     AORL   = mean longitude L (radians)
*     DM     = daily motion (radians)
*
*     Option JFORM=2, suitable for minor planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e (range 0 to <1)
*     AORL   = mean anomaly M (radians)
*
*     Option JFORM=3, suitable for comets:
*
*     EPOCH  = epoch of perihelion (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = perihelion distance, q (AU)
*     E      = eccentricity, e (range 0 to 10)
*
*  7  Unused elements (DM for JFORM=2, AORL and DM for JFORM=3) are
*     not accessed.
*
*  8  The algorithm was originally adapted from the EPHSLA program of
*     D.H.P.Jones (private communication, 1996).  The method is based
*     on Stumpff's Universal Variables.
*
*  Reference:  Everhart & Pitkin, Am.J.Phys. 51, 712 (1983).
*
*  Last revision:   8 September 2005
*
*  Copyright P.T.Wallace.  All rights reserved.
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

      DOUBLE PRECISION DATE
      INTEGER JFORM
      DOUBLE PRECISION EPOCH,ORBINC,ANODE,PERIH,AORQ,E,AORL,DM,U(13)
      INTEGER JSTAT

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Sin and cos of J2000 mean obliquity (IAU 1976)
      DOUBLE PRECISION SE,CE
      PARAMETER (SE=0.3977771559319137D0,
     :           CE=0.9174820620691818D0)

      INTEGER J

      DOUBLE PRECISION PHT,ARGPH,Q,W,CM,ALPHA,PHS,SW,CW,SI,CI,SO,CO,
     :                 X,Y,Z,PX,PY,PZ,VX,VY,VZ,DT,FC,FP,PSI,
     :                 UL(13),PV(6)



*  Validate arguments.
      IF (JFORM.LT.1.OR.JFORM.GT.3) THEN
         JSTAT = -1
         GO TO 9999
      END IF
      IF (E.LT.0D0.OR.E.GT.10D0.OR.(E.GE.1D0.AND.JFORM.NE.3)) THEN
         JSTAT = -2
         GO TO 9999
      END IF
      IF (AORQ.LE.0D0) THEN
         JSTAT = -3
         GO TO 9999
      END IF
      IF (JFORM.EQ.1.AND.DM.LE.0D0) THEN
         JSTAT = -4
         GO TO 9999
      END IF

*
*  Transform elements into standard form:
*
*  PHT   = epoch of perihelion passage
*  ARGPH = argument of perihelion (little omega)
*  Q     = perihelion distance (q)
*  CM    = combined mass, M+m (mu)

      IF (JFORM.EQ.1) THEN

*     Major planet.
         PHT = EPOCH-(AORL-PERIH)/DM
         ARGPH = PERIH-ANODE
         Q = AORQ*(1D0-E)
         W = DM/GCON
         CM = W*W*AORQ*AORQ*AORQ

      ELSE IF (JFORM.EQ.2) THEN

*     Minor planet.
         PHT = EPOCH-AORL*SQRT(AORQ*AORQ*AORQ)/GCON
         ARGPH = PERIH
         Q = AORQ*(1D0-E)
         CM = 1D0

      ELSE

*     Comet.
         PHT = EPOCH
         ARGPH = PERIH
         Q = AORQ
         CM = 1D0

      END IF

*  The universal variable alpha.  This is proportional to the total
*  energy of the orbit:  -ve for an ellipse, zero for a parabola,
*  +ve for a hyperbola.

      ALPHA = CM*(E-1D0)/Q

*  Speed at perihelion.

      PHS = SQRT(ALPHA+2D0*CM/Q)

*  In a Cartesian coordinate system which has the x-axis pointing
*  to perihelion and the z-axis normal to the orbit (such that the
*  object orbits counter-clockwise as seen from +ve z), the
*  perihelion position and velocity vectors are:
*
*    position   [Q,0,0]
*    velocity   [0,PHS,0]
*
*  To express the results in J2000 equatorial coordinates we make a
*  series of four rotations of the Cartesian axes:
*
*           axis      Euler angle
*
*     1      z        argument of perihelion (little omega)
*     2      x        inclination (i)
*     3      z        longitude of the ascending node (big omega)
*     4      x        J2000 obliquity (epsilon)
*
*  In each case the rotation is clockwise as seen from the +ve end of
*  the axis concerned.

*  Functions of the Euler angles.
      SW = SIN(ARGPH)
      CW = COS(ARGPH)
      SI = SIN(ORBINC)
      CI = COS(ORBINC)
      SO = SIN(ANODE)
      CO = COS(ANODE)

*  Position at perihelion (AU).
      X = Q*CW
      Y = Q*SW
      Z = Y*SI
      Y = Y*CI
      PX = X*CO-Y*SO
      Y = X*SO+Y*CO
      PY = Y*CE-Z*SE
      PZ = Y*SE+Z*CE

*  Velocity at perihelion (AU per canonical day).
      X = -PHS*SW
      Y = PHS*CW
      Z = Y*SI
      Y = Y*CI
      VX = X*CO-Y*SO
      Y = X*SO+Y*CO
      VY = Y*CE-Z*SE
      VZ = Y*SE+Z*CE

*  Time from perihelion to date (in Canonical Days: a canonical day
*  is 58.1324409... days, defined as 1/GCON).

      DT = (DATE-PHT)*GCON

*  First approximation to the Universal Eccentric Anomaly, PSI,
*  based on the circle (FC) and parabola (FP) values.

      FC = DT/Q
      W = (3D0*DT+SQRT(9D0*DT*DT+8D0*Q*Q*Q))**(1D0/3D0)
      FP = W-2D0*Q/W
      PSI = (1D0-E)*FC+E*FP

*  Assemble local copy of element set.
      UL(1) = CM
      UL(2) = ALPHA
      UL(3) = PHT
      UL(4) = PX
      UL(5) = PY
      UL(6) = PZ
      UL(7) = VX
      UL(8) = VY
      UL(9) = VZ
      UL(10) = Q
      UL(11) = 0D0
      UL(12) = DATE
      UL(13) = PSI

*  Predict position+velocity at epoch of osculation.
      CALL sla_UE2PV(DATE,UL,PV,J)
      IF (J.NE.0) GO TO 9010

*  Convert back to universal elements.
      CALL sla_PV2UE(PV,DATE,CM-1D0,U,J)
      IF (J.NE.0) GO TO 9010

*  OK exit.
      JSTAT = 0
      GO TO 9999

*  Quasi-impossible numerical errors.
 9010 CONTINUE
      JSTAT = -5

 9999 CONTINUE
      END
