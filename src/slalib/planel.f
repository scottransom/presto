      SUBROUTINE sla_PLANEL (DATE, JFORM, EPOCH, ORBINC, ANODE, PERIH,
     :                       AORQ, E, AORL, DM, PV, JSTAT)
*+
*     - - - - - - -
*      P L A N E L
*     - - - - - - -
*
*  Heliocentric position and velocity of a planet, asteroid or comet,
*  starting from orbital elements.
*
*  Given:
*     DATE     d     date, Modified Julian Date (JD - 2400000.5, Note 1)
*     JFORM    i     choice of element set (1-3; Note 3)
*     EPOCH    d     epoch of elements (TT MJD, Note 4)
*     ORBINC   d     inclination (radians)
*     ANODE    d     longitude of the ascending node (radians)
*     PERIH    d     longitude or argument of perihelion (radians)
*     AORQ     d     mean distance or perihelion distance (AU)
*     E        d     eccentricity
*     AORL     d     mean anomaly or longitude (radians, JFORM=1,2 only)
*     DM       d     daily motion (radians, JFORM=1 only)
*
*  Returned:
*     PV       d(6)  heliocentric x,y,z,xdot,ydot,zdot of date,
*                                     J2000 equatorial triad (AU,AU/s)
*     JSTAT    i     status:  0 = OK
*                            -1 = illegal JFORM
*                            -2 = illegal E
*                            -3 = illegal AORQ
*                            -4 = illegal DM
*                            -5 = numerical error
*
*  Called:  sla_EL2UE, sla_UE2PV
*
*  Notes
*
*  1  DATE is the instant for which the prediction is required.  It is
*     in the TT timescale (formerly Ephemeris Time, ET) and is a
*     Modified Julian Date (JD-2400000.5).
*
*  2  The elements are with respect to the J2000 ecliptic and equinox.
*
*  3  A choice of three different element-set options is available:
*
*     Option JFORM = 1, suitable for the major planets:
*
*       EPOCH  = epoch of elements (TT MJD)
*       ORBINC = inclination i (radians)
*       ANODE  = longitude of the ascending node, big omega (radians)
*       PERIH  = longitude of perihelion, curly pi (radians)
*       AORQ   = mean distance, a (AU)
*       E      = eccentricity, e (range 0 to <1)
*       AORL   = mean longitude L (radians)
*       DM     = daily motion (radians)
*
*     Option JFORM = 2, suitable for minor planets:
*
*       EPOCH  = epoch of elements (TT MJD)
*       ORBINC = inclination i (radians)
*       ANODE  = longitude of the ascending node, big omega (radians)
*       PERIH  = argument of perihelion, little omega (radians)
*       AORQ   = mean distance, a (AU)
*       E      = eccentricity, e (range 0 to <1)
*       AORL   = mean anomaly M (radians)
*
*     Option JFORM = 3, suitable for comets:
*
*       EPOCH  = epoch of elements and perihelion (TT MJD)
*       ORBINC = inclination i (radians)
*       ANODE  = longitude of the ascending node, big omega (radians)
*       PERIH  = argument of perihelion, little omega (radians)
*       AORQ   = perihelion distance, q (AU)
*       E      = eccentricity, e (range 0 to 10)
*
*     Unused arguments (DM for JFORM=2, AORL and DM for JFORM=3) are not
*     accessed.
*
*  4  Each of the three element sets defines an unperturbed heliocentric
*     orbit.  For a given epoch of observation, the position of the body
*     in its orbit can be predicted from these elements, which are
*     called "osculating elements", using standard two-body analytical
*     solutions.  However, due to planetary perturbations, a given set
*     of osculating elements remains usable for only as long as the
*     unperturbed orbit that it describes is an adequate approximation
*     to reality.  Attached to such a set of elements is a date called
*     the "osculating epoch", at which the elements are, momentarily,
*     a perfect representation of the instantaneous position and
*     velocity of the body.
*
*     Therefore, for any given problem there are up to three different
*     epochs in play, and it is vital to distinguish clearly between
*     them:
*
*     . The epoch of observation:  the moment in time for which the
*       position of the body is to be predicted.
*
*     . The epoch defining the position of the body:  the moment in time
*       at which, in the absence of purturbations, the specified
*       position (mean longitude, mean anomaly, or perihelion) is
*       reached.
*
*     . The osculating epoch:  the moment in time at which the given
*       elements are correct.
*
*     For the major-planet and minor-planet cases it is usual to make
*     the epoch that defines the position of the body the same as the
*     epoch of osculation.  Thus, only two different epochs are
*     involved:  the epoch of the elements and the epoch of observation.
*
*     For comets, the epoch of perihelion fixes the position in the
*     orbit and in general a different epoch of osculation will be
*     chosen.  Thus, all three types of epoch are involved.
*
*     For the present routine:
*
*     . The epoch of observation is the argument DATE.
*
*     . The epoch defining the position of the body is the argument
*       EPOCH.
*
*     . The osculating epoch is not used and is assumed to be close
*       enough to the epoch of observation to deliver adequate accuracy.
*       If not, a preliminary call to sla_PERTEL may be used to update
*       the element-set (and its associated osculating epoch) by
*       applying planetary perturbations.
*
*  5  The reference frame for the result is with respect to the mean
*     equator and equinox of epoch J2000.
*
*  6  The algorithm was originally adapted from the EPHSLA program of
*     D.H.P.Jones (private communication, 1996).  The method is based
*     on Stumpff's Universal Variables.
*
*  Reference:  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
*
*  P.T.Wallace   Starlink   31 December 2002
*
*  Copyright (C) 2002 Rutherford Appleton Laboratory
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
      DOUBLE PRECISION EPOCH,ORBINC,ANODE,PERIH,AORQ,E,AORL,DM,PV(6)
      INTEGER JSTAT

      DOUBLE PRECISION U(13)
      INTEGER J



*  Validate elements and convert to "universal variables" parameters.
      CALL sla_EL2UE(DATE,JFORM,
     :               EPOCH,ORBINC,ANODE,PERIH,AORQ,E,AORL,DM,U,J)

*  Determine the position and velocity.
      IF (J.EQ.0) THEN
         CALL sla_UE2PV(DATE,U,PV,J)
         IF (J.NE.0) J=-5
      END IF

*  Wrap up.
      JSTAT = J

      END
