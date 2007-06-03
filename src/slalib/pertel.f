      SUBROUTINE sla_PERTEL (JFORM, DATE0, DATE1,
     :                 EPOCH0, ORBI0, ANODE0, PERIH0, AORQ0, E0, AM0,
     :                 EPOCH1, ORBI1, ANODE1, PERIH1, AORQ1, E1, AM1,
     :                       JSTAT)
*+
*     - - - - - - -
*      P E R T E L
*     - - - - - - -
*
*  Update the osculating orbital elements of an asteroid or comet by
*  applying planetary perturbations.
*
*  Given (format and dates):
*     JFORM   i    choice of element set (2 or 3; Note 1)
*     DATE0   d    date of osculation (TT MJD) for the given elements
*     DATE1   d    date of osculation (TT MJD) for the updated elements
*
*  Given (the unperturbed elements):
*     EPOCH0  d    epoch (TT MJD) of the given element set (Note 2)
*     ORBI0   d    inclination (radians)
*     ANODE0  d    longitude of the ascending node (radians)
*     PERIH0  d    argument of perihelion (radians)
*     AORQ0   d    mean distance or perihelion distance (AU)
*     E0      d    eccentricity
*     AM0     d    mean anomaly (radians, JFORM=2 only)
*
*  Returned (the updated elements):
*     EPOCH1  d    epoch (TT MJD) of the updated element set (Note 2)
*     ORBI1   d    inclination (radians)
*     ANODE1  d    longitude of the ascending node (radians)
*     PERIH1  d    argument of perihelion (radians)
*     AORQ1   d    mean distance or perihelion distance (AU)
*     E1      d    eccentricity
*     AM1     d    mean anomaly (radians, JFORM=2 only)
*
*  Returned (status flag):
*     JSTAT   i    status: +102 = warning, distant epoch
*                          +101 = warning, large timespan ( > 100 years)
*                     +1 to +10 = coincident with planet (Note 6)
*                             0 = OK
*                            -1 = illegal JFORM
*                            -2 = illegal E0
*                            -3 = illegal AORQ0
*                            -4 = internal error
*                            -5 = numerical error
*
*  Notes:
*
*  1  Two different element-format options are available:
*
*     Option JFORM=2, suitable for minor planets:
*
*     EPOCH   = epoch of elements (TT MJD)
*     ORBI    = inclination i (radians)
*     ANODE   = longitude of the ascending node, big omega (radians)
*     PERIH   = argument of perihelion, little omega (radians)
*     AORQ    = mean distance, a (AU)
*     E       = eccentricity, e
*     AM      = mean anomaly M (radians)
*
*     Option JFORM=3, suitable for comets:
*
*     EPOCH   = epoch of perihelion (TT MJD)
*     ORBI    = inclination i (radians)
*     ANODE   = longitude of the ascending node, big omega (radians)
*     PERIH   = argument of perihelion, little omega (radians)
*     AORQ    = perihelion distance, q (AU)
*     E       = eccentricity, e
*
*  2  DATE0, DATE1, EPOCH0 and EPOCH1 are all instants of time in
*     the TT timescale (formerly Ephemeris Time, ET), expressed
*     as Modified Julian Dates (JD-2400000.5).
*
*     DATE0 is the instant at which the given (i.e. unperturbed)
*     osculating elements are correct.
*
*     DATE1 is the specified instant at which the updated osculating
*     elements are correct.
*
*     EPOCH0 and EPOCH1 will be the same as DATE0 and DATE1
*     (respectively) for the JFORM=2 case, normally used for minor
*     planets.  For the JFORM=3 case, the two epochs will refer to
*     perihelion passage and so will not, in general, be the same as
*     DATE0 and/or DATE1 though they may be similar to one another.
*
*  3  The elements are with respect to the J2000 ecliptic and equinox.
*
*  4  Unused elements (AM0 and AM1 for JFORM=3) are not accessed.
*
*  5  See the sla_PERTUE routine for details of the algorithm used.
*
*  6  This routine is not intended to be used for major planets, which
*     is why JFORM=1 is not available and why there is no opportunity
*     to specify either the longitude of perihelion or the daily
*     motion.  However, if JFORM=2 elements are somehow obtained for a
*     major planet and supplied to the routine, sensible results will,
*     in fact, be produced.  This happens because the sla_PERTUE routine
*     that is called to perform the calculations checks the separation
*     between the body and each of the planets and interprets a
*     suspiciously small value (0.001 AU) as an attempt to apply it to
*     the planet concerned.  If this condition is detected, the
*     contribution from that planet is ignored, and the status is set to
*     the planet number (1-10 = Mercury, Venus, EMB, Mars, Jupiter,
*     Saturn, Uranus, Neptune, Earth, Moon) as a warning.
*
*  Reference:
*
*     Sterne, Theodore E., "An Introduction to Celestial Mechanics",
*     Interscience Publishers Inc., 1960.  Section 6.7, p199.
*
*  Called:  sla_EL2UE, sla_PERTUE, sla_UE2EL
*
*  This revision:   19 June 2004
*
*  Copyright (C) 2004 P.T.Wallace.  All rights reserved.
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
      INTEGER JFORM
      DOUBLE PRECISION DATE0,DATE1,
     :                 EPOCH0,ORBI0,ANODE0,PERIH0,AORQ0,E0,AM0,
     :                 EPOCH1,ORBI1,ANODE1,PERIH1,AORQ1,E1,AM1
      INTEGER JSTAT

      DOUBLE PRECISION U(13),DM
      INTEGER J,JF



*  Check that the elements are either minor-planet or comet format.
      IF (JFORM.LT.2.OR.JFORM.GT.3) THEN
         JSTAT = -1
         GO TO 9999
      ELSE

*     Provisionally set the status to OK.
         JSTAT = 0
      END IF

*  Transform the elements from conventional to universal form.
      CALL sla_EL2UE(DATE0,JFORM,EPOCH0,ORBI0,ANODE0,PERIH0,
     :               AORQ0,E0,AM0,0D0,U,J)
      IF (J.NE.0) THEN
         JSTAT = J
         GO TO 9999
      END IF

*  Update the universal elements.
      CALL sla_PERTUE(DATE1,U,J)
      IF (J.GT.0) THEN
         JSTAT = J
      ELSE IF (J.LT.0) THEN
         JSTAT = -5
         GO TO 9999
      END IF

*  Transform from universal to conventional elements.
      CALL sla_UE2EL(U,JFORM,
     :               JF, EPOCH1, ORBI1, ANODE1, PERIH1,
     :               AORQ1, E1, AM1, DM, J)
      IF (JF.NE.JFORM.OR.J.NE.0) JSTAT=-5

 9999 CONTINUE
      END
