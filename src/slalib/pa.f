      DOUBLE PRECISION FUNCTION sla_PA (HA, DEC, PHI)
*+
*     - - -
*      P A
*     - - -
*
*  HA, Dec to Parallactic Angle (double precision)
*
*  Given:
*     HA     d     hour angle in radians (geocentric apparent)
*     DEC    d     declination in radians (geocentric apparent)
*     PHI    d     observatory latitude in radians (geodetic)
*
*  The result is in the range -pi to +pi
*
*  Notes:
*
*  1)  The parallactic angle at a point in the sky is the position
*      angle of the vertical, i.e. the angle between the direction to
*      the pole and to the zenith.  In precise applications care must
*      be taken only to use geocentric apparent HA,Dec and to consider
*      separately the effects of atmospheric refraction and telescope
*      mount errors.
*
*  2)  At the pole a zero result is returned.
*
*  P.T.Wallace   Starlink   16 August 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION HA,DEC,PHI

      DOUBLE PRECISION CP,SQSZ,CQSZ



      CP=COS(PHI)
      SQSZ=CP*SIN(HA)
      CQSZ=SIN(PHI)*COS(DEC)-CP*SIN(DEC)*COS(HA)
      IF (SQSZ.EQ.0D0.AND.CQSZ.EQ.0D0) CQSZ=1D0
      sla_PA=ATAN2(SQSZ,CQSZ)

      END
