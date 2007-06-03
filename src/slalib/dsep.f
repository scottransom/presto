      DOUBLE PRECISION FUNCTION sla_DSEP (A1, B1, A2, B2)
*+
*     - - - - -
*      D S E P
*     - - - - -
*
*  Angle between two points on a sphere.
*
*  (double precision)
*
*  Given:
*     A1,B1    d     spherical coordinates of one point
*     A2,B2    d     spherical coordinates of the other point
*
*  (The spherical coordinates are [RA,Dec], [Long,Lat] etc, in radians.)
*
*  The result is the angle, in radians, between the two points.  It
*  is always positive.
*
*  Called:  sla_DCS2C, sla_DSEPV
*
*  Last revision:   7 May 2000
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

      DOUBLE PRECISION A1,B1,A2,B2

      DOUBLE PRECISION V1(3),V2(3)
      DOUBLE PRECISION sla_DSEPV



*  Convert coordinates from spherical to Cartesian.
      CALL sla_DCS2C(A1,B1,V1)
      CALL sla_DCS2C(A2,B2,V2)

*  Angle between the vectors.
      sla_DSEP = sla_DSEPV(V1,V2)

      END
