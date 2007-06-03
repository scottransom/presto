      REAL FUNCTION sla_SEP (A1, B1, A2, B2)
*+
*     - - - -
*      S E P
*     - - - -
*
*  Angle between two points on a sphere.
*
*  (single precision)
*
*  Given:
*     A1,B1    r     spherical coordinates of one point
*     A2,B2    r     spherical coordinates of the other point
*
*  (The spherical coordinates are [RA,Dec], [Long,Lat] etc, in radians.)
*
*  The result is the angle, in radians, between the two points.  It
*  is always positive.
*
*  Called:  sla_DSEP
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

      REAL A1,B1,A2,B2

      DOUBLE PRECISION sla_DSEP



*  Use double precision version.
      sla_SEP = REAL(sla_DSEP(DBLE(A1),DBLE(B1),DBLE(A2),DBLE(B2)))

      END
