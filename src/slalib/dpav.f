      DOUBLE PRECISION FUNCTION sla_DPAV ( V1, V2 )
*+
*     - - - - -
*      D P A V
*     - - - - -
*
*  Position angle of one celestial direction with respect to another.
*
*  (double precision)
*
*  Given:
*     V1    d(3)    direction cosines of one point
*     V2    d(3)    direction cosines of the other point
*
*  (The coordinate frames correspond to RA,Dec, Long,Lat etc.)
*
*  The result is the bearing (position angle), in radians, of point
*  V2 with respect to point V1.  It is in the range +/- pi.  The
*  sense is such that if V2 is a small distance east of V1, the
*  bearing is about +pi/2.  Zero is returned if the two points
*  are coincident.
*
*  V1 and V2 need not be unit vectors.
*
*  The routine sla_DBEAR performs an equivalent function except
*  that the points are specified in the form of spherical
*  coordinates.
*
*  Last revision:   16 March 2005
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

      DOUBLE PRECISION V1(3),V2(3)

      DOUBLE PRECISION X1,Y1,Z1,W,X2,Y2,Z2,SQ,CQ



*  The unit vector to point 1.
      X1 = V1(1)
      Y1 = V1(2)
      Z1 = V1(3)
      W = SQRT(X1*X1+Y1*Y1+Z1*Z1)
      IF (W.NE.0D0) THEN
         X1 = X1/W
         Y1 = Y1/W
         Z1 = Z1/W
      END IF

*  The vector to point 2.
      X2 = V2(1)
      Y2 = V2(2)
      Z2 = V2(3)

*  Position angle.
      SQ = Y2*X1-X2*Y1
      CQ = Z2*(X1*X1+Y1*Y1)-Z1*(X2*X1+Y2*Y1)
      IF (SQ.EQ.0D0.AND.CQ.EQ.0D0) CQ=1D0
      sla_DPAV = ATAN2(SQ,CQ)

      END
