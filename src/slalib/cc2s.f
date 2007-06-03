      SUBROUTINE sla_CC2S (V, A, B)
*+
*     - - - - -
*      C C 2 S
*     - - - - -
*
*  Cartesian to spherical coordinates (single precision)
*
*  Given:
*     V     r(3)   x,y,z vector
*
*  Returned:
*     A,B   r      spherical coordinates in radians
*
*  The spherical coordinates are longitude (+ve anticlockwise looking
*  from the +ve latitude pole) and latitude.  The Cartesian coordinates
*  are right handed, with the x axis at zero longitude and latitude, and
*  the z axis at the +ve latitude pole.
*
*  If V is null, zero A and B are returned.  At either pole, zero A is
*  returned.
*
*  Last revision:   22 July 2004
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

      REAL V(3),A,B

      REAL X,Y,Z,R


      X = V(1)
      Y = V(2)
      Z = V(3)
      R = SQRT(X*X+Y*Y)

      IF (R.EQ.0.0) THEN
         A = 0.0
      ELSE
         A = ATAN2(Y,X)
      END IF

      IF (Z.EQ.0.0) THEN
         B = 0.0
      ELSE
         B = ATAN2(Z,R)
      END IF

      END
