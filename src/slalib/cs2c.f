      SUBROUTINE sla_CS2C (A, B, V)
*+
*     - - - - -
*      C S 2 C
*     - - - - -
*
*  Spherical coordinates to direction cosines (single precision)
*
*  Given:
*     A,B      real      spherical coordinates in radians
*                           (RA,Dec), (long,lat) etc.
*
*  Returned:
*     V        real(3)   x,y,z unit vector
*
*  The spherical coordinates are longitude (+ve anticlockwise looking
*  from the +ve latitude pole) and latitude.  The Cartesian coordinates
*  are right handed, with the x axis at zero longitude and latitude, and
*  the z axis at the +ve latitude pole.
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

      REAL A,B,V(3)

      REAL COSB



      COSB = COS(B)

      V(1) = COS(A)*COSB
      V(2) = SIN(A)*COSB
      V(3) = SIN(B)

      END
