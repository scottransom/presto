      SUBROUTINE sla_CS2C6 ( A, B, R, AD, BD, RD, V )
*+
*     - - - - - -
*      C S 2 C 6
*     - - - - - -
*
*  Conversion of position & velocity in spherical coordinates
*  to Cartesian coordinates (single precision)
*
*  Given:
*     A     r      longitude (radians)
*     B     r      latitude (radians)
*     R     r      radial coordinate
*     AD    r      longitude derivative (radians per unit time)
*     BD    r      latitude derivative (radians per unit time)
*     RD    r      radial derivative
*
*  Returned:
*     V     r(6)   Cartesian position & velocity vector
*
*  Last revision:   11 September 2005
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

      REAL A, B, R, AD, BD, RD, V(6)

      REAL SA, CA, SB, CB, RCB, X, Y, RBD, W



*  Useful functions.
      SA = SIN(A)
      CA = COS(A)
      SB = SIN(B)
      CB = COS(B)
      RCB = R*CB
      X = RCB*CA
      Y = RCB*SA
      RBD = R*BD
      W = RBD*SB-CB*RD

*  Position.
      V(1) = X
      V(2) = Y
      V(3) = R*SB

*  Velocity.
      V(4) = -Y*AD-W*CA
      V(5) = X*AD-W*SA
      V(6) = RBD*CB+SB*RD

      END
