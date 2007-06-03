      SUBROUTINE sla_DAV2M (AXVEC, RMAT)
*+
*     - - - - - -
*      D A V 2 M
*     - - - - - -
*
*  Form the rotation matrix corresponding to a given axial vector.
*  (double precision)
*
*  A rotation matrix describes a rotation about some arbitrary axis,
*  called the Euler axis.  The "axial vector" supplied to this routine
*  has the same direction as the Euler axis, and its magnitude is the
*  amount of rotation in radians.
*
*  Given:
*    AXVEC  d(3)     axial vector (radians)
*
*  Returned:
*    RMAT   d(3,3)   rotation matrix
*
*  If AXVEC is null, the unit matrix is returned.
*
*  The reference frame rotates clockwise as seen looking along
*  the axial vector from the origin.
*
*  Last revision:   26 November 2005
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

      DOUBLE PRECISION AXVEC(3),RMAT(3,3)

      DOUBLE PRECISION X,Y,Z,PHI,S,C,W



*  Rotation angle - magnitude of axial vector - and functions
      X = AXVEC(1)
      Y = AXVEC(2)
      Z = AXVEC(3)
      PHI = SQRT(X*X+Y*Y+Z*Z)
      S = SIN(PHI)
      C = COS(PHI)
      W = 1D0-C

*  Euler axis - direction of axial vector (perhaps null)
      IF (PHI.NE.0D0) THEN
         X = X/PHI
         Y = Y/PHI
         Z = Z/PHI
      END IF

*  Compute the rotation matrix
      RMAT(1,1) = X*X*W+C
      RMAT(1,2) = X*Y*W+Z*S
      RMAT(1,3) = X*Z*W-Y*S
      RMAT(2,1) = X*Y*W-Z*S
      RMAT(2,2) = Y*Y*W+C
      RMAT(2,3) = Y*Z*W+X*S
      RMAT(3,1) = X*Z*W+Y*S
      RMAT(3,2) = Y*Z*W-X*S
      RMAT(3,3) = Z*Z*W+C

      END
