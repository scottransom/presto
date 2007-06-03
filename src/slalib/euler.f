      SUBROUTINE sla_EULER (ORDER, PHI, THETA, PSI, RMAT)
*+
*     - - - - - -
*      E U L E R
*     - - - - - -
*
*  Form a rotation matrix from the Euler angles - three successive
*  rotations about specified Cartesian axes (single precision)
*
*  Given:
*    ORDER  c*(*)    specifies about which axes the rotations occur
*    PHI    r        1st rotation (radians)
*    THETA  r        2nd rotation (   "   )
*    PSI    r        3rd rotation (   "   )
*
*  Returned:
*    RMAT   r(3,3)   rotation matrix
*
*  A rotation is positive when the reference frame rotates
*  anticlockwise as seen looking towards the origin from the
*  positive region of the specified axis.
*
*  The characters of ORDER define which axes the three successive
*  rotations are about.  A typical value is 'ZXZ', indicating that
*  RMAT is to become the direction cosine matrix corresponding to
*  rotations of the reference frame through PHI radians about the
*  old Z-axis, followed by THETA radians about the resulting X-axis,
*  then PSI radians about the resulting Z-axis.
*
*  The axis names can be any of the following, in any order or
*  combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
*  axis labelling/numbering conventions apply;  the xyz (=123)
*  triad is right-handed.  Thus, the 'ZXZ' example given above
*  could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
*  is terminated by length or by the first unrecognized character.
*
*  Fewer than three rotations are acceptable, in which case the later
*  angle arguments are ignored.  If all rotations are zero, the
*  identity matrix is produced.
*
*  Called:  sla_DEULER
*
*  P.T.Wallace   Starlink   23 May 1997
*
*  Copyright (C) 1997 Rutherford Appleton Laboratory
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

      CHARACTER*(*) ORDER
      REAL PHI,THETA,PSI,RMAT(3,3)

      INTEGER J,I
      DOUBLE PRECISION W(3,3)



*  Compute matrix in double precision
      CALL sla_DEULER(ORDER,DBLE(PHI),DBLE(THETA),DBLE(PSI),W)

*  Copy the result
      DO J=1,3
         DO I=1,3
            RMAT(I,J) = REAL(W(I,J))
         END DO
      END DO

      END
