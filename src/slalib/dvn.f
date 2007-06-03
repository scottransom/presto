      SUBROUTINE sla_DVN (V, UV, VM)
*+
*     - - - -
*      D V N
*     - - - -
*
*  Normalizes a 3-vector also giving the modulus (double precision)
*
*  Given:
*     V       d(3)      vector
*
*  Returned:
*     UV      d(3)      unit vector in direction of V
*     VM      d         modulus of V
*
*  Notes:
*
*  1  If the modulus of V is zero, UV is set to zero as well.
*
*  2  To comply with the ANSI Fortran 77 standard, V and UV must be
*     different arrays.  However, the routine is coded so as to work
*     properly on most platforms even if this rule is violated.
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

      DOUBLE PRECISION V(3),UV(3),VM

      INTEGER I
      DOUBLE PRECISION W1,W2


*  Modulus.
      W1 = 0D0
      DO I=1,3
         W2 = V(I)
         W1 = W1+W2*W2
      END DO
      W1 = SQRT(W1)
      VM = W1

*  Normalize the vector.
      IF (W1.LE.0D0) W1 = 1D0
      DO I=1,3
         UV(I) = V(I)/W1
      END DO

      END
