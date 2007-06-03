      SUBROUTINE sla_VN (V, UV, VM)
*+
*     - - -
*      V N
*     - - -
*
*  Normalizes a 3-vector also giving the modulus (single precision)
*
*  Given:
*     V       real(3)      vector
*
*  Returned:
*     UV      real(3)      unit vector in direction of V
*     VM      real         modulus of V
*
*  If the modulus of V is zero, UV is set to zero as well
*
*  P.T.Wallace   Starlink   23 November 1995
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

      REAL V(3),UV(3),VM

      INTEGER I
      REAL W1,W2


*  Modulus
      W1=0.0
      DO I=1,3
         W2=V(I)
         W1=W1+W2*W2
      END DO
      W1=SQRT(W1)
      VM=W1

*  Normalize the vector
      IF (W1.LE.0.0) W1=1.0
      DO I=1,3
         UV(I)=V(I)/W1
      END DO

      END
