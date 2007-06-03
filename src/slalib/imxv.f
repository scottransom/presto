      SUBROUTINE sla_IMXV (RM, VA, VB)
*+
*     - - - - -
*      I M X V
*     - - - - -
*
*  Performs the 3-D backward unitary transformation:
*
*     vector VB = (inverse of matrix RM) * vector VA
*
*  (single precision)
*
*  (n.b.  the matrix must be unitary, as this routine assumes that
*   the inverse and transpose are identical)
*
*  Given:
*     RM       real(3,3)    matrix
*     VA       real(3)      vector
*
*  Returned:
*     VB       real(3)      result vector
*
*  P.T.Wallace   Starlink   November 1984
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

      REAL RM(3,3),VA(3),VB(3)

      INTEGER I,J
      REAL W,VW(3)



*  Inverse of matrix RM * vector VA -> vector VW
      DO J=1,3
         W=0.0
         DO I=1,3
            W=W+RM(I,J)*VA(I)
         END DO
         VW(J)=W
      END DO

*  Vector VW -> vector VB
      DO J=1,3
         VB(J)=VW(J)
      END DO

      END
