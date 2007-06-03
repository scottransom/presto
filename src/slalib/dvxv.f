      SUBROUTINE sla_DVXV (VA, VB, VC)
*+
*     - - - - -
*      D V X V
*     - - - - -
*
*  Vector product of two 3-vectors  (double precision)
*
*  Given:
*      VA      dp(3)     first vector
*      VB      dp(3)     second vector
*
*  Returned:
*      VC      dp(3)     vector result
*
*  P.T.Wallace   Starlink   March 1986
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

      DOUBLE PRECISION VA(3),VB(3),VC(3)

      DOUBLE PRECISION VW(3)
      INTEGER I


*  Form the vector product VA cross VB
      VW(1)=VA(2)*VB(3)-VA(3)*VB(2)
      VW(2)=VA(3)*VB(1)-VA(1)*VB(3)
      VW(3)=VA(1)*VB(2)-VA(2)*VB(1)

*  Return the result
      DO I=1,3
         VC(I)=VW(I)
      END DO

      END
