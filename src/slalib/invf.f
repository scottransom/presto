      SUBROUTINE sla_INVF (FWDS,BKWDS,J)
*+
*     - - - - -
*      I N V F
*     - - - - -
*
*  Invert a linear model of the type produced by the sla_FITXY routine.
*
*  Given:
*     FWDS    d(6)      model coefficients
*
*  Returned:
*     BKWDS   d(6)      inverse model
*     J        i        status:  0 = OK, -1 = no inverse
*
*  The models relate two sets of [X,Y] coordinates as follows.
*  Naming the elements of FWDS:
*
*     FWDS(1) = A
*     FWDS(2) = B
*     FWDS(3) = C
*     FWDS(4) = D
*     FWDS(5) = E
*     FWDS(6) = F
*
*  where two sets of coordinates [X1,Y1] and [X2,Y1] are related
*  thus:
*
*     X2 = A + B*X1 + C*Y1
*     Y2 = D + E*X1 + F*Y1
*
*  the present routine generates a new set of coefficients:
*
*     BKWDS(1) = P
*     BKWDS(2) = Q
*     BKWDS(3) = R
*     BKWDS(4) = S
*     BKWDS(5) = T
*     BKWDS(6) = U
*
*  such that:
*
*     X1 = P + Q*X2 + R*Y2
*     Y1 = S + T*X2 + U*Y2
*
*  Two successive calls to sla_INVF will thus deliver a set
*  of coefficients equal to the starting values.
*
*  To comply with the ANSI Fortran standard, FWDS and BKWDS must
*  not be the same array, even though the routine is coded to
*  work on many platforms even if this rule is violated.
*
*  See also sla_FITXY, sla_PXY, sla_XY2XY, sla_DCMPF
*
*  Last revision:   26 December 2004
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

      DOUBLE PRECISION FWDS(6),BKWDS(6)
      INTEGER J

      DOUBLE PRECISION A,B,C,D,E,F,DET



      A=FWDS(1)
      B=FWDS(2)
      C=FWDS(3)
      D=FWDS(4)
      E=FWDS(5)
      F=FWDS(6)
      DET=B*F-C*E
      IF (DET.NE.0D0) THEN
         BKWDS(1)=(C*D-A*F)/DET
         BKWDS(2)=F/DET
         BKWDS(3)=-C/DET
         BKWDS(4)=(A*E-B*D)/DET
         BKWDS(5)=-E/DET
         BKWDS(6)=B/DET
         J=0
      ELSE
         J=-1
      END IF

      END
