      SUBROUTINE sla_XY2XY (X1,Y1,COEFFS,X2,Y2)
*+
*     - - - - - -
*      X Y 2 X Y
*     - - - - - -
*
*  Transform one [X,Y] into another using a linear model of the type
*  produced by the sla_FITXY routine.
*
*  Given:
*     X1       d        x-coordinate
*     Y1       d        y-coordinate
*     COEFFS  d(6)      transformation coefficients (see note)
*
*  Returned:
*     X2       d        x-coordinate
*     Y2       d        y-coordinate
*
*  The model relates two sets of [X,Y] coordinates as follows.
*  Naming the elements of COEFFS:
*
*     COEFFS(1) = A
*     COEFFS(2) = B
*     COEFFS(3) = C
*     COEFFS(4) = D
*     COEFFS(5) = E
*     COEFFS(6) = F
*
*  the present routine performs the transformation:
*
*     X2 = A + B*X1 + C*Y1
*     Y2 = D + E*X1 + F*Y1
*
*  See also sla_FITXY, sla_PXY, sla_INVF, sla_DCMPF
*
*  P.T.Wallace   Starlink   5 December 1994
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

      DOUBLE PRECISION X1,Y1,COEFFS(6),X2,Y2


      X2=COEFFS(1)+COEFFS(2)*X1+COEFFS(3)*Y1
      Y2=COEFFS(4)+COEFFS(5)*X1+COEFFS(6)*Y1

      END
