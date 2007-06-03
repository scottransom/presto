      SUBROUTINE sla_CC62S (V, A, B, R, AD, BD, RD)
*+
*     - - - - - -
*      C C 6 2 S
*     - - - - - -
*
*  Conversion of position & velocity in Cartesian coordinates
*  to spherical coordinates (single precision)
*
*  Given:
*     V      r(6)   Cartesian position & velocity vector
*
*  Returned:
*     A      r      longitude (radians)
*     B      r      latitude (radians)
*     R      r      radial coordinate
*     AD     r      longitude derivative (radians per unit time)
*     BD     r      latitude derivative (radians per unit time)
*     RD     r      radial derivative
*
*  P.T.Wallace   Starlink   28 April 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
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

      REAL V(6),A,B,R,AD,BD,RD

      REAL X,Y,Z,XD,YD,ZD,RXY2,RXY,R2,XYP



*  Components of position/velocity vector
      X=V(1)
      Y=V(2)
      Z=V(3)
      XD=V(4)
      YD=V(5)
      ZD=V(6)

*  Component of R in XY plane squared
      RXY2=X*X+Y*Y

*  Modulus squared
      R2=RXY2+Z*Z

*  Protection against null vector
      IF (R2.EQ.0.0) THEN
         X=XD
         Y=YD
         Z=ZD
         RXY2=X*X+Y*Y
         R2=RXY2+Z*Z
      END IF

*  Position and velocity in spherical coordinates
      RXY=SQRT(RXY2)
      XYP=X*XD+Y*YD
      IF (RXY2.NE.0.0) THEN
         A=ATAN2(Y,X)
         B=ATAN2(Z,RXY)
         AD=(X*YD-Y*XD)/RXY2
         BD=(ZD*RXY2-Z*XYP)/(R2*RXY)
      ELSE
         A=0.0
         IF (Z.NE.0.0) THEN
            B=ATAN2(Z,RXY)
         ELSE
            B=0.0
         END IF
         AD=0.0
         BD=0.0
      END IF
      R=SQRT(R2)
      IF (R.NE.0.0) THEN
         RD=(XYP+Z*ZD)/R
      ELSE
         RD=0.0
      END IF

      END
