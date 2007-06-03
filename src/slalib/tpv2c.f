      SUBROUTINE sla_TPV2C (XI, ETA, V, V01, V02, N)
*+
*     - - - - - -
*      T P V 2 C
*     - - - - - -
*
*  Given the tangent-plane coordinates of a star and its direction
*  cosines, determine the direction cosines of the tangent-point.
*
*  (single precision)
*
*  Given:
*     XI,ETA    r       tangent plane coordinates of star
*     V         r(3)    direction cosines of star
*
*  Returned:
*     V01       r(3)    direction cosines of tangent point, solution 1
*     V02       r(3)    direction cosines of tangent point, solution 2
*     N         i       number of solutions:
*                         0 = no solutions returned (note 2)
*                         1 = only the first solution is useful (note 3)
*                         2 = both solutions are useful (note 3)
*
*  Notes:
*
*  1  The vector V must be of unit length or the result will be wrong.
*
*  2  Cases where there is no solution can only arise near the poles.
*     For example, it is clearly impossible for a star at the pole
*     itself to have a non-zero XI value, and hence it is meaningless
*     to ask where the tangent point would have to be.
*
*  3  Also near the poles, cases can arise where there are two useful
*     solutions.  The argument N indicates whether the second of the
*     two solutions returned is useful.  N=1 indicates only one useful
*     solution, the usual case;  under these circumstances, the second
*     solution can be regarded as valid if the vector V02 is interpreted
*     as the "over-the-pole" case.
*
*  4  This routine is the Cartesian equivalent of the routine sla_TPS2C.
*
*  P.T.Wallace   Starlink   5 June 1995
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

      REAL XI,ETA,V(3),V01(3),V02(3)
      INTEGER N

      REAL X,Y,Z,RXY2,XI2,ETA2P1,SDF,R2,R,C


      X=V(1)
      Y=V(2)
      Z=V(3)
      RXY2=X*X+Y*Y
      XI2=XI*XI
      ETA2P1=ETA*ETA+1.0
      SDF=Z*SQRT(XI2+ETA2P1)
      R2=RXY2*ETA2P1-Z*Z*XI2
      IF (R2.GT.0.0) THEN
         R=SQRT(R2)
         C=(SDF*ETA+R)/(ETA2P1*SQRT(RXY2*(R2+XI2)))
         V01(1)=C*(X*R+Y*XI)
         V01(2)=C*(Y*R-X*XI)
         V01(3)=(SDF-ETA*R)/ETA2P1
         R=-R
         C=(SDF*ETA+R)/(ETA2P1*SQRT(RXY2*(R2+XI2)))
         V02(1)=C*(X*R+Y*XI)
         V02(2)=C*(Y*R-X*XI)
         V02(3)=(SDF-ETA*R)/ETA2P1
         IF (ABS(SDF).LT.1.0) THEN
            N=1
         ELSE
            N=2
         END IF
      ELSE
         N=0
      END IF

      END
