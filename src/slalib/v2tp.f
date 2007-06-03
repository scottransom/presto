      SUBROUTINE sla_V2TP (V, V0, XI, ETA, J)
*+
*     - - - - -
*      V 2 T P
*     - - - - -
*
*  Given the direction cosines of a star and of the tangent point,
*  determine the star's tangent-plane coordinates.
*
*  (single precision)
*
*  Given:
*     V         r(3)    direction cosines of star
*     V0        r(3)    direction cosines of tangent point
*
*  Returned:
*     XI,ETA    r       tangent plane coordinates of star
*     J         i       status:   0 = OK
*                                 1 = error, star too far from axis
*                                 2 = error, antistar on tangent plane
*                                 3 = error, antistar too far from axis
*
*  Notes:
*
*  1  If vector V0 is not of unit length, or if vector V is of zero
*     length, the results will be wrong.
*
*  2  If V0 points at a pole, the returned XI,ETA will be based on the
*     arbitrary assumption that the RA of the tangent point is zero.
*
*  3  This routine is the Cartesian equivalent of the routine sla_S2TP.
*
*  P.T.Wallace   Starlink   27 November 1996
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

      REAL V(3),V0(3),XI,ETA
      INTEGER J

      REAL X,Y,Z,X0,Y0,Z0,R2,R,W,D

      REAL TINY
      PARAMETER (TINY=1E-6)


      X=V(1)
      Y=V(2)
      Z=V(3)
      X0=V0(1)
      Y0=V0(2)
      Z0=V0(3)
      R2=X0*X0+Y0*Y0
      R=SQRT(R2)
      IF (R.EQ.0.0) THEN
         R=1E-20
         X0=R
      END IF
      W=X*X0+Y*Y0
      D=W+Z*Z0
      IF (D.GT.TINY) THEN
         J=0
      ELSE IF (D.GE.0.0) THEN
         J=1
         D=TINY
      ELSE IF (D.GT.-TINY) THEN
         J=2
         D=-TINY
      ELSE
         J=3
      END IF
      D=D*R
      XI=(Y*X0-X*Y0)/D
      ETA=(Z*R2-Z0*W)/D

      END
