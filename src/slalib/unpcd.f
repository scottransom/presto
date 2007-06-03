      SUBROUTINE sla_UNPCD ( DISCO, X, Y )
*+
*     - - - - - -
*      U N P C D
*     - - - - - -
*
*  Remove pincushion/barrel distortion from a distorted [x,y] to give
*  tangent-plane [x,y].
*
*  Given:
*     DISCO    d      pincushion/barrel distortion coefficient
*     X,Y      d      distorted coordinates
*
*  Returned:
*     X,Y      d      tangent-plane coordinates
*
*  Notes:
*
*  1)  The distortion is of the form RP = R*(1+C*R^2), where R is
*      the radial distance from the tangent point, C is the DISCO
*      argument, and RP is the radial distance in the presence of
*      the distortion.
*
*  2)  For pincushion distortion, C is +ve;  for barrel distortion,
*      C is -ve.
*
*  3)  For X,Y in "radians" - units of one projection radius,
*      which in the case of a photograph is the focal length of
*      the camera - the following DISCO values apply:
*
*          Geometry          DISCO
*
*          astrograph         0.0
*          Schmidt           -0.3333
*          AAT PF doublet  +147.069
*          AAT PF triplet  +178.585
*          AAT f/8          +21.20
*          JKT f/8          +13.32
*
*  4)  The present routine is a rigorous inverse of the companion
*      routine sla_PCD.  The expression for RP in Note 1 is rewritten
*      in the form x^3+a*x+b=0 and solved by standard techniques.
*
*  5)  Cases where the cubic has multiple real roots can sometimes
*      occur, corresponding to extreme instances of barrel distortion
*      where up to three different undistorted [X,Y]s all produce the
*      same distorted [X,Y].  However, only one solution is returned,
*      the one that produces the smallest change in [X,Y].
*
*  P.T.Wallace   Starlink   3 September 2000
*
*  Copyright (C) 2000 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION DISCO,X,Y

      DOUBLE PRECISION THIRD
      PARAMETER (THIRD=1D0/3D0)
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925286766559D0)

      DOUBLE PRECISION RP,Q,R,D,W,S,T,F,C,T3,F1,F2,F3,W1,W2,W3



*  Distance of the point from the origin.
      RP = SQRT(X*X+Y*Y)

*  If zero, or if no distortion, no action is necessary.
      IF (RP.NE.0D0.AND.DISCO.NE.0D0) THEN

*     Begin algebraic solution.
         Q = 1D0/(3D0*DISCO)
         R = RP/(2D0*DISCO)
         W = Q*Q*Q+R*R

*     Continue if one real root, or three of which only one is positive.
         IF (W.GE.0D0) THEN
            D = SQRT(W)
            W = R+D
            S = SIGN(ABS(W)**THIRD,W)
            W = R-D
            T = SIGN((ABS(W))**THIRD,W)
            F = S+T
         ELSE

*        Three different real roots:  use geometrical method instead.
            W = 2D0/SQRT(-3D0*DISCO)
            C = 4D0*RP/(DISCO*W*W*W)
            S = SQRT(1D0-MIN(C*C,1D0))
            T3 = ATAN2(S,C)

*        The three solutions.
            F1 = W*COS((D2PI-T3)/3D0)
            F2 = W*COS((T3)/3D0)
            F3 = W*COS((D2PI+T3)/3D0)

*        Pick the one that moves [X,Y] least.
            W1 = ABS(F1-RP)
            W2 = ABS(F2-RP)
            W3 = ABS(F3-RP)
            IF (W1.LT.W2) THEN
               IF (W1.LT.W3) THEN
                  F = F1
               ELSE
                  F = F3
               END IF
            ELSE
               IF (W2.LT.W3) THEN
                  F = F2
               ELSE
                  F = F3
               END IF
            END IF

         END IF

*     Remove the distortion.
         F = F/RP
         X = F*X
         Y = F*Y

      END IF

      END
