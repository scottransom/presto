      SUBROUTINE sla_DCMPF (COEFFS,XZ,YZ,XS,YS,PERP,ORIENT)
*+
*     - - - - - -
*      D C M P F
*     - - - - - -
*
*  Decompose an [X,Y] linear fit into its constituent parameters:
*  zero points, scales, nonperpendicularity and orientation.
*
*  Given:
*     COEFFS  d(6)      transformation coefficients (see note)
*
*  Returned:
*     XZ       d        x zero point
*     YZ       d        y zero point
*     XS       d        x scale
*     YS       d        y scale
*     PERP     d        nonperpendicularity (radians)
*     ORIENT   d        orientation (radians)
*
*  Called:  sla_DRANGE
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
*  the model transforms coordinates [X1,Y1] into coordinates
*  [X2,Y2] as follows:
*
*     X2 = A + B*X1 + C*Y1
*     Y2 = D + E*X1 + F*Y1
*
*  The transformation can be decomposed into four steps:
*
*     1)  Zero points:
*
*             x' = XZ + X1
*             y' = YZ + Y1
*
*     2)  Scales:
*
*             x'' = XS*x'
*             y'' = YS*y'
*
*     3)  Nonperpendicularity:
*
*             x''' = cos(PERP/2)*x'' + sin(PERP/2)*y''
*             y''' = sin(PERP/2)*x'' + cos(PERP/2)*y''
*
*     4)  Orientation:
*
*             X2 = cos(ORIENT)*x''' + sin(ORIENT)*y'''
*             Y2 =-sin(ORIENT)*y''' + cos(ORIENT)*y'''
*
*  See also sla_FITXY, sla_PXY, sla_INVF, sla_XY2XY
*
*  P.T.Wallace   Starlink   19 December 2001
*
*  Copyright (C) 2001 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION COEFFS(6),XZ,YZ,XS,YS,PERP,ORIENT

      DOUBLE PRECISION A,B,C,D,E,F,RB2E2,RC2F2,XSC,YSC,P1,P2,P,WS,WC,
     :                 OR,HP,SHP,CHP,SOR,COR,DET,X0,Y0,sla_DRANGE



*  Copy the six coefficients.
      A = COEFFS(1)
      B = COEFFS(2)
      C = COEFFS(3)
      D = COEFFS(4)
      E = COEFFS(5)
      F = COEFFS(6)

*  Scales.
      RB2E2 = SQRT(B*B+E*E)
      RC2F2 = SQRT(C*C+F*F)
      IF (B*F-C*E.GE.0D0) THEN
         XSC = RB2E2
      ELSE
         B = -B
         E = -E
         XSC = -RB2E2
      END IF
      YSC = RC2F2

*  Non-perpendicularity.
      IF (C.NE.0D0.OR.F.NE.0D0) THEN
         P1 = ATAN2(C,F)
      ELSE
         P1 = 0D0
      END IF
      IF (E.NE.0D0.OR.B.NE.0D0) THEN
         P2 = ATAN2(E,B)
      ELSE
         P2 = 0D0
      END IF
      P = sla_DRANGE(P1+P2)

*  Orientation.
      WS = C*RB2E2-E*RC2F2
      WC = B*RC2F2+F*RB2E2
      IF (WS.NE.0D0.OR.WC.NE.0D0) THEN
         OR = ATAN2(WS,WC)
      ELSE
         OR = 0D0
      END IF

*  Zero points.
      HP = P/2D0
      SHP = SIN(HP)
      CHP = COS(HP)
      SOR = SIN(OR)
      COR = COS(OR)
      DET = XSC*YSC*(CHP+SHP)*(CHP-SHP)
      IF (ABS(DET).GT.0D0) THEN
         X0 = YSC*(A*(CHP*COR-SHP*SOR)-D*(CHP*SOR+SHP*COR))/DET
         Y0 = XSC*(A*(CHP*SOR-SHP*COR)+D*(CHP*COR+SHP*SOR))/DET
      ELSE
         X0 = 0D0
         Y0 = 0D0
      END IF

*  Results.
      XZ = X0
      YZ = Y0
      XS = XSC
      YS = YSC
      PERP = P
      ORIENT = OR

      END
