      SUBROUTINE sla_PXY (NP,XYE,XYM,COEFFS,XYP,XRMS,YRMS,RRMS)
*+
*     - - - -
*      P X Y
*     - - - -
*
*  Given arrays of "expected" and "measured" [X,Y] coordinates, and a
*  linear model relating them (as produced by sla_FITXY), compute
*  the array of "predicted" coordinates and the RMS residuals.
*
*  Given:
*     NP       i        number of samples
*     XYE     d(2,np)   expected [X,Y] for each sample
*     XYM     d(2,np)   measured [X,Y] for each sample
*     COEFFS  d(6)      coefficients of model (see below)
*
*  Returned:
*     XYP     d(2,np)   predicted [X,Y] for each sample
*     XRMS     d        RMS in X
*     YRMS     d        RMS in Y
*     RRMS     d        total RMS (vector sum of XRMS and YRMS)
*
*  The model is supplied in the array COEFFS.  Naming the
*  elements of COEFF as follows:
*
*     COEFFS(1) = A
*     COEFFS(2) = B
*     COEFFS(3) = C
*     COEFFS(4) = D
*     COEFFS(5) = E
*     COEFFS(6) = F
*
*  the model is applied thus:
*
*     XP = A + B*XM + C*YM
*     YP = D + E*XM + F*YM
*
*  The residuals are (XP-XE) and (YP-YE).
*
*  If NP is less than or equal to zero, no coordinates are
*  transformed, and the RMS residuals are all zero.
*
*  See also sla_FITXY, sla_INVF, sla_XY2XY, sla_DCMPF
*
*  Called:  sla_XY2XY
*
*  P.T.Wallace   Starlink   22 May 1996
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

      INTEGER NP
      DOUBLE PRECISION XYE(2,NP),XYM(2,NP),COEFFS(6),
     :                 XYP(2,NP),XRMS,YRMS,RRMS

      INTEGER I
      DOUBLE PRECISION SDX2,SDY2,XP,YP,DX,DY,DX2,DY2,P



*  Initialize summations
      SDX2=0D0
      SDY2=0D0

*  Loop by sample
      DO I=1,NP

*     Transform "measured" [X,Y] to "predicted" [X,Y]
         CALL sla_XY2XY(XYM(1,I),XYM(2,I),COEFFS,XP,YP)
         XYP(1,I)=XP
         XYP(2,I)=YP

*     Compute residuals in X and Y, and update summations
         DX=XYE(1,I)-XP
         DY=XYE(2,I)-YP
         DX2=DX*DX
         DY2=DY*DY
         SDX2=SDX2+DX2
         SDY2=SDY2+DY2

*     Next sample
      END DO

*  Compute RMS values
      P=MAX(1D0,DBLE(NP))
      XRMS=SQRT(SDX2/P)
      YRMS=SQRT(SDY2/P)
      RRMS=SQRT(XRMS*XRMS+YRMS*YRMS)

      END
