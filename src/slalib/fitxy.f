      SUBROUTINE sla_FITXY (ITYPE,NP,XYE,XYM,COEFFS,J)
*+
*     - - - - - -
*      F I T X Y
*     - - - - - -
*
*  Fit a linear model to relate two sets of [X,Y] coordinates.
*
*  Given:
*     ITYPE    i        type of model: 4 or 6 (note 1)
*     NP       i        number of samples (note 2)
*     XYE     d(2,np)   expected [X,Y] for each sample
*     XYM     d(2,np)   measured [X,Y] for each sample
*
*  Returned:
*     COEFFS  d(6)      coefficients of model (note 3)
*     J        i        status:  0 = OK
*                               -1 = illegal ITYPE
*                               -2 = insufficient data
*                               -3 = no solution
*
*  Notes:
*
*  1)  ITYPE, which must be either 4 or 6, selects the type of model
*      fitted.  Both allowed ITYPE values produce a model COEFFS which
*      consists of six coefficients, namely the zero points and, for
*      each of XE and YE, the coefficient of XM and YM.  For ITYPE=6,
*      all six coefficients are independent, modelling squash and shear
*      as well as origin, scale, and orientation.  However, ITYPE=4
*      selects the "solid body rotation" option;  the model COEFFS
*      still consists of the same six coefficients, but now two of
*      them are used twice (appropriately signed).  Origin, scale
*      and orientation are still modelled, but not squash or shear -
*      the units of X and Y have to be the same.
*
*  2)  For NC=4, NP must be at least 2.  For NC=6, NP must be at
*      least 3.
*
*  3)  The model is returned in the array COEFFS.  Naming the
*      elements of COEFFS as follows:
*
*                  COEFFS(1) = A
*                  COEFFS(2) = B
*                  COEFFS(3) = C
*                  COEFFS(4) = D
*                  COEFFS(5) = E
*                  COEFFS(6) = F
*
*      the model is:
*
*            XE = A + B*XM + C*YM
*            YE = D + E*XM + F*YM
*
*      For the "solid body rotation" option (ITYPE=4), the
*      magnitudes of B and F, and of C and E, are equal.  The
*      signs of these coefficients depend on whether there is a
*      sign reversal between XE,YE and XM,YM;  fits are performed
*      with and without a sign reversal and the best one chosen.
*
*  4)  Error status values J=-1 and -2 leave COEFFS unchanged;
*      if J=-3 COEFFS may have been changed.
*
*  See also sla_PXY, sla_INVF, sla_XY2XY, sla_DCMPF
*
*  Called:  sla_DMAT, sla_DMXV
*
*  Last revision:   8 September 2005
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

      INTEGER ITYPE,NP
      DOUBLE PRECISION XYE(2,NP),XYM(2,NP),COEFFS(6)
      INTEGER J

      INTEGER I,JSTAT,IW(4),NSOL
      DOUBLE PRECISION A,B,C,D,AOLD,BOLD,COLD,DOLD,SOLD,
     :                 P,SXE,SXEXM,SXEYM,SYE,SYEYM,SYEXM,SXM,
     :                 SYM,SXMXM,SXMYM,SYMYM,XE,YE,
     :                 XM,YM,V(4),DM3(3,3),DM4(4,4),DET,
     :                 SGN,SXXYY,SXYYX,SX2Y2,SDR2,XR,YR



*  Preset the status
      J=0

*  Variable initializations to avoid compiler warnings
      A = 0D0
      B = 0D0
      C = 0D0
      D = 0D0
      AOLD = 0D0
      BOLD = 0D0
      COLD = 0D0
      DOLD = 0D0
      SOLD = 0D0

*  Float the number of samples
      P=DBLE(NP)

*  Check ITYPE
      IF (ITYPE.EQ.6) THEN

*
*     Six-coefficient linear model
*     ----------------------------

*     Check enough samples
         IF (NP.GE.3) THEN

*     Form summations
         SXE=0D0
         SXEXM=0D0
         SXEYM=0D0
         SYE=0D0
         SYEYM=0D0
         SYEXM=0D0
         SXM=0D0
         SYM=0D0
         SXMXM=0D0
         SXMYM=0D0
         SYMYM=0D0
         DO I=1,NP
            XE=XYE(1,I)
            YE=XYE(2,I)
            XM=XYM(1,I)
            YM=XYM(2,I)
            SXE=SXE+XE
            SXEXM=SXEXM+XE*XM
            SXEYM=SXEYM+XE*YM
            SYE=SYE+YE
            SYEYM=SYEYM+YE*YM
            SYEXM=SYEXM+YE*XM
            SXM=SXM+XM
            SYM=SYM+YM
            SXMXM=SXMXM+XM*XM
            SXMYM=SXMYM+XM*YM
            SYMYM=SYMYM+YM*YM
         END DO

*        Solve for A,B,C in  XE = A + B*XM + C*YM
            V(1)=SXE
            V(2)=SXEXM
            V(3)=SXEYM
            DM3(1,1)=P
            DM3(1,2)=SXM
            DM3(1,3)=SYM
            DM3(2,1)=SXM
            DM3(2,2)=SXMXM
            DM3(2,3)=SXMYM
            DM3(3,1)=SYM
            DM3(3,2)=SXMYM
            DM3(3,3)=SYMYM
            CALL sla_DMAT(3,DM3,V,DET,JSTAT,IW)
            IF (JSTAT.EQ.0) THEN
               DO I=1,3
                  COEFFS(I)=V(I)
               END DO

*           Solve for D,E,F in  YE = D + E*XM + F*YM
               V(1)=SYE
               V(2)=SYEXM
               V(3)=SYEYM
               CALL sla_DMXV(DM3,V,COEFFS(4))

            ELSE

*           No 6-coefficient solution possible
               J=-3

            END IF

         ELSE

*        Insufficient data for 6-coefficient fit
            J=-2

         END IF

      ELSE IF (ITYPE.EQ.4) THEN

*
*     Four-coefficient solid body rotation model
*     ------------------------------------------

*     Check enough samples
         IF (NP.GE.2) THEN

*        Try two solutions, first without then with flip in X
            DO NSOL=1,2
               IF (NSOL.EQ.1) THEN
                  SGN=1D0
               ELSE
                  SGN=-1D0
               END IF

*           Form summations
               SXE=0D0
               SXXYY=0D0
               SXYYX=0D0
               SYE=0D0
               SXM=0D0
               SYM=0D0
               SX2Y2=0D0
               DO I=1,NP
                  XE=XYE(1,I)*SGN
                  YE=XYE(2,I)
                  XM=XYM(1,I)
                  YM=XYM(2,I)
                  SXE=SXE+XE
                  SXXYY=SXXYY+XE*XM+YE*YM
                  SXYYX=SXYYX+XE*YM-YE*XM
                  SYE=SYE+YE
                  SXM=SXM+XM
                  SYM=SYM+YM
                  SX2Y2=SX2Y2+XM*XM+YM*YM
               END DO

*
*           Solve for A,B,C,D in:  +/- XE = A + B*XM - C*YM
*                                    + YE = D + C*XM + B*YM
               V(1)=SXE
               V(2)=SXXYY
               V(3)=SXYYX
               V(4)=SYE
               DM4(1,1)=P
               DM4(1,2)=SXM
               DM4(1,3)=-SYM
               DM4(1,4)=0D0
               DM4(2,1)=SXM
               DM4(2,2)=SX2Y2
               DM4(2,3)=0D0
               DM4(2,4)=SYM
               DM4(3,1)=SYM
               DM4(3,2)=0D0
               DM4(3,3)=-SX2Y2
               DM4(3,4)=-SXM
               DM4(4,1)=0D0
               DM4(4,2)=SYM
               DM4(4,3)=SXM
               DM4(4,4)=P
               CALL sla_DMAT(4,DM4,V,DET,JSTAT,IW)
               IF (JSTAT.EQ.0) THEN
                  A=V(1)
                  B=V(2)
                  C=V(3)
                  D=V(4)

*              Determine sum of radial errors squared
                  SDR2=0D0
                  DO I=1,NP
                     XM=XYM(1,I)
                     YM=XYM(2,I)
                     XR=A+B*XM-C*YM-XYE(1,I)*SGN
                     YR=D+C*XM+B*YM-XYE(2,I)
                     SDR2=SDR2+XR*XR+YR*YR
                  END DO

               ELSE

*              Singular: set flag
                  SDR2=-1D0

               END IF

*           If first pass and non-singular, save variables
               IF (NSOL.EQ.1.AND.JSTAT.EQ.0) THEN
                  AOLD=A
                  BOLD=B
                  COLD=C
                  DOLD=D
                  SOLD=SDR2
               END IF

            END DO

*        Pick the best of the two solutions
            IF (SOLD.GE.0D0.AND.(SOLD.LE.SDR2.OR.NP.EQ.2)) THEN
               COEFFS(1)=AOLD
               COEFFS(2)=BOLD
               COEFFS(3)=-COLD
               COEFFS(4)=DOLD
               COEFFS(5)=COLD
               COEFFS(6)=BOLD
            ELSE IF (JSTAT.EQ.0) THEN
               COEFFS(1)=-A
               COEFFS(2)=-B
               COEFFS(3)=C
               COEFFS(4)=D
               COEFFS(5)=C
               COEFFS(6)=B
            ELSE

*           No 4-coefficient fit possible
               J=-3
            END IF
         ELSE

*        Insufficient data for 4-coefficient fit
            J=-2
         END IF
      ELSE

*     Illegal ITYPE - not 4 or 6
         J=-1
      END IF

      END
