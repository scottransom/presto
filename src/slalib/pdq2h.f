      SUBROUTINE sla_PDQ2H (P, D, Q, H1, J1, H2, J2)
*+
*     - - - - - -
*      P D Q 2 H
*     - - - - - -
*
*  Hour Angle corresponding to a given parallactic angle
*
*  (double precision)
*
*  Given:
*     P       d        latitude
*     D       d        declination
*     Q       d        parallactic angle
*
*  Returned:
*     H1      d        hour angle:  first solution if any
*     J1      i        flag: 0 = solution 1 is valid
*     H2      d        hour angle:  second solution if any
*     J2      i        flag: 0 = solution 2 is valid
*
*  Called:  sla_DRANGE
*
*  P.T.Wallace   Starlink   6 October 1994
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

      DOUBLE PRECISION P,D,Q,H1
      INTEGER J1
      DOUBLE PRECISION H2
      INTEGER J2

      DOUBLE PRECISION DPI
      PARAMETER (DPI=3.141592653589793238462643D0)
      DOUBLE PRECISION D90
      PARAMETER (D90=DPI/2D0)
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-12)
      DOUBLE PRECISION PN,QN,DN,SQ,CQ,SQSD,QT,QB,HPT,T
      DOUBLE PRECISION sla_DRANGE


*  Preset status flags to OK
      J1=0
      J2=0

*  Adjust latitude, declination, parallactic angle to avoid critical values
      PN=sla_DRANGE(P)
      IF (ABS(ABS(PN)-D90).LT.TINY) THEN
         PN=PN-SIGN(TINY,PN)
      ELSE IF (ABS(PN).LT.TINY) THEN
         PN=TINY
      END IF
      QN=sla_DRANGE(Q)
      IF (ABS(ABS(QN)-DPI).LT.TINY) THEN
         QN=QN-SIGN(TINY,QN)
      ELSE IF (ABS(QN).LT.TINY) THEN
         QN=TINY
      END IF
      DN=sla_DRANGE(D)
      IF (ABS(ABS(D)-ABS(P)).LT.TINY) THEN
         DN=DN-SIGN(TINY,DN)
      ELSE IF (ABS(ABS(D)-D90).LT.TINY) THEN
         DN=DN-SIGN(TINY,DN)
      END IF

*  Useful functions
      SQ=SIN(QN)
      CQ=COS(QN)
      SQSD=SQ*SIN(DN)

*  Quotient giving sin(h+t)
      QT=SIN(PN)*SQ*COS(DN)
      QB=COS(PN)*SQRT(CQ*CQ+SQSD*SQSD)

*  Any solutions?
      IF (ABS(QT).LE.QB) THEN

*     Yes: find h+t and t
         HPT=ASIN(QT/QB)
         T=ATAN2(SQSD,CQ)

*     The two solutions
         H1=sla_DRANGE(HPT-T)
         H2=sla_DRANGE(-HPT-(T+DPI))

*     Reject if h and Q different signs
         IF (H1*QN.LT.0D0) J1=-1
         IF (H2*QN.LT.0D0) J2=-1
      ELSE
         J1=-1
         J2=-1
      END IF

      END
