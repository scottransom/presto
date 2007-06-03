      SUBROUTINE sla_AMPQK (RA, DA, AMPRMS, RM, DM)
*+
*     - - - - - -
*      A M P Q K
*     - - - - - -
*
*  Convert star RA,Dec from geocentric apparent to mean place
*
*  The mean coordinate system is the post IAU 1976 system,
*  loosely called FK5.
*
*  Use of this routine is appropriate when efficiency is important
*  and where many star positions are all to be transformed for
*  one epoch and equinox.  The star-independent parameters can be
*  obtained by calling the sla_MAPPA routine.
*
*  Given:
*     RA       d      apparent RA (radians)
*     DA       d      apparent Dec (radians)
*
*     AMPRMS   d(21)  star-independent mean-to-apparent parameters:
*
*       (1)      time interval for proper motion (Julian years)
*       (2-4)    barycentric position of the Earth (AU)
*       (5-7)    heliocentric direction of the Earth (unit vector)
*       (8)      (grav rad Sun)*2/(Sun-Earth distance)
*       (9-11)   ABV: barycentric Earth velocity in units of c
*       (12)     sqrt(1-v**2) where v=modulus(ABV)
*       (13-21)  precession/nutation (3,3) matrix
*
*  Returned:
*     RM       d      mean RA (radians)
*     DM       d      mean Dec (radians)
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Note:
*
*     Iterative techniques are used for the aberration and
*     light deflection corrections so that the routines
*     sla_AMP (or sla_AMPQK) and sla_MAP (or sla_MAPQK) are
*     accurate inverses;  even at the edge of the Sun's disc
*     the discrepancy is only about 1 nanoarcsecond.
*
*  Called:  sla_DCS2C, sla_DIMXV, sla_DVDV, sla_DVN, sla_DCC2S,
*           sla_DRANRM
*
*  P.T.Wallace   Starlink   7 May 2000
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

      DOUBLE PRECISION RA,DA,AMPRMS(21),RM,DM

      INTEGER I,J

      DOUBLE PRECISION GR2E,AB1,EHN(3),ABV(3),P3(3),P2(3),
     :                 AB1P1,P1DV,P1DVP1,P1(3),W,PDE,PDEP1,P(3)

      DOUBLE PRECISION sla_DVDV,sla_DRANRM



*  Unpack scalar and vector parameters
      GR2E = AMPRMS(8)
      AB1 = AMPRMS(12)
      DO I=1,3
         EHN(I) = AMPRMS(I+4)
         ABV(I) = AMPRMS(I+8)
      END DO

*  Apparent RA,Dec to Cartesian
      CALL sla_DCS2C(RA,DA,P3)

*  Precession and nutation
      CALL sla_DIMXV(AMPRMS(13),P3,P2)

*  Aberration
      AB1P1 = AB1+1D0
      DO I=1,3
         P1(I) = P2(I)
      END DO
      DO J=1,2
         P1DV = sla_DVDV(P1,ABV)
         P1DVP1 = 1D0+P1DV
         W = 1D0+P1DV/AB1P1
         DO I=1,3
            P1(I) = (P1DVP1*P2(I)-W*ABV(I))/AB1
         END DO
         CALL sla_DVN(P1,P3,W)
         DO I=1,3
            P1(I) = P3(I)
         END DO
      END DO

*  Light deflection
      DO I=1,3
         P(I) = P1(I)
      END DO
      DO J=1,5
         PDE = sla_DVDV(P,EHN)
         PDEP1 = 1D0+PDE
         W = PDEP1-GR2E*PDE
         DO I=1,3
            P(I) = (PDEP1*P1(I)-GR2E*EHN(I))/W
         END DO
         CALL sla_DVN(P,P2,W)
         DO I=1,3
            P(I) = P2(I)
         END DO
      END DO

*  Mean RA,Dec
      CALL sla_DCC2S(P,RM,DM)
      RM = sla_DRANRM(RM)

      END
