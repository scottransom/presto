      SUBROUTINE sla_MAPQK (RM, DM, PR, PD, PX, RV, AMPRMS, RA, DA)
*+
*     - - - - - -
*      M A P Q K
*     - - - - - -
*
*  Quick mean to apparent place:  transform a star RA,Dec from
*  mean place to geocentric apparent place, given the
*  star-independent parameters.
*
*  Use of this routine is appropriate when efficiency is important
*  and where many star positions, all referred to the same equator
*  and equinox, are to be transformed for one epoch.  The
*  star-independent parameters can be obtained by calling the
*  sla_MAPPA routine.
*
*  If the parallax and proper motions are zero the sla_MAPQKZ
*  routine can be used instead.
*
*  The reference frames and timescales used are post IAU 1976.
*
*  Given:
*     RM,DM    d      mean RA,Dec (rad)
*     PR,PD    d      proper motions:  RA,Dec changes per Julian year
*     PX       d      parallax (arcsec)
*     RV       d      radial velocity (km/sec, +ve if receding)
*
*     AMPRMS   d(21)  star-independent mean-to-apparent parameters:
*
*       (1)      time interval for proper motion (Julian years)
*       (2-4)    barycentric position of the Earth (AU)
*       (5-7)    heliocentric direction of the Earth (unit vector)
*       (8)      (grav rad Sun)*2/(Sun-Earth distance)
*       (9-11)   barycentric Earth velocity in units of c
*       (12)     sqrt(1-v**2) where v=modulus(ABV)
*       (13-21)  precession/nutation (3,3) matrix
*
*  Returned:
*     RA,DA    d      apparent RA,Dec (rad)
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Notes:
*
*  1)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
*      the mean equinox and equator of epoch EQ.
*
*  2)  Strictly speaking, the routine is not valid for solar-system
*      sources, though the error will usually be extremely small.
*      However, to prevent gross errors in the case where the
*      position of the Sun is specified, the gravitational
*      deflection term is restrained within about 920 arcsec of the
*      centre of the Sun's disc.  The term has a maximum value of
*      about 1.85 arcsec at this radius, and decreases to zero as
*      the centre of the disc is approached.
*
*  Called:
*     sla_DCS2C       spherical to Cartesian
*     sla_DVDV        dot product
*     sla_DMXV        matrix x vector
*     sla_DCC2S       Cartesian to spherical
*     sla_DRANRM      normalize angle 0-2Pi
*
*  P.T.Wallace   Starlink   15 January 2000
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

      DOUBLE PRECISION RM,DM,PR,PD,PX,RV,AMPRMS(21),RA,DA

*  Arc seconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=0.484813681109535994D-5)

*  Km/s to AU/year
      DOUBLE PRECISION VF
      PARAMETER (VF=0.21094502D0)

      INTEGER I

      DOUBLE PRECISION PMT,GR2E,AB1,EB(3),EHN(3),ABV(3),
     :                 Q(3),PXR,W,EM(3),P(3),PN(3),PDE,PDEP1,
     :                 P1(3),P1DV,P2(3),P3(3)

      DOUBLE PRECISION sla_DVDV,sla_DRANRM



*  Unpack scalar and vector parameters
      PMT = AMPRMS(1)
      GR2E = AMPRMS(8)
      AB1 = AMPRMS(12)
      DO I=1,3
         EB(I) = AMPRMS(I+1)
         EHN(I) = AMPRMS(I+4)
         ABV(I) = AMPRMS(I+8)
      END DO

*  Spherical to x,y,z
      CALL sla_DCS2C(RM,DM,Q)

*  Space motion (radians per year)
      PXR = PX*AS2R
      W = VF*RV*PXR
      EM(1) = -PR*Q(2)-PD*COS(RM)*SIN(DM)+W*Q(1)
      EM(2) =  PR*Q(1)-PD*SIN(RM)*SIN(DM)+W*Q(2)
      EM(3) =          PD*COS(DM)        +W*Q(3)

*  Geocentric direction of star (normalized)
      DO I=1,3
         P(I) = Q(I)+PMT*EM(I)-PXR*EB(I)
      END DO
      CALL sla_DVN(P,PN,W)

*  Light deflection (restrained within the Sun's disc)
      PDE = sla_DVDV(PN,EHN)
      PDEP1 = PDE+1D0
      W = GR2E/MAX(PDEP1,1D-5)
      DO I=1,3
         P1(I) = PN(I)+W*(EHN(I)-PDE*PN(I))
      END DO

*  Aberration (normalization omitted)
      P1DV = sla_DVDV(P1,ABV)
      W = 1D0+P1DV/(AB1+1D0)
      DO I=1,3
         P2(I) = AB1*P1(I)+W*ABV(I)
      END DO

*  Precession and nutation
      CALL sla_DMXV(AMPRMS(13),P2,P3)

*  Geocentric apparent RA,Dec
      CALL sla_DCC2S(P3,RA,DA)
      RA = sla_DRANRM(RA)

      END
