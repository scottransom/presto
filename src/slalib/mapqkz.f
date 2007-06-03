      SUBROUTINE sla_MAPQKZ (RM, DM, AMPRMS, RA, DA)
*+
*     - - - - - - -
*      M A P Q K Z
*     - - - - - - -
*
*  Quick mean to apparent place:  transform a star RA,Dec from
*  mean place to geocentric apparent place, given the
*  star-independent parameters, and assuming zero parallax
*  and proper motion.
*
*  Use of this routine is appropriate when efficiency is important
*  and where many star positions, all with parallax and proper
*  motion either zero or already allowed for, and all referred to
*  the same equator and equinox, are to be transformed for one
*  epoch.  The star-independent parameters can be obtained by
*  calling the sla_MAPPA routine.
*
*  The corresponding routine for the case of non-zero parallax
*  and proper motion is sla_MAPQK.
*
*  The reference frames and timescales used are post IAU 1976.
*
*  Given:
*     RM,DM    d      mean RA,Dec (rad)
*     AMPRMS   d(21)  star-independent mean-to-apparent parameters:
*
*       (1-4)    not used
*       (5-7)    heliocentric direction of the Earth (unit vector)
*       (8)      (grav rad Sun)*2/(Sun-Earth distance)
*       (9-11)   ABV: barycentric Earth velocity in units of c
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
*  1)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to the
*      mean equinox and equator of epoch EQ.
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
*  Called:  sla_DCS2C, sla_DVDV, sla_DMXV, sla_DCC2S, sla_DRANRM
*
*  P.T.Wallace   Starlink   18 March 1999
*
*  Copyright (C) 1999 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION RM,DM,AMPRMS(21),RA,DA

      INTEGER I

      DOUBLE PRECISION GR2E,AB1,EHN(3),ABV(3),
     :                 P(3),PDE,PDEP1,W,P1(3),P1DV,
     :                 P1DVP1,P2(3),P3(3)

      DOUBLE PRECISION sla_DVDV,sla_DRANRM




*  Unpack scalar and vector parameters
      GR2E = AMPRMS(8)
      AB1 = AMPRMS(12)
      DO I=1,3
         EHN(I) = AMPRMS(I+4)
         ABV(I) = AMPRMS(I+8)
      END DO

*  Spherical to x,y,z
      CALL sla_DCS2C(RM,DM,P)

*  Light deflection
      PDE = sla_DVDV(P,EHN)
      PDEP1 = PDE+1D0
      W = GR2E/MAX(PDEP1,1D-5)
      DO I=1,3
         P1(I) = P(I)+W*(EHN(I)-PDE*P(I))
      END DO

*  Aberration
      P1DV = sla_DVDV(P1,ABV)
      P1DVP1 = P1DV+1D0
      W = 1D0+P1DV/(AB1+1D0)
      DO I=1,3
         P2(I) = (AB1*P1(I)+W*ABV(I))/P1DVP1
      END DO

*  Precession and nutation
      CALL sla_DMXV(AMPRMS(13),P2,P3)

*  Geocentric apparent RA,Dec
      CALL sla_DCC2S(P3,RA,DA)
      RA = sla_DRANRM(RA)

      END
