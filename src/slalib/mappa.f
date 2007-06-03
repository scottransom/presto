      SUBROUTINE sla_MAPPA (EQ, DATE, AMPRMS)
*+
*     - - - - - -
*      M A P P A
*     - - - - - -
*
*  Compute star-independent parameters in preparation for
*  conversions between mean place and geocentric apparent place.
*
*  The parameters produced by this routine are required in the
*  parallax, light deflection, aberration, and precession/nutation
*  parts of the mean/apparent transformations.
*
*  The reference frames and timescales used are post IAU 1976.
*
*  Given:
*     EQ       d      epoch of mean equinox to be used (Julian)
*     DATE     d      TDB (JD-2400000.5)
*
*  Returned:
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
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Notes:
*
*  1)  For DATE, the distinction between the required TDB and TT
*      is always negligible.  Moreover, for all but the most
*      critical applications UTC is adequate.
*
*  2)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
*      the mean equinox and equator of epoch EQ.
*
*  3)  The parameters AMPRMS produced by this routine are used by
*      sla_AMPQK, sla_MAPQK and sla_MAPQKZ.
*
*  4)  The accuracy is sub-milliarcsecond, limited by the
*      precession-nutation model (IAU 1976 precession, Shirai &
*      Fukushima 2001 forced nutation and precession corrections).
*
*  5)  A further limit to the accuracy of routines using the parameter
*      array AMPRMS is imposed by the routine sla_EVP, used here to
*      compute the Earth position and velocity by the methods of
*      Stumpff.  The maximum error in the resulting aberration
*      corrections is about 0.3 milliarcsecond.
*
*  Called:
*     sla_EPJ         MDJ to Julian epoch
*     sla_EVP         earth position & velocity
*     sla_DVN         normalize vector
*     sla_PRENUT      precession/nutation matrix
*
*  P.T.Wallace   Starlink   24 October 2003
*
*  Copyright (C) 2003 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION EQ,DATE,AMPRMS(21)

*  Light time for 1 AU (sec)
      DOUBLE PRECISION CR
      PARAMETER (CR=499.004782D0)

*  Gravitational radius of the Sun x 2 (2*mu/c**2, AU)
      DOUBLE PRECISION GR2
      PARAMETER (GR2=2D0*9.87063D-9)

      INTEGER I

      DOUBLE PRECISION EBD(3),EHD(3),EH(3),E,VN(3),VM

      DOUBLE PRECISION sla_EPJ



*  Time interval for proper motion correction
      AMPRMS(1) = sla_EPJ(DATE)-EQ

*  Get Earth barycentric and heliocentric position and velocity
      CALL sla_EVP(DATE,EQ,EBD,AMPRMS(2),EHD,EH)

*  Heliocentric direction of earth (normalized) and modulus
      CALL sla_DVN(EH,AMPRMS(5),E)

*  Light deflection parameter
      AMPRMS(8) = GR2/E

*  Aberration parameters
      DO I=1,3
         AMPRMS(I+8) = EBD(I)*CR
      END DO
      CALL sla_DVN(AMPRMS(9),VN,VM)
      AMPRMS(12) = SQRT(1D0-VM*VM)

*  Precession/nutation matrix
      CALL sla_PRENUT(EQ,DATE,AMPRMS(13))

      END
