      SUBROUTINE sla_NUT (DATE, RMATN)
*+
*     - - - -
*      N U T
*     - - - -
*
*  Form the matrix of nutation for a given date - Shirai & Fukushima
*  2001 theory (double precision)
*
*  Reference:
*     Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
*
*  Given:
*     DATE    d          TDB (loosely ET) as Modified Julian Date
*                                           (=JD-2400000.5)
*  Returned:
*     RMATN   d(3,3)     nutation matrix
*
*  Notes:
*
*  1  The matrix is in the sense  v(true) = rmatn * v(mean) .
*     where v(true) is the star vector relative to the true equator and
*     equinox of date and v(mean) is the star vector relative to the
*     mean equator and equinox of date.
*
*  2  The matrix represents forced nutation (but not free core
*     nutation) plus corrections to the IAU~1976 precession model.
*
*  3  Earth attitude predictions made by combining the present nutation
*     matrix with IAU~1976 precession are accurate to 1~mas (with
*     respect to the ICRS) for a few decades around 2000.
*
*  4  The distinction between the required TDB and TT is always
*     negligible.  Moreover, for all but the most critical applications
*     UTC is adequate.
*
*  Called:   sla_NUTC, sla_DEULER
*
*  Last revision:   1 December 2005
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

      DOUBLE PRECISION DATE,RMATN(3,3)

      DOUBLE PRECISION DPSI,DEPS,EPS0



*  Nutation components and mean obliquity
      CALL sla_NUTC(DATE,DPSI,DEPS,EPS0)

*  Rotation matrix
      CALL sla_DEULER('XZX',EPS0,-DPSI,-(EPS0+DEPS),RMATN)

      END
