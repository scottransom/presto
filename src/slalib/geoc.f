      SUBROUTINE sla_GEOC (P, H, R, Z)
*+
*     - - - - -
*      G E O C
*     - - - - -
*
*  Convert geodetic position to geocentric (double precision)
*
*  Given:
*     P     dp     latitude (geodetic, radians)
*     H     dp     height above reference spheroid (geodetic, metres)
*
*  Returned:
*     R     dp     distance from Earth axis (AU)
*     Z     dp     distance from plane of Earth equator (AU)
*
*  Notes:
*
*  1  Geocentric latitude can be obtained by evaluating ATAN2(Z,R).
*
*  2  IAU 1976 constants are used.
*
*  Reference:
*
*     Green,R.M., Spherical Astronomy, CUP 1985, p98.
*
*  Last revision:   22 July 2004
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

      DOUBLE PRECISION P,H,R,Z

*  Earth equatorial radius (metres)
      DOUBLE PRECISION A0
      PARAMETER (A0=6378140D0)

*  Reference spheroid flattening factor and useful function
      DOUBLE PRECISION F,B
      PARAMETER (F=1D0/298.257D0,B=(1D0-F)**2)

*  Astronomical unit in metres
      DOUBLE PRECISION AU
      PARAMETER (AU=1.49597870D11)

      DOUBLE PRECISION SP,CP,C,S



*  Geodetic to geocentric conversion
      SP = SIN(P)
      CP = COS(P)
      C = 1D0/SQRT(CP*CP+B*SP*SP)
      S = B*C
      R = (A0*C+H)*CP/AU
      Z = (A0*S+H)*SP/AU

      END
