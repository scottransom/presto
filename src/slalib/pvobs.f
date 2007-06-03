      SUBROUTINE sla_PVOBS (P, H, STL, PV)
*+
*     - - - - - -
*      P V O B S
*     - - - - - -
*
*  Position and velocity of an observing station (double precision)
*
*  Given:
*     P     dp     latitude (geodetic, radians)
*     H     dp     height above reference spheroid (geodetic, metres)
*     STL   dp     local apparent sidereal time (radians)
*
*  Returned:
*     PV    dp(6)  position/velocity 6-vector (AU, AU/s, true equator
*                                              and equinox of date)
*
*  Called:  sla_GEOC
*
*  IAU 1976 constants are used.
*
*  P.T.Wallace   Starlink   14 November 1994
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

      DOUBLE PRECISION P,H,STL,PV(6)

      DOUBLE PRECISION R,Z,S,C,V

*  Mean sidereal rate (at J2000) in radians per (UT1) second
      DOUBLE PRECISION SR
      PARAMETER (SR=7.292115855306589D-5)



*  Geodetic to geocentric conversion
      CALL sla_GEOC(P,H,R,Z)

*  Functions of ST
      S=SIN(STL)
      C=COS(STL)

*  Speed
      V=SR*R

*  Position
      PV(1)=R*C
      PV(2)=R*S
      PV(3)=Z

*  Velocity
      PV(4)=-V*S
      PV(5)=V*C
      PV(6)=0D0

      END
