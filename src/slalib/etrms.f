      SUBROUTINE sla_ETRMS (EP, EV)
*+
*     - - - - - -
*      E T R M S
*     - - - - - -
*
*  Compute the E-terms (elliptic component of annual aberration)
*  vector (double precision)
*
*  Given:
*     EP      dp      Besselian epoch
*
*  Returned:
*     EV      dp(3)   E-terms as (dx,dy,dz)
*
*  Note the use of the J2000 aberration constant (20.49552 arcsec).
*  This is a reflection of the fact that the E-terms embodied in
*  existing star catalogues were computed from a variety of
*  aberration constants.  Rather than adopting one of the old
*  constants the latest value is used here.
*
*  References:
*     1  Smith, C.A. et al., 1989.  Astr.J. 97, 265.
*     2  Yallop, B.D. et al., 1989.  Astr.J. 97, 274.
*
*  P.T.Wallace   Starlink   23 August 1996
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

      DOUBLE PRECISION EP,EV(3)

*  Arcseconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=0.484813681109535994D-5)

      DOUBLE PRECISION T,E,E0,P,EK,CP



*  Julian centuries since B1950
      T=(EP-1950D0)*1.00002135903D-2

*  Eccentricity
      E=0.01673011D0-(0.00004193D0+0.000000126D0*T)*T

*  Mean obliquity
      E0=(84404.836D0-(46.8495D0+(0.00319D0+0.00181D0*T)*T)*T)*AS2R

*  Mean longitude of perihelion
      P=(1015489.951D0+(6190.67D0+(1.65D0+0.012D0*T)*T)*T)*AS2R

*  E-terms
      EK=E*20.49552D0*AS2R
      CP=COS(P)
      EV(1)= EK*SIN(P)
      EV(2)=-EK*CP*COS(E0)
      EV(3)=-EK*CP*SIN(E0)

      END
