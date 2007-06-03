      DOUBLE PRECISION FUNCTION sla_EQEQX (DATE)
*+
*     - - - - - -
*      E Q E Q X
*     - - - - - -
*
*  Equation of the equinoxes  (IAU 1994, double precision)
*
*  Given:
*     DATE    dp      TDB (loosely ET) as Modified Julian Date
*                                          (JD-2400000.5)
*
*  The result is the equation of the equinoxes (double precision)
*  in radians:
*
*     Greenwich apparent ST = GMST + sla_EQEQX
*
*  References:  IAU Resolution C7, Recommendation 3 (1994)
*               Capitaine, N. & Gontier, A.-M., Astron. Astrophys.,
*               275, 645-650 (1993)
*
*  Called:  sla_NUTC
*
*  Patrick Wallace   Starlink   23 August 1996
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

      DOUBLE PRECISION DATE

*  Turns to arc seconds and arc seconds to radians
      DOUBLE PRECISION T2AS,AS2R
      PARAMETER (T2AS=1296000D0,
     :           AS2R=0.484813681109535994D-5)

      DOUBLE PRECISION T,OM,DPSI,DEPS,EPS0



*  Interval between basic epoch J2000.0 and current epoch (JC)
      T=(DATE-51544.5D0)/36525D0

*  Longitude of the mean ascending node of the lunar orbit on the
*   ecliptic, measured from the mean equinox of date
      OM=AS2R*(450160.280D0+(-5D0*T2AS-482890.539D0
     :         +(7.455D0+0.008D0*T)*T)*T)

*  Nutation
      CALL sla_NUTC(DATE,DPSI,DEPS,EPS0)

*  Equation of the equinoxes
      sla_EQEQX=DPSI*COS(EPS0)+AS2R*(0.00264D0*SIN(OM)+
     :                               0.000063D0*SIN(OM+OM))

      END
