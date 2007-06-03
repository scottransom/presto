      DOUBLE PRECISION FUNCTION sla_AIRMAS (ZD)
*+
*     - - - - - - -
*      A I R M A S
*     - - - - - - -
*
*  Air mass at given zenith distance (double precision)
*
*  Given:
*     ZD     d     Observed zenith distance (radians)
*
*  The result is an estimate of the air mass, in units of that
*  at the zenith.
*
*  Notes:
*
*  1)  The "observed" zenith distance referred to above means "as
*      affected by refraction".
*
*  2)  Uses Hardie's (1962) polynomial fit to Bemporad's data for
*      the relative air mass, X, in units of thickness at the zenith
*      as tabulated by Schoenberg (1929). This is adequate for all
*      normal needs as it is accurate to better than 0.1% up to X =
*      6.8 and better than 1% up to X = 10. Bemporad's tabulated
*      values are unlikely to be trustworthy to such accuracy
*      because of variations in density, pressure and other
*      conditions in the atmosphere from those assumed in his work.
*
*  3)  The sign of the ZD is ignored.
*
*  4)  At zenith distances greater than about ZD = 87 degrees the
*      air mass is held constant to avoid arithmetic overflows.
*
*  References:
*     Hardie, R.H., 1962, in "Astronomical Techniques"
*        ed. W.A. Hiltner, University of Chicago Press, p180.
*     Schoenberg, E., 1929, Hdb. d. Ap.,
*        Berlin, Julius Springer, 2, 268.
*
*  Original code by P.W.Hill, St Andrews
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

      DOUBLE PRECISION ZD

      DOUBLE PRECISION SECZM1


      SECZM1 = 1D0/(COS(MIN(1.52D0,ABS(ZD))))-1D0
      sla_AIRMAS = 1D0 + SECZM1*(0.9981833D0
     :             - SECZM1*(0.002875D0 + 0.0008083D0*SECZM1))

      END
