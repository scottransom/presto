      DOUBLE PRECISION FUNCTION sla_DRANGE (ANGLE)
*+
*     - - - - - - -
*      D R A N G E
*     - - - - - - -
*
*  Normalize angle into range +/- pi  (double precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result (double precision) is ANGLE expressed in the range +/- pi.
*
*  P.T.Wallace   Starlink   23 November 1995
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

      DOUBLE PRECISION ANGLE

      DOUBLE PRECISION DPI,D2PI
      PARAMETER (DPI=3.141592653589793238462643D0)
      PARAMETER (D2PI=6.283185307179586476925287D0)


      sla_DRANGE=MOD(ANGLE,D2PI)
      IF (ABS(sla_DRANGE).GE.DPI)
     :          sla_DRANGE=sla_DRANGE-SIGN(D2PI,ANGLE)

      END
