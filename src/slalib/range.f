      REAL FUNCTION sla_RANGE (ANGLE)
*+
*     - - - - - -
*      R A N G E
*     - - - - - -
*
*  Normalize angle into range +/- pi  (single precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result is ANGLE expressed in the +/- pi (single
*  precision).
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

      REAL ANGLE

      REAL API,A2PI
      PARAMETER (API=3.141592653589793238462643)
      PARAMETER (A2PI=6.283185307179586476925287)


      sla_RANGE=MOD(ANGLE,A2PI)
      IF (ABS(sla_RANGE).GE.API)
     :          sla_RANGE=sla_RANGE-SIGN(A2PI,ANGLE)

      END
