      REAL FUNCTION sla_RANORM (ANGLE)
*+
*     - - - - - - -
*      R A N O R M
*     - - - - - - -
*
*  Normalize angle into range 0-2 pi  (single precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result is ANGLE expressed in the range 0-2 pi (single
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

      REAL A2PI
      PARAMETER (A2PI=6.283185307179586476925287)


      sla_RANORM=MOD(ANGLE,A2PI)
      IF (sla_RANORM.LT.0.0) sla_RANORM=sla_RANORM+A2PI

      END
