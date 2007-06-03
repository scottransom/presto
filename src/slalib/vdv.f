      REAL FUNCTION sla_VDV (VA, VB)
*+
*     - - - -
*      V D V
*     - - - -
*
*  Scalar product of two 3-vectors  (single precision)
*
*  Given:
*      VA      real(3)     first vector
*      VB      real(3)     second vector
*
*  The result is the scalar product VA.VB (single precision)
*
*  P.T.Wallace   Starlink   November 1984
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

      REAL VA(3),VB(3)


      sla_VDV=VA(1)*VB(1)+VA(2)*VB(2)+VA(3)*VB(3)

      END
