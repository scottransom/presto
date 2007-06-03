      REAL FUNCTION sla_SEPV (V1, V2)
*+
*     - - - - -
*      S E P V
*     - - - - -
*
*  Angle between two vectors.
*
*  (single precision)
*
*  Given:
*     V1      r(3)    first vector
*     V2      r(3)    second vector
*
*  The result is the angle, in radians, between the two vectors.  It
*  is always positive.
*
*  Notes:
*
*  1  There is no requirement for the vectors to be unit length.
*
*  2  If either vector is null, zero is returned.
*
*  3  The simplest formulation would use dot product alone.  However,
*     this would reduce the accuracy for angles near zero and pi.  The
*     algorithm uses both cross product and dot product, which maintains
*     accuracy for all sizes of angle.
*
*  Called:  sla_DSEPV
*
*  Last revision:   7 May 2000
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

      REAL V1(3),V2(3)

      INTEGER I
      DOUBLE PRECISION DV1(3),DV2(3)
      DOUBLE PRECISION sla_DSEPV



*  Use double precision version.
      DO I=1,3
         DV1(I) = DBLE(V1(I))
         DV2(I) = DBLE(V2(I))
      END DO
      sla_SEPV = REAL(sla_DSEPV(DV1,DV2))

      END
