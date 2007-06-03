      SUBROUTINE sla_KBJ (JB, E, K, J)
*+
*     - - - -
*      K B J
*     - - - -
*
*  Select epoch prefix 'B' or 'J'
*
*  Given:
*     JB     int     sla_DBJIN prefix status:  0=none, 1='B', 2='J'
*     E      dp      epoch - Besselian or Julian
*
*  Returned:
*     K      char    'B' or 'J'
*     J      int     status:  0=OK
*
*  If JB=0, B is assumed for E < 1984D0, otherwise J.
*
*  P.T.Wallace   Starlink   31 July 1989
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

      INTEGER JB
      DOUBLE PRECISION E
      CHARACTER K*(*)
      INTEGER J

*  Preset status
      J=0

*  If prefix given expressly, use it
      IF (JB.EQ.1) THEN
         K='B'
      ELSE IF (JB.EQ.2) THEN
         K='J'

*  If no prefix, examine the epoch
      ELSE IF (JB.EQ.0) THEN

*     If epoch is pre-1984.0, assume Besselian;  otherwise Julian
         IF (E.LT.1984D0) THEN
            K='B'
         ELSE
            K='J'
         END IF

*  If illegal prefix, return error status
      ELSE
         K=' '
         J=1
      END IF

      END
