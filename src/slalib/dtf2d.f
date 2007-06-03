      SUBROUTINE sla_DTF2D (IHOUR, IMIN, SEC, DAYS, J)
*+
*     - - - - - -
*      D T F 2 D
*     - - - - - -
*
*  Convert hours, minutes, seconds to days (double precision)
*
*  Given:
*     IHOUR       int       hours
*     IMIN        int       minutes
*     SEC         dp        seconds
*
*  Returned:
*     DAYS        dp        interval in days
*     J           int       status:  0 = OK
*                                    1 = IHOUR outside range 0-23
*                                    2 = IMIN outside range 0-59
*                                    3 = SEC outside range 0-59.999...
*
*  Notes:
*
*     1)  The result is computed even if any of the range checks fail.
*
*     2)  The sign must be dealt with outside this routine.
*
*  P.T.Wallace   Starlink   July 1984
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

      INTEGER IHOUR,IMIN
      DOUBLE PRECISION SEC,DAYS
      INTEGER J

*  Seconds per day
      DOUBLE PRECISION D2S
      PARAMETER (D2S=86400D0)



*  Preset status
      J=0

*  Validate sec, min, hour
      IF (SEC.LT.0D0.OR.SEC.GE.60D0) J=3
      IF (IMIN.LT.0.OR.IMIN.GT.59) J=2
      IF (IHOUR.LT.0.OR.IHOUR.GT.23) J=1

*  Compute interval
      DAYS=(60D0*(60D0*DBLE(IHOUR)+DBLE(IMIN))+SEC)/D2S

      END
