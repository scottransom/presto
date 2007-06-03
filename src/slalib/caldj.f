      SUBROUTINE sla_CALDJ (IY, IM, ID, DJM, J)
*+
*     - - - - - -
*      C A L D J
*     - - - - - -
*
*  Gregorian Calendar to Modified Julian Date
*
*  (Includes century default feature:  use sla_CLDJ for years
*   before 100AD.)
*
*  Given:
*     IY,IM,ID     int    year, month, day in Gregorian calendar
*
*  Returned:
*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
*     J            int    status:
*                           0 = OK
*                           1 = bad year   (MJD not computed)
*                           2 = bad month  (MJD not computed)
*                           3 = bad day    (MJD computed)
*
*  Acceptable years are 00-49, interpreted as 2000-2049,
*                       50-99,     "       "  1950-1999,
*                       100 upwards, interpreted literally.
*
*  Called:  sla_CLDJ
*
*  P.T.Wallace   Starlink   November 1985
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

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

      INTEGER NY




*  Default century if appropriate
      IF (IY.GE.0.AND.IY.LE.49) THEN
         NY=IY+2000
      ELSE IF (IY.GE.50.AND.IY.LE.99) THEN
         NY=IY+1900
      ELSE
         NY=IY
      END IF

*  Modified Julian Date
      CALL sla_CLDJ(NY,IM,ID,DJM,J)

      END
