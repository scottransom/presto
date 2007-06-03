      SUBROUTINE sla_CALYD (IY, IM, ID, NY, ND, J)
*+
*     - - - - - -
*      C A L Y D
*     - - - - - -
*
*  Gregorian calendar date to year and day in year (in a Julian
*  calendar aligned to the 20th/21st century Gregorian calendar).
*
*  (Includes century default feature:  use sla_CLYD for years
*   before 100AD.)
*
*  Given:
*     IY,IM,ID   int    year, month, day in Gregorian calendar
*                       (year may optionally omit the century)
*  Returned:
*     NY         int    year (re-aligned Julian calendar)
*     ND         int    day in year (1 = January 1st)
*     J          int    status:
*                         0 = OK
*                         1 = bad year (before -4711)
*                         2 = bad month
*                         3 = bad day (but conversion performed)
*
*  Notes:
*
*  1  This routine exists to support the low-precision routines
*     sla_EARTH, sla_MOON and sla_ECOR.
*
*  2  Between 1900 March 1 and 2100 February 28 it returns answers
*     which are consistent with the ordinary Gregorian calendar.
*     Outside this range there will be a discrepancy which increases
*     by one day for every non-leap century year.
*
*  3  Years in the range 50-99 are interpreted as 1950-1999, and
*     years in the range 00-49 are interpreted as 2000-2049.
*
*  Called:  sla_CLYD
*
*  P.T.Wallace   Starlink   23 November 1994
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

      INTEGER IY,IM,ID,NY,ND,J

      INTEGER I



*  Default century if appropriate
      IF (IY.GE.0.AND.IY.LE.49) THEN
         I=IY+2000
      ELSE IF (IY.GE.50.AND.IY.LE.99) THEN
         I=IY+1900
      ELSE
         I=IY
      END IF

*  Perform the conversion
      CALL sla_CLYD(I,IM,ID,NY,ND,J)

      END
