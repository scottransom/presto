      SUBROUTINE sla_CLYD (IY, IM, ID, NY, ND, JSTAT)
*+
*     - - - - -
*      C L Y D
*     - - - - -
*
*  Gregorian calendar to year and day in year (in a Julian calendar
*  aligned to the 20th/21st century Gregorian calendar).
*
*  Given:
*     IY,IM,ID   i    year, month, day in Gregorian calendar
*
*  Returned:
*     NY         i    year (re-aligned Julian calendar)
*     ND         i    day in year (1 = January 1st)
*     JSTAT      i    status:
*                       0 = OK
*                       1 = bad year (before -4711)
*                       2 = bad month
*                       3 = bad day (but conversion performed)
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
*  3  The essence of the algorithm is first to express the Gregorian
*     date as a Julian Day Number and then to convert this back to
*     a Julian calendar date, with day-in-year instead of month and
*     day.  See 12.92-1 and 12.95-1 in the reference.
*
*  Reference:  Explanatory Supplement to the Astronomical Almanac,
*              ed P.K.Seidelmann, University Science Books (1992),
*              p604-606.
*
*  P.T.Wallace   Starlink   26 November 1994
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

      INTEGER IY,IM,ID,NY,ND,JSTAT

      INTEGER I,J,K,L,N

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB/31,28,31,30,31,30,31,31,30,31,30,31/



*  Preset status
      JSTAT=0

*  Validate year
      IF (IY.GE.-4711) THEN

*     Validate month
         IF (IM.GE.1.AND.IM.LE.12) THEN

*        Allow for (Gregorian) leap year
            IF (MOD(IY,4).EQ.0.AND.
     :         (MOD(IY,100).NE.0.OR.MOD(IY,400).EQ.0)) THEN
               MTAB(2)=29
            ELSE
               MTAB(2)=28
            END IF

*        Validate day
            IF (ID.LT.1.OR.ID.GT.MTAB(IM)) JSTAT=3

*        Perform the conversion
            I=(14-IM)/12
            K=IY-I
            J=(1461*(K+4800))/4+(367*(IM-2+12*I))/12
     :        -(3*((K+4900)/100))/4+ID-30660
            K=(J-1)/1461
            L=J-1461*K
            N=(L-1)/365-L/1461
            J=((80*(L-365*N+30))/2447)/11
            I=N+J
            ND=59+L-365*I+((4-N)/4)*(1-J)
            NY=4*K+I-4716

*        Bad month
         ELSE
            JSTAT=2
         END IF
      ELSE

*     Bad year
         JSTAT=1
      END IF

      END
