      SUBROUTINE sla_DJCAL (NDP, DJM, IYMDF, J)
*+
*     - - - - - -
*      D J C A L
*     - - - - - -
*
*  Modified Julian Date to Gregorian Calendar, expressed
*  in a form convenient for formatting messages (namely
*  rounded to a specified precision, and with the fields
*  stored in a single array)
*
*  Given:
*     NDP      i      number of decimal places of days in fraction
*     DJM      d      modified Julian Date (JD-2400000.5)
*
*  Returned:
*     IYMDF    i(4)   year, month, day, fraction in Gregorian
*                     calendar
*     J        i      status:  nonzero = out of range
*
*  Any date after 4701BC March 1 is accepted.
*
*  NDP should be 4 or less if internal overflows are to be avoided
*  on machines which use 32-bit integers.
*
*  The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
*
*  Last revision:   22 July 2004
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

      INTEGER NDP
      DOUBLE PRECISION DJM
      INTEGER IYMDF(4),J

      INTEGER NFD
      DOUBLE PRECISION FD,DF,F,D
      INTEGER JD,N4,ND10


*  Validate.
      IF ( DJM.LE.-2395520D0 .OR. DJM.GE.1D9 ) THEN
         J = -1
      ELSE
         J = 0

*     Denominator of fraction.
         NFD = 10**MAX(NDP,0)
         FD = DBLE(NFD)

*     Round date and express in units of fraction.
         DF = ANINT(DJM*FD)

*     Separate day and fraction.
         F = MOD(DF,FD)
         IF (F.LT.0D0) F = F+FD
         D = (DF-F)/FD

*     Express day in Gregorian calendar.
         JD = NINT(D)+2400001

         N4 = 4*(JD+((2*((4*JD-17918)/146097)*3)/4+1)/2-37)
         ND10 = 10*(MOD(N4-237,1461)/4)+5

         IYMDF(1) = N4/1461-4712
         IYMDF(2) = MOD(ND10/306+2,12)+1
         IYMDF(3) = MOD(ND10,306)/10+1
         IYMDF(4) = NINT(F)

      END IF

      END
