      SUBROUTINE sla_CD2TF (NDP, DAYS, SIGN, IHMSF)
*+
*     - - - - - -
*      C D 2 T F
*     - - - - - -
*
*  Convert an interval in days into hours, minutes, seconds
*
*  (single precision)
*
*  Given:
*     NDP       int      number of decimal places of seconds
*     DAYS      real     interval in days
*
*  Returned:
*     SIGN      char     '+' or '-'
*     IHMSF     int(4)   hours, minutes, seconds, fraction
*
*  Notes:
*
*     1)  NDP less than zero is interpreted as zero.
*
*     2)  The largest useful value for NDP is determined by the size of
*         DAYS, the format of REAL floating-point numbers on the target
*         machine, and the risk of overflowing IHMSF(4).  On some
*         architectures, for DAYS up to 1.0, the available floating-
*         point precision corresponds roughly to NDP=3.  This is well
*         below the ultimate limit of NDP=9 set by the capacity of a
*         typical 32-bit IHMSF(4).
*
*     3)  The absolute value of DAYS may exceed 1.0.  In cases where it
*         does not, it is up to the caller to test for and handle the
*         case where DAYS is very nearly 1.0 and rounds up to 24 hours,
*         by testing for IHMSF(1)=24 and setting IHMSF(1-4) to zero.
*
*  Called:  sla_DD2TF
*
*  Last revision:   26 December 2004
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
      REAL DAYS
      CHARACTER SIGN*(*)
      INTEGER IHMSF(4)



*  Call double precision version
      CALL sla_DD2TF(NDP,DBLE(DAYS),SIGN,IHMSF)

      END
