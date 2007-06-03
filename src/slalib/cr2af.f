      SUBROUTINE sla_CR2AF (NDP, ANGLE, SIGN, IDMSF)
*+
*     - - - - - -
*      C R 2 A F
*     - - - - - -
*
*  Convert an angle in radians into degrees, arcminutes, arcseconds
*  (single precision)
*
*  Given:
*     NDP       int      number of decimal places of arcseconds
*     ANGLE     real     angle in radians
*
*  Returned:
*     SIGN      char     '+' or '-'
*     IDMSF     int(4)   degrees, arcminutes, arcseconds, fraction
*
*  Notes:
*
*     1)  NDP less than zero is interpreted as zero.
*
*     2)  The largest useful value for NDP is determined by the size of
*         ANGLE, the format of REAL floating-point numbers on the target
*         machine, and the risk of overflowing IDMSF(4).  On some
*         architectures, for ANGLE up to 2pi, the available floating-
*         point precision corresponds roughly to NDP=3.  This is well
*         below the ultimate limit of NDP=9 set by the capacity of a
*         typical 32-bit IDMSF(4).
*
*     3)  The absolute value of ANGLE may exceed 2pi.  In cases where it
*         does not, it is up to the caller to test for and handle the
*         case where ANGLE is very nearly 2pi and rounds up to 360 deg,
*         by testing for IDMSF(1)=360 and setting IDMSF(1-4) to zero.
*
*  Called:  sla_CD2TF
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
      REAL ANGLE
      CHARACTER SIGN*(*)
      INTEGER IDMSF(4)

*  Hours to degrees * radians to turns
      REAL F
      PARAMETER (F=15.0/6.283185307179586476925287)



*  Scale then use days to h,m,s routine
      CALL sla_CD2TF(NDP,ANGLE*F,SIGN,IDMSF)

      END
