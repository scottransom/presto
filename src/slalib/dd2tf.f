      SUBROUTINE sla_DD2TF (NDP, DAYS, SIGN, IHMSF)
*+
*     - - - - - -
*      D D 2 T F
*     - - - - - -
*
*  Convert an interval in days into hours, minutes, seconds
*  (double precision)
*
*  Given:
*     NDP      i      number of decimal places of seconds
*     DAYS     d      interval in days
*
*  Returned:
*     SIGN     c      '+' or '-'
*     IHMSF    i(4)   hours, minutes, seconds, fraction
*
*  Notes:
*
*     1)  NDP less than zero is interpreted as zero.
*
*     2)  The largest useful value for NDP is determined by the size
*         of DAYS, the format of DOUBLE PRECISION floating-point numbers
*         on the target machine, and the risk of overflowing IHMSF(4).
*         On some architectures, for DAYS up to 1D0, the available
*         floating-point precision corresponds roughly to NDP=12.
*         However, the practical limit is NDP=9, set by the capacity of
*         a typical 32-bit IHMSF(4).
*
*     3)  The absolute value of DAYS may exceed 1D0.  In cases where it
*         does not, it is up to the caller to test for and handle the
*         case where DAYS is very nearly 1D0 and rounds up to 24 hours,
*         by testing for IHMSF(1)=24 and setting IHMSF(1-4) to zero.
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
      DOUBLE PRECISION DAYS
      CHARACTER SIGN*(*)
      INTEGER IHMSF(4)

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER (D2S=86400D0)

      INTEGER NRS,N
      DOUBLE PRECISION RS,RM,RH,A,AH,AM,AS,AF



*  Handle sign
      IF (DAYS.GE.0D0) THEN
         SIGN='+'
      ELSE
         SIGN='-'
      END IF

*  Field units in terms of least significant figure
      NRS=1
      DO N=1,NDP
         NRS=NRS*10
      END DO
      RS=DBLE(NRS)
      RM=RS*60D0
      RH=RM*60D0

*  Round interval and express in smallest units required
      A=ANINT(RS*D2S*ABS(DAYS))

*  Separate into fields
      AH=AINT(A/RH)
      A=A-AH*RH
      AM=AINT(A/RM)
      A=A-AM*RM
      AS=AINT(A/RS)
      AF=A-AS*RS

*  Return results
      IHMSF(1)=MAX(NINT(AH),0)
      IHMSF(2)=MAX(MIN(NINT(AM),59),0)
      IHMSF(3)=MAX(MIN(NINT(AS),59),0)
      IHMSF(4)=MAX(NINT(MIN(AF,RS-1D0)),0)

      END
