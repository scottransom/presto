      SUBROUTINE sla_INTIN (STRING, NSTRT, IRESLT, JFLAG)
*+
*     - - - - - -
*      I N T I N
*     - - - - - -
*
*  Convert free-format input into an integer
*
*  Given:
*     STRING     c     string containing number to be decoded
*     NSTRT      i     pointer to where decoding is to start
*     IRESLT     i     current value of result
*
*  Returned:
*     NSTRT      i      advanced to next number
*     IRESLT     i      result
*     JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
*
*  Called:  sla__IDCHI
*
*  Notes:
*
*     1     The reason INTIN has separate OK status values for +
*           and - is to enable minus zero to be detected.   This is
*           of crucial importance when decoding mixed-radix numbers.
*           For example, an angle expressed as deg, arcmin, arcsec
*           may have a leading minus sign but a zero degrees field.
*
*     2     A TAB is interpreted as a space.
*
*     3     The basic format is the sequence of fields #^, where
*           # is a sign character + or -, and ^ means a string of
*           decimal digits.
*
*     4     Spaces:
*
*             .  Leading spaces are ignored.
*
*             .  Spaces between the sign and the number are allowed.
*
*             .  Trailing spaces are ignored;  the first signifies
*                end of decoding and subsequent ones are skipped.
*
*     5     Delimiters:
*
*             .  Any character other than +,-,0-9 or space may be
*                used to signal the end of the number and terminate
*                decoding.
*
*             .  Comma is recognized by INTIN as a special case;  it
*                is skipped, leaving the pointer on the next character.
*                See 9, below.
*
*     6     The sign is optional.  The default is +.
*
*     7     A "null result" occurs when the string of characters being
*           decoded does not begin with +,- or 0-9, or consists
*           entirely of spaces.  When this condition is detected, JFLAG
*           is set to 1 and IRESLT is left untouched.
*
*     8     NSTRT = 1 for the first character in the string.
*
*     9     On return from INTIN, NSTRT is set ready for the next
*           decode - following trailing blanks and any comma.  If a
*           delimiter other than comma is being used, NSTRT must be
*           incremented before the next call to INTIN, otherwise
*           all subsequent calls will return a null result.
*
*     10    Errors (JFLAG=2) occur when:
*
*             .  there is a + or - but no number;  or
*
*             .  the number is greater than BIG (defined below).
*
*     11    When an error has been detected, NSTRT is left
*           pointing to the character following the last
*           one used before the error came to light.
*
*     12    See also FLOTIN and DFLTIN.
*
*  P.T.Wallace   Starlink   27 April 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
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

      CHARACTER*(*) STRING
      INTEGER NSTRT,IRESLT,JFLAG

*  Maximum allowed value
      DOUBLE PRECISION BIG
      PARAMETER (BIG=2147483647D0)

      INTEGER JPTR,MSIGN,NVEC,J
      DOUBLE PRECISION DRES,DIGIT



*  Current character
      JPTR=NSTRT

*  Set defaults
      DRES=0D0
      MSIGN=1

*  Look for sign
 100  CONTINUE
      CALL sla__IDCHI(STRING,JPTR,NVEC,DIGIT)
      GO TO ( 400, 100,  300,  200, 9110, 9100, 9110),NVEC
*             0-9   SP     +     -     ,   ELSE   END

*  Negative
 200  CONTINUE
      MSIGN=-1

*  Look for first decimal
 300  CONTINUE
      CALL sla__IDCHI(STRING,JPTR,NVEC,DIGIT)
      GO TO ( 400, 300, 9200, 9200, 9200, 9200, 9210),NVEC
*             0-9   SP     +     -     ,   ELSE   END

*  Accept decimals
 400  CONTINUE
      DRES=DRES*1D1+DIGIT

*  Test for overflow
      IF (DRES.GT.BIG) GO TO 9200

*  Look for subsequent decimals
      CALL sla__IDCHI(STRING,JPTR,NVEC,DIGIT)
      GO TO ( 400, 1610, 1600, 1600, 1600, 1600, 1610),NVEC
*             0-9   SP     +     -     ,   ELSE   END

*  Get result & status
 1600 CONTINUE
      JPTR=JPTR-1
 1610 CONTINUE
      J=0
      IF (MSIGN.EQ.1) GO TO 1620
      J=-1
      DRES=-DRES
 1620 CONTINUE
      IRESLT=NINT(DRES)

*  Skip to end of field
 1630 CONTINUE
      CALL sla__IDCHI(STRING,JPTR,NVEC,DIGIT)
      GO TO (1720, 1630, 1720, 1720, 9900, 1720, 9900),NVEC
*             0-9   SP     +     -     ,   ELSE   END

 1720 CONTINUE
      JPTR=JPTR-1
      GO TO 9900

*  Exits

*  Null field
 9100 CONTINUE
      JPTR=JPTR-1
 9110 CONTINUE
      J=1
      GO TO 9900

*  Errors
 9200 CONTINUE
      JPTR=JPTR-1
 9210 CONTINUE
      J=2

*  Return
 9900 CONTINUE
      NSTRT=JPTR
      JFLAG=J

      END
