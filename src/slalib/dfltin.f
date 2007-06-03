      SUBROUTINE sla_DFLTIN (STRING, NSTRT, DRESLT, JFLAG)
*+
*     - - - - - - -
*      D F L T I N
*     - - - - - - -
*
*  Convert free-format input into double precision floating point
*
*  Given:
*     STRING     c     string containing number to be decoded
*     NSTRT      i     pointer to where decoding is to start
*     DRESLT     d     current value of result
*
*  Returned:
*     NSTRT      i      advanced to next number
*     DRESLT     d      result
*     JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
*
*  Notes:
*
*     1     The reason DFLTIN has separate OK status values for +
*           and - is to enable minus zero to be detected.   This is
*           of crucial importance when decoding mixed-radix numbers.
*           For example, an angle expressed as deg, arcmin, arcsec
*           may have a leading minus sign but a zero degrees field.
*
*     2     A TAB is interpreted as a space, and lowercase characters
*           are interpreted as uppercase.
*
*     3     The basic format is the sequence of fields #^.^@#^, where
*           # is a sign character + or -, ^ means a string of decimal
*           digits, and @, which indicates an exponent, means D or E.
*           Various combinations of these fields can be omitted, and
*           embedded blanks are permissible in certain places.
*
*     4     Spaces:
*
*             .  Leading spaces are ignored.
*
*             .  Embedded spaces are allowed only after +, -, D or E,
*                and after the decomal point if the first sequence of
*                digits is absent.
*
*             .  Trailing spaces are ignored;  the first signifies
*                end of decoding and subsequent ones are skipped.
*
*     5     Delimiters:
*
*             .  Any character other than +,-,0-9,.,D,E or space may be
*                used to signal the end of the number and terminate
*                decoding.
*
*             .  Comma is recognized by DFLTIN as a special case;  it
*                is skipped, leaving the pointer on the next character.
*                See 13, below.
*
*     6     Both signs are optional.  The default is +.
*
*     7     The mantissa ^.^ defaults to 1.
*
*     8     The exponent @#^ defaults to D0.
*
*     9     The strings of decimal digits may be of any length.
*
*     10    The decimal point is optional for whole numbers.
*
*     11    A "null result" occurs when the string of characters being
*           decoded does not begin with +,-,0-9,.,D or E, or consists
*           entirely of spaces.  When this condition is detected, JFLAG
*           is set to 1 and DRESLT is left untouched.
*
*     12    NSTRT = 1 for the first character in the string.
*
*     13    On return from DFLTIN, NSTRT is set ready for the next
*           decode - following trailing blanks and any comma.  If a
*           delimiter other than comma is being used, NSTRT must be
*           incremented before the next call to DFLTIN, otherwise
*           all subsequent calls will return a null result.
*
*     14    Errors (JFLAG=2) occur when:
*
*             .  a +, -, D or E is left unsatisfied;  or
*
*             .  the decimal point is present without at least
*                one decimal digit before or after it;  or
*
*             .  an exponent more than 100 has been presented.
*
*     15    When an error has been detected, NSTRT is left
*           pointing to the character following the last
*           one used before the error came to light.  This
*           may be after the point at which a more sophisticated
*           program could have detected the error.  For example,
*           DFLTIN does not detect that '1D999' is unacceptable
*           (on a computer where this is so) until the entire number
*           has been decoded.
*
*     16    Certain highly unlikely combinations of mantissa &
*           exponent can cause arithmetic faults during the
*           decode, in some cases despite the fact that they
*           together could be construed as a valid number.
*
*     17    Decoding is left to right, one pass.
*
*     18    See also FLOTIN and INTIN
*
*  Called:  sla__IDCHF
*
*  P.T.Wallace   Starlink   18 March 1999
*
*  Copyright (C) 1999 Rutherford Appleton Laboratory
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
      INTEGER NSTRT
      DOUBLE PRECISION DRESLT
      INTEGER JFLAG

      INTEGER NPTR,MSIGN,NEXP,NDP,NVEC,NDIGIT,ISIGNX,J
      DOUBLE PRECISION DMANT,DIGIT



*  Current character
      NPTR=NSTRT

*  Set defaults: mantissa & sign, exponent & sign, decimal place count
      DMANT=0D0
      MSIGN=1
      NEXP=0
      ISIGNX=1
      NDP=0

*  Look for sign
 100  CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO ( 400,  100,  800,  500,  300,  200, 9110, 9100, 9110),NVEC
*             0-9    SP   D/E    .     +     -     ,   ELSE   END

*  Negative
 200  CONTINUE
      MSIGN=-1

*  Look for first leading decimal
 300  CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO ( 400, 300,  800,  500, 9200, 9200, 9200, 9200, 9210),NVEC
*             0-9   SP   D/E    .     +     -     ,   ELSE   END

*  Accept leading decimals
 400  CONTINUE
      DMANT=DMANT*1D1+DIGIT
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO ( 400, 1310,  900,  600, 1300, 1300, 1300, 1300, 1310),NVEC
*             0-9   SP    D/E    .     +     -     ,   ELSE   END

*  Look for decimal when none preceded the point
 500  CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO ( 700, 500, 9200, 9200, 9200, 9200, 9200, 9200, 9210),NVEC
*             0-9   SP   D/E    .     +     -     ,   ELSE   END

*  Look for trailing decimals
 600  CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO ( 700, 1310,  900, 1300, 1300, 1300, 1300, 1300, 1310),NVEC
*             0-9   SP    D/E    .     +     -     ,   ELSE   END

*  Accept trailing decimals
 700  CONTINUE
      NDP=NDP+1
      DMANT=DMANT*1D1+DIGIT
      GO TO 600

*  Exponent symbol first in field: default mantissa to 1
 800  CONTINUE
      DMANT=1D0

*  Look for sign of exponent
 900  CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO (1200, 900, 9200, 9200, 1100, 1000, 9200, 9200, 9210),NVEC
*             0-9   SP   D/E    .     +     -     ,   ELSE   END

*  Exponent negative
 1000 CONTINUE
      ISIGNX=-1

*  Look for first digit of exponent
 1100 CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO (1200, 1100, 9200, 9200, 9200, 9200, 9200, 9200, 9210),NVEC
*             0-9   SP    D/E    .     +     -     ,   ELSE   END

*  Use exponent digit
 1200 CONTINUE
      NEXP=NEXP*10+NDIGIT
      IF (NEXP.GT.100) GO TO 9200

*  Look for subsequent digits of exponent
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO (1200, 1310, 1300, 1300, 1300, 1300, 1300, 1300, 1310),NVEC
*             0-9   SP    D/E    .     +     -     ,   ELSE   END

*  Combine exponent and decimal place count
 1300 CONTINUE
      NPTR=NPTR-1
 1310 CONTINUE
      NEXP=NEXP*ISIGNX-NDP

*  Skip if net exponent negative
      IF (NEXP.LT.0) GO TO 1500

*  Positive exponent: scale up
 1400 CONTINUE
      IF (NEXP.LT.10) GO TO 1410
      DMANT=DMANT*1D10
      NEXP=NEXP-10
      GO TO 1400
 1410 CONTINUE
      IF (NEXP.LT.1) GO TO 1600
      DMANT=DMANT*1D1
      NEXP=NEXP-1
      GO TO 1410

*  Negative exponent: scale down
 1500 CONTINUE
      IF (NEXP.GT.-10) GO TO 1510
      DMANT=DMANT/1D10
      NEXP=NEXP+10
      GO TO 1500
 1510 CONTINUE
      IF (NEXP.GT.-1) GO TO 1600
      DMANT=DMANT/1D1
      NEXP=NEXP+1
      GO TO 1510

*  Get result & status
 1600 CONTINUE
      J=0
      IF (MSIGN.EQ.1) GO TO 1610
      J=-1
      DMANT=-DMANT
 1610 CONTINUE
      DRESLT=DMANT

*  Skip to end of field
 1620 CONTINUE
      CALL sla__IDCHF(STRING,NPTR,NVEC,NDIGIT,DIGIT)
      GO TO (1720, 1620, 1720, 1720, 1720, 1720, 9900, 1720, 9900),NVEC
*             0-9   SP    D/E    .     +     -     ,   ELSE   END

 1720 CONTINUE
      NPTR=NPTR-1
      GO TO 9900


*  Exits

*  Null field
 9100 CONTINUE
      NPTR=NPTR-1
 9110 CONTINUE
      J=1
      GO TO 9900

*  Errors
 9200 CONTINUE
      NPTR=NPTR-1
 9210 CONTINUE
      J=2

*  Return
 9900 CONTINUE
      NSTRT=NPTR
      JFLAG=J

      END
