      SUBROUTINE sla_FLOTIN (STRING, NSTRT, RESLT, JFLAG)
*+
*     - - - - - - -
*      F L O T I N
*     - - - - - - -
*
*  Convert free-format input into single precision floating point
*
*  Given:
*     STRING     c     string containing number to be decoded
*     NSTRT      i     pointer to where decoding is to start
*     RESLT      r     current value of result
*
*  Returned:
*     NSTRT      i      advanced to next number
*     RESLT      r      result
*     JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
*
*  Called:  sla_DFLTIN
*
*  Notes:
*
*     1     The reason FLOTIN has separate OK status values for +
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
*             .  Comma is recognized by FLOTIN as a special case;  it
*                is skipped, leaving the pointer on the next character.
*                See 13, below.
*
*     6     Both signs are optional.  The default is +.
*
*     7     The mantissa ^.^ defaults to 1.
*
*     8     The exponent @#^ defaults to E0.
*
*     9     The strings of decimal digits may be of any length.
*
*     10    The decimal point is optional for whole numbers.
*
*     11    A "null result" occurs when the string of characters being
*           decoded does not begin with +,-,0-9,.,D or E, or consists
*           entirely of spaces.  When this condition is detected, JFLAG
*           is set to 1 and RESLT is left untouched.
*
*     12    NSTRT = 1 for the first character in the string.
*
*     13    On return from FLOTIN, NSTRT is set ready for the next
*           decode - following trailing blanks and any comma.  If a
*           delimiter other than comma is being used, NSTRT must be
*           incremented before the next call to FLOTIN, otherwise
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
*           FLOTIN does not detect that '1E999' is unacceptable
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
*     18    See also DFLTIN and INTIN
*
*  P.T.Wallace   Starlink   23 November 1995
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

      CHARACTER*(*) STRING
      INTEGER NSTRT
      REAL RESLT
      INTEGER JFLAG

      DOUBLE PRECISION DRESLT


*  Call the double precision version
      CALL sla_DFLTIN(STRING,NSTRT,DRESLT,JFLAG)
      IF (JFLAG.LE.0) RESLT=REAL(DRESLT)

      END
