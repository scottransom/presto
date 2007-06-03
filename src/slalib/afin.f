      SUBROUTINE sla_AFIN (STRING, IPTR, A, J)
*+
*     - - - - -
*      A F I N
*     - - - - -
*
*  Sexagesimal character string to angle (single precision)
*
*  Given:
*     STRING  c*(*)   string containing deg, arcmin, arcsec fields
*     IPTR      i     pointer to start of decode (1st = 1)
*
*  Returned:
*     IPTR      i     advanced past the decoded angle
*     A         r     angle in radians
*     J         i     status:  0 = OK
*                             +1 = default, A unchanged
*                             -1 = bad degrees      )
*                             -2 = bad arcminutes   )  (note 3)
*                             -3 = bad arcseconds   )
*
*  Example:
*
*    argument    before                           after
*
*    STRING      '-57 17 44.806  12 34 56.7'      unchanged
*    IPTR        1                                16 (points to 12...)
*    A           ?                                -1.00000
*    J           ?                                0
*
*    A further call to sla_AFIN, without adjustment of IPTR, will
*    decode the second angle, 12deg 34min 56.7sec.
*
*  Notes:
*
*     1)  The first three "fields" in STRING are degrees, arcminutes,
*         arcseconds, separated by spaces or commas.  The degrees field
*         may be signed, but not the others.  The decoding is carried
*         out by the DFLTIN routine and is free-format.
*
*     2)  Successive fields may be absent, defaulting to zero.  For
*         zero status, the only combinations allowed are degrees alone,
*         degrees and arcminutes, and all three fields present.  If all
*         three fields are omitted, a status of +1 is returned and A is
*         unchanged.  In all other cases A is changed.
*
*     3)  Range checking:
*
*           The degrees field is not range checked.  However, it is
*           expected to be integral unless the other two fields are
*           absent.
*
*           The arcminutes field is expected to be 0-59, and integral if
*           the arcseconds field is present.  If the arcseconds field
*           is absent, the arcminutes is expected to be 0-59.9999...
*
*           The arcseconds field is expected to be 0-59.9999...
*
*     4)  Decoding continues even when a check has failed.  Under these
*         circumstances the field takes the supplied value, defaulting
*         to zero, and the result A is computed and returned.
*
*     5)  Further fields after the three expected ones are not treated
*         as an error.  The pointer IPTR is left in the correct state
*         for further decoding with the present routine or with DFLTIN
*         etc.  See the example, above.
*
*     6)  If STRING contains hours, minutes, seconds instead of degrees
*         etc, or if the required units are turns (or days) instead of
*         radians, the result A should be multiplied as follows:
*
*           for        to obtain    multiply
*           STRING     A in         A by
*
*           d ' "      radians      1       =  1.0
*           d ' "      turns        1/2pi   =  0.1591549430918953358
*           h m s      radians      15      =  15.0
*           h m s      days         15/2pi  =  2.3873241463784300365
*
*  Called:  sla_DAFIN
*
*  P.T.Wallace   Starlink   13 September 1990
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
      INTEGER IPTR
      REAL A
      INTEGER J

      DOUBLE PRECISION AD



*  Call the double precision version
      CALL sla_DAFIN(STRING,IPTR,AD,J)
      IF (J.LE.0) A=REAL(AD)

      END
