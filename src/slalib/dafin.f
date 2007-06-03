      SUBROUTINE sla_DAFIN (STRING, IPTR, A, J)
*+
*     - - - - - -
*      D A F I N
*     - - - - - -
*
*  Sexagesimal character string to angle (double precision)
*
*  Given:
*     STRING  c*(*)   string containing deg, arcmin, arcsec fields
*     IPTR      i     pointer to start of decode (1st = 1)
*
*  Returned:
*     IPTR      i     advanced past the decoded angle
*     A         d     angle in radians
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
*    A           ?                                -1.00000D0
*    J           ?                                0
*
*    A further call to sla_DAFIN, without adjustment of IPTR, will
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
*           expected to be integral unless the other two fields are absent.
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
*         etc. See the example, above.
*
*     6)  If STRING contains hours, minutes, seconds instead of degrees
*         etc, or if the required units are turns (or days) instead of
*         radians, the result A should be multiplied as follows:
*
*           for        to obtain    multiply
*           STRING     A in         A by
*
*           d ' "      radians      1       =  1D0
*           d ' "      turns        1/2pi   =  0.1591549430918953358D0
*           h m s      radians      15      =  15D0
*           h m s      days         15/2pi  =  2.3873241463784300365D0
*
*  Called:  sla_DFLTIN
*
*  P.T.Wallace   Starlink   1 August 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
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
      DOUBLE PRECISION A
      INTEGER J

      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=4.84813681109535993589914102358D-6)
      INTEGER JF,JD,JM,JS
      DOUBLE PRECISION DEG,ARCMIN,ARCSEC



*  Preset the status to OK
      JF=0

*  Defaults
      DEG=0D0
      ARCMIN=0D0
      ARCSEC=0D0

*  Decode degrees, arcminutes, arcseconds
      CALL sla_DFLTIN(STRING,IPTR,DEG,JD)
      IF (JD.GT.1) THEN
         JF=-1
      ELSE
         CALL sla_DFLTIN(STRING,IPTR,ARCMIN,JM)
         IF (JM.LT.0.OR.JM.GT.1) THEN
            JF=-2
         ELSE
            CALL sla_DFLTIN(STRING,IPTR,ARCSEC,JS)
            IF (JS.LT.0.OR.JS.GT.1) THEN
               JF=-3

*        See if the combination of fields is credible
            ELSE IF (JD.GT.0) THEN
*           No degrees:  arcmin, arcsec ought also to be absent
               IF (JM.EQ.0) THEN
*              Suspect arcmin
                  JF=-2
               ELSE IF (JS.EQ.0) THEN
*              Suspect arcsec
                  JF=-3
               ELSE
*              All three fields absent
                  JF=1
               END IF
*        Degrees present:  if arcsec present so ought arcmin to be
            ELSE IF (JM.NE.0.AND.JS.EQ.0) THEN
               JF=-3

*        Tests for range and integrality

*        Degrees
            ELSE IF (JM.EQ.0.AND.DINT(DEG).NE.DEG) THEN
               JF=-1
*        Arcminutes
            ELSE IF ((JS.EQ.0.AND.DINT(ARCMIN).NE.ARCMIN).OR.
     :               ARCMIN.GE.60D0) THEN
               JF=-2
*        Arcseconds
            ELSE IF (ARCSEC.GE.60D0) THEN
               JF=-3
            END IF
         END IF
      END IF

*  Unless all three fields absent, compute angle value
      IF (JF.LE.0) THEN
         A=AS2R*(60D0*(60D0*ABS(DEG)+ARCMIN)+ARCSEC)
         IF (JD.LT.0) A=-A
      END IF

*  Return the status
      J=JF

      END
