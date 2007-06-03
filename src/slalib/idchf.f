      SUBROUTINE sla__IDCHF (STRING, NPTR, NVEC, NDIGIT, DIGIT)
*+
*     - - - - - -
*      I D C H F
*     - - - - - -
*
*  Internal routine used by DFLTIN
*
*  Identify next character in string
*
*  Given:
*     STRING      char        string
*     NPTR        int         pointer to character to be identified
*
*  Returned:
*     NPTR        int         incremented unless end of field
*     NVEC        int         vector for identified character
*     NDIGIT      int         0-9 if character was a numeral
*     DIGIT       double      equivalent of NDIGIT
*
*     NVEC takes the following values:
*
*      1     0-9
*      2     space or TAB   !!! n.b. ASCII TAB assumed !!!
*      3     D,d,E or e
*      4     .
*      5     +
*      6     -
*      7     ,
*      8     else
*      9     outside field
*
*  If the character is not 0-9, NDIGIT and DIGIT are either not
*  altered or are set to arbitrary values.
*
*  P.T.Wallace   Starlink   22 December 1992
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
      INTEGER NPTR,NVEC,NDIGIT
      DOUBLE PRECISION DIGIT

      CHARACTER K
      INTEGER NCHAR

*  Character/vector tables
      INTEGER NCREC
      PARAMETER (NCREC=19)
      CHARACTER KCTAB(NCREC)
      INTEGER KVTAB(NCREC)
      DATA KCTAB/'0','1','2','3','4','5','6','7','8','9',
     :           ' ','D','d','E','e','.','+','-',','/
      DATA KVTAB/10*1,2,4*3,4,5,6,7/


*  Handle pointer outside field
      IF (NPTR.LT.1.OR.NPTR.GT.LEN(STRING)) THEN
         NVEC=9
      ELSE

*     Not end of field: identify the character
         K=STRING(NPTR:NPTR)
         DO NCHAR=1,NCREC
            IF (K.EQ.KCTAB(NCHAR)) THEN

*           Recognized
               NVEC=KVTAB(NCHAR)
               NDIGIT=NCHAR-1
               DIGIT=DBLE(NDIGIT)
               GO TO 2300
            END IF
         END DO

*     Not recognized: check for TAB    !!! n.b. ASCII assumed !!!
         IF (K.EQ.CHAR(9)) THEN

*        TAB: treat as space
            NVEC=2
         ELSE

*        Unrecognized
            NVEC=8
         END IF

*     Increment pointer
 2300    CONTINUE
         NPTR=NPTR+1
      END IF

      END
