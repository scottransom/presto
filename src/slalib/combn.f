      SUBROUTINE sla_COMBN ( NSEL, NCAND, LIST, J )
*+
*     - - - - - -
*      C O M B N
*     - - - - - -
*
*  Generate the next combination, a subset of a specified size chosen
*  from a specified number of items.
*
*  Given:
*     NSEL     i        number of items (subset size)
*     NCAND    i        number of candidates (set size)
*
*  Given and returned:
*     LIST     i(NSEL)  latest combination, LIST(1)=0 to initialize
*
*  Returned:
*     J        i        status: -1 = illegal NSEL or NCAND
*                                0 = OK
*                               +1 = no more combinations available
*
*  Notes:
*
*  1) NSEL and NCAND must both be at least 1, and NSEL must be less
*     than or equal to NCAND.
*
*  2) This routine returns, in the LIST array, a subset of NSEL integers
*     chosen from the range 1 to NCAND inclusive, in ascending order.
*     Before calling the routine for the first time, the caller must set
*     the first element of the LIST array to zero (any value less than 1
*     will do) to cause initialization.
*
*  2) The first combination to be generated is:
*
*        LIST(1)=1, LIST(2)=2, ..., LIST(NSEL)=NSEL
*
*     This is also the combination returned for the "finished" (J=1)
*     case.
*
*     The final permutation to be generated is:
*
*        LIST(1)=NCAND, LIST(2)=NCAND-1, ..., LIST(NSEL)=NCAND-NSEL+1
*
*  3) If the "finished" (J=1) status is ignored, the routine
*     continues to deliver combinations, the pattern repeating
*     every NCAND!/(NSEL!*(NCAND-NSEL)!) calls.
*
*  4) The algorithm is by R.F.Warren-Smith (private communication).
*
*  P.T.Wallace   Starlink   25 August 1999
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

      INTEGER NSEL,NCAND,LIST(NSEL),J

      INTEGER I,LISTI,NMAX,M
      LOGICAL MORE


*  Validate, and set status.
      IF (NSEL.LT.1.OR.NCAND.LT.1.OR.NSEL.GT.NCAND) THEN
         J = -1
         GO TO 9999
      ELSE
         J = 0
      END IF

*  Just starting?
      IF (LIST(1).LT.1) THEN

*     Yes: return 1,2,3...
         DO I=1,NSEL
            LIST(I) = I
         END DO

      ELSE

*     No: find the first selection that we can increment.

*     Start with the first list item.
         I = 1

*     Loop.
         MORE = .TRUE.
         DO WHILE (MORE)

*        Current list item.
            LISTI = LIST(I)

*        Is this the final list item?
            IF (I.GE.NSEL) THEN

*           Yes:  comparison value is number of candidates plus one.
               NMAX = NCAND+1
            ELSE

*           No:  comparison value is next list item.
               NMAX = LIST(I+1)
            END IF

*        Can the current item be incremented?
            IF (NMAX-LISTI.GT.1) THEN

*           Yes:  increment it.
               LIST(I) = LISTI+1

*           Reinitialize the preceding items.
               DO M=1,I-1
                  LIST(M) = M
               END DO

*           Break.
               MORE = .FALSE.
            ELSE

*           Can't increment the current item:  is it the final one?
               IF (I.GE.NSEL) THEN

*              Yes:  set the status.
                  J = 1

*              Restart the sequence.
                  DO I=1,NSEL
                     LIST(I) = I
                  END DO

*              Break.
                  MORE = .FALSE.
               ELSE

*              No:  next list item.
                  I = I+1
               END IF
            END IF
         END DO
      END IF
 9999 CONTINUE

      END
