#include <sla_config.h>
      REAL FUNCTION sla_RANDOM (SEED)
*+
*     - - - - - - -
*      R A N D O M
*     - - - - - - -
*
*  Generate pseudo-random real number in the range 0 <= X < 1.
*  (single precision)
*
*
*  Given:
*     SEED     real     an arbitrary real number
*
*  Notes:
*
*  1)  The result is a pseudo-random REAL number in the range
*      0 <= sla_RANDOM < 1.
*
*  2)  SEED is used first time through only.
*
*  Called:  RAN or RAND (a REAL function returning a random variate --
*           the precise function which is called depends on which functions
*           are available when the library is built).  If neither of these
*           is available, we use the local substitute RANDOM defined
*           in rtl_random.c
*
*  P.T.Wallace   Starlink   14 October 1991
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
*-

      IMPLICIT NONE

      REAL SEED

#if HAVE_RAND
      REAL RAND
#elif HAVE_RANDOM
      REAL RANDOM
#else
 error "Can't find random-number function"
#endif

      REAL AS
      INTEGER ISEED
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./



*  If first time, turn SEED into a large, odd integer
      IF (FIRST) THEN
         AS=ABS(SEED)+1.0
         ISEED=NINT(AS/10.0**(NINT(ALOG10(AS))-6))
         IF (MOD(ISEED,2).EQ.0) ISEED=ISEED+1
         FIRST=.FALSE.
#if HAVE_RAND
         AS = RAND(ISEED)
#endif
      ELSE
         ISEED=0
      END IF

*  Next pseudo-random number
#if HAVE_RAND
      sla_RANDOM=RAND(0)
#elif HAVE_RANDOM
      sla_RANDOM=RANDOM(ISEED)
#endif

      END
