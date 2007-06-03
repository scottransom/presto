#include <sla_config.h>
      REAL FUNCTION sla_GRESID (S)
*+
*     - - - - - - -
*      G R E S I D
*     - - - - - - -
*
*  Generate pseudo-random normal deviate ( = 'Gaussian residual')
*  (single precision)
*
*  Given:
*     S      real     standard deviation
*
*  The results of many calls to this routine will be
*  normally distributed with mean zero and standard deviation S.
*
*  The Box-Muller algorithm is used.  This is described in
*  Numerical Recipes, section 7.2.
*
*  Called:  RAN or RAND (a REAL function returning a random variate --
*           the precise function which is called depends on which functions
*           are available when the library is built).  If neither of these
*           is available, we use the local substitute RANDOM defined
*           in rtl_random.c
*
*  P.T.Wallace   Starlink   14 October 1991
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
*-

      IMPLICIT NONE

      REAL S

      REAL X,Y,R,W,GNEXT,G
      LOGICAL FTF,FIRST

#if HAVE_RAND
      REAL RAND
#elif HAVE_RANDOM
      REAL RANDOM
#else
 error "Can't find random-number function"
#endif

      SAVE GNEXT,FTF,FIRST
      DATA FTF,FIRST / .TRUE.,.TRUE. /

      X = 0.0
      Y = 0.0

*  First time through, initialise the random-number generator
#if HAVE_RAND
      IF (FTF) THEN
         X = RAND(123456789)
         FTF = .FALSE.
      END IF
#endif

*  Second normal deviate of the pair available?
      IF (FIRST) THEN

*     No - generate two random numbers inside unit circle
         R = 2.0
         DO WHILE (R.GE.1.0)

*        Generate two random numbers in range +/- 1
#if HAVE_RAND
            X = 2.0*RAND(0)-1.0
            Y = 2.0*RAND(0)-1.0
#elif HAVE_RANDOM
            X = 2.0*RAN(ISEED)-1.0
            Y = 2.0*RAN(ISEED)-1.0
#endif

*        Try again if not in unit circle
            R = X*X+Y*Y
         END DO

*     Box-Muller transformation, generating two deviates
         W = SQRT(-2.0*LOG(R)/MAX(R,1E-20))
         GNEXT = X*W
         G = Y*W

*     Set flag to indicate availability of next deviate
         FIRST = .FALSE.
      ELSE

*     Return second deviate of the pair & reset flag
         G = GNEXT
         FIRST = .TRUE.
      END IF

*  Scale the deviate by the required standard deviation
      sla_GRESID = G*S

      END
