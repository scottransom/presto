#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slamac.h"

void slaDjcl(double djm, int *iy, int *im, int *id, double *fd, int *j);

void slaDjcl(double djm, int *iy, int *im, int *id, double *fd, int *j)
/*
   **  - - - - - - - -
   **   s l a D j c l
   **  - - - - - - - -
   **
   **  Modified Julian Date to Gregorian year, month, day,
   **  and fraction of a day.
   **
   **  Given:
   **     djm      double     Modified Julian Date (JD-2400000.5)
   **
   **  Returned:
   **     *iy      int        year
   **     *im      int        month
   **     *id      int        day
   **     *fd      double     fraction of day
   **     *j       int        status:
   **                      -1 = unacceptable date (before 4701BC March 1)
   **
   **  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
   **
   **  Defined in slamac.h:  dmod
   **
   **  Last revision:   20 April 1996
   **
   **  Copyright P.T.Wallace.  All rights reserved.
 */
{
    double f, d;
    long jd, n4, nd10;

/* Check if date is acceptable */
    if ((djm <= -2395522.0) || (djm >= 1e9)) {
        *j = -1;
        return;
    } else {
        *j = 0;

        /* Separate day and fraction */
        f = dmod(djm, 1.0);
        if (f < 0.0)
            f += 1.0;
        d = djm - f;
        d = dnint(d);

        /* Express day in Gregorian calendar */
        jd = (long) dnint(d) + 2400001;
        n4 = 4L * (jd + ((6L * ((4L * jd - 17918L) / 146097L)) / 4L + 1L) / 2L -
                   37L);
        nd10 = 10L * (((n4 - 237L) % 1461L) / 4L) + 5L;
        *iy = (int) (n4 / 1461L - 4712L);
        *im = (int) (((nd10 / 306L + 2L) % 12L) + 1L);
        *id = (int) ((nd10 % 306L) / 10L + 1L);
        *fd = f;
        *j = 0;
    }
}
