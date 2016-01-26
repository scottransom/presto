/*
 *     - - - - - 
 *      C L D J 
 *     - - - - - 
 *
 *  Gregorian Calendar to Modified Julian Date
 *
 *  Given: 
 *     IY,IM,ID     int    year, month, day in Gregorian calendar 
 *
 *  Returned: 
 *     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs 
 *     J            int    status: 
 *                           0 = OK 
 *                           1 = bad year   (MJD not computed) 
 *                           2 = bad month  (MJD not computed) 
 *                           3 = bad day    (MJD computed) 
 *
 *  The year must be -4699 (i.e. 4700BC) or later. 
 *
 *  The algorithm is derived from that of Hatcher 1984 
 *  (QJRAS 25, 53-55). 
 *
 *  P.T.Wallace   Starlink   11 March 1998 
 *
 *  Copyright (C) 1998 Rutherford Appleton Laboratory 
 *
 */

static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

double slaCldj(int iy, int im, int id, int *j)
{
    double mjd = 0.0;
    /*  Month lengths in days */

    /* Preset status */
    *j = 0;
    /* Validate year */
    if (iy < -4699) {
        *j = 1;
    } else {
        /* Validate month */
        if (im >= 1 && im <= 12) {
            /* Allow for leap year */
            if (iy % 4 == 0) {
                mtab[1] = 29;
            } else {
                mtab[1] = 28;
            }
            if (iy % 100 == 0 && iy % 400 != 0) {
                mtab[1] = 28;
            }
            /* Validate day */
            if (id < 1 || id > mtab[im - 1]) {
                *j = 3;
            }
            /* Modified Julian Date */
            mjd = (double) ((iy - (12 - im) / 10 + 4712) * 1461 / 4 +
                            ((im + 9) % 12 * 306 + 5) / 10 -
                            (iy - (12 - im) / 10 + 4900) / 100 * 3 / 4 + id -
                            2399904);
        } else {                /* Bad month */
            *j = 2;
        }
    }
    return mjd;
}
