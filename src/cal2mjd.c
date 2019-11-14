#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double slaCldj(int iy, int im, int id, int *j);

/* Does not include jump discontinuities from leap seconds.  */
/* I don't really know if it should !                        */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    int year, month = 1, day = 1, hour = 0, min = 0, err;
    double MJD, fracday, sec = 0.0;

    if (argc < 2) {
        printf("\nUsage:  'cal2mjd YYYY MM DD HH MM SS.SSSSS'\n\n");
        exit(1);
    }
    year = strtol(argv[1], NULL, 10);
    if (argc > 2) {
        month = strtol(argv[2], NULL, 10);
        if (month < 1 || month > 12) {
            printf("\nmonth = %d is out-of-range.\n", month);
            exit(1);
        }
    }
    if (argc > 3) {
        day = strtol(argv[3], NULL, 10);
        if (day < 1 || day > 31) {
            printf("\nday = %d is out-of-range.\n", month);
            exit(1);
        }
    }
    if (argc > 4) {
        hour = strtol(argv[4], NULL, 10);
        if (hour == 24)
            hour = 0;
        if (hour < 0 || hour > 23) {
            printf("\nhour = %d is out-of-range.\n", month);
            exit(1);
        }
    }
    if (argc > 5) {
        min = strtol(argv[5], NULL, 10);
        if (min < 0 || min > 59) {
            printf("\nmin = %d is out-of-range.\n", month);
            exit(1);
        }
    }
    if (argc > 6) {
        sec = strtod(argv[6], NULL);
        if (sec < 0.0 || sec >= 60.0) {
            printf("\nsec = %d is out-of-range.\n", month);
            exit(1);
        }
    }
    fracday = (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
    MJD = slaCldj(year, month, day, &err);
    MJD += fracday;
    if (err == 1) {
        printf("\nTry again.  Bad year.\n\n");
        exit(1);
    }
    if (err == 2) {
        printf("\nTry again.  Bad month.\n\n");
        exit(1);
    }
    if (err == 3) {
        printf("\nTry again.  Bad day.\n\n");
        exit(1);
    }
    printf("\nMJD is %17.11f\n\n", MJD);
    exit(0);
}
