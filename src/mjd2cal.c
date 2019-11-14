#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "slamac.h"

char months[12][4] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

void slaDjcl(double djm, int *iy, int *im, int *id, double *fd, int *j);

/* Does not include jump discontinuities from leap seconds.  */
/* I don't really know if it should !                        */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    int year, month, day, hour, min, err;
    double MJD, fracday, dh, dm, sec;

    if (argc != 2) {
        printf("\nUsage:  'mjd2cal MJD'\n\n");
        exit(1);
    }
    MJD = atof(argv[1]);
    slaDjcl(MJD, &year, &month, &day, &fracday, &err);

    if (err == -1) {
        printf("\nTry again.  Bad MJD.\n\n");
        exit(1);
    }
    dh = fracday * 24.0;
    hour = (int) dh;
    dm = (dh - hour) * 60.0;
    min = (int) dm;
    sec = (dm - min) * 60.0;

    printf("\nDate is %2d %s %4d at %2d hours %2d minutes and %.8g seconds\n\n",
           day, months[month - 1], year, hour, min, sec);
    exit(0);
}
