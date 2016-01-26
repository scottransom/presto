#include <stdio.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
    fitsfile *infptr;           /* FITS file pointers defined in fitsio.h */
    int status = 0;             /* status must always be initialized = 0  */

    if (argc == 1) {
        printf("usage:  fitsdelrow fitsfile.fits[3] 1 7\n"
               "        where the optional [#] means the HDU number, and 1 and 7\n"
               "        are the start row and number of rows to delete\n\n");
        exit(0);
    }

    /* Open the input file */
    if (!fits_open_file(&infptr, argv[1], READWRITE, &status)) {
        /* Delete the first 7 rows */
        fits_delete_rows(infptr, atoi(argv[2]), atoi(argv[3]), &status);
        /* Close the input file */
        fits_close_file(infptr, &status);
    }

    /* if error occured, print out error message */
    if (status)
        fits_report_error(stderr, status);
    return (status);
}
