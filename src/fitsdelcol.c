#include <stdio.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
    fitsfile *infptr;           /* FITS file pointers defined in fitsio.h */
    int status = 0;             /* status must always be initialized = 0  */
    int errors = 0;             /* 0 means no errors were encountered     */
    int colnum;                 /* column number                          */
    int ii;                     /* dummy index                            */

    if (argc < 3) {
        printf("usage:  fitsdelcol fitsfile.fits[HDU] COLNUM1 [COLNUM2 ...]\n"
               "        where HDU is the HDU number or name to use, \n"
               "        and COLNUM is the name or index of the column to remove\n\n");
        exit(0);
    }

    /* Open the input file */
    if (!fits_open_file(&infptr, argv[1], READWRITE, &status)) {
        for (ii = 2; ii < argc; ii++) {
            fits_get_colnum(infptr, CASEINSEN, argv[ii], &colnum, &status);
            if (status) {
                /* if error occured, print out error message */
                fits_report_error(stderr, status);
                status = 0;     /* Reset to properly register new errors */
                errors = 1;
            } else {
                fits_delete_col(infptr, colnum, &status);
                /* if error occured, print out error message */
                if (status) {
                    fits_report_error(stderr, status);
                    status = 0; /* Reset to properly register new errors */
                    errors = 1;
                }
            }
        }
        /* Close the input file */
        fits_close_file(infptr, &status);
        if (status) {
            fits_report_error(stderr, status);
            errors = 1;
        }
    }

    return (errors);
}
