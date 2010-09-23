#include <stdio.h>
#include "fitsio.h"
// v2
int main(int argc, char *argv[])
{
    fitsfile *infptr;   /* FITS file pointers defined in fitsio.h */
    int status = 0, ii = 1;       /* status must always be initialized = 0  */

    /* Open the input file */
    if ( !fits_open_file(&infptr, argv[1], READWRITE, &status) )
    {
        /* Delete the first 7 rows */
      fits_delete_rows(infptr, atoi(argv[2]), atoi(argv[3]), &status);
	/* Close the input file */
        fits_close_file(infptr, &status);
    }

    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);
    return(status);
}
