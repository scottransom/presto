#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fitsio.h"

int get_colnum(fitsfile * infile, char *colname, int *status)
{
    int colnum;

    fits_get_colnum(infile, 0, colname, &colnum, status);
    if (*status) {
        fprintf(stderr, "  Error in 'get_colnum' looking for '%s'!\n", colname);
        fits_report_error(stderr, *status);
        exit(1);
    }
    return colnum;
}


void get_float_column(fitsfile * infile, long long rownum, char *colname,
                      float *data, int N, int *status)
{
    int colnum, anynull;

    colnum = get_colnum(infile, colname, status);
    fits_read_col(infile, TFLOAT, colnum, rownum, 1L, N, 0, data, &anynull, status);
    if (*status) {
        fprintf(stderr, "  Error in 'get_float_column' reading '%s'!\n", colname);
        fits_report_error(stderr, *status);
        exit(1);
    }
}


int main(int argc, char *argv[])
{
    fitsfile *infile;
    int ii, jj, status = 0;
    int nchan, npol;
    long long nrows, rownum = 1;
    float *freqs, *weights, *offsets, *scales;
    char comment[120];

    if (argc == 1) {
        printf("usage:  psrfits_dumparrays PSRFITS_filename [row=1]\n");
        exit(0);
    }
    if (argc > 2) {
        rownum = atoll(argv[2]);        // defaults to 1
    }
    // Open the PSRFITS file
    fits_open_file(&infile, argv[1], READONLY, &status);
    if (status) {
        fprintf(stderr, "  Error!  Cannot open '%s'!\n", argv[1]);
        fits_report_error(stderr, status);
        exit(1);
    }
    // Is the file a PSRFITS file?
    {
        char ctmp[80];
        fits_read_key(infile, TSTRING, "FITSTYPE", ctmp, comment, &status);
        if (status || strcmp(ctmp, "PSRFITS")) {
            fprintf(stderr,
                    "  Error!  '%s' does not appear to be PSRFITS!\n", argv[1]);
            if (status)
                fits_report_error(stderr, status);
            exit(1);
        }
    }

    // Move to the SUBINT HDU
    fits_movnam_hdu(infile, BINARY_TBL, "SUBINT", 0, &status);
    if (status) {
        fits_report_error(stderr, status);
        exit(1);
    }
    // Read the number of channels and polarizations
    fits_read_key(infile, TINT, "NCHAN", &nchan, comment, &status);
    if (status) {
        fprintf(stderr, "  Error!  Cannot find NCHAN in '%s'!\n", argv[1]);
        fits_report_error(stderr, status);
        exit(1);
    }
    fits_read_key(infile, TINT, "NPOL", &npol, comment, &status);
    if (status) {
        fprintf(stderr, "  Error!  Cannot find NPOL in '%s'!\n", argv[1]);
        fits_report_error(stderr, status);
        exit(1);
    }
    // How many rows are there?
    fits_get_num_rowsll(infile, &nrows, &status);
    if (status) {
        fprintf(stderr,
                "  Error!  Cannot read the number of rows in '%s'!\n", argv[1]);
        fits_report_error(stderr, status);
        exit(1);
    }
    if (rownum > nrows) {
        fprintf(stderr,
                "  Error!  Requested row %lld is greater than the number of rows %lld!\n",
                rownum, nrows);
        exit(1);
    }
    // Allocate the arrays
    freqs = (float *) malloc(sizeof(float) * nchan);
    weights = (float *) malloc(sizeof(float) * nchan);
    scales = (float *) malloc(sizeof(float) * nchan * npol);
    offsets = (float *) malloc(sizeof(float) * nchan * npol);

    // Read the columns
    get_float_column(infile, rownum, "DAT_FREQ", freqs, nchan, &status);
    get_float_column(infile, rownum, "DAT_WTS", weights, nchan, &status);
    get_float_column(infile, rownum, "DAT_SCL", scales, nchan * npol, &status);
    get_float_column(infile, rownum, "DAT_OFFS", offsets, nchan * npol, &status);

    // Now print them to STDOUT
    printf
        ("#  Chan            Freq(MHz)         Weights         Scales          Offsets\n"
         "#-----------------------------------------------------------------------------\n");
    for (ii = 0; ii < npol; ii++) {
        for (jj = 0; jj < nchan; jj++) {
            if (ii > 0) {
                printf("%6d (pol%d)  %15.7f %15.7f %15.7f %15.7f\n",
                       jj, ii + 1, freqs[jj], weights[jj],
                       scales[ii * nchan + jj], offsets[ii * nchan + jj]);
            } else {
                printf("%6d         %15.7f %15.7f %15.7f %15.7f\n",
                       jj, freqs[jj], weights[jj], scales[jj], offsets[jj]);
            }
        }
    }

    // Close the file and cleanup
    fits_close_file(infile, &status);
    if (status) {
        fprintf(stderr, "  Error!  Cannot properly close '%s'!\n", argv[1]);
        fits_report_error(stderr, status);
        exit(1);
    }
    free(freqs);
    free(weights);
    free(scales);
    free(offsets);
    exit(0);
}
