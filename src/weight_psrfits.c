#include "presto.h"
#include "mask.h"
#include "psrfits.h"

void read_wgts_and_offs(char *filenm, int *numchan, float **weights, float **offsets)
{
    FILE *infile;
    int N, chan;
    float wgt, offs;
    char line[80];

    infile = chkfopen(filenm, "r");

    // Read the input file once to count the lines
    N = 0;
    while (!feof(infile)) {
        fgets(line, 80, infile);
        if (line[0] != '#') {
            sscanf(line, "%d %f %f\n", &chan, &wgt, &offs);
            N++;
        }
    }
    N--;
    *numchan = N;

    // Allocate the output arrays
    *weights = (float *) malloc(N * sizeof(float));
    *offsets = (float *) malloc(N * sizeof(float));

    // Rewind and read the EVENTs for real
    rewind(infile);
    N = 0;
    while (!feof(infile)) {
        fgets(line, 80, infile);
        if (line[0] != '#') {
            sscanf(line, "%d %f %f\n", &chan, *weights + N, *offsets + N);
            N++;
        }
    }
    fclose(infile);
}


int main(int argc, char *argv[])
{
    fitsfile *infile;
    int ii, jj, kk, status = 0;
    int nchan, nchan2, npol, wgts_col, offs_col;
    long nrows;
    float *weights, *offsets;
    char comment[120];

    // Read the weights and offsets
    read_wgts_and_offs(argv[1], &nchan, &weights, &offsets);
    printf("Read in %d channels of weights and offsets from\n\t'%s'\n",
           nchan, argv[1]);

    // Step through the FITS files
    for (ii = 0; ii < argc - 2; ii++) {
        printf("Updating '%s'\n", argv[ii + 2]);

        // Is the file a PSRFITS file?
        if (!is_PSRFITS(argv[ii + 2])) {
            fprintf(stderr,
                    "  Error!  '%s' does not appear to be PSRFITS!\n", argv[ii + 2]);
            exit(1);
        }
        // Open the PSRFITS file
        fits_open_file(&infile, argv[ii + 2], READWRITE, &status);
        if (status) {
            printf("  Error!  Cannot open '%s'!\n", argv[ii + 2]);
            exit(1);
        }
        // Move to the SUBINT HDU
        fits_movnam_hdu(infile, BINARY_TBL, "SUBINT", 0, &status);
        if (status) {
            printf("  Warning!  Cannot find NPOL in '%s'!  Assuming NPOL=1\n",
                   argv[ii + 2]);
            status = 0;
        }
        // Read the number of channels and polarizations
        fits_read_key(infile, TINT, "NCHAN", &nchan2, comment, &status);
        if (status) {
            printf("  Warning!  Cannot find NCHAN in '%s'!\n", argv[ii + 2]);
            status = 0;
        } else if (nchan != nchan2) {
            printf("  Error!  The number of channels in '%s'\n", argv[1]);
            printf("          and in '%s' do not match!\n", argv[ii + 2]);
            exit(1);
        }
        fits_read_key(infile, TINT, "NPOL", &npol, comment, &status);
        if (status) {
            printf("  Warning!  Cannot find NPOL in '%s'!  Assuming NPOL=1\n",
                   argv[ii + 2]);
            npol = 1;
            status = 0;
        }
        // How many rows are there?
        fits_get_num_rows(infile, &nrows, &status);
        if (status) {
            printf("  Error!  Cannot read the number of rows in '%s'!\n",
                   argv[ii + 2]);
            exit(1);
        }
        // Get the column numbers for the weights
        fits_get_colnum(infile, 0, "DAT_WTS", &wgts_col, &status);
        if (status == COL_NOT_FOUND) {
            printf("  Warning!:  Can't find the channel weights!\n");
            status = 0;
        } else {
            // update the weights, row by row
            for (jj = 1; jj < nrows + 1; jj++)
                fits_write_col(infile, TFLOAT, wgts_col, jj,
                               1L, nchan, weights, &status);
        }

        // Get the column numbers for the offsets
        if (0) {
            fits_get_colnum(infile, 0, "DAT_OFFS", &offs_col, &status);
            if (status == COL_NOT_FOUND) {
                printf("  Warning!:  Can't find the channel offsets!\n");
                status = 0;
            } else {
                // update the offsets, row by row
                for (jj = 1; jj < nrows + 1; jj++)
                    for (kk = 0; kk < npol; kk++)
                        fits_write_col(infile, TFLOAT, offs_col, jj,
                                       kk * nchan + 1L, nchan, offsets, &status);
            }
        }
        // Close the file
        fits_close_file(infile, &status);
        if (status) {
            printf("  Warning!:  Cannot properly close '%s' (status=%d)!\n",
                   argv[ii + 2], status);
            status = 0;
        }
    }
    free(weights);
    free(offsets);
    printf("Finished.\n");
    exit(0);
}
