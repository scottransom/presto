#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "presto.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define BUFFLEN 65536

int main(int argc, char *argv[])
{
    FILE *infile, *outfile;
    int ii, offset, numread;
    long long N = 0;
    float min = 1e99, max = -1e99, fbuffer[sizeof(float) * BUFFLEN];
    double avg = 0.0, var = 0.0;
    short sbuffer[sizeof(short) * BUFFLEN];
    char *outname, *statsname;

    if (argc != 2) {
        printf("\nUsage:  dat2sdat datfilename\n\n");
        exit(1);
    }
    printf("\n   Floats to Shorts Data Conversion Program\n\n");

    /* Get the root of the input filename and */
    /* generate the output filenames from it. */
    {
        char *rootnm, *suffix;

        split_root_suffix(argv[1], &rootnm, &suffix);
        outname = (char *) calloc(strlen(rootnm) + 6, sizeof(char));
        statsname = (char *) calloc(strlen(rootnm) + 7, sizeof(char));
        sprintf(outname, "%s.sdat", rootnm);
        sprintf(statsname, "%s.stats", rootnm);
        free(rootnm);
        free(suffix);
    }

    /* Determine time series statistics using one-pass technique */
    printf(" Determining simple statistics of '%s'...", argv[1]);
    fflush(NULL);
    infile = chkfopen(argv[1], "rb");
    while ((numread = chkfread(fbuffer, sizeof(float), BUFFLEN, infile))) {
        float xx;
        double dev;
        for (ii = 0; ii < numread; ii++) {
            xx = fbuffer[ii];
            /* Check the max and min values */
            if (xx > max)
                max = xx;
            if (xx < min)
                min = xx;
            /* Use clever single pass mean and variance calculation */
            dev = xx - avg;
            avg += dev / (N + ii + 1.0);
            var += dev * (xx - avg);
        }
        N += numread;
    }
    var /= (N - 1.0);
    offset = (int) (floor(avg));

    printf("done.\n\n");
    printf("    Number:  %lld\n", N);
    printf(" Max Value:  %.2f\n", max);
    printf(" Min Value:  %.2f\n", min);
    printf("   Average:  %.2f\n", avg);
    printf("   Std Dev:  %.2f\n", sqrt(var));
    printf("  Variance:  %.2f\n\n", var);

    if ((max - min) > (SHRT_MAX - SHRT_MIN)) {
        if ((max - min) < 1.5 * (SHRT_MAX - SHRT_MIN)) {
            printf("Warning:  There is more dynamic range in the data\n"
                   "          than can be handled perfectly:\n"
                   "               max - min = %.2f - %.2f = %.2f\n"
                   "          Clipping the low values...\n\n", max, min, max - min);
            offset = max - SHRT_MAX;
        } else {
            printf("Error:  There is way too much dynamic range in the data:\n"
                   "               max - min = %.2f - %.2f = %.2f\n"
                   "        Exiting.\n\n", max, min, max - min);
            exit(1);
        }
    }

    printf(" Value to be subtracted from each point is:  %d\n\n", offset);

    outfile = chkfopen(statsname, "w");
    fprintf(outfile, "Simple stats for '%s':\n", argv[1]);
    fprintf(outfile, "   Number of points:  %lld\n", N);
    fprintf(outfile, "      Maximum Value:  %.2f\n", max);
    fprintf(outfile, "      Minimum Value:  %.2f\n", min);
    fprintf(outfile, "      Average Value:  %.2f\n", avg);
    fprintf(outfile, " Standard Deviation:  %.2f\n", sqrt(var));
    fprintf(outfile, "           Variance:  %.2f\n", var);
    fprintf(outfile, "     Offset applied:  %d\n", -offset);
    fclose(outfile);

    printf(" Writing the new file...");
    fflush(NULL);

    rewind(infile);
    outfile = chkfopen(outname, "wb");
    while ((numread = chkfread(fbuffer, sizeof(float), BUFFLEN, infile))) {
        for (ii = 0; ii < numread; ii++)
            sbuffer[ii] = (short) (fbuffer[ii] + 1e-20 - offset);
        fwrite(sbuffer, sizeof(short), numread, outfile);
    }
    printf("done.\n\n");
    fclose(infile);
    fclose(outfile);
    free(outname);
    free(statsname);
    exit(0);
}
