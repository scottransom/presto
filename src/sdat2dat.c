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
    int ii, numread;
    long long N = 0;
    float offset = 0.0, fbuffer[sizeof(float) * BUFFLEN];
    short sbuffer[sizeof(short) * BUFFLEN];
    char *outname;

    if (argc != 2 && argc != 3) {
        printf("\nUsage:  sdat2dat sdatfilename [offset to add]\n\n");
        exit(1);
    }
    printf("\n   Shorts to Floats Data Conversion Program\n\n");

    /* Get the root of the input filename and */
    /* generate the output filenames from it. */
    {
        char *rootnm, *suffix;

        split_root_suffix(argv[1], &rootnm, &suffix);
        outname = (char *) calloc(strlen(rootnm) + 5, sizeof(char));
        sprintf(outname, "%s.dat", rootnm);
        free(rootnm);
        free(suffix);
    }

    if (argc == 3) {
        offset = strtod(argv[2], NULL);
        printf(" Converting, adding %g, and writing to '%s'...", offset, outname);
    } else {
        printf(" Converting the data and writing to '%s'...", outname);
    }
    fflush(NULL);

    infile = chkfopen(argv[1], "rb");
    outfile = chkfopen(outname, "wb");
    while ((numread = chkfread(sbuffer, sizeof(short), BUFFLEN, infile))) {
        N += numread;
        for (ii = 0; ii < numread; ii++)
            fbuffer[ii] = (float) (sbuffer[ii]) + offset;
        fwrite(fbuffer, sizeof(float), numread, outfile);
    }
    printf("done.\n Wrote %lld points.\n\n", N);
    fclose(infile);
    fclose(outfile);
    free(outname);
    exit(0);
}
