#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{

    unsigned char buffer[4];
    int tmp;
    long ct = 0;
    char innm[40], outnm[40];
    FILE *infile, *outfile;

    /*  Get the filename of the datafile and open it */

    printf("\n\nSwap a Big-Endian Data file to Little Endian\n");
    printf("          or vise-versa.  (32 bit data only)\n");
    if (argc > 1) {
        sprintf(innm, "%s", argv[1]);
        printf("\nReading data from \"%s\".\n", innm);
        if (NULL == (infile = fopen(innm, "r+b"))) {
            printf("\nCan't open \"%s\", exiting.\n", innm);
            exit(1);
        }
        sprintf(outnm, "%s_swapped", argv[1]);
        printf("Writing data to \"%s\".\n\n", outnm);
        if (NULL == (outfile = fopen(outnm, "w+b"))) {
            printf("\nCan't open %s, exiting.\n", outnm);
            exit(1);
        }
    } else {
        printf("\n\nPlease enter a filename for the input data file");
        printf(" after \"readdata\".\n");
        exit(1);
    }

    ct = 0;
    while (fread(buffer, sizeof(unsigned char), 4, infile)) {
        tmp = buffer[0];
        buffer[0] = buffer[3];
        buffer[3] = tmp;
        tmp = buffer[1];
        buffer[1] = buffer[2];
        buffer[2] = tmp;
        fwrite(buffer, sizeof(unsigned char), 4, outfile);
        ct++;
    }
    printf("Converted %ld points.\n\n", ct);
    fclose(infile);
    fclose(outfile);
    exit(0);
}
