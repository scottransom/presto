#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HDRLEN 640
#define DATLEN 49152

int main(int argc, char *argv[])
{
    FILE *datfile, *rawfile;
    char hdr[HDRLEN + 8], filenm[100];
    unsigned char dat[DATLEN + 8];
    int ii = 0, numread;

    if (argc != 2) {
        printf("usage:  un_sc_td basefilename\n");
        exit(0);
    }

    sprintf(filenm, "%s.hdr", argv[1]);
    datfile = fopen(filenm, "r");
    if ((numread = fread(hdr, sizeof(char), HDRLEN + 8, datfile)) != HDRLEN + 8) {
        printf("  Problem reading the header file '%s'\n", argv[1]);
    } else {
        printf("Successfully read the header file '%s'.\n"
               "  Now merging '.hdr' and '.dat' files...\n", argv[1]);
    }
    fclose(datfile);

    sprintf(filenm, "%s.dat", argv[1]);
    datfile = fopen(filenm, "r");
    sprintf(filenm, "%s.pkmb", argv[1]);
    rawfile = fopen(filenm, "w");
    while (fread(dat, sizeof(unsigned char), DATLEN + 8, datfile)) {
        fwrite(hdr + 4, sizeof(char), HDRLEN, rawfile);
        fwrite(dat + 4, sizeof(unsigned char), DATLEN, rawfile);
        ii++;
    }
    fclose(datfile);
    fclose(rawfile);
    printf("Done.  Wrote %d blocks.\n", ii);
    return 0;
}
