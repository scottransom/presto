#include "presto.h"
#include "mask.h"
#include "multibeam.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
/* This routine checks a Parkes Multibeam raw data file for */
/* obvious errors.                                          */
{
    FILE *infile;
    long numrec = 1, currentctr = 1, lastctr = 1;
    int h = 0, m = 0, padding = 0;
    double currentut = 0.0, lastut = 0.0, s = 0.0;
    PKMB_tapehdr hdr;
    unsigned char data[RECLEN];
    char tmp[100];

    if (argc != 2) {
        printf("\n  Usage:  check_parkes_raw filename\n\n");
        printf("  This routine checks a Parkes Multibeam raw data file for ");
        printf("obvious errors.\n\n");
        exit(0);
    }

    printf("\n\n     Parkes Multibeam Raw Data File Checker\n");
    printf("                by Scott M. Ransom\n");
    printf("                    20 July 1998\n\n");

    infile = chkfopen(argv[1], "rb");
    while (read_PKMB_rawblock(&infile, 1, &hdr, data, &padding)) {
        sprintf(tmp, " %.8s ", hdr.blk_cntr);
        currentctr = strtol(tmp, NULL, 10);
        sprintf(tmp, " %.16s ", hdr.ut_blk);
        ra_dec_from_string(tmp, &h, &m, &s);
        currentut = hms2hours(h, m, s);
        if (numrec == 1) {
            sprintf(tmp, " %.16s ", hdr.ut_start);
            printf("Observation start time (UT) = %s\n", tmp);
            sprintf(tmp, " %.16s ", hdr.ut_blk);
            printf("First block time (UT)       = %s\n", tmp);
            printf("The first block number is %ld.\n\n", currentctr);
            if (currentctr != 1) {
                printf("\nWarning:  Possibly missing first %ld blocks.\n",
                       currentctr - 1);
            }
        } else {
            if ((currentut - lastut) > 0.001) {
                printf("\nWarning:  UT times make a (larger than normal ;) jump ");
                printf("between blocks %ld and %ld\n", lastctr, currentctr);
            }
            if (currentctr - lastctr != 1) {
                printf("\nWarning:  Record(s) missing between blocks %ld and %ld\n",
                       lastctr, currentctr);
            }
        }
        lastctr = currentctr;
        lastut = currentut;
        numrec++;
    }
    printf("\nThe file contains %ld records.\n", numrec - 1);
    printf("Looks OK.\n\n");
    fclose(infile);
    exit(0);
}
