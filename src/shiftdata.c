/*                                    */
/* shiftdata:  Data shifting routine. */
/*        by Scott M. Ransom          */
/*            1 Mar 1999              */
/*                                    */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define WORKLEN 16384

void fshift(float *indata, float *outdata, int arrlen, double shift, double overlap);

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    int numread;
    double shift;
    float indata[WORKLEN], outdata[WORKLEN + 1], overlap = 0.0;
    char *infilenm = NULL, *outfilenm = NULL;
    FILE *infile = NULL, *outfile = NULL, *textout = NULL;

    /*  Get the filename of the datafile and open it */

    textout = (argc == 3) ? stderr : stdout;
    if (argc == 3 || argc == 4) {
        shift = strtod(argv[1], (char **) NULL);
        if (shift <= -1.0 || shift >= 1.0) {
            fprintf(textout, "\nbins_to_shift must be between -1.0 and 1.0.");
            fprintf(textout, "  Exiting.\n\n");
            exit(-1);
        }
        infilenm = (char *) malloc(strlen(argv[2]) + 1);
        strcpy(infilenm, argv[2]);
        fprintf(textout, "\nShifting data from \"%s\" by %f bin.\n\n",
                infilenm, shift);
        if (NULL == (infile = fopen(infilenm, "rb"))) {
            fprintf(textout, "\nCan't open \"%s\".  Exiting.\n\n", infilenm);
            exit(-1);
        }

        /* Open an output file if needed, else STDOUT */

        if (argc == 4) {
            outfilenm = (char *) malloc(strlen(argv[3]) + 1);
            strcpy(outfilenm, argv[3]);
            fprintf(textout, "Shifted data will be placed in \"%s\".\n\n",
                    outfilenm);
            if (NULL == (outfile = fopen(outfilenm, "wb"))) {
                fprintf(textout, "\nCan't open \"%s\".  Exiting.\n\n", outfilenm);
                exit(-1);
            }
        } else {
            outfile = stdout;
        }

    } else {

        printf("\nUsage:  shiftdata bins_to_shift infile [outfile]\n\n");
        printf("   This routine is used to shift a file of single\n");
        printf("   precision floating point data to the left or right\n");
        printf("   by up to 1 bin.  This allows you to patch two time\n");
        printf("   series together with the correct phase, assuming\n");
        printf("   that the data point durations are the same.\n");
        printf("   The arguments are:\n");
        printf("      bins_to_shift (double):  The amount to shift the\n");
        printf("         data.  Negative shifts left, positive shifts\n");
        printf("         right.  For example, if I have the time series\n");
        printf("         {1.0, 2,0, 1.0} and I give a bins_to_shift of\n");
        printf("         -0.3, my new time series would look like\n");
        printf("         {0.3, 1.3, 1.7, 0.7}.  Notice that my time\n");
        printf("         series increased in length by one and that I\n");
        printf("         added some noise due to the fact that I am\n");
        printf("         performing a type of two-bin average.\n");
        printf("      infile (string):  The data file to shift.\n");
        printf("      outfile (string):  An optional output file.  The\n");
        printf("         default behavior sends the raw output to STDOUT\n\n");
        exit(0);
    }

    /* Read the data in WORKLEN sized chunks */

    do {
        numread = fread(indata, sizeof(float), WORKLEN, infile);
        fshift(indata, outdata, numread, shift, overlap);
        overlap = outdata[numread];
        fwrite(outdata, sizeof(float), (unsigned long) numread, outfile);
    } while (numread == WORKLEN);

    /* Write the last data point if we didn't shift at all */

    if (fabs(shift) > 1.0e-7)
        fwrite(&overlap, sizeof(float), 1, outfile);

    /* Cleanup */

    free(infilenm);
    fclose(infile);
    if (argc == 4) {
        fclose(outfile);
        free(outfilenm);
    }
    exit(0);
}


void fshift(float *indata, float *outdata, int arrlen, double shift, double overlap)
{
    static double lopart, hipart;
    static int firsttime = 1;
    int ii;

    if (firsttime) {
        if (shift < 0.0) {
            lopart = -shift;
            hipart = 1.0 + shift;
        } else {
            lopart = 1.0 - shift;
            hipart = shift;
        }
        firsttime = 0;
    }
    outdata[0] = overlap + lopart * indata[0];
    outdata[arrlen] = hipart * indata[arrlen - 1];
    for (ii = 1; ii < arrlen; ii++)
        outdata[ii] = lopart * indata[ii] + hipart * indata[ii - 1];
}
