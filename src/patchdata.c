/*                                            */
/* patchdata:  Data patching/padding routine. */
/*           by Scott M. Ransom               */
/*               1 Mar 1999                   */
/*                                            */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WORKLEN 16384

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    long ii, bins_to_patch;
    float *a, patchvalue;
    char *outfilenm = NULL;
    FILE *outfile, *textout = NULL;

    /*  Get the filename of the datafile and open it */

    textout = (argc == 3) ? stderr : stdout;
    if (argc == 3 || argc == 4) {
        bins_to_patch = strtol(argv[1], (char **) NULL, 10);
        if (bins_to_patch <= 0) {
            fprintf(textout, "\nbins_to_patch must be greater than 0.");
            fprintf(textout, "  Exiting.\n\n");
            exit(-1);
        }
        patchvalue = strtod(argv[2], (char **) NULL);

        /* Open an output file if needed, else STDOUT */

        if (argc == 4) {
            outfilenm = (char *) malloc(strlen(argv[3]) + 1);
            strcpy(outfilenm, argv[3]);
            fprintf(textout,
                    "\nPlacing %ld bins of value %f at the end of \"%s\".\n\n",
                    bins_to_patch, patchvalue, outfilenm);
            if (NULL == (outfile = fopen(outfilenm, "ab"))) {
                fprintf(textout, "\nCan't open \"%s\".  Exiting.\n\n", outfilenm);
                exit(-1);
            }
        } else {
            outfile = stdout;
        }

    } else {

        printf("\nUsage:  patchdata bins_to_patch patch_value [outfile]\n\n");
        printf("   This routine is used to generate floating point data\n");
        printf("   of a single value -- 'patch_value' -- to be used in\n");
        printf("   patching together or padding a data file.\n");
        printf("   The arguments are:\n");
        printf("      bins_to_patch (long):  The number of points to\n");
        printf("         generate.  This value must be greater than 0.\n");
        printf("      patch_value (float):  The value of each point.\n");
        printf("      outfile (string):  An optional output file to append\n");
        printf("         the padding to.  The default behavior sends the raw\n");
        printf("         output to STDOUT\n\n");
        exit(0);
    }

    /* Setup the data buffer */

    a = (float *) malloc(sizeof(float) * WORKLEN);
    for (ii = 0; ii < WORKLEN; ii++)
        a[ii] = patchvalue;

    /* Write the file */

    for (ii = 0; ii < bins_to_patch / WORKLEN; ii++)
        fwrite(a, sizeof(float), WORKLEN, outfile);

    /* Write the remnants */

    fwrite(a, sizeof(float), (unsigned long) (bins_to_patch % WORKLEN), outfile);

    /* Cleanup */

    free(a);
    if (argc == 4) {
        fclose(outfile);
        free(outfilenm);
    }
    exit(0);
}
