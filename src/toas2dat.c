#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "toas2dat_cmd.h"

#define WORKLEN 65536
#define SECPERDAY 86400
#define MAXREAD 32768
/* #define DEBUG */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

unsigned long getfilelen(FILE * file, size_t size)
{
    int filenum, rt;
    struct stat buf;

    filenum = fileno(file);
    rt = fstat(filenum, &buf);
    if (rt == -1) {
        perror("\nError in getfilelen()");
        printf("\n");
        exit(-1);
    }
    return (unsigned long) (buf.st_size / size);
}

int compare_doubles(const void *a, const void *b)
/* qsort comparison function for doubles */
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}

int read_toas(FILE *infile, double **toas)
/* Read a text file containing ASCII text TOAs. */
/* The number of TOAs read is returned.         */
/* Lines beginning with '#' are ignored.        */
{
    double dtmp;
    char line[80], *sptr = NULL;
    int ii=0, numtoa=0;

    /* Read the input file once to count TOAs */
    while (1) {
        sptr = fgets(line, 80, infile);
        if (!feof(infile) && sptr != NULL && sptr[0] != '\n') {
            if (line[0] != '#' && sscanf(line, "%lf", &dtmp) == 1)
                numtoa++;
        } else {
            break;
        }
    }

    *toas = (double *) malloc(sizeof(double) * numtoa);

    /* Rewind and read the TOAs for real */
    rewind(infile);
    while (1) {
        sptr = fgets(line, 80, infile);
        if (!feof(infile) && sptr != NULL && sptr[0] != '\n') {
            if (line[0] != '#' && sscanf(line, "%lf", &(*toas)[ii]) == 1)
                ii++;
        } else {
            break;
        }
    }
    return numtoa;
}


int main(int argc, char *argv[])
/* Convert a file of TOAs in either text or binary format   */
/* into a floating point time series.  The time series will */
/* have 'cmd->numout' points with each bin of length        */
/* 'cmd->dt' seconds.                                       */
{
    long ii, jj, ntoas, numwrites, numtowrite, numplaced = 0;
    double To, toa, *toaptr, *ddata, lotime, hitime, dtfract, blockt;
    float *fdata;
    FILE *infile, *outfile;
    Cmdline *cmd;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

#ifdef DEBUG
    showOptionValues();
#endif

    fprintf(stderr, "\n\n  TOA to Time Series Converter\n");
    fprintf(stderr, "      by Scott M. Ransom\n");
    fprintf(stderr, "        17 October 2000\n\n");

    /* Open our files and read the TOAs */

    printf("\nReading TOAs from '%s'.\n", cmd->argv[0]);
    if (cmd->textP) {           /* Text data */
        infile = fopen(cmd->argv[0], "r");
        ntoas = read_toas(infile, &ddata);
        printf("   Found %ld TOAs.\n", ntoas);
    } else {                    /* Binary data */
        infile = fopen(cmd->argv[0], "rb");
        if (cmd->floatP) {      /* Floating point data */
            ntoas = getfilelen(infile, sizeof(float));
            printf("   Found %ld TOAs.\n", ntoas);
            ddata = (double *) malloc(sizeof(double) * ntoas);
            fdata = (float *) malloc(sizeof(float) * ntoas);
            jj = fread(fdata, sizeof(float), ntoas, infile);
            if (jj != ntoas) {
                printf("\nError reading TOA file.  Only %ld points read.\n\n", jj);
                exit(-1);
            }
            for (jj = 0; jj < ntoas; jj++)
                ddata[jj] = (double) fdata[jj];
            free(fdata);
        } else {                /* Double precision data */
            ntoas = getfilelen(infile, sizeof(double));
            printf("   Found %ld TOAs.\n", ntoas);
            ddata = (double *) malloc(sizeof(double) * ntoas);
            jj = fread(ddata, sizeof(double), ntoas, infile);
            if (jj != ntoas) {
                printf("\nError reading TOA file.  Only %ld points read.\n\n", jj);
                exit(-1);
            }
        }
    }
    fclose(infile);
    outfile = fopen(cmd->outfile, "wb");

    /* Allocate our output array */

    fdata = (float *) malloc(sizeof(float) * WORKLEN);
    printf("\nWriting time series of %ld points of\n", cmd->numout);
    printf("length %f seconds to '%s'.\n\n", cmd->dt, cmd->outfile);

    /* Sort the TOAs */

    qsort(ddata, ntoas, sizeof(double), compare_doubles);

    /* Convert the TOAs to seconds offset from the first TOA */

    if (cmd->t0P)
        To = cmd->t0;
    else
        To = ddata[0];
    if (cmd->secP)
        for (ii = 0; ii < ntoas; ii++)
            ddata[ii] = (ddata[ii] - To);
    else
        for (ii = 0; ii < ntoas; ii++)
            ddata[ii] = (ddata[ii] - To) * SECPERDAY;
    toaptr = ddata;
    toa = *toaptr;

    /* Determine the number of output writes we need */

    numwrites = (cmd->numout % WORKLEN) ?
        cmd->numout / WORKLEN + 1 : cmd->numout / WORKLEN;
    dtfract = 1.0 / cmd->dt;
    blockt = WORKLEN * cmd->dt;

    /* Loop over the number of writes */

    for (ii = 0; ii < numwrites; ii++) {

        /* Determine the beginning and ending times of the output array */

        lotime = ii * blockt;
        hitime = (ii + 1) * blockt;
        numtowrite = ((cmd->numout % WORKLEN) && (ii == (numwrites - 1))) ?
            cmd->numout % WORKLEN : WORKLEN;

        /* Initialize the output data array to all zeros */

        for (jj = 0; jj < WORKLEN; jj++)
            fdata[jj] = 0.0;

        /* Place any TOAs we need to in the current output array */

        while ((toaptr - ddata) < ntoas) {
            if (toa >= hitime)
                break;
            else if (toa >= lotime) {
                fdata[(int) ((toa - lotime) * dtfract)] += 1.0;
                numplaced++;
            }
            toaptr++;
            toa = *toaptr;
        }

        /* Write the output data */

        fwrite(fdata, sizeof(float), numtowrite, outfile);
    }

    /* Cleanup */

    printf("Done.\n   Placed %ld TOAs.\n\n", numplaced);
    fclose(outfile);
    free(fdata);
    free(ddata);
    exit(0);
}
