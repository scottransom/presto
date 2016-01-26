#include "presto.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    int ii, numfiles, filelen, datalen;
    multifile *mfile;
    float *data, *scratch;
    char **filenames;

    srand(3);
    if (argc < 4) {
        printf("usage:  test_multifiles numfiles filelen datalen\n");
        exit(0);
    }
    numfiles = atoi(argv[1]);
    filenames = (char **) malloc(sizeof(char *) * numfiles);
    filelen = atoi(argv[2]);
    for (ii = 0; ii < numfiles; ii++) {
        filenames[ii] = calloc(10, 1);
        sprintf(filenames[ii], "test%d.dat", ii);
    }

    /* Generate the data and scratch spaces */

    datalen = atoi(argv[3]);
    data = gen_fvect(datalen);
    scratch = gen_fvect(datalen);
    for (ii = 0; ii < datalen; ii++)
        data[ii] = (float) 2.0 *(rand() / ((double) RAND_MAX)) - 1.0;

    /* Write the data */

    printf("\nWrite test...\n");
    mfile = fopen_multifile(numfiles, filenames, "w", filelen);
    print_multifile(mfile, 1);
    fwrite_multifile(data, sizeof(float), datalen, mfile);
    print_multifile(mfile, 0);
    fclose_multifile(mfile);

    /* Read the data */

    printf("\nRead test...\n");
    mfile = fopen_multifile(numfiles, filenames, "r", 0);
    print_multifile(mfile, 1);
    fread_multifile(scratch, sizeof(float), datalen, mfile);
    print_multifile(mfile, 0);
    fclose_multifile(mfile);

    /* Compare data and scratch */

    for (ii = 0; ii < datalen; ii++) {
        if (scratch[ii] != data[ii])
            printf("Data point %d doesn't match:  %f, %f\n", ii, data[ii],
                   scratch[ii]);
    }

    /* Seek to a random place in the file and check that data */

    printf("\nSeek test...\n");
    {
        long long byteindex;

        mfile = fopen_multifile(numfiles, filenames, "r", 0);
        print_multifile(mfile, 1);
        byteindex = mfile->filelens[0] + mfile->filelens[1]
            - sizeof(float) * 10;
        print_multifile(mfile, 0);
        fseek_multifile(mfile, byteindex, SEEK_SET);
        print_multifile(mfile, 0);
        fread_multifile(scratch + byteindex / sizeof(float), sizeof(float), 30,
                        mfile);
        print_multifile(mfile, 0);
        fclose_multifile(mfile);
    }

    /* Compare data and scratch */

    for (ii = 0; ii < datalen; ii++) {
        if (scratch[ii] != data[ii])
            printf("Data point %d doesn't match:  %f, %f\n", ii, data[ii],
                   scratch[ii]);
    }
    /* Free stuff up */

    vect_free(data);
    vect_free(scratch);
    for (ii = 0; ii < numfiles; ii++)
        free(filenames[ii]);
    free(filenames);
    exit(0);
}
