#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RECLEN 49792
#define NUMBEAMS 13
#define ALLBEAMS RECLEN*NUMBEAMS

int split_root_suffix(char *input, char **root, char **suffix)
/* This routine splits an input string into a root name */
/* + suffix.  Since it allocates the memory for the     */
/* root and suffix dynamically, the calling program     */
/* must free both "root" and "suffix".                  */
/* If the routine finds a suffix, it returns 1, else 0. */
{
    char *sptr = NULL;
    unsigned int len, rootlen = 0, suffixlen = 0;

    len = strlen(input);
    sptr = strrchr(input, '.');
    if (sptr == NULL) {
        *root = (char *) calloc(len + 1, sizeof(char));
        strncpy(*root, input, len);
        return 0;
    } else {
        rootlen = sptr - input;
        *root = (char *) calloc(rootlen + 1, sizeof(char));
        strncpy(*root, input, rootlen);
        suffixlen = len - rootlen - 1;
        *suffix = (char *) calloc(suffixlen + 1, sizeof(char));
        strncpy(*suffix, sptr + 1, suffixlen);
        return 1;
    }
}

int main(int argc, char *argv[])
{
    FILE *infile, *outfiles[NUMBEAMS];
    int ii, numrec = 0;
    char *outname, *root = NULL, *suffix = NULL;
    unsigned char buffer[ALLBEAMS];

    if (argc != 2) {
        printf("\n  Usage:  split_parkes_beams infilename\n\n");
        printf
            ("  This routine splits all 13 beams from a Parkes Multibeam raw data file.\n\n");
        exit(0);
    }

    split_root_suffix(argv[1], &root, &suffix);
    if (suffix)
        outname = calloc(strlen(root) + strlen(suffix) + 5, 1);
    else
        outname = calloc(strlen(root) + 5, 1);
    for (ii = 0; ii < NUMBEAMS; ii++) {
        if (suffix)
            sprintf(outname, "%s_%02d.%s", root, ii + 1, suffix);
        else
            sprintf(outname, "%s_%02d", root, ii + 1);
        outfiles[ii] = fopen(outname, "wb");
    }
    free(outname);
    free(root);
    free(suffix);

    infile = fopen(argv[1], "rb");
    while (fread(buffer, ALLBEAMS, 1, infile)) {
        for (ii = 0; ii < NUMBEAMS; ii++)
            fwrite(buffer + ii * RECLEN, RECLEN, 1, outfiles[ii]);
        numrec++;
    }
    fclose(infile);

    for (ii = 0; ii < NUMBEAMS; ii++)
        fclose(outfiles[ii]);
    printf("\nDone.  Split %d records.\n", numrec);
    exit(0);
}
