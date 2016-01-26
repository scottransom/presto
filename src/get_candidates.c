#include "presto.h"

int read_rzw_cand(FILE * file, fourierprops * cands)
/* Read the next rzw candidate from the file */
/* If successful, return 1, else 0           */
{
    return chkfread(cands, sizeof(fourierprops), 1, file);
}


int read_bin_cand(FILE * file, binaryprops * cands)
/* Read the next binary candidate from the file */
/* If successful, return 1, else 0              */
{
    return chkfread(cands, sizeof(binaryprops), 1, file);
}


int read_rawbin_cand(FILE * file, rawbincand * cands)
/* Read the next rawbin candidate from the file */
/* If successful, return 1, else 0              */
{
    return chkfread(cands, sizeof(rawbincand), 1, file);
}


void get_rzw_cand(char *filenm, int candnum, fourierprops * cand)
/*  Read the rzw candidate file 'filenm' and return a        */
/*  pointer to the fourierprops that describes it.           */
{
    FILE *candfile;

    candfile = chkfopen(filenm, "rb");
    chkfileseek(candfile, candnum - 1, sizeof(fourierprops), SEEK_SET);
    chkfread(cand, sizeof(fourierprops), 1, candfile);
    fclose(candfile);
}


void get_bin_cand(char *filenm, int candnum, binaryprops * cand)
/*  Read the bin candidate file 'filenm' and return a        */
/*  pointer to the binaryprops that describes it.            */
{
    FILE *candfile;

    candfile = chkfopen(filenm, "rb");
    chkfileseek(candfile, candnum - 1, sizeof(binaryprops), SEEK_SET);
    chkfread(cand, sizeof(binaryprops), 1, candfile);
    fclose(candfile);
}


void get_rawbin_cand(char *filenm, int candnum, rawbincand * cand)
/*  Read the rawbin candidate file 'filenm' and return a     */
/*  pointer to the rawbincand that describe it.              */
{
    FILE *candfile;

    candfile = chkfopen(filenm, "rb");
    chkfileseek(candfile, candnum - 1, sizeof(rawbincand), SEEK_SET);
    chkfread(cand, sizeof(rawbincand), 1, candfile);
    fclose(candfile);
}
