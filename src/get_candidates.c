#include "presto.h"

void get_rzw_cand(char *filenm, int candnum, fourierprops * cand)
/*  Read the rzw candidate file 'filenm' and return a        */
/*  pointer to the fourierprops that describes it.           */
{
  FILE *candfile;

  candfile = chkfopen(filenm, "rb");
  chkfileseek(candfile, candnum-1, sizeof(fourierprops), SEEK_SET);
  chkfread(cand, sizeof(fourierprops), 1, candfile);
  fclose(candfile);
}


void get_bin_cand(char *filenm, int candnum, binaryprops * cand)
/*  Read the bin candidate file 'filenm' and return a        */
/*  pointer to the binaryprops that describes it.            */
{
  FILE *candfile;

  candfile = chkfopen(filenm, "rb");
  chkfileseek(candfile, candnum-1, sizeof(binaryprops), SEEK_SET);
  chkfread(cand, sizeof(binaryprops), 1, candfile);
  fclose(candfile);
}



