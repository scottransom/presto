#include "presto.h"
#include "mask.h"

void fill_mask(double sigma, double mjd, double dtint, double lofreq, 
	       double dfreq, int numchan, int numint, int ptsperint, 
	       int num_user_chans, int *user_chans, int num_user_ints, 
	       int *user_ints, unsigned char **bytemask, mask *obsmask)
/* Fill a mask structure with the appropriate values */
{
  int ii, jj, count;

  obsmask->sigma = sigma;
  obsmask->mjd = mjd;
  obsmask->dtint = dtint;
  obsmask->lofreq = lofreq;
  obsmask->dfreq = dfreq;
  obsmask->numchan = numchan;
  obsmask->numint = numint;
  obsmask->ptsperint = ptsperint;
  obsmask->num_user_chans = num_user_chans;
  if (obsmask->num_user_chans){
    obsmask->user_chans = gen_ivect(obsmask->num_user_chans);
    for (ii=0; ii<obsmask->num_user_chans; ii++)
      obsmask->user_chans[ii] = user_chans[ii];
  }
  obsmask->num_user_ints = num_user_ints;
  if (obsmask->num_user_ints){
    obsmask->user_ints = gen_ivect(obsmask->num_user_ints);
    for (ii=0; ii<obsmask->num_user_ints; ii++)
      obsmask->user_ints[ii] = user_ints[ii];
  }
  obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
  obsmask->chans = (int **)malloc(obsmask->numint * sizeof(int *));
  for (ii=0; ii<obsmask->numint; ii++){
    /* Count the bad channels first */
    count = 0;
    for (jj=0; jj<obsmask->numchan; jj++)
      if (bytemask[ii][jj]) count++;
    obsmask->num_chans_per_int[ii] = count;
    if (count){
      /* Now determine which channels */
      count = 0;
      obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
      for (jj=0; jj<obsmask->numchan; jj++){
	if (bytemask[ii][jj]){
	  obsmask->chans[ii][count] = jj;
	  count++;
	}
      }
    }
  }
}


void free_mask(mask obsmask)
/* Free the contents of an mask structure */
{
  int ii;

  for (ii=0; ii<obsmask.numint; ii++){
    if (obsmask.num_chans_per_int[ii] > 0 &&
	obsmask.num_chans_per_int[ii] < obsmask.numchan)
      free(obsmask.chans[ii]);
  }
  free(obsmask.chans);
  free(obsmask.num_chans_per_int);
  if (obsmask.num_user_chans)
    free(obsmask.user_chans);
  if (obsmask.num_user_ints)
    free(obsmask.user_ints);
}

void read_mask(char *maskfilenm, mask *obsmask)
/* Read the contents of a mask structure from a file */
{
  FILE *infile;
  int ii;

  infile = chkfopen(maskfilenm, "rb");
  chkfread(&(obsmask->sigma), sizeof(double), 1, infile);
  chkfread(&(obsmask->mjd), sizeof(double), 1, infile);
  chkfread(&(obsmask->dtint), sizeof(double), 1, infile);
  chkfread(&(obsmask->lofreq), sizeof(double), 1, infile);
  chkfread(&(obsmask->dfreq), sizeof(double), 1, infile);
  chkfread(&(obsmask->numchan), sizeof(int), 1, infile);
  chkfread(&(obsmask->numint), sizeof(int), 1, infile);
  chkfread(&(obsmask->ptsperint), sizeof(int), 1, infile);
  chkfread(&(obsmask->num_user_chans), sizeof(int), 1, infile);
  if (obsmask->num_user_chans){
    obsmask->user_chans = gen_ivect(obsmask->num_user_chans);
    chkfread(obsmask->user_chans, sizeof(int), 
	     obsmask->num_user_chans, infile);
  }
  chkfread(&(obsmask->num_user_ints), sizeof(int), 1, infile);
  if (obsmask->num_user_ints){
    obsmask->user_ints = gen_ivect(obsmask->num_user_ints);
    chkfread(obsmask->user_ints, sizeof(int), 
	     obsmask->num_user_ints, infile);
  }
  obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
  chkfread(obsmask->num_chans_per_int, sizeof(int), 
	   obsmask->numint, infile);
  obsmask->chans = (int **)malloc(obsmask->numint * sizeof(int *));
  for (ii=0; ii<obsmask->numint; ii++){
    if (obsmask->num_chans_per_int[ii] > 0 &&
	obsmask->num_chans_per_int[ii] < obsmask->numchan){
      obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
      chkfread(obsmask->chans[ii], sizeof(int), 
	       obsmask->num_chans_per_int[ii], infile);
    }
  }
  fclose(infile);
}

void write_mask(char *maskfilenm, mask *obsmask)
/* Write the contents of an mask structure to a file */
{
  FILE *outfile;
  int ii;

  outfile = chkfopen(maskfilenm, "wb");
  chkfwrite(&(obsmask->sigma), sizeof(double), 1, outfile);
  chkfwrite(&(obsmask->mjd), sizeof(double), 1, outfile);
  chkfwrite(&(obsmask->dtint), sizeof(double), 1, outfile);
  chkfwrite(&(obsmask->lofreq), sizeof(double), 1, outfile);
  chkfwrite(&(obsmask->dfreq), sizeof(double), 1, outfile);
  chkfwrite(&(obsmask->numchan), sizeof(int), 1, outfile);
  chkfwrite(&(obsmask->numint), sizeof(int), 1, outfile);
  chkfwrite(&(obsmask->ptsperint), sizeof(int), 1, outfile);
  chkfwrite(&(obsmask->num_user_chans), sizeof(int), 1, outfile);
  if (obsmask->num_user_chans)
    chkfwrite(obsmask->user_chans, sizeof(int), 
	      obsmask->num_user_chans, outfile);
  chkfwrite(&(obsmask->num_user_ints), sizeof(int), 1, outfile);
  if (obsmask->num_user_ints)
    chkfwrite(obsmask->user_ints, sizeof(int), 
	      obsmask->num_user_ints, outfile);
  chkfwrite(obsmask->num_chans_per_int, sizeof(int), 
	    obsmask->numint, outfile);
  for (ii=0; ii<obsmask->numint; ii++){
    if (obsmask->num_chans_per_int[ii] > 0 &&
	obsmask->num_chans_per_int[ii] < obsmask->numchan){
      chkfwrite(obsmask->chans[ii], sizeof(int), 
		obsmask->num_chans_per_int[ii], outfile);
    }
  }
  fclose(outfile);
}
