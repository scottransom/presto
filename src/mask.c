#include "presto.h"
#include "mask.h"

extern int compare_ints(const void *a, const void *b);

static find_num(int num, int *arr, int arrlen);
static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, 
			  int *merged);

void fill_mask(double sigma, double mjd, double dtint, double lofreq, 
	       double dfreq, int numchan, int numint, int ptsperint, 
	       int num_zap_chans, int *zap_chans, int num_zap_ints, 
	       int *zap_ints, unsigned char **bytemask, mask *obsmask)
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
  obsmask->num_zap_chans = num_zap_chans;
  if (obsmask->num_zap_chans){
    obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
    for (ii=0; ii<obsmask->num_zap_chans; ii++)
      obsmask->zap_chans[ii] = zap_chans[ii];
  }
  obsmask->num_zap_ints = num_zap_ints;
  if (obsmask->num_zap_ints){
    obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
    for (ii=0; ii<obsmask->num_zap_ints; ii++)
      obsmask->zap_ints[ii] = zap_ints[ii];
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
	if (bytemask[ii][jj])
	  obsmask->chans[ii][count++] = jj;
      }
    }
  }
}


void set_bytes_from_mask(mask *obsmask, unsigned char **bytematrix,
			 unsigned char fillval)
/* Inverse of fill_mask() */
{
  int ii, jj;

  for (ii=0; ii<obsmask->numint; ii++)
    for (jj=0; jj<obsmask->num_chans_per_int[ii]; jj++)
      bytarr[ii][obsmask->chans[ii][jj]] = fillval;
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
  if (obsmask.num_zap_chans)
    free(obsmask.zap_chans);
  if (obsmask.num_zap_ints)
    free(obsmask.zap_ints);
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
  chkfread(&(obsmask->num_zap_chans), sizeof(int), 1, infile);
  if (obsmask->num_zap_chans){
    obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
    chkfread(obsmask->zap_chans, sizeof(int), 
	     obsmask->num_zap_chans, infile);
  }
  chkfread(&(obsmask->num_zap_ints), sizeof(int), 1, infile);
  if (obsmask->num_zap_ints){
    obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
    chkfread(obsmask->zap_ints, sizeof(int), 
	     obsmask->num_zap_ints, infile);
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
  chkfwrite(&(obsmask->num_zap_chans), sizeof(int), 1, outfile);
  if (obsmask->num_zap_chans)
    chkfwrite(obsmask->zap_chans, sizeof(int), 
	      obsmask->num_zap_chans, outfile);
  chkfwrite(&(obsmask->num_zap_ints), sizeof(int), 1, outfile);
  if (obsmask->num_zap_ints)
    chkfwrite(obsmask->zap_ints, sizeof(int), 
	      obsmask->num_zap_ints, outfile);
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


int check_mask(double starttime, double duration, mask *obsmask, 
	       int *maskchans)
/* Return value is the number of channels to mask.  The */
/* channel numbers are placed in maskchans (which must  */
/* have a length of numchan).  If -1 is returned, all   */
/* channels should be masked.                           */
{
  int ii, loint, hiint, numchan;
  double endtime;
  static int old_loint=-1, old_hiint=-1, old_numchan=0;
  
  endtime = starttime + duration;
  loint = (int)(starttime / dtint);
  hiint = (int)(endtime / dtint);
  if (loint == old_loint && hiint == old_hiint)
    /* Notice we don't mess with the maskchans array! */
    return old_numchan;
  if (loint == hiint){
    old_loint = old_hiint = loint;
    if (obsmask->num_zap_ints)
      if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)){
	old_numchan = -1;
	return old_numchan;
      }
    old_numchan = merge_no_dupes(obsmask->zap_chans, 
				 obsmask->num_zap_chans,
				 obsmask->chans[loint], 
				 obsmask->num_chans_per_int[loint],
				 maskchans);
  } else {
    int *tmpchans;

    old_loint = loint;
    old_hiint = hiint;
    if (obsmask->num_zap_ints){
      if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)){
	old_numchan = -1;
	return old_numchan;
      }
      if (find_num(hiint, obsmask->zap_ints, obsmask->num_zap_ints)){
	old_numchan = -1;
	return old_numchan;
      }
      tmpchans = gen_ivect(obsmask->numchan);
      old_numchan = merge_no_dupes(obsmask->zap_chans, 
				   obsmask->num_zap_chans,
				   obsmask->chans[loint], 
				   obsmask->num_chans_per_int[loint],
				   tmpchans);
    } else {
      tmpchans = obsmask->chans[loint];
      old_numchan = obsmask->num_chans_per_int[loint];
    }
    old_numchan = merge_no_dupes(tmpchans, 
				 old_numchan,
				 obsmask->chans[hiint], 
				 obsmask->num_chans_per_int[hiint],
				 maskchans);
    if (obsmask->num_zap_ints)
      free(tmpchans);
  }
  return old_numchan;
}


static find_num(int num, int *arr, int arrlen)
{
  int ii;

  for (ii=0; ii<arrlen; ii++)
    if (arr[ii] >= num)
      return (arr[ii] == num) ? 1 : 0;
  return 0;
}


static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, 
			  int *merged)
{
  int ptr1=0, ptr2=0, count=0;
    
  while (1){
    if (arr1[ptr1] < arr2[ptr2])
      out[count++] = arr1[ptr1++];
    else if (arr1[ptr1] > arr2[ptr2])
      out[count++] = arr2[ptr2++];
    else {
      out[count++] = arr1[ptr1];
      ptr1++;
      ptr2++;
    }
    if (ptr1 == len1){
      while (ptr2 < len2)
	out[count++] = arr2[ptr2++];
      break;
    } else if (ptr2 == len2){
      while (ptr1 < len1)
	out[count++] = arr1[ptr1++];
      break;
    }
  }
  return count;
}
