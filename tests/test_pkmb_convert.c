#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

extern long int random(void);
double time_diff(struct timeval t1, struct timeval t2);
int check_results(int nn, unsigned char *gooddata, unsigned char *data);
void convert1(unsigned char *rec, unsigned char *data,
	      int numchan, int decreasing_f);
void convert2(unsigned char *rec, unsigned char *data,
	      int numchan, int decreasing_f);

#define NUMRAW  32000000
#define NUMDATA 256000000

double time_diff(struct timeval tend, struct timeval tbeg)
{
  struct timeval diff;
  
  diff.tv_sec = tend.tv_sec - tbeg.tv_sec;
  diff.tv_usec = tend.tv_usec - tbeg.tv_usec;
  /* normalize */
  while (diff.tv_usec < 0) {
    diff.tv_usec += 1000000L;
    diff.tv_sec -= 1;
  }
  return diff.tv_sec + (double) diff.tv_usec / 1000000.0;
}

#define GET_BIT(c, n) (*(c+(n>>3)) >> (7-(n&7)) & 1)
#define CHARBIT(c, n) (*(c) >> (7-(n&7)) & 1)

int main(void)
{
  int ii, jj;
  long *rawptr;
  unsigned char *rawdata;
  unsigned char *gooddata;
  unsigned char *data;
  struct timeval t1, t2;
  
  rawdata = (unsigned char *)malloc(NUMRAW);
  gooddata = (unsigned char *)malloc(NUMDATA);
  data = (unsigned char *)malloc(NUMDATA);

  /* Set the raw data */

  rawptr = (long *) rawdata;
  for (ii = 0; ii < NUMRAW/4; ii++)
    rawptr[ii] = random();

  printf("Increasing freqs...\n");

  /* Do the 'correct' conversion */

  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, gooddata+jj, 256, 0);

  /* Do the 'correct' conversion */

  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, data+jj, 256, 0);

  /* Time the first one */

  gettimeofday(&t1, 0);
  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, data+jj, 256, 0);
  gettimeofday(&t2, 0);
  check_results(NUMDATA, gooddata, data);
  printf("     Standard Time = %.4f s\n", time_diff(t2, t1));

  /* Time the second one */
			     
  gettimeofday(&t1, 0);
  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert2(rawdata+ii, data+jj, 256, 0);
  gettimeofday(&t2, 0);
  check_results(NUMDATA, gooddata, data);
  printf("   New Method Time = %.4f s\n", time_diff(t2, t1));

  printf("Decreasing freqs...\n");

  /* Do the 'correct' conversion */

  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, gooddata+jj, 256, 1);

  /* Do the 'correct' conversion */

  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, data+jj, 256, 1);

  /* Time the first one */

  gettimeofday(&t1, 0);
  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert1(rawdata+ii, data+jj, 256, 1);
  gettimeofday(&t2, 0);
  check_results(NUMDATA, gooddata, data);
  printf("     Standard Time = %.4f s\n", time_diff(t2, t1));

  /* Time the second one */
			     
  gettimeofday(&t1, 0);
  for (ii = 0, jj = 0; ii < NUMRAW/32; ii+=32, jj+=256)
    convert2(rawdata+ii, data+jj, 256, 1);
  gettimeofday(&t2, 0);
  check_results(NUMDATA, gooddata, data);
  printf("   New Method Time = %.4f s\n", time_diff(t2, t1));

  free(rawdata);
  free(gooddata);
  free(data);

  return 0;
}


int check_results(int nn, unsigned char *gooddata, unsigned char *data)
/* Return value of 1 is bad, 0 is good (they match)  */
{
  int ii, flag=0;

  for (ii = 0; ii < nn; ii++)
    if (gooddata[ii] != data[ii]){
      flag = 1;
      printf("gooddata, data = %d, %d at %d\n", 
	     gooddata[ii], data[ii], ii);
    }
  return flag;
}


void convert_old(unsigned char *rec, unsigned char *data, \
	      int numchan, int decreasing_f)
/* This routine converts 1 bit digitized data with 'numchan' */
/* channels to an array of 'numchan' floats.                 */
{
  int ii, jj;

  if (decreasing_f){
    for(ii = 0, jj = numchan-8; ii < numchan / 8; ii++, jj-=8){
      data[jj] = rec[ii] & 0x80 ? 1 : 0;
      data[jj+1] = rec[ii] & 0x40 ? 1 : 0;
      data[jj+2] = rec[ii] & 0x20 ? 1 : 0;
      data[jj+3] = rec[ii] & 0x10 ? 1 : 0;
      data[jj+4] = rec[ii] & 0x08 ? 1 : 0;
      data[jj+5] = rec[ii] & 0x04 ? 1 : 0;
      data[jj+6] = rec[ii] & 0x02 ? 1 : 0;
      data[jj+7] = rec[ii] & 0x01 ? 1 : 0;
    }
  } else {
    for(ii = 0, jj = 0; ii < numchan / 8; ii++, jj+=8){
      data[jj] = rec[ii] & 0x01 ? 1 : 0;
      data[jj+1] = rec[ii] & 0x02 ? 1 : 0;
      data[jj+2] = rec[ii] & 0x04 ? 1 : 0;
      data[jj+3] = rec[ii] & 0x08 ? 1 : 0;
      data[jj+4] = rec[ii] & 0x10 ? 1 : 0;
      data[jj+5] = rec[ii] & 0x20 ? 1 : 0;
      data[jj+6] = rec[ii] & 0x40 ? 1 : 0;
      data[jj+7] = rec[ii] & 0x80 ? 1 : 0;
    }
  }
}

void convert1(unsigned char *rec, unsigned char *data, \
	      int numchan, int decreasing_f)
/* This routine converts 1 bit digitized data with 'numchan' */
/* channels to an array of 'numchan' floats.                 */
{
  int ii, jj;

  if (decreasing_f){
    for(ii = numchan/8-1, jj = 0; ii >= 0; ii--, jj+=8){
      data[jj] = (rec[ii] >> 7) & 1;
      data[jj+1] = (rec[ii] >> 6) & 1;
      data[jj+2] = (rec[ii] >> 5) & 1;
      data[jj+3] = (rec[ii] >> 4) & 1;
      data[jj+4] = (rec[ii] >> 3) & 1;
      data[jj+5] = (rec[ii] >> 2) & 1;
      data[jj+6] = (rec[ii] >> 1) & 1;
      data[jj+7] = rec[ii] & 1;
    }
  } else {
    for(ii = 0, jj = 0; ii < numchan/8; ii++, jj+=8){
      data[jj] = rec[ii] & 1;
      data[jj+1] = (rec[ii] >> 1) & 1;
      data[jj+2] = (rec[ii] >> 2) & 1;
      data[jj+3] = (rec[ii] >> 3) & 1;
      data[jj+4] = (rec[ii] >> 4) & 1;
      data[jj+5] = (rec[ii] >> 5) & 1;
      data[jj+6] = (rec[ii] >> 6) & 1;
      data[jj+7] = (rec[ii] >> 7) & 1;
    }
  }
}

void convert2(unsigned char *rec, unsigned char *data, \
	      int numchan, int decreasing_f)
/* This routine converts 1 bit digitized data with 'numchan' */
/* channels to an array of 'numchan' floats.                 */
{
  register unsigned char one = 1;
  register unsigned char two = 2;
  register unsigned char thr = 3;
  register unsigned char fou = 4;
  register unsigned char fiv = 5;
  register unsigned char six = 6;
  register unsigned char sev = 7;
  int ii, jj;

  if (decreasing_f){
    for(ii = numchan/8-1, jj = 0; ii >= 0; ii--, jj+=8){
      data[jj] = (rec[ii] >> sev) & one;
      data[jj+1] = (rec[ii] >> six) & one;
      data[jj+2] = (rec[ii] >> fiv) & one;
      data[jj+3] = (rec[ii] >> fou) & one;
      data[jj+4] = (rec[ii] >> thr) & one;
      data[jj+5] = (rec[ii] >> two) & one;
      data[jj+6] = (rec[ii] >> one) & one;
      data[jj+7] = rec[ii] & one;
    }
  } else {
    for(ii = 0, jj = 0; ii < numchan/8; ii++, jj+=8){
      data[jj] = rec[ii] & one;
      data[jj+1] = (rec[ii] >> 1) & one;
      data[jj+2] = (rec[ii] >> 2) & one;
      data[jj+3] = (rec[ii] >> 3) & one;
      data[jj+4] = (rec[ii] >> 4) & one;
      data[jj+5] = (rec[ii] >> 5) & one;
      data[jj+6] = (rec[ii] >> 6) & one;
      data[jj+7] = (rec[ii] >> 7) & one;
    }
  }
}

