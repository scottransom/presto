#include "presto.h"

#define PRUNELEV 25
#define NEWLEV   5

/* Function declarations */

int compare_powindex(const void *ca, const void *cb);
   /*  Used as compare function for qsort() */


int prune_powers(float *arr, int n, int numsumpows)
/* Sets powers that are more than approx PRUNELEV standard */
/* devs above the median value to NEWLEV times the median  */
/* value.  This helps 'clean' the spectrum of high power   */
/* signals that probably have nothing to do with a phase   */
/* modulation spectrum (i.e. they are RF noise or strong   */
/* solitary pulsars.                                       */
{
  int ii, ct = 0;
  float median, cutoff, *tmparr;

  /* Determine the median power */

  tmparr = gen_fvect(n);
  memcpy(tmparr, arr, sizeof(float) * n);
  median = selectkth(n / 2, n, tmparr - 1);
  free(tmparr);

  /* Throw away powers that are bigger that PRUNELEV * median */

  cutoff = median * PRUNELEV / sqrt((float) numsumpows);
  for (ii = 0; ii < n; ii++) {
    if (arr[ii] > cutoff) {
      arr[ii] = NEWLEV * median;
      ct++;
    }
  }
  return ct;
}


float selectkth(long k, long n, float arr[])
/* Selects the kth largest value from the array arr */
{
  long i, ir, j, l, mid;
  float a, tempzz = 0.0;

  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
	SWAP(arr[l], arr[ir]);
      }
      return arr[k];
    } else {
      mid = (l + ir) >> 1;
      SWAP(arr[mid], arr[l + 1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l], arr[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
	SWAP(arr[l + 1], arr[ir]);
      }
      if (arr[l] > arr[l + 1]) {
	SWAP(arr[l], arr[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      for (;;) {
	do
	  i++;
	while (arr[i] < a);
	do
	  j--;
	while (arr[j] > a);
	if (j < i)
	  break;
	SWAP(arr[i], arr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      if (j >= k)
	ir = j - 1;
      if (j <= k)
	l = i;
    }
  }
}


void hpselect(unsigned long m, unsigned long n, \
	      float arr[], powindex heap[])
/* Selects the m largest values from the array arr          */
/* and stores them and their indices in heap and heapindex. */
{
  unsigned long i, j, k;
  powindex tempzz;

  if (m > n / 2 || m < 1) {
    printf("Probable misuse of hpselect.\n");
    printf("Number to select is out of range.  Exiting.\n");
    exit(1);
  }
  for (i = 1; i <= m; i++) {
    heap[i].pow = arr[i];
    heap[i].ind = i - 1;
  }
  qsort(heap + 1, m, sizeof(powindex), compare_powindex);
  for (i = m + 1; i <= n; i++) {
    if (arr[i] > heap[1].pow) {
      heap[1].pow = arr[i];
      heap[1].ind = i - 1;
      for (j = 1;;) {
	k = j << 1;
	if (k > m)
	  break;
	if (k != m && heap[k].pow > heap[k + 1].pow)
	  k++;
	if (heap[j].pow <= heap[k].pow)
	  break;
	SWAP(heap[k], heap[j]);
	j = k;
      }
    }
  }
}
