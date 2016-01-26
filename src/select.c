#include "presto.h"

#define PRUNELEV 25
#define NEWLEV   5

/* Function declarations */

int compare_powindex(const void *ca, const void *cb);
   /*  Used as compare function for qsort() */
float median(float arr[], int n);
/* Finds the median (but messes up the array order) */

int prune_powers(float *arr, int n, int numsumpows)
/* Sets powers that are more than approx PRUNELEV standard */
/* devs above the median value to NEWLEV times the median  */
/* value.  This helps 'clean' the spectrum of high power   */
/* signals that probably have nothing to do with a phase   */
/* modulation spectrum (i.e. they are RF noise or strong   */
/* solitary pulsars.                                       */
{
    int ii, ct = 0;
    float med, cutoff, *tmparr;

    /* Determine the median power */

    tmparr = gen_fvect(n);
    memcpy(tmparr, arr, sizeof(float) * n);
    med = median(tmparr, n);
    vect_free(tmparr);

    /* Throw away powers that are bigger that PRUNELEV * median */

    cutoff = med * PRUNELEV / sqrt((float) numsumpows);
    for (ii = 0; ii < n; ii++) {
        if (arr[ii] > cutoff) {
            arr[ii] = NEWLEV * med;
            ct++;
        }
    }
    return ct;
}


void hpselect(unsigned long m, unsigned long n, float arr[], powindex heap[])
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
