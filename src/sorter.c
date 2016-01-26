#include "presto.h"

/* Function declarations */

int compare_powindex(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int compare_positions(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int compare_fourierprops(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int compare_birds(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int comp_bin_pow(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int comp_bin_nfftbins(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int remove_dupes(position * list, int nlist);
   /*  Removes list values that are 1 unit of search away from a higher */
   /*  power candidate (dr = 0.5 or dz = 2.0).  Returns # removed.      */
int remove_dupes2(fourierprops * list, int nlist);
   /*  Removes list values that are within measurement error away from  */
   /*  a higher power candidate.  Returns # removed.                    */
int remove_dupes_bin(binaryprops * list, int nlist);
   /*  Removes list values that are within 1 Fourier bin of the PSR freq */
   /*  from a higher power candidate. Returns # removed.                 */
int remove_other(fourierprops * list, int nlist, long rlo,
                 long rhi, double locpow, char zapfile,
                 double *lobins, double *hibins, int numzap);
/*  Removes list values whose frequencies fall outside rlo and rhi, */
/*  candidates whose local power levels are below locpow, and       */
/*  candidates close to known birdies.  Returns # removed.          */
int remove_other_bin(binaryprops * list, int nlist);
   /*  Removes list values whose binary parameters are unrealistic.   */
   /*  For example, orbital periods under 200 sec.                    */
   /*  Returns # removed.                                             */


/* int compare_positions(const void *ca, const void *cb) */
int compare_positions(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    position *a, *b;

    a = (position *) ca;
    b = (position *) cb;
    if ((b->pow - a->pow) < 0.0)
        return -1;
    if ((b->pow - a->pow) > 0.0)
        return 1;
    return 0;
}

int compare_floats(const void *a, const void *b)
/* qsort comparison function for floats */
{
    const float *da = (const float *) a;
    const float *db = (const float *) b;

    return (*da > *db) - (*da < *db);
}

int compare_doubles(const void *a, const void *b)
/* qsort comparison function for doubles */
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}

int compare_ints(const void *a, const void *b)
/* qsort comparison function for ints */
{
    const int *da = (const int *) a;
    const int *db = (const int *) b;

    return (*da > *db) - (*da < *db);
}

/* int compare_powindex(const void *ca, const void *cb) */
int compare_powindex(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    powindex *a, *b;

    a = (powindex *) ca;
    b = (powindex *) cb;
    if ((b->pow - a->pow) < 0.0)
        return -1;
    if ((b->pow - a->pow) > 0.0)
        return 1;
    return 0;
}

/* int comp_bin_pow(const void *ca, const void *cb) */
int comp_bin_pow(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    binaryprops *a, *b;

    a = (binaryprops *) ca;
    b = (binaryprops *) cb;
    if ((b->pow - a->pow) < 0.0)
        return -1;
    if ((b->pow - a->pow) > 0.0)
        return 1;
    return 0;
}

/* int comp_bin_nfftbins(const void *ca, const void *cb) */
int comp_bin_nfftbins(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    binaryprops *a, *b;

    a = (binaryprops *) ca;
    b = (binaryprops *) cb;
    if ((b->nfftbins - a->nfftbins) < 0.0)
        return -1;
    if ((b->nfftbins - a->nfftbins) > 0.0)
        return 1;
    return 0;
}

/* int compare_fourierprops(const void *ca, const void *cb) */
int compare_fourierprops(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    fourierprops *a, *b;

    a = (fourierprops *) ca;
    b = (fourierprops *) cb;
    if ((b->pow - a->pow) < 0.0)
        return -1;
    if ((b->pow - a->pow) > 0.0)
        return 1;
    return 0;
}

/* int compare_birds(const void *ca, const void *cb) */
int compare_birds(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
    bird *a, *b;

    a = (bird *) ca;
    b = (bird *) cb;
    if ((b->lobin - a->lobin) < 0.0)
        return 1;
    if ((b->lobin - a->lobin) > 0.0)
        return -1;
    return 0;
}

float percolate(position * list, int nlist, int spot)
/*  Pushes a position structure as far up a sorted list of positions */
/*  as it needs to go to keep the list sorted.  Returns the new low  */
/*  power in the list.                                               */
{
    int ii;
    position tempzz;

    ii = spot;
    while (list[ii - 1].pow < list[ii].pow && ii > 0) {
        SWAP(list[ii - 1], list[ii]);
        ii--;
    }
    return list[nlist - 1].pow;
}


float percolate_bin(binaryprops * list, int nlist)
/*  Pushes a binaryprops structure as far up a sorted list of structs */
/*  as it needs to go to keep the list sorted.  Returns the new low  */
/*  power in the list.                                               */
{
    int ct;
    binaryprops tempzz;

    for (ct = nlist - 2; ct >= 0; ct--) {
        if (list[ct].pow < list[ct + 1].pow) {
            SWAP(list[ct], list[ct + 1]);
        } else {
            break;
        }
    }
    return list[nlist - 1].pow;
}


int remove_dupes(position * list, int nlist)
/*  Removes list values that are 1 unit of search away from a higher */
/*  power candidate (dr = 0.5 or dz = 2.0).  Returns # removed.      */
{
    int i, j, k, ct = 0;
    position tempzz;

    for (i = 0; i < nlist - 1; i++) {
        if (list[i].pow == 0.0)
            break;
        j = i + 1;
        while (j < nlist) {
            if (list[j].pow == 0.0)
                break;
            if ((fabs(list[j].p1 - list[i].p1) < 0.51) &&
                (fabs(list[j].p2 - list[i].p2) < 2.01)) {
                if (j < nlist - 1) {
                    for (k = j; k < nlist - 1; k++) {
                        SWAP(list[k], list[k + 1]);
                    }
                }
                list[nlist - 1].pow = 0.0;
                ct++;
            } else
                j++;
        }
    }
    printf("Removed %d duplicate candidates.\n", ct);
    return ct;
}


int remove_dupes_bin(binaryprops * list, int nlist)
/*  Removes list values that are within 1 Fourier bin of the PSR freq */
/*  from a higher power candidate. Returns # removed.                 */
{
    int i, j, k, ct = 0;
    binaryprops tempzz;

    for (i = 0; i < nlist - 1; i++) {
        if (list[i].pow == 0.0)
            break;
        j = i + 1;
        while (j < nlist) {
            if (list[j].pow == 0.0)
                break;
            if ((fabs(list[j].rdetect - list[i].rdetect) < 0.6) &&
                (fabs(list[j].rpsr - list[i].rpsr) < list[j].nfftbins / 2) &&
                (list[j].nfftbins == list[i].nfftbins)) {
                if (j < nlist - 1) {
                    for (k = j; k < nlist - 1; k++) {
                        SWAP(list[k], list[k + 1]);
                    }
                }
                list[nlist - 1].pow = 0.0;
                ct++;
            } else
                j++;
        }
    }
    return ct;
}


int remove_dupes2(fourierprops * list, int nlist)
/*  Removes list values that are within measurement error away from  */
/*  a higher power candidate.  Returns # removed.                    */
{
    int i, j, k, ct = 0;
    fourierprops tempzz;

    for (i = 0; i < nlist - 1; i++) {
        if (list[i].pow == 0.0)
            break;
        j = i + 1;
        while (j < nlist) {
            if (list[j].pow == 0.0)
                break;
            if ((fabs(list[j].r - list[i].r) < list[i].rerr) &&
                (fabs(list[j].z - list[i].z) < list[i].zerr)) {
                if (j < nlist - 1) {
                    for (k = j; k < nlist - 1; k++) {
                        SWAP(list[k], list[k + 1]);
                    }
                }
                list[nlist - 1].pow = 0.0;
                ct++;
            } else
                j++;
        }
    }
    printf("Removed %d duplicate optimized candidates.\n", ct);
    return ct;
}


int remove_other(fourierprops * list, int nlist, long rlo,
                 long rhi, double locpow, char zapfile,
                 double *lobins, double *hibins, int numzap)
/*  Removes list values whose frequencies fall outside rlo and rhi, */
/*  candidates whose local power levels are below locpow, and       */
/*  candidates close to known birdies.  Returns # removed.          */
{
    int i = 0, j, ct = 0;
    fourierprops tempzz;

    while (i < nlist) {
        if (list[i].pow == 0.0)
            break;
        if (list[i].r < rlo || list[i].r > rhi ||
            list[i].pow < locpow ||
            (zapfile && check_to_zap(list[i].r, lobins, hibins, numzap))) {
            if (i < nlist - 1) {
                for (j = i; j < nlist - 1; j++) {
                    SWAP(list[j], list[j + 1]);
                }
            }
            list[nlist - 1].pow = 0.0;
            ct++;
        } else
            i++;
    }
    printf("Removed %d candidates for various reasons.\n", ct);
    return ct;
}


int remove_other_bin(binaryprops * list, int nlist)
  /*  Removes list values whose binary parameters are unrealistic.   */
  /*  For example, orbital periods under 300 sec.                    */
  /*  Returns # removed.                                             */
{
    int i = 0, j, ct = 0;
    float cutoff = 300.0;
    binaryprops tempzz;

    while (i < nlist) {
        if (list[i].pow == 0.0)
            break;
        if (list[i].pbin < cutoff) {
            if (i < nlist - 1) {
                for (j = i; j < nlist - 1; j++) {
                    SWAP(list[j], list[j + 1]);
                }
            }
            list[nlist - 1].pow = 0.0;
            ct++;
        } else
            i++;
    }
    return ct;
}
