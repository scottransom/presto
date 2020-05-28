#include "presto.h"
#include "mask.h"

extern int compare_ints(const void *a, const void *b);
extern int compare_floats(const void *a, const void *b);

static int find_num(int num, int *arr, int arrlen);
static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, int *merged, int lenout);

void fill_mask(double timesigma, double freqsigma, double mjd,
               double dtint, double lofreq, double dfreq,
               int numchan, int numint, int ptsperint,
               int num_zap_chans, int *zap_chans, int num_zap_ints,
               int *zap_ints, unsigned char **bytemask, mask * obsmask)
/* Fill a mask structure with the appropriate values */
{
    int ii, jj, count;

    obsmask->timesigma = timesigma;
    obsmask->freqsigma = freqsigma;
    obsmask->mjd = mjd;
    obsmask->dtint = dtint;
    obsmask->lofreq = lofreq;
    obsmask->dfreq = dfreq;
    obsmask->numchan = numchan;
    obsmask->numint = numint;
    obsmask->ptsperint = ptsperint;
    obsmask->num_zap_chans = num_zap_chans;
    if (obsmask->num_zap_chans) {
        obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
        for (ii = 0; ii < obsmask->num_zap_chans; ii++)
            obsmask->zap_chans[ii] = zap_chans[ii];
    }
    obsmask->num_zap_ints = num_zap_ints;
    if (obsmask->num_zap_ints) {
        obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
        for (ii = 0; ii < obsmask->num_zap_ints; ii++)
            obsmask->zap_ints[ii] = zap_ints[ii];
    }
    obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
    obsmask->chans = (int **) malloc(obsmask->numint * sizeof(int *));
    for (ii = 0; ii < obsmask->numint; ii++) {
        count = 0;
        /* Count the bad channels first */
        for (jj = 0; jj < obsmask->numchan; jj++)
            if ((bytemask[ii][jj] & BADDATA) | (bytemask[ii][jj] & USERZAP))
                count++;
        obsmask->num_chans_per_int[ii] = count;
        if (count) {
            /* Now determine which channels */
            count = 0;
            obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
            for (jj = 0; jj < obsmask->numchan; jj++) {
                if ((bytemask[ii][jj] & BADDATA) | (bytemask[ii][jj] & USERZAP))
                    obsmask->chans[ii][count++] = jj;
            }
        }
    }
}


void set_oldmask_bits(mask * oldmask, unsigned char **bytemask)
/* Sets the oldmask bit in the appropriate bytes in bytemask */
{
    int ii, jj;

    for (ii = 0; ii < oldmask->numint; ii++)
        for (jj = 0; jj < oldmask->num_chans_per_int[ii]; jj++)
            bytemask[ii][oldmask->chans[ii][jj]] |= OLDMASK;
}


void unset_oldmask_bits(mask * oldmask, unsigned char **bytemask)
/* Unsets the oldmask bits in bytemask */
{
    int ii, jj;

    for (ii = 0; ii < oldmask->numint; ii++)
        for (jj = 0; jj < oldmask->numchan; jj++)
            bytemask[ii][jj] &= ~OLDMASK;
}


void free_mask(mask obsmask)
/* Free the contents of an mask structure */
{
    int ii;

    for (ii = 0; ii < obsmask.numint; ii++) {
        if (obsmask.num_chans_per_int[ii] > 0 &&
            obsmask.num_chans_per_int[ii] <= obsmask.numchan)
            vect_free(obsmask.chans[ii]);
    }
    free(obsmask.chans);
    vect_free(obsmask.num_chans_per_int);
    if (obsmask.num_zap_chans)
        vect_free(obsmask.zap_chans);
    if (obsmask.num_zap_ints)
        vect_free(obsmask.zap_ints);
}


void read_mask(char *maskfilenm, mask * obsmask)
/* Read the contents of a mask structure from a file */
{
    FILE *infile;
    int ii;

    infile = chkfopen(maskfilenm, "rb");
    chkfread(&(obsmask->timesigma), sizeof(double), 1, infile);
    chkfread(&(obsmask->freqsigma), sizeof(double), 1, infile);
    chkfread(&(obsmask->mjd), sizeof(double), 1, infile);
    chkfread(&(obsmask->dtint), sizeof(double), 1, infile);
    chkfread(&(obsmask->lofreq), sizeof(double), 1, infile);
    chkfread(&(obsmask->dfreq), sizeof(double), 1, infile);
    chkfread(&(obsmask->numchan), sizeof(int), 1, infile);
    chkfread(&(obsmask->numint), sizeof(int), 1, infile);
    chkfread(&(obsmask->ptsperint), sizeof(int), 1, infile);
    chkfread(&(obsmask->num_zap_chans), sizeof(int), 1, infile);
    if (obsmask->num_zap_chans) {
        obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
        chkfread(obsmask->zap_chans, sizeof(int), obsmask->num_zap_chans, infile);
    }
    chkfread(&(obsmask->num_zap_ints), sizeof(int), 1, infile);
    if (obsmask->num_zap_ints) {
        obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
        chkfread(obsmask->zap_ints, sizeof(int), obsmask->num_zap_ints, infile);
    }
    obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
    chkfread(obsmask->num_chans_per_int, sizeof(int), obsmask->numint, infile);
    obsmask->chans = (int **) malloc(obsmask->numint * sizeof(int *));
    for (ii = 0; ii < obsmask->numint; ii++) {
        if (obsmask->num_chans_per_int[ii] > 0 &&
            obsmask->num_chans_per_int[ii] < obsmask->numchan) {
            obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
            chkfread(obsmask->chans[ii], sizeof(int),
                     obsmask->num_chans_per_int[ii], infile);
        } else if (obsmask->num_chans_per_int[ii] == obsmask->numchan) {
            int jj;
            obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
            for (jj = 0; jj < obsmask->numchan; jj++)
                obsmask->chans[ii][jj] = jj;
        }
    }
    fclose(infile);
}


void calc_avgmedstd(float *arr, int numarr, float fraction,
                    int step, float *avg, float *med, float *std)
/* Calculates the median and middle-'fraction' std deviation  */
/* and average of the array 'arr'.  Values are returned in    */
/* 'avg', 'med' and 'std'.  The array is not modified.        */
{
    int ii, jj, len, start;
    float *tmparr;
    double davg, dstd;

    len = (int) (numarr * fraction + 0.5);
    if (len > numarr || len < 0) {
        printf("fraction (%g) out-of-bounds in calc_avgmedstd()\n", fraction);
        exit(1);
    }
    start = (numarr - len) / 2;
    tmparr = gen_fvect(numarr);
    for (ii = 0, jj = 0; ii < numarr; ii++, jj += step)
        tmparr[ii] = arr[jj];
    qsort(tmparr, numarr, sizeof(float), compare_floats);
    avg_var(tmparr + start, len, &davg, &dstd);
    *avg = (float) davg;
    *med = tmparr[numarr / 2];
    *std = sqrt(dstd);
    vect_free(tmparr);
}


int determine_padvals(char *maskfilenm, mask * obsmask, float *padvals)
// Determine reasonable padding values from the rfifind produced
// *.stats file if it is available.  The pre-allocated vector (of
// length numchan) is in padvals.  Return a '1' if the routine used
// the stats file, return 0 if the padding was set to zeros.
{
    FILE *statsfile;
    int ii, numchan, numint, ptsperint, lobin, numbetween;
    float **dataavg, tmp1, tmp2;
    char *statsfilenm, *root, *suffix;

    if (split_root_suffix(maskfilenm, &root, &suffix) == 0) {
        printf("\nThe mask filename (%s) must have a suffix!\n\n", maskfilenm);
        exit(1);
    } else {
        /* Determine the stats file name */
        statsfilenm = calloc(strlen(maskfilenm) + 2, sizeof(char));
        sprintf(statsfilenm, "%s.stats", root);
        free(root);
        free(suffix);
        /* Check to see if the file exists */
        printf("Attempting to read the data statistics from '%s'...\n", statsfilenm);
        statsfile = chkfopen(statsfilenm, "rb");
        free(statsfilenm);
        if (statsfile) {        /* Read the stats */
            chkfread(&numchan, sizeof(int), 1, statsfile);
            chkfread(&numint, sizeof(int), 1, statsfile);
            chkfread(&ptsperint, sizeof(int), 1, statsfile);
            chkfread(&lobin, sizeof(int), 1, statsfile);
            chkfread(&numbetween, sizeof(int), 1, statsfile);
            dataavg = gen_fmatrix(numint, numchan);
            /* These are the powers */
            chkfread(dataavg[0], sizeof(float), numchan * numint, statsfile);
            /* These are the averages */
            chkfread(dataavg[0], sizeof(float), numchan * numint, statsfile);
            /* Set the padding values equal to the mid-80% channel averages */
            for (ii = 0; ii < numchan; ii++)
                calc_avgmedstd(dataavg[0] + ii, numint, 0.8, numchan,
                               padvals + ii, &tmp1, &tmp2);
            printf
                ("...succeded.  Set the padding values equal to the mid-80%% channel averages.\n");
            vect_free(dataavg[0]);
            vect_free(dataavg);
            fclose(statsfile);
            return 1;
        } else {
            /* This is a temporary solution */
            for (ii = 0; ii < obsmask->numchan; ii++)
                padvals[ii] = 0.0;
            printf("...failed.\n  Set the padding values to 0.\n");
            return 0;
        }
    }
}


void write_mask(char *maskfilenm, mask * obsmask)
/* Write the contents of an mask structure to a file */
{
    FILE *outfile;
    int ii;

    outfile = chkfopen(maskfilenm, "wb");
    chkfwrite(&(obsmask->timesigma), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->freqsigma), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->mjd), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->dtint), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->lofreq), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->dfreq), sizeof(double), 1, outfile);
    chkfwrite(&(obsmask->numchan), sizeof(int), 1, outfile);
    chkfwrite(&(obsmask->numint), sizeof(int), 1, outfile);
    chkfwrite(&(obsmask->ptsperint), sizeof(int), 1, outfile);
    chkfwrite(&(obsmask->num_zap_chans), sizeof(int), 1, outfile);
    if (obsmask->num_zap_chans)
        chkfwrite(obsmask->zap_chans, sizeof(int), obsmask->num_zap_chans, outfile);
    chkfwrite(&(obsmask->num_zap_ints), sizeof(int), 1, outfile);
    if (obsmask->num_zap_ints)
        chkfwrite(obsmask->zap_ints, sizeof(int), obsmask->num_zap_ints, outfile);
    chkfwrite(obsmask->num_chans_per_int, sizeof(int), obsmask->numint, outfile);
    for (ii = 0; ii < obsmask->numint; ii++) {
        if (obsmask->num_chans_per_int[ii] > 0 &&
            obsmask->num_chans_per_int[ii] < obsmask->numchan) {
            chkfwrite(obsmask->chans[ii], sizeof(int),
                      obsmask->num_chans_per_int[ii], outfile);
        }
    }
    fclose(outfile);
}


int check_mask(double starttime, double duration, mask * obsmask, int *maskchans)
/* Return value is the number of channels to mask.  The */
/* channel numbers are placed in maskchans (which must  */
/* have a length of numchan).  If -1 is returned, all   */
/* channels should be masked.                           */
{
    int loint, hiint;
    double endtime;
    static int old_loint = -1, old_hiint = -1, old_numchan = 0;

    /*
       static int firsttime = 1;
       if (firsttime){
       int ii;
       printf("\n\n numzapints = %d\n : ", obsmask->num_zap_ints);
       for (ii=0; ii<obsmask->num_zap_ints; ii++)
       printf("%d ", obsmask->zap_ints[ii]);
       printf("\n\n numzapchans = %d\n : ", obsmask->num_zap_chans);
       for (ii=0; ii<obsmask->num_zap_chans; ii++)
       printf("%d ", obsmask->zap_chans[ii]);
       printf("\n\n");
       firsttime = 0;
       }
     */

    endtime = starttime + duration;
    loint = (int) (starttime / obsmask->dtint);
    hiint = (int) (endtime / obsmask->dtint);

    /* Mask the same channels as for the last call */
    if (loint == old_loint && hiint == old_hiint)
        return old_numchan;

    /* Make sure that we aren't past the last interval */
    if (loint >= obsmask->numint)
        loint = obsmask->numint - 1;
    if (hiint >= obsmask->numint)
        hiint = loint;

    /* Determine new channels to mask */
    if (loint == hiint) {
        old_loint = old_hiint = loint;
        /* Check to see if this is an interval where we zap all the channels */
        if (obsmask->num_zap_ints) {
            if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)) {
                old_numchan = -1;
                return old_numchan;
            }
        }
        /* Merge the overall channels to zap with the local channels to zap */
        old_numchan = merge_no_dupes(obsmask->zap_chans,
                                     obsmask->num_zap_chans,
                                     obsmask->chans[loint],
                                     obsmask->num_chans_per_int[loint],
                                     maskchans, obsmask->numchan);
    } else {                    /* We are straddling a rfifind interval boundary */
        int *tmpchans;

        old_loint = loint;
        old_hiint = hiint;
        /* Check to see if this is an interval where we zap all the channels */
        if (obsmask->num_zap_ints) {
            if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)) {
                old_numchan = -1;
                return old_numchan;
            }
            if (find_num(hiint, obsmask->zap_ints, obsmask->num_zap_ints)) {
                old_numchan = -1;
                return old_numchan;
            }
        }
        /* Merge the overall channels to zap with the loint channels to zap */
        if (obsmask->num_zap_chans) {
            tmpchans = gen_ivect(obsmask->numchan);
            old_numchan = merge_no_dupes(obsmask->zap_chans,
                                         obsmask->num_zap_chans,
                                         obsmask->chans[loint],
                                         obsmask->num_chans_per_int[loint],
                                         tmpchans, obsmask->numchan);
        } else {
            tmpchans = obsmask->zap_chans;
            old_numchan = obsmask->num_zap_chans;
        }
        /* Merge the loint+overall channels to zap with the hiint channels to zap */
        old_numchan = merge_no_dupes(tmpchans,
                                     old_numchan,
                                     obsmask->chans[hiint],
                                     obsmask->num_chans_per_int[hiint],
                                     maskchans, obsmask->numchan);
        if (obsmask->num_zap_chans)
            vect_free(tmpchans);
    }
    return old_numchan;
}


static int find_num(int num, int *arr, int arrlen)
{
    int ii;

    /* Note:  I should make sure the array is sorted and do a binary search */
    for (ii = 0; ii < arrlen; ii++)
        if (arr[ii] == num)
            return 1;
    return 0;
}


static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, int *merged, int lenout)
{
    int ptr1 = 0, ptr2 = 0, count = 0;

    while (1) {
        if (ptr1 == len1) {
            while ((ptr2 < len2) && (merged[count-1] < lenout - 1)) {
                if (arr2[ptr2] > merged[count-1])
                    merged[count++] = arr2[ptr2++];
                else ptr2++;
            }
            break;
        } else if (ptr2 == len2) {
            while ((ptr1 < len1) && (merged[count-1] < lenout - 1)) {
                if (arr1[ptr1] > merged[count-1])
                    merged[count++] = arr1[ptr1++];
                else ptr1++;
            }
            break;
        }
        if (arr1[ptr1] < arr2[ptr2])
            merged[count++] = arr1[ptr1++];
        else if (arr1[ptr1] > arr2[ptr2])
            merged[count++] = arr2[ptr2++];
        else {
            merged[count++] = arr1[ptr1];
            ptr1++;
            ptr2++;
        }
    }
    return count;
}
