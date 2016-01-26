#include "presto.h"
#include "rfifind.h"

rfi *new_rfi(int numchan, int numint)
/* Create an rfi structure */
{
    int num;
    rfi *newrfi;

    newrfi = (rfi *) malloc(sizeof(rfi));
    newrfi->freq_avg = 0.0;
    newrfi->freq_var = 0.0;
    newrfi->sigma_avg = 0.0;
    newrfi->numobs = 0;
    num = (numint % 8) ? numint / 8 + 1 : numint / 8;
    newrfi->times = (unsigned char *) calloc(1, num);
    num = (numchan % 8) ? numchan / 8 + 1 : numchan / 8;
    newrfi->chans = (unsigned char *) calloc(1, num);
    return newrfi;
}

void write_rfi(FILE * outfile, rfi * outrfi, int numchan, int numint)
/* Write the contents of an rfi structure to a file */
{
    int num;

    chkfwrite(&(outrfi->freq_avg), sizeof(float), 1, outfile);
    chkfwrite(&(outrfi->freq_var), sizeof(float), 1, outfile);
    chkfwrite(&(outrfi->sigma_avg), sizeof(float), 1, outfile);
    chkfwrite(&(outrfi->numobs), sizeof(int), 1, outfile);
    num = (numint % 8) ? numint / 8 + 1 : numint / 8;
    chkfwrite(outrfi->times, 1, num, outfile);
    num = (numchan % 8) ? numchan / 8 + 1 : numchan / 8;
    chkfwrite(outrfi->chans, 1, num, outfile);
}

void read_rfi(FILE * infile, rfi * inrfi, int numchan, int numint)
/* Read the contents of an rfi structure from a file */
{
    int num;

    chkfread(&(inrfi->freq_avg), sizeof(float), 1, infile);
    chkfread(&(inrfi->freq_var), sizeof(float), 1, infile);
    chkfread(&(inrfi->sigma_avg), sizeof(float), 1, infile);
    chkfread(&(inrfi->numobs), sizeof(int), 1, infile);
    num = (numint % 8) ? numint / 8 + 1 : numint / 8;
    chkfread(inrfi->times, 1, num, infile);
    num = (numchan % 8) ? numchan / 8 + 1 : numchan / 8;
    chkfread(inrfi->chans, 1, num, infile);
}

void free_rfi(rfi oldrfi)
/* Free the contents of an rfi structure */
{
    free(oldrfi.times);
    free(oldrfi.chans);
}

rfi *rfi_vector(rfi * rfivect, int numchan, int numint, int oldnum, int newnum)
/* Create or reallocate an rfi_vector */
{
    int ii, numt, numc;

    numt = (numint % 8) ? numint / 8 + 1 : numint / 8;
    numc = (numchan % 8) ? numchan / 8 + 1 : numchan / 8;
    rfivect = (rfi *) realloc(rfivect, sizeof(rfi) * newnum);
    for (ii = oldnum; ii < newnum; ii++) {
        rfivect[ii].freq_avg = 0.0;
        rfivect[ii].freq_var = 0.0;
        rfivect[ii].sigma_avg = 0.0;
        rfivect[ii].numobs = 0;
        rfivect[ii].times = (unsigned char *) calloc(1, numt);
        rfivect[ii].chans = (unsigned char *) calloc(1, numc);
    }
    return rfivect;
}

void free_rfi_vector(rfi * rfivect, int numrfi)
/* Free an rfi vector and its contents */
{
    int ii;

    for (ii = 0; ii < numrfi; ii++)
        free_rfi(rfivect[ii]);
    free(rfivect);
}

void update_rfi(rfi * oldrfi, float freq, float sigma, int channel, int interval)
/* Updates an rfi structure with a new detection */
{
    oldrfi->numobs++;
    if (oldrfi->numobs == 1) {
        oldrfi->freq_avg = freq;
        oldrfi->freq_var = 0.0;
        oldrfi->sigma_avg = sigma;
    } else {
        double dx, an, an1;

        an = (double) (oldrfi->numobs + 1);
        an1 = (double) (oldrfi->numobs);
        dx = (freq - oldrfi->freq_avg) / an;
        oldrfi->freq_var *= (an1 - 1.0);
        oldrfi->freq_var += an * an1 * dx * dx;
        oldrfi->freq_avg += dx;
        oldrfi->freq_var /= an1;
        dx = (sigma - oldrfi->sigma_avg) / an;
        oldrfi->sigma_avg += dx;
    }
    SET_BIT(oldrfi->times, interval);
    SET_BIT(oldrfi->chans, channel);
}

int find_rfi(rfi * rfivect, int numrfi, double freq, double fract_error)
/* Try to find a birdie in an rfi ector.  Compare all */
/* currently known birdies with the new freq.  If it  */
/* finds one with a freq within fractional error, it  */
/* returns the number of the birdie -- otherwise, -1. */
{
    float err;
    int ii = 0;

    for (ii = 0; ii < numrfi; ii++) {
        err = fabs(rfivect[ii].freq_avg - freq) / freq;
        if (err < fract_error)
            return ii;
    }
    return -1;
}

int compare_rfi_freq(const void *ca, const void *cb)
/*  Used as compare function for qsort() */
{
    rfi *a, *b;

    a = (rfi *) ca;
    b = (rfi *) cb;
    if ((b->freq_avg - a->freq_avg) < 0.0)
        return 1;
    if ((b->freq_avg - a->freq_avg) > 0.0)
        return -1;
    return 0;
}

int compare_rfi_sigma(const void *ca, const void *cb)
/*  Used as compare function for qsort() */
{
    rfi *a, *b;

    a = (rfi *) ca;
    b = (rfi *) cb;
    if ((b->sigma_avg - a->sigma_avg) < 0.0)
        return -1;
    if ((b->sigma_avg - a->sigma_avg) > 0.0)
        return 1;
    if ((b->freq_avg - a->freq_avg) < 0.0)
        return -1;
    if ((b->freq_avg - a->freq_avg) > 0.0)
        return 1;
    return 0;
}

int compare_rfi_numobs(const void *ca, const void *cb)
/*  Used as compare function for qsort() */
{
    rfi *a, *b;

    a = (rfi *) ca;
    b = (rfi *) cb;
    if ((b->numobs - a->numobs) < 0.0)
        return -1;
    if ((b->numobs - a->numobs) > 0.0)
        return 1;
    return 0;
}
