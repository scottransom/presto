#include "presto.h"

/* Number of bins on each side of a freq to use for interpolation */
#define INTERPBINS 5

/* Minimum sigma detection to return.            */
/* (Including corrections for independant freqs) */
#define MINRETURNSIG 1.5

/* Routines defined at the bottom */

float percolate_rawbincands(rawbincand * cands, int numcands);
float percolate_fftcands(fftcand * cands, int numcands);
void print_rawbincand(rawbincand cand);
int comp_rawbin_to_cand(rawbincand * cand, infodata * idata, char *output, int full);

fftcand *search_fft(fcomplex * fft, int numfft, int lobin, int hibin,
                    int numharmsum, int numbetween,
                    presto_interptype interptype,
                    float norm, float sigmacutoff, int *numcands,
                    float *powavg, float *powvar, float *powmax)
/* This routine searches a short FFT of 'numfft' complex freqs      */
/* and returns a candidate vector of fftcand structures containing  */
/* information about the best candidates found.                     */
/* The routine uses either interbinning or interpolation as well    */
/* as harmonic summing during the search.                           */
/* The number of candidates returned is either 'numcands' if != 0,  */
/* or is determined automatically by 'sigmacutoff' -- which         */
/* takes into account the number of bins searched.                  */
/* The returned vector is sorted in order of decreasing power.      */
/* Arguments:                                                       */
/*   'fft' is the FFT to search (complex valued)                    */
/*   'numfft' is the number of complex points in 'fft'              */
/*   'lobin' is the lowest Fourier freq to search                   */
/*   'hibin' is the highest Fourier freq to search                  */
/*   'numharmsum' the number of harmonics to sum during the search  */
/*   'numbetween' the points to interpolate per bin                 */
/*   'interptype' is either INTERBIN or INTERPOLATE.                */
/*      INTERBIN = (interbinning) is fast but less sensitive.       */
/*         NOTE:  INTERBINNING is conducted by this routine!        */
/*      INTERPOLATE = (Fourier interpolation) is slower but more    */
/*        sensitive.                                                */
/*         NOTE:  The interpolation is assumed to ALREADY have been */
/*                completed by the calling function!  The easiest   */
/*                way is by zero-padding to 2*numfft and FFTing.    */
/*                If you use this method, make sure numfft is the   */
/*                original length rather than the interpolated      */
/*                length and also make sure numbetween is correct.  */
/*   'norm' is the normalization constant to multiply each power by */
/*   'sigmacutoff' if the number of candidates will be determined   */
/*      automatically, is the minimum Gaussian significance of      */
/*      candidates to keep -- taking into account the number of     */
/*      bins searched                                               */
/*   'numcands' if !0, is the number of candates to return.         */
/*      if 0, is a return value giving the number of candidates.    */
/*   'powavg' is a return value giving the average power level      */
/*   'powvar' is a return value giving the power level variance     */
/*   'powmax' is a return value giving the maximum power            */
{
    int ii, jj, offset, numtosearch, dynamic = 0;
    int numspread = 0, nc = 0, startnc = 10;
    float powargr, powargi, *fullpows = NULL, *sumpows, ftmp;
    double twobypi, minpow = 0.0, tmpminsig = 0.0, dr, davg, dvar;
    fftcand *cands, newcand;
    fcomplex *spread;

    /* Override the value of numbetween if interbinning */

    if (interptype == INTERBIN)
        numbetween = 2;
    norm = 1.0 / norm;
    *powmax = 0.0;

    /* Decide if we will manage the number of candidates */

    if (*numcands > 0)
        startnc = *numcands;
    else {
        dynamic = 1;
        minpow = power_for_sigma(sigmacutoff, 1, hibin - lobin);
    }
    cands = (fftcand *) malloc(startnc * sizeof(fftcand));
    for (ii = 0; ii < startnc; ii++)
        cands[ii].sig = 0.0;

    /* Prep some other values we will need */

    dr = 1.0 / (double) numbetween;
    twobypi = 2.0 / PI;
    numtosearch = numfft * numbetween;

    /* Spread and interpolate the fft */

    numspread = numfft * numbetween + 1;
    if (interptype == INTERPOLATE) {    /* INTERPOLATE */
        spread = fft;
    } else {                    /* INTERBIN */
        spread = gen_cvect(numspread);
        spread_with_pad(fft, numfft, spread, numspread, numbetween, 0);
        for (ii = 1; ii < numtosearch; ii += 2) {
            spread[ii].r = twobypi * (spread[ii - 1].r - spread[ii + 1].r);
            spread[ii].i = twobypi * (spread[ii - 1].i - spread[ii + 1].i);
        }
    }
    spread[0].r = spread[numtosearch].r = 1.0;
    spread[0].i = spread[numtosearch].i = 0.0;

    /* First generate the original powers in order to         */
    /* calculate the statistics.  Yes, this is inefficient... */

    fullpows = gen_fvect(numtosearch);
    for (ii = lobin, jj = 0; ii < hibin; ii++, jj++) {
        ftmp = POWER(fft[ii].r, fft[ii].i) * norm;
        fullpows[jj] = ftmp;
        if (ftmp > *powmax)
            *powmax = ftmp;
    }
    avg_var(fullpows, hibin - lobin, &davg, &dvar);
    *powavg = davg;
    *powvar = dvar;
    fullpows[0] = 1.0;
    for (ii = 1; ii < numtosearch; ii++)
        fullpows[ii] = POWER(spread[ii].r, spread[ii].i) * norm;
    if (interptype == INTERBIN)
        vect_free(spread);

    /* Search the raw powers */

    for (ii = lobin * numbetween; ii < hibin * numbetween; ii++) {
        if (fullpows[ii] > minpow) {
            newcand.r = dr * (double) ii;
            newcand.p = fullpows[ii];
            newcand.sig = candidate_sigma(fullpows[ii], 1, hibin - lobin);
            newcand.nsum = 1;
            cands[startnc - 1] = newcand;
            tmpminsig = percolate_fftcands(cands, startnc);
            if (dynamic) {
                nc++;
                if (nc == startnc) {
                    startnc *= 2;
                    cands = (fftcand *) realloc(cands, startnc * sizeof(fftcand));
                    for (jj = nc; jj < startnc; jj++)
                        cands[jj].sig = 0.0;
                }
            } else {
                minpow = cands[startnc - 1].p;
                if (nc < startnc)
                    nc++;
            }
        }
    }

    /* If needed, sum and search the harmonics */

    if (numharmsum > 1) {
        sumpows = gen_fvect(numtosearch);
        memcpy(sumpows, fullpows, sizeof(float) * numtosearch);
        for (ii = 2; ii <= numharmsum; ii++) {
            offset = ii / 2;
            if (dynamic)
                minpow = power_for_sigma(sigmacutoff, ii, hibin - lobin);
            else
                minpow = power_for_sigma(tmpminsig, ii, hibin - lobin);
            for (jj = lobin * numbetween; jj < numtosearch; jj++) {
                sumpows[jj] += fullpows[(jj + offset) / ii];
                if (sumpows[jj] > minpow) {
                    newcand.r = dr * (double) jj;
                    newcand.p = sumpows[jj];
                    newcand.sig = candidate_sigma(sumpows[jj], ii, hibin - lobin);
                    newcand.nsum = ii;
                    cands[startnc - 1] = newcand;
                    tmpminsig = percolate_fftcands(cands, startnc);
                    if (dynamic) {
                        nc++;
                        if (nc == startnc) {
                            startnc *= 2;
                            cands =
                                (fftcand *) realloc(cands,
                                                    startnc * sizeof(fftcand));
                            for (jj = nc; jj < startnc; jj++)
                                cands[jj].sig = 0.0;
                        }
                    } else {
                        minpow = power_for_sigma(tmpminsig, ii, hibin - lobin);
                        if (nc < startnc)
                            nc++;
                    }
                }
            }
        }
        vect_free(sumpows);
    }
    vect_free(fullpows);

    /* Chop off the unused parts of the dynamic array */

    if (dynamic)
        cands = (fftcand *) realloc(cands, nc * sizeof(fftcand));
    *numcands = nc;
    return cands;
}


void search_minifft(fcomplex * minifft, int numminifft,
                    double min_orb_p, double max_orb_p,
                    rawbincand * cands, int numcands, int numharmsum,
                    int numbetween, double numfullfft, double timefullfft,
                    double lorfullfft, presto_interptype interptype,
                    presto_checkaliased checkaliased)
  /* This routine searches a short FFT (usually produced using the   */
  /* MiniFFT binary search method) and returns a candidte vector     */
  /* containing information about the best binary candidates found.  */
  /* The routine uses either interbinning or interpolation as well   */
  /* as harmonic summing during the search.                          */
  /* Arguments:                                                      */
  /*   'minifft' is the FFT to search (complex valued)               */
  /*   'numminifft' is the number of complex points in 'minifft'     */
  /*   'min_orb_p' is the minimum orbital period (s) to search       */
  /*   'max_orb_p' is the maximum orbital period (s) to search       */
  /*   'cands' is a pre-allocated vector of rawbincand type in which */
  /*      the sorted (in decreasing sigma) candidates are returned   */
  /*   'numcands' is the length of the 'cands' vector                */
  /*   'numharmsum' the number of harmonics to sum during the search */
  /*   'numbetween' the points to interpolate per bin                */
  /*   'numfullfft' the number of points in the original long FFT    */
  /*   'timefullfft' the duration of the original time series (s)    */
  /*   'lorfullfft' the 1st bin of the long FFT that was miniFFT'd   */
  /*   'interptype' is either INTERBIN or INTERPOLATE.               */
/*      INTERBIN = (interbinning) is fast but less sensitive.        */
/*         NOTE:  INTERBINNING is conducted by this routine!         */
/*      INTERPOLATE = (Fourier interpolation) is slower but more     */
/*        sensitive.                                                 */
/*         NOTE:  The interpolation is assumed to ALREADY have been  */
/*                completed by the calling function!  The easiest    */
/*                way is by zero-padding to 2*numminifft and FFTing. */
/*                If you use this method, make sure numminifft is the*/
/*                original length rather than the interpolated       */
/*                length and also make sure numbetween is correct.   */
  /*   'checkaliased' is either CHECK_ALIASED or NO_CHECK_ALIASED.   */
  /*      NO_CHECK_ALIASED = harmonic summing does not include       */
  /*        aliased freqs making it faster but less sensitive.       */
  /*      CHECK_ALIASED = harmonic summing includes aliased freqs    */
  /*        making it slower but more sensitive.                     */
{
    int ii, jj, fftlen, offset, numtosearch = 0, lobin, hibin, numspread = 0;
    float powargr, powargi, *fullpows = NULL, *sumpows;
    double twobypi, minpow, minsig, dr, numindep;
    fcomplex *spread;

    /* Override the value of numbetween if interbinning */

    if (interptype == INTERBIN)
        numbetween = 2;

    /* Prep some other values we will need */

    dr = 1.0 / (double) numbetween;
    twobypi = 2.0 / PI;
    fftlen = numminifft * numbetween;
    for (ii = 0; ii < numcands; ii++) {
        cands[ii].mini_sigma = 0.0;
        cands[ii].mini_power = 0.0;
    }
    lobin = ceil(2 * numminifft * min_orb_p / timefullfft);
    if (lobin <= 0)
        lobin = 1;
    hibin = floor(2 * numminifft * max_orb_p / timefullfft);
    if (hibin >= 2 * numminifft)
        hibin = 2 * numminifft - 1;
    lobin *= numbetween;
    hibin *= numbetween;

    /* Spread and interpolate the fft */

    numtosearch = (checkaliased == CHECK_ALIASED) ? 2 * fftlen : fftlen;
    numspread = numminifft * numbetween + 1;
    if (interptype == INTERPOLATE) {    /* INTERPOLATE */
        spread = minifft;
    } else {                    /* INTERBIN */
        spread = gen_cvect(numspread);
        spread_with_pad(minifft, numminifft, spread, numspread, numbetween, 0);
        for (ii = 1; ii < fftlen; ii += 2) {
            spread[ii].r = twobypi * (spread[ii - 1].r - spread[ii + 1].r);
            spread[ii].i = twobypi * (spread[ii - 1].i - spread[ii + 1].i);
        }
    }
    spread[0].r = spread[fftlen].r = 1.0;
    spread[0].i = spread[fftlen].i = 0.0;

    fullpows = gen_fvect(numtosearch);
    fullpows[0] = 1.0;
    if (checkaliased == CHECK_ALIASED)
        fullpows[fftlen] = 1.0; /* used to be nyquist^2 */

    /* The following wraps the data around the Nyquist freq such that */
    /* we consider aliased frequencies as well (If CHECK_ALIASED).    */

    if (checkaliased == CHECK_ALIASED)
        for (ii = 1, jj = numtosearch - 1; ii < fftlen; ii++, jj--)
            fullpows[ii] = fullpows[jj] = POWER(spread[ii].r, spread[ii].i);
    else
        for (ii = 1; ii < numtosearch; ii++)
            fullpows[ii] = POWER(spread[ii].r, spread[ii].i);
    if (interptype == INTERBIN)
        vect_free(spread);

    /* Search the raw powers */

    numindep = hibin - lobin + 1.0;
    minpow = power_for_sigma(MINRETURNSIG, 1, numindep);
    for (ii = lobin; ii < hibin; ii++) {
        if (fullpows[ii] > minpow) {
            cands[numcands - 1].mini_r = dr * (double) ii;
            cands[numcands - 1].mini_power = fullpows[ii];
            cands[numcands - 1].mini_numsum = 1.0;
            cands[numcands - 1].mini_sigma =
                candidate_sigma(fullpows[ii], 1, numindep);
            minsig = percolate_rawbincands(cands, numcands);
            if (cands[numcands - 1].mini_power > minpow)
                minpow = cands[numcands - 1].mini_power;
        }
    }

    /* If needed, sum and search the harmonics */

    if (numharmsum > 1) {
        sumpows = gen_fvect(numtosearch);
        memcpy(sumpows, fullpows, sizeof(float) * numtosearch);
        for (ii = 2; ii <= numharmsum; ii++) {
            offset = ii / 2;
            numindep = (hibin - lobin + 1.0) / (double) ii;
            if (cands[numcands - 1].mini_sigma < MINRETURNSIG)
                minsig = MINRETURNSIG;
            else
                minsig = cands[numcands - 1].mini_sigma;
            minpow = power_for_sigma(minsig, ii, numindep);
            for (jj = lobin * ii; jj < hibin; jj++) {
                sumpows[jj] += fullpows[(jj + offset) / ii];
                if (sumpows[jj] > minpow) {
                    cands[numcands - 1].mini_r = (dr * (double) jj) / ii;
                    cands[numcands - 1].mini_power = sumpows[jj];
                    cands[numcands - 1].mini_numsum = (double) ii;
                    cands[numcands - 1].mini_sigma =
                        candidate_sigma(sumpows[jj], ii, numindep);
                    minsig = percolate_rawbincands(cands, numcands);
                    if (minsig > MINRETURNSIG)
                        minpow = power_for_sigma(minsig, ii, numindep);
                }
            }
        }
        vect_free(sumpows);
    }
    vect_free(fullpows);

    /* Add the rest of the rawbincand data to the candidate array */

    for (ii = 0; ii < numcands; ii++) {
        cands[ii].full_N = numfullfft;
        cands[ii].full_T = timefullfft;
        cands[ii].full_lo_r = lorfullfft;
        cands[ii].mini_N = 2 * numminifft;      /* # of real points */
        cands[ii].psr_p = timefullfft / (lorfullfft + numminifft);
        cands[ii].orb_p = timefullfft * cands[ii].mini_r / cands[ii].mini_N;
    }
}


void print_rawbincand(rawbincand cand)
{
    printf("  Sigma       =  %-7.3f\n", cand.mini_sigma);
    printf("  Orbit p     =  %-8.2f\n", cand.orb_p);
    if (cand.psr_p < 0.001)
        printf("  Pulsar p    =  %-12.5e\n", cand.psr_p);
    else
        printf("  Pulsar p    =  %-12.9f\n", cand.psr_p);
    printf("  rlo (full)  =  %-10.0f\n", cand.full_lo_r);
    printf("  N (mini)    =  %-6.0f\n", cand.mini_N);
    printf("  r (detect)  =  %-9.3f\n", cand.mini_r);
    printf("  Power       =  %-8.3f\n", cand.mini_power);
    printf("  Numsum      =  %-2.0f\n", cand.mini_numsum);
    printf("  N (full)    =  %-10.0f\n", cand.full_N);
    printf("  T (full)    =  %-13.6f\n\n", cand.full_T);
}

float percolate_fftcands(fftcand * cands, int numcands)
  /*  Pushes a fftcand candidate as far up the array of   */
  /*  candidates as it shoud go to keep the array sorted  */
  /*  in indecreasing sigmas.  Returns the new lowest     */
  /*  sigma in the array.                                 */
{
    int ii;
    fftcand tempzz;

    for (ii = numcands - 2; ii >= 0; ii--) {
        if (cands[ii].sig < cands[ii + 1].sig) {
            SWAP(cands[ii], cands[ii + 1]);
        } else {
            break;
        }
    }
    return cands[numcands - 1].sig;
}


float percolate_rawbincands(rawbincand * cands, int numcands)
  /*  Pushes a rawbincand candidate as far up the array of   */
  /*  candidates as it shoud go to keep the array sorted in  */
  /*  indecreasing significance.  Returns the new lowest     */
  /*  sigma in the array.                                    */
{
    int ii;
    rawbincand tempzz;

    for (ii = numcands - 2; ii >= 0; ii--) {
        if (cands[ii].mini_sigma < cands[ii + 1].mini_sigma) {
            SWAP(cands[ii], cands[ii + 1]);
        } else {
            break;
        }
    }
    return cands[numcands - 1].mini_sigma;
}


int not_already_there_rawbin(rawbincand newcand, rawbincand * list, int nlist)
{
    int ii;

    /* Loop through the candidates already in the list */

    for (ii = 0; ii < nlist; ii++) {
        if (list[ii].mini_sigma == 0.0)
            break;

        /* Do not add the candidate to the list if it is a lower power */
        /* version of an already listed candidate.                     */

        if (list[ii].mini_N == newcand.mini_N) {
            if (fabs(list[ii].mini_r - newcand.mini_r) < 0.6) {
                if (list[ii].mini_sigma > newcand.mini_sigma) {
                    return 0;
                }
            }
        }
    }
    return 1;
}


void compare_rawbin_cands(rawbincand * list, int nlist, char *notes)
{
    double perr;
    int ii, jj, kk, ll;
    char tmp[30];

    /* Loop through the candidates (reference cands) */

    for (ii = 0; ii < nlist; ii++) {

        /* Loop through the candidates (referenced cands) */

        for (jj = 0; jj < nlist; jj++) {
            if (ii == jj)
                continue;
            perr = 0.5 * list[jj].full_T / list[jj].mini_N;

            /* Loop through the possible PSR period harmonics */

            for (kk = 1; kk < 41; kk++) {

                /* Check if the PSR Fourier freqs are close enough */

                if (fabs(list[ii].full_lo_r - list[jj].full_lo_r / kk) <
                    list[ii].mini_N) {

                    /* Loop through the possible binary period harmonics */

                    for (ll = 1; ll < 10; ll++) {

                        /* Check if the binary Fourier freqs are close enough */

                        if (fabs(list[ii].orb_p - list[jj].orb_p / ll) < perr) {

                            /* Check if the note has already been written */

                            sprintf(tmp, "%.18s", notes + jj * 18);
                            if (!strcmp("                  ", tmp)) {

                                /* Write the note */

                                if (ll == 1 && kk == 1)
                                    sprintf(notes + jj * 18, "Same as #%d?", ii + 1);
                                else
                                    sprintf(notes + jj * 18, "MH=%d H=%d of #%d", ll,
                                            kk, ii + 1);

                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}


void file_rawbin_candidates(rawbincand * cand, char *notes,
                            int numcands, int numharm, char name[])
/* Outputs a .ps file describing all the binary candidates from a    */
/*   binary search. */
{
    FILE *fname;
    int i, j, k = 0;
    int nlines = 87, pages, extralines, linestoprint;
    char *filenm, command[200];
    double orbperr, psrperr;

    filenm = (char *) calloc(strlen(name) + 10, sizeof(char));
    sprintf(filenm, "%s_bin%d", name, numharm);
    fname = chkfopen(filenm, "w");

    if (numcands <= 0) {
        printf(" Must have at least 1 candidate in ");
        printf("file_bin_candidates().\n\n");
        exit(1);
    }
    pages = numcands / nlines + 1;
    extralines = numcands % nlines;

    for (i = 1; i <= pages; i++) {

        /*                       1         2         3         4         5         6         7         8         9         0         1    */
        /*              123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234 */
        fprintf(fname,
                "#               P_orbit +/- Error   P_pulsar +/- Error   FullFFT   MiniFFT   MiniFFT  Num   Sum                   \n");
        fprintf(fname,
                "# Cand  Sigma         (sec)                (sec)         Low Bin   Length      Bin    Sum  Power  Notes           \n");
        fprintf(fname,
                "#------------------------------------------------------------------------------------------------------------------\n");

        if (i == pages) {
            linestoprint = extralines;
        } else {
            linestoprint = nlines;
        }

        for (j = 0; j < linestoprint; j++, k++) {

            /* Calculate the approximate error in our value of orbital period */
            orbperr = 0.5 * cand[k].full_T / cand[k].mini_N;

            /* Calculate the approximate error in our value of spin period */

            if (cand[k].full_lo_r == 0.0)
                psrperr = cand[k].psr_p;
            else
                psrperr = fabs(cand[k].full_T / (cand[k].full_lo_r +
                                                 0.5 * cand[k].mini_N) -
                               cand[k].full_T / cand[k].full_lo_r);

            /*  Now output it... */

            fprintf(fname, "%-5d %7.3f  ", k + 1, cand[k].mini_sigma);
            fprintf(fname, " %8.2f", cand[k].orb_p);
            fprintf(fname, " %-7.2g ", orbperr);
            if (cand[k].psr_p < 0.001)
                fprintf(fname, " %12.5e", cand[k].psr_p);
            else
                fprintf(fname, " %12.9f", cand[k].psr_p);
            fprintf(fname, " %-7.2g ", psrperr);
            fprintf(fname, " %9.0f  ", cand[k].full_lo_r);
            fprintf(fname, " %6.0f ", cand[k].mini_N);
            fprintf(fname, " %8.1f ", cand[k].mini_r);
            fprintf(fname, " %2.0f ", cand[k].mini_numsum);
            fprintf(fname, "%7.2f ", cand[k].mini_power);
            fprintf(fname, " %.18s\n", notes + k * 18);
            fflush(fname);
        }
    }
    fprintf(fname, "\n Notes:  MH = Modulation harmonic.  ");
    fprintf(fname, "H = Pulsar harmonic.  # indicates the candidate number.\n\n");
    fclose(fname);
    sprintf(command, "cat %s.inf >> %s", name, filenm);
    system(command);
    /* This is not necessary
       sprintf(command, \
       "$PRESTO/bin/a2x -c1 -n90 -title -date -num %s > %s.ps", \
       filenm, filenm);
       system(command);
     */
    free(filenm);
}
