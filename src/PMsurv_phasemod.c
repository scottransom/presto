#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "tape_hdr.h"

/* If the following is defined, the program will print        */
/* debugging info to stdout or the results to the output dir. */
/* #define DEBUGOUT */

#ifndef TWOPI
#define TWOPI 6.2831853071795864769252867665590057683943387987502
#endif

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 100.0

/* Minimum miniFFT bin number to accept as 'real' */
#define MINMINIBIN 3

/* Minimum pulsar frequency (Hz) to search */
#define MINFREQ 1.0

/* Maximum pulsar frequency (Hz) to search */
#define MAXFREQ 900.0

/* Factor to overlap the miniFFTs (power-of-two only) */
/*    1 = no overlap of successive miniFFTs           */
/*    2 = overlap 1/2 of each miniFFT                 */
/*    4 = overlap 3/4 of each miniFFT                 */
#define OVERLAPFACT 2

/* Blocks of length maxfft to work with at a time */
#define WORKBLOCK 4

/* The power level at which to prune the initial power spectrum */
/* in order to remove some of the basic RFI etc.                */
#define PRUNELEV 25

/* The power level to reduce pruned powers to */
#define NEWLEV   5

/* Candidate file directory */
/* #define OUTDIR "/home/ransom" */
#define OUTDIR "/tmp_mnt/usr/users/sransom/results"

/* The number of miniffts to use  */
#define NUMMINIFFTS 7

/* Define USEFFTW in order to use FFTW */
/* #define USEFFTW */
#ifdef USEFFTW
#include "fftw.h"
/* FFTW plans for the FFTs */
static fftw_plan fftwplan[NUMMINIFFTS];
static fftw_plan plan;
#else
typedef struct fftw_plan {
    int decoy;
} fftw_plan;
static double *ffttables[NUMMINIFFTS], *table;
#endif

/* The lengths of the miniffts to use */
static int minifftlen[NUMMINIFFTS] = { 32, 64, 128, 256, 512, 1024, 2048 };

/* Output thresholds */
/* If any candidates have a power level greater than the   */
/* following (which corresponds to a sigma of ~6.5 when    */
/* you take into account the number of bins searched in    */
/* each minifft), then print the candidate and observation */
/* info to the output candfile.                            */

static double threshold[NUMMINIFFTS] =
    { 26.71, 27.40, 28.09, 28.79, 29.48, 30.17, 30.86 };
/* The following is ~7.5 sigma if needed
  {33.85, 34.54, 35.23, 35.93, 36.62, 37.31, 38.01};
*/

/* A single precision floating point complex structure */

typedef struct FCOMPLEX {
    float r;                    /* real */
    float i;                    /* imag */
} fcomplex;

/* The binary candidate structure to write if we find candidates */

typedef struct RAWBINCAND {
    double full_N;              /* Number of points in original time series  */
    double full_T;              /* Length (s) of original time series        */
    double full_lo_r;           /* Lowest Fourier bin that was miniFFTd      */
    double mini_N;              /* Number of points in short (mini) FFT      */
    double mini_r;              /* Candidate Fourier bin in miniFFT          */
    double mini_power;          /* Candidate normalized power in miniFFT     */
    double mini_numsum;         /* Number of powers summed to get candidate  */
    double mini_sigma;          /* Equivalent candidate sigma (for sum pow)  */
    double psr_p;               /* Approx PSR period (miniFFT center bin)    */
    double orb_p;               /* Approx orbital period (s)                 */
} rawbincand;

/* A short structure containing useful info */

typedef struct INFO {
    double dt;                  /* Sample time (s)      */
    double dm;                  /* Dispersion measure   */
    double T;                   /* Integration time     */
    int N;                      /* Number of samples    */
    int file_cntr;              /* File number on tape  */
    int ibeam;                  /* Multibeam number     */
    char tape_lbl[7];           /* Tape label           */
    char pname[17];             /* Pointing ID          */
    char outfilebase[50];       /* Base output filename */
    char rawheader[640];        /* Multibeam header     */
} info;

static float median(float arr[], int n);
#ifdef USEFFTW
static void make_plans(void);
static void destroy_plans(void);
#endif
static void info_from_header(struct tphdr *header, int N,
                             double dt, double dm, info * idata);
static int prune_powers(float *arr, int n, int numsumpows);
static void realfft(float idata[], long n, int isign);
static void cand_output(double pow, double bin, int numminifft,
                        int lorfullfft, info * idata, int candnum);
static void search_phasemod(fcomplex * minifft, int numminifft,
                            int lorfullfft, double threshold,
                            info * idata, int *numcands);
static void tablesplitfft(float data[], long nn, int isign);
static void tablefft(float data[], long nn, int isign);
static double *maketable(long nn, int isign);
static void fft_scramble(float data[], long nn);
static void tablefftraw(float data[], double table[], long n);
static void tablesplitfftraw(float data[], double table[], long n, int isign);

/******************************************************************/

int PMsurv_phasemod_search(struct tphdr *header, int N,
                           fcomplex * bigfft, double dt, double dm)
/*
 * This routine searches an FFT (assumed to be single precision        
 * complex values in increasing order by Fourier Freq) from the        
 * Parkes Multibeam Pulsar Survey for binary pulsars in very short     
 * orbits.  The method used looks for periodic phase-modulation        
 * induced sidelobes around a pulsar spin freq in the power spectrum   
 * (see http://cfa160.harvard.edu/~ransom/readable_poster.ps.gz for    
 * more details).  The inputs are:                                     
 *   *header - the 640 byte raw header for the observation             
 *   N - the number of frequencies in the input FFT 
 *       (should be 2**22 for a full length pointing)  
 *   *bigfft - the single precision full-length FFT to be searched     
 *   dt - the sample time for the FFT (should be 0.00025 for
 *       a non-decimated pointing)
 *   dm - the DM used to prepare the data in the bigfft                
 * The routine as now stands takes approximately 10-15 sec to run
 * on a full-length FFT (i.e. a 2**22 point pointing).  
 * We should be sensitive to pulsars in very-low to very-high mass     
 * binaries with orbital periods <~ 20min and flux densities of ~1 mJy  
 * or a little less (if the duty cycle of the pulsar is short enough). 
 * Return value is 0 if the code didn't find any candidates. If the 
 * return value is positive, then a significant candidate (or more) 
 * was found in the FFT (this could be the trigger to save the bigFFT 
 * to a file for a more in depth analysis...)
 */
{
    int ii, jj, worklen, fftindex, fftlen, binsleft, overlaplen;
    int bigfft_offset, max_offset, powers_offset, wrkblk = WORKBLOCK;
    int minfft, maxfft, numcands = 0;
    float *powers, *minifft, *powers_pos;
    double norm, thresh;
    info idata;

    /* Convert the Header into usable info... */
    info_from_header(header, 2 * N, dt, dm, &idata);

    /* Make the FFTW plans or twiddle tables */
#ifdef USEFFTW
    make_plans();
#else
    for (ii = 0; ii < NUMMINIFFTS; ii++)
        ffttables[ii] = maketable(minifftlen[ii] / 2, 1);
#endif

    /* Set some key variables */
    fftindex = NUMMINIFFTS - 1;
    minfft = minifftlen[0];
    maxfft = minifftlen[fftindex];
    worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
    bigfft_offset = (int) ceil(MINFREQ * idata.T);
    max_offset = (int) floor(MAXFREQ * idata.T);
    if (max_offset > N - 1)
        max_offset = N - 1;

    /* Allocate the arrays that will store the powers from */
    /* the bigFFT as well as the miniFFTs.                 */
    powers = (float *) malloc(sizeof(float) * worklen);
    minifft = (float *) malloc(sizeof(float) * maxfft);

    while (bigfft_offset < max_offset) {        /* Loop through the bigFFT */
#ifdef DEBUGOUT
        printf("\nbigfft_offset = %d\n", bigfft_offset);
#endif

        /* How close are we to the end of the freqs we want to search? */
        binsleft = max_offset - bigfft_offset;

        /* make the maximum fft size the default */
        fftindex = NUMMINIFFTS - 1;

        /* Adjust our search parameters if close to end of zone to search */
        if (binsleft < worklen) {
            wrkblk = 1;
            worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
            while (binsleft < worklen && fftindex > 0) {
                fftindex--;
                maxfft = minifftlen[fftindex];
                worklen = (wrkblk + 1) * maxfft - (maxfft / OVERLAPFACT);
            }
            if (worklen < minfft)
                break;
        }

        /* Get the powers from the bigFFT */
        for (ii = 0, jj = bigfft_offset; ii < worklen; ii++, jj++)
            powers[ii] = bigfft[jj].r * bigfft[jj].r + bigfft[jj].i * bigfft[jj].i;

        /* Chop the powers that are way above the median.  */
        /* This is a crude way of removing strong coherent */
        /* pulsations or RFI from the power spectrum.      */
        prune_powers(powers, worklen, 1);

        /* Loop through the different small FFT sizes */
        while (fftindex >= 0) {
            thresh = threshold[fftindex];
            fftlen = minifftlen[fftindex];
#ifdef USEFFTW
            plan = fftwplan[fftindex];
#else
            table = ffttables[fftindex];
#endif
#ifdef DEBUGOUT
            printf("\n%d:\n", fftlen);
#endif
            powers_pos = powers;
            powers_offset = 0;
            overlaplen = fftlen / OVERLAPFACT;

            /* Perform miniffts at each section of the powers array */
            while (powers_offset < wrkblk * maxfft) {
#ifdef DEBUGOUT
                printf("%d ", bigfft_offset + powers_offset);
#endif
                /* Copy the proper amount and portion of powers into minifft */
                memcpy(minifft, powers_pos, fftlen * sizeof(float));

                /* Perform the minifft */
                realfft(minifft, fftlen, -1);

                /* Normalize and search the miniFFT */
                norm = sqrt((double) fftlen) / minifft[0];
                for (ii = 0; ii < fftlen; ii++)
                    minifft[ii] *= norm;
                search_phasemod((fcomplex *) minifft, fftlen / 2,
                                bigfft_offset + powers_offset,
                                thresh, &idata, &numcands);

                /* Increment our data pointers */
                powers_pos += overlaplen;
                powers_offset += overlaplen;
            }

            /* Switch to the next smaller miniFFT size */
            fftindex--;
        }

        bigfft_offset += wrkblk * maxfft;
    }

    /* Free up our data arrays */

    ii = -1;
    search_phasemod(NULL, 0, 0, 0.0, &idata, &ii);
    free(powers);
    free(minifft);
    /* Free the FFTW plans or twiddle tables */
#ifdef USEFFTW
    destroy_plans();
#else
    for (ii = 0; ii < NUMMINIFFTS; ii++)
        free(ffttables[ii]);
#endif
    return numcands;
}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipies in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 */

/* Fast computation of the median of an array. */
/* Note:  It messes up the order!              */

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

static float median(float arr[], int n)
{
    int low, high;
    int median;
    int middle, ll, hh;

    low = 0;
    high = n - 1;
    median = (low + high) / 2;
    for (;;) {
        if (high <= low)        /* One element only */
            return arr[median];

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]);
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            ELEM_SWAP(arr[middle], arr[high]);
        if (arr[low] > arr[high])
            ELEM_SWAP(arr[low], arr[high]);
        if (arr[middle] > arr[low])
            ELEM_SWAP(arr[middle], arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low + 1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do
                ll++;
            while (arr[low] > arr[ll]);
            do
                hh--;
            while (arr[hh] > arr[low]);

            if (hh < ll)
                break;

            ELEM_SWAP(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]);

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}


#ifdef USEFFTW
static void make_plans(void)
{
    int ii;
    for (ii = 0; ii < NUMMINIFFTS; ii++)
        fftwplan[ii] = fftw_create_plan(minifftlen[ii] / 2, -1,
                                        FFTW_MEASURE | FFTW_USE_WISDOM |
                                        FFTW_IN_PLACE);
}

static void destroy_plans(void)
{
    int ii;
    for (ii = 0; ii < NUMMINIFFTS; ii++)
        fftw_destroy_plan(fftwplan[ii]);
}
#endif

static char *rmtrail(char *str)
/* Removes trailing space from a string */
{
    int i;

    if (str && 0 != (i = strlen(str))) {
        while (--i >= 0) {
            if (!isspace(str[i]))
                break;
        }
        str[++i] = '\0';
    }
    return str;
}


static void info_from_header(struct tphdr *header, int N,
                             double dt, double dm, info * idata)
{
    idata->dt = dt;
    idata->dm = dm;
    idata->N = N;
    idata->T = N * dt;
    strncpy(idata->tape_lbl, header->tape_lbl, 6);
    idata->tape_lbl[6] = '\0';
    rmtrail(idata->tape_lbl);
    strncpy(idata->pname, header->pname, 16);
    idata->pname[16] = '\0';
    rmtrail(idata->pname);
    idata->ibeam = strtol(header->ibeam, NULL, 10);
    idata->file_cntr = strtol(header->file_cntr, NULL, 10);
    sprintf(idata->outfilebase, "%s_File%d_Beam%02d_DM%.2f_%s",
            idata->tape_lbl, idata->file_cntr, idata->ibeam, idata->dm,
            idata->pname);
    memcpy(idata->rawheader, (char *) header, 640);
}


static int prune_powers(float *arr, int n, int numsumpows)
{
    int ii, ct = 0;
    float med, cutoff, *tmparr;

    /* Determine the median power */

    tmparr = (float *) malloc(sizeof(float) * n);
    memcpy(tmparr, arr, sizeof(float) * n);
    med = median(tmparr, n);
    free(tmparr);

    /* Throw away powers that are bigger than PRUNELEV * median */

    cutoff = med * PRUNELEV / sqrt((float) numsumpows);
    for (ii = 0; ii < n; ii++) {
        if (arr[ii] > cutoff) {
            arr[ii] = NEWLEV * med;
            ct++;
        }
    }
    return ct;
}


static void realfft(float idata[], long n, int isign)
/*  This is a modified version of the NR routine with correct (-)  */
/*  exponent.  The forward transform (i.e. normal FFT) is isign=-1 */
{
    long nby2, il, ih;
    double cc, h1r, h1i, h2r, h2i, h2rwr, h2iwr, h2rwi, h2iwi;
    double wr, wi, wpr, wpi, tmp1, theta;
    fcomplex *data;

    if (n % 2) {
        printf("\nrealfft() arrays lengths must be evenly divisible by 2.\n\n");
        exit(-1);
    }
    nby2 = n >> 1;
    data = (fcomplex *) idata;
    if (isign == -1) {
        cc = -0.5;
        theta = -TWOPI / (double) n;
#ifdef USEFFTW
        fftw_one(plan, (FFTW_COMPLEX *) idata, NULL);
#else
        tablesplitfftraw(idata, table, n / 2, isign);
        /* tablefftraw(idata, table, n/2); */
        fft_scramble(idata, n / 2);
#endif
    } else {
        cc = 0.5;
        theta = TWOPI / (double) n;
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
    }
    /* Prep the trig recursion */
    wr = cos(theta);
    wi = sin(theta);
    tmp1 = sin(0.5 * theta);
    wpr = -2.0 * tmp1 * tmp1;
    wpi = wi;
    il = 1;                     /* n     */
    ih = nby2 - il;             /* N/2-n */
    for (; il <= (n >> 2); il++, ih--) {
        h1r = 0.5 * (data[il].r + data[ih].r);
        h1i = 0.5 * (data[il].i - data[ih].i);
        h2r = -cc * (data[il].i + data[ih].i);
        h2i = cc * (data[il].r - data[ih].r);
        h2rwr = h2r * wr;
        h2rwi = h2r * wi;
        h2iwr = h2i * wr;
        h2iwi = h2i * wi;
        data[il].r = h1r + h2rwr - h2iwi;
        data[il].i = h1i + h2iwr + h2rwi;
        data[ih].r = h1r - h2rwr + h2iwi;
        data[ih].i = -h1i + h2iwr + h2rwi;
        tmp1 = wr;
        wr = tmp1 * wpr - wi * wpi + wr;
        wi = wi * wpr + tmp1 * wpi + wi;
    }
    if (isign == -1) {
        /* Set data[0].r to Freq 0 value  */
        /* Set data[0].i to Nyquist value */
        tmp1 = data[0].r;
        data[0].r = tmp1 + data[0].i;
        data[0].i = tmp1 - data[0].i;
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
    } else {
        tmp1 = data[0].r;
        data[0].r = 0.5 * (tmp1 + data[0].i);
        data[0].i = 0.5 * (tmp1 - data[0].i);
#ifdef USEFFTW
        fftw_one(plan, (FFTW_COMPLEX *) idata, NULL);
#else
        tablesplitfftraw(idata, table, n / 2, isign);
        /* tablefftraw(idata, table, n/2); */
        fft_scramble(idata, n / 2);
#endif
        tmp1 = 2.0 / (double) n;
        for (il = 0; il < n; il++)
            idata[il] *= tmp1;
    }
}


static void cand_output(double pow, double bin, int numminifft,
                        int lorfullfft, info * idata, int candnum)
/* Write the candidate */
{
    static FILE *out_binfile = NULL, *out_txtfile = NULL;
    static int filesopen = 0;
    double orbperr, psrperr;
    char filenm[100];
    rawbincand cand;

    if (candnum == 1) {
        /* Write the header */
        sprintf(filenm, "%s/%s.hdr", OUTDIR, idata->outfilebase);
        out_binfile = fopen(filenm, "wb");
        fwrite(idata->rawheader, 640, 1, out_binfile);
        fclose(out_binfile);
        /* Open the binary file */
        sprintf(filenm, "%s/%s.cand", OUTDIR, idata->outfilebase);
        out_binfile = fopen(filenm, "wb");
        /* Open the text file */
        sprintf(filenm, "%s/%s.phsout", OUTDIR, idata->outfilebase);
        out_txtfile = fopen(filenm, "w");
        fprintf(out_txtfile, "#     Tape: %s\n", idata->tape_lbl);
        fprintf(out_txtfile, "#   File #: %d\n", idata->file_cntr);
        fprintf(out_txtfile, "#     Beam: %d\n", idata->ibeam);
        fprintf(out_txtfile, "# Pointing: %s\n", idata->pname);
        fprintf(out_txtfile, "#       DM: %.2f\n", idata->dm);
        fprintf(out_txtfile, "#       dt: %.7g\n", idata->dt);
        fprintf(out_txtfile, "#        N: %d\n", idata->N);
        fprintf(out_txtfile,
                "#-----------------------------------------------------------------------------------\n");
        fprintf(out_txtfile,
                "#       Power/   P_orbit+/-Err       P_pulsar+/-Err       FullFFT  MiniFFT   MiniFFT\n");
        fprintf(out_txtfile,
                "# Cand  LocPow       (sec)                 (sec)          Low Bin   Length     Bin  \n");
        fprintf(out_txtfile,
                "#-----------------------------------------------------------------------------------\n");
        filesopen = 1;
    } else if (candnum == 0) {
        if (filesopen) {
            fclose(out_binfile);
            fclose(out_txtfile);
            filesopen = 0;
        }
        return;
    }

    /* Calculate the candidate */

    cand.full_N = idata->N;
    cand.full_T = idata->T;
    cand.full_lo_r = lorfullfft;
    cand.mini_N = 2 * numminifft;
    cand.mini_r = bin;
    cand.mini_sigma = 7.5;
    cand.mini_power = pow;
    cand.mini_numsum = 1;
    cand.psr_p = idata->T / (lorfullfft + numminifft);
    cand.orb_p = idata->T * bin / cand.mini_N;
    orbperr = 0.5 * cand.full_T / cand.mini_N;
    if (cand.full_lo_r == 0.0)
        psrperr = cand.psr_p;
    else
        psrperr = fabs(cand.full_T / (cand.full_lo_r +
                                      0.5 * cand.mini_N) -
                       cand.full_T / cand.full_lo_r);

    /*  Now output it... */
    fwrite(&cand, sizeof(rawbincand), 1, out_binfile);
    fprintf(out_txtfile, " %4d ", candnum);
    fprintf(out_txtfile, " %7.2f ", cand.mini_power);
    fprintf(out_txtfile, " %8.2f", cand.orb_p);
    fprintf(out_txtfile, " %-7.2g ", orbperr);
    if (cand.psr_p < 0.001)
        fprintf(out_txtfile, " %12.5e", cand.psr_p);
    else
        fprintf(out_txtfile, " %12.9f", cand.psr_p);
    fprintf(out_txtfile, " %-7.2g ", psrperr);
    fprintf(out_txtfile, " %9.0f  ", cand.full_lo_r);
    fprintf(out_txtfile, " %6.0f ", cand.mini_N);
    fprintf(out_txtfile, " %8.1f\n", cand.mini_r);
    fflush(out_txtfile);
    return;
}


static void search_phasemod(fcomplex * minifft, int numminifft,
                            int lorfullfft, double threshold,
                            info * idata, int *numcands)
/* Search the short FFT of a group of powers from a long FFT for */
/* pulsations that have been phase-modulated.                    */
{
    int ii, lobin, hibin, numpowers;
    double twobypi = 0.6366197723675813430755351, pow;
    fcomplex interbin;

    /* Close the outputfiles */
    if (*numcands < 0) {
        cand_output(0.0, 0.0, 0, 0, idata, 0);
        return;
    }

    /* Prep some other values we will need */
    numpowers = 2 * numminifft;
    lobin = (int) ceil(numpowers * MINORBP / idata->T);
    if (lobin < MINMINIBIN)
        lobin = MINMINIBIN;
    hibin = numminifft - 2;

    /* Loop through the minifft */
    for (ii = lobin; ii < hibin; ii++) {
        pow = minifft[ii].r * minifft[ii].r + minifft[ii].i * minifft[ii].i;
        if (pow > threshold) {
            (*numcands)++;
            cand_output(pow, ii, numminifft, lorfullfft, idata, *numcands);
        }
        interbin.r = twobypi * (minifft[ii].r - minifft[ii + 1].r);
        interbin.i = twobypi * (minifft[ii].i - minifft[ii + 1].i);
        pow = interbin.r * interbin.r + interbin.i * interbin.i;
        if (pow > threshold) {
            (*numcands)++;
            cand_output(pow, ii + 0.5, numminifft, lorfullfft, idata, *numcands);
        }
    }
}


#undef ELEM_SWAP
#undef MINORBP
#undef MINMINIBIN
#undef MINFREQ
#undef MAXFREQ
#undef OVERLAPFACT
#undef WORKBLOCK
#undef PRUNELEV
#undef NEWLEV
#undef OUTDIR
#undef NUMMINIFFTS
#ifdef DEBUGOUT
#undef DEBUGOUT
#endif

/* Various FFT routines and aux. routines */

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

static void tablesplitfft(float data[], long nn, int isign)
/*  This is a split-radix Decimation in Frequency FFT */
{
    double *table;

    table = maketable(nn, 1);
    tablesplitfftraw(data, table, nn, isign);
    fft_scramble(data, nn);
    free(table);
}


static void tablefft(float data[], long nn, int isign)
/*  This is a radix-2 Gentleman-Sande or Decimation in Frequency FFT */
{
    double *table;

    table = maketable(nn, isign);
    tablefftraw(data, table, nn);
    fft_scramble(data, nn);
    free(table);
}


static double *maketable(long nn, int isign)
{
    long i, n;
    double tmp1, wr, wpr, wpi, wi, theta;
    double *table;

    n = (nn << 1);
    table = (double *) malloc(sizeof(double) * n);
    table[0] = 1.0;
    table[1] = 0.0;
    theta = isign * (TWOPI / nn);
    wr = cos(theta);
    wi = sin(theta);
    tmp1 = sin(0.5 * theta);
    wpr = -2.0 * tmp1 * tmp1;
    wpi = wi;
    for (i = 2; i < n; i += 2) {
        table[i] = wr;
        table[i + 1] = wi;
        wr = (tmp1 = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + tmp1 * wpi + wi;
    }
    /* To check trig recursion above...
       for (i = 0; i < n; i += 2) {
       theta = isign*i*(PI/nn);
       table[i] = cos(theta);
       table[i + 1] = sin(theta);
       }
     */
    return table;
}


static void fft_scramble(float data[], long nn)
{
    long i, j, m, n;
    float tempzz;

    data--;
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}


static void tablefftraw(float data[], double table[], long n)
/*  This is a radix-2 Gentleman-Sande or Decimation in Frequency FFT */
{
    register long j, m = 2, p, q, k, n2 = n, n1, nn;
    register double c, s, rtmp, itmp;

    nn = (n << 1);
    while (m < nn) {
        n1 = n2;
        n2 >>= 1;
        for (j = 0, q = 0; j < n1; j += 2) {
            c = table[q];
            s = table[q + 1];
            q += m;
            for (k = j; k < nn; k += n1 * 2) {
                p = (k + n1);
                rtmp = data[k] - data[p];
                itmp = data[k + 1] - data[p + 1];
                data[k] += data[p];
                data[k + 1] += data[p + 1];
                data[p] = c * rtmp - s * itmp;
                data[p + 1] = c * itmp + s * rtmp;
            }
        }
        m <<= 1;
    }
}


static void tablesplitfftraw(float data[], double table[], long n, int isign)
/*  This is a split-radix Decimation in Frequency FFT */
{
    int m, n2, j, is, id;
    register int i0, n4, n3;
    register int i0i, i1i, i2i, i3i;
    double r1, r2, s1, s2, s3, cc1, ss1, cc3, ss3;
    int a, a3, ai, a3i, ndec = n - 1;
    float *x;

    /* The following is a total HACK.  See below also. */
    if (isign == 1)
        for (j = 1; j < n * 2; j += 2)
            data[j] = -data[j];
    x = data - 2;
    n2 = n << 1;
    m = 1;
    while (m < n / 2) {
        n2 >>= 1;
        n4 = n2 >> 2;
        n3 = n2 >> 1;
        a = 0;
        for (j = 1; j <= n4; j++) {
            ai = a << 1;
            a3 = (a + (a << 1)) & ndec;
            a3i = a3 << 1;
            cc1 = table[ai];
            ss1 = table[ai + 1];
            cc3 = table[a3i];
            ss3 = table[a3i + 1];
            a = (a + m) & ndec;
            is = j;
            id = n2 << 1;
            do {
                for (i0 = is; i0 <= n - 1; i0 += id) {
                    i0i = i0 << 1;
                    i1i = i0i + n3;
                    i2i = i1i + n3;
                    i3i = i2i + n3;
                    r1 = x[i0i] - x[i2i];
                    x[i0i] += x[i2i];
                    r2 = x[i1i] - x[i3i];
                    x[i1i] += x[i3i];
                    s1 = x[i0i + 1] - x[i2i + 1];
                    x[i0i + 1] += x[i2i + 1];
                    s2 = x[i1i + 1] - x[i3i + 1];
                    x[i1i + 1] += x[i3i + 1];
                    s3 = r1 - s2;
                    r1 += s2;
                    s2 = r2 - s1;
                    r2 += s1;
                    x[i2i] = r1 * cc1 - s2 * ss1;
                    x[i2i + 1] = -s2 * cc1 - r1 * ss1;
                    x[i3i] = s3 * cc3 + r2 * ss3;
                    x[i3i + 1] = r2 * cc3 - s3 * ss3;
                }
                is = (id << 1) - n2 + j;
                id <<= 2;
            }
            while (is < n);
        }
        m <<= 1;
    }
    is = 1;
    id = 4;
    do {
        for (i0 = is; i0 <= n; i0 += id) {
            i0i = i0 << 1;
            i1i = i0i + 2;
            r1 = x[i0i];
            x[i0i] = r1 + x[i1i];
            x[i1i] = r1 - x[i1i];
            r1 = x[i0i + 1];
            x[i0i + 1] = r1 + x[i1i + 1];
            x[i1i + 1] = r1 - x[i1i + 1];
        }
        is = (id << 1) - 1;
        id <<= 2;
    }
    while (is < n);
    /* The following is a total HACK. */
    if (isign == 1)
        for (j = 1; j < n * 2; j += 2)
            data[j] = -data[j];
}

#undef TWOPI
#ifdef SWAP
#undef SWAP
#endif
