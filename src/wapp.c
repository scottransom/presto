#include "presto.h"
#include "mask.h"
#include "wapp.h"

/*  NOTES:
bytesperblk_st is the number of bytes in the RAW LAGS for 
               a _single_ WAPP
sampperblk_st  is the number of samples (i.e. lags) in a block 
               for all WAPPs together
*/

/* All of the following have an _st to indicate static */
static long long numpts_st[MAXPATCHFILES], padpts_st[MAXPATCHFILES], N_st;
static long long filedatalen_st[MAXPATCHFILES];
static int numblks_st[MAXPATCHFILES], corr_level_st, decreasing_freqs_st = 0;
static int bytesperpt_st, bytesperblk_st, bits_per_samp_st;
static int numwapps_st = 0, numwappchan_st, numchan_st, numifs_st, ptsperblk_st;
static int need_byteswap_st, sampperblk_st, usewindow_st = 0;
static double times_st[MAXPATCHFILES], mjds_st[MAXPATCHFILES];
static double elapsed_st[MAXPATCHFILES], T_st, dt_st, dtus_st;
static double startblk_st[MAXPATCHFILES], endblk_st[MAXPATCHFILES];
static double corr_rate_st, corr_scale_st;
static double center_freqs_st[WAPP_MAXNUMWAPPS], *window_st = NULL;
static infodata idata_st[MAXPATCHFILES];
//static unsigned char padvals[WAPP_MAXLAGLEN], newpadvals[WAPP_MAXLAGLEN];
//static unsigned char databuffer[2 * WAPP_MAXDATLEN], padval = 128;
//static unsigned char lagbuffer[WAPP_MAXLAGLEN];
static int currentfile, currentblock;
static int header_version_st, header_size_st;
//static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
static fftwf_plan fftplan;
static float clip_sigma_st = 0.0, *lags = NULL;
//static int using_MPI = 0;

double slaCldj(int iy, int im, int id, int *j);
extern short transpose_bytes(unsigned char *a, int nx, int ny, unsigned char *move,
                             int move_size);

static double *hanning_window(int numlags)
{
    double *win;
    int ii;

    win = gen_dvect(numlags);
    for (ii = 0; ii < numlags; ii++)
        win[ii] = 0.5 + 0.5 * cos(PI * ii / (numlags - 1));
    return win;
}

char *get_hdr_string(struct HEADERP *h, char *name, int *slen)
{
    struct HEADERVAL val;
    char *cptr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    cptr = (char *) val.value;
    *slen = strlen(cptr);
    return cptr;
}

double get_hdr_double(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    double dval;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dval = *((double *) val.value);
    if (need_byteswap_st) {
        dval = swap_double(dval);
    }
    return dval;
}

int get_hdr_int(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    int ival;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    ival = *((int *) val.value);
    if (need_byteswap_st) {
        ival = swap_int(ival);
    }
    return ival;
}

long long get_hdr_longlong(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    long long llval;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    llval = *((long long *) val.value);
    if (need_byteswap_st) {
        llval = swap_longlong(llval);
    }
    return llval;
}

double *get_hdr_double_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii;
    double *dptr, *darr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dptr = (double *) val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    darr = gen_dvect(*len);
    if (need_byteswap_st) {
        for (ii = 0; ii < *len; ii++, dptr++) {
            darr[ii] = swap_double(*dptr);
        }
    } else {
        for (ii = 0; ii < *len; ii++, dptr++) {
            darr[ii] = *dptr;
        }
    }
    return darr;
}

int *get_hdr_int_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii, *iptr, *iarr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    iptr = (int *) val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    iarr = gen_ivect(*len);
    if (need_byteswap_st) {
        for (ii = 0; ii < *len; ii++, iptr++) {
            iarr[ii] = swap_int(*iptr);
        }
    } else {
        for (ii = 0; ii < *len; ii++, iptr++) {
            iarr[ii] = *iptr;
        }
    }
    return iarr;
}

void set_WAPP_HEADER_version(struct HEADERP *hdr)
{
    int bin_hdrlen, ascii_hdrlen;

    bin_hdrlen = hdr->headlen;
    ascii_hdrlen = hdr->offset;
    header_version_st = get_hdr_int(hdr, "header_version");
    header_size_st = bin_hdrlen + ascii_hdrlen;

    /* The following tries to determine if we need to byteswap */
    if ((header_version_st < 1 || header_version_st > 15) &&
        ((bin_hdrlen < 1000 || bin_hdrlen > 4000) ||
         (ascii_hdrlen < 1000 || ascii_hdrlen > 40000))) {
        header_version_st = swap_int(header_version_st);
        bin_hdrlen = swap_int(bin_hdrlen);
        ascii_hdrlen = swap_int(ascii_hdrlen);
        header_size_st = bin_hdrlen + ascii_hdrlen;
        need_byteswap_st = 1;
    }
#if 0
    printf("Header version:  %d\n", header_version_st);
    printf("Header  length:  %d\n", header_size_st);
#endif
}

static double UT_strings_to_MJD(char *obs_date, char *start_time,
                                int *mjd_day, double *mjd_fracday)
{
    int year, month, day, hour, min, sec, err;

    sscanf(obs_date, "%4d%2d%2d", &year, &month, &day);
    sscanf(start_time, "%2d:%2d:%2d", &hour, &min, &sec);
    *mjd_fracday = (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
    *mjd_day = slaCldj(year, month, day, &err);
    return *mjd_day + *mjd_fracday;
}

static double wappcorrect(double mjd)
/* subroutine to return correction to wapp_time (us) based on mjd */
{
    double correction;

    /* assume no correction initially */
    correction = 0.0;

    if ((mjd >= 51829.0) && (mjd < 51834.0))
        correction = -0.08;
    if ((mjd >= 51834.0) && (mjd < 51854.0))
        correction = -0.68;
    if ((mjd >= 51854.0) && (mjd < 51969.0))
        correction = +0.04;

    if (correction != 0.0) {
        fprintf(stderr, "WARNING: correction %f us applied for MJD %.1f\n",
                correction, mjd);
        fflush(stderr);
    }

    return (correction);
}

static void WAPP_hdr_to_inf(struct HEADERP *h, infodata * idata)
/* Convert WAPP header into an infodata structure */
{
    int len;
    double MJD, dval;
    char ctmp[80], *cptr1, *cptr2;

    cptr1 = get_hdr_string(h, "src_name", &len);
    strncpy(idata->object, cptr1, 24);
    dval = get_hdr_double(h, "src_ra");
    idata->ra_h = (int) floor(dval / 10000.0);
    idata->ra_m = (int) floor((dval - idata->ra_h * 10000) / 100.0);
    idata->ra_s = dval - idata->ra_h * 10000 - idata->ra_m * 100;
    dval = get_hdr_double(h, "src_dec");
    idata->dec_d = (int) floor(fabs(dval) / 10000.0);
    idata->dec_m = (int) floor((fabs(dval) - idata->dec_d * 10000) / 100.0);
    idata->dec_s = fabs(dval) - idata->dec_d * 10000 - idata->dec_m * 100;
    if (dval < 0.0)
        idata->dec_d = -idata->dec_d;
    strcpy(idata->telescope, "Arecibo");
    strcpy(idata->instrument, "WAPP");
    idata->num_chan = get_hdr_int(h, "num_lags");
    cptr1 = get_hdr_string(h, "obs_date", &len);
    cptr2 = get_hdr_string(h, "start_time", &len);
    MJD = UT_strings_to_MJD(cptr1, cptr2, &(idata->mjd_i), &(idata->mjd_f));
    idata->dt = (wappcorrect(MJD) + get_hdr_double(h, "wapp_time")) / 1000000.0;
    /* This is to allow folding starting from files that aren't the first of an obs */
    idata->mjd_f += get_hdr_longlong(h, "timeoff") * idata->dt / SECPERDAY;
    if (idata->mjd_f > 1.0) {
        idata->mjd_f -= 1.0;
        idata->mjd_i += 1;
    }
    idata->N = get_hdr_double(h, "obs_time") / idata->dt;
    idata->freqband = get_hdr_double(h, "bandwidth");
    idata->chan_wid = fabs(idata->freqband / idata->num_chan);
    idata->freq =
        get_hdr_double(h,
                       "cent_freq") - 0.5 * idata->freqband + 0.5 * idata->chan_wid;
    idata->fov = 1.2 * SOL * 3600.0 / (1000000.0 * idata->freq * 300.0 * DEGTORAD);
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    strcpy(idata->analyzer, "Scott Ransom");
    cptr1 = get_hdr_string(h, "observers", &len);
    strncpy(idata->observer, cptr1, 24);
    if (get_hdr_int(h, "sum"))
        sprintf(ctmp,
                "%d %d-level IF(s) were summed.  Lags are %d bit ints.",
                numifs_st, corr_level_st, bits_per_samp_st);
    else
        sprintf(ctmp, "%d %d-level IF(s) were not summed.  Lags are %d bit ints.",
                numifs_st, corr_level_st, bits_per_samp_st);
    sprintf(idata->notes,
            "Starting Azimuth (deg) = %.15g,  Zenith angle (deg) = %.15g\n"
            "Project ID %s, Scan number %d, Date: %s %s.\n    %s\n",
            get_hdr_double(h, "start_az"),
            get_hdr_double(h, "start_za"),
            get_hdr_string(h, "project_id", &len),
            get_hdr_int(h, "scan_number"),
            get_hdr_string(h, "obs_date", &len),
            get_hdr_string(h, "start_time", &len), ctmp);
}


void get_WAPP_file_info(FILE * files[], int numwapps, int numfiles, int usewindow,
                        float clipsig, long long *N, int *ptsperblock,
                        int *numchan, double *dt, double *T,
                        infodata * idata, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of WAPP rawfiles that you want to patch      */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
    int ii, jj, ival;
    double dval;
    struct HEADERP *hdr, *hdr2;

    if (numfiles > MAXPATCHFILES) {
        printf("\nThe number of input files (%d) is greater than \n", numfiles);
        printf("   MAXPATCHFILES=%d.  Exiting.\n\n", MAXPATCHFILES);
        exit(0);
    }
    numwapps_st = numwapps;

    /* Read the header with the yacc/lex generated tools */
    hdr = head_parse(files[0]);
    /* Check the header version and find out the header offsets */
    set_WAPP_HEADER_version(hdr);

    /* Skip the ASCII and binary headers of all the WAPP files */
    for (ii = 0; ii < numwapps_st; ii++)
        chkfseek(files[ii], header_size_st, SEEK_SET);

    /* Now start getting the technical info */
    numifs_st = get_hdr_int(hdr, "nifs");

    if (get_hdr_int(hdr, "freqinversion")) {
        decreasing_freqs_st = (decreasing_freqs_st == 1 ? 0 : 1);
    }
    if (get_hdr_int(hdr, "iflo_flip")) {
        decreasing_freqs_st = (decreasing_freqs_st == 1 ? 0 : 1);
    }
    printf("freqinversion = %d   iflo_flip = %d : using decreasing_freqs = %d\n",
           get_hdr_int(hdr, "freqinversion"),
           get_hdr_int(hdr, "iflo_flip"), decreasing_freqs_st);

    ival = get_hdr_int(hdr, "level");
    if (ival == 1)
        corr_level_st = 3;
    else if (ival == 2)
        corr_level_st = 9;
    else
        printf("\nERROR:  Unrecognized level setting!\n\n");

    ival = get_hdr_int(hdr, "lagformat");
    if (ival == 0)
        bits_per_samp_st = 16;
    else if (ival == 1)
        bits_per_samp_st = 32;
    else
        printf("\nERROR:  Unrecognized number of bits per sample!\n\n");

    /* Quick hack to allow offsets of the WAPP center freq without recompiling */
    dval = get_hdr_double(hdr, "cent_freq");
    {
        char *envval = getenv("WAPP_FREQ_ADJ");
        if (envval != NULL) {
            double dblval = strtod(envval, NULL);
            if (dblval) {
                printf
                    ("Offsetting band by %.4g MHz as per WAPP_FREQ_ADJ env variable.\n",
                     dblval);
                dval += dblval;
            }
        }
    }
    center_freqs_st[0] = dval;

    if (numifs_st == 2)
        printf("Both IFs are present.\n");
    WAPP_hdr_to_inf(hdr, &idata_st[0]);
    WAPP_hdr_to_inf(hdr, idata);
    for (ii = 1; ii < numwapps_st; ii++)
        center_freqs_st[ii] = center_freqs_st[0] + ii * idata->freqband;
    /* Hack to invert band when being used for very low frequencies */
    if (center_freqs_st[0] < 400.0) {
        decreasing_freqs_st = 1;
        printf("Inverting the band since the center frequency is < 400MHz...\n");
    }
    /* Hack to invert band when using the 12.5 MHz mode */
    if (center_freqs_st[0] > 400.0 && idata_st[0].freqband < 13.0) {
        decreasing_freqs_st = 1;
        printf("Inverting the band since the BW < 12.5 MHz...\n");
    }
    /* Are we going to clip the data? */
    if (clipsig > 0.0)
        clip_sigma_st = clipsig;
    numwappchan_st = idata_st[0].num_chan;
    *numchan = numchan_st = numwapps_st * numwappchan_st;
    /* The following should be freed sometime... */
    lags = (float *) fftwf_malloc((numwappchan_st + 1) * sizeof(float));
    fftplan =
        fftwf_plan_r2r_1d(numwappchan_st + 1, lags, lags, FFTW_REDFT00,
                          FFTW_PATIENT);
    if (usewindow) {
        usewindow_st = 1;
        printf("Calculated Hanning window for use.\n");
        /* Note:  Since the lags we get are only half of the lags that   */
        /* we really need to FFT in order to get spectra (i.e. the       */
        /* transform that we compute is real and even so we comute a     */
        /* DCT-I instead of an FFT), we will multiply the lags by the    */
        /* second half of the window.  The other half of the data (which */
        /* we don't store since it is redundant)  gets the 2st half of   */
        /* the window implicitly since the data wraps around.            */
        window_st = hanning_window(numwappchan_st);
    }
    /* Calculate the maximum number of points we can have in a */
    /* block (power of two), based on the number of samples in */
    /* each file.                                              */
    filedatalen_st[0] = chkfilelen(files[0], 1) - header_size_st;
    bytesperpt_st = (numwappchan_st * numifs_st * bits_per_samp_st) / 8;
    numpts_st[0] = filedatalen_st[0] / bytesperpt_st;
    if (filedatalen_st[0] % bytesperpt_st)
        printf("\nWARNING!!!:\n\t"
               "File 0 has a non-integer number of complete samples!\n\n");
    ptsperblk_st = WAPP_MAXPTSPERBLOCK;
    while (numpts_st[0] % ptsperblk_st)
        ptsperblk_st /= 2;
    bytesperblk_st = ptsperblk_st * bytesperpt_st;
    if (filedatalen_st[0] % bytesperblk_st)
        printf("\nWARNING!!!:\n\t"
               "File 0 has a non-integer number of complete blocks!\n\n");
    *ptsperblock = ptsperblk_st;
    sampperblk_st = ptsperblk_st * numchan_st;
    numblks_st[0] = filedatalen_st[0] / bytesperblk_st;
    N_st = numpts_st[0];
    dtus_st = idata_st[0].dt * 1000000.0;
    corr_rate_st = 1.0 / (dtus_st - WAPP_DEADTIME);
    corr_scale_st = corr_rate_st / idata->freqband;
    /* Correction for narrow band use */
    if (idata->freqband < 50.0)
        corr_scale_st = corr_rate_st / 50.0;
    if (corr_level_st == 9)     /* 9-level sampling */
        corr_scale_st /= 16.0;
    if (get_hdr_int(hdr, "sum"))        /* summed IFs (search mode) */
        corr_scale_st /= 2.0;
    corr_scale_st *= pow(2.0, (double) get_hdr_int(hdr, "lagtrunc"));
    idata->freqband *= numwapps_st;
    idata->num_chan *= numwapps_st;
    dt_st = *dt = idata_st[0].dt;
    times_st[0] = numpts_st[0] * dt_st;
    mjds_st[0] = idata_st[0].mjd_i + idata_st[0].mjd_f;
    elapsed_st[0] = 0.0;
    startblk_st[0] = 1;
    endblk_st[0] = (double) numpts_st[0] / ptsperblk_st;
    padpts_st[0] = padpts_st[numfiles / numwapps_st - 1] = 0;
    for (ii = 1; ii < numfiles / numwapps_st; ii++) {
        /* Read the header with the yacc/lex generated tools */
        hdr2 = head_parse(files[ii * numwapps_st]);
        /* Skip the ASCII and binary headers of the other WAPP files */
        for (jj = 0; jj < numwapps_st; jj++)
            chkfseek(files[ii * numwapps_st + jj], header_size_st, SEEK_SET);
        WAPP_hdr_to_inf(hdr2, &idata_st[ii]);
        close_parse(hdr2);
        if (idata_st[ii].num_chan != numwappchan_st) {
            printf("Number of channels (file %d) is not the same!\n\n", ii + 1);
        }
        if (idata_st[ii].dt != dt_st) {
            printf("Sample time (file %d) is not the same!\n\n", ii + 1);
        }
        filedatalen_st[ii] = chkfilelen(files[ii * numwapps_st], 1) - header_size_st;
        numblks_st[ii] = filedatalen_st[ii] / bytesperblk_st;
        numpts_st[ii] = numblks_st[ii] * ptsperblk_st;
        times_st[ii] = numpts_st[ii] * dt_st;
        /* If the MJDs are equal, then this is a continuation */
        /* file.  In that case, calculate the _real_ time     */
        /* length of the previous file and add it to the      */
        /* previous files MJD to get the current MJD.         */
        mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        if (fabs(mjds_st[ii] - mjds_st[0]) < 1.0e-6 / SECPERDAY) {
            elapsed_st[ii] = (filedatalen_st[ii - 1] / bytesperpt_st) * dt_st;
            idata_st[ii].mjd_f = idata_st[ii - 1].mjd_f + elapsed_st[ii] / SECPERDAY;
            idata_st[ii].mjd_i = idata_st[ii - 1].mjd_i;
            if (idata_st[ii].mjd_f >= 1.0) {
                idata_st[ii].mjd_f -= 1.0;
                idata_st[ii].mjd_i++;
            }
            mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        } else {
            elapsed_st[ii] = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                                          idata_st[ii - 1].mjd_i,
                                          idata_st[ii - 1].mjd_f);
        }
        padpts_st[ii - 1] =
            (long long) ((elapsed_st[ii] - times_st[ii - 1]) / dt_st + 0.5);
        elapsed_st[ii] += elapsed_st[ii - 1];
        N_st += numpts_st[ii] + padpts_st[ii - 1];
        startblk_st[ii] = (double) (N_st - numpts_st[ii]) / ptsperblk_st + 1;
        endblk_st[ii] = (double) (N_st) / ptsperblk_st;
    }
    padpts_st[numfiles / numwapps_st - 1] =
        ((long long) ceil(endblk_st[numfiles / numwapps_st - 1]) * ptsperblk_st -
         N_st);
    N_st += padpts_st[numfiles / numwapps_st - 1];
    *N = N_st;
    *T = T_st = N_st * dt_st;
    currentfile = currentblock = 0;
    if (output) {
        printf("   Number of files = %d\n", numfiles);
        printf(" Num of files/WAPP = %d\n", numfiles / numwapps_st);
        printf("      Points/block = %d\n", ptsperblk_st);
        printf("     Channels/WAPP = %d\n", numwappchan_st);
        printf("   Num of channels = %d\n", numchan_st);
        for (ii = 0; ii < numwapps_st; ii++)
            printf(" WAPP%d center freq = %.8g\n", ii + 1, center_freqs_st[ii]);
        printf("  Total points (N) = %lld\n", N_st);
        printf("  Sample time (dt) = %-14.14g\n", dt_st);
        printf("    Total time (s) = %-14.14g\n", T_st);
        printf("  ASCII Header (B) = %d\n", header_size_st -
               get_hdr_int(hdr, "header_size"));
        printf(" Binary Header (B) = %d\n\n", get_hdr_int(hdr, "header_size"));
        printf
            ("File  Start Block    Last Block     Points      Elapsed (s)      Time (s)            MJD           Padding\n");
        printf
            ("----  ------------  ------------  ----------  --------------  --------------  ------------------  ----------\n");
        for (ii = 0; ii < numfiles / numwapps_st; ii++)
            printf
                ("%2d    %12.11g  %12.11g  %10lld  %14.13g  %14.13g  %17.12f  %10lld\n",
                 ii + 1, startblk_st[ii], endblk_st[ii], numpts_st[ii],
                 elapsed_st[ii], times_st[ii], mjds_st[ii], padpts_st[ii]);
        printf("\n");
    }
    close_parse(hdr);
}


void WAPP_update_infodata(int numfiles, infodata * idata)
/* Update the onoff bins section in case we used multiple files */
{
    int ii, index = 2;

    idata->N = N_st;
    if (numfiles / numwapps_st == 1 && padpts_st[0] == 0) {
        idata->numonoff = 0;
        return;
    }
    /* Determine the topocentric onoff bins */
    idata->numonoff = 1;
    idata->onoff[0] = 0.0;
    idata->onoff[1] = numpts_st[0] - 1.0;
    for (ii = 1; ii < numfiles / numwapps_st; ii++) {
        if (padpts_st[ii - 1]) {
            idata->onoff[index] = idata->onoff[index - 1] + padpts_st[ii - 1];
            idata->onoff[index + 1] = idata->onoff[index] + numpts_st[ii];
            idata->numonoff++;
            index += 2;
        } else {
            idata->onoff[index - 1] += numpts_st[ii];
        }
    }
    if (padpts_st[numfiles / numwapps_st - 1]) {
        idata->onoff[index] =
            idata->onoff[index - 1] + padpts_st[numfiles / numwapps_st - 1];
        idata->onoff[index + 1] = idata->onoff[index];
        idata->numonoff++;
    }
}


void print_WAPP_hdr(struct HEADERP *hdr)
/* Output a WAPP header in human readable form */
{
    int mjd_i, len;
    double mjd_d;

    printf("\n             Header version = %d\n",
           get_hdr_int(hdr, "header_version"));
    printf("  ASCII Header size (bytes) = %d\n",
           header_size_st - get_hdr_int(hdr, "header_size"));
    printf(" Binary Header size (bytes) = %d\n", get_hdr_int(hdr, "header_size"));
    printf("                Source Name = %s\n",
           get_hdr_string(hdr, "src_name", &len));
    printf("           Observation Type = %s\n",
           get_hdr_string(hdr, "obs_type", &len));
    printf(" Observation Date (YYYMMDD) = %s\n",
           get_hdr_string(hdr, "obs_date", &len));
    printf("    Obs Start UT (HH:MM:SS) = %s\n",
           get_hdr_string(hdr, "start_time", &len));
    printf("             MJD start time = %.12f\n",
           UT_strings_to_MJD(get_hdr_string(hdr, "obs_date", &len),
                             get_hdr_string(hdr, "start_time", &len), &mjd_i,
                             &mjd_d));
    printf("                 Project ID = %s\n",
           get_hdr_string(hdr, "project_id", &len));
    printf("                  Observers = %s\n",
           get_hdr_string(hdr, "observers", &len));
    printf("                Scan Number = %d\n", get_hdr_int(hdr, "scan_number"));
    printf("    RA (J2000, HHMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_ra"));
    printf("   DEC (J2000, DDMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_dec"));
    printf("        Start Azimuth (deg) = %-17.15g\n",
           get_hdr_double(hdr, "start_az"));
    printf("     Start Zenith Ang (deg) = %-17.15g\n",
           get_hdr_double(hdr, "start_za"));
    printf("            Start AST (sec) = %-17.15g\n",
           get_hdr_double(hdr, "start_ast"));
    printf("            Start LST (sec) = %-17.15g\n",
           get_hdr_double(hdr, "start_lst"));
    printf("           Obs Length (sec) = %-17.15g\n",
           get_hdr_double(hdr, "obs_time"));
    printf("      Requested T_samp (us) = %-17.15g\n",
           get_hdr_double(hdr, "samp_time"));
    printf("         Actual T_samp (us) = %-17.15g\n",
           get_hdr_double(hdr, "wapp_time"));
    printf("         Central freq (MHz) = %-17.15g\n",
           get_hdr_double(hdr, "cent_freq"));
    printf("      Total Bandwidth (MHz) = %-17.15g\n",
           get_hdr_double(hdr, "bandwidth"));
    printf("             Number of lags = %d\n", get_hdr_int(hdr, "num_lags"));
    printf("              Number of IFs = %d\n", get_hdr_int(hdr, "nifs"));
    printf("    Samples since obs start = %lld\n", get_hdr_longlong(hdr, "timeoff"));
    printf("   Other information:\n");
    if (get_hdr_int(hdr, "sum") == 1)
        printf("      IFs are summed.\n");
    if (header_version_st < 7) {
        if (get_hdr_int(hdr, "freqinversion")) {
            decreasing_freqs_st = 1;
        }
    } else {
        if (get_hdr_int(hdr, "iflo_flip")) {
            decreasing_freqs_st = 1;
        }
        if (get_hdr_int(hdr, "isalfa")) {
            printf("      Data is from the ALFA receiver.\n");
            decreasing_freqs_st = 1;
        }
    }
    if (decreasing_freqs_st)
        printf("      Frequency band is inverted.\n");
    if (get_hdr_int(hdr, "lagformat") == 0)
        printf("      Lags are 16 bit integers.\n\n");
    else
        printf("      Lags are 32 bit integers.\n\n");
}
