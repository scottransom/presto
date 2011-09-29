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
static unsigned char padvals[WAPP_MAXLAGLEN], newpadvals[WAPP_MAXLAGLEN];
static unsigned char databuffer[2 * WAPP_MAXDATLEN], padval = 128;
static unsigned char lagbuffer[WAPP_MAXLAGLEN];
static int currentfile, currentblock;
static int header_version_st, header_size_st;
static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
static fftwf_plan fftplan;
static float clip_sigma_st = 0.0, *lags = NULL;
static int using_MPI = 0;

double slaCldj(int iy, int im, int id, int *j);
static double inv_cerf(double input);
static void vanvleck3lev(float *rho, int npts);
static void vanvleck9lev(float *rho, int npts);
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

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    cptr = (char *)val.value;
    *slen = strlen(cptr);
    return cptr;
}

double get_hdr_double(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    double dval;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dval = *((double *)val.value);
    if (need_byteswap_st) {
        dval = swap_double(dval);
    }
    return dval;
}

int get_hdr_int(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    int ival;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    ival = *((int *)val.value);
    if (need_byteswap_st) {
        ival = swap_int(ival);
    }
    return ival;
}

long long get_hdr_longlong(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    long long llval;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    llval = *((long long *)val.value);
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

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dptr = (double *)val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    darr = gen_dvect(*len);
    if (need_byteswap_st) {
        for (ii=0; ii<*len; ii++, dptr++) {
            darr[ii] = swap_double(*dptr);
        } 
    } else {
        for (ii=0; ii<*len; ii++, dptr++) {
            darr[ii] = *dptr;
        }
    }
    return darr;
}

int *get_hdr_int_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii, *iptr, *iarr;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    iptr = (int *)val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    iarr = gen_ivect(*len);
    if (need_byteswap_st) {
        for (ii=0; ii<*len; ii++, iptr++) {
            iarr[ii] = swap_int(*iptr);
        } 
    } else {
        for (ii=0; ii<*len; ii++, iptr++) {
            iarr[ii] = *iptr;
        }
    }
    return iarr;
}

void get_WAPP_static(int *bytesperpt, int *bytesperblk, int *numifs,
                     float *clip_sigma)
{
   *bytesperpt = bytesperpt_st;
   *bytesperblk = bytesperblk_st;
   *numifs = numifs_st;
   *clip_sigma = clip_sigma_st;
}

void set_WAPP_static(int ptsperblk, int bytesperpt, int bytesperblk,
                     int numchan, int numifs, float clip_sigma, double dt)
{
   using_MPI = 1;
   currentblock = 0;
   ptsperblk_st = ptsperblk;
   bytesperpt_st = bytesperpt;
   bytesperblk_st = bytesperblk;
   numchan_st = numchan;
   numifs_st = numifs;
   sampperblk_st = ptsperblk_st * numchan_st;
   clip_sigma_st = clip_sigma;
   dt_st = dt;
}

void set_WAPP_padvals(float *fpadvals, int good_padvals)
{
   int ii;
   float sum_padvals = 0.0;

   if (good_padvals) {
      for (ii = 0; ii < numchan_st; ii++) {
         padvals[ii] = newpadvals[ii] = (unsigned char) (fpadvals[ii] + 0.5);
         sum_padvals += fpadvals[ii];
      }
      padval = (unsigned char) (sum_padvals / numchan_st + 0.5);
   } else {
      for (ii = 0; ii < numchan_st; ii++)
         padvals[ii] = newpadvals[ii] = padval;
   }
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
         (ascii_hdrlen < 1000 || ascii_hdrlen > 40000))){
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
        get_hdr_double(h, "cent_freq") - 0.5 * idata->freqband + 0.5 * idata->chan_wid;
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
        decreasing_freqs_st = (decreasing_freqs_st==1 ? 0 : 1);
    }
    if (get_hdr_int(hdr, "iflo_flip")) {
        decreasing_freqs_st = (decreasing_freqs_st==1 ? 0 : 1);
    }
    printf("freqinversion = %d   iflo_flip = %d : using decreasing_freqs = %d\n", 
           get_hdr_int(hdr, "freqinversion"),
           get_hdr_int(hdr, "iflo_flip"),
           decreasing_freqs_st);

    ival = get_hdr_int(hdr, "level");
    if (ival==1)
        corr_level_st = 3;
    else if (ival==2)
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
        fftwf_plan_r2r_1d(numwappchan_st + 1, lags, lags, FFTW_REDFT00, FFTW_PATIENT);
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
    if (corr_level_st == 9)      /* 9-level sampling */
        corr_scale_st /= 16.0;
    if (get_hdr_int(hdr, "sum")) /* summed IFs (search mode) */
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
        filedatalen_st[ii] = chkfilelen(files[ii * numwapps_st], 1) -
            header_size_st;
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
                 ii + 1, startblk_st[ii], endblk_st[ii], numpts_st[ii], elapsed_st[ii],
                 times_st[ii], mjds_st[ii], padpts_st[ii]);
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


int skip_to_WAPP_rec(FILE * infiles[], int numfiles, int rec)
/* This routine skips to the record 'rec' in the input files   */
/* *infiles.  *infiles contain data from the WAPP at Arecibo   */
/* Returns the record skipped to.                              */
{
   double floor_blk;
   int filenum = 0, ii;

   if (rec < startblk_st[0])
      rec += (startblk_st[0] - 1);
   if (rec > 0 && rec < endblk_st[numfiles / numwapps_st - 1]) {

      /* Find which file we need */
      while (rec > endblk_st[filenum])
         filenum++;

      currentblock = rec - 1;
      shiftbuffer = 1;
      floor_blk = floor(startblk_st[filenum]);

      /* Set the data buffer to all padding just in case */
      memset(databuffer, padval, 2 * WAPP_MAXDATLEN);

      /* Warning:  I'm not sure if the following is correct. */
      /*   If really needs accurate testing to see if my     */
      /*   offsets are correct.  Bottom line, don't trust    */
      /*   a TOA determined using the following!             */

      if (rec < startblk_st[filenum]) { /* Padding region */
         currentfile = filenum - 1;
         for (ii = 0; ii < numwapps_st; ii++)
            chkfileseek(infiles[currentfile * numwapps_st + ii], 0, 1, SEEK_END);
         bufferpts = padpts_st[currentfile] % ptsperblk_st;
         padnum = ptsperblk_st * (rec - endblk_st[currentfile] - 1);
         /*
            printf("Padding:  currentfile = %d  bufferpts = %d  padnum = %d\n",
            currentfile, bufferpts, padnum);
          */
      } else {                  /* Data region */
         currentfile = filenum;
         for (ii = 0; ii < numwapps_st; ii++)
            chkfileseek(infiles[currentfile * numwapps_st + ii],
                        rec - startblk_st[filenum], bytesperblk_st, SEEK_CUR);
         bufferpts = (int) ((startblk_st[filenum] - floor_blk) * ptsperblk_st + 0.5);
         padnum = 0;
         /*
            printf("Data:  currentfile = %d  bufferpts = %d  padnum = %d\n",
            currentfile, bufferpts, padnum);
          */
      }

   } else {
      printf("\n rec = %d out of range in skip_to_WAPP_rec()\n", rec);
      exit(1);
   }
   return rec;
}


void print_WAPP_hdr(struct HEADERP *hdr)
/* Output a WAPP header in human readable form */
{
    int mjd_i, len;
    double mjd_d;

    printf("\n             Header version = %d\n", get_hdr_int(hdr, "header_version"));
    printf("  ASCII Header size (bytes) = %d\n", header_size_st -
           get_hdr_int(hdr, "header_size"));
    printf(" Binary Header size (bytes) = %d\n", get_hdr_int(hdr, "header_size"));
    printf("                Source Name = %s\n", get_hdr_string(hdr, "src_name", &len));
    printf("           Observation Type = %s\n", get_hdr_string(hdr, "obs_type", &len));
    printf(" Observation Date (YYYMMDD) = %s\n", get_hdr_string(hdr, "obs_date", &len));
    printf("    Obs Start UT (HH:MM:SS) = %s\n", get_hdr_string(hdr, "start_time", &len));
    printf("             MJD start time = %.12f\n",
           UT_strings_to_MJD(get_hdr_string(hdr, "obs_date", &len), 
                             get_hdr_string(hdr, "start_time", &len), &mjd_i, &mjd_d));
    printf("                 Project ID = %s\n", get_hdr_string(hdr, "project_id", &len));
    printf("                  Observers = %s\n", get_hdr_string(hdr, "observers", &len));
    printf("                Scan Number = %d\n", get_hdr_int(hdr, "scan_number"));
    printf("    RA (J2000, HHMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_ra"));
    printf("   DEC (J2000, DDMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_dec"));
    printf("        Start Azimuth (deg) = %-17.15g\n", get_hdr_double(hdr, "start_az"));
    printf("     Start Zenith Ang (deg) = %-17.15g\n", get_hdr_double(hdr, "start_za"));
    printf("            Start AST (sec) = %-17.15g\n", get_hdr_double(hdr, "start_ast"));
    printf("            Start LST (sec) = %-17.15g\n", get_hdr_double(hdr, "start_lst"));
    printf("           Obs Length (sec) = %-17.15g\n", get_hdr_double(hdr, "obs_time"));
    printf("      Requested T_samp (us) = %-17.15g\n", get_hdr_double(hdr, "samp_time"));
    printf("         Actual T_samp (us) = %-17.15g\n", get_hdr_double(hdr, "wapp_time"));
    printf("         Central freq (MHz) = %-17.15g\n", get_hdr_double(hdr, "cent_freq"));
    printf("      Total Bandwidth (MHz) = %-17.15g\n", get_hdr_double(hdr, "bandwidth"));
    printf("             Number of lags = %d\n", get_hdr_int(hdr, "num_lags"));
    printf("              Number of IFs = %d\n", get_hdr_int(hdr, "nifs"));
    printf("    Samples since obs start = %lld\n", get_hdr_longlong(hdr, "timeoff"));
    printf("   Other information:\n");
    if (get_hdr_int(hdr, "sum") == 1) printf("      IFs are summed.\n");
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
    if (decreasing_freqs_st) printf("      Frequency band is inverted.\n");
    if (get_hdr_int(hdr, "lagformat") == 0)
        printf("      Lags are 16 bit integers.\n\n");
    else
        printf("      Lags are 32 bit integers.\n\n");
}

int read_WAPP_rawblock(FILE * infiles[], int numfiles,
                       unsigned char *data, int *padding, IFs ifs)
/* This routine reads a single record from the          */
/* input files *infiles which contain 16 or 32 bit lags */
/* data from the WAPP correlator at Arecibo.            */
/* A WAPP record is ptsperblk_st*numchan_st*#bits long. */
/* *data must be sampperblk_st bytes long.              */
/* If padding is returned as 1, then padding was added  */
/* and statistics should not be calculated.             */
{
   int offset = 0, numtopad = 0, ii, jj;
   unsigned char *dataptr = data;

   /* If our buffer array is offset from last time */
   /* copy the second part into the first.         */

   if (bufferpts) {
      offset = bufferpts * numchan_st;
      dataptr = databuffer + offset;
      if (shiftbuffer)
         memcpy(databuffer, databuffer + sampperblk_st, offset);
   }
   shiftbuffer = 1;

   /* Make sure our current file number is valid */

   if (currentfile >= numfiles / numwapps_st)
      return 0;

   /* First, attempt to read data from the current file */

   if (chkfread(lagbuffer, bytesperblk_st, 1, infiles[currentfile * numwapps_st])) {    /* Got data */
      for (ii = 1; ii < numwapps_st; ii++)      /* Get data from other WAPPs */
         chkfread(lagbuffer + ii * bytesperblk_st, bytesperblk_st, 1,
                  infiles[currentfile * numwapps_st + ii]);
      /* See if we need to byte-swap and if so, doit */
      if (need_byteswap_st) {
         if (bits_per_samp_st == 16) {
            unsigned short *sptr = (unsigned short *) lagbuffer;
            for (ii = 0; ii < sampperblk_st; ii++, sptr++)
               *sptr = swap_ushort(*sptr);
         }
         if (bits_per_samp_st == 32) {
            unsigned int *iptr = (unsigned int *) lagbuffer;
            for (ii = 0; ii < sampperblk_st; ii++, iptr++)
               *iptr = swap_uint(*iptr);
         }
      }
      /* Convert from Correlator Lags to Filterbank Powers */
      for (ii = 0; ii < ptsperblk_st; ii++)
         for (jj = 0; jj < numwapps_st; jj++)
            convert_WAPP_point(lagbuffer + jj * bytesperblk_st + ii * bytesperpt_st,
                               dataptr + jj * numwappchan_st + ii * numchan_st, ifs);

      /* Clip nasty RFI if requested */
      if (clip_sigma_st > 0.0)
         clip_times(dataptr, ptsperblk_st, numchan_st, clip_sigma_st, newpadvals);
      *padding = 0;

      /* Put the new data into the databuffer if needed */
      if (bufferpts) {
         memcpy(data, databuffer, sampperblk_st);
      }
      currentblock++;
      return 1;
   } else {                     /* Didn't get data */
      if (feof(infiles[currentfile * numwapps_st])) {   /* End of file? */
         numtopad = padpts_st[currentfile] - padnum;
         if (numtopad) {        /* Pad the data? */
            *padding = 1;
            if (numtopad >= ptsperblk_st - bufferpts) { /* Lots of padding */
               if (bufferpts) { /* Buffer the padding? */
                  /* Add the amount of padding we need to */
                  /* make our buffer offset = 0           */
                  numtopad = ptsperblk_st - bufferpts;
                  for (ii = 0; ii < numtopad; ii++)
                     memcpy(dataptr + ii * numchan_st, newpadvals, numchan_st);
                  /* Copy the new data/padding into the output array */
                  memcpy(data, databuffer, sampperblk_st);
                  bufferpts = 0;
               } else {         /* Add a full record of padding */
                  numtopad = ptsperblk_st;
                  for (ii = 0; ii < numtopad; ii++)
                     memcpy(data + ii * numchan_st, newpadvals, numchan_st);
               }
               padnum += numtopad;
               currentblock++;
               /* If done with padding reset padding variables */
               if (padnum == padpts_st[currentfile]) {
                  padnum = 0;
                  currentfile++;
               }
               return 1;
            } else {            /* Need < 1 block (or remaining block) of padding */
               int pad;
               /* Add the remainder of the padding and */
               /* then get a block from the next file. */
               for (ii = 0; ii < numtopad; ii++)
                  memcpy(databuffer + bufferpts * numchan_st + ii * numchan_st,
                         newpadvals, numchan_st);
               bufferpts += numtopad;
               padnum = 0;
               shiftbuffer = 0;
               currentfile++;
               return read_WAPP_rawblock(infiles, numfiles, data, &pad, ifs);
            }
         } else {               /* No padding needed.  Try reading the next file */
            currentfile++;
            shiftbuffer = 0;
            return read_WAPP_rawblock(infiles, numfiles, data, padding, ifs);
         }
      } else {
         printf("\nProblem reading record from WAPP data file:\n");
         printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                currentfile, currentblock);
         exit(1);
      }
   }
}


int read_WAPP_rawblocks(FILE * infiles[], int numfiles,
                        unsigned char rawdata[], int numblocks,
                        int *padding, IFs ifs)
/* This routine reads numblocks WAPP records from the input */
/* files *infiles.  The 8-bit filterbank data is returned   */
/* in rawdata which must have a size of numblocks*          */
/* sampperblk_st.  The number  of blocks read is returned.  */
/* If padding is returned as 1, then padding was added      */
/* and statistics should not be calculated                  */
{
   int ii, retval = 0, pad, numpad = 0;

   *padding = 0;
   for (ii = 0; ii < numblocks; ii++) {
      retval += read_WAPP_rawblock(infiles, numfiles,
                                   rawdata + ii * sampperblk_st, &pad, ifs);
      if (pad)
         numpad++;
   }
   /* Return padding 'true' if more than */
   /* half of the blocks are padding.    */
   /* 
      if (numpad > numblocks / 2)
      *padding = 1;
    */
   /* Return padding 'true' if any block was padding */
   if (numpad)
      *padding = 1;
   return retval;
}


int read_WAPP(FILE * infiles[], int numfiles, float *data,
              int numpts, double *dispdelays, int *padding,
              int *maskchans, int *nummasked, mask * obsmask, IFs ifs)
/* This routine reads numpts from the WAPP raw input   */
/* files *infiles.  These files contain raw correlator */
/* data from WAPP at Arecibo.  Time delays             */
/* and a mask are applied to each channel.  It returns */
/* the # of points read if successful, 0 otherwise.    */
/* If padding is returned as 1, then padding was       */
/* added and statistics should not be calculated.      */
/* maskchans is an array of length numchans contains   */
/* a list of the number of channels that were masked.  */
{
   int ii, jj, numread = 0, offset;
   double starttime = 0.0;
   static unsigned char *tempzz, *rawdata1, *rawdata2;
   static unsigned char *currentdata, *lastdata;
   static int firsttime = 1, numblocks = 1, allocd = 0, mask = 0;
   static double duration = 0.0, timeperblk = 0.0;

   *nummasked = 0;
   if (firsttime) {
      if (numpts % ptsperblk_st) {
         printf("numpts must be a multiple of %d in read_WAPP()!\n", ptsperblk_st);
         exit(1);
      } else
         numblocks = numpts / ptsperblk_st;

      if (obsmask->numchan)
         mask = 1;
      rawdata1 = gen_bvect(numblocks * sampperblk_st);
      rawdata2 = gen_bvect(numblocks * sampperblk_st);
      allocd = 1;
      timeperblk = ptsperblk_st * dt_st;
      duration = numblocks * timeperblk;
      currentdata = rawdata1;
      lastdata = rawdata2;
   }

   /* Read and de-disperse */

   if (allocd) {
      while (1) {
         numread = read_WAPP_rawblocks(infiles, numfiles, currentdata,
                                       numblocks, padding, ifs);

         if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
         }

         /* Only use the recently measured padding if all the channels aren't masked */
         if ((clip_sigma_st > 0.0) && 
             !(mask && (*nummasked == -1)) &&
             (padvals != newpadvals))
             memcpy(padvals, newpadvals, WAPP_MAXLAGLEN);
         
         if (mask) {
            if (*nummasked == -1) {     /* If all channels are masked */
               for (ii = 0; ii < numpts; ii++)
                  memcpy(currentdata + ii * numchan_st, padvals, numchan_st);
            } else if (*nummasked > 0) {        /* Only some of the channels are masked */
               int channum;
               for (ii = 0; ii < numpts; ii++) {
                  offset = ii * numchan_st;
                  for (jj = 0; jj < *nummasked; jj++) {
                     channum = maskchans[jj];
                     currentdata[offset + channum] = padvals[channum];
                  }
               }
            }
         }

         if (!firsttime)
            dedisp(currentdata, lastdata, numpts, numchan_st, dispdelays, data);
         SWAP(currentdata, lastdata);
         if (numread != numblocks) {
            vect_free(rawdata1);
            vect_free(rawdata2);
            fftwf_destroy_plan(fftplan);
            fftwf_free(lags);
            vect_free(window_st);
            allocd = 0;
         }
         if (firsttime)
            firsttime = 0;
         else
            break;
      }
      return numblocks * ptsperblk_st;
   } else {
      return 0;
   }
}

void get_WAPP_channel(int channum, float chandat[],
                      unsigned char rawdata[], int numblocks)
/* Return the values for channel 'channum' of a block of       */
/* 'numblocks' raw WAPP data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                */
/* read_WAPP_rawblocks(), and 'chandat' must have at least     */
/* 'numblocks' * 'ptsperblk_st' spaces.                        */
/* Channel 0 is assumed to be the lowest freq channel.         */
{
   int ii, jj, ptsperchan;

   if (channum > numchan_st * numifs_st || channum < 0) {
      printf("\nchannum = %d is out of range in get_WAPP_channel()!\n\n", channum);
      exit(1);
   }
   ptsperchan = ptsperblk_st * numblocks;

   /* Since the following is only called from rfifind, we know that the */
   /* channel accesses will be in order from 0 to the numchan-1         */
   if (channum == 0) {          /* Transpose the data */
      short trtn;
      int move_size;
      unsigned char *move;

      move_size = (ptsperchan + numchan_st) / 2;
      move = gen_bvect(move_size);
      if ((trtn = transpose_bytes(rawdata, ptsperchan, numchan_st,
                                  move, move_size)) < 0)
         printf("Error %d in transpose_bytes().\n", trtn);
      vect_free(move);
   }

   /* Select the correct channel */
   for (ii = 0, jj = ptsperchan * channum; ii < ptsperchan; ii++, jj++)
      chandat[ii] = (float) rawdata[jj];

   /* Select the correct channel */
/*   for (ii=0, jj=channum;  */
/*        ii<numblocks*ptsperblk_st;  */
/*        ii++, jj+=numchan_st) */
/*     chandat[ii] = (float)rawdata[jj]; */
}


int prep_WAPP_subbands(unsigned char *rawdata, float *data,
                       double *dispdelays, int numsubbands,
                       int transpose, int *maskchans, int *nummasked, mask * obsmask)
/* This routine preps a block from the WAPP system.  The routine uses     */
/* dispersion delays in 'dispdelays' to de-disperse the data into         */
/* 'numsubbands' subbands.  It stores the resulting data in vector 'data' */
/* of length 'numsubbands' * 'ptsperblk_st'.  The low freq subband is     */
/* stored first, then the next highest subband etc, with 'ptsperblk_st'   */
/* floating points per subband.  It returns the # of points read if       */
/* succesful, 0 otherwise.  'maskchans' is an array of length numchans    */
/* which contains a list of the number of channels that were masked.  The */
/* # of channels masked is returned in 'nummasked'.  'obsmask' is the     */
/* mask structure to use for masking.  If 'transpose'==0, the data will   */
/* be kept in time order instead of arranged by subband as above.         */
{
   int ii, jj, trtn, offset;
   double starttime = 0.0;
   static unsigned char *tempzz;
   static unsigned char rawdata1[WAPP_MAXDATLEN], rawdata2[WAPP_MAXDATLEN];
   static unsigned char *currentdata, *lastdata, *move;
   static int firsttime = 1, move_size = 0, mask = 0;
   static double timeperblk = 0.0;

   *nummasked = 0;
   if (firsttime) {
      if (obsmask->numchan)
         mask = 1;
      move_size = (ptsperblk_st + numsubbands) / 2;
      move = gen_bvect(move_size);
      currentdata = rawdata1;
      lastdata = rawdata2;
      memcpy(currentdata, rawdata, sampperblk_st);
      timeperblk = ptsperblk_st * dt_st;
   }

   /* Read and de-disperse */

   memcpy(currentdata, rawdata, sampperblk_st);
   if (mask) {
      starttime = currentblock * timeperblk;
      *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
   }

   /* Only use the recently measured padding if all the channels aren't masked */
   if ((clip_sigma_st > 0.0) && 
       !(mask && (*nummasked == -1)) &&
       (padvals != newpadvals))
       memcpy(padvals, newpadvals, WAPP_MAXLAGLEN);
         
   if (mask) {
      if (*nummasked == -1) {   /* If all channels are masked */
         for (ii = 0; ii < ptsperblk_st; ii++)
            memcpy(currentdata + ii * numchan_st, padvals, numchan_st);
      } else if (*nummasked > 0) {      /* Only some of the channels are masked */
         int channum;
         for (ii = 0; ii < ptsperblk_st; ii++) {
            offset = ii * numchan_st;
            for (jj = 0; jj < *nummasked; jj++) {
               channum = maskchans[jj];
               currentdata[offset + channum] = padvals[channum];
            }
         }
      }
   }

   /* In mpiprepsubband, the nodes do not call read_*_rawblock() */
   /* where currentblock gets incremented.                       */
   if (using_MPI) currentblock++;

   if (firsttime) {
      SWAP(currentdata, lastdata);
      firsttime = 0;
      return 0;
   } else {
      dedisp_subbands(currentdata, lastdata, ptsperblk_st, numchan_st,
                      dispdelays, numsubbands, data);
      SWAP(currentdata, lastdata);
      /* Transpose the data into vectors in the result array */
      if (transpose) {
         if ((trtn = transpose_float(data, ptsperblk_st, numsubbands,
                                     move, move_size)) < 0)
            printf("Error %d in transpose_float().\n", trtn);
      }
      return ptsperblk_st;
   }
}


int read_WAPP_subbands(FILE * infiles[], int numfiles, float *data,
                       double *dispdelays, int numsubbands,
                       int transpose, int *padding,
                       int *maskchans, int *nummasked, mask * obsmask, IFs ifs)
/* This routine reads a record from the input files *infiles[]   */
/* which contain data from the WAPP system.  The routine uses    */
/* dispersion delays in 'dispdelays' to de-disperse the data     */
/* into 'numsubbands' subbands.  It stores the resulting data    */
/* in vector 'data' of length 'numsubbands' * 'ptsperblk_st'.    */
/* The low freq subband is stored first, then the next highest   */
/* subband etc, with 'ptsperblk_st' floating points per subband. */
/* It returns the # of points read if succesful, 0 otherwise.    */
/* If padding is returned as 1, then padding was added and       */
/* statistics should not be calculated.  'maskchans' is an array */
/* of length numchans which contains a list of the number of     */
/* channels that were masked.  The # of channels masked is       */
/* returned in 'nummasked'.  'obsmask' is the mask structure     */
/* to use for masking.  If 'transpose'==0, the data will be kept */
/* in time order instead of arranged by subband as above.        */
{
   static int firsttime = 1, memfreed = 0;
   static unsigned char rawdata[WAPP_MAXDATLEN];

   if (firsttime) {
      if (!read_WAPP_rawblock(infiles, numfiles, rawdata, padding, ifs)) {
         printf("Problem reading the raw WAPP data file.\n\n");
         return 0;
      }
      if (0 != prep_WAPP_subbands(rawdata, data, dispdelays, numsubbands,
                                  transpose, maskchans, nummasked, obsmask)) {
         printf("Problem initializing prep_WAPP_subbands()\n\n");
         return 0;
      }
      firsttime = 0;
   }
   if (!read_WAPP_rawblock(infiles, numfiles, rawdata, padding, ifs)) {
      /* printf("Problem reading the raw WAPP data file.\n\n"); */
      if (!memfreed) {
         fftwf_destroy_plan(fftplan);
         fftwf_free(lags);
         free(window_st);
         memfreed = 1;
      }
      return 0;
   }
   return prep_WAPP_subbands(rawdata, data, dispdelays, numsubbands,
                             transpose, maskchans, nummasked, obsmask);
}


void convert_WAPP_point(void *rawdata, unsigned char *bytes, IFs ifs)
/* This routine converts a single point of WAPP lags   */
/* into a filterbank style array of bytes.             */
/* Van Vleck corrections are applied but no window     */
/* functions can be applied as of yet...               */
{
   int ii, ifnum = 0, index = 0;
   float *templags = NULL;
   double power, pfact;
   double scale_min_st = 0.0, scale_max_st = 3.0;

   if (ifs == IF0) {
      ifnum = 1;
      index = 0;
   } else if (ifs == IF1) {
      ifnum = 1;
      index = numwappchan_st;
   } else {                     /* Sum the IFs (or they were already summed) */
      if (numifs_st == 2) {
         scale_min_st *= 2.0;
         scale_max_st *= 2.0;
      }
   }

   /* Loop over the IFs */
   for (ifnum = 0; ifnum < numifs_st; ifnum++, index += numwappchan_st) {

      /* Fill lag array with scaled CFs */
      if (bits_per_samp_st == 16) {
         unsigned short *sdata = (unsigned short *) rawdata;
         for (ii = 0; ii < numwappchan_st; ii++)
            lags[ii] = corr_scale_st * sdata[ii + index] - 1.0;
      } else {
         unsigned int *idata = (unsigned int *) rawdata;
         for (ii = 0; ii < numwappchan_st; ii++)
            lags[ii] = corr_scale_st * idata[ii + index] - 1.0;
      }

      /* Calculate power */
      power = inv_cerf(lags[0]);
      power = 0.1872721836 / (power * power);

      /* Apply Van Vleck Corrections to the Lags */
      if (corr_level_st == 3)
         vanvleck3lev(lags, numwappchan_st);
      else if (corr_level_st == 9)
         vanvleck9lev(lags, numwappchan_st);
      else
         printf("\nError:  corr_level_st (%d) does not equal 3 or 9!\n\n",
                corr_level_st);

      for (ii = 0; ii < numwappchan_st; ii++)
         lags[ii] *= power;

      if (usewindow_st)
         for (ii = 0; ii < numwappchan_st; ii++)
            lags[ii] *= window_st[ii];

      /* FFT the ACF lags (which are real and even) -> real and even FFT */
      lags[numwappchan_st] = 0.0;
      fftwf_execute(fftplan);

#if 0
      printf("\n");
      for (ii = 0; ii < numwappchan_st; ii++)
         printf("%d  %.7g\n", ii, lags[ii]);
      printf("\n");
      exit(0);
#endif

      /* Reverse band if it needs it */
      if (decreasing_freqs_st) {
         float tempzz = 0.0, *loptr, *hiptr;
         loptr = lags + 0;
         hiptr = lags + numwappchan_st - 1;
         for (ii = 0; ii < numwappchan_st / 2; ii++, loptr++, hiptr--) {
            SWAP(*loptr, *hiptr);
         }
      }

      if (numifs_st == 2 && ifs == SUMIFS) {
         if (ifnum == 0) {
            templags = gen_fvect(numwappchan_st);
            /* Copy the unscaled values to the templag array */
            for (ii = 0; ii < numwappchan_st; ii++)
               templags[ii] = lags[ii];
         } else {
            /* Sum the unscaled IFs */
            for (ii = 0; ii < numwappchan_st; ii++)
               lags[ii] += templags[ii];
            vect_free(templags);
         }
      }
   }
   /* Scale and pack the powers */
   pfact = 255.0 / (scale_max_st - scale_min_st);
   for (ii = 0; ii < numwappchan_st; ii++) {
      double templag;
      templag = (lags[ii] > scale_max_st) ? scale_max_st : lags[ii];
      templag = (templag < scale_min_st) ? scale_min_st : templag;
      bytes[ii] = (unsigned char) ((templag - scale_min_st) * pfact + 0.5);
   }

#if 0
   {                            /* Show what the raw powers are (avg ~1.05, var ~0.2) */
      double avg, var;
      avg_var(lags, numwappchan_st, &avg, &var);
      printf("avg = %f    var = %f\n", avg, var);
      exit(0);
   }
#endif
}


static double inv_cerf(double input)
/* Approximation for Inverse Complementary Error Function */
{
   static double numerator_const[3] = {
      1.591863138, -2.442326820, 0.37153461
   };
   static double denominator_const[3] = {
      1.467751692, -3.013136362, 1.0
   };
   double num, denom, temp_data, temp_data_srq, erf_data;

   erf_data = 1.0 - input;
   temp_data = erf_data * erf_data - 0.5625;
   temp_data_srq = temp_data * temp_data;
   num = erf_data * (numerator_const[0] +
                     (temp_data * numerator_const[1]) +
                     (temp_data_srq * numerator_const[2]));
   denom = denominator_const[0] + temp_data * denominator_const[1] +
       temp_data_srq * denominator_const[2];
   return num / denom;
}


#define NO    0
#define YES   1
/*------------------------------------------------------------------------*
 * Van Vleck Correction for 9-level sampling/correlation
 *  Samples {-4,-3,-2,-1,0,1,2,3,4}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into ACF of at least 'count' points
 * This routine takes the first value as the zerolag and corrects the
 * remaining count-1 points.  Zerolag is set to a normalized 1
 * NOTE - The available routine works on lags normaized to -16<rho<16, so
 *  I need to adjust the values before/after the fit
 * Coefficent ranges
 *   c1
 *     all
 *   c2
 *     r0 > 4.5
 *     r0 < 2.1
 *     rest
 * NOTE - correction is done INPLACE ! Original values are destroyed
 * As reported by M. Lewis -> polynomial fits are OK, but could be improved
 *------------------------------------------------------------------------*/
static void vanvleck9lev(float *rho, int npts)
{
   double acoef[5], dtmp, zl;
   int i;
   static double coef1[5] =
       { 1.105842267, -0.053258115, 0.011830276, -0.000916417, 0.000033479 };
   static double coef2rg4p5[5] =
       { 0.111705575, -0.066425925, 0.014844439, -0.001369796, 0.000044119 };
   static double coef2rl2p1[5] =
       { 1.285303775, -1.472216011, 0.640885537, -0.123486209, 0.008817175 };
   static double coef2rother[5] =
       { 0.519701391, -0.451046837, 0.149153116, -0.021957940, 0.001212970 };
   static double coef3rg2p0[5] =
       { 1.244495105, -0.274900651, 0.022660239, -0.000760938, -1.993790548 };
   static double coef3rother[5] =
       { 1.249032787, 0.101951346, -0.126743165, 0.015221707, -2.625961708 };
   static double coef4rg3p15[5] =
       { 0.664003237, -0.403651682, 0.093057131, -0.008831547, 0.000291295 };
   static double coef4rother[5] =
       { 9.866677289, -12.858153787, 6.556692205, -1.519871179, 0.133591758 };
   static double coef5rg4p0[4] =
       { 0.033076469, -0.020621902, 0.001428681, 0.000033733 };
   static double coef5rg2p2[4] =
       { 5.284269565, 6.571535249, -2.897741312, 0.443156543 };
   static double coef5rother[4] =
       { -1.475903733, 1.158114934, -0.311659264, 0.028185170 };

   zl = rho[0] * 16;
   /* ro = *rho;                 */
   /*  for(i=0; i<npts; i++)     */
   /*    (rho+i) *= *(rho+i)/ro; */
   acoef[0] =
       ((((coef1[4] * zl + coef1[3]) * zl + coef1[2]) * zl +
         coef1[1]) * zl + coef1[0]);
   if (zl > 4.50)
      acoef[1] =
          ((((coef2rg4p5[4] * zl + coef2rg4p5[3]) * zl +
             coef2rg4p5[2]) * zl + coef2rg4p5[1]) * zl + coef2rg4p5[0]);
   else if (zl < 2.10)
      acoef[1] =
          ((((coef2rl2p1[4] * zl + coef2rl2p1[3]) * zl +
             coef2rl2p1[2]) * zl + coef2rl2p1[1]) * zl + coef2rl2p1[0]);
   else
      acoef[1] =
          ((((coef2rother[4] * zl + coef2rother[3]) * zl +
             coef2rother[2]) * zl + coef2rother[1]) * zl + coef2rother[0]);
   if (zl > 2.00)
      acoef[2] =
          coef3rg2p0[4] / zl +
          (((coef3rg2p0[3] * zl + coef3rg2p0[2]) * zl +
            coef3rg2p0[1]) * zl + coef3rg2p0[0]);
   else
      acoef[2] =
          coef3rother[4] / zl +
          (((coef3rother[3] * zl + coef3rother[2]) * zl +
            coef3rother[1]) * zl + coef3rother[0]);
   if (zl > 3.15)
      acoef[3] =
          ((((coef4rg3p15[4] * zl + coef4rg3p15[3]) * zl +
             coef4rg3p15[2]) * zl + coef4rg3p15[1]) * zl + coef4rg3p15[0]);
   else
      acoef[3] =
          ((((coef4rg3p15[4] * zl + coef4rother[3]) * zl +
             coef4rother[2]) * zl + coef4rother[1]) * zl + coef4rother[0]);
   if (zl > 4.00)
      acoef[4] =
          (((coef5rg4p0[3] * zl + coef5rg4p0[2]) * zl +
            coef5rg4p0[1]) * zl + coef5rg4p0[0]);
   else if (zl < 2.2)
      acoef[4] =
          (((coef5rg2p2[3] * zl + coef5rg2p2[2]) * zl +
            coef5rg2p2[1]) * zl + coef5rg2p2[0]);
   else
      acoef[4] =
          (((coef5rother[3] * zl + coef5rother[2]) * zl +
            coef5rother[1]) * zl + coef5rother[0]);
   for (i = 1; i < npts; i++) {
      dtmp = rho[i];
      rho[i] =
          ((((acoef[4] * dtmp + acoef[3]) * dtmp + acoef[2]) * dtmp +
            acoef[1]) * dtmp + acoef[0]) * dtmp;
   }
   rho[0] = 1.0;
   return;
}

/*------------------------------------------------------------------------*
 * Van Vleck Correction for 3-level sampling/correlation
 *  Samples {-1,0,1}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into ACF of at least 'count' points
 * This routine takes the first value as the zerolag and corrects the
 * remaining count-1 points.  Zerolag is set to a normalized 1
 *
 * NOTE - correction is done INPLACE ! Original values are destroyed
 *------------------------------------------------------------------------*/
static void vanvleck3lev(float *rho, int npts)
{
   double lo_u[3], lo_h[3];
   double high_u[5], high_h[5];
   double lo_coefficient[3];
   double high_coefficient[5];
   double zho, zho_3;
   double temp_data;
   double temp_data_1;
   int ichan, ico, flag_any_high;
   static double lo_const[3][4] = {
      {0.939134371719, -0.567722496249, 1.02542540932, 0.130740914912},
      {-0.369374472755, -0.430065136734, -0.06309459132, -0.00253019992917},
      {0.888607422108, -0.230608118885, 0.0586846424223, 0.002012775510695}
   };
   static double high_const[5][4] = {
      {-1.83332160595, 0.719551585882, 1.214003774444, 7.15276068378e-5},
      {1.28629698818, -1.45854382672, -0.239102591283, -0.00555197725185},
      {-7.93388279993, 1.91497870485, 0.351469403030, 0.00224706453982},
      {8.04241371651, -1.51590759772, -0.18532022393, -0.00342644824947},
      {-13.076435520, 0.769752851477, 0.396594438775, 0.0164354218208}
   };

   /* Perform Lo correction on All data that is not flaged 
      for high correction  */
   zho = (double) rho[0];
   zho_3 = zho * zho * zho;
   lo_u[0] = zho;
   lo_u[1] = zho_3 - (61.0 / 512.0);
   lo_u[2] = zho - (63.0 / 128.0);
   lo_h[0] = zho * zho;
   lo_h[2] = zho_3 * zho_3 * zho;       /* zlag ^7 */
   lo_h[1] = zho * lo_h[2];     /* zlag ^8 */
   /* determine lo-correct coefficents - */
   for (ico = 0; ico < 3; ico++) {
      lo_coefficient[ico] =
          (lo_u[ico] *
           (lo_u[ico] *
            (lo_u[ico] * lo_const[ico][0] + lo_const[ico][1]) +
            lo_const[ico][2]) + lo_const[ico][3]) / lo_h[ico];
   }
   /* perform correction -- */
   for (ichan = 1, flag_any_high = NO; ichan < npts; ichan++) {
      temp_data = (double) rho[ichan];
      if (fabs(temp_data) > 0.199) {
         if (flag_any_high == NO) {
            high_u[0] = lo_h[2];        /* zlag ^7 */
            high_u[1] = zho - (63.0 / 128.0);
            high_u[2] = zho * zho - (31.0 / 128.0);
            high_u[3] = zho_3 - (61.0 / 512.0);
            high_u[4] = zho - (63.0 / 128.0);
            high_h[0] = lo_h[1];        /* zlag ^8 */
            high_h[1] = lo_h[1];        /* zlag ^8 */
            high_h[2] = lo_h[1] * zho_3 * zho;  /* zlag ^12 */
            high_h[3] = lo_h[1] * lo_h[1] * zho;        /* zlag ^17 */
            high_h[4] = high_h[3];      /* zlag ^17 */
            for (ico = 0; ico < 5; ico++) {
               high_coefficient[ico] =
                   (high_u[ico] *
                    (high_u[ico] *
                     (high_u[ico] * high_const[ico][0] +
                      high_const[ico][1]) + high_const[ico][2]) +
                    high_const[ico][3]) / high_h[ico];
            }
            flag_any_high = YES;
         }
         temp_data_1 = fabs(temp_data * temp_data * temp_data);
         rho[ichan] =
             (temp_data *
              (temp_data_1 *
               (temp_data_1 *
                (temp_data_1 *
                 (temp_data_1 * high_coefficient[4] +
                  high_coefficient[3]) + high_coefficient[2]) +
                high_coefficient[1]) + high_coefficient[0]));
      } else {
         temp_data_1 = temp_data * temp_data;
         rho[ichan] =
             (temp_data *
              (temp_data_1 *
               (temp_data_1 * lo_coefficient[2] + lo_coefficient[1]) +
               lo_coefficient[0]));
      }
   }
   rho[0] = 1.0;
}
