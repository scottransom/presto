#include "presto.h"
#include "mask.h"
#include "bpp.h"

/* All of the following have an _st to indicate static */
static long long numpts_st[MAXPATCHFILES], padpts_st[MAXPATCHFILES], N_st;
static long long filedatalen_st[MAXPATCHFILES];
static int numblks_st[MAXPATCHFILES];
static int bytesperpt_st, bytesperblk_st;
static int numchan_st, numifs_st, ptsperblk_st = PTSPERBLOCK;
static double times_st[MAXPATCHFILES], mjds_st[MAXPATCHFILES];
static double elapsed_st[MAXPATCHFILES], T_st, dt_st;
static double startblk_st[MAXPATCHFILES], endblk_st[MAXPATCHFILES];
static infodata idata_st[MAXPATCHFILES];
//static float padvals[MAXNUMCHAN], padval = 4;
//static unsigned char databuffer[2 * MAXDATLEN];
static unsigned char *splitbytes_buffer[MAXPATCHFILES];
static int currentfile, currentblock, both_IFs_present = 0;
//static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
/* To fix the fact that the BCPM sometimes splits a sample over 2 files. */
static int splitbytes_st[MAXPATCHFILES];
static double mid_freq_st, ch1_freq_st, delta_freq_st;
static double chan_freqs[2 * MAXNUMCHAN];
static int chan_index[2 * MAXNUMCHAN], chan_mapping[2 * MAXNUMCHAN];
static float clip_sigma_st = 0.0;
//static int using_MPI = 0;

#if 0
static double chan_freqs2[2 * MAXNUMCHAN];
static int chan_index2[2 * MAXNUMCHAN];
#endif
static int dfb_chan_lookup[MAXREGS][NIBPERREG] = {
    {4, 0, 4, 0},
    {5, 1, 5, 1},
    {6, 2, 6, 2},
    {7, 3, 7, 3},
    {4, 0, 4, 0},
    {5, 1, 5, 1},
    {6, 2, 6, 2},
    {7, 3, 7, 3}
};

/* This takes care of byte swap in outreg_b */
static float sideband_lookup[MAXREGS][NIBPERREG] = {
    {-1.0, -1.0, +1.0, +1.0},
    {-1.0, -1.0, +1.0, +1.0},
    {-1.0, -1.0, +1.0, +1.0},
    {-1.0, -1.0, +1.0, +1.0},
    {+1.0, +1.0, -1.0, -1.0},
    {+1.0, +1.0, -1.0, -1.0},
    {+1.0, +1.0, -1.0, -1.0},
    {+1.0, +1.0, -1.0, -1.0}
};

typedef struct chanmap {
    double freq;
    int index;
    int mapping;
} chanmap;

double slaCldj(int iy, int im, int id, int *j);


int compare_freq(const void *ca, const void *cb)
/* qsort comparison function for chanmap */
{
    chanmap *a, *b;

    a = (chanmap *) ca;
    b = (chanmap *) cb;
    if ((b->freq - a->freq) < 0.0)
        return 1;
    if ((b->freq - a->freq) > 0.0)
        return -1;
    /* This keeps IF0 freqs before IF1 freqs */
    /*
       if ((b->index - a->index) < 0)
       return 1;
       if ((b->index - a->index) > 0)
       return -1;
     */
    return 0;
}


int compare_index(const void *ca, const void *cb)
/* qsort comparison function for chanmap */
{
    chanmap *a, *b;

    a = (chanmap *) ca;
    b = (chanmap *) cb;
    if ((b->index - a->index) < 0)
        return 1;
    if ((b->index - a->index) > 0)
        return -1;
    return 0;
}


static double UT_strings_to_MJD(char *date, char *start_time,
                                int *mjd_day, double *mjd_fracday)
{
    int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    int julday, year, month = 1, day, hour, min, sec, err;

    sscanf(date, "%3d:%4d", &julday, &year);
    julday++;                   /* UT date starts from 0, I believe */
    /* Allow for leap years */
    if (year % 4 == 0)
        mtab[1] = 29;
    if (year % 100 == 0 && year % 400 != 0)
        mtab[1] = 28;
    /* Convert Julian day to day and month */
    while (julday - mtab[month - 1] > 0) {
        julday -= mtab[month - 1];
        month++;
    }
    day = julday;
    sscanf(start_time, "%2d:%2d:%2d", &hour, &min, &sec);
    /* Recorded BCPM start time is 1 second late */
    *mjd_fracday = (hour + (min + ((sec + 1) / 60.0)) / 60.0) / 24.0;
    *mjd_day = slaCldj(year, month, day, &err);
    return *mjd_day + *mjd_fracday;
}


static void swapendian_BPP_header(BPP_SEARCH_HEADER * hdr)
/* This is required since it is a binary header */
{
    int ii;

    hdr->header_version = swap_int(hdr->header_version);
    hdr->bit_mode = swap_int(hdr->bit_mode);
    hdr->num_chans = swap_int(hdr->num_chans);
    hdr->lmst = swap_int(hdr->lmst);
    hdr->scan_file_number = swap_int(hdr->scan_file_number);
    hdr->file_size = swap_int(hdr->file_size);
    hdr->tape_num = swap_int(hdr->tape_num);
    hdr->tape_file_number = swap_int(hdr->tape_file_number);
    hdr->enabled_CBs = swap_int(hdr->enabled_CBs);
    hdr->mb_start_address = swap_int(hdr->mb_start_address);
    hdr->mb_end_address = swap_int(hdr->mb_end_address);
    hdr->mb_start_board = swap_int(hdr->mb_start_board);
    hdr->mb_end_board = swap_int(hdr->mb_end_board);
    hdr->mb_vme_mid_address = swap_int(hdr->mb_vme_mid_address);
    hdr->mb_ack_enabled = swap_int(hdr->mb_ack_enabled);
    hdr->start_from_ste = swap_int(hdr->start_from_ste);
    hdr->cb_sum_polarizations = swap_int(hdr->cb_sum_polarizations);
    hdr->cb_direct_mode = swap_int(hdr->cb_direct_mode);
    hdr->cb_accum_length = swap_int(hdr->cb_accum_length);
    hdr->tb_outs_reg = swap_int(hdr->tb_outs_reg);
    hdr->tb_ste = swap_int(hdr->tb_ste);
    hdr->tb_stc = swap_int(hdr->tb_stc);
    hdr->H_deci_factor = swap_int(hdr->H_deci_factor);
    hdr->GenStat0 = swap_int(hdr->GenStat0);
    hdr->GenStat1 = swap_int(hdr->GenStat1);
    hdr->Ack_Reg = swap_int(hdr->Ack_Reg);
    hdr->dfb_sram_length = swap_int(hdr->dfb_sram_length);
    hdr->ASYMMETRIC = swap_int(hdr->ASYMMETRIC);
    hdr->mb_long_ds0 = swap_int(hdr->mb_long_ds0);
    hdr->aib_serial = swap_int(hdr->aib_serial);
    hdr->aib_rev = swap_int(hdr->aib_rev);
    hdr->BACKEND_TYPE = swap_int(hdr->BACKEND_TYPE);
    hdr->UPDATE_DONE = swap_int(hdr->UPDATE_DONE);
    hdr->HEADER_TYPE = swap_int(hdr->HEADER_TYPE);
    hdr->tb_id = swap_int(hdr->tb_id);
    hdr->aib_if_switch = swap_int(hdr->aib_if_switch);
    hdr->mb_rev = swap_int(hdr->mb_rev);
    hdr->mb_serial = swap_int(hdr->mb_serial);
    hdr->tb_rev = swap_int(hdr->tb_rev);
    hdr->tb_serial = swap_int(hdr->tb_serial);
    hdr->mb_xtal_freq = swap_int(hdr->mb_xtal_freq);
    hdr->scan_num = swap_uint(hdr->scan_num);
    hdr->ll_file_offset = swap_longlong(hdr->ll_file_offset);
    hdr->ll_file_size = swap_longlong(hdr->ll_file_size);
    hdr->length_of_integration = swap_double(hdr->length_of_integration);
    hdr->samp_rate = swap_double(hdr->samp_rate);
    hdr->ra_2000 = swap_double(hdr->ra_2000);
    hdr->dec_2000 = swap_double(hdr->dec_2000);
    hdr->tele_x = swap_double(hdr->tele_x);
    hdr->tele_y = swap_double(hdr->tele_y);
    hdr->tele_z = swap_double(hdr->tele_z);
    hdr->tele_inc = swap_double(hdr->tele_inc);
    hdr->Fclk = swap_double(hdr->Fclk);
    hdr->Har_Clk = swap_double(hdr->Har_Clk);
    hdr->bandwidth = swap_double(hdr->bandwidth);
    hdr->rf_lo = swap_double(hdr->rf_lo);
    hdr->max_dfb_freq = swap_double(hdr->max_dfb_freq);
    hdr->mjd_start = swap_longdouble(hdr->mjd_start);
    for (ii = 0; ii < FB_CHAN_PER_BRD; ii++) {
        hdr->dfb_sram_freqs[ii] = swap_float(hdr->dfb_sram_freqs[ii]);
    }
    for (ii = 0; ii < MAX_HARRIS_TAPS; ii++) {
        hdr->i_hcoef[ii] = swap_int(hdr->i_hcoef[ii]);
        hdr->q_hcoef[ii] = swap_int(hdr->q_hcoef[ii]);
    }
    for (ii = 0; ii < MAXNUMCB; ii++) {
        hdr->cb_id[ii] = swap_int(hdr->cb_id[ii]);
        hdr->cb_rev[ii] = swap_int(hdr->cb_rev[ii]);
        hdr->cb_serial[ii] = swap_int(hdr->cb_serial[ii]);
        hdr->cb_eprom_mode[ii] = swap_int(hdr->cb_eprom_mode[ii]);
    }
    for (ii = 0; ii < MAX_NUM_MF_BOARDS; ii++) {
        hdr->mf_rev[ii] = swap_int(hdr->mf_rev[ii]);
        hdr->mf_serial[ii] = swap_int(hdr->mf_serial[ii]);
        hdr->mf_filt_width[ii] = swap_double(hdr->mf_filt_width[ii]);
        hdr->mf_atten[ii] = swap_double(hdr->mf_atten[ii]);
    }
    for (ii = 0; ii < MAX_NUM_LO_BOARDS; ii++) {
        hdr->lo_rev[ii] = swap_int(hdr->lo_rev[ii]);
        hdr->lo_serial[ii] = swap_int(hdr->lo_serial[ii]);
        hdr->aib_los[ii] = swap_double(hdr->aib_los[ii]);
    }
    for (ii = 0; ii < MAXNUMDFB; ii++) {
        hdr->dfb_mixer_reg[ii] = swap_int(hdr->dfb_mixer_reg[ii]);
        hdr->dfb_conf_reg[ii] = swap_int(hdr->dfb_conf_reg[ii]);
        hdr->dfb_sram_addr_msb[ii] = swap_int(hdr->dfb_sram_addr_msb[ii]);
        hdr->dfb_id[ii] = swap_int(hdr->dfb_id[ii]);
        hdr->dfb_rev[ii] = swap_int(hdr->dfb_rev[ii]);
        hdr->dfb_serial[ii] = swap_int(hdr->dfb_serial[ii]);
        hdr->dfb_sun_program[ii] = swap_int(hdr->dfb_sun_program[ii]);
        hdr->dfb_eprom[ii] = swap_int(hdr->dfb_eprom[ii]);
        hdr->dfb_sram_addr[ii] = swap_int(hdr->dfb_sram_addr[ii]);
        hdr->dfb_har_addr[ii] = swap_int(hdr->dfb_har_addr[ii]);
        hdr->dfb_clip_adc_neg8[ii] = swap_int(hdr->dfb_clip_adc_neg8[ii]);
        hdr->dfb_shften_[ii] = swap_int(hdr->dfb_shften_[ii]);
        hdr->dfb_fwd_[ii] = swap_int(hdr->dfb_fwd_[ii]);
        hdr->dfb_rvrs_[ii] = swap_int(hdr->dfb_rvrs_[ii]);
        hdr->dfb_asymmetric[ii] = swap_int(hdr->dfb_asymmetric[ii]);
        hdr->dfb_i_dc[ii] = swap_double(hdr->dfb_i_dc[ii]);
        hdr->dfb_q_dc[ii] = swap_double(hdr->dfb_q_dc[ii]);
        hdr->dfb_gain[ii] = swap_double(hdr->dfb_gain[ii]);
    }
}


void calc_BPP_chans(BPP_SEARCH_HEADER * hdr)
/* Calculates freqs and ordering index table for BPP channels */
{
    int ii, n = 0, dfb_chan, logical_board, regid, bid, nibble, nchans;
    double f_aib, u_or_l, f_sram, fc;
    chanmap *chanmaparr;

    /* The following is probably a bad way to see if */
    /* we need to swap the endianess of the header.  */
    if (hdr->num_chans < 0 || hdr->num_chans > 2 * MAXNUMCHAN) {
        swapendian_BPP_header(hdr);
    }
    nchans = (hdr->mb_end_address / 2 - hdr->mb_start_address / 2 + 1) *
        (hdr->mb_end_board - hdr->mb_start_board + 1) * 4;
    if (nchans > 2 * MAXNUMCHAN) {
        printf("Error:  nchans (%d) > 2*MAXNUMCHAN (%d) in bpp_calc_chans()\n\n",
               nchans, 2 * MAXNUMCHAN);
        exit(1);
    } else if (nchans == 2 * MAXNUMCHAN) {
        both_IFs_present = 1;
    }
    if (fabs(hdr->rf_lo) >= 1.e6)
        hdr->rf_lo /= 1.e6;

    /* Loop over (16-bit) regs per board. divide by 2's are to make them   */
    /* word addresses instead of byte addresses so we can index with them. */
    /* Normal modes will be regid = 0..3, 0..7, or 4..7                    */
    for (regid = hdr->mb_start_address / 2; regid <= hdr->mb_end_address / 2;
         regid++) {
        /* Loop over each board */
        for (bid = hdr->mb_start_board; bid <= hdr->mb_end_board; bid++) {
            /* Now find which LOGICAL CB we are reading */
            logical_board = -1;
            for (ii = 0; ii < MAXNUMCB; ii++) {
                if (bid == hdr->cb_id[ii]) {
                    logical_board = ii;
                    break;
                }
            }
            if (logical_board == -1) {
                printf("calc_BPP_chans() - logical_board not found");
                exit(1);
            }
            /* Assumes cabling so that LO0 feeds MF0,1 which feeds leftmost CB! */
            f_aib = hdr->aib_los[logical_board];
            /* Loop over 4 nibbles per reg */
            for (nibble = 0; nibble < 4; nibble++) {
                dfb_chan = dfb_chan_lookup[regid][nibble];
                u_or_l = sideband_lookup[regid][nibble];
                f_sram = hdr->dfb_sram_freqs[dfb_chan];
                fc = f_aib + f_sram + u_or_l * hdr->bandwidth / 4.0;
                /* obs below 10 GHz are LSB; above 10 GHz are USB */
                if (hdr->rf_lo < 1.e4)
                    chan_freqs[n++] = hdr->rf_lo + 800 - fc / 1.0e6;
                else
                    chan_freqs[n++] = hdr->rf_lo + fc / 1.0e6;
                /* The following is a hack in order to get 350MHz working with the right freqs */
                {
                    char *envval = getenv("BCPM_FREQ_ADJ");
                    if (envval != NULL) {
                        double dblval = strtod(envval, NULL);
                        if (dblval) {
                            /* For 350MHz, this dblval should be -1470 ; */
                            chan_freqs[n - 1] += dblval;
                        }
                    }
                }
            }
        }
    }

    /* Make a lookup table which gives chans in order of increasing freq */

    numchan_st = MAXNUMCHAN;
    chanmaparr = (chanmap *) malloc(sizeof(chanmap) * nchans);
    for (ii = 0; ii < nchans; ii++) {
        chanmaparr[ii].freq = chan_freqs[ii];
        chanmaparr[ii].index = ii;
    }
    /* Sort IF1 or summed IFs according to freq */
    qsort(chanmaparr, MAXNUMCHAN, sizeof(chanmap), compare_freq);
    for (ii = 0; ii < MAXNUMCHAN; ii++) {
        chan_index[ii] = chanmaparr[ii].index;
        chanmaparr[ii].mapping = ii;
    }
    /* Sort IF1 or summed IFs according to index */
    qsort(chanmaparr, MAXNUMCHAN, sizeof(chanmap), compare_index);

    if (both_IFs_present) {
        /* Sort IF2 according to freq */
        qsort(chanmaparr + MAXNUMCHAN, MAXNUMCHAN, sizeof(chanmap), compare_freq);
        for (ii = MAXNUMCHAN; ii < 2 * MAXNUMCHAN; ii++) {
            chan_index[ii] = chanmaparr[ii].index;
            chanmaparr[ii].mapping = ii - MAXNUMCHAN;
        }
        /* Sort IF2 according to index */
        qsort(chanmaparr + MAXNUMCHAN, MAXNUMCHAN, sizeof(chanmap), compare_index);
    }
    for (ii = 0; ii < nchans; ii++)
        chan_mapping[ii] = chanmaparr[ii].mapping;
    /* 
       for (ii=0; ii<nchans; ii++){
       chan_mapping[ii] = chanmaparr[ii].mapping;
       printf("%3d %3d %3d %10.4f %10.4f\n", ii, chan_index[ii], 
       chan_mapping[ii], chan_freqs[chan_index[ii]], chan_freqs[ii]);
       }
     */
    free(chanmaparr);

    /* Set the static variables */

    n = nchans / 2;
    mid_freq_st = 0.5 * (chan_freqs[chan_index[n]] + chan_freqs[chan_index[n - 1]]);
    ch1_freq_st = chan_freqs[chan_index[0]];
    delta_freq_st = chan_freqs[chan_index[1]] - chan_freqs[chan_index[0]];
}


void BPP_hdr_to_inf(BPP_SEARCH_HEADER * hdr, infodata * idata)
/* Convert BPP header into an infodata structure */
{
    double MJD;
    char ctmp[80];

    strncpy(idata->object, hdr->target_name, 32);
    idata->ra_h = (int) floor(hdr->ra_2000 / 10000.0);
    idata->ra_m = (int) floor((hdr->ra_2000 - idata->ra_h * 10000) / 100.0);
    idata->ra_s = hdr->ra_2000 - idata->ra_h * 10000 - idata->ra_m * 100;
    idata->dec_d = (int) floor(fabs(hdr->dec_2000) / 10000.0);
    idata->dec_m = (int) floor((fabs(hdr->dec_2000) - idata->dec_d * 10000) / 100.0);
    idata->dec_s = fabs(hdr->dec_2000) - idata->dec_d * 10000 - idata->dec_m * 100;
    if (hdr->dec_2000 < 0.0)
        idata->dec_d = -idata->dec_d;
    strcpy(idata->telescope, "GBT");
    strcpy(idata->instrument, "BCPM1");
    idata->num_chan = numchan_st;
    idata->dt = hdr->samp_rate / 1000000.0;
    idata->N = hdr->length_of_integration / idata->dt;
    idata->chan_wid = fabs(delta_freq_st);
    idata->freqband = idata->num_chan * idata->chan_wid;
    idata->freq = ch1_freq_st;
    idata->fov = 1.2 * SOL / (mid_freq_st * 1e6) / 100.0 * RADTODEG * 3600.0;
    MJD = UT_strings_to_MJD(hdr->date, hdr->start_time,
                            &(idata->mjd_i), &(idata->mjd_f));
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    strcpy(idata->analyzer, "Scott Ransom");
    strcpy(idata->observer, "--");
    if (hdr->cb_sum_polarizations)
        sprintf(ctmp, "Polarizations were summed in hardware.");
    else
        sprintf(ctmp, "Polarizations were not summed in hardware.");
    sprintf(idata->notes,
            "Scan number %010u from tape %d, file %d.\n    Topo UTC Date (DDD:YYYY) & Time at file start = %s, %s\n    %s\n",
            hdr->scan_num, hdr->tape_num, hdr->tape_file_number, hdr->date,
            hdr->start_time, ctmp);
}


void get_BPP_file_info(FILE * files[], int numfiles, float clipsig,
                       long long *N, int *ptsperblock, int *numchan, double *dt,
                       double *T, infodata * idata, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of BPP rawfiles that you want to patch      */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
    int ii, newscan = 0;
    char rawhdr[BPP_HEADER_SIZE];
    double last_file_epoch = 0.0;
    BPP_SEARCH_HEADER *header;

    if (numfiles > MAXPATCHFILES) {
        printf("\nThe number of input files (%d) is greater than \n", numfiles);
        printf("   MAXPATCHFILES=%d.  Exiting.\n\n", MAXPATCHFILES);
        exit(0);
    }
    chkfread(rawhdr, BPP_HEADER_SIZE, 1, files[0]);
    header = (BPP_SEARCH_HEADER *) rawhdr;
    calc_BPP_chans(header);
    BPP_hdr_to_inf(header, &idata_st[0]);
    BPP_hdr_to_inf(header, idata);
    /* Are we going to clip the data? */
    if (clipsig > 0.0)
        clip_sigma_st = clipsig;
    *numchan = numchan_st = idata_st[0].num_chan;
    *ptsperblock = ptsperblk_st;
    if (both_IFs_present) {
        printf("  (Note:  Both IFs are present.)\n");
        numifs_st = 2;
    } else
        numifs_st = 1;
    /* The following is the number of bytes in the _raw_ data */
    bytesperpt_st = (numchan_st * numifs_st * 4) / 8;
    bytesperblk_st = ptsperblk_st * bytesperpt_st;
    filedatalen_st[0] = chkfilelen(files[0], 1) - BPP_HEADER_SIZE;
    numblks_st[0] = filedatalen_st[0] / bytesperblk_st;
    splitbytes_st[0] = filedatalen_st[0] % bytesperblk_st;
    if (splitbytes_st[0]) {
        splitbytes_st[0] = bytesperblk_st - filedatalen_st[0] % bytesperblk_st;
        if (numfiles > 1) {
            printf
                ("  File %2d has a non-integer number of complete blocks and/or samples!\n"
                 "\tApplying work-around.  (bytes split = %d)\n", 1,
                 splitbytes_st[0]);
            numblks_st[0]++;
            /* This makes a memory leak as it is never freed (although the OS */
            /* does free it when the program exits -- which is usually when   */
            /* we want it freed anyways...)                                   */
            splitbytes_buffer[0] = (unsigned char *) malloc(splitbytes_st[0]);
        }
    }
    numpts_st[0] = numblks_st[0] * ptsperblk_st;
    N_st = numpts_st[0];
    dt_st = *dt = idata_st[0].dt;
    times_st[0] = numpts_st[0] * dt_st;
    mjds_st[0] = idata_st[0].mjd_i + idata_st[0].mjd_f;
    last_file_epoch = mjds_st[0];
    elapsed_st[0] = 0.0;
    startblk_st[0] = 1;
    endblk_st[0] = (double) numpts_st[0] / ptsperblk_st;
    padpts_st[0] = padpts_st[numfiles - 1] = 0;
    for (ii = 1; ii < numfiles; ii++) {
        chkfread(rawhdr, BPP_HEADER_SIZE, 1, files[ii]);
        header = (BPP_SEARCH_HEADER *) rawhdr;
        calc_BPP_chans(header);
        BPP_hdr_to_inf(header, &idata_st[ii]);
        if (idata_st[ii].num_chan != numchan_st) {
            printf("Number of channels (file %d) is not the same!\n\n", ii + 1);
        }
        if (idata_st[ii].dt != dt_st) {
            printf("Sample time (file %d) is not the same!\n\n", ii + 1);
        }
        /* If the MJDs are equal, then this is a continuation    */
        /* of the same scan.  In that case, calculate the _real_ */
        /* duration of the previous file and add it to the       */
        /* previous files MJD to get the current MJD.            */
        mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        if (fabs(mjds_st[ii] - last_file_epoch) * SECPERDAY > 1.0e-5)
            newscan = 1;
        else
            newscan = 0;
        last_file_epoch = mjds_st[ii];
        filedatalen_st[ii] = chkfilelen(files[ii], 1) - BPP_HEADER_SIZE;
        if (splitbytes_st[ii - 1]) {
            if (newscan) {      /* Fill the buffer with padding */
                //memset(splitbytes_buffer[ii - 1], (padval << 4) + padval,
                memset(splitbytes_buffer[ii - 1], (8 << 4) + 8,
                       splitbytes_st[ii - 1]);
            } else {            /* Fill the buffer with the data from the missing block part */
                chkfread(splitbytes_buffer[ii - 1], splitbytes_st[ii - 1], 1,
                         files[ii]);
                filedatalen_st[ii] -= splitbytes_st[ii - 1];
            }
        }
        numblks_st[ii] = filedatalen_st[ii] / bytesperblk_st;
        splitbytes_st[ii] = filedatalen_st[ii] % bytesperblk_st;
        if (splitbytes_st[ii]) {
            splitbytes_st[ii] = bytesperblk_st - filedatalen_st[ii] % bytesperblk_st;
            printf
                ("  File %2d has a non-integer number of complete blocks and/or samples!\n"
                 "\tApplying work-around.  (bytes split = %d)\n", ii + 1,
                 splitbytes_st[ii]);
            if (numfiles > ii + 1) {
                numblks_st[ii]++;
                /* This makes a memory leak as it is never freed (although the OS */
                /* does free it when the program exits -- which is usually when   */
                /* we want it freed anyways...)                                   */
                splitbytes_buffer[ii] = (unsigned char *) malloc(splitbytes_st[ii]);
            }
        }
        numpts_st[ii] = numblks_st[ii] * ptsperblk_st;
        times_st[ii] = numpts_st[ii] * dt_st;
        if (newscan) {
            elapsed_st[ii] = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                                          idata_st[ii - 1].mjd_i,
                                          idata_st[ii - 1].mjd_f);
        } else {
            elapsed_st[ii] = times_st[ii - 1];
            idata_st[ii].mjd_f = idata_st[ii - 1].mjd_f + elapsed_st[ii] / SECPERDAY;
            idata_st[ii].mjd_i = idata_st[ii - 1].mjd_i;
            if (idata_st[ii].mjd_f >= 1.0) {
                idata_st[ii].mjd_f -= 1.0;
                idata_st[ii].mjd_i++;
            }
            mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        }
        padpts_st[ii - 1] = (long long) ((elapsed_st[ii] - times_st[ii - 1]) /
                                         dt_st + 0.5);
        elapsed_st[ii] += elapsed_st[ii - 1];
        N_st += numpts_st[ii] + padpts_st[ii - 1];
        startblk_st[ii] = (double) (N_st - numpts_st[ii]) / ptsperblk_st + 1;
        endblk_st[ii] = (double) (N_st) / ptsperblk_st;
    }
    padpts_st[numfiles - 1] = ((long long) ceil(endblk_st[numfiles - 1]) *
                               ptsperblk_st - N_st);
    N_st += padpts_st[numfiles - 1];
    *N = N_st;
    *T = T_st = N_st * dt_st;
    currentfile = currentblock = 0;
#if 0
    {
        double *freq;
        int *nibble;

        for (ii = 0; ii < 96; ii++) {
            gen_channel_mapping(header, &nibble, &freq, NULL);
            chan_index2[ii] = nibble[ii];
            chan_freqs2[ii] = freq[ii] / 1000000.0;
        }
        for (ii = 0; ii < numchan_st; ii++) {
            printf("%3d  %3d  %3d  %10.3f  %3d  %10.3f\n", ii,
                   chan_mapping[ii], chan_index[ii], chan_freqs[chan_index[ii]],
                   chan_index2[ii], chan_freqs2[ii]);
        }
    }
#endif
    if (output) {
        printf("  Number of files = %d\n", numfiles);
        printf("     Points/block = %d\n", ptsperblk_st);
        printf("  Num of channels = %d\n", numchan_st);
        printf(" Total points (N) = %lld\n", N_st);
        printf(" Sample time (dt) = %-14.14g\n", dt_st);
        printf("   Total time (s) = %-14.14g\n\n", T_st);
        printf
            ("File  Start Block    Last Block     Points      Elapsed (s)      Time (s)            MJD           Padding\n");
        printf
            ("----  ------------  ------------  ----------  --------------  --------------  ------------------  ----------\n");
        for (ii = 0; ii < numfiles; ii++)
            printf
                ("%2d    %12.11g  %12.11g  %10lld  %14.13g  %14.13g  %17.12f  %10lld\n",
                 ii + 1, startblk_st[ii], endblk_st[ii], numpts_st[ii],
                 elapsed_st[ii], times_st[ii], mjds_st[ii], padpts_st[ii]);
        printf("\n");
    }
}


void BPP_update_infodata(int numfiles, infodata * idata)
/* Update the onoff bins section in case we used multiple files */
{

    int ii, index = 2;

    idata->N = N_st;
    if (numfiles == 1 && padpts_st[0] == 0) {
        idata->numonoff = 0;
        return;
    }
    /* Determine the topocentric onoff bins */
    idata->numonoff = 1;
    idata->onoff[0] = 0.0;
    idata->onoff[1] = numpts_st[0] - 1.0;
    for (ii = 1; ii < numfiles; ii++) {
        if (padpts_st[ii - 1]) {
            idata->onoff[index] = idata->onoff[index - 1] + padpts_st[ii - 1];
            idata->onoff[index + 1] = idata->onoff[index] + numpts_st[ii];
            idata->numonoff++;
            index += 2;
        } else {
            idata->onoff[index - 1] += numpts_st[ii];
        }
    }
    if (padpts_st[numfiles - 1]) {
        idata->onoff[index] = idata->onoff[index - 1] + padpts_st[numfiles - 1];
        idata->onoff[index + 1] = idata->onoff[index];
        idata->numonoff++;
    }
}


void print_BPP_hdr(BPP_SEARCH_HEADER * hdr)
/* Output a BPP header in human readable form */
{
    int mjd_i;
    double mjd_d;

    calc_BPP_chans(hdr);
    printf("\n'%s' (version %d)\n", hdr->head, hdr->header_version);
    printf("                     Target = %s\n", hdr->target_name);
    if (strlen(hdr->obs_group))
        printf("                Observed by = %s\n", hdr->obs_group);
    printf("   Scan number (DDDYYYY###) = %010u\n", hdr->scan_num);
    printf("           Scan file number = %d\n", hdr->scan_file_number);
    /*
       printf("          File size (bytes) = %d\n", hdr->file_size);
     */
    printf("                Tape number = %d\n", hdr->tape_num);
    printf("           Tape file number = %d\n", hdr->tape_file_number);
    printf("       LMST in sec since 0h = %d\n", hdr->lmst);
    printf("        UTC date (DDD:YYYY) = %s\n", hdr->date);
    printf("        UTC time (HH:MM:SS) = %s\n", hdr->start_time);
    printf("             MJD start time = %.11f\n",
           UT_strings_to_MJD(hdr->date, hdr->start_time, &mjd_i, &mjd_d));
    printf("    RA (J2000, HHMMSS.SSSS) = %.4f\n", hdr->ra_2000);
    printf("   DEC (J2000, DDMMSS.SSSS) = %.4f\n", hdr->dec_2000);
    printf("     Integration length (s) = %.17g\n", hdr->length_of_integration);
    printf("           Sample time (us) = %.17g\n", hdr->samp_rate);
    if (hdr->bit_mode == 4)
        printf("               Channel mode = Powers\n");
    else if (hdr->bit_mode == -4)
        printf("               Channel mode = Direct voltages\n");
    else
        printf("        UNKNOWN CHANNEL MODE!!\n");
    if (hdr->cb_sum_polarizations) {
        printf("      Polarizations summed? = Yes\n");
        printf("         Number of channels = %d\n", hdr->num_chans);
        printf("    Overall bandwidth (MHz) = %.17g\n",
               hdr->num_chans * delta_freq_st);
    } else {
        printf("      Polarizations summed? = No\n");
        printf("         Number of channels = %d x 2 IFs\n", hdr->num_chans / 2);
        printf("    Overall bandwidth (MHz) = %.17g\n",
               0.5 * hdr->num_chans * delta_freq_st);
    }
    printf("    Channel bandwidth (MHz) = %.17g\n", delta_freq_st);
    printf("  Lowest channel freq (MHz) = %.17g\n", ch1_freq_st);
    printf("          Middle freq (MHz) = %.17g\n", mid_freq_st);
    printf("  LO freq used for IF (MHz) = %.17g\n", hdr->rf_lo);
    printf("\n");
}
