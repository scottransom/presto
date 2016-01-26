#include "presto.h"
#include "mask.h"
#include "multibeam.h"

#define MAXNUMCHAN 1024

/* All of the following have an _st to indicate static */
static long long numpts_st[MAXPATCHFILES], padpts_st[MAXPATCHFILES], N_st;
static int numblks_st[MAXPATCHFILES], nfilter_st, filter_st = 0, offsetbytes_st = 0;
static int decreasing_freqs_st, numchan_st, ptsperblk_st, bytesperpt_st;
static double times_st[MAXPATCHFILES], mjds_st[MAXPATCHFILES];
static double elapsed_st[MAXPATCHFILES], T_st, dt_st;
static double startblk_st[MAXPATCHFILES], endblk_st[MAXPATCHFILES];
static infodata idata_st[MAXPATCHFILES];
static const unsigned char padval = 0x55;       /* 01010101 (even channels are 0) */
static unsigned char chanmask[MAXNUMCHAN];
static unsigned char databuffer[DATLEN * 2];
static int currentfile, currentblock;
static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
static int using_MPI = 0;
static float clip_sigma_st = 0.0;

static int clip_PKMB_times(unsigned char *rawdata, int ptsperblk, int numchan,
                           float clip_sigma);

void get_PKMB_static(int *decreasing_freqs, float *clip_sigma,
                     int *bytesperpt, int *offsetbytes)
{
    *decreasing_freqs = decreasing_freqs_st;
    *clip_sigma = clip_sigma_st;
    *bytesperpt = bytesperpt_st;
    *offsetbytes = offsetbytes_st;
}

void set_PKMB_static(int ptsperblk, int bytesperpt,
                     int numchan, int decreasing_freqs,
                     int offsetbytes, float clip_sigma, double dt)
{
    using_MPI = 1;
    currentblock = 0;
    ptsperblk_st = ptsperblk;
    bytesperpt_st = bytesperpt;
    numchan_st = numchan;
    decreasing_freqs_st = decreasing_freqs;
    offsetbytes_st = offsetbytes;
    clip_sigma_st = clip_sigma;
    dt_st = dt;
}


void get_PKMB_file_info(FILE * files[], int numfiles, float clipsig,
                        long long *N, int *ptsperblock, int *numchan,
                        double *dt, double *T, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of PKMB rawfiles that you want to patch     */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
    int ii;
    double block_offset;
    char ctmp[12];
    PKMB_tapehdr header;

    if (numfiles > MAXPATCHFILES) {
        printf("\nThe number of input files (%d) is greater than \n", numfiles);
        printf("   MAXPATCHFILES=%d.  Exiting.\n\n", MAXPATCHFILES);
        exit(0);
    }
    chkfread(&header, 1, HDRLEN, files[0]);
    rewind(files[0]);
    nfilter_st = strtol(header.nfilter, NULL, 10);
    /* Choose the higher filterbank for 10/50 data */
    if (nfilter_st > 1) {
        char *envval = getenv("PARKES_FILTER");
        filter_st = 1;
        if (envval != NULL) {
            strlower(envval);
            if (envval[0] == 'h' && envval[1] == 'i') {
                filter_st = 0;
                printf("Looks like this is 10cm/50cm data.  Using 10cm part.\n");
            } else {
                printf("Looks like this is 10cm/50cm data.  Using 50cm part.\n");
            }
        } else {
            printf("Looks like this is 10cm/50cm data.  Using 50cm part.\n");
        }
    }
    /* Are we going to clip the data? */
    if (clipsig > 0.0)
        clip_sigma_st = clipsig;
    PKMB_hdr_to_inf(&header, &idata_st[0]);
    numchan_st = *numchan = idata_st[0].num_chan;
    bytesperpt_st = numchan_st / 8;
    if (nfilter_st == 1) {
        ptsperblk_st = *ptsperblock = DATLEN * 8 / numchan_st;
        offsetbytes_st = 0;
        /* New 288 channel wide-band mode */
        if ((idata_st[0].num_chan == 288) && (idata_st[0].chan_wid == 3.0)) {
            ptsperblk_st = *ptsperblock = DATLEN * 8 / 384;
            bytesperpt_st = 384 / 8;
        }
    } else {                    /* 10cm/50cm data */
        ptsperblk_st = *ptsperblock = DATLEN * 8 / 512;
        bytesperpt_st = 512 / 8;
        if (filter_st == 0) {   /* 10cm */
            offsetbytes_st = 0;
        } else {                /* 50cm */
            offsetbytes_st = 256 / 8;
        }
    }
    decreasing_freqs_st = (strtod(header.chanbw[filter_st], NULL) > 0.0) ? 0 : 1;
    numblks_st[0] = chkfilelen(files[0], RECLEN);
    numpts_st[0] = numblks_st[0] * ptsperblk_st;
    N_st = numpts_st[0];
    dt_st = *dt = idata_st[0].dt;
    times_st[0] = numpts_st[0] * dt_st;
    mjds_st[0] = idata_st[0].mjd_i + idata_st[0].mjd_f;
    elapsed_st[0] = 0.0;
    sprintf(ctmp, " %.8s ", header.blk_cntr);
    sscanf(ctmp, "%8lf", &block_offset);
    startblk_st[0] = block_offset;
    endblk_st[0] = (double) numpts_st[0] / ptsperblk_st + block_offset - 1;
    padpts_st[0] = padpts_st[numfiles - 1] = 0;
    for (ii = 1; ii < numfiles; ii++) {
        chkfread(&header, 1, HDRLEN, files[ii]);
        rewind(files[ii]);
        PKMB_hdr_to_inf(&header, &idata_st[ii]);
        if (idata_st[ii].num_chan != numchan_st) {
            printf("Number of channels (file %d) is not the same!\n\n", ii + 1);
        }
        if (idata_st[ii].dt != dt_st) {
            printf("Sample time (file %d) is not the same!\n\n", ii + 1);
        }
        numblks_st[ii] = chkfilelen(files[ii], RECLEN);
        numpts_st[ii] = numblks_st[ii] * ptsperblk_st;
        times_st[ii] = numpts_st[ii] * dt_st;
        mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        elapsed_st[ii] = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                                      idata_st[ii - 1].mjd_i,
                                      idata_st[ii - 1].mjd_f);
        padpts_st[ii - 1] =
            (long long) ((elapsed_st[ii] - times_st[ii - 1]) / dt_st + 0.5);
        elapsed_st[ii] += elapsed_st[ii - 1];
        N_st += numpts_st[ii] + padpts_st[ii - 1];
        startblk_st[ii] = (double) (N_st - numpts_st[ii]) /
            ptsperblk_st + block_offset;
        endblk_st[ii] = (double) (N_st) / ptsperblk_st + block_offset - 1;
    }
    padpts_st[numfiles - 1] = ((long long) ceil(endblk_st[numfiles - 1] -
                                                block_offset + 1.0) *
                               ptsperblk_st - N_st);
    N_st += padpts_st[numfiles - 1];
    *N = N_st;
    *T = T_st = N_st * dt_st;
    currentfile = currentblock = 0;
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

void PKMB_update_infodata(int numfiles, infodata * idata)
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


int skip_to_PKMB_rec(FILE * infiles[], int numfiles, int rec)
/* This routine skips to the record 'rec' in the input files   */
/* *infiles.  *infiles contains 1 bit digitized data from the  */
/* PKMB backend at Parkes.  Returns the record skipped to.     */
{
    double floor_blk;
    int filenum = 0;

    if (rec < startblk_st[0])
        rec += (startblk_st[0] - 1);
    if (rec > 0 && rec < endblk_st[numfiles - 1]) {

        /* Find which file we need */
        while (rec > endblk_st[filenum])
            filenum++;

        currentblock = rec - 1;
        shiftbuffer = 1;
        floor_blk = floor(startblk_st[filenum]);

        /* Set the data buffer to all padding just in case */
        memset(databuffer, padval, DATLEN * 2);

        /* Warning:  I'm not sure if the following is correct. */
        /*   If really needs accurate testing to see if my     */
        /*   offsets are correct.  Bottom line, don't trust    */
        /*   a TOA determined using the following!             */

        if (rec < startblk_st[filenum]) {       /* Padding region */
            currentfile = filenum - 1;
            chkfileseek(infiles[currentfile], 0, 1, SEEK_END);
            bufferpts = padpts_st[currentfile] % ptsperblk_st;
            padnum = ptsperblk_st * (rec - endblk_st[currentfile] - 1);
            /*
               printf("Padding:  currentfile = %d  bufferpts = %d  padnum = %d\n", 
               currentfile, bufferpts, padnum);
             */
        } else {                /* Data region */
            currentfile = filenum;
            chkfileseek(infiles[currentfile], rec - startblk_st[filenum],
                        RECLEN, SEEK_SET);
            bufferpts =
                (int) ((startblk_st[filenum] - floor_blk) * ptsperblk_st + 0.5);
            padnum = 0;
            /*
               printf("Data:  currentfile = %d  bufferpts = %d  padnum = %d\n", 
               currentfile, bufferpts, padnum);
             */
        }

    } else {
        printf("\n rec = %d out of range in skip_to_PKMB_rec()\n", rec);
        exit(1);
    }
    return rec;
}


int read_PKMB_rawblock(FILE * infiles[], int numfiles,
                       PKMB_tapehdr * hdr, unsigned char *data, int *padding)
/* This routine reads a single record from the         */
/* input files *infiles which contain 1 bit digitized  */
/* data from the PKMB pulsar backend at Parkes.        */
/* Length of a PKMB record is 640 bytes for the header */
/* plus 48k of data = 49792 bytes.                     */
/* The header of the record read is placed in hdr.     */
/* *data must be pre-allocated with a size of 48k.     */
/* If padding is returned as 1, then padding was       */
/* added and statistics should not be calculated       */
{
    int offset, numtopad = 0;
    unsigned char record[RECLEN];

    /* If our buffer array is offset from last time */
    /* copy the second part into the first.         */

    if (bufferpts && shiftbuffer)
        memcpy(databuffer, databuffer + DATLEN, bufferpts * bytesperpt_st);
    shiftbuffer = 1;

    /* Make sure our current file number is valid */

    if (currentfile >= numfiles)
        return 0;

    /* First, attempt to read data from the current file */

    if (fread(record, RECLEN, 1, infiles[currentfile])) {       /* Got Data */
        *padding = 0;
        /* Put the new header into the header array */
        memcpy(hdr, record, HDRLEN);
        /* Put the new data into the databuffer or directly */
        /* into the return array if the bufferoffset is 0.  */
        if (bufferpts) {
            offset = bufferpts * bytesperpt_st;
            memcpy(databuffer + offset, record + HDRLEN, DATLEN);
            memcpy(data, databuffer, DATLEN);
        } else {
            memcpy(data, record + HDRLEN, DATLEN);
        }
        currentblock++;
        return 1;
    } else {                    /* Didn't get data */
        if (feof(infiles[currentfile])) {       /* End of file? */
            numtopad = padpts_st[currentfile] - padnum;
            if (numtopad) {     /* Pad the data? */
                *padding = 1;
                if (numtopad >= ptsperblk_st - bufferpts) {     /* Lots of padding */
                    if (bufferpts) {    /* Buffer the padding? */
                        /* Add the amount of padding we need to */
                        /* make our buffer offset = 0           */
                        numtopad = ptsperblk_st - bufferpts;
                        memset(databuffer + bufferpts * bytesperpt_st,
                               padval, numtopad * bytesperpt_st);
                        /* Copy the new data/padding into the output array */
                        memcpy(data, databuffer, DATLEN);
                        bufferpts = 0;
                    } else {    /* Add a full record of padding */
                        numtopad = ptsperblk_st;
                        memset(data, padval, DATLEN);
                    }
                    padnum += numtopad;
                    currentblock++;
                    /* If done with padding reset padding variables */
                    if (padnum == padpts_st[currentfile]) {
                        padnum = 0;
                        currentfile++;
                    }
                    return 1;
                } else {        /* Need < 1 block (or remaining block) of padding */
                    int pad;
                    /* Add the remainder of the padding and */
                    /* then get a block from the next file. */
                    memset(databuffer + bufferpts * bytesperpt_st,
                           padval, numtopad * bytesperpt_st);
                    padnum = 0;
                    currentfile++;
                    shiftbuffer = 0;
                    bufferpts += numtopad;
                    return read_PKMB_rawblock(infiles, numfiles, hdr, data, &pad);
                }
            } else {            /* No padding needed.  Try reading the next file */
                currentfile++;
                shiftbuffer = 0;
                return read_PKMB_rawblock(infiles, numfiles, hdr, data, padding);
            }
        } else {
            printf("\nProblem reading record from PKMB data file:\n");
            printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                   currentfile, currentblock);
            exit(1);
        }
    }
}


int read_PKMB_rawblocks(FILE * infiles[], int numfiles,
                        unsigned char rawdata[], int numblocks, int *padding)
/* This routine reads numblocks PKMB records from the input */
/* files *infiles.  The raw bit data is returned in rawdata */
/* which must have a size of numblocks*DATLEN.  The number  */
/* of blocks read is returned.                              */
/* If padding is returned as 1, then padding was added      */
/* and statistics should not be calculated                  */
{
    int ii, retval = 0, pad, numpad = 0;
    PKMB_tapehdr hdr;

    *padding = 0;
    for (ii = 0; ii < numblocks; ii++) {
        retval += read_PKMB_rawblock(infiles, numfiles, &hdr,
                                     rawdata + ii * DATLEN, &pad);
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


int read_PKMB(FILE * infiles[], int numfiles, float *data,
              int numpts, double *dispdelays, int *padding,
              int *maskchans, int *nummasked, mask * obsmask)
/* This routine reads numpts from the PKMB raw input   */
/* files *infiles.  These files contain 1 bit data     */
/* from the PKMB backend at Parkes.  Time delays and   */
/* a mask are applied to each channel.  It returns     */
/* the # of points read if successful, 0 otherwise.    */
/* If padding is returned as 1, then padding was       */
/* added and statistics should not be calculated.      */
/* maskchans is an array of length numchans contains   */
/* a list of the number of channels that were masked.  */
/* The # of channels masked is returned in nummasked.  */
/* obsmask is the mask structure to use for masking.   */
{
    int ii, jj, numread = 0, offset;
    double starttime = 0.0;
    static unsigned char *tempzz, *raw, *rawdata1, *rawdata2;
    static unsigned char *currentdata, *lastdata;
    static int firsttime = 1, numblocks = 1, allocd = 0, mask = 0;
    static double duration = 0.0, timeperblk = 0.0;

    *nummasked = 0;
    if (firsttime) {
        if (numpts % ptsperblk_st) {
            printf("numpts must be a multiple of %d in read_PKMB()!\n",
                   ptsperblk_st);
            exit(1);
        } else
            numblocks = numpts / ptsperblk_st;

        if (obsmask->numchan)
            mask = 1;
        raw = gen_bvect(numblocks * DATLEN);
        rawdata1 = gen_bvect(numblocks * SAMPPERBLK);
        rawdata2 = gen_bvect(numblocks * SAMPPERBLK);
        allocd = 1;
        timeperblk = ptsperblk_st * dt_st;
        duration = numblocks * timeperblk;

        numread = read_PKMB_rawblocks(infiles, numfiles, raw, numblocks, padding);
        if (numread != numblocks && allocd) {
            printf("Problem reading the raw PKMB data file.\n\n");
            vect_free(raw);
            vect_free(rawdata1);
            vect_free(rawdata2);
            allocd = 0;
            return 0;
        }

        currentdata = rawdata1;
        lastdata = rawdata2;

        if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
            if (*nummasked == -1)       /* If all channels are masked */
                memset(raw, padval, numblocks * DATLEN);
        }

        for (ii = 0; ii < numpts; ii++)
            convert_PKMB_point(raw + ii * bytesperpt_st,
                               currentdata + ii * numchan_st);

        /* Clip nasty RFI if requested */
        if (clip_sigma_st > 0.0)
            clip_PKMB_times(currentdata, numpts, numchan_st, clip_sigma_st);

        if (*nummasked > 0) {   /* Only some of the channels are masked */
            for (ii = 0; ii < *nummasked; ii++)
                chanmask[ii] = maskchans[ii] & 0x01;
            for (ii = 0; ii < numpts; ii++) {
                offset = ii * numchan_st;
                for (jj = 0; jj < *nummasked; jj++)
                    currentdata[offset + maskchans[jj]] = chanmask[jj];
            }
        }

        SWAP(currentdata, lastdata);
        firsttime = 0;
    }

    /* Read, convert and de-disperse */

    if (allocd) {
        numread = read_PKMB_rawblocks(infiles, numfiles, raw, numblocks, padding);

        if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
            if (*nummasked == -1)       /* If all channels are masked */
                memset(raw, padval, numblocks * DATLEN);
        }

        for (ii = 0; ii < numpts; ii++)
            convert_PKMB_point(raw + ii * bytesperpt_st,
                               currentdata + ii * numchan_st);

        /* Clip nasty RFI if requested */
        if (clip_sigma_st > 0.0)
            clip_PKMB_times(currentdata, numpts, numchan_st, clip_sigma_st);

        if (*nummasked > 0) {   /* Only some of the channels are masked */
            for (ii = 0; ii < *nummasked; ii++)
                chanmask[ii] = maskchans[ii] & 0x01;
            for (ii = 0; ii < numpts; ii++) {
                offset = ii * numchan_st;
                for (jj = 0; jj < *nummasked; jj++)
                    currentdata[offset + maskchans[jj]] = chanmask[jj];
            }
        }

        dedisp(currentdata, lastdata, numpts, numchan_st, dispdelays, data);
        SWAP(currentdata, lastdata);

        if (numread != numblocks) {
            vect_free(raw);
            vect_free(rawdata1);
            vect_free(rawdata2);
            allocd = 0;
        }
        return numread * ptsperblk_st;
    } else {
        return 0;
    }
}


void get_PKMB_channel(int channum, float chandat[],
                      unsigned char rawdata[], int numblocks)
/* Return the values for channel 'channum' of a block of       */
/* 'numblocks' raw PKMB data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                */
/* read_PKMB_rawblocks(), and 'chandat' must have at least     */
/* 'numblocks' * 'ptsperblk_st' spaces.                        */
/* Channel 0 is assumed to be the lowest freq channel.         */
{
    int ii, bit;

    if (channum > numchan_st || channum < 0) {
        printf("\nchannum = %d is out of range in get_PKMB_channel()!\n\n", channum);
        exit(1);
    }
    bit = (decreasing_freqs_st) ? numchan_st - 1 - channum : channum;
    bit += offsetbytes_st * 8;
    for (ii = 0; ii < numblocks * ptsperblk_st; ii++)
        chandat[ii] = (float) GET_BIT(rawdata + ii * bytesperpt_st, bit);
}


void PKMB_hdr_to_inf(PKMB_tapehdr * hdr, infodata * idata)
/* Convert PKMB header into an infodata structure */
{
    double tmp1, tmp2;
    char ctmp[100];
    int ptsperrec = 0;
    int itmp1, itmp2;

    sprintf(idata->object, "%.16s", hdr->pname);
    sscanf(hdr->ra_start, "%3d:%2d:%9lf", &idata->ra_h, &idata->ra_m, &idata->ra_s);
    sscanf(hdr->dec_start,
           "%3d:%2d:%9lf", &idata->dec_d, &idata->dec_m, &idata->dec_s);
    /* make sure to catch minus sign for "+/-00:XX:XX.XX" Dec string in header */
    if (idata->dec_d == 0 && (hdr->dec_start)[0] == '-') {
        idata->dec_m = -idata->dec_m;
        idata->dec_s = -idata->dec_s;
    }
    if (!strncmp(hdr->telid, "PARKES", 6))
        strcpy(idata->telescope, "Parkes");
    else if (strncmp(hdr->telid, "  0.", 4) == 0 ||
             strncmp(hdr->telid, "Lovell", 6) == 0 ||
             strncmp(hdr->telid, "LOVELL", 6) == 0)
        strcpy(idata->telescope, "Jodrell");
    itmp1 = strtol(hdr->ibeam, NULL, 10);
    itmp2 = strtol(hdr->nbeam, NULL, 10);
    tmp1 = strtod(hdr->chanbw[filter_st], NULL);
    idata->chan_wid = fabs(tmp1);
    sscanf(hdr->nchan[filter_st], "%4d", &idata->num_chan);
    ptsperrec = DATLEN / (idata->num_chan / 8);
    if (nfilter_st == 2) {
        sprintf(idata->instrument, "10cm+50cm Receiver");
    } else {
        if ((idata->num_chan == 384) && (idata->chan_wid == 3.0)) {
            sprintf(idata->instrument, "288 channel wideband");
            idata->num_chan = 288;
        } else {
            sprintf(idata->instrument, "Multibeam (Beam %d of %d)", itmp1, itmp2);
        }
    }
    /* throw away the top 64 channels of the 10cm filterbank */
    if (nfilter_st == 2 && filter_st == 0)
        idata->num_chan = 192;
    sscanf(hdr->samp_int[filter_st], "%12lf", &idata->dt);
    idata->dt *= 0.001;
    sscanf(hdr->nblk_read, "%8lf", &tmp1);
    sprintf(ctmp, " %.8s ", hdr->obs_time);
    tmp2 = strtod(ctmp, NULL);
    if (tmp1) {
        idata->N = ptsperrec * tmp1;
    } else {
        idata->N = tmp2 / idata->dt;
    }
    tmp1 = strtod(hdr->chanbw[filter_st], NULL);
    idata->chan_wid = fabs(tmp1);
    idata->freqband = idata->num_chan * idata->chan_wid;
    idata->freq = strtod(hdr->freq_ch1[filter_st], NULL);
    if (tmp1 < 0.0) {
        idata->freq -= (idata->freqband - idata->chan_wid);
    }
    idata->fov = 1.2 * SOL * 3600.0 / (1000000.0 * idata->freq * 64 * DEGTORAD);
    idata->mjd_i = strtol(hdr->mjd, NULL, 10);
    sscanf(hdr->ut_start, "%lf:%lf:%lf", &idata->mjd_f, &tmp1, &tmp2);
    idata->mjd_f = (idata->mjd_f + (tmp1 + tmp2 / 60.0) / 60.0) / 24.0;
    sprintf(ctmp, " %.8s ", hdr->blk_cntr);
    tmp1 = strtod(ctmp, NULL) - 1.0;
    idata->mjd_f += (tmp1 * ptsperrec * idata->dt) / SECPERDAY;
    if (idata->mjd_f >= 1.0) {
        idata->mjd_f -= 1.0;
        idata->mjd_i += 1;
    }
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    strcpy(idata->analyzer, "Scott Ransom");
    strcpy(idata->observer, "--");
    sprintf(idata->notes,
            "Topo UT Date & Time at file start = %.8s, %.16s\n    From tape %.6s  file #%.4s  block #%.8s\n    Comment: %.64s\n",
            hdr->date, hdr->ut_start, hdr->tape_lbl, hdr->file_cntr, hdr->blk_cntr,
            hdr->comment);

}


void print_PKMB_hdr(PKMB_tapehdr * hdr)
/* Output a PKMB header in human readable form */
{
    printf("\n");
    printf(" Program name                                      = %.6s\n",
           hdr->program);
    printf(" Version number 'vv.v'                             = %.4s\n",
           hdr->version);
    printf(" Tape label                                        = %.6s\n",
           hdr->tape_lbl);
    printf(" File number on tape                               = %.4s\n",
           hdr->file_cntr);
    printf(" Block counter per beam in file                    = %.8s\n",
           hdr->blk_cntr);
    printf(" UT Date 'yymmdd  '                                = %.8s\n", hdr->date);
    printf(" MJD at 00h UT                                     = %.8s\n", hdr->mjd);
    printf(" UT at file start  ' hh:mm:ss.ssssss'              = %.16s\n",
           hdr->ut_start);
    printf(" LST at file start ' hh:mm:ss.ss'                  = %.12s\n",
           hdr->lst_start);
    printf(" Coord system ' 4' = Galactic, ' 5' = J2000        = %.2s\n",
           hdr->coord_type);
    printf(" RA at file start  ' hh:mm:ss.ssss'                = %.16s\n",
           hdr->ra_start);
    printf(" Dec at file start '-dd:mm:ss.sss'                 = %.16s\n",
           hdr->dec_start);
    printf(" Galactic long at file start (deg) 'ddd.dddd'      = %.8s\n",
           hdr->gl_start);
    printf(" Galactic lat. at file start (deg) '-dd.dddd'      = %.8s\n",
           hdr->gb_start);
    printf(" FA, PA or GPA (deg) '-ddd.ddd'                    = %.8s\n",
           hdr->fangle);
    printf(" Obs length 'ssss.sss'                             = %.8s\n",
           hdr->obs_time);
    printf(" Comment                                           = %.64s\n",
           hdr->comment);
    printf(" Number of filter systems, ' 1' or ' 2'            = %.2s\n",
           hdr->nfilter);
    {
        int ii;
        int nfilter = strtol(hdr->nfilter, NULL, 10);

        for (ii = 0; ii < nfilter; ii++) {
            printf(" Channel incrt (MHz, -ve if inverted) '-b.bbbbb'   = %.8s\n",
                   hdr->chanbw[ii]);
            printf(" Nr of channels in each filter system              = %.4s\n",
                   hdr->nchan[ii]);
            printf(" RF of first channel centre 'ffffff.fffff'         = %.12s\n",
                   hdr->freq_ch1[ii]);
            printf(" Samp intval in ms (1 = 2 always) 'mm.mmmmmmmmm'   = %.12s\n",
                   hdr->samp_int[ii]);
            printf(" Samp per group (For backwards compat., now = 1)   = %.4s\n",
                   hdr->samp_grp[ii]);
        }
    }
    printf(" Samples per block (Was groups per block)          = %.8s\n",
           hdr->samp_blk);
    printf(" Seconds per tape block per beam 'ss.sssss'        = %.8s\n",
           hdr->blk_sec);
    printf(" 0 = none, 1 = fixed FA, 2 = fixed PA or GPA       = %.2s\n",
           hdr->fdctrl);
    printf(" Data type code ('2'=pol, '3'=normal, '5'=dedisp)  = %.1s\n",
           hdr->data[0]);
    printf(" UT at blk start ' hh:mm:ss.ssssss'                = %.16s\n",
           hdr->ut_blk);
    printf(" LST at blk start ' hh:mm:ss.ss'                   = %.12s\n",
           hdr->lst_blk);
    printf(" RA at blk start ' hh:mm:ss.ssss'                  = %.16s\n",
           hdr->ra_blk);
    printf(" Dec at blk start '-dd:mm:ss.sss'                  = %.16s\n",
           hdr->dec_blk);
    printf(" Gal long. at blk start 'ddd.dddd'                 = %.8s\n",
           hdr->gl_blk);
    printf(" Gal lat at blk start '-dd.dddd'                   = %.8s\n",
           hdr->gb_blk);
    printf(" Zenith angle at blk start (deg) 'ddd.dddd'        = %.8s\n",
           hdr->zen_blk);
    printf(" Azimuth angle at blk start (deg) 'ddd.dddd'       = %.8s\n",
           hdr->az_blk);
    printf(" Attenuator settings (A,B) for Sys 1,2 'dd.d'      = %.4s\n",
           hdr->atten[0]);
    printf(" Total powers for 32-channel groups 'iiii'         = %.4s\n",
           hdr->tpower[0]);
    printf(" Nr of tape blocks in disk file                    = %.8s\n",
           hdr->nblk_read);
    printf(" Scan rate in `long' (deg/min) '-r.rrrrr'          = %.8s\n",
           hdr->scan_ratel);
    printf(" Scan rate in `lat'  (deg/min) '-r.rrrrr'          = %.8s\n",
           hdr->scan_rateb);
    printf(" Total number of beams                             = %.4s\n",
           hdr->nbeam);
    printf(" Beam number                                       = %.4s\n",
           hdr->ibeam);
    printf(" Psrname (B or J) or pointing name                 = %.16s\n",
           hdr->pname);
    printf(" Config file name                                  = %.8s\n",
           hdr->cfg_file);
    printf(" Schedule file name                                = %.8s\n",
           hdr->schd_file);
    printf(" Nr of bits per sample                             = %.2s\n",
           hdr->nbits[0]);
    printf(" DM for dedispersion 'ddddd.dd'                    = %.8s\n", hdr->ddm);
    printf(" Nr of channels per dedispersed band               = %.4s\n",
           hdr->nddch[0]);
    printf(" Tel status 'k0': On posn k(0-9)                   = %.2s\n", hdr->move);
    printf("            '01': On grid posn\n");
    printf("            '02': Off position\n");
    printf("            '03': GC scan\n");
    printf("            '04': SC scan\n");
    printf(" Pointing error (arcmin) 'mmm.mm'                  = %.6s\n",
           hdr->pnterr);
    printf(" Dedisp type: '  ': Norm, ' T': Tree, ' D': PDM    = %.2s\n", hdr->tree);
    printf(" Filter system for disk file                       = %.2s\n",
           hdr->ifilter);
    printf(" Telescope ID (as in tel.def)                      = %.10s\n",
           hdr->telid);
    printf(" Parallactic angle (deg) '-ddd.ddd'                = %.8s\n",
           hdr->pangle);
    printf(" Nr of samples per beam switch                     = %.8s\n", hdr->bsw);
    printf(" Cal cycle period in samples                       = %.4s\n",
           hdr->cal_cycle);
    printf(" Spare                                             = %.22s\n",
           hdr->spare);
    printf("\n");
}


void convert_PKMB_point(unsigned char *bits, unsigned char *bytes)
/* This routine converts 1 bit digitized data */
/* into an array of 'numchan' bytes.          */
{
    int ii, jj;
    int true_bytesperpt = numchan_st / 8;

    bits += offsetbytes_st;     /* for 10cm/50cm data */
    if (decreasing_freqs_st) {
        for (ii = true_bytesperpt - 1, jj = 0; ii >= 0; ii--, jj += 8) {
            bytes[jj] = (bits[ii] >> 0x07) & 0x01;
            bytes[jj + 1] = (bits[ii] >> 0x06) & 0x01;
            bytes[jj + 2] = (bits[ii] >> 0x05) & 0x01;
            bytes[jj + 3] = (bits[ii] >> 0x04) & 0x01;
            bytes[jj + 4] = (bits[ii] >> 0x03) & 0x01;
            bytes[jj + 5] = (bits[ii] >> 0x02) & 0x01;
            bytes[jj + 6] = (bits[ii] >> 0x01) & 0x01;
            bytes[jj + 7] = bits[ii] & 0x01;
        }
    } else {
        for (ii = 0, jj = 0; ii < true_bytesperpt; ii++, jj += 8) {
            bytes[jj] = bits[ii] & 0x01;
            bytes[jj + 1] = (bits[ii] >> 0x01) & 0x01;
            bytes[jj + 2] = (bits[ii] >> 0x02) & 0x01;
            bytes[jj + 3] = (bits[ii] >> 0x03) & 0x01;
            bytes[jj + 4] = (bits[ii] >> 0x04) & 0x01;
            bytes[jj + 5] = (bits[ii] >> 0x05) & 0x01;
            bytes[jj + 6] = (bits[ii] >> 0x06) & 0x01;
            bytes[jj + 7] = (bits[ii] >> 0x07) & 0x01;
        }
    }
}

int clip_PKMB_times(unsigned char *rawdata, int ptsperblk, int numchan,
                    float clip_sigma)
/* Perform time-domain clipping of rawdata.   This is a 2D   */
/* array with ptsperblk*numchan points, each of which is an  */
/* unsigned char.  The clipping is done at clip_sigma sigma  */
/* above/below the running mean.                             */
{
    static float running_avg = 0.0, running_std = 0.0;
    static int blocksread = 0;
    float *zero_dm_block, *median_temp;
    float current_med, trigger, running_wgt = 0.05;
    double current_avg = 0.0, current_std = 0.0;
    unsigned char *powptr;
    int ii, jj, clipit = 0, clipped = 0;

    zero_dm_block = gen_fvect(ptsperblk);
    median_temp = gen_fvect(ptsperblk);

    /* Calculate the zero DM time series */
    for (ii = 0; ii < ptsperblk; ii++) {
        zero_dm_block[ii] = 0.0;
        powptr = rawdata + ii * numchan;
        for (jj = 0; jj < numchan; jj++)
            zero_dm_block[ii] += *powptr++;
        median_temp[ii] = zero_dm_block[ii];
    }
    current_med = median(median_temp, ptsperblk);

    /* Calculate the current standard deviation and mean  */
    /* but only for data points that are within a certain */
    /* fraction of the median value.  This removes the    */
    /* really strong RFI from the calculation.            */
    {
        float lo_cutoff, hi_cutoff;
        int numgoodpts = 0;

        lo_cutoff = 0.7 * current_med;
        hi_cutoff = 1.3 * current_med;
        /* Find the "good" points */
        for (ii = 0; ii < ptsperblk; ii++) {
            if (zero_dm_block[ii] > lo_cutoff && zero_dm_block[ii] < hi_cutoff) {
                median_temp[numgoodpts] = zero_dm_block[ii];
                numgoodpts++;
            }
        }
        /* Calculate the current average and stddev */
        if (numgoodpts < 1) {
            current_avg = running_avg;
            current_std = running_std;
        } else {
            avg_var(median_temp, numgoodpts, &current_avg, &current_std);
            current_std = sqrt(current_std);
        }
    }

    /* Update a pseudo running average and stdev */
    if (blocksread) {
        running_avg =
            (running_avg * (1.0 - running_wgt) + running_wgt * current_avg);
        running_std =
            (running_std * (1.0 - running_wgt) + running_wgt * current_std);
    } else {
        running_avg = current_avg;
        running_std = current_std;
        if (running_avg == 0.0 || current_avg == 0.0)
            printf("BAD RFI IN BLOCK#1!!!\n\n");
    }

    /* See if any points need clipping */
    trigger = clip_sigma * running_std;
    for (ii = 0; ii < ptsperblk; ii++) {
        if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
            clipit = 1;
            break;
        }
    }

    /* Replace the bad channel data with alternating zeros and ones */
    if (clipit) {
        for (ii = 0; ii < ptsperblk; ii++) {
            if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
                powptr = rawdata + ii * numchan;
                for (jj = 0; jj < numchan / 2; jj++) {
                    *powptr++ = 0;
                     /**/ *powptr++ = 1;
                }
                clipped++;
            }
        }
    }
    blocksread++;

    vect_free(zero_dm_block);
    vect_free(median_temp);

    return clipped;
}
