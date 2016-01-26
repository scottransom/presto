#include "presto.h"
#include "mask.h"
#include "spigot.h"
#include "fitsfile.h"
#include "fitshead.h"
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#define DEBUGOUT 0

/*  NOTES:
bytesperblk_st is the number of bytes in the RAW LAGS for a SPIGOT
sampperblk_st  is the number of lags in a block 
*/

/* All of the following have an _st to indicate static */
static SPIGOT_INFO *spigot;
static infodata *idata_st;
static fftwf_plan fftplan;
static long long N_st;
static int decreasing_freqs_st = 0, bytesperpt_st, bits_per_lag_st = 0;
static int numchan_st = 0, numifs_st, ptsperblk_st, bytesperblk_st;
static int sampperblk_st, usewindow_st = 0;     //, vanvleck_st = 0;
static int currentfile, currentblock;
//static int currentfile, currentblock, bufferpts = 0, padnum = 0, shiftbuffer = 1;
static double T_st, dt_st, center_freq_st, *window_st = NULL;
static float clip_sigma_st = 0.0, *lags = NULL;
static float lag_factor[SPIGOT_MAXLAGLEN], lag_offset[SPIGOT_MAXLAGLEN];
//static unsigned char padvals[SPIGOT_MAXLAGLEN], newpadvals[SPIGOT_MAXLAGLEN];
//static unsigned char databuffer[2 * SPIGOT_MAXDATLEN], padval = 128;
//static unsigned char lagbuffer[SPIGOT_MAXLAGLEN];
//static int using_MPI = 0;
static float lag_scale_env = 1.0;

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

static double UT_strings_to_MJD(char *obs_date, char *start_time,
                                int *mjd_day, double *mjd_fracday)
/* Convert date and time strings to a double precision MJD */
{
    int year, month, day, hour, min, sec, err;

    sscanf(remove_whitespace(obs_date), "%4d-%2d-%2d", &year, &month, &day);
    sscanf(remove_whitespace(start_time), "%2d:%2d:%2d", &hour, &min, &sec);
    *mjd_fracday = (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
    *mjd_day = slaCldj(year, month, day, &err);
    return *mjd_day + *mjd_fracday;
}


int read_SPIGOT_header(char *filename, SPIGOT_INFO * spigot)
/* Read and convert SPIGOT header information and place it into */
/* a SPIGOT_INFO structure.  Return 1 if successful, 0 if not.  */
{
    int hdrlen, data_offset, itmp;
    double dtmp1, dtmp2;
    char *hdr;
    static int firsttime = 1;

    hdr = fitsrhead(filename, &hdrlen, &data_offset);
    if (!hdr) {
        printf("\n  Error!  Could not read '%s'!\n\n", filename);
        return (0);
    }
    if (DEBUGOUT) {
        printf("Read %d bytes containing the header of '%s'\n", hdrlen, filename);
        printf("\tTrue header length is %d bytes.\n", data_offset);
    }
    /* Device or program of origin */
    hgets(hdr, "INSTRUME", 40, spigot->instrument);
    remove_whitespace(spigot->instrument);
    /* Version of observing software */
    hgets(hdr, "SOFTVER", 40, spigot->software_vers);
    remove_whitespace(spigot->software_vers);
    /* Telescope used */
    hgets(hdr, "TELESCOP", 40, spigot->telescope);
    remove_whitespace(spigot->telescope);
    /* Source observed */
    hgets(hdr, "OBJECT", 40, spigot->object);
    remove_whitespace(spigot->object);
    /* Name(s) of observer(s) */
    hgets(hdr, "OBSERVER", 40, spigot->observer);
    remove_whitespace(spigot->observer);
    /* Project identifier */
    hgets(hdr, "PROJID", 16, spigot->project_id);
    remove_whitespace(spigot->project_id);
    /* Observation identifier */
    hgets(hdr, "OBSID", 16, spigot->obs_id);
    remove_whitespace(spigot->obs_id);
    /* Start of observation (YYYY-MM-DD) */
    hgets(hdr, "DATE-OBS", 16, spigot->date_obs);
    remove_whitespace(spigot->date_obs);
    /* Start of observation (HH:MM:SS) */
    hgets(hdr, "TIME-OBS", 16, spigot->time_obs);
    remove_whitespace(spigot->time_obs);
    /* Time scale specification (usually UTC) */
    hgets(hdr, "TIMESYS", 16, spigot->time_sys);
    remove_whitespace(spigot->time_sys);
    if (strcmp(spigot->time_sys, "UTC") != 0) {
        printf("\n  Warning!  '%s' is an unrecognized time system in '%s'!\n\n",
               spigot->time_sys, filename);
    }
    /* Offset coordinate mode of GBT */
    hgets(hdr, "INDICSYS", 16, spigot->coord_sys);
    remove_whitespace(spigot->coord_sys);
    /* RA and DEC of observation (deg, J2000) */
    hgetr8(hdr, "RA", &(spigot->ra));
    hgetr8(hdr, "DEC", &(spigot->dec));
    /* Al and El of GBT in degrees */
    hgetr8(hdr, "AZIMUTH", &(spigot->az));
    hgetr8(hdr, "ELEVATN", &(spigot->el));
    /* Polarization recorded (L or C) */
    hgets(hdr, "POL-TYPE", 8, spigot->pol_type);
    remove_whitespace(spigot->pol_type);
    /* Spigot correlator mode */
    hgets(hdr, "MODE", 8, spigot->corr_mode);
    remove_whitespace(spigot->corr_mode);
    /* Sampling time (us) */
    hgetr8(hdr, "TSAMP", &(spigot->dt_us));
    /* Sky center freq (MHz) for sampler 1 */
    hgetr8(hdr, "CENTFREQ", &(spigot->freq_ctr));
    /* Bandwidth (MHz) for sampler 1 */
    hgetr8(hdr, "SAMP-BW", &(spigot->bandwidth));
    spigot->MJD_obs = UT_strings_to_MJD(spigot->date_obs, spigot->time_obs,
                                        &itmp, &dtmp1);
    /* Start time in sec since MJD=40587.0 */
    hgeti4(hdr, "SEC-OBS", &(spigot->sec_obs));
    /* Calculate the MJD a different way to check it */
    dtmp1 = 40587.0 + spigot->sec_obs / SECPERDAY;
    dtmp2 = (dtmp1 - spigot->MJD_obs) * SECPERDAY;
    if (fabs(dtmp2) > 1e-6) {
        printf
            ("\n  Warning!  The quoted start times for '%s' disagree by %.15g s!\n\n",
             filename, dtmp2);
    }
    /* Tracking (T) or drift scan (F) */
    hgetl(hdr, "TRACK", &(spigot->tracking));
    /* Number of scan */
    hgeti4(hdr, "SCAN", &(spigot->scan_number));
    /* Size of header (bytes) */
    hgeti4(hdr, "HEADSIZE", &(spigot->header_len));
    if (spigot->header_len != data_offset) {
        printf("\n  Warning!  Inconsistent header length calculation in '%s'!\n\n",
               filename);
    }
    /* Number of samplers */
    hgeti4(hdr, "SAMPLERS", &(spigot->num_samplers));
    /* Are polarizations summed? */
    hgetl(hdr, "SUMPOL", &(spigot->summed_pols));
    if (spigot->summed_pols)
        numifs_st = 1;
    /* Upper sideband? */
    hgetl(hdr, "UPPERSB", &(spigot->upper_sideband));
    if (!spigot->upper_sideband)
        decreasing_freqs_st = 1;
    /* For now, ignore BITPIX/NAXIS1 in favor of BITS/NLAGS 
       This is because FITS wants files to have BITPIX >= 8, 
       so for 4-bit data BITPIX says 8 while BITS says 4.
       Should come up with a better solution for this.
       DLK
     */
    /* Bits/lag */
    hgeti4(hdr, "BITPIX", &itmp);
    /* Bits/lag */
    hgeti4(hdr, "BITS", &(spigot->bits_per_lag));
    if (0 && spigot->bits_per_lag != itmp) {
        printf("\n  Warning!  '%s' claims both %d and %d bits/lag!\n\n",
               filename, spigot->bits_per_lag, itmp);
    }
    if (bits_per_lag_st == 0)
        bits_per_lag_st = spigot->bits_per_lag;
    /* Number of lags/sample */
    hgeti4(hdr, "NLAGS", &(spigot->lags_per_sample));
    if (numchan_st == 0)
        numchan_st = spigot->lags_per_sample;
    /* Number of lags/sample */
    hgeti4(hdr, "NAXIS1", &itmp);
    if (0 && spigot->lags_per_sample != itmp) {
        printf("\n  Warning!  '%s' claims both %d and %d lags/sample!\n\n",
               filename, spigot->lags_per_sample, itmp);
    }
    /* Number of spectra in this file */
    hgeti4(hdr, "NAXIS2", &(spigot->samples_per_file));
    /* Total duration of the file in sec */
    spigot->file_duration = spigot->samples_per_file * 1e-6 * spigot->dt_us;
    /* Total (planned) number of spectra */
    hgeti4(hdr, "SPECTRA", &(spigot->tot_num_samples));
    /* The initial spectrum number in the file  */
    hgeti4(hdr, "FRSTSPEC", &(spigot->first_spectrum));
    /* Update the MJD based on the number of spectra */
    spigot->MJD_obs += (spigot->first_spectrum * 1e-6 * spigot->dt_us) / SECPERDAY;
    /* Set the lag scaling and offset values if this is the firsttime */
    if (firsttime) {
        int NomOffset = 0, ii;
        char keyword[10], ctmp[100];
        float ftmp1, ftmp2;

        for (ii = 0; ii < spigot->lags_per_sample; ii++) {
            sprintf(keyword, "LAGD%04d", ii);
            hgets(hdr, keyword, 100, ctmp);
            sscanf(ctmp, "%f %f %f %f", lag_factor + ii, lag_offset + ii, &ftmp1,
                   &ftmp2);
        }
        firsttime = 0;
        switch ((int) (spigot->dt_us / 2.56)) {
            /* Offset depends on sampling time & BW */
        case 32:
            if (spigot->bandwidth >= 200)
                NomOffset = 32780;
            else {
                NomOffset = 61444;
            }
            break;
        case 16:
            if (spigot->bandwidth >= 200)
                NomOffset = 32792;
            else
                NomOffset = 63492;
            break;
        case 8:
            if (spigot->bandwidth >= 200)
                NomOffset = 49176;
            else
                NomOffset = 64516;
            break;
        case 4:
            if (spigot->bandwidth >= 200)
                NomOffset = 57368;
            else
                NomOffset = 65028;
            break;
        case 2:
            if (spigot->bandwidth >= 200)
                NomOffset = 61464;
            else
                NomOffset = 65284;
            break;
        case 1:
            if (spigot->bandwidth >= 200)
                NomOffset = 63512;
            else
                NomOffset = 65412;
            break;
        }
        {
            float other_fact = 0.0;

            if (spigot->bits_per_lag == 16) {
                other_fact = 16.0;
            } else if (spigot->bits_per_lag == 8) {
                other_fact = 16.0 * 256.0;
            } else if (spigot->bits_per_lag == 4) {
                other_fact = 16.0 * 4096.0;
            } else if (spigot->bits_per_lag == 2) {
                other_fact = 16.0 * 65536.0;
            }
            /* Do this once so we don't have to do it later in the decoding loop */
            for (ii = 0; ii < spigot->lags_per_sample; ii++) {
                lag_offset[ii] = NomOffset - lag_offset[ii];
                lag_factor[ii] = other_fact / lag_factor[ii];
            }
        }
    }
    free(hdr);
    if (lags == NULL) {
        /* The following should be freed sometime... */
        lags = (float *) fftwf_malloc((numchan_st + 1) * sizeof(float));
        /* Generate the FFTW plan */
        fftplan =
            fftwf_plan_r2r_1d(numchan_st + 1, lags, lags, FFTW_REDFT00,
                              FFTW_PATIENT);
    }
    return (1);
}


void SPIGOT_INFO_to_inf(SPIGOT_INFO * spigot, infodata * idata)
/* Convert a SPIGOT_INFO structure into an infodata structure */
{
    double MJD;
    char ctmp[100];
    struct passwd *pwd;

    strncpy(idata->object, spigot->object, 24);
    hours2hms(spigot->ra / 15.0, &(idata->ra_h), &(idata->ra_m), &(idata->ra_s));
    deg2dms(spigot->dec, &(idata->dec_d), &(idata->dec_m), &(idata->dec_s));
    strcpy(idata->telescope, "GBT");
    strcpy(idata->instrument, spigot->instrument);
    idata->num_chan = spigot->lags_per_sample;
    idata->dt = spigot->dt_us * 1e-6;
    MJD = UT_strings_to_MJD(spigot->date_obs, spigot->time_obs,
                            &(idata->mjd_i), &(idata->mjd_f));
    idata->mjd_f += spigot->first_spectrum * idata->dt / SECPERDAY;
    if (idata->mjd_f >= 1.0) {
        idata->mjd_f -= 1.0;
        idata->mjd_i++;
    }
    idata->N = spigot->samples_per_file;
    idata->freqband = spigot->bandwidth;
    idata->chan_wid = fabs(idata->freqband / idata->num_chan);
    idata->freq = spigot->freq_ctr - 0.5 * idata->freqband + 0.5 * idata->chan_wid;
    idata->fov = 1.2 * SOL * 3600.0 / (idata->freq * 1.0e6) / 100.0 * RADTODEG;
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    // strcpy(idata->analyzer, getlogin());
    pwd = getpwuid(geteuid());
    strcpy(idata->analyzer, pwd->pw_name);
    strncpy(idata->observer, spigot->observer, 24);
    if (spigot->summed_pols)
        sprintf(ctmp,
                "%d IF(s) were summed.  Lags are %d bit ints.",
                spigot->num_samplers, spigot->bits_per_lag);
    else
        sprintf(ctmp, "%d IF(s) were not summed.  Lags are %d bit ints.",
                spigot->num_samplers, spigot->bits_per_lag);
    sprintf(idata->notes,
            "Project ID %s, ObsID %s, Scan #%d, Date: %s %s.\n    %s\n",
            spigot->project_id, spigot->obs_id, spigot->scan_number,
            spigot->date_obs, spigot->time_obs, ctmp);
}


void print_SPIGOT_header(SPIGOT_INFO * spigot)
/* Output a SPIGOT header in human readable form */
{
    char pos_str[20];
    int h_or_d, m;
    double s;

    printf("\n           Software version = '%s'\n", spigot->software_vers);
    printf("        Header size (bytes) = %d\n", spigot->header_len);
    printf("                  Telescope = %s\n", spigot->telescope);
    printf("                 Instrument = %s\n", spigot->instrument);
    printf("                Source Name = %s\n", spigot->object);
    if (spigot->tracking)
        printf("                  Tracking? = True\n");
    else
        printf("                  Tracking? = False\n");
    printf("      Obs Date (YYYY-MM-DD) = %s\n", spigot->date_obs);
    printf("    Obs Start UT (HH:MM:SS) = %s\n", spigot->time_obs);
    printf("             MJD start time = %.12f\n", spigot->MJD_obs);
    printf("                 Project ID = %s\n", spigot->project_id);
    printf("                      ObsID = %s\n", spigot->obs_id);
    printf("                   Observer = %s\n", spigot->observer);
    printf("                Scan Number = %d\n", spigot->scan_number);
    hours2hms(spigot->ra / 15.0, &h_or_d, &m, &s);
    ra_dec_to_string(pos_str, h_or_d, m, s);
    printf("                 RA (J2000) = %s\n", pos_str);
    deg2dms(spigot->dec, &h_or_d, &m, &s);
    ra_dec_to_string(pos_str, h_or_d, m, s);
    printf("                DEC (J2000) = %s\n", pos_str);
    printf("              Azimuth (deg) = %-17.15g\n", spigot->az);
    printf("            Elevation (deg) = %-17.15g\n", spigot->el);
    printf(" Planned Obs Duration (sec) = %-17.15g\n",
           spigot->tot_num_samples * spigot->dt_us / 1e6);
    printf("                T_samp (us) = %-17.15g\n", spigot->dt_us);
    printf("         Central freq (MHz) = %-17.15g\n", spigot->freq_ctr);
    printf("      Total Bandwidth (MHz) = %-17.15g\n", spigot->bandwidth);
    printf("              Number of IFs = %d\n", spigot->num_samplers);
    printf("          Polarization type = '%s'\n", spigot->pol_type);
    printf("               Bits per lag = %d\n", spigot->bits_per_lag);
    printf("            Lags per sample = %d\n", spigot->lags_per_sample);
    printf("      First spectra in file = %d\n", spigot->first_spectrum);
    printf("           Spectra per file = %d\n", spigot->samples_per_file);
    printf("            Correlator mode = '%s'\n", spigot->corr_mode);
    printf("   Other information:\n");
    if (spigot->summed_pols)
        printf("      IFs were summed in hardware.\n");
    if (spigot->upper_sideband)
        printf("      Lags are upper sideband.\n");
    else
        printf("      Lags are lower sideband.\n");
}


void get_SPIGOT_file_info(FILE * files[], SPIGOT_INFO * spigot_files,
                          int numfiles, int usewindow,
                          float clipsig, long long *N, int *ptsperblock,
                          int *numchan, double *dt, double *T,
                          infodata * idata, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of SPIGOT rawfiles that you want to patch   */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
    long long filedatalen, calc_filedatalen, numpts;
    int ii;

    /* Allocate memory for our information structures */
    spigot = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
    idata_st = (infodata *) malloc(sizeof(infodata) * numfiles);
    /* Copy the SPIGOT_INFO structures into the static versions */
    for (ii = 0; ii < numfiles; ii++)
        spigot[ii] = spigot_files[ii];

    /* Quick hack to allow offsets of the SPIGOT center freq without recompiling */
    {
        char *envval = getenv("SPIGOT_FREQ_ADJ");
        if (envval != NULL) {
            double dblval = strtod(envval, NULL);
            if (dblval) {
                if (output)
                    printf
                        ("Offsetting band by %.4g MHz as per SPIGOT_FREQ_ADJ env variable.\n",
                         dblval);
                spigot[0].freq_ctr += dblval;
            }
        }
    }

    /* Quick hack to allow scaling of lags without recompiling */
    {
        char *envval = getenv("SPIGOT_LAG_SCALE");
        if (envval != NULL) {
            double dblval = strtod(envval, NULL);
            if (dblval) {
                if (output)
                    printf
                        ("Scaling lags by %f as per SPIGOT_LAG_SCALE env variable.\n",
                         dblval);
            }
            lag_scale_env = dblval;
        }
    }



    /* Convert the SPIGOT_INFO structures into infodata structures */
    SPIGOT_INFO_to_inf(spigot, &idata_st[0]);
    SPIGOT_INFO_to_inf(spigot, idata);

    /* Determine some important static variable values */
    bits_per_lag_st = spigot[0].bits_per_lag;
    center_freq_st = spigot[0].freq_ctr;
    *numchan = numchan_st = idata_st[0].num_chan;
    numifs_st = spigot[0].num_samplers;
    /* Override the numifs_st if the polarizations are summed in hardware */
    if (spigot[0].summed_pols) {
        numifs_st = 1;
        if (numifs_st == 2 && output)
            printf("Both IFs are present and summed.\n");
    } else {
        if (numifs_st == 1 && output)
            printf("A single IF is present.\n");
        if (numifs_st == 2 && output)
            printf("Both IFs are present.\n");
    }
    /* Flip the band if required */
    if (!spigot[0].upper_sideband) {
        decreasing_freqs_st = 1;
        if (output)
            printf("Flipping the band.");
    }
    /* We currently can't do full stokes */
    if (numifs_st > 2) {
        fprintf(stderr,
                "\n  Error:  There are more than 2 IFs present!  We can't handle this yet!\n\n");
        exit(0);
    }
    /* Are we going to clip the data? */
    if (clipsig > 0.0)
        clip_sigma_st = clipsig;
    if (usewindow) {
        usewindow_st = 1;
        if (output)
            printf("Calculated Hanning window for use.\n");
        /* Note:  Since the lags we get are only half of the lags that   */
        /* we really need to FFT in order to get spectra (i.e. the       */
        /* transform that we compute is real and even so we comute a     */
        /* DCT-I instead of an FFT), we will multiply the lags by the    */
        /* second half of the window.  The other half of the data (which */
        /* we don't store since it is redundant)  gets the 1st half of   */
        /* the window implicitly since the data wraps around.            */
        window_st = hanning_window(numchan_st);
    }
    /* Calculate the maximum number of points we can have in a */
    /* block (power of two), based on the number of samples in */
    /* each file.                                              */
    bytesperpt_st = (numchan_st * numifs_st * bits_per_lag_st) / 8;
    filedatalen = chkfilelen(files[0], 1) - spigot[0].header_len;
    calc_filedatalen =
        (spigot[0].lags_per_sample * (long long) spigot[0].samples_per_file *
         bits_per_lag_st) / 8;
    if (filedatalen != calc_filedatalen)
        fprintf(stderr,
                "\n  Warning!  The calculated (%lld) and measured (%lld) data lengths of file %d are different!\n\n",
                calc_filedatalen, filedatalen, 0);
    numpts = filedatalen / bytesperpt_st;
    if (filedatalen % bytesperpt_st)
        fprintf(stderr,
                "\n  Warning!  File %d has a non-integer number of complete samples!\n\n",
                0);
    if (numpts != spigot[0].samples_per_file) {
        fprintf(stderr,
                "\n  Warning!  The calculated (%lld) and reported (%d) number of samples in file %d are different!\n\n",
                numpts, spigot[0].samples_per_file, 0);
        spigot[0].samples_per_file = numpts;
    }
    /* Calculate the largest block size the fits evenly in each file */
    ptsperblk_st = SPIGOT_MAXPTSPERBLOCK;
    while (numpts % ptsperblk_st)
        ptsperblk_st--;
    bytesperblk_st = ptsperblk_st * bytesperpt_st;
    if (filedatalen % bytesperblk_st)
        fprintf(stderr,
                "\n  Warning!  File %d has a non-integer number of complete blocks!\n\n",
                0);
    *ptsperblock = ptsperblk_st;
    sampperblk_st = ptsperblk_st * numchan_st;
    spigot[0].num_blocks = filedatalen / bytesperblk_st;
    N_st = numpts;
    *dt = dt_st = idata_st[0].dt;
    spigot[0].file_duration = numpts * dt_st;
    spigot[0].elapsed_time = 0.0;
    spigot[0].start_block = 1.0;
    spigot[0].end_block = (double) numpts / ptsperblk_st;
    spigot[0].padding_samples = spigot[numfiles - 1].padding_samples = 0;
    /* Position the file stream at the beginning of the data */
    chkfseek(files[0], spigot[0].header_len, SEEK_SET);
    /* Now step through the other files and determine the important values */
    for (ii = 1; ii < numfiles; ii++) {
        chkfseek(files[ii], spigot[ii].header_len, SEEK_SET);
        SPIGOT_INFO_to_inf(spigot + ii, &idata_st[ii]);
        if (idata_st[ii].num_chan != numchan_st) {
            fprintf(stderr,
                    "\n  Warning!  Number of channels (file %d) is not the same!\n\n",
                    ii + 1);
        }
        if (idata_st[ii].dt != dt_st) {
            fprintf(stderr,
                    "\n  Warning!  Sample time (file %d) is not the same!\n\n",
                    ii + 1);
        }
        filedatalen = chkfilelen(files[ii], 1) - spigot[ii].header_len;
        spigot[ii].num_blocks = filedatalen / bytesperblk_st;
        numpts = spigot[ii].num_blocks * ptsperblk_st;
        if (numpts != spigot[ii].samples_per_file) {
            fprintf(stderr,
                    "\n  Warning!  The calculated and reported number of samples in file %d are different!\n\n",
                    ii);
            spigot[ii].samples_per_file = numpts;
        }
        spigot[ii].file_duration = numpts * dt_st;
        spigot[ii].MJD_obs = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
        spigot[ii].elapsed_time =
            mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                         idata_st[ii - 1].mjd_i, idata_st[ii - 1].mjd_f);
        spigot[ii - 1].padding_samples =
            (int) ((spigot[ii].elapsed_time - spigot[ii - 1].file_duration) / dt_st +
                   0.5);
        spigot[ii].elapsed_time += spigot[ii - 1].elapsed_time;
        N_st += numpts + spigot[ii - 1].padding_samples;
        spigot[ii].start_block = (double) (N_st - numpts) / ptsperblk_st + 1;
        spigot[ii].end_block = (double) (N_st) / ptsperblk_st;
    }
    spigot[numfiles - 1].padding_samples =
        ((int) ceil(spigot[numfiles - 1].end_block) * ptsperblk_st - N_st);
    N_st += spigot[numfiles - 1].padding_samples;
    *N = N_st;
    *T = T_st = N_st * dt_st;
    currentfile = currentblock = 0;
    if (output) {
        printf("   Number of files = %d\n", numfiles);
        printf("      Points/block = %d\n", ptsperblk_st);
        printf("   Num of channels = %d\n", numchan_st);
        printf("       Center freq = %.8g\n", center_freq_st);
        printf("  Total points (N) = %lld\n", N_st);
        printf("  Sample time (dt) = %-14.14g\n", dt_st);
        printf("    Total time (s) = %-14.14g\n", T_st);
        printf(" Header length (B) = %d\n", spigot[0].header_len);
        printf
            ("File  Start Block    Last Block     Points      Elapsed (s)      Time (s)            MJD           Padding\n");
        printf
            ("----  ------------  ------------  ----------  --------------  --------------  ------------------  ----------\n");
        for (ii = 0; ii < numfiles; ii++)
            printf
                ("%2d    %12.11g  %12.11g  %10d  %14.13g  %14.13g  %17.12f  %10d\n",
                 ii + 1, spigot[ii].start_block, spigot[ii].end_block,
                 spigot[ii].samples_per_file, spigot[ii].elapsed_time,
                 spigot[ii].file_duration, spigot[ii].MJD_obs,
                 spigot[ii].padding_samples);
        printf("\n");
    }
}


void SPIGOT_update_infodata(int numfiles, infodata * idata)
/* Update the onoff bins section in case we used multiple files */
{
    int ii, index = 2;

    idata->N = N_st;
    if (numfiles == 1 && spigot[0].padding_samples == 0) {
        idata->numonoff = 0;
        return;
    }
    /* Determine the topocentric onoff bins */
    idata->numonoff = 1;
    idata->onoff[0] = 0.0;
    idata->onoff[1] = spigot[0].samples_per_file - 1.0;
    for (ii = 1; ii < numfiles; ii++) {
        if (spigot[ii - 1].padding_samples) {
            idata->onoff[index] =
                idata->onoff[index - 1] + spigot[ii - 1].padding_samples;
            idata->onoff[index + 1] =
                idata->onoff[index] + spigot[ii].samples_per_file;
            idata->numonoff++;
            index += 2;
        } else {
            idata->onoff[index - 1] += spigot[ii].samples_per_file;
        }
    }
    if (spigot[numfiles - 1].padding_samples) {
        idata->onoff[index] =
            idata->onoff[index - 1] + spigot[numfiles - 1].padding_samples;
        idata->onoff[index + 1] = idata->onoff[index];
        idata->numonoff++;
    }
}
