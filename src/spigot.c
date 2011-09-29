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
static int need_byteswap_st = 0, sampperblk_st, usewindow_st = 0, vanvleck_st = 0;
static int currentfile, currentblock, bufferpts = 0, padnum = 0, shiftbuffer = 1;
static double T_st, dt_st, center_freq_st, *window_st = NULL;
static float clip_sigma_st = 0.0, *lags = NULL;
static float lag_factor[SPIGOT_MAXLAGLEN], lag_offset[SPIGOT_MAXLAGLEN];
static unsigned char padvals[SPIGOT_MAXLAGLEN], newpadvals[SPIGOT_MAXLAGLEN];
static unsigned char databuffer[2 * SPIGOT_MAXDATLEN], padval = 128;
static unsigned char lagbuffer[SPIGOT_MAXLAGLEN];
static int using_MPI = 0;
static float lag_scale_env = 1.0;

double slaCldj(int iy, int im, int id, int *j);
static double inv_cerf(double input);
static void vanvleck3lev(float *rho, int npts);
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

void get_SPIGOT_static(int *bytesperpt, int *bytesperblk, int *numifs,
                       float *clip_sigma)
/* This is used by the MPI-based mpiprepsubband code to insure that each */
/* processor can access (correctly) the required static variables.       */
{
   *bytesperpt = bytesperpt_st;
   *bytesperblk = bytesperblk_st;
   *numifs = numifs_st;
   *clip_sigma = clip_sigma_st;
}

void set_SPIGOT_static(int ptsperblk, int bytesperpt, int bytesperblk,
                       int numchan, int numifs, float clip_sigma, double dt)
/* This is used by the MPI-based mpiprepsubband code to insure that each */
/* processor can access (correctly) the required static variables.       */
{
   using_MPI = 1;
   currentblock = 0;
   ptsperblk_st = ptsperblk;
   bytesperpt_st = bytesperpt;
   bytesperblk_st = bytesperblk;
   bits_per_lag_st = (bytesperpt * 8) / numchan;
   numchan_st = numchan;
   numifs_st = numifs;
   sampperblk_st = ptsperblk_st * numchan_st;
   clip_sigma_st = clip_sigma;
   dt_st = dt;
}

void set_SPIGOT_padvals(float *fpadvals, int good_padvals)
/* If reasonable padding values are available from rifind,       */
/* use these as the padding values, otherwise, set them all      */
/* to the value of padval defined in the static variables above. */
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
   if (fabs(dtmp2 > 1e-6)) {
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
          fftwf_plan_r2r_1d(numchan_st + 1, lags, lags, FFTW_REDFT00, FFTW_PATIENT);
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
   sprintf(idata->notes, "Project ID %s, ObsID %s, Scan #%d, Date: %s %s.\n    %s\n",
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
         fprintf(stderr, "\n  Warning!  Sample time (file %d) is not the same!\n\n",
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
      spigot[ii].elapsed_time = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                                             idata_st[ii - 1].mjd_i,
                                             idata_st[ii - 1].mjd_f);
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
         printf("%2d    %12.11g  %12.11g  %10d  %14.13g  %14.13g  %17.12f  %10d\n",
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
         idata->onoff[index + 1] = idata->onoff[index] + spigot[ii].samples_per_file;
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


int skip_to_SPIGOT_rec(FILE * infiles[], int numfiles, int rec)
/* This routine skips to the record 'rec' in the input files   */
/* *infiles.  *infiles contain data from the SPIGOT at GBT     */
/* Returns the record skipped to.                              */
{
   double floor_blk;
   int filenum = 0;

   if (rec < spigot[0].start_block)
      rec += (spigot[0].start_block - 1);
   if (rec > 0 && rec < spigot[numfiles - 1].end_block) {

      /* Find which file we need */
      while (rec > spigot[filenum].end_block)
         filenum++;

      currentblock = rec - 1;
      shiftbuffer = 1;
      floor_blk = floor(spigot[filenum].start_block);

      /* Set the data buffer to all padding just in case */
      memset(databuffer, padval, 2 * SPIGOT_MAXDATLEN);

      /* Warning:  I'm not sure if the following is correct. */
      /*   If really needs accurate testing to see if my     */
      /*   offsets are correct.  Bottom line, don't trust    */
      /*   a TOA determined using the following!             */

      if (rec < spigot[filenum].start_block) {  /* Padding region */
         currentfile = filenum - 1;
         chkfileseek(infiles[currentfile], 0, 1, SEEK_END);
         bufferpts = spigot[currentfile].padding_samples % ptsperblk_st;
         padnum = ptsperblk_st * (rec - spigot[currentfile].end_block - 1);
         /*
            printf("Padding:  currentfile = %d  bufferpts = %d  padnum = %d\n",
            currentfile, bufferpts, padnum);
          */
      } else {                  /* Data region */
         currentfile = filenum;
         chkfileseek(infiles[currentfile],
                     rec - spigot[filenum].start_block, bytesperblk_st, SEEK_CUR);
         bufferpts = (int) ((spigot[filenum].start_block - floor_blk) *
                            ptsperblk_st + 0.5);
         padnum = 0;
         /*
            printf("Data:  currentfile = %d  bufferpts = %d  padnum = %d\n",
            currentfile, bufferpts, padnum);
          */
      }

   } else {
      printf("\n rec = %d out of range in skip_to_SPIGOT_rec()\n", rec);
      exit(1);
   }
   return rec;
}


int read_SPIGOT_rawblock(FILE * infiles[], int numfiles,
                         unsigned char *data, int *padding, IFs ifs)
/* This routine reads a single record from the input      */
/* files *infiles which contain lag data from the SPIGOT. */
/* A SPIGOT record is ptsperblk_st*numchan_st*#bits long. */
/* *data must be sampperblk_st bytes long.                */
/* If padding is returned as 1, then padding was added    */
/* and statistics should not be calculated.               */
{
   int offset = 0, numtopad = 0, ii;
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

   if (currentfile >= numfiles)
      return 0;

   /* First, attempt to read data from the current file */

   if (chkfread(lagbuffer, bytesperblk_st, 1, infiles[currentfile])) {  /* Got data */
      /* See if we need to byte-swap and if so, doit */
      if (need_byteswap_st) {
         if (bits_per_lag_st == 16) {
            unsigned short *sptr = (unsigned short *) lagbuffer;
            for (ii = 0; ii < sampperblk_st; ii++, sptr++)
               *sptr = swap_ushort(*sptr);
         }
      }
      /* Convert from Correlator Lags to Filterbank Powers */
      for (ii = 0; ii < ptsperblk_st; ii++)
         convert_SPIGOT_point(lagbuffer + ii * bytesperpt_st,
                              dataptr + ii * numchan_st, ifs, 1.0);

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
      if (feof(infiles[currentfile])) { /* End of file? */
         numtopad = spigot[currentfile].padding_samples - padnum;
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
               if (padnum == spigot[currentfile].padding_samples) {
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
               return read_SPIGOT_rawblock(infiles, numfiles, data, &pad, ifs);
            }
         } else {               /* No padding needed.  Try reading the next file */
            currentfile++;
            shiftbuffer = 0;
            return read_SPIGOT_rawblock(infiles, numfiles, data, padding, ifs);
         }
      } else {
         printf("\nProblem reading record from SPIGOT data file:\n");
         printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                currentfile, currentblock);
         exit(1);
      }
   }
}


int read_SPIGOT_rawblocks(FILE * infiles[], int numfiles,
                          unsigned char rawdata[], int numblocks,
                          int *padding, IFs ifs)
/* This routine reads numblocks SPIGOT records from the input */
/* files *infiles.  The 8-bit filterbank data is returned     */
/* in rawdata which must have a size of numblocks*            */
/* sampperblk_st.  The number  of blocks read is returned.    */
/* If padding is returned as 1, then padding was added        */
/* and statistics should not be calculated                    */
{
   int ii, retval = 0, pad, numpad = 0;

   *padding = 0;
   for (ii = 0; ii < numblocks; ii++) {
      retval += read_SPIGOT_rawblock(infiles, numfiles,
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


int read_SPIGOT(FILE * infiles[], int numfiles, float *data,
                int numpts, double *dispdelays, int *padding,
                int *maskchans, int *nummasked, mask * obsmask, IFs ifs)
/* This routine reads numpts from the SPIGOT raw input */
/* files *infiles.  These files contain raw correlator */
/* data from the SPIGOT at GBT.  Time delays           */
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
         printf("numpts must be a multiple of %d in read_SPIGOT()!\n", ptsperblk_st);
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
         numread = read_SPIGOT_rawblocks(infiles, numfiles, currentdata,
                                         numblocks, padding, ifs);
         if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
         }

         /* Only use the recently measured padding if all the channels aren't masked */
         if ((clip_sigma_st > 0.0) && 
             !(mask && (*nummasked == -1)) &&
             (padvals != newpadvals))
             memcpy(padvals, newpadvals, SPIGOT_MAXLAGLEN);

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


void get_SPIGOT_channel(int channum, float chandat[],
                        unsigned char rawdata[], int numblocks)
/* Return the values for channel 'channum' of a block of         */
/* 'numblocks' raw SPIGOT data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                  */
/* read_SPIGOT_rawblocks(), and 'chandat' must have at least     */
/* 'numblocks' * 'ptsperblk_st' spaces.                          */
/* Channel 0 is assumed to be the lowest freq channel.          */
{
   int ii, jj, ptsperchan;

   if (channum > numchan_st * numifs_st || channum < 0) {
      printf("\nchannum = %d is out of range in get_SPIGOT_channel()!\n\n", channum);
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


int prep_SPIGOT_subbands(unsigned char *rawdata, float *data,
                         double *dispdelays, int numsubbands,
                         int transpose, int *maskchans,
                         int *nummasked, mask * obsmask)
/* This routine preps a block from the SPIGOT system.  The routine uses   */
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
   static unsigned char rawdata1[SPIGOT_MAXDATLEN], rawdata2[SPIGOT_MAXDATLEN];
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
       memcpy(padvals, newpadvals, SPIGOT_MAXLAGLEN);

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


int read_SPIGOT_subbands(FILE * infiles[], int numfiles, float *data,
                         double *dispdelays, int numsubbands,
                         int transpose, int *padding,
                         int *maskchans, int *nummasked, mask * obsmask, IFs ifs)
/* This routine reads a record from the input files *infiles[]   */
/* which contain data from the SPIGOT system.  The routine uses  */
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
   static unsigned char rawdata[SPIGOT_MAXDATLEN];

   if (firsttime) {
      if (!read_SPIGOT_rawblock(infiles, numfiles, rawdata, padding, ifs)) {
         printf("Problem reading the raw SPIGOT data file.\n\n");
         return 0;
      }
      if (0 != prep_SPIGOT_subbands(rawdata, data, dispdelays, numsubbands,
                                    transpose, maskchans, nummasked, obsmask)) {
         printf("Problem initializing prep_SPIGOT_subbands()\n\n");
         return 0;
      }
      firsttime = 0;
   }
   if (!read_SPIGOT_rawblock(infiles, numfiles, rawdata, padding, ifs)) {
      /* printf("Problem reading the raw SPIGOT data file.\n\n"); */
      if (!memfreed) {
         fftwf_destroy_plan(fftplan);
         fftwf_free(lags);
         free(window_st);
         memfreed = 1;
      }
      return 0;
   }
   return prep_SPIGOT_subbands(rawdata, data, dispdelays, numsubbands,
                               transpose, maskchans, nummasked, obsmask);
}

void scale_rawlags(void *rawdata, int index)
/* Scale the raw lags so that they are "calibrated" */
{
   int ii, lag;
   //  To turn off the wrap correction, set hilag = -1
   static int hilag = -1, last_lags[20], wrapvals[20];
   static double ravg = 0.0, firsttime = 1;

   if (firsttime) {
      if (bits_per_lag_st==16 && hilag > 0)
         hilag = 1;  // For 16-bit data only check the zero lag for wrap issues
      for (ii = 0; ii < hilag; ii++) {
         last_lags[ii] = 0;
         wrapvals[ii] = 0;
      }
      firsttime = 0;
   }
   /* Fill lag array with scaled CFs */
   /* Note:  The 2 and 4 bit options will almost certainly need to be fixed */
   /*        since I don't know the precise data layout. SMR 15Aug2003      */
   /*        Also, I have no idea how multiple IFs/RFs are stored...        */
   if (bits_per_lag_st == 16) {
      short *data = (short *) rawdata;
      for (ii = 0; ii < numchan_st; ii++) {
         lag = (int) data[ii + index];
         /* Attempt to lags that have either over-flowed or under-flowed */
         if (ii < hilag) {
            lag += wrapvals[ii];
            if (abs(lag - last_lags[ii]) > 40000) {   // 40K is quite arbitrary...
               if (lag < last_lags[ii]) {
                  // Prevent the wrapping from skrocketting due to RFI
                  if (wrapvals[ii] < 2*65536) {
                     wrapvals[ii] += 65536;
                     lag += 65536;
                  }
               } else {
                  // Prevent the wrapping from skrocketting due to RFI
                  if (wrapvals[ii] > -2*65536) {
                     wrapvals[ii] -= 65536;
                     lag -= 65536;
                  }
               }
            }
            last_lags[ii] = lag;
         }
         lags[ii] = (lag * lag_factor[ii] + lag_offset[ii]) / lag_scale_env;
      }
   } else if (bits_per_lag_st == 8) {
      char *data = (char *) rawdata;
      for (ii = 0; ii < numchan_st; ii++) {
         lag = (int) data[ii + index];
         /* Attempt to lags that have either over-flowed or under-flowed */
         if (ii < hilag) {
            lag += wrapvals[ii];
            if (abs(lag - last_lags[ii]) > 170) {   // 170 is quite arbitrary...
		 if (lag < last_lags[ii]) {
                  // Prevent the wrapping from skrocketting due to RFI
                  if (wrapvals[ii] < 1*256) {
                     wrapvals[ii] += 256;
                     lag += 256;
                  }
               } else {
                  // Prevent the wrapping from skrocketting due to RFI
                  if (wrapvals[ii] > -1*256) {
                     wrapvals[ii] -= 256;
                     lag -= 256;
                  }
               }
            }
            last_lags[ii] = lag;
         }
         lags[ii] = (lag * lag_factor[ii] + lag_offset[ii]) / lag_scale_env;
      }
   } else if (bits_per_lag_st == 4) {
      static int wrapval;
      int jj, lag0;
      char tmplag;
      unsigned char *data = (unsigned char *) rawdata, byte;
      /* Attempt to fix zerolags that have either over-flowed or under-flowed */
      {
         byte = data[index / 2];
         lag0 = ((byte & 0x70) >> 4) - ((byte & 0x80) >> 4);
         lag0 += wrapval;

         if (fabs(lag0 - ravg) > 12) {  // 12 is quite arbitrary
            if (lag0 < ravg) {
               if (wrapval <= 32)
                  wrapval += 16;        // Prevent the wrapval from skyrocketting
               lag0 += 16;
            } else {
               if (wrapval >= -32)
                  wrapval -= 16;        // Prevent the wrapval from skyrocketting
               lag0 -= 16;
            }
            // printf("wrapping at %d.  wrapval = %d\n", ii, wrapval);
         }
         ravg = 0.1 * ravg + 0.9 * lag0;
         lags[0] = lag0 * lag_factor[0] + lag_offset[0];
         tmplag = (byte & 0x07) - (byte & 0x08);
         lags[1] = tmplag * lag_factor[1] + lag_offset[1];
      }
      for (ii = 1, jj = 2; ii < numchan_st / 2; ii++) {
         byte = data[ii + index / 2];
         /* Treat the 4 msbs as a twos-complement 4-bit int */
         tmplag = ((byte & 0x70) >> 4) - ((byte & 0x80) >> 4);
         lags[jj] = tmplag * lag_factor[jj] + lag_offset[jj];
         jj++;
         /* Treat the 4 lsbs as a twos-complement 4-bit int */
         tmplag = (byte & 0x07) - (byte & 0x08);
         lags[jj] = tmplag * lag_factor[jj] + lag_offset[jj];
         jj++;
      }
   } else if (bits_per_lag_st == 2) {   /* The following is _not_ correct yet ! */
      int jj;
      unsigned char *data = (unsigned char *) rawdata, byte;
      for (ii = 0, jj = 0; ii < numchan_st / 4; ii++) {
         byte = data[ii + index / 4];
         lags[jj] = (byte >> 6) * lag_factor[jj] + lag_offset[jj];
         jj++;
         lags[jj] = ((byte >> 4) & 0x03) * lag_factor[jj] + lag_offset[jj];
         jj++;
         lags[jj] = ((byte >> 2) & 0x03) * lag_factor[jj] + lag_offset[jj];
         jj++;
         lags[jj] = (byte & 0x03) * lag_factor[jj] + lag_offset[jj];
         jj++;
      }
   }
}


void get_calibrated_lags(void *rawlags, float *calibrated_lags)
/* Return the calibrated lags for a single sample of raw lags */
{
    int ii;
    
    /* Scale the raw lags */
    scale_rawlags(rawlags, 0);
    
    /* Apply Van Vleck Corrections to the Lags */
    if (vanvleck_st) {
        double power;
        
        /* Calculate power */
        //power = inv_cerf(lags[0]);
        //power = 0.1872721836 / (power * power);
        power = sqrt(2.0)*inv_cerf(lags[0]);
        power = 1.0/(power*power);
        vanvleck3lev(lags, numchan_st);
        for (ii=0; ii<numchan_st; ii++)
            lags[ii] *= power;
    }

    /* Window the Lags */
    if (usewindow_st)
        for (ii = 0; ii < numchan_st; ii++)
            lags[ii] *= window_st[ii];
    
    /* Copy the lags into the output array */
    for (ii = 0; ii < numchan_st; ii++)
        calibrated_lags[ii] = lags[ii];
}


void convert_SPIGOT_point(void *rawdata, unsigned char *bytes,
                          IFs ifs, float lag_scaling)
/* This routine converts a single point of SPIGOT lags   */
/* into a filterbank style array of bytes.               */
{
   int ii, hiclipcount = 0, loclipcount = 0, ifnum = 0, index = 0;
   float *templags = NULL;
   double pfact;
   //double power;
   static int counter = 0;
   static double scale_min = 9e19, scale_max = -9e19;

   if (ifs == IF0)
      index = 0;
   else if (ifs == IF1)
      index = numchan_st;

   /* Loop over the IFs */
   for (ifnum = 0; ifnum < numifs_st; ifnum++, index += numchan_st) {

      /* Scale the raw lags */
      scale_rawlags(rawdata, index);

      /* Apply the user-specified scaling */
      if (lag_scaling != 1.0) {
          if (counter==0)
              printf("Scaling the lags by %f\n", lag_scaling);
          //lags[0] += lag_scaling;
          // Should we scale all the lags?
          for (ii = 0; ii < numchan_st; ii++)
              lags[ii] *= lag_scaling;
      }

#if 0
      printf("\n");
      for (ii = 0; ii < numchan_st; ii++)
         printf("%d  %.7g\n", ii, lags[ii]);
      printf("\n");
      exit(0);
#endif

      /* Apply Van Vleck Corrections to the Lags */
      if (vanvleck_st) {
          double power;
          
          if ((lag_scale_env != 1.0) && (counter % 10000 == 0))
              printf("For Van Vleck correction:  power = %.4f\n", lags[0]);
          /* Calculate power */
          //power = inv_cerf(lags[0]);
          //power = 0.1872721836 / (power * power);
          power = sqrt(2.0)*inv_cerf(lags[0]);
          power = 1.0/(power*power);
          vanvleck3lev(lags, numchan_st);
          for (ii=0; ii<numchan_st; ii++)
              lags[ii] *= power;
      }

      /* Window the Lags */
      if (usewindow_st)
         for (ii = 0; ii < numchan_st; ii++)
            lags[ii] *= window_st[ii];

#if 0
      printf("\n");
      for (ii = 0; ii < numchan_st; ii++)
         printf("%d  %.7g\n", ii, lags[ii]);
      printf("\n");
#endif

      /* FFT the ACF lags (which are real and even) -> real and even FFT */
      /* Set the missing lag as per Carl Heiles, PASP paper */
      lags[numchan_st] = lags[numchan_st - 1];

# if 0
      if (counter > 200000) { 
         static int sctr=0, nsum=1000;
         static float spectrum[1024];
         
         if (sctr==0)
            for (ii = 0; ii < numchan_st; ii++) 
               spectrum[ii] = 0.0;
         if (sctr<nsum){
            for (ii = 0; ii < numchan_st; ii++)
               spectrum[ii] += lags[ii];
         } else {
            for (ii = 0; ii < numchan_st; ii++)
               printf("%.7g\n", spectrum[ii]/nsum);
            exit(0);
         }
         sctr++;
      }
#endif

      fftwf_execute(fftplan);

#if 0
      printf("\n");
      for (ii = 0; ii < numchan_st; ii++)
         printf("%d  %.7g\n", ii, lags[ii]);
      printf("\n");
      exit(0);
#endif

      /* Determine some simple statistics of the spectra */
      if (counter % 10240 == 0) {
         static double min = 9e9, max = -9e9;
         double avg = 0, var = 0, max2 = -9e9, max3 = -9e9;

         avg_var(lags, numchan_st, &avg, &var);
         for (ii = 0; ii < numchan_st; ii++) {
            if (lags[ii] < min)
               min = lags[ii];
            if (lags[ii] > max) {
               max3 = max2;  //  3nd highest value for RFI avoidance
               max2 = max;   //  2nd highest value for RFI avoidance
               max = lags[ii];
            }
         }
         if (counter == 0) {
            double range;

            // 2048-lag 8-bit data with bad calibration hack
            if (numchan_st==2048) {
               double left_lo = 0, right_lo = 0;
               
               // printf("Using modified spectra scaling for 8-bit 2048-lag data...\n");
               avg_var(lags+1, 5, &left_lo, &var);
               avg_var(lags+numchan_st-7, 5, &right_lo, &var);
               avg_var(lags+1, numchan_st-1, &avg, &var);
               if (left_lo < right_lo) min = left_lo;
               else  min = right_lo;
               //fprintf(stderr, "  freq[0] = %f  'lo'-freq avg = %f  'hi'-freq avg = %f\n", lags[0], left_lo, right_lo);
               //fprintf(stderr, "  freqs_mean = %f  freqs_std = %f  freq_max1,2,3 = %f,%f,%f\n", avg, sqrt(var), max, max2, max3);
            }
            range = max3 - min;
            //  The maximum power tends to vary more (due to RFI) than the
            //  minimum power.  (I think)
            scale_max = max3 + 1.5 * range;
            scale_min = min  - 0.5 * range;
            //if (numchan_st==2048)
            // fprintf(stderr, "  scale_min = %f  scale_max = %f\n", scale_min, scale_max);
            /* Check to see if we are overriding the scaling values */
            {
               char *envval;
               envval = getenv("SPIGOT_SCALING_MIN");
               if (envval != NULL) {
                  double dblval = strtod(envval, NULL);
                  if (dblval) {
                     fprintf(stderr, "Measured scaling min = %20.15f\n", scale_min);
                     scale_min = dblval;
                     fprintf(stderr, "    User scaling min = %20.15f\n", scale_min);
                  }
               }
               envval = getenv("SPIGOT_SCALING_MAX");
               if (envval != NULL) {
                  double dblval = strtod(envval, NULL);
                  if (dblval) {
                     fprintf(stderr, "Measured scaling max = %20.15f\n", scale_max);
                     scale_max = dblval;
                     fprintf(stderr, "    User scaling max = %20.15f\n", scale_max);
                  }
               }
            }
         }
#if 0
         printf("counter = %d  avg = %.3g  stdev = %.3g  min = %.3g  max = %.3g\n",
                counter, avg, sqrt(var), min, max);
#endif
      }
      counter++;

      /* Reverse band if it needs it */
      if (decreasing_freqs_st) {
         float tempzz = 0.0, *loptr, *hiptr;
         loptr = lags + 0;
         hiptr = lags + numchan_st - 1;
         for (ii = 0; ii < numchan_st / 2; ii++, loptr++, hiptr--) {
            SWAP(*loptr, *hiptr);
         }
      }

      if (numifs_st == 2 && ifs == SUMIFS) {
         if (ifnum == 0) {
            templags = gen_fvect(numchan_st);
            /* Copy the unscaled values to the templag array */
            for (ii = 0; ii < numchan_st; ii++)
               templags[ii] = lags[ii];
         } else {
            /* Sum the unscaled IFs */
            for (ii = 0; ii < numchan_st; ii++)
               lags[ii] += templags[ii];
            vect_free(templags);
         }
      }
   }

   /* Scale and pack the powers into unsigned chars */
   pfact = 255.0 / (scale_max - scale_min);
   for (ii = 0; ii < numchan_st; ii++) {
      double templag;
      //templag = (lags[ii] > scale_max) ? scale_max : lags[ii];
      //templag = (templag < scale_min) ? scale_min : templag;
      if (lags[ii] > scale_max) {
         templag = scale_max;
         hiclipcount++;
      } else {
         templag = lags[ii];
      }
      if (templag < scale_min) {
         templag = scale_min;
         loclipcount++;
      }
      bytes[ii] = (unsigned char) ((templag - scale_min) * pfact + 0.5);
   }
   if (0 && (hiclipcount > 100 || loclipcount > 100))
      fprintf(stderr, "Warning:  For sample %d, (lo, hi) clipcounts = %d, %d\n",
              counter, loclipcount, hiclipcount);
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
