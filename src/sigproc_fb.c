#include "presto.h"
#include "mask.h"
#include "sigproc_fb.h"

#define MAXNUMCHAN 8192
#define BLOCKLEN   512

/* All of the following have an _st to indicate static */
static sigprocfb *fb_st;
static infodata *idata_st;
static long long *numpts_st, *padpts_st, N_st;
static int *numblks_st, need_byteswap_st = 0, sampperblk_st;
static int numchan_st, ptsperblk_st, bytesperpt_st = 1, bytesperblk_st;
static double *times_st, *mjds_st, Tdiam=100.0;
static double *elapsed_st, T_st, dt_st;
static double *startblk_st, *endblk_st;
static unsigned char padvals[MAXNUMCHAN], padval = 128;
static int rawdatabuffer[MAXNUMCHAN * BLOCKLEN];
static unsigned char databuffer[MAXNUMCHAN * BLOCKLEN];
static int currentfile, currentblock;
static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
static float clip_sigma_st = 0.0;
static int using_MPI = 0;
static int ptsperbyte_st = 1;

/* Note:  Much of this has been ripped out of SIGPROC      */
/* and then slightly modified.  Thanks Dunc!               */

static void send_string(char *string, FILE * outfile)
{
   int len;
   len = strlen(string);
   chkfwrite(&len, sizeof(int), 1, outfile);
   chkfwrite(string, sizeof(char), len, outfile);
}

static void get_string(FILE * inputfile, int *nbytes, char string[])
{
   int nchar;
   strcpy(string, "ERROR");
   chkfread(&nchar, sizeof(int), 1, inputfile);
   if (nchar > 80 || nchar < 1)
      return;
   *nbytes = sizeof(int);
   chkfread(string, nchar, 1, inputfile);
   string[nchar] = '\0';
   *nbytes += nchar;
}

static int strings_equal(char *string1, char *string2)
{
   if (!strcmp(string1, string2)) {
      return 1;
   } else {
      return 0;
   }
}

static void send_double(char *name, double double_precision, FILE * outfile)
{
   send_string(name, outfile);
   chkfwrite(&double_precision, sizeof(double), 1, outfile);
}

static void send_int(char *name, int integer, FILE * outfile)
{
   send_string(name, outfile);
   chkfwrite(&integer, sizeof(int), 1, outfile);
}

static void send_coords(double raj, double dej, double az, double za, FILE * outfile)
{
   if ((raj != 0.0) || (raj != -1.0))
      send_double("src_raj", raj, outfile);
   if ((dej != 0.0) || (dej != -1.0))
      send_double("src_dej", dej, outfile);
   if ((az != 0.0) || (az != -1.0))
      send_double("az_start", az, outfile);
   if ((za != 0.0) || (za != -1.0))
      send_double("za_start", za, outfile);
}

static char *telescope_name(int telescope_id)
{
   char *telescope, string[80];
   switch (telescope_id) {
   case 0:
      strcpy(string, "Fake");
      break;
   case 1:
      strcpy(string, "Arecibo");
      Tdiam = 305.0;
      break;
   case 2:
      strcpy(string, "Ooty");
      break;
   case 3:
      strcpy(string, "Nancay");
      break;
   case 4:
      strcpy(string, "Parkes");
      Tdiam = 64.0;
      break;
   case 5:
      strcpy(string, "Jodrell");
      Tdiam = 76.0;
      break;
   case 6:
      strcpy(string, "GBT");
      Tdiam = 100.0;
      break;
   case 7:
      strcpy(string, "GMRT");
      Tdiam = 45.0;  // possibly not right if using phased array
      break;
   case 8:
      strcpy(string, "Effelsberg");
      Tdiam = 100.0;
      break;
   case 9:
      strcpy(string, "ATA");
      break;
   case 10:
      strcpy(string, "UTR-2");
      break;
   case 11:
      strcpy(string, "LOFAR");
      Tdiam = 999.0; // depends on configuration
      break;
   default:
      strcpy(string, "???????");
      break;
   }
   telescope = (char *) calloc(strlen(string) + 1, 1);
   strcpy(telescope, string);
   return telescope;
}

static char *backend_name(int machine_id)
{
   char *backend, string[80];
   switch (machine_id) {
   case 0:
      strcpy(string, "FAKE");
      break;
   case 1:
      strcpy(string, "PSPM");
      break;
   case 2:
      strcpy(string, "WAPP");
      break;
   case 3:
      strcpy(string, "AOFTM");
      break;
   case 4:
      strcpy(string, "BPP");
      break;
   case 5:
      strcpy(string, "OOTY");
      break;
   case 6:
      strcpy(string, "SCAMP");
      break;
   case 7:
      strcpy(string, "SPIGOT");
      break;
   case 11:
      strcpy(string, "BG/P");
      break;
   default:
      strcpy(string, "????");
      break;
   }
   backend = (char *) calloc(strlen(string) + 1, 1);
   strcpy(backend, string);
   return backend;
}

/*
static char *data_category(int data_type)
{
  char *datatype, string[80];
  switch (data_type) {
  case 0:
    strcpy(string,"raw data");
    break;
  case 1:
    strcpy(string,"filterbank");
    break;
  case 2:
    strcpy(string,"time series");
    break;
  case 3:
    strcpy(string,"pulse profiles");
    break;
  case 4:
    strcpy(string,"amplitude spectrum");
    break;
  case 5:
    strcpy(string,"complex spectrum");
    break;
  case 6:
    strcpy(string,"dedispersed subbands");
    break;
  default:
    strcpy(string,"unknown!");
    break;
  }
  datatype=(char *)calloc(strlen(string)+1, 1);
  strcpy(datatype,string);
  return datatype;
}
*/

void write_filterbank_header(sigprocfb * fb, FILE * outfile)
{
   int ii, jj;

   if (fb->machine_id != 0) {
      send_string("HEADER_START", outfile);
      send_string("rawdatafile", outfile);
      send_string(fb->inpfile, outfile);
      if (!strings_equal(fb->source_name, "")) {
         send_string("source_name", outfile);
         send_string(fb->source_name, outfile);
      }
      send_int("machine_id", fb->machine_id, outfile);
      send_int("telescope_id", fb->telescope_id, outfile);
      send_coords(fb->src_raj, fb->src_dej, fb->az_start, fb->za_start, outfile);
      send_int("data_type", 1, outfile);        /* filterbank data */
      send_double("fch1", fb->fch1, outfile);
      send_double("foff", fb->foff, outfile);
      send_int("nchans", fb->nchans, outfile);
      send_int("nbits", fb->nbits, outfile);
      send_double("tstart", fb->tstart, outfile);
      send_double("tsamp", fb->tsamp, outfile);
      if (fb->sumifs) {
         send_int("nifs", 1, outfile);
      } else {
         jj = 0;
         for (ii = 1; ii <= fb->nifs; ii++)
            if (fb->ifstream[ii - 1] == 'Y')
               jj++;
         if (jj == 0)
            printf("\nNo valid IF streams selected!\n\n");
         send_int("nifs", jj, outfile);
      }
      send_string("HEADER_END", outfile);
   }
}

/* attempt to read in the general header info from a pulsar data file */
int read_filterbank_header(sigprocfb * fb, FILE * inputfile)
{
   char string[80], message[80];
   int itmp, nbytes = 0, totalbytes;
   int expecting_rawdatafile = 0, expecting_source_name = 0;
   int barycentric,pulsarcentric;
   /* try to read in the first line of the header */
   get_string(inputfile, &nbytes, string);
   if (!strings_equal(string, "HEADER_START")) {
      /* the data file is not in standard format, rewind and return */
      rewind(inputfile);
      return 0;
   }
   /* store total number of bytes read so far */
   totalbytes = nbytes;

   /* loop over and read remaining header lines until HEADER_END reached */
   while (1) {
      get_string(inputfile, &nbytes, string);
      if (strings_equal(string, "HEADER_END"))
         break;
      totalbytes += nbytes;
      if (strings_equal(string, "rawdatafile")) {
         expecting_rawdatafile = 1;
      } else if (strings_equal(string, "source_name")) {
         expecting_source_name = 1;
      } else if (strings_equal(string, "az_start")) {
         chkfread(&(fb->az_start), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "za_start")) {
         chkfread(&(fb->za_start), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "src_raj")) {
         chkfread(&(fb->src_raj), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "src_dej")) {
         chkfread(&(fb->src_dej), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "tstart")) {
         chkfread(&(fb->tstart), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "tsamp")) {
         chkfread(&(fb->tsamp), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "fch1")) {
         chkfread(&(fb->fch1), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "foff")) {
         chkfread(&(fb->foff), sizeof(double), 1, inputfile);
         totalbytes += sizeof(double);
      } else if (strings_equal(string, "nchans")) {
         chkfread(&(fb->nchans), sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "telescope_id")) {
         chkfread(&(fb->telescope_id), sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "machine_id")) {
         chkfread(&(fb->machine_id), sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "data_type")) {
         chkfread(&itmp, sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "nbits")) {
         chkfread(&(fb->nbits), sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string,"barycentric")) {
         chkfread(&barycentric,sizeof(int),1,inputfile);
	 totalbytes+=sizeof(int);
      } else if (strings_equal(string,"pulsarcentric")) {
         chkfread(&pulsarcentric,sizeof(int),1,inputfile);
         totalbytes+=sizeof(int);
      } else if (strings_equal(string, "nsamples")) {
         /* read this one only for backwards compatibility */
         chkfread(&itmp, sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "nifs")) {
         chkfread(&(fb->nifs), sizeof(int), 1, inputfile);
         totalbytes += sizeof(int);
      } else if (strings_equal(string, "nbeams")) {
         chkfread(&(fb->nbeams),sizeof(int),1,inputfile);
         totalbytes += sizeof(int);
      } else if(strings_equal(string, "ibeam")) {
         chkfread(&(fb->ibeam),sizeof(int),1,inputfile);
         totalbytes += sizeof(int);
      } else if (expecting_rawdatafile) {
         strcpy(fb->inpfile, string);
         expecting_rawdatafile = 0;
      } else if (expecting_source_name) {
         strcpy(fb->source_name, string);
         expecting_source_name = 0;
      } else {
         sprintf(message, "read_filterbank_header - unknown parameter: %s\n",
                 string);
         fprintf(stderr, "ERROR: %s\n", message);
         exit(1);
      }
   }
   /* add on last header string */
   totalbytes += nbytes;
   /* return total number of bytes read */
   fb->headerlen = totalbytes;
   /* Calculate the number of samples in the file */
   fb->N = (chkfilelen(inputfile, 1) - fb->headerlen) / fb->nchans * 8 / fb->nbits;
   return totalbytes;
}


void print_filterbank_header(sigprocfb * fb)
/* Output a SIGPROC filterbank header in human readable form */
{
   char *ctmp;

   printf("\n        Header size (bytes) = %d\n", fb->headerlen);
   ctmp = telescope_name(fb->telescope_id);
   printf("                  Telescope = %s (%d)\n", ctmp, fb->telescope_id);
   free(ctmp);
   ctmp = backend_name(fb->machine_id);
   printf("                 Instrument = %s (%d)\n", ctmp, fb->machine_id);
   free(ctmp);
   printf("                Source Name = %s\n", fb->source_name);
   printf("   Original input file name = '%s'\n", fb->inpfile);
   printf("             MJD start time = %.12f\n", fb->tstart);
   printf("    RA (J2000, HHMMSS.SSSS) = %.4f\n", fb->src_raj);
   printf("   DEC (J2000, DDMMSS.SSSS) = %.4f\n", fb->src_dej);
   printf("          Number of samples = %lld\n", fb->N);
   printf("        File duration (sec) = %-17.15g\n", fb->N * fb->tsamp);
   printf("                T_samp (us) = %-17.6g\n", fb->tsamp * 1e6);
   printf("    High channel freq (MHz) = %-17.15g\n", fb->fch1);
   printf("         Central freq (MHz) = %-17.15g\n",
          fb->fch1 + (fb->nchans / 2 - 0.5) * fb->foff);
   printf("         Number of channels = %d\n", fb->nchans);
   printf("    Channel Bandwidth (MHz) = %-17.15g\n", fb->foff);
   printf("      Total Bandwidth (MHz) = %-17.15g\n", fabs(fb->foff * fb->nchans));
   printf("      Number of IFs present = %d\n", fb->nifs);
   printf("            Bits per sample = %d\n", fb->nbits);
   if (fb->az_start != 0.0 || fb->za_start != 0.0){
      printf("        Start azimuth (deg) = %.3f\n", fb->az_start);
      printf("   Start zenith angle (deg) = %.3f\n", fb->za_start);
   }
   if (fb->sumifs)
      printf(" Note:  IFs were summed in hardware.\n");
}

void get_filterbank_static(int *ptsperbyte, int *bytesperpt, 
                           int *bytesperblk, float *clip_sigma)
{
   *ptsperbyte = ptsperbyte_st;
   *bytesperpt = bytesperpt_st;
   *bytesperblk = bytesperblk_st;
   *clip_sigma = clip_sigma_st;
}


void set_filterbank_static(int ptsperbyte, int ptsperblk, 
                           int bytesperpt, int bytesperblk,
                           int numchan, float clip_sigma, double dt)
{
   using_MPI = 1;
   currentblock = 0;
   ptsperblk_st = ptsperblk;
   bytesperpt_st = bytesperpt;
   bytesperblk_st = bytesperblk;
   ptsperbyte_st = ptsperbyte;
   numchan_st = numchan;
   sampperblk_st = ptsperblk_st * numchan_st;
   clip_sigma_st = clip_sigma;
   dt_st = dt;
}


void set_filterbank_padvals(float *fpadvals, int good_padvals)
{
   int ii;
   float sum_padvals = 0.0;

   if (good_padvals) {
      for (ii = 0; ii < numchan_st; ii++) {
         padvals[ii] = (unsigned char) (fpadvals[ii] + 0.5);
         sum_padvals += fpadvals[ii];
      }
      padval = (unsigned char) (sum_padvals / numchan_st + 0.5);
   } else {
      for (ii = 0; ii < numchan_st; ii++)
         padvals[ii] = padval;
   }
}


void sigprocfb_to_inf(sigprocfb * fb, infodata * idata)
/* Convert filterbank header structure into an infodata structure */
{
   char *tmpstr;
   strncpy(idata->object, fb->source_name, 80);
   idata->ra_h = (int) floor(fb->src_raj / 10000.0);
   idata->ra_m = (int) floor((fb->src_raj - idata->ra_h * 10000) / 100.0);
   idata->ra_s = fb->src_raj - idata->ra_h * 10000 - idata->ra_m * 100;
   idata->dec_d = (int) floor(fabs(fb->src_dej) / 10000.0);
   idata->dec_m = (int) floor((fabs(fb->src_dej) - idata->dec_d * 10000) / 100.0);
   idata->dec_s = fabs(fb->src_dej) - idata->dec_d * 10000 - idata->dec_m * 100;
   if (fb->src_dej < 0.0)
      idata->dec_d = -idata->dec_d;
   tmpstr = telescope_name(fb->telescope_id);
   strcpy(idata->telescope, tmpstr);
   free(tmpstr);
   tmpstr = backend_name(fb->machine_id);
   strcpy(idata->instrument, tmpstr);
   free(tmpstr);
   idata->mjd_i = (int) floor(fb->tstart);
   idata->mjd_f = fb->tstart - idata->mjd_i;
   idata->dt = fb->tsamp;
   idata->N = (double) fb->N;
   idata->chan_wid = fabs(fb->foff);
   idata->num_chan = fb->nchans;
   idata->freqband = idata->chan_wid * idata->num_chan;
   idata->freq = fb->fch1 - (idata->num_chan - 1) * idata->chan_wid;
   idata->fov = 1.2 * SOL * 3600.0 / (idata->freq * 1.0e6) / Tdiam * RADTODEG;
   idata->bary = 0;
   idata->numonoff = 0;
   idata->dm = 0.0;
   strcpy(idata->band, "Radio");
   strcpy(idata->analyzer, "Scott Ransom");
   strcpy(idata->observer, "Unknown");
   sprintf(idata->notes, "Input filterbank samples have %d bits.\n", fb->nbits);
}


void get_filterbank_file_info(FILE * files[], int numfiles, float clipsig,
                              long long *N, int *ptsperblock, int *numchan,
                              double *dt, double *T, int output)
/* Read basic information into static variables and make padding        */
/* calculations for a set of filterbank rawfiles that you want to patch */
/* together.  N, numchan, dt, and T are return values and include all   */
/* the files with the required padding.  If output is true, prints      */
/* a table showing a summary of the values.                             */
{
   int ii, headerlen;

   /* Allocate memory for our information structures */
   fb_st = (sigprocfb *) malloc(sizeof(sigprocfb) * numfiles);
   idata_st = (infodata *) malloc(sizeof(infodata) * numfiles);
   numpts_st = (long long *) malloc(sizeof(long long) * numfiles);
   padpts_st = (long long *) malloc(sizeof(long long) * numfiles);
   numblks_st = (int *) malloc(sizeof(int) * numfiles);
   times_st = (double *) malloc(sizeof(double) * numfiles);
   mjds_st = (double *) malloc(sizeof(double) * numfiles);
   elapsed_st = (double *) malloc(sizeof(double) * numfiles);
   startblk_st = (double *) malloc(sizeof(double) * numfiles);
   endblk_st = (double *) malloc(sizeof(double) * numfiles);

   /* Now read the first header... */
   headerlen = read_filterbank_header(&(fb_st[0]), files[0]);
   if (fb_st[0].nbits == 8) {
       ptsperbyte_st = 1;
   } else if (fb_st[0].nbits == 4) {
       ptsperbyte_st = 2;
   } else if (fb_st[0].nbits == 2) {
       ptsperbyte_st = 4;
   } else if (fb_st[0].nbits == 1) {
       ptsperbyte_st = 8;
   }
   sigprocfb_to_inf(fb_st + 0, idata_st + 0);
   chkfseek(files[0], fb_st[0].headerlen, SEEK_SET);
   if (clipsig > 0.0)
      clip_sigma_st = clipsig;
   *numchan = numchan_st = idata_st[0].num_chan;
   *ptsperblock = ptsperblk_st = BLOCKLEN;
   sampperblk_st = ptsperblk_st * numchan_st;
   bytesperblk_st = sampperblk_st / ptsperbyte_st;
   numblks_st[0] = chkfilelen(files[0], bytesperblk_st);
   numpts_st[0] = numblks_st[0] * ptsperblk_st;
   N_st = numpts_st[0];
   dt_st = *dt = idata_st[0].dt;
   times_st[0] = numpts_st[0] * dt_st;
   mjds_st[0] = idata_st[0].mjd_i + idata_st[0].mjd_f;
   elapsed_st[0] = 0.0;
   startblk_st[0] = 1;
   endblk_st[0] = numblks_st[0] - 1;
   padpts_st[0] = padpts_st[numfiles - 1] = 0;
   /* And read the rest of the headers */
   for (ii = 1; ii < numfiles; ii++) {
      headerlen = read_filterbank_header(&(fb_st[ii]), files[ii]);
      sigprocfb_to_inf(fb_st + ii, idata_st + ii);
      chkfseek(files[ii], fb_st[ii].headerlen, SEEK_SET);
      if (idata_st[ii].num_chan != numchan_st) {
         printf("Number of channels (file %d) is not the same!\n\n", ii + 1);
      }
      if (idata_st[ii].dt != dt_st) {
         printf("Sample time (file %d) is not the same!\n\n", ii + 1);
      }
      numblks_st[ii] = chkfilelen(files[ii], bytesperblk_st);
      numpts_st[ii] = numblks_st[ii] * ptsperblk_st;
      times_st[ii] = numpts_st[ii] * dt_st;
      mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
      elapsed_st[ii] = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
                                    idata_st[ii - 1].mjd_i, idata_st[ii - 1].mjd_f);
      padpts_st[ii - 1] = (long long) ((elapsed_st[ii] - times_st[ii - 1]) /
                                       dt_st + 0.5);
      elapsed_st[ii] += elapsed_st[ii - 1];
      N_st += numpts_st[ii] + padpts_st[ii - 1];
      startblk_st[ii] = (double) (N_st - numpts_st[ii]) / ptsperblk_st + 1;
      endblk_st[ii] = (double) (N_st) / ptsperblk_st - 1;
   }
   padpts_st[numfiles - 1] = ((long long) ceil(endblk_st[numfiles - 1] + 1.0) *
                              ptsperblk_st - N_st);
   N_st += padpts_st[numfiles - 1];
   *N = N_st;
   *T = T_st = N_st * dt_st;
   currentfile = currentblock = 0;
   if (output) {
      char *ctmp;
      ctmp = telescope_name(fb_st[0].telescope_id);
      printf("         Telescope = %s (%d)\n", ctmp, fb_st[0].telescope_id);
      free(ctmp);
      ctmp = backend_name(fb_st[0].machine_id);
      printf("        Instrument = %s (%d)\n", ctmp, fb_st[0].machine_id);
      free(ctmp);
      printf("   Number of files = %d\n", numfiles);
      printf("    Bits per point = %d\n", fb_st[0].nbits);
      printf("      Points/block = %d\n", ptsperblk_st);
      printf("   Num of channels = %d\n", numchan_st);
      printf("  Total points (N) = %lld\n", N_st);
      printf("  Sample time (dt) = %-14.14g\n", dt_st);
      printf("    Total time (s) = %-14.14g\n", T_st);
      printf(" Header length (B) = %d\n\n", headerlen);
      printf
          ("File  Start Block    Last Block     Points      Elapsed (s)      Time (s)            MJD           Padding\n");
      printf
          ("----  ------------  ------------  ----------  --------------  --------------  ------------------  ----------\n");
      for (ii = 0; ii < numfiles; ii++)
         printf
             ("%2d    %12.11g  %12.11g  %10lld  %14.13g  %14.13g  %17.12f  %10lld\n",
              ii + 1, startblk_st[ii], endblk_st[ii], numpts_st[ii], elapsed_st[ii],
              times_st[ii], mjds_st[ii], padpts_st[ii]);
      printf("\n");
   }
}

void filterbank_update_infodata(int numfiles, infodata * idata)
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


int skip_to_filterbank_rec(FILE * infiles[], int numfiles, int rec)
/* This routine skips to the record 'rec' in the input files     */
/* *infiles.  *infiles contains SIGPROC format data.             */
/* Returns the record skipped to.                                */
{
   double floor_blk;
   int filenum = 0, ii;

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
      for (ii = 0; ii < MAXNUMCHAN * BLOCKLEN; ii++)
         databuffer[ii] = padval;

      /* Warning:  I'm not sure if the following is correct. */
      /*   If really needs accurate testing to see if my     */
      /*   offsets are correct.  Bottom line, don't trust    */
      /*   a TOA determined using the following!             */

      if (rec < startblk_st[filenum]) { /* Padding region */
         currentfile = filenum - 1;
         chkfileseek(infiles[currentfile], 0, 1, SEEK_END);
         bufferpts = padpts_st[currentfile] % ptsperblk_st;
         padnum = ptsperblk_st * (rec - endblk_st[currentfile] - 1);
         /*
            printf("Padding:  currentfile = %d  bufferpts = %d  padnum = %d\n", 
            currentfile, bufferpts, padnum);
          */
      } else {                  /* Data region */
         currentfile = filenum;
         chkfileseek(infiles[currentfile], rec - startblk_st[filenum],
                     bytesperblk_st, SEEK_CUR);
         bufferpts = (int) ((startblk_st[filenum] - floor_blk) * ptsperblk_st + 0.5);
         padnum = 0;
         /*
            printf("Data:  currentfile = %d  bufferpts = %d  padnum = %d\n", 
            currentfile, bufferpts, padnum);
          */
      }

   } else {
      printf("\n rec = %d out of range in skip_to_filterbank_rec()\n", rec);
      exit(1);
   }
   return rec;
}


int read_filterbank_rawblock(FILE * infiles[], int numfiles,
                             unsigned char *data, int *padding)
/* This routine reads a single record from the   */
/* input files *infiles which contain 1 byte     */
/* data in SIGPROC filterbank format.            */
/* If padding is returned as 1, then padding was */
/* added and statistics should not be calculated */
{
   int offset, numtopad = 0;
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

   if (chkfread((unsigned char *) rawdatabuffer, bytesperblk_st, 1, infiles[currentfile])) {    /* Got Data */
      /* See if we need to byte-swap and if so, doit */
      if (need_byteswap_st) {
         /* Need to add this later */
      }
      convert_filterbank_block(rawdatabuffer, dataptr);

      *padding = 0;
      /* Put the new data into the databuffer if needed */
      if (bufferpts)
         memcpy(data, databuffer, sampperblk_st);
      currentblock++;
      return 1;
   } else {                     /* Didn't get data */
      if (feof(infiles[currentfile])) { /* End of file? */
         numtopad = padpts_st[currentfile] - padnum;
         if (numtopad) {        /* Pad the data? */
            *padding = 1;
            if (numtopad >= ptsperblk_st - bufferpts) { /* Lots of padding */
               if (bufferpts) { /* Buffer the padding? */
                  /* Add the amount of padding we need to */
                  /* make our buffer offset = 0           */
                  numtopad = ptsperblk_st - bufferpts;
                  memset(dataptr, padval, numtopad * numchan_st);
                  /* Copy the new data/padding into the output array */
                  memcpy(data, databuffer, sampperblk_st);
                  bufferpts = 0;
               } else {         /* Add a full record of padding */
                  numtopad = ptsperblk_st;
                  memset(data, padval, sampperblk_st);
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
               memset(databuffer + bufferpts * numchan_st,
                      padval, numtopad * numchan_st);
               bufferpts += numtopad;
               padnum = 0;
               shiftbuffer = 0;
               currentfile++;
               return read_filterbank_rawblock(infiles, numfiles, data, &pad);
            }
         } else {               /* No padding needed.  Try reading the next file */
            currentfile++;
            shiftbuffer = 0;
            return read_filterbank_rawblock(infiles, numfiles, data, padding);
         }
      } else {
         printf("\nProblem reading record from filterbank data file:\n");
         printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                currentfile, currentblock);
         exit(1);
      }
   }
}


int read_filterbank_rawblocks(FILE * infiles[], int numfiles,
                              unsigned char rawdata[], int numblocks, int *padding)
/* This routine reads numblocks filterbank records from the input */
/* files *infiles.  The 8-bit filterbank data is returned   */
/* in rawdata which must have a size of numblocks*          */
/* sampperblk_st.  The number  of blocks read is returned.  */
/* If padding is returned as 1, then padding was added      */
/* and statistics should not be calculated                  */
{
   int ii, retval = 0, pad, numpad = 0;

   *padding = 0;
   for (ii = 0; ii < numblocks; ii++) {
      retval += read_filterbank_rawblock(infiles, numfiles,
                                         rawdata + ii * sampperblk_st, &pad);
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


int read_filterbank(FILE * infiles[], int numfiles, float *data,
                    int numpts, double *dispdelays, int *padding,
                    int *maskchans, int *nummasked, mask * obsmask)
/* This routine reads numpts from the filterbank raw   */
/* files *infiles.  These files contain 16 bit data    */
/* from the filterbank backend.  Time delays and       */
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
   static unsigned char *tempzz, *rawdata1, *rawdata2;
   static unsigned char *currentdata, *lastdata;
   static int firsttime = 1, numblocks = 1, allocd = 0, mask = 0;
   static double duration = 0.0, timeperblk = 0.0;

   *nummasked = 0;
   if (firsttime) {
      if (numpts % ptsperblk_st) {
         printf("numpts must be a multiple of %d in read_filterbank()!\n",
                ptsperblk_st);
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

   /* Read, convert and de-disperse */
   if (allocd) {
      while (1) {
         numread = read_filterbank_rawblocks(infiles, numfiles, currentdata,
                                             numblocks, padding);
         if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
         }

         /* Clip nasty RFI if requested and we're not masking all the channels */
         if ((clip_sigma_st > 0.0) && !(mask && (*nummasked == -1)))
            clip_times(currentdata, numpts, numchan_st, clip_sigma_st, padvals);

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


void get_filterbank_channel(int channum, float chandat[],
                            unsigned char rawdata[], int numblocks)
/* Return the values for channel 'channum' of a block of       */
/* 'numblocks' raw filterbank data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                */
/* read_filterbank_rawblocks(), and 'chandat' must have at least     */
/* 'numblocks' * 'ptsperblk_st' spaces.                        */
/* Channel 0 is assumed to be the lowest freq channel.         */
{
   int ii, jj;

   if (channum > numchan_st || channum < 0) {
      printf("\nchannum = %d is out of range in get_GMR_channel()!\n\n", channum);
      exit(1);
   }
   /* Select the correct channel */
   for (ii = 0, jj = channum; ii < numblocks * ptsperblk_st; ii++, jj += numchan_st)
      chandat[ii] = rawdata[jj];
}


int prep_filterbank_subbands(unsigned char *rawdata, float *data,
                             double *dispdelays, int numsubbands,
                             int transpose, int *maskchans,
                             int *nummasked, mask * obsmask)
/* This routine preps a block from from the filterbank system.  It uses         */
/* dispersion delays in 'dispdelays' to de-disperse the data into         */
/* 'numsubbands' subbands.  It stores the resulting data in vector 'data' */
/* of length 'numsubbands' * 'ptsperblk_st'.  The low freq subband is     */
/* stored first, then the next highest subband etc, with 'ptsperblk_st'   */
/* floating points per subband. It returns the # of points read if        */
/* succesful, 0 otherwise.  'maskchans' is an array of length numchans    */
/* which contains a list of the number of channels that were masked.  The */
/* # of channels masked is returned in 'nummasked'.  'obsmask' is the     */
/* mask structure to use for masking.  If 'transpose'==0, the data will   */
/* be kept in time order instead of arranged by subband as above.         */
{
   int ii, jj, trtn, offset;
   double starttime = 0.0;
   static unsigned char *tempzz;
   static unsigned char rawdata1[MAXNUMCHAN * BLOCKLEN],
       rawdata2[MAXNUMCHAN * BLOCKLEN];
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
      timeperblk = ptsperblk_st * dt_st;
   }

   /* Read and de-disperse */
   memcpy(currentdata, rawdata, sampperblk_st);
   if (mask) {
      starttime = currentblock * timeperblk;
      *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
   }

   /* Clip nasty RFI if requested and we're not masking all the channels*/
   if ((clip_sigma_st > 0.0) && !(mask && (*nummasked == -1)))
      clip_times(currentdata, ptsperblk_st, numchan_st, clip_sigma_st, padvals);

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


int read_filterbank_subbands(FILE * infiles[], int numfiles, float *data,
                             double *dispdelays, int numsubbands,
                             int transpose, int *padding,
                             int *maskchans, int *nummasked, mask * obsmask)
/* This routine reads a record from the input files *infiles[]   */
/* in SIGPROC filterbank format.  The routine uses               */
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
   static int firsttime = 1;
   static unsigned char rawdata[MAXNUMCHAN * BLOCKLEN];

   if (firsttime) {
      if (!read_filterbank_rawblock(infiles, numfiles, rawdata, padding)) {
         printf("Problem reading the raw filterbank data file.\n\n");
         return 0;
      }
      if (0 != prep_filterbank_subbands(rawdata, data, dispdelays, numsubbands,
                                        transpose, maskchans, nummasked, obsmask)) {
         printf("Problem initializing prep_filterbank_subbands()\n\n");
         return 0;
      }
      firsttime = 0;
   }
   if (!read_filterbank_rawblock(infiles, numfiles, rawdata, padding)) {
      /* printf("Problem reading the raw filterbank data file.\n\n"); */
      return 0;
   }
   return prep_filterbank_subbands(rawdata, data, dispdelays, numsubbands,
                                   transpose, maskchans, nummasked, obsmask);
}


void convert_filterbank_block(int *indata, unsigned char *outdata)
/* This routine converts SIGPROC filterbank-format data into PRESTO format */
{
   int ii, jj, samp_ct, offset;

   if (ptsperbyte_st == 1) {
       unsigned char *chardata = (unsigned char *) indata;
       for (samp_ct = 0; samp_ct < ptsperblk_st; samp_ct++) {
           offset = samp_ct * numchan_st;
           for (ii = 0, jj = numchan_st - 1; ii < numchan_st; ii++, jj--)
               outdata[ii + offset] = chardata[jj + offset];
       }
   } else if (ptsperbyte_st == 2) {
       unsigned char c, *chardata = (unsigned char *) indata;
       for (samp_ct = 0; samp_ct < ptsperblk_st; samp_ct++) {
           offset = samp_ct * numchan_st;
           for (ii = 0, jj = numchan_st/2 - 1; ii < numchan_st; ii+=2, jj--) {
               c = chardata[(jj + offset/2)];
               outdata[(ii+1) + offset] = (unsigned char ) (c & 15);
               outdata[ii + offset] = (unsigned char ) ( (c & 240) >> 4 );
           }
       }
   } else if (ptsperbyte_st == 4) {
       unsigned char c, *chardata = (unsigned char *) indata;
       unsigned char *outdataptr;
       for (samp_ct = 0; samp_ct < ptsperblk_st; samp_ct++) {
           offset = samp_ct * numchan_st;
           outdataptr = outdata + offset;
           for (ii = 0, jj = numchan_st/4 - 1; ii < numchan_st; ii+=4, jj--) {
               c = chardata[(jj + offset/4)];
               *outdataptr++ = (c >> 0x06) & 0x03;
               *outdataptr++ = (c >> 0x04) & 0x03;
               *outdataptr++ = (c >> 0x02) & 0x03;
               *outdataptr++ = c & 0x03;
           }
       }
   } else if (ptsperbyte_st == 8) {
       unsigned char c, *chardata = (unsigned char *) indata;
       unsigned char *outdataptr;
       for (samp_ct = 0; samp_ct < ptsperblk_st; samp_ct++) {
           offset = samp_ct * numchan_st;
           outdataptr = outdata + offset;
           for (ii = 0, jj = numchan_st/8 - 1; ii < numchan_st; ii+=8, jj--) {
               c = chardata[(jj + offset/8)];
               *outdataptr++ = (c >> 0x07) & 0x01;
               *outdataptr++ = (c >> 0x06) & 0x01;
               *outdataptr++ = (c >> 0x05) & 0x01;
               *outdataptr++ = (c >> 0x04) & 0x01;
               *outdataptr++ = (c >> 0x03) & 0x01;
               *outdataptr++ = (c >> 0x02) & 0x01;
               *outdataptr++ = (c >> 0x01) & 0x01;
               *outdataptr++ = c & 0x01;
           }
       }
   } else {
       printf("\nYikes!!! Not supposed to be here in convert_filterbank_block()\n\n");
   }
}
