#include <sys/types.h>
#include <pwd.h>
#include "presto.h"
#include "mask.h"
#include "psrfits.h"

#define DEBUGOUT 1

static struct spectra_info S;
static unsigned char *rawbuffer, *ringbuffer, *tmpbuffer;
static float *offsets, *scales, global_scale = 1.0;
static unsigned char *padvals=NULL, *newpadvals=NULL;
static int cur_file = 0, cur_subint = 1, cur_specoffs = 0, padval = 0;
static int bufferspec = 0, padnum = 0, shiftbuffer = 1, missing_blocks = 0;
static int using_MPI = 0, default_poln = 0, user_poln = 0;
static double last_offs_sub = 0.0;

extern double slaCldj(int iy, int im, int id, int *j);
extern short transpose_bytes(unsigned char *a, int nx, int ny, unsigned char *move,
                             int move_size);

// This tests to see if 2 times are within 100ns of each other
#define TEST_CLOSE(a, b) (fabs((a)-(b)) <= 1e-7 ? 1 : 0)

void get_PSRFITS_static(int *bytesperpt, int *bytesperblk, int *numifs,
                       float *clip_sigma)
// This is used by the MPI-based mpiprepsubband code to insure that each
// processor can access (correctly) the required static variables.
{
   *bytesperpt = S.bytes_per_spectra;
   *bytesperblk = S.bytes_per_subint;
   *numifs = S.num_polns;
   *clip_sigma = S.clip_sigma;
}

void set_PSRFITS_static(int ptsperblk, int bytesperpt, int bytesperblk,
                       int numchan, int numifs, float clip_sigma, double dt)
// This is used by the MPI-based mpiprepsubband code to insure that each
// processor can access (correctly) the required static variables.
{
   using_MPI = 1;
   cur_subint = 1;
   S.spectra_per_subint = ptsperblk;
   S.bytes_per_spectra = bytesperpt;
   S.bytes_per_subint = bytesperblk;
   S.bits_per_sample = (bytesperpt * 8) / numchan / numifs;
   S.num_channels = numchan;
   S.num_polns = numifs;
   S.summed_polns = (numifs==1) ? 1 : 0;
   S.clip_sigma = clip_sigma;
   S.dt = dt;
   S.time_per_subint = S.dt * S.spectra_per_subint;
   // Note: the following only makes this correct for starting at 
   // the beginning of the data
   if (last_offs_sub==0.0)
       last_offs_sub = - S.time_per_subint;
}

void set_PSRFITS_padvals(float *fpadvals, int good_padvals)
// If reasonable padding values are available from rifind,
// use these as the padding values, otherwise, set them all
// to the value of padval defined in the static variables above.
{
   int ii;
   float sum_padvals = 0.0;

   if (padvals==NULL) { // This is for the clients in mpiprepsubband
       padvals = gen_bvect(S.num_channels);
       newpadvals = padvals;
   }

   if (good_padvals) {
      for (ii = 0; ii < S.num_channels; ii++) {
         padvals[ii] = newpadvals[ii] = (unsigned char) (fpadvals[ii] + 0.5);
         sum_padvals += fpadvals[ii];
      }
      padval = (unsigned char) (sum_padvals / S.num_channels + 0.5);
   } else {
      for (ii = 0; ii < S.num_channels; ii++)
         padvals[ii] = newpadvals[ii] = padval;
   }
}

static double DATEOBS_to_MJD(char *dateobs, int *mjd_day, double *mjd_fracday)
// Convert DATE-OBS string from PSRFITS primary HDU to a MJD
{
   int year, month, day, hour, min, err;
   double sec;

   sscanf(dateobs, "%4d-%2d-%2dT%2d:%2d:%lf", 
          &year, &month, &day, &hour, &min, &sec);
   *mjd_fracday = (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
   *mjd_day = slaCldj(year, month, day, &err);
   return *mjd_day + *mjd_fracday;
}

int is_PSRFITS(char *filename)
// Return 1 if the file described by filename is a PSRFITS file
// Return 0 otherwise.                                         
{
    fitsfile *fptr;
    int status=0;
    char ctmp[80], comment[120];
    
    // Read the primary HDU
    fits_open_file(&fptr, filename, READONLY, &status);
    if (status) return 0;
    
    // Make the easy check first
    fits_read_key(fptr, TSTRING, "FITSTYPE", ctmp, comment, &status);
    if (status || strcmp(ctmp, "PSRFITS")) return 0;

    // See if the data are search-mode
    fits_read_key(fptr, TSTRING, "OBS_MODE", ctmp, comment, &status);
    if (status || (strcmp(ctmp, "SEARCH") && 
                   strcmp(ctmp, "SRCH"))) return 0;

    fits_close_file(fptr, &status);
    return 1;  // it is search-mode  PSRFITS
}

#define get_hdr_string(name, param) {                                   \
        fits_read_key(s->files[ii], TSTRING, (name), ctmp, comment, &status); \
        if (status) {\
            printf("Error %d reading key %s\n", status, name); \
            if (ii==0) param[0]='\0'; \
            if (status==KEY_NO_EXIST) status=0;                      \
        } else {                                                     \
            if (ii==0) strncpy((param), ctmp, 40);                          \
            else if (strcmp((param), ctmp)!=0)                              \
                printf("Warning!:  %s values don't match for files 0 and %d!\n", \
                       (name), ii);                                         \
        }                                                               \
    }

#define get_hdr_int(name, param) {                                      \
        fits_read_key(s->files[ii], TINT, (name), &itmp, comment, &status); \
        if (status) {\
            printf("Error %d reading key %s\n", status, name); \
            if (ii==0) param=0; \
            if (status==KEY_NO_EXIST) status=0;\
        } else {                                                          \
            if (ii==0) param = itmp;                                        \
            else if (param != itmp)                                         \
                printf("Warning!:  %s values don't match for files 0 and %d!\n", \
                       (name), ii);                                         \
        }                                                               \
    }

#define get_hdr_double(name, param) {                                   \
        fits_read_key(s->files[ii], TDOUBLE, (name), &dtmp, comment, &status); \
        if (status) {\
            printf("Error %d reading key %s\n", status, name); \
            if (ii==0.0) param=0.0; \
            if (status==KEY_NO_EXIST) status=0;\
        } else {                                                          \
            if (ii==0) param = dtmp;                                        \
            else if (param != dtmp)                                         \
                printf("Warning!:  %s values don't match for files 0 and %d!\n", \
                       (name), ii);                                         \
        }                                                               \
    }

int read_PSRFITS_files(char **filenames, int numfiles, struct spectra_info *s)
// Read and convert PSRFITS information from a group of files 
// and place the resulting info into a spectra_info structure.
// Return 1 if successful, 0 if not.
{
    int IMJD, SMJD, itmp, ii, status = 0;
    double OFFS, dtmp;
    long double MJDf;
    char ctmp[80], comment[120];
    
    if (numfiles > MAXPFITSFILES) {
        printf("Error!: There are more than %d input files!\n", MAXPFITSFILES);
        exit(1);
    }
    s->num_files = numfiles;
    s->N = 0;

    // The apply_{scale,offset,weight} flags should be -1 if 
    // we are going to let the data decide if we need to apply
    // them.  If they come preset to 0 (False), then we will
    // not apply them even if the data suggests we should.
    // However, by default, apply_flipband should start at 0
    s->apply_flipband = 0;

    // Step through the other files
    for (ii = 0 ; ii < numfiles ; ii++) {

#if DEBUG_OUT       
        printf("Reading '%s'\n", filenames[ii]);
#endif
        // Is the file a PSRFITS file?
        if (!is_PSRFITS(filenames[ii])) {
            fprintf(stderr, 
                    "\nError!  File '%s' does not appear to be PSRFITS!\n", 
                    filenames[ii]);
            return 0;
        }
        
        // Open the PSRFITS file
        fits_open_file(&(s->files[ii]), filenames[ii], READONLY, &status);

        // Is the data in search mode?
        fits_read_key(s->files[ii], TSTRING, "OBS_MODE", ctmp, comment, &status);
        // Quick fix for Parkes DFB data (SRCH?  why????)...
        if (strcmp("SRCH", ctmp)==0) {
            strncpy(ctmp, "SEARCH", 40);
        }
        if (strcmp(ctmp, "SEARCH")) {
            fprintf(stderr, 
                    "\nError!  File '%s' does not contain SEARCH-mode data!\n", 
                    filenames[ii]);
            return 0;
        }

        // Now get the stuff we need from the primary HDU header
        fits_read_key(s->files[ii], TSTRING, "TELESCOP", ctmp, comment, &status); \
        // Quick fix for MockSpec data...
        if (strcmp("ARECIBO 305m", ctmp)==0) {
            strncpy(ctmp, "Arecibo", 40);
        }
        // Quick fix for Parkes DFB data...
        {
            char newctmp[80];

            // Copy ctmp first since strlower() is in-place
            strcpy(newctmp, ctmp);
            if (strcmp("parkes", strlower(remove_whitespace(newctmp)))==0) {
                strncpy(ctmp, "Parkes", 40);
            }
        }
        if (status) {
            printf("Error %d reading key %s\n", status, "TELESCOP");
            if (ii==0) s->telescope[0]='\0';
            if (status==KEY_NO_EXIST) status=0;
        } else {
            if (ii==0) strncpy(s->telescope, ctmp, 40);
            else if (strcmp(s->telescope, ctmp)!=0)
                printf("Warning!:  %s values don't match for files 0 and %d!\n",
                       "TELESCOP", ii);
        }

        get_hdr_string("OBSERVER", s->observer);
        get_hdr_string("SRC_NAME", s->source);
        get_hdr_string("FRONTEND", s->frontend);
        get_hdr_string("BACKEND", s->backend);
        get_hdr_string("PROJID", s->project_id);
        get_hdr_string("DATE-OBS", s->date_obs);
        get_hdr_string("FD_POLN", s->poln_type);
        get_hdr_string("RA", s->ra_str);
        get_hdr_string("DEC", s->dec_str);
        get_hdr_double("OBSFREQ", s->fctr);
        get_hdr_int("OBSNCHAN", s->orig_num_chan);
        get_hdr_double("OBSBW", s->orig_df);
        //get_hdr_double("CHAN_DM", s->chan_dm);
        get_hdr_double("BMIN", s->beam_FWHM);

        /* This is likely not in earlier versions of PSRFITS so */
        /* treat it a bit differently                           */
        fits_read_key(s->files[ii], TDOUBLE, "CHAN_DM", 
                      &(s->chan_dm), comment, &status);
        if (status==KEY_NO_EXIST) {
            status = 0;
            s->chan_dm = 0.0;
        }

        // Don't use the macros unless you are using the struct!
        fits_read_key(s->files[ii], TINT, "STT_IMJD", &IMJD, comment, &status);
        s->start_MJD[ii] = (long double) IMJD;
        fits_read_key(s->files[ii], TINT, "STT_SMJD", &SMJD, comment, &status);
        fits_read_key(s->files[ii], TDOUBLE, "STT_OFFS", &OFFS, comment, &status);
        s->start_MJD[ii] += ((long double) SMJD + (long double) OFFS) / SECPERDAY;

        // Are we tracking?
        fits_read_key(s->files[ii], TSTRING, "TRK_MODE", ctmp, comment, &status);
        itmp = (strcmp("TRACK", ctmp)==0) ? 1 : 0;
        if (ii==0) s->tracking = itmp;
        else if (s->tracking != itmp)
            printf("Warning!:  TRK_MODE values don't match for files 0 and %d!\n", ii);

        // Now switch to the SUBINT HDU header
        fits_movnam_hdu(s->files[ii], BINARY_TBL, "SUBINT", 0, &status);
        get_hdr_double("TBIN", s->dt);
        get_hdr_int("NCHAN", s->num_channels);
        get_hdr_int("NPOL", s->num_polns);
        {
            char *envval = getenv("PSRFITS_POLN");
            if (envval != NULL) {
                int ival = strtol(envval, NULL, 10);
                if ((ival > -1) && (ival < s->num_polns)) {
                    printf
                        ("Using polarization %d (from 0-%d) from PSRFITS_POLN.\n",
                         ival, s->num_polns-1);
                    default_poln = ival;
                    user_poln = 1;
                }
            }
        }
        get_hdr_string("POL_TYPE", s->poln_order);
        fits_read_key(s->files[ii], TINT, "NCHNOFFS", &itmp, comment, &status);
        if (itmp > 0)
            printf("Warning!:  First freq channel is not 0 in file %d!\n", ii);
        get_hdr_int("NSBLK", s->spectra_per_subint);
        get_hdr_int("NBITS", s->bits_per_sample);
        fits_read_key(s->files[ii], TINT, "NAXIS2", 
                      &(s->num_subint[ii]), comment, &status);
        fits_read_key(s->files[ii], TINT, "NSUBOFFS", 
                      &(s->start_subint[ii]), comment, &status);
        s->time_per_subint = s->dt * s->spectra_per_subint;

        /* This is likely not in earlier versions of PSRFITS so */
        /* treat it a bit differently                           */
        fits_read_key(s->files[ii], TFLOAT, "ZERO_OFF", 
                      &(s->zero_offset), comment, &status);
        if (status==KEY_NO_EXIST) {
            status = 0;
            s->zero_offset = 0.0;
        }
        s->zero_offset = fabs(s->zero_offset);

        // Get the time offset column info and the offset for the 1st row
        {
            double offs_sub;
            int colnum, anynull, numrows;

            // Identify the OFFS_SUB column number
            fits_get_colnum(s->files[ii], 0, "OFFS_SUB", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the OFFS_SUB column!\n");
                status = 0; // Reset status
            } else {
                if (ii==0) {
                    s->offs_sub_col = colnum;
                } else if (colnum != s->offs_sub_col) {
                    printf("Warning!:  OFFS_SUB column changes between files!\n");
                }
            }

            // Read the OFFS_SUB column value for the 1st row
            fits_read_col(s->files[ii], TDOUBLE,
                          s->offs_sub_col, 1L, 1L, 1L,
                          0, &offs_sub, &anynull, &status);

            numrows = (int)((offs_sub - 0.5 * s->time_per_subint) /
                            s->time_per_subint + 1e-7);
            // Check to see if any rows have been deleted or are missing
            if (numrows > s->start_subint[ii]) {
                printf("Warning: NSUBOFFS reports %d previous rows\n"
                       "         but OFFS_SUB implies %d.  Using OFFS_SUB.\n"
                       "         Will likely be able to correct for this.\n",
                       s->start_subint[ii], numrows);
            }
            s->start_subint[ii] = numrows;
        }

        // This is the MJD offset based on the starting subint number
        MJDf = (s->time_per_subint * s->start_subint[ii]) / SECPERDAY;
        // The start_MJD values should always be correct
        s->start_MJD[ii] += MJDf;

        // Compute the starting spectra from the times
        MJDf = s->start_MJD[ii] - s->start_MJD[0];
        if (MJDf < 0.0) {
            printf("Error!: File %d seems to be from before file 0!\n", ii); 
            exit(1);
        }
        s->start_spec[ii] = (long long)(MJDf * SECPERDAY / s->dt + 0.5);

        // Now pull stuff from the other columns
        {
            float ftmp;
            long repeat, width;
            int colnum, anynull;
            
            // Identify the data column and the data type
            fits_get_colnum(s->files[ii], 0, "DATA", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the DATA column!\n");
                status = 0; // Reset status
            } else {
                if (ii==0) {
                    s->data_col = colnum;
                    fits_get_coltype(s->files[ii], colnum, &(s->FITS_typecode), 
                                     &repeat, &width, &status);
                    // This makes CFITSIO treat 1-bit data as written in 'B' mode
                    // even if it was written in 'X' mode originally.  This means
                    // that we unpack it ourselves.
                    if (s->bits_per_sample==1 && s->FITS_typecode==1) {
                        s->FITS_typecode = 11;
                    }
                } else if (colnum != s->data_col) {
                    printf("Warning!:  DATA column changes between files!\n");
                }
            }
            
            // Telescope azimuth
            fits_get_colnum(s->files[ii], 0, "TEL_AZ", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                s->azimuth = 0.0;
                status = 0; // Reset status
            } else {
                fits_read_col(s->files[ii], TFLOAT, colnum, 
                              1L, 1L, 1L, 0, &ftmp, &anynull, &status);
                if (ii==0) s->azimuth = (double) ftmp;
            }
            
            // Telescope zenith angle
            fits_get_colnum(s->files[ii], 0, "TEL_ZEN", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                s->zenith_ang = 0.0;
                status = 0; // Reset status
            } else {
                fits_read_col(s->files[ii], TFLOAT, colnum, 
                              1L, 1L, 1L, 0, &ftmp, &anynull, &status);
                if (ii==0) s->zenith_ang = (double) ftmp;
            }
            
            // Observing frequencies
            fits_get_colnum(s->files[ii], 0, "DAT_FREQ", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel freq column!\n");
                status = 0; // Reset status
            } else {
                int jj;
                float *freqs = (float *)malloc(sizeof(float) * s->num_channels);
                fits_read_col(s->files[ii], TFLOAT, colnum, 1L, 1L, 
                              s->num_channels, 0, freqs, &anynull, &status);
                
                if (ii==0) {
		  s->df = ((double)freqs[s->num_channels-1]-(double)freqs[0])/(double)(s->num_channels-1);
                    s->lo_freq = freqs[0];
                    s->hi_freq = freqs[s->num_channels-1];
                    // Now check that the channel spacing is the same throughout
                    for (jj = 0 ; jj < s->num_channels - 1 ; jj++) {
                        ftmp = freqs[jj+1] - freqs[jj];
                        if (fabs(ftmp - s->df) > 1e-7)
                            printf("Warning!:  Channel spacing changes in file %d!\n", ii);
                    }
                } else {
                    ftmp = fabs(s->df-(freqs[1]-freqs[0]));
                    if (ftmp > 1e-7)
                        printf("Warning!:  Channel spacing changes between files!\n");
                    ftmp = fabs(s->lo_freq-freqs[0]);
                    if (ftmp > 1e-7)
                        printf("Warning!:  Low channel changes between files!\n");
                    ftmp = fabs(s->hi_freq-freqs[s->num_channels-1]);
                    if (ftmp > 1e-7)
                        printf("Warning!:  High channel changes between files!\n");
                }
                free(freqs);
            }
            
            // Data weights
            fits_get_colnum(s->files[ii], 0, "DAT_WTS", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel weights!\n");
                status = 0; // Reset status
            } else {
                if (s->apply_weight < 0) { // Use the data to decide
                    int jj;
                    if (ii==0) {
                        s->dat_wts_col = colnum;
                    } else if (colnum != s->dat_wts_col) {
                        printf("Warning!:  DAT_WTS column changes between files!\n");
                    }
                    float *fvec = (float *)malloc(sizeof(float) * s->num_channels);
                    fits_read_col(s->files[ii], TFLOAT, s->dat_wts_col, 1L, 1L, 
                                  s->num_channels, 0, fvec, &anynull, &status);
                    for (jj = 0 ; jj < s->num_channels ; jj++) {
                        // If the weights are not 1, apply them
                        if (fvec[jj] != 1.0) {
                            s->apply_weight = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_weight < 0) s->apply_weight = 0;  // not needed
            }
            
            // Data offsets
            fits_get_colnum(s->files[ii], 0, "DAT_OFFS", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel offsets!\n");
                status = 0; // Reset status
            } else {
                if (s->apply_offset < 0) { // Use the data to decide
                    int jj;
                    if (ii==0) {
                        s->dat_offs_col = colnum;
                    } else if (colnum != s->dat_offs_col) {
                        printf("Warning!:  DAT_OFFS column changes between files!\n");
                    }
                    float *fvec = (float *)malloc(sizeof(float) * 
                                                  s->num_channels * s->num_polns);
                    fits_read_col(s->files[ii], TFLOAT, s->dat_offs_col, 1L, 1L, 
                                  s->num_channels * s->num_polns, 
                                  0, fvec, &anynull, &status);
                    for (jj = 0 ; jj < s->num_channels * s->num_polns ; jj++) {
                        // If the offsets are not 0, apply them
                        if (fvec[jj] != 0.0) {
                            s->apply_offset = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_offset < 0) s->apply_offset = 0; // not needed
            }
            
            // Data scalings
            fits_get_colnum(s->files[ii], 0, "DAT_SCL", &colnum, &status);
            if (status==COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel scalings!\n");
                status = 0; // Reset status
            } else {
                if (s->apply_scale < 0) { // Use the data to decide
                    int jj;
                    if (ii==0) {
                        s->dat_scl_col = colnum;
                    } else if (colnum != s->dat_scl_col) {
                        printf("Warning!:  DAT_SCL column changes between files!\n");
                    }
                    float *fvec = (float *)malloc(sizeof(float) * 
                                                  s->num_channels * s->num_polns);
                    fits_read_col(s->files[ii], TFLOAT, colnum, 1L, 1L, 
                                  s->num_channels * s->num_polns, 
                                  0, fvec, &anynull, &status);
                    for (jj = 0 ; jj < s->num_channels * s->num_polns ; jj++) {
                        // If the scales are not 1, apply them
                        if (fvec[jj] != 1.0) {
                            s->apply_scale = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_scale < 0) s->apply_scale = 0; // not needed
            }
        }
        
        // Compute the samples per file and the amount of padding
        // that the _previous_ file has
        s->num_pad[ii] = 0;
        s->num_spec[ii] = s->spectra_per_subint * s->num_subint[ii];
        if (ii > 0) {
            if (s->start_spec[ii] > s->N) { // Need padding
                s->num_pad[ii-1] = s->start_spec[ii] - s->N;
                s->N += s->num_pad[ii-1];
            }
        }
        s->N += s->num_spec[ii];
    }

    // Convert the position strings into degrees
    {
        int d, h, m;
        double sec;
        ra_dec_from_string(s->ra_str, &h, &m, &sec);
        s->ra2000 = hms2rad(h, m, sec) * RADTODEG;
        ra_dec_from_string(s->dec_str, &d, &m, &sec);
        s->dec2000 = dms2rad(d, m, sec) * RADTODEG;
    }

    // Are the polarizations summed?
    if ((strncmp("AA+BB", s->poln_order, 5)==0) ||
        (strncmp("INTEN", s->poln_order, 5)==0))
        s->summed_polns = 1;
    else
        s->summed_polns = 0;

    // Calculate some others
    s->T = s->N * s->dt;
    s->orig_df /= (double) s->orig_num_chan;
    s->samples_per_spectra = s->num_polns * s->num_channels;
    // Note:  the following is the number of bytes that will be in
    //        the returned array from CFITSIO, possibly after processing by
    //        us given that we turn 1-, 2-, and 4-bit data into bytes
    //        immediately if bits_per_sample < 8
    if (s->bits_per_sample < 8)
        s->bytes_per_spectra = s->samples_per_spectra;
    else
        s->bytes_per_spectra = (s->bits_per_sample * s->samples_per_spectra) / 8;
    s->samples_per_subint = s->samples_per_spectra * s->spectra_per_subint;
    s->bytes_per_subint = s->bytes_per_spectra * s->spectra_per_subint;
    
    // Flip the band?
    if (s->hi_freq < s->lo_freq) {
        float ftmp = s->hi_freq;
        s->hi_freq = s->lo_freq;
        s->lo_freq = ftmp;
        s->df *= -1.0;
        s->apply_flipband = 1;
    }
    // Compute the bandwidth
    s->BW = s->num_channels * s->df;

    // Flip the bytes for Parkes FB_1BIT data
    if (s->bits_per_sample==1 &&
        strcmp(s->telescope, "Parkes")==0 &&
        strcmp(s->backend, "FB_1BIT")==0) {
        printf("Flipping bit ordering since Parkes FB_1BIT data.\n");
        s->flip_bytes = 1;
    } else {
        s->flip_bytes = 0;
    }

    // Copy the structures and return success
    S = *s;

    // Allocate the buffers
    rawbuffer = gen_bvect(s->bytes_per_subint/s->num_polns);
    ringbuffer = gen_bvect(2 * s->bytes_per_subint/s->num_polns);
    if (s->num_polns > 1)
        tmpbuffer = gen_bvect(s->samples_per_subint);
    // TODO:  The following is temporary, until I fix the padding
    padvals = gen_bvect(s->samples_per_spectra/s->num_polns);
    offsets = gen_fvect(s->samples_per_spectra);
    scales = gen_fvect(s->samples_per_spectra);
    for (ii = 0 ; ii < s->samples_per_spectra ; ii++) {
        offsets[ii] = 0.0;
        scales[ii] = 1.0;
    }
    for (ii = 0 ; ii < s->samples_per_spectra/s->num_polns ; ii++)
        padvals[ii] = 0;
    newpadvals = padvals;
    return 1;
}

void print_PSRFITS_info(struct spectra_info *s)
// Output a spectra_info structure in human readable form
{
    int ii, numhdus, hdutype, status = 0, itmp;
    char ctmp[40], comment[120];
    double dtmp;

    printf("From the PSRFITS file '%s':\n", s->files[0]->Fptr->filename);
    fits_get_num_hdus(s->files[0], &numhdus, &status);
    printf("                       HDUs = primary, ");
    for (ii = 2 ; ii < numhdus + 1 ; ii++) {
        fits_movabs_hdu(s->files[0], ii, &hdutype, &status);
        fits_read_key(s->files[0], TSTRING, "EXTNAME", ctmp, comment, &status);
        printf("%s%s", ctmp, (ii < numhdus) ? ", " : "\n");
    }
    printf("                  Telescope = %s\n", s->telescope);
    printf("                   Observer = %s\n", s->observer);
    printf("                Source Name = %s\n", s->source);
    printf("                   Frontend = %s\n", s->frontend);
    printf("                    Backend = %s\n", s->backend);
    printf("                 Project ID = %s\n", s->project_id);
    // printf("                Scan Number = %d\n", s->scan_number);
    printf("            Obs Date String = %s\n", s->date_obs);
    DATEOBS_to_MJD(s->date_obs, &itmp, &dtmp);
    sprintf(ctmp, "%.14f", dtmp);
    printf("  MJD start time (DATE-OBS) = %5i.%14s\n", itmp, ctmp+2);
    printf("     MJD start time (STT_*) = %19.14Lf\n", s->start_MJD[0]);
    printf("                   RA J2000 = %s\n", s->ra_str);
    printf("             RA J2000 (deg) = %-17.15g\n", s->ra2000);
    printf("                  Dec J2000 = %s\n", s->dec_str);
    printf("            Dec J2000 (deg) = %-17.15g\n", s->dec2000);
    printf("                  Tracking? = %s\n", s->tracking ? "True" : "False");
    printf("              Azimuth (deg) = %-.7g\n", s->azimuth);
    printf("           Zenith Ang (deg) = %-.7g\n", s->zenith_ang);
    printf("          Polarization type = %s\n", s->poln_type);
    if (s->num_polns>=2 && !s->summed_polns)
        printf("            Number of polns = %d\n", s->num_polns);
    else if (s->summed_polns)
        printf("            Number of polns = 2 (summed)\n");
    else 
        printf("            Number of polns = 1\n");
    printf("         Polarization order = %s\n", s->poln_order);
    printf("           Sample time (us) = %-17.15g\n", s->dt * 1e6);
    printf("         Central freq (MHz) = %-17.15g\n", s->fctr);
    printf("          Low channel (MHz) = %-17.15g\n", s->lo_freq);
    printf("         High channel (MHz) = %-17.15g\n", s->hi_freq);
    printf("        Channel width (MHz) = %-17.15g\n", s->df);
    printf("         Number of channels = %d\n", s->num_channels);
    if (s->chan_dm != 0.0) {
        printf("   Orig Channel width (MHz) = %-17.15g\n", s->orig_df);
        printf("    Orig Number of channels = %d\n", s->orig_num_chan);
        printf("    DM used for chan dedisp = %-17.15g\n", s->chan_dm);
    }
    printf("      Total Bandwidth (MHz) = %-17.15g\n", s->BW);
    printf("         Spectra per subint = %d\n", s->spectra_per_subint);
    printf("            Starting subint = %d\n", s->start_subint[0]);
    printf("           Subints per file = %d\n", s->num_subint[0]);
    printf("           Spectra per file = %lld\n", s->num_spec[0]);
    printf("      Time per subint (sec) = %-.12g\n", s->time_per_subint);
    printf("        Time per file (sec) = %-.12g\n", s->num_spec[0]*s->dt);
    printf("              FITS typecode = %d\n", s->FITS_typecode);
    printf("                DATA column = %d\n", s->data_col);
    printf("            bits per sample = %d\n", s->bits_per_sample);
    if (s->bits_per_sample < 8)
        itmp = (s->bytes_per_spectra * s->bits_per_sample) / 8;
    else
        itmp = s->bytes_per_spectra;
    printf("          bytes per spectra = %d\n", itmp);
    printf("        samples per spectra = %d\n", s->samples_per_spectra);
    if (s->bits_per_sample < 8)
        itmp = (s->bytes_per_subint * s->bits_per_sample) / 8;
    else
        itmp = s->bytes_per_subint;
    printf("           bytes per subint = %d\n", itmp);
    printf("         samples per subint = %d\n", s->samples_per_subint);
    printf("                zero offset = %-17.15g\n", s->zero_offset);
    printf("             Apply scaling? = %s\n", s->apply_scale ? "True" : "False");
    printf("             Apply offsets? = %s\n", s->apply_offset ? "True" : "False");
    printf("             Apply weights? = %s\n", s->apply_weight ? "True" : "False");
    printf("           Invert the band? = %s\n", s->apply_flipband ? "True" : "False");
}


void spectra_info_to_inf(struct spectra_info * s, infodata * idata)
// Convert a spectra_info structure into an infodata structure
{
    char ctmp[100];
    struct passwd *pwd;
    
    strncpy(idata->object, s->source, 80);
    hours2hms(s->ra2000 / 15.0, &(idata->ra_h), &(idata->ra_m), &(idata->ra_s));
    deg2dms(s->dec2000, &(idata->dec_d), &(idata->dec_m), &(idata->dec_s));
    strcpy(idata->telescope, s->telescope);
    strcpy(idata->instrument, s->backend);
    idata->num_chan = s->num_channels;
    idata->dt = s->dt;
    // DATEOBS_to_MJD(s->date_obs, &(idata->mjd_i), &(idata->mjd_f));
    idata->mjd_i = (int)(s->start_MJD[0]);
    idata->mjd_f = s->start_MJD[0] - idata->mjd_i;
    idata->N = s->N;
    idata->freqband = s->BW;
    idata->chan_wid = s->df;
    idata->freq = s->lo_freq;
    idata->fov = s->beam_FWHM * 3600.0; // in arcsec
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    pwd = getpwuid(geteuid());
    strcpy(idata->analyzer, pwd->pw_name);
    strncpy(idata->observer, s->observer, 80);
    if (s->summed_polns)
        sprintf(ctmp,
                "2 polns were summed.  Samples have %d bits.",
                s->bits_per_sample);
    else
        sprintf(ctmp, "%d polns were not summed.  Samples have %d bits.",
                s->num_polns, s->bits_per_sample);
    sprintf(idata->notes, "Project ID %s, Date: %s.\n    %s\n",
            s->project_id, s->date_obs, ctmp);
}


void get_PSRFITS_file_info(char **filenames, int numfiles, float clipsig, 
                           struct spectra_info * s, 
                           infodata * idata, int output)
// Read basic information into static variables and make padding
// calculations for a set of PSRFITS rawfiles that you want to patch
// together.  If output is true, prints a table showing a summary
// of the values.
{
    int ii, status;
    
    // Are we going to clip the data?
    if (clipsig > 0.0) s->clip_sigma = clipsig;
    else s->clip_sigma = 0.0;
    
    // Read the information from the input files
    status = read_PSRFITS_files(filenames, numfiles, s);

    // Make an infodata struct
    spectra_info_to_inf(s, idata);
    
    if (status && output) {
        printf("    Number of files = %d\n", numfiles);
        printf("     Spectra/subint = %d\n", s->spectra_per_subint);
        if (s->num_polns>=2 && !s->summed_polns)
            printf("       Num of polns = %d\n", s->num_polns);
        else if (s->summed_polns)
            printf("       Num of polns = 2 (summed)\n");
        else 
            printf("       Num of polns = 1\n");
        printf("    Num of channels = %d\n", s->num_channels);
        printf("  Center freq (MHz) = %.8g\n", s->fctr);
        printf("   Total points (N) = %lld\n", s->N);
        printf("   Sample time (dt) = %-14.14g\n", s->dt);
        printf("     Total time (s) = %-14.14g\n\n", s->T);
        printf("File   Samples      Padding        Start MJD\n");
        printf("----  ----------  ----------  --------------------\n");
        for (ii = 0; ii < numfiles; ii++)
            printf("%-4d  %10lld  %10lld  %19.14Lf\n", ii + 1,
                   s->num_spec[ii], s->num_pad[ii], s->start_MJD[ii]);
        printf("\n");
    }
}


void PSRFITS_update_infodata(infodata * idata)
// Update the onoff bins section in case we used multiple files
{
   int ii, index = 2;
   
   idata->N = S.N;
   if (S.num_files == 1 && S.num_pad[0] == 0) {
       idata->numonoff = 0;
       return;
   }
   // Determine the topocentric onoff bins
   idata->numonoff = 1;
   idata->onoff[0] = 0.0;
   idata->onoff[1] = S.num_spec[0] - 1.0;
   for (ii = 1; ii < S.num_files; ii++) {
       if (S.num_pad[ii - 1]) {
           idata->onoff[index] =
               idata->onoff[index - 1] + S.num_pad[ii - 1];
           idata->onoff[index + 1] = idata->onoff[index] + S.num_spec[ii];
           idata->numonoff++;
           index += 2;
       } else {
           idata->onoff[index - 1] += S.num_spec[ii];
       }
   }
   if (S.num_pad[S.num_files - 1]) {
       idata->onoff[index] =
           idata->onoff[index - 1] + S.num_pad[S.num_files - 1];
       idata->onoff[index + 1] = idata->onoff[index];
       idata->numonoff++;
   }
}


int skip_to_PSRFITS_samp(long long samplenum)
// This routine skips to the sample samplenum in the input
// PSRFITS files that are stored in S.  It returns 
// 1 if successful, 0 if unsuccessful.
{
    long long dsamp = 0;
    
    if (samplenum > 0)
        printf("\nWarning:  Offsetting by %lld > 0 samples.  This needs work!\n", 
               samplenum);
    else 
        return 0;

    // Find which file we should be in
    cur_file = 0;
    while (cur_file < S.num_files &&
           samplenum > S.start_spec[cur_file])
        cur_file++;
    if (cur_file == S.num_files) {
        printf("Error!:  Sample %lld does not exist in the input files!\n", 
               samplenum);
        exit(1);
    }
    
    shiftbuffer = 1;
    // If we are in the "data" part of the file
    if (samplenum < S.start_spec[cur_file] + S.num_spec[cur_file]) {
        dsamp = samplenum - S.start_spec[cur_file];
        // Remember that the subints are ones-offset
        cur_subint = dsamp / S.spectra_per_subint + 1;
        // Compute the number of spectra we need to offset in the
        // current subint to get to samplenum
        cur_specoffs = dsamp % S.spectra_per_subint;
        bufferspec = 0;
        padnum = 0;
#if DEBUGOUT        
        printf("Data:  currentfile = %d  bufferspec = %d  padnum = %d\n",
               cur_file, bufferspec, padnum);
#endif
    } else {  // We are in the padding part
        // Fill the buffer with padding (use zeros for now)
        memset(rawbuffer, 0, S.samples_per_subint/S.num_polns);
        cur_subint = 1; // placeholder
        bufferspec = 0;
        padnum = 0;
#if DEBUGOUT        
        printf("Padding:  currentfile = %d  bufferspec = %d  padnum = %d\n",
               cur_file, bufferspec, padnum);
#endif
    }
    return 1;
}


int read_PSRFITS_rawblock(unsigned char *data, int *padding)
// This routine reads a single subint (or equivalent) from the PSRFITs
// files in the static structure S.  *data must be at least
// S.bytes_per_subint bytes.  If padding is returned as 1, then
// padding was added and statistics should not be calculated.
{
    int offset = 0, numtopad = 0, ii, anynull, status = 0;
    unsigned char *dataptr = data;
    double offs_sub;
    static int firsttime = 1;
    static float *weights, *offsets;
    
    // The data will go directly into *data unless there is buffered
    // data in the ring buffer (this really only happens if patching
    // non-contiguous files together

    if (bufferspec) {
        offset = bufferspec * S.bytes_per_spectra/S.num_polns;
        dataptr = ringbuffer + offset;
        // If our buffer array is offset from last time, 
        // copy the previously offset part into the beginning.
        // New data comes after the old data in the buffer.
        if (shiftbuffer)
            memcpy(ringbuffer, ringbuffer + S.bytes_per_subint/S.num_polns, offset);
    }
    shiftbuffer = 1;
    
    // Make sure our current file number is valid
    if (cur_file >= S.num_files) return 0;

    // Read a subint of data from the DATA col
    if (cur_subint <= S.num_subint[cur_file]) {
        if (S.num_polns==1) tmpbuffer = dataptr;

        // Read the OFFS_SUB column value in case there were dropped blocks
        fits_read_col(S.files[cur_file], TDOUBLE, 
                      S.offs_sub_col, cur_subint, 1L, 1L, 
                      0, &offs_sub, &anynull, &status);

        if (firsttime) {
            if (S.apply_weight)
                weights = gen_fvect(S.num_channels);
            if (S.apply_offset)
                offsets = gen_fvect(S.num_channels*S.num_polns);
            if (S.apply_scale)
                scales = gen_fvect(S.num_channels*S.num_polns);
            last_offs_sub = offs_sub - S.time_per_subint;
        }

        // Read the weights, offsets, and scales if required
        if (S.apply_weight)
            fits_read_col(S.files[cur_file], TFLOAT, S.dat_wts_col, cur_subint, 1L, 
                          S.num_channels, 0, weights, &anynull, &status);
        if (S.apply_offset)
            fits_read_col(S.files[cur_file], TFLOAT, S.dat_offs_col, cur_subint, 1L, 
                          S.num_channels*S.num_polns, 0, offsets, &anynull, &status);
        if (S.apply_scale)
            fits_read_col(S.files[cur_file], TFLOAT, S.dat_scl_col, cur_subint, 1L, 
                          S.num_channels*S.num_polns, 0, scales, &anynull, &status);

        // The following determines if there were lost blocks
        if TEST_CLOSE(offs_sub-last_offs_sub, S.time_per_subint) {
            // if so, read the data from the column
            {
                int num_to_read = S.samples_per_subint;
                // The following allows us to read byte-packed data
                if (S.bits_per_sample < 8)
                    num_to_read = S.samples_per_subint * S.bits_per_sample / 8;
                fits_read_col(S.files[cur_file], S.FITS_typecode, 
                              S.data_col, cur_subint, 1L, num_to_read, 
                              0, tmpbuffer, &anynull, &status);
                // The following converts that byte-packed data into 
                // bytes, in place
                if (S.bits_per_sample == 4) {
                    int ii, jj;
                    unsigned char uctmp;
                    for (ii = num_to_read - 1, jj = 2 * num_to_read - 1 ; 
                         ii >= 0 ; ii--, jj -= 2) {
                        uctmp = (unsigned char)tmpbuffer[ii];
                        tmpbuffer[jj] = uctmp & 0x0F;
                        tmpbuffer[jj-1] = uctmp >> 4;
                    }
                } else if (S.bits_per_sample == 2) {
                    int ii, jj;
                    unsigned char uctmp;
                    for (ii = num_to_read - 1, jj = 4 * num_to_read - 1 ; 
                         ii >= 0 ; ii--, jj -= 4) {
                        uctmp = (unsigned char)tmpbuffer[ii];
                        tmpbuffer[jj] = (uctmp & 0x03);
                        tmpbuffer[jj-1] = ((uctmp >> 0x02) & 0x03);
                        tmpbuffer[jj-2] = ((uctmp >> 0x04) & 0x03);
                        tmpbuffer[jj-3] = ((uctmp >> 0x06) & 0x03);
                    }
                } else if (S.bits_per_sample == 1) {
                    int ii, jj;
                    unsigned char uctmp;
                    for (ii = num_to_read - 1, jj = 8 * num_to_read - 1 ; 
                         ii >= 0 ; ii--, jj -= 8) {
                        uctmp = (unsigned char)tmpbuffer[ii];
                        tmpbuffer[jj] = (uctmp & 0x01);
                        tmpbuffer[jj-1] = ((uctmp >> 0x01) & 0x01);
                        tmpbuffer[jj-2] = ((uctmp >> 0x02) & 0x01);
                        tmpbuffer[jj-3] = ((uctmp >> 0x03) & 0x01);
                        tmpbuffer[jj-4] = ((uctmp >> 0x04) & 0x01);
                        tmpbuffer[jj-5] = ((uctmp >> 0x05) & 0x01);
                        tmpbuffer[jj-6] = ((uctmp >> 0x06) & 0x01);
                        tmpbuffer[jj-7] = ((uctmp >> 0x07) & 0x01);
                    }
                }
            }
            last_offs_sub = offs_sub;
        } else {
            // if not, use padding instead
            double dnumblocks;
            int numblocks;
            dnumblocks = (offs_sub-last_offs_sub)/S.time_per_subint;
            numblocks = (int) round(dnumblocks);
            //printf("\n%d  %20.15f  %20.15f  %20.15f  %20.15f  %20.15f  %d \n", 
            //       cur_subint, last_offs_sub, offs_sub, 
            //       offs_sub-last_offs_sub, S.time_per_subint, dnumblocks, numblocks);
            missing_blocks++;
            if (fabs(dnumblocks - (double)numblocks) > 1e-6) {
                printf("\nYikes!  We missed a fraction (%.20f) of a subint!\n", 
                       fabs(dnumblocks - (double)numblocks));
            }
            printf("At subint %d found %d dropped subints (%d total), adding 1.\n", 
                   cur_subint, numblocks-1, missing_blocks);
            // Add a full block of padding
            numtopad = S.spectra_per_subint * S.num_polns;
            for (ii = 0; ii < numtopad; ii++)
                memcpy(tmpbuffer + ii * S.bytes_per_spectra/S.num_polns, 
                       newpadvals, S.bytes_per_spectra/S.num_polns);
            // Update the time of the last subs based on this padding
            last_offs_sub += S.time_per_subint;
            // Set the padding flag
            *padding = 1;
            // Decrement the subint counter since it gets incremented later
            // but we don't want to move to the next subint yet
            cur_subint--;
        }

        // This loop allows us to work with single polns out of many
        // or to sum polarizations if required
        if (S.num_polns > 1) {
            int sum_polns = 0;

            if ((0==strncmp(S.poln_order, "AABB", 4)) || (S.num_polns == 2)) sum_polns = 1;
            if (user_poln || ((S.num_polns > 2) && !sum_polns)) {  // The user chose the poln
                int ii, offset;
                unsigned char *tmpptr = dataptr;
                for (ii = 0 ; ii < S.spectra_per_subint ; ii++) {
                    offset = ii * S.samples_per_spectra + default_poln * S.num_channels;
                    memcpy(tmpptr+ii*S.num_channels, tmpbuffer+offset, S.num_channels);
                }
            } else if (sum_polns) { // sum the polns if there are 2 by default
                int ii, jj, itmp, offset0, offset1;
                unsigned char *tmpptr = dataptr;
                for (ii = 0 ; ii < S.spectra_per_subint ; ii++) {
                    offset0 = ii * S.samples_per_spectra;
                    offset1 = offset0 + S.num_channels;
                    for (jj = 0 ; jj < S.num_channels ; jj++) {
                        itmp = (int) tmpbuffer[offset0+jj] + (int) tmpbuffer[offset1+jj];
                        tmpptr[ii*S.num_channels+jj] = (unsigned char) (itmp >> 1);
                    }
                }
            }
        }

        if (firsttime) {
            // Determine overall scaling of the data so it will fit
            // into an unsigned char.  If we used _floats_ we wouldn't
            // need to do that!  TODO:  fix this!
            if (S.apply_offset || S.apply_scale) {
                int ii, jj, d_idx, os_idx;
                unsigned char *tmpptr;
                float *fvec, offset, scale, fmed;
                fvec = gen_fvect(S.spectra_per_subint*S.num_channels);
                for (ii = 0 ; ii < S.spectra_per_subint ; ii++) {
                    d_idx = ii * S.num_channels;
                    os_idx = default_poln * S.num_channels;
                    tmpptr = dataptr + d_idx;
                    for (jj = 0 ; jj < S.num_channels ; jj++, os_idx++) {
                        offset = (S.apply_offset) ? offsets[os_idx] : 0.0;
                        scale = (S.apply_scale) ? scales[os_idx] : 1.0;
                        fvec[d_idx+jj] = (((float)tmpptr[jj] - S.zero_offset) * scale) + offset;
                    }
                }
                // Now determine the median of the data...
                fmed = median(fvec, S.spectra_per_subint*S.num_channels);
                // Set the scale so that the median is at about 1/3 of the
                // dynamic range of an unsigned char.  Note that this assumes
                // that the data are properly offset so that the min values
                // are at values of zero...
                global_scale = (256.0/3.0) / fmed;
                printf("\nSetting PSRFITS global scale to %f\n", global_scale);
                vect_free(fvec);

            }
            firsttime = 0;
        }

        // Apply offsets and scales if needed
        if (S.apply_offset || S.apply_scale) {
            int ii, jj, d_idx, os_idx;
            unsigned char *tmpptr;
            float ftmp, offset, scale;
            for (ii = 0 ; ii < S.spectra_per_subint ; ii++) {
                d_idx = ii * S.num_channels;
                os_idx = default_poln * S.num_channels;
                tmpptr = dataptr + d_idx;
                for (jj = 0 ; jj < S.num_channels ; jj++, os_idx++) {
                    offset = (S.apply_offset) ? offsets[os_idx] : 0.0;
                    scale = (S.apply_scale) ? scales[os_idx] : 1.0;
                    ftmp = (((float)tmpptr[jj] - S.zero_offset) * scale) + offset;
                    ftmp = (ftmp < 0.0) ? 0.0 : ftmp;
                    tmpptr[jj] = (unsigned char)(ftmp * global_scale + 0.5);
                }
            }
        }

        // Apply weights if needed
        if (S.apply_weight) {
            int ii, jj, offset;
            unsigned char *tmpptr;
            float ftmp;
            for (ii = 0 ; ii < S.spectra_per_subint ; ii++) {
                offset = ii * S.num_channels;
                tmpptr = dataptr + offset;
                for (jj = 0 ; jj < S.num_channels ; jj++) {
                    ftmp = (float)tmpptr[jj] * weights[jj] + 0.5;
                    tmpptr[jj] = (unsigned char)ftmp;
                }
            }
        }

        // Flip the band if needed
        if (S.apply_flipband) {
            unsigned char uctmp;
            int ii, jj, offset;
            for (jj = 0 ; jj < S.spectra_per_subint ; jj++) {
                offset = jj * S.num_channels;
                for (ii = 0 ; ii < S.num_channels/2 ; ii++) {
                    uctmp = dataptr[offset+ii];
                    dataptr[offset+ii] = dataptr[offset+S.num_channels-ii-1];
                    dataptr[offset+S.num_channels-ii-1] = uctmp;
                }
            }
        }

        // Hack to flip each byte of data if needed
        if (S.flip_bytes) {
            unsigned char uctmp;
            int ii, jj;
            for (ii = 0 ; ii < S.bytes_per_subint/8 ; ii++) {
                offset = ii * 8;
                for (jj = 0 ; jj < 4 ; jj++) {
                    uctmp = dataptr[offset+jj];
                    dataptr[offset+jj] = dataptr[offset+8-1-jj];
                    dataptr[offset+8-1-jj] = uctmp;
                }
            }
        }

        if (status) {
            printf("\nProblem reading record from PSRFITS data file:\n");
            printf("   currentfile = %d, currentsubint = %d.  Exiting.\n",
                   cur_file, cur_subint);
            exit(1);
        }
        // Clip nasty RFI if requested
        if (S.clip_sigma > 0.0)  {
            clip_times(dataptr, S.spectra_per_subint, S.num_channels, 
                       S.clip_sigma, newpadvals);
        }
        *padding = 0;
        // Pull the new data from the ringbuffer if there is an offset
        if (bufferspec)
            memcpy(data, ringbuffer, S.bytes_per_subint/S.num_polns);
        cur_subint++;
        return 1;
        
    
    } else { // We can't read anymore...  so read OFFS_SUB
        // for the last row of the current file to see about padding
        fits_read_col(S.files[cur_file], TDOUBLE, 
                      S.offs_sub_col, S.num_subint[cur_file], 1L, 1L, 
                      0, &offs_sub, &anynull, &status);
    }

    if (S.num_pad[cur_file]==0 ||
        TEST_CLOSE(last_offs_sub, offs_sub)) {  // No padding is necessary
        // The TEST_CLOSE check means that the lack of data we noticed
        // upon reading the file was due to dropped data in the
        // middle of the file that we already fixed.  So no
        // padding is really necessary.
        cur_file++;
        cur_subint = 1;
        shiftbuffer = 0;  // Since recursively calling, don't shift again
        return read_PSRFITS_rawblock(data, padding);
    } else { // add padding
        numtopad = S.num_pad[cur_file] - padnum;
        // The amount of padding still to be sent for this file
        if (numtopad) {
            *padding = 1;
            if (numtopad >= S.spectra_per_subint - bufferspec) {
                if (bufferspec) {  // Buffer the padding?
                    // Add the amount of padding we need to
                    // make our buffer offset = 0
                    numtopad = S.spectra_per_subint - bufferspec;
                    for (ii = 0; ii < numtopad; ii++)
                        memcpy(dataptr + ii * S.bytes_per_spectra/S.num_polns, 
                               newpadvals, S.bytes_per_spectra/S.num_polns);
                    // Copy the new data/padding into the output array
                    memcpy(data, ringbuffer, S.bytes_per_subint/S.num_polns);
                    bufferspec = 0;
                } else {  // Add a full record of padding
                    numtopad = S.spectra_per_subint;
                    for (ii = 0; ii < numtopad; ii++)
                        memcpy(data + ii * S.bytes_per_spectra/S.num_polns, 
                               newpadvals, S.bytes_per_spectra/S.num_polns);
                }
                padnum += numtopad;
                cur_subint++;
                // If done with padding reset padding variables and go to next file
                if (padnum == S.num_pad[cur_file]) {
                    padnum = 0;
                    cur_file++;
                    cur_subint = 1;
                }
                // Update the time of the last subs based on this padding
                last_offs_sub += S.dt * numtopad;
                return 1;
            } else {  // Need < 1 block (or remaining block) of padding
                int pad;
                // Add the remainder of the padding and then get a
                // block from the next file.
                for (ii = 0; ii < numtopad; ii++) {
                    offset = bufferspec * S.bytes_per_spectra/S.num_polns;
                    memcpy(ringbuffer + offset + ii * S.bytes_per_spectra/S.num_polns,
                           newpadvals, S.bytes_per_spectra/S.num_polns);
                }
                bufferspec += numtopad;
                padnum = 0;
                cur_file++;
                cur_subint = 1;
                shiftbuffer = 0;  // Since recursively calling, don't shift again
                // Update the time of the last subs based on this padding
                last_offs_sub += S.dt * numtopad;
                return read_PSRFITS_rawblock(data, &pad);
            }
        }
    }
    return 0;  // Should never get here
}


int read_PSRFITS_rawblocks(unsigned char rawdata[], int numblocks,
                           int *padding)
// This routine reads numblocks PSRFITS blocks (i.e. subints) from the
// input fils.  The raw data is returned in rawdata which must have a
// size of numblocks* S.bytes_per_subint.  The number of blocks read
// is returned.  If padding is returned as 1, then padding was added
// and statistics should not be calculated.
{
   int ii, retval = 0, pad, numpad = 0;
   
   *padding = 0;
   for (ii = 0; ii < numblocks; ii++) {
       pad = 0;
       retval += read_PSRFITS_rawblock(rawdata + ii * S.samples_per_subint/S.num_polns, &pad);
       if (pad) numpad++;
   }
   if (numpad)
       *padding = 1;
   return retval;
}


int read_PSRFITS(float *data, int numspec, double *dispdelays, int *padding,
                 int *maskchans, int *nummasked, mask * obsmask)
// This routine reads numspec from the PSRFITS raw input files.  These
// files contain raw data in PSRFITS format.  Time delays and a mask
// are applied to each channel.  It returns the # of points read if
// successful, 0 otherwise.  If padding is returned as 1, then padding
// was added and statistics should not be calculated.  maskchans is an
// array of length numchans contains a list of the number of channels
// that were masked.
{
   int ii, jj, numread = 0, offset;
   static unsigned char *tempzz, *rawdata1, *rawdata2;
   static unsigned char *currentdata, *lastdata;
   static int firsttime = 1, numblocks = 1, allocd = 0, mask = 0;
   static double duration = 0.0;

   *nummasked = 0;
   if (firsttime) {
       if (numspec % S.spectra_per_subint) {
           printf("numspec must be a multiple of %d in read_PSRFITS()!\n", 
                  S.spectra_per_subint);
           exit(1);
       } else
           numblocks = numspec / S.spectra_per_subint;
       
       if (obsmask->numchan)
           mask = 1;
       rawdata1 = gen_bvect(numblocks * S.bytes_per_subint/S.num_polns);
       rawdata2 = gen_bvect(numblocks * S.bytes_per_subint/S.num_polns);
       allocd = 1;
       duration = numblocks * S.time_per_subint;
       currentdata = rawdata1;
       lastdata = rawdata2;
   }

   // Read and de-disperse
   if (allocd) {
       while (1) {
           numread = read_PSRFITS_rawblocks(currentdata, numblocks, padding);
           if (mask) {
               // Remember that last_offs_sub gets updated before 
               // read_PSRFITS_rawblock returns.  And also, offs_sub is the
               // midpoint of each subint.  (note:  this is only correct 
               // if numblocks is 1, which it should be, I think)
               *nummasked = check_mask(last_offs_sub - 0.5 * duration - \
                                       S.start_subint[0] * S.time_per_subint, 
                                       duration, obsmask, maskchans);
           }
           // Only use the recently measured padding if all the channels aren't masked
           if ((S.clip_sigma > 0.0) && 
               !(mask && (*nummasked == -1)) &&
               (padvals != newpadvals))
               memcpy(padvals, newpadvals, S.bytes_per_spectra/S.num_polns);
           if (mask) {
               //if (S.num_polns > 1 && !S.summed_polns) {
               //    printf("WARNING!:  masking does not currently work with multiple polns!\n");
               //}
               if (*nummasked == -1) {  // If all channels are masked
                   for (ii = 0; ii < numspec; ii++)
                       memcpy(currentdata + ii * S.bytes_per_spectra/S.num_polns, 
                              padvals, S.bytes_per_spectra/S.num_polns);
               } else if (*nummasked > 0) { // Only some of the channels are masked
                   int channum;
                   for (ii = 0; ii < numspec; ii++) {
                       offset = ii * S.num_channels;
                       for (jj = 0; jj < *nummasked; jj++) {
                           channum = maskchans[jj];
                           currentdata[offset + channum] = padvals[channum];
                       }
                   }
               }
           }

           if (!firsttime)
               dedisp(currentdata, lastdata, numspec, S.num_channels, dispdelays, data);
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
       return numblocks * S.spectra_per_subint;
   } else {
       return 0;
   }
}


void get_PSRFITS_channel(int channum, float chandat[],
                         unsigned char rawdata[], int numblocks)
// Return the values for channel 'channum' of a block of 'numblocks'
// raw PSRFITS data stored in 'rawdata' in 'chandat'.  'rawdata'
// should have been initialized using read_PSRFITS_rawblocks(), and
// 'chandat' must have at least 'numblocks' * S.samples_per_subint
// spaces.  Channel 0 is assumed to be the lowest freq channel.
{
    int ii, jj, ptsperchan;
    
    if (channum > S.num_channels || channum < 0) {
        printf("\nchannum = %d is out of range in get_PSRFITS_channel()!\n\n", 
               channum);
        exit(1);
    }
    ptsperchan = S.spectra_per_subint * numblocks;
    
    // Since the following is only called from rfifind, we know that the
    // channel accesses will be in order from 0 to numchan-1
    if (channum == 0) {          // Transpose the data
        short trtn;
        int move_size;
        unsigned char *move;
        
        move_size = (ptsperchan + S.num_channels) / 2;
        move = gen_bvect(move_size);
        if ((trtn = transpose_bytes(rawdata, ptsperchan, S.num_channels,
                                    move, move_size)) < 0)
            printf("Error %d in transpose_bytes().\n", trtn);
        vect_free(move);
    }

    // Select the correct channel
    for (ii = 0, jj = ptsperchan * channum; ii < ptsperchan; ii++, jj++)
        chandat[ii] = (float) rawdata[jj];
}


int prep_PSRFITS_subbands(unsigned char *rawdata, float *data,
                          double *dispdelays, int numsubbands,
                          int transpose, int *maskchans,
                          int *nummasked, mask * obsmask)
// This routine preps a block from PSRFITS files.  The routine uses
// dispersion delays in 'dispdelays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'data' of length 'numsubbands' * S.spectra_per_subint.  The low
// freq subband is stored first, then the next highest subband etc,
// with S.spectra_per_subint floating points per subband.  It returns
// the # of points read if succesful, 0 otherwise.  'maskchans' is an
// array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    int ii, jj, trtn, offset;
    static unsigned char *tempzz, *rawdata1, *rawdata2;
    static unsigned char *currentdata, *lastdata, *move;
    static int firsttime = 1, move_size = 0, mask = 0;
    
    *nummasked = 0;
    if (firsttime) {
        if (obsmask->numchan)
            mask = 1;
        rawdata1 = gen_bvect(S.bytes_per_subint/S.num_polns);
        rawdata2 = gen_bvect(S.bytes_per_subint/S.num_polns);
        move_size = (S.spectra_per_subint + numsubbands) / 2;
        move = gen_bvect(move_size);
        currentdata = rawdata1;
        lastdata = rawdata2;
        memcpy(currentdata, rawdata, S.bytes_per_subint/S.num_polns);
    }
    
    // Read and de-disperse
    memcpy(currentdata, rawdata, S.bytes_per_subint/S.num_polns);
    if (mask) {
        // Remember that last_offs_sub gets updated before 
        // read_PSRFITS_rawblock returns.  And also, offs_sub is the
        // midpoint of each subint.
        *nummasked = check_mask(last_offs_sub - 0.5 * S.time_per_subint - \
                                S.start_subint[0] * S.time_per_subint, 
                                S.time_per_subint, obsmask, maskchans);
    }
    
    // Only use the recently measured padding if all the channels aren't masked
    if ((S.clip_sigma > 0.0) && 
        !(mask && (*nummasked == -1)) &&
        (padvals != newpadvals))
        memcpy(padvals, newpadvals, S.bytes_per_spectra/S.num_polns);
    
    if (mask) {
        if (S.num_polns > 1 && !S.summed_polns) {
            printf("WARNING!:  masking does not currently work with multiple polns!\n");
        }
        
        if (*nummasked == -1) {  // If all channels are masked
            for (ii = 0; ii < S.spectra_per_subint; ii++)
                memcpy(currentdata + ii * S.num_channels, padvals, S.num_channels);
        } else if (*nummasked > 0) {  // Only some of the channels are masked
            int channum;
            for (ii = 0; ii < S.spectra_per_subint; ii++) {
                offset = ii * S.num_channels;
                for (jj = 0; jj < *nummasked; jj++) {
                    channum = maskchans[jj];
                    currentdata[offset + channum] = padvals[channum];
                }
            }
        }
    }
    
    // In mpiprepsubband, the nodes do not call read_*_rawblock()
    // where cur_subint gets incremented.
    if (using_MPI) {
        cur_subint++;
        last_offs_sub += S.time_per_subint;
    }
    
    if (firsttime) {
        SWAP(currentdata, lastdata);
        firsttime = 0;
        return 0;
    } else {
        dedisp_subbands(currentdata, lastdata, S.spectra_per_subint, S.num_channels,
                        dispdelays, numsubbands, data);
        SWAP(currentdata, lastdata);
        /* Transpose the data into vectors in the result array */
        if (transpose) {
            if ((trtn = transpose_float(data, S.spectra_per_subint, numsubbands,
                                        move, move_size)) < 0)
                printf("Error %d in transpose_float().\n", trtn);
        }
        return S.spectra_per_subint;
    }
}


int read_PSRFITS_subbands(float *data, double *dispdelays, int numsubbands,
                          int transpose, int *padding,
                          int *maskchans, int *nummasked, mask * obsmask)
// This routine reads a record from the input files which contain
// PSRFITS data.  The routine uses dispersion delays in 'dispdelays'
// to de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'data' of length 'numsubbands' *
// S.spectra_per_subint.  The low freq subband is stored first, then
// the next highest subband etc, with S.spectra_per_subint floating
// points per subband. It returns the # of points read if succesful, 0
// otherwise.  If padding is returned as 1, then padding was added and
// statistics should not be calculated.  'maskchans' is an array of
// length numchans which contains a list of the number of channels
// that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static unsigned char *rawdata;
    
    if (firsttime) {
        rawdata = gen_bvect(S.bytes_per_subint/S.num_polns);
        if (!read_PSRFITS_rawblock(rawdata, padding)) {
            printf("Problem reading the raw PSRFITS data.\n\n");
            return 0;
        }
        if (0 != prep_PSRFITS_subbands(rawdata, data, dispdelays, numsubbands,
                                       transpose, maskchans, nummasked, obsmask)) {
            printf("Problem initializing prep_PSRFITS_subbands()\n\n");
            return 0;
        }
        firsttime = 0;
    }
    if (!read_PSRFITS_rawblock(rawdata, padding)) {
        return 0;
    }
    return prep_PSRFITS_subbands(rawdata, data, dispdelays, numsubbands,
                                 transpose, maskchans, nummasked, obsmask);
}
