#include <sys/types.h>
#include <pwd.h>
#include "presto.h"
#include "mask.h"
#include "psrfits.h"

#define DEBUG_OUT 1

static unsigned char *cdatabuffer;
static float *fdatabuffer, *offsets, *scales, *weights;
static int cur_file = 0, cur_subint = 1, numbuffered = 0, offs_sub_are_zero = 0;
static long long cur_spec = 0, new_spec = 0;

extern double slaCldj(int iy, int im, int id, int *j);
extern void add_padding(float *fdata, float *padding, int numchan, int numtopad);

void get_PSRFITS_subint(float *fdata, unsigned char *cdata, struct spectra_info *s);

double DATEOBS_to_MJD(char *dateobs, int *mjd_day, double *mjd_fracday)
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
    int status = 0;
    char ctmp[80], comment[120], err_text[81];

    // Read the primary HDU
    fits_open_file(&fptr, filename, READONLY, &status);
    if (status) {
        fits_get_errstatus(status, err_text);
        printf("Error %d opening %s : %s\n", status, filename, err_text);
        fits_close_file(fptr, &status);
        return 0;
    }

    // Make the easy check first
    fits_read_key(fptr, TSTRING, "FITSTYPE", ctmp, comment, &status);
    if (status) {
        fits_get_errstatus(status, err_text);
        printf("Error %d reading 'FITSTYPE' from %s : %s\n",
               status, filename, err_text);
        fits_close_file(fptr, &status);
        return 0;
    } else {
        if (strcmp(ctmp, "PSRFITS")) {
            printf("Error 'FITSTYPE' is not 'PSRFITS' in %s\n", filename);
            fits_close_file(fptr, &status);
            return 0;
        }
    }

    // See if the data are search-mode
    fits_read_key(fptr, TSTRING, "OBS_MODE", ctmp, comment, &status);
    if (status) {
        fits_get_errstatus(status, err_text);
        printf("Error %d reading 'OBS_MODE' from %s : %s\n",
               status, filename, err_text);
        fits_close_file(fptr, &status);
        return 0;
    } else {
        if ((strcmp(ctmp, "SEARCH") && strcmp(ctmp, "SRCH"))) {
            printf("Error 'OBS_MODE' is not 'SEARCH' in %s\n", filename);
            fits_close_file(fptr, &status);
            return 0;
        }
    }
    fits_close_file(fptr, &status);
    return 1;                   // it is search-mode  PSRFITS
}

#define check_read_status(name) {                                   \
        if (status) {\
            fits_get_errstatus(status, err_text); \
            printf("Error %d reading %s : %s\n", status, name, err_text); \
            status=0;      \
        }                                                               \
    }

#define get_hdr_string(name, param) {                                   \
        fits_read_key(s->fitsfiles[ii], TSTRING, (name), ctmp, comment, &status); \
        if (status) {\
            fits_get_errstatus(status, err_text); \
            printf("Error %d reading %s : %s\n", status, name, err_text); \
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
        fits_read_key(s->fitsfiles[ii], TINT, (name), &itmp, comment, &status); \
        if (status) {\
            fits_get_errstatus(status, err_text); \
            printf("Error %d reading %s : %s\n", status, name, err_text); \
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
        fits_read_key(s->fitsfiles[ii], TDOUBLE, (name), &dtmp, comment, &status); \
        if (status) {\
            fits_get_errstatus(status, err_text); \
            printf("Error %d reading %s : %s\n", status, name, err_text); \
            if (ii==0.0) param=0.0; \
            if (status==KEY_NO_EXIST) status=0;\
        } else {                                                          \
            if (ii==0) param = dtmp;                                        \
            else if (param != dtmp)                                         \
                printf("Warning!:  %s values don't match for files 0 and %d!\n", \
                       (name), ii);                                         \
        }                                                               \
    }



void read_PSRFITS_files(struct spectra_info *s)
// Read and convert PSRFITS information from a group of files
// and place the resulting info into a spectra_info structure.
{
    int IMJD, SMJD, itmp, ii, status = 0;
    double OFFS, BE_DELAY, dtmp;
    long double MJDf;
    char ctmp[80], comment[120], err_text[81];

    s->datatype = PSRFITS;
    s->fitsfiles = (fitsfile **) malloc(sizeof(fitsfile *) * s->num_files);
    s->start_subint = gen_ivect(s->num_files);
    s->num_subint = gen_ivect(s->num_files);
    s->start_spec = (long long *) malloc(sizeof(long long) * s->num_files);
    s->num_spec = (long long *) malloc(sizeof(long long) * s->num_files);
    s->num_pad = (long long *) malloc(sizeof(long long) * s->num_files);
    s->start_MJD = (long double *) malloc(sizeof(long double) * s->num_files);
    s->N = 0;
    s->num_beams = 1;
    s->get_rawblock = &get_PSRFITS_rawblock;
    s->offset_to_spectra = &offset_to_PSRFITS_spectra;

    // By default, don't flip the band.  But don't change
    // the input value if it is aleady set to flip the band always
    if (s->apply_flipband == -1)
        s->apply_flipband = 0;

    // Step through the other files
    for (ii = 0; ii < s->num_files; ii++) {

        // Is the file a PSRFITS file?
        if (!is_PSRFITS(s->filenames[ii])) {
            fprintf(stderr,
                    "\nError!  File '%s' does not appear to be PSRFITS!\n",
                    s->filenames[ii]);
            exit(1);
        }
        // Open the PSRFITS file
        fits_open_file(&(s->fitsfiles[ii]), s->filenames[ii], READONLY, &status);

        // Is the data in search mode?
        fits_read_key(s->fitsfiles[ii], TSTRING, "OBS_MODE", ctmp, comment, &status);
        check_read_status("OBS_MODE");
        // Quick fix for Parkes DFB data (SRCH?  why????)...
        if (strcmp("SRCH", ctmp) == 0) {
            strncpy(ctmp, "SEARCH", 40);
        }
        if (strcmp(ctmp, "SEARCH")) {
            fprintf(stderr,
                    "\nError!  File '%s' does not contain SEARCH-mode data!\n",
                    s->filenames[ii]);
            exit(1);
        }
        // Now get the stuff we need from the primary HDU header
        fits_read_key(s->fitsfiles[ii], TSTRING, "TELESCOP", ctmp, comment, &status);
        // Quick fix for MockSpec data...
        if (strcmp("ARECIBO 305m", ctmp) == 0) {
            strncpy(ctmp, "Arecibo", 40);
        }
        // Quick fix for Parkes DFB data...
        {
            char newctmp[80];

            // Copy ctmp first since strlower() is in-place
            strcpy(newctmp, ctmp);
            if (strcmp("parkes", strlower(remove_whitespace(newctmp))) == 0) {
                strncpy(ctmp, "Parkes", 40);
            }
        }
        if (status) {
            printf("Error %d reading key %s\n", status, "TELESCOP");
            if (ii == 0)
                s->telescope[0] = '\0';
            if (status == KEY_NO_EXIST)
                status = 0;
        } else {
            if (ii == 0)
                strncpy(s->telescope, ctmp, 40);
            else if (strcmp(s->telescope, ctmp) != 0)
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

        /* This is likely not in earlier versions of PSRFITS */
        s->chan_dm = 0.0;
        fits_read_key(s->fitsfiles[ii], TDOUBLE, "CHAN_DM",
                      &(s->chan_dm), comment, &status);
        if (status==KEY_NO_EXIST) status=0; // Prevents error messages on old files
        check_read_status("CHAN_DM");
        // Don't use the macros unless you are using the struct!
        fits_read_key(s->fitsfiles[ii], TINT, "STT_IMJD", &IMJD, comment, &status);
        check_read_status("STT_IMJD");
        s->start_MJD[ii] = (long double) IMJD;
        fits_read_key(s->fitsfiles[ii], TINT, "STT_SMJD", &SMJD, comment, &status);
        check_read_status("STT_SMJD");
        fits_read_key(s->fitsfiles[ii], TDOUBLE, "STT_OFFS", &OFFS, comment,
                      &status);
        check_read_status("STT_OFFS");
        BE_DELAY = 0.0; // Back-end delay.  Will only be applied to STT*-based times
        fits_read_key(s->fitsfiles[ii], TDOUBLE, "BE_DELAY", &BE_DELAY, comment,
                      &status);
        if (status==KEY_NO_EXIST) status=0; // Prevents error messages on old files
        check_read_status("BE_DELAY");
        s->start_MJD[ii] += ((long double) SMJD +
                             (long double) OFFS +
                             (long double) BE_DELAY) / SECPERDAY;

        // Are we tracking?
        fits_read_key(s->fitsfiles[ii], TSTRING, "TRK_MODE", ctmp, comment, &status);
        check_read_status("TRK_MODE");
        itmp = (strcmp("TRACK", ctmp) == 0) ? 1 : 0;
        if (ii == 0)
            s->tracking = itmp;
        else if (s->tracking != itmp)
            printf("Warning!:  TRK_MODE values don't match for files 0 and %d!\n",
                   ii);

        // Now switch to the SUBINT HDU header
        fits_movnam_hdu(s->fitsfiles[ii], BINARY_TBL, "SUBINT", 0, &status);
        check_read_status("SUBINT");
        get_hdr_double("TBIN", s->dt);
        get_hdr_int("NCHAN", s->num_channels);
        get_hdr_int("NPOL", s->num_polns);
        get_hdr_string("POL_TYPE", s->poln_order);
        fits_read_key(s->fitsfiles[ii], TINT, "NCHNOFFS", &itmp, comment, &status);
        check_read_status("NCHNOFFS");
        if (itmp > 0)
            printf("Warning!:  First freq channel is not 0 in file %d!\n", ii);
        get_hdr_int("NSBLK", s->spectra_per_subint);
        get_hdr_int("NBITS", s->bits_per_sample);
        fits_read_key(s->fitsfiles[ii], TINT, "NAXIS2",
                      &(s->num_subint[ii]), comment, &status);
        check_read_status("NAXIS2");
        fits_read_key(s->fitsfiles[ii], TINT, "NSUBOFFS",
                      &(s->start_subint[ii]), comment, &status);
        check_read_status("NSUBOFFS");
        s->time_per_subint = s->dt * s->spectra_per_subint;

        /* This is likely not in earlier versions of PSRFITS */
        s->zero_offset = 0.0;
        fits_read_key(s->fitsfiles[ii], TFLOAT, "ZERO_OFF",
                      &(s->zero_offset), comment, &status);
        if (status==KEY_NO_EXIST) status=0; // Prevents error messages on old files
        check_read_status("ZERO_OFF");
        s->zero_offset = fabs(s->zero_offset);

        // Get the time offset column info and the offset for the 1st row
        {
            double offs_sub;
            int colnum, anynull, numrows;

            // Identify the OFFS_SUB column number
            fits_get_colnum(s->fitsfiles[ii], 0, "OFFS_SUB", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the OFFS_SUB column!\n");
                status = 0;     // Reset status
            } else {
                if (ii == 0) {
                    s->offs_sub_col = colnum;
                } else if (colnum != s->offs_sub_col) {
                    printf("Warning!:  OFFS_SUB column changes between files!\n");
                }
            }

            // Read the OFFS_SUB column value for the 1st row
            fits_read_col(s->fitsfiles[ii], TDOUBLE,
                          s->offs_sub_col, 1L, 1L, 1L,
                          0, &offs_sub, &anynull, &status);

            if (offs_sub != 0.0) {
                numrows = (int) ((offs_sub - 0.5 * s->time_per_subint) /
                                 s->time_per_subint + 1e-7);
                // Check to see if any rows have been deleted or are missing
                if (numrows > s->start_subint[ii]) {
                    printf("Warning!:  NSUBOFFS reports %d previous rows\n"
                           "           but OFFS_SUB implies %d.  Using OFFS_SUB.\n"
                           "           Will likely be able to correct for this.\n",
                           s->start_subint[ii], numrows);
                }
                s->start_subint[ii] = numrows;
            } else {
                int indexval_col, jj;
                offs_sub_are_zero = 1;

                // If OFFS_SUB are all 0.0, then we will assume that there are
                // no gaps in the file.  This isn't truly proper PSRFITS, but
                // we should still be able to handle it
                for (jj = 1; jj <= s->num_subint[ii]; jj++) {
                    fits_read_col(s->fitsfiles[ii], TDOUBLE,
                                  s->offs_sub_col, jj, 1L, 1L,
                                  0, &offs_sub, &anynull, &status);
                    if (offs_sub != 0.0) {
                        offs_sub_are_zero = 0;
                        perror
                            ("Error!:  Some, but not all OFFS_SUB are 0.0.  Not good PSRFITS.\n");
                        exit(EXIT_FAILURE);
                    }
                }
                if (offs_sub_are_zero) {
                    printf
                        ("Warning!:  All OFFS_SUB are 0.0.  Assuming no missing rows.\n");
                }
                // Check to see if there is an INDEXVAL column.  That should tell
                // us if we are missing any subints.  Use it in lieu of OFFS_SUB
                fits_get_colnum(s->fitsfiles[ii], 0, "INDEXVAL", &indexval_col,
                                &status);
                if (status == COL_NOT_FOUND) {
                    printf
                        ("Warning!:  No INDEXVAL column, either.  This is not proper PSRFITS.\n");
                    status = 0; // Reset status
                    s->start_subint[ii] = 0;
                } else {
                    double subint_index;
                    // Read INDEXVAL
                    fits_read_col(s->fitsfiles[ii], TDOUBLE,
                                  indexval_col, 1L, 1L, 1L,
                                  0, &subint_index, &anynull, &status);
                    s->start_subint[ii] = (int) (subint_index + 1e-7 - 1.0);
                }
            }
        }

        // This is the MJD offset based on the starting subint number
        MJDf = (s->time_per_subint * s->start_subint[ii]) / SECPERDAY;
        // The start_MJD values should always be correct
        s->start_MJD[ii] += MJDf;

        // Compute the starting spectra from the times
        MJDf = s->start_MJD[ii] - s->start_MJD[0];
        if (MJDf < 0.0) {
            fprintf(stderr, "Error!: File %d seems to be from before file 0!\n", ii);
            exit(1);
        }
        s->start_spec[ii] = (long long) (MJDf * SECPERDAY / s->dt + 0.5);

        // Now pull stuff from the other columns
        {
            float ftmp;
            long repeat, width;
            int colnum, anynull;

            // Identify the data column and the data type
            fits_get_colnum(s->fitsfiles[ii], 0, "DATA", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the DATA column!\n");
                status = 0;     // Reset status
            } else {
                if (ii == 0) {
                    s->data_col = colnum;
                    fits_get_coltype(s->fitsfiles[ii], colnum, &(s->FITS_typecode),
                                     &repeat, &width, &status);
                    // This makes CFITSIO treat 1-bit data as written in 'B' mode
                    // even if it was written in 'X' mode originally.  This means
                    // that we unpack it ourselves.
                    if (s->bits_per_sample < 8 && s->FITS_typecode == 1) {
                        s->FITS_typecode = 11;
                    }
                } else if (colnum != s->data_col) {
                    printf("Warning!:  DATA column changes between files!\n");
                }
            }

            // Telescope azimuth
            fits_get_colnum(s->fitsfiles[ii], 0, "TEL_AZ", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                s->azimuth = 0.0;
                status = 0;     // Reset status
            } else {
                fits_read_col(s->fitsfiles[ii], TFLOAT, colnum,
                              1L, 1L, 1L, 0, &ftmp, &anynull, &status);
                if (ii == 0)
                    s->azimuth = (double) ftmp;
            }

            // Telescope zenith angle
            fits_get_colnum(s->fitsfiles[ii], 0, "TEL_ZEN", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                s->zenith_ang = 0.0;
                status = 0;     // Reset status
            } else {
                fits_read_col(s->fitsfiles[ii], TFLOAT, colnum,
                              1L, 1L, 1L, 0, &ftmp, &anynull, &status);
                if (ii == 0)
                    s->zenith_ang = (double) ftmp;
            }

            // Observing frequencies
            fits_get_colnum(s->fitsfiles[ii], 0, "DAT_FREQ", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel freq column!\n");
                status = 0;     // Reset status
            } else {
                int jj;
                double *freqs = (double *) malloc(sizeof(double) * s->num_channels);
                fits_read_col(s->fitsfiles[ii], TDOUBLE, colnum, 1L, 1L,
                              s->num_channels, 0, freqs, &anynull, &status);

                if (ii == 0) {
                    int trigger = 0;
                    s->df = ((double) freqs[s->num_channels - 1] -
                             (double) freqs[0]) / (double) (s->num_channels - 1);
                    s->lo_freq = freqs[0];
                    s->hi_freq = freqs[s->num_channels - 1];
                    // Now check that the channel spacing is the same throughout
                    for (jj = 0; jj < s->num_channels - 1; jj++) {
                        ftmp = freqs[jj + 1] - freqs[jj];
                        if ((fabs(ftmp - s->df) > 1e-7) && !trigger) {
                            trigger = 1;
                            printf
                                ("Warning!:  Channel spacing changes in file %d!\n",
                                 ii);
                        }
                    }
                } else {
                    ftmp = fabs(s->df - (freqs[1] - freqs[0]));
                    if (ftmp > 1e-7)
                        printf
                            ("Warning!:  Channel spacing changes between files!\n");
                    ftmp = fabs(s->lo_freq - freqs[0]);
                    if (ftmp > 1e-7)
                        printf("Warning!:  Low channel changes between files!\n");
                    ftmp = fabs(s->hi_freq - freqs[s->num_channels - 1]);
                    if (ftmp > 1e-7)
                        printf("Warning!:  High channel changes between files!\n");
                }
                free(freqs);
            }

            // Data weights
            fits_get_colnum(s->fitsfiles[ii], 0, "DAT_WTS", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel weights!\n");
                status = 0;     // Reset status
            } else {
                if (s->apply_weight < 0) {      // Use the data to decide
                    int jj;
                    if (ii == 0) {
                        s->dat_wts_col = colnum;
                    } else if (colnum != s->dat_wts_col) {
                        printf("Warning!:  DAT_WTS column changes between files!\n");
                    }
                    float *fvec = (float *) malloc(sizeof(float) * s->num_channels);
                    fits_read_col(s->fitsfiles[ii], TFLOAT, s->dat_wts_col, 1L, 1L,
                                  s->num_channels, 0, fvec, &anynull, &status);
                    for (jj = 0; jj < s->num_channels; jj++) {
                        // If the weights are not 1, apply them
                        if (fvec[jj] != 1.0) {
                            s->apply_weight = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_weight < 0)
                    s->apply_weight = 0;        // not needed
            }

            // Data offsets
            fits_get_colnum(s->fitsfiles[ii], 0, "DAT_OFFS", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel offsets!\n");
                status = 0;     // Reset status
            } else {
                if (s->apply_offset < 0) {      // Use the data to decide
                    int jj;
                    if (ii == 0) {
                        s->dat_offs_col = colnum;
                    } else if (colnum != s->dat_offs_col) {
                        printf
                            ("Warning!:  DAT_OFFS column changes between files!\n");
                    }
                    float *fvec = (float *) malloc(sizeof(float) *
                                                   s->num_channels * s->num_polns);
                    fits_read_col(s->fitsfiles[ii], TFLOAT, s->dat_offs_col, 1L, 1L,
                                  s->num_channels * s->num_polns,
                                  0, fvec, &anynull, &status);
                    for (jj = 0; jj < s->num_channels * s->num_polns; jj++) {
                        // If the offsets are not 0, apply them
                        if (fvec[jj] != 0.0) {
                            s->apply_offset = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_offset < 0)
                    s->apply_offset = 0;        // not needed
            }

            // Data scalings
            fits_get_colnum(s->fitsfiles[ii], 0, "DAT_SCL", &colnum, &status);
            if (status == COL_NOT_FOUND) {
                printf("Warning!:  Can't find the channel scalings!\n");
                status = 0;     // Reset status
            } else {
                if (s->apply_scale < 0) {       // Use the data to decide
                    int jj;
                    if (ii == 0) {
                        s->dat_scl_col = colnum;
                    } else if (colnum != s->dat_scl_col) {
                        printf("Warning!:  DAT_SCL column changes between files!\n");
                    }
                    float *fvec = (float *) malloc(sizeof(float) *
                                                   s->num_channels * s->num_polns);
                    fits_read_col(s->fitsfiles[ii], TFLOAT, colnum, 1L, 1L,
                                  s->num_channels * s->num_polns,
                                  0, fvec, &anynull, &status);
                    for (jj = 0; jj < s->num_channels * s->num_polns; jj++) {
                        // If the scales are not 1, apply them
                        if (fvec[jj] != 1.0) {
                            s->apply_scale = 1;
                            break;
                        }
                    }
                    free(fvec);
                }
                if (s->apply_scale < 0)
                    s->apply_scale = 0; // not needed
            }
        }

        // Compute the samples per file and the amount of padding
        // that the _previous_ file has
        s->num_pad[ii] = 0;
        s->num_spec[ii] = s->spectra_per_subint * s->num_subint[ii];
        if (ii > 0) {
            if (s->start_spec[ii] > s->N) {     // Need padding
                s->num_pad[ii - 1] = s->start_spec[ii] - s->N;
                s->N += s->num_pad[ii - 1];
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
    if ((strncmp("AA+BB", s->poln_order, 5) == 0) ||
        (strncmp("INTEN", s->poln_order, 5) == 0))
        s->summed_polns = 1;
    else
        s->summed_polns = 0;

    // Is the data IQUV and the user poln is not set?
    if ((strncmp("IQUV", s->poln_order, 4) == 0) && (s->use_poln == 0))
        s->use_poln = 1;        // 1st poln = I

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
    if (s->bits_per_sample == 1 &&
        strcmp(s->telescope, "Parkes") == 0 &&
        strcmp(s->backend, "FB_1BIT") == 0) {
        printf("Flipping bit ordering since Parkes FB_1BIT data.\n");
        s->flip_bytes = 1;
    } else {
        s->flip_bytes = 0;
    }

    // Allocate the buffers
    cdatabuffer = gen_bvect(s->bytes_per_subint);
    // Following is twice as big because we use it as a ringbuffer too
    fdatabuffer = gen_fvect(2 * s->spectra_per_subint * s->num_channels);
    s->padvals = gen_fvect(s->num_channels);
    for (ii = 0; ii < s->num_channels; ii++)
        s->padvals[ii] = 0.0;
    offsets = gen_fvect(s->num_channels * s->num_polns);
    scales = gen_fvect(s->num_channels * s->num_polns);
    weights = gen_fvect(s->num_channels);
    // Initialize these if we won't be reading them from the file
    if (s->apply_offset == 0)
        for (ii = 0; ii < s->num_channels * s->num_polns; ii++)
            offsets[ii] = 0.0;
    if (s->apply_scale == 0)
        for (ii = 0; ii < s->num_channels * s->num_polns; ii++)
            scales[ii] = 1.0;
    if (s->apply_weight == 0)
        for (ii = 0; ii < s->num_channels; ii++)
            weights[ii] = 1.0;
}


long long offset_to_PSRFITS_spectra(long long specnum, struct spectra_info *s)
// This routine offsets into the PSRFITS files to the spectra
// 'specnum'.  It returns the current spectra number.
{
    int filenum = 0;

    if (specnum > s->N) {
        fprintf(stderr, "Error:  offset spectra %lld is > total spectra %lld\n\n",
                specnum, s->N);
        exit(1);
    }
    // Find which file we need
    while (filenum + 1 < s->num_files && specnum > s->start_spec[filenum + 1])
        filenum++;

    // Shift to that file
    cur_spec = specnum;
    new_spec = specnum;
    cur_file = filenum;
    numbuffered = 0;

    // Are we in a padding zone?
    if (specnum > (s->start_spec[cur_file] + s->num_spec[cur_file])) {
        // "Seek" to the end of the file
        cur_subint = s->num_subint[cur_file] + 1;
        new_spec = s->start_spec[cur_file + 1];
        return specnum;
    }
    // Otherwise, "seek" to the spectra (really a whole subint)
    // Check to make sure that specnum is the start of a subint
    if ((specnum - s->start_spec[cur_file]) % s->spectra_per_subint) {
        fprintf(stderr,
                "Error:  requested spectra %lld is not the start of a PSRFITS subint\n\n",
                specnum);
        exit(1);
    }
    // Remember zero offset for CFITSIO...
    cur_subint = (specnum - s->start_spec[cur_file]) / s->spectra_per_subint + 1;
    // printf("Current spectra = %lld, subint = %d, in file %d\n", cur_spec, cur_subint, cur_file);
    return specnum;
}


int get_PSRFITS_rawblock(float *fdata, struct spectra_info *s, int *padding)
// This routine reads a single block (i.e subint) from the input files
// which contain raw data in PSRFITS format.  If padding is
// returned as 1, then padding was added and statistics should not be
// calculated.  Return 1 on success.
{
    int numtopad = 0, numtoread, status = 0, anynull;
    float *fdataptr = fdata;

    fdataptr = fdata + numbuffered * s->num_channels;
    // numtoread is always this size since we need to read
    // full PSRFITS subints...
    numtoread = s->spectra_per_subint;

    // If our buffer array is offset from last time,
    // copy the previously offset part into the beginning.
    // New data comes after the old data in the buffer.
    if (numbuffered)
        memcpy((char *) fdata, (char *) (fdata + numtoread * s->num_channels),
               numbuffered * s->num_channels * sizeof(float));

    // Make sure our current file number is valid
    if (cur_file >= s->num_files)
        return 0;

    // Read a subint of data from the DATA col
    if (cur_subint <= s->num_subint[cur_file]) {
        double offs_sub = 0.0;
        if (!offs_sub_are_zero) {
            // Read the OFFS_SUB column value in case there were dropped blocks
            fits_read_col(s->fitsfiles[cur_file], TDOUBLE,
                          s->offs_sub_col, cur_subint, 1L, 1L,
                          0, &offs_sub, &anynull, &status);
            // Set new_spec to proper value, accounting for possibly
            // missing initial rows of data and/or combining observations
            // Note: need to remove start_subint because that was already put
            // into start_spec.  This is important if initial rows are gone.
            new_spec = s->start_spec[cur_file] +
                roundl((offs_sub - (s->start_subint[cur_file] + 0.5)
                        * s->time_per_subint) / s->dt);
        } else {
            new_spec = s->start_spec[cur_file] +
                (cur_subint - 1) * s->spectra_per_subint;
        }

        //printf("cur/new_spec = %lld, %lld  s->start_spec[cur_file] = %lld\n",
        //       cur_spec, new_spec, s->start_spec[cur_file]);

        // The following determines if there were lost blocks, or if
        // we are putting different observations together so that
        // the blocks are not aligned
        if (new_spec == cur_spec + numbuffered) {
            // if things look good, with no missing blocks, read the data
            get_PSRFITS_subint(fdataptr, cdatabuffer, s);
            cur_subint++;
            goto return_block;
        } else {
            goto padding_block;
        }
    } else {
        // We are going to move to the next file, so update
        // new_spec to be the starting spectra from the next file
        // so we can see if any padding is necessary
        if (cur_file < s->num_files - 1)
            new_spec = s->start_spec[cur_file + 1];
        else
            new_spec = cur_spec + numbuffered;
    }

    if (new_spec == cur_spec + numbuffered) {
        // No padding is necessary, so switch files
        cur_file++;
        cur_subint = 1;
        return get_PSRFITS_rawblock(fdata, s, padding);
    } else {                    // add padding
        goto padding_block;
    }

  padding_block:
    if (new_spec < cur_spec) {
        // Files out of order?  Shouldn't get here.
        fprintf(stderr, "Error!:  Current subint has earlier time than previous!\n\n"
                "\tfilename = '%s', subint = %d\n"
                "\tcur_spec = %lld  new_spec = %lld\n",
                s->filenames[cur_file], cur_subint, cur_spec, new_spec);
        exit(1);
    }
    numtopad = new_spec - cur_spec;
    // Don't add more than 1 block and if buffered, then realign the buffer
    if (numtopad > (s->spectra_per_subint - numbuffered))
        numtopad = s->spectra_per_subint - numbuffered;
    add_padding(fdataptr, s->padvals, s->num_channels, numtopad);
    // Update pointer into the buffer
    numbuffered = (numbuffered + numtopad) % s->spectra_per_subint;
    // Set the padding flag
    *padding = 1;
    // If we haven't gotten a full block, or completed the buffered one
    // then recursively call get_PSRFITS_rawblock()
    if (numbuffered) {
        printf("Adding %d spectra of padding to buffer at subint %d\n",
               numtopad, cur_subint);
        return get_PSRFITS_rawblock(fdata, s, padding);
    } else {
        printf("Adding %d spectra of padding at subint %d\n", numtopad, cur_subint);
        goto return_block;
    }

  return_block:
    // Apply the corrections that need a full block

    // Invert the band if needed
    if (s->apply_flipband)
        flip_band(fdata, s);

    // Perform Zero-DMing if requested
    if (s->remove_zerodm)
        remove_zerodm(fdata, s);

    // Increment our static counter (to determine how much data we
    // have written on the fly).
    cur_spec += s->spectra_per_subint;

    return 1;
}


void get_PSRFITS_subint(float *fdata, unsigned char *cdata, struct spectra_info *s)
{
    unsigned char *ctmp = cdata;
    int ii, status = 0, anynull;
    int numtoread = s->samples_per_subint;

    // The following allows us to read byte-packed data
    if (s->bits_per_sample < 8) {
        numtoread = (s->samples_per_subint * s->bits_per_sample) / 8;
        ctmp = gen_bvect(numtoread);
    }
    // or 16-bit data that is listed as being bytes
    if (s->bits_per_sample == 16 && s->FITS_typecode == 11)
        numtoread = s->samples_per_subint * 2;

    // Read the weights, offsets, and scales if required
    if (s->apply_weight)
        fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_wts_col, cur_subint, 1L,
                      s->num_channels, 0, weights, &anynull, &status);
    if (s->apply_offset)
        fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_offs_col, cur_subint,
                      1L, s->num_channels * s->num_polns, 0, offsets, &anynull,
                      &status);
    if (s->apply_scale)
        fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_scl_col, cur_subint, 1L,
                      s->num_channels * s->num_polns, 0, scales, &anynull, &status);

    // Now actually read the subint into the temporary buffer
    fits_read_col(s->fitsfiles[cur_file], s->FITS_typecode,
                  s->data_col, cur_subint, 1L, numtoread,
                  0, ctmp, &anynull, &status);

    if (status) {
        fprintf(stderr, "Error!:  Problem reading record from PSRFITS data file\n"
                "\tfilename = '%s', subint = %d.  FITS status = %d.  Exiting.\n",
                s->filenames[cur_file], cur_subint, status);
        exit(1);
    }
    // The following converts that byte-packed data into bytes
    if (s->bits_per_sample == 4) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 2 * ii;
            cdata[jj] = uctmp >> 4;
            cdata[jj + 1] = uctmp & 0x0F;
        }
    } else if (s->bits_per_sample == 2) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 4 * ii;
            cdata[jj] = ((uctmp >> 0x06) & 0x03);
            cdata[jj + 1] = ((uctmp >> 0x04) & 0x03);
            cdata[jj + 2] = ((uctmp >> 0x02) & 0x03);
            cdata[jj + 3] = (uctmp & 0x03);
        }
    } else if (s->bits_per_sample == 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 8 * ii;
            cdata[jj] = ((uctmp >> 0x07) & 0x01);
            cdata[jj + 1] = ((uctmp >> 0x06) & 0x01);
            cdata[jj + 2] = ((uctmp >> 0x05) & 0x01);
            cdata[jj + 3] = ((uctmp >> 0x04) & 0x01);
            cdata[jj + 4] = ((uctmp >> 0x03) & 0x01);
            cdata[jj + 5] = ((uctmp >> 0x02) & 0x01);
            cdata[jj + 6] = ((uctmp >> 0x01) & 0x01);
            cdata[jj + 7] = (uctmp & 0x01);
        }
    }

    if (s->bits_per_sample < 8)
        vect_free(ctmp);

    if (s->bits_per_sample == 1 && s->flip_bytes) {
        // Hack to flip each byte of data if needed
        for (ii = 0; ii < s->bytes_per_subint / 8; ii++) {
            int jj;
            const int offset = ii * 8;
            for (jj = 0; jj < 4; jj++) {
                unsigned char uctmp = cdata[offset + jj];
                cdata[offset + jj] = cdata[offset + 8 - 1 - jj];
                cdata[offset + 8 - 1 - jj] = uctmp;
            }
        }
    }
    // Now convert all of the data into floats

    // The following allows us to work with single polns out of many
    // or to sum polarizations if required
    if (s->num_polns > 1) {
        int sum_polns = 0;

        if ((0 == strncmp(s->poln_order, "AABB", 4)) || (s->num_polns == 2))
            sum_polns = 1;
        // User chose which poln to use
        if (s->use_poln > 0 || ((s->num_polns > 2) && !sum_polns)) {
            const int idx = (s->use_poln - 1) * s->num_channels;
            if (s->bits_per_sample == 16) {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const short *sptr =
                        (short *) cdata + ii * s->samples_per_spectra + idx;
                    for (jj = 0; jj < s->num_channels; jj++)
                        fptr[jj] =
                            (((float) sptr[jj] - s->zero_offset) * scales[idx + jj] +
                             offsets[idx + jj]) * weights[jj];
                }
            } else if (s->bits_per_sample == 32) {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const float *ftptr =
                        (float *) cdata + ii * s->samples_per_spectra + idx;
                    for (jj = 0; jj < s->num_channels; jj++)
                        fptr[jj] =
                            (((float) ftptr[jj] - s->zero_offset) * scales[idx +
                                                                           jj] +
                             offsets[idx + jj]) * weights[jj];
                }
            } else {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const unsigned char *cptr =
                        cdata + ii * s->samples_per_spectra + idx;
                    for (jj = 0; jj < s->num_channels; jj++)
                        fptr[jj] =
                            (((float) cptr[jj] - s->zero_offset) * scales[idx + jj] +
                             offsets[idx + jj]) * weights[jj];
                }
            }
        } else if (sum_polns) { // sum the polns if there are 2 by default
            const int idx = s->num_channels;
            if (s->bits_per_sample == 16) {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const short *sptr =
                        (short *) cdata + ii * s->samples_per_spectra;
                    for (jj = 0; jj < s->num_channels; jj++) {
                        fptr[jj] =
                            (((float) sptr[jj] - s->zero_offset) * scales[jj] +
                             offsets[jj]) * weights[jj];
                        fptr[jj] +=
                            (((float) sptr[jj + idx] - s->zero_offset) * scales[idx +
                                                                                jj] +
                             offsets[idx + jj]) * weights[jj];
                    }
                }
            } else if (s->bits_per_sample == 32) {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const float *ftptr =
                        (float *) cdata + ii * s->samples_per_spectra;
                    for (jj = 0; jj < s->num_channels; jj++) {
                        fptr[jj] =
                            (((float) ftptr[jj] - s->zero_offset) * scales[jj] +
                             offsets[jj]) * weights[jj];
                        fptr[jj] +=
                            (((float) ftptr[jj + idx] -
                              s->zero_offset) * scales[idx + jj] + offsets[idx +
                                                                           jj]) *
                            weights[jj];
                    }
                }
            } else {
#ifdef _OPENMP
#pragma omp parallel for shared(s,cdata,fdata,scales,offsets,weights)
#endif
                for (ii = 0; ii < s->spectra_per_subint; ii++) {
                    int jj;
                    float *fptr = fdata + ii * s->num_channels;
                    const unsigned char *cptr = cdata + ii * s->samples_per_spectra;
                    for (jj = 0; jj < s->num_channels; jj++) {
                        fptr[jj] =
                            (((float) cptr[jj] - s->zero_offset) * scales[jj] +
                             offsets[jj]) * weights[jj];
                        fptr[jj] +=
                            (((float) cptr[jj + idx] - s->zero_offset) * scales[idx +
                                                                                jj] +
                             offsets[idx + jj]) * weights[jj];
                    }
                }
            }
        }
    } else {                    // This is for normal single-polarization data
        if (s->bits_per_sample == 16) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < s->spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * s->num_channels;
                const short *sptr = (short *) cdata + ii * s->samples_per_spectra;
                for (jj = 0; jj < s->num_channels; jj++)
                    fptr[jj] = (((float) sptr[jj] - s->zero_offset) * scales[jj] +
                                offsets[jj]) * weights[jj];
            }
        } else if (s->bits_per_sample == 32) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < s->spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * s->num_channels;
                const float *ftptr = (float *) cdata + ii * s->samples_per_spectra;
                for (jj = 0; jj < s->num_channels; jj++)
                    fptr[jj] = (((float) ftptr[jj] - s->zero_offset) * scales[jj] +
                                offsets[jj]) * weights[jj];
            }
        } else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < s->spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * s->num_channels;
                const unsigned char *cptr = cdata + ii * s->samples_per_spectra;
                for (jj = 0; jj < s->num_channels; jj++)
                    fptr[jj] = (((float) cptr[jj] - s->zero_offset) * scales[jj] +
                                offsets[jj]) * weights[jj];
            }
        }
    }
}
