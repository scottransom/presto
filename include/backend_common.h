#include "fitsio.h"
#include "mask.h"
#include "makeinf.h"


typedef enum {
    SIGPROCFB, PSRFITS, SCAMP, BPP, WAPP, SPIGOT, \
    SUBBAND, DAT, SDAT, EVENTS, UNSET
} psrdatatype;


struct spectra_info {
    char telescope[40];     // Telescope used
    char observer[100];     // Observer's name
    char source[100];       // Source name
    char frontend[100];     // Frontend used
    char backend[100];      // Backend or instrument used
    char project_id[40];    // Project identifier
    char date_obs[40];      // Start of observation (YYYY-MM-DDTHH:MM:SS.SSS)
    char ra_str[40];        // Right Ascension string (HH:MM:SS.SSSS)
    char dec_str[40];       // Declination string (DD:MM:SS.SSSS)
    char poln_type[40];     // Polarization recorded (LIN or CIRC)
    char poln_order[40];    // Order of polarizations (i.e. XXYYXYYX)
    long long N;            // Total number of spectra in observation
    double T;               // Total observation time (N * dt)
    double dt;              // Sample duration (s)
    double fctr;            // Center frequency of the observing band (MHz)
    double lo_freq;         // Center freq of the lowest channel (MHz)
    double hi_freq;         // Center freq of the highest channel (MHz)
    double orig_df;         // Original frequency spacing between the channels (MHz)
    double chan_dm;         // DM used to de-disperse the freq channels
    double df;              // Frequency spacing between the channels (MHz)
    double BW;              // Bandwidth of the observing band (MHz)
    double ra2000;          // RA  of observation (deg, J2000) at obs start
    double dec2000;         // Dec of observation (deg, J2000) at obs start
    double azimuth;         // Azimuth (commanded) at the start of the obs (deg)
    double zenith_ang;      // Zenith angle (commanded) at the start of the obs (deg)
    double beam_FWHM;       // Beam FWHM (deg)
    double time_per_subint; // Duration (in sec) of a full SUBINT entry
    psrdatatype datatype;   // The data format being used
    int scan_number;        // Number of scan
    int tracking;           // Tracking (1) or drift scan (0)
    int orig_num_chan;      // Number of spectral channels per sample
    int num_channels;       // Number of spectral channels per sample
    int num_polns;          // 1, 2, or 4 (for Full Stokes)
    int num_beams;          // Number of beams in the data
    int beamnum;            // Beam number if num_beams=1 or beam to process if num_beams > 1 (0 offset)
    int summed_polns;       // Are polarizations summed? (1=Yes, 0=No)
    int FITS_typecode;      // FITS data typecode as per CFITSIO
    int bits_per_sample;    // Number of bits per sample
    int bytes_per_spectra;  // Number of bytes in a spectra (inluding all polns)
    int samples_per_spectra;// Number of samples in a spectra (inluding all polns)
    int bytes_per_subint;   // Number of bytes per SUBINT entry
    int spectra_per_subint; // Number of spectra per SUBINT entry
    int samples_per_subint; // Number of samples per SUBINT entry
    int min_spect_per_read; // The minimum number of spectra you can read at once
    int num_files;          // Number of files in the observation
    int offs_sub_col;       // The number of the OFFS_SUB column in the SUBINT HDU
    int dat_wts_col;        // The number of the DAT_WTS column in the SUBINT HDU
    int dat_offs_col;       // The number of the DAT_OFFS column in the SUBINT HDU
    int dat_scl_col;        // The number of the DAT_SCL column in the SUBINT HDU
    int data_col;           // The number of the DATA column in the SUBINT HDU
    int apply_scale;        // Do we apply the scales to the data? (1=Yes, 0=No)
    int apply_offset;       // Do we apply the offsets to the data? (1=Yes, 0=No)
    int apply_weight;       // Do we apply the weights to the data? (1=Yes, 0=No)
    int apply_flipband;     // Do we invert the band?
    int signedints;         // Used signed bytes rather than default unsigned bytes
    int remove_zerodm;      // Do zero-DM substraction?
    int use_poln;           // The number of the specific polarization to use 0-num_polns-1
    int flip_bytes;         // Hack to flip the order of the bits in a byte of data
    int num_ignorechans;    // Number of channels to explicitly ignore (set to zero)
    float zero_offset;      // A DC zero-offset value to apply to all the data
    float clip_sigma;       // Clipping value in standard deviations to use
    long double *start_MJD; // Array of long double MJDs for the file starts
    char **filenames;       // Array of the input file names
    FILE **files;           // Array of normal file pointers if needed
    fitsfile **fitsfiles;    // Array of FITS file pointers if needed
    float *padvals;         // Array of length num_channels for current padding values
    int *header_offset;     // Number of bytes to skip in each file to get to spectra
    int *start_subint;      // Starting ISUBINT in each file (for PSRFITS)
    int *num_subint;        // Number of subints per file
    int *ignorechans;       // Array of channels to explicitly ignore (set to zero)
    char *ignorechans_str; // String describing the ignorechans
    long long *start_spec; // Starting spectra for each file (zero offset)
    long long *num_spec;   // Number of spectra per file
    long long *num_pad;    // Number of padding samples after each file
    int (*get_rawblock)(float *, struct spectra_info *, int *);  // Raw data block function pointer
    long long (*offset_to_spectra)(long long, struct spectra_info *);  // Shift into file(s) function pointer
};


/* backend_common.c */
void psrdatatype_description(char *outstr, psrdatatype ptype);
void identify_psrdatatype(struct spectra_info *s, int output);
void close_rawfiles(struct spectra_info *s);
void spectra_info_set_defaults(struct spectra_info *s);
void read_rawdata_files(struct spectra_info *s);
void print_spectra_info(struct spectra_info *s);
void print_spectra_info_summary(struct spectra_info *s);
void spectra_info_to_inf(struct spectra_info *s, infodata *idata);
long long offset_to_spectra(long long specnum, struct spectra_info *s);
void set_currentspectra(long long specnum);
int read_rawblocks(float *fdata, int numsubints, struct spectra_info *s, int *padding);
int read_psrdata(float *fdata, int numspect, struct spectra_info *s, int *delays, int *padding, int *maskchans, int *nummasked, mask *obsmask);
void get_channel(float chandat[], int channum, int numsubints, float rawdata[], struct spectra_info *s);
int prep_subbands(float *fdata, float *rawdata, int *delays, int numsubbands, struct spectra_info *s, int transpose, int *maskchans, int *nummasked, mask *obsmask);
int read_subbands(float *fdata, int *delays, int numsubbands, struct spectra_info *s, int transpose, int *padding, int *maskchans, int *nummasked, mask *obsmask);
void flip_band(float *fdata, struct spectra_info *s);
int *get_ignorechans(char *ignorechans_str, int minchan, int maxchan, int *num_ignorechans, char **filestr);

/* zerodm.c */
void remove_zerodm(float *fdata, struct spectra_info *s);
