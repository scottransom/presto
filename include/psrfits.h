#include "backend_common.h"

/* psrfits.c  (generated automatically by cproto) */
void get_PSRFITS_static(int *bytesperpt, int *bytesperblk, int *numifs, float *clip_sigma);
void set_PSRFITS_static(int ptsperblk, int bytesperpt, int bytesperblk, int numchan, int numifs, float clip_sigma, double dt);
void set_PSRFITS_padvals(float *fpadvals, int good_padvals);
int is_PSRFITS(char *filename);
int read_PSRFITS_files(char **filenames, int numfiles, struct spectra_info *s);
void print_PSRFITS_info(struct spectra_info *s);
void spectra_info_to_inf(struct spectra_info *s, infodata *idata);
void get_PSRFITS_file_info(char **filenames, int numfiles, float clipsig, struct spectra_info *s, infodata *idata, int output);
void PSRFITS_update_infodata(infodata *idata);
int skip_to_PSRFITS_samp(long long samplenum);
int read_PSRFITS_rawblock(unsigned char *data, int *padding);
int read_PSRFITS_rawblocks(unsigned char rawdata[], int numblocks, int *padding);
int read_PSRFITS(float *data, int numspec, double *dispdelays, int *padding, int *maskchans, int *nummasked, mask *obsmask);
void get_PSRFITS_channel(int channum, float chandat[], unsigned char rawdata[], int numblocks);
int prep_PSRFITS_subbands(unsigned char *rawdata, float *data, double *dispdelays, int numsubbands, int transpose, int *maskchans, int *nummasked, mask *obsmask);
int read_PSRFITS_subbands(float *data, double *dispdelays, int numsubbands, int transpose, int *padding, int *maskchans, int *nummasked, mask *obsmask);



