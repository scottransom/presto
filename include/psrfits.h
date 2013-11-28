#include "backend_common.h"

/* psrfits.c  (generated automatically by cproto) */
int is_PSRFITS(char *filename);
void read_PSRFITS_files(struct spectra_info *s);
long long offset_to_PSRFITS_spectra(long long specnum, struct spectra_info *s);
int get_PSRFITS_rawblock(float *fdata, struct spectra_info *s, int *padding);
