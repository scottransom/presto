#include "backend_common.h"
#include "sigproc_fb.h"

void print_psrdatatype(psrdatatype type)
    switch(psrdatatype) {
    case SIGPROCFB:
        printf("Data are SIGPROC filterbank format.\n");
        break;
    case PSRFITS:
        printf("Data are PSRFITS searchmode format.\n");
        break;
    case SCAMP:
        printf("Data are 'SCAMP' 1-bit format (un_sc_td'd).\n");
        break;
    case BPP:
        printf("Data are BCPM (i.e. BPP) format.\n");
        break;
    case GMRT:
        printf("Data are old-type GMRT format.\n");
        break;
    case SPIGOT:
        printf("Data are GBT SPIGOT format.\n");
        break;
    case SUBBAND:
        printf("Data are PRESTO subband format.\n");
        break;
    default:
        printf("Data are in an unknown format!!!\n");
    }
}


void spectra_info_set_defaults(struct spectra_info *s) {
    strcpy(s->telescope, "unset");
    strcpy(s->observer, "unset");
    strcpy(s->source, "unset");
    strcpy(s->frontend, "unset");
    strcpy(s->backend, "unset");
    strcpy(s->project_id, "unset");
    strcpy(s->date_obs, "unset");
    strcpy(s->ra_str, "unset");
    strcpy(s->dec_str, "unset");
    strcpy(s->poln_type, "unset");
    strcpy(s->poln_order, "unset");
    s->datatype = UNSET;
    s->N = 0;
    s->T = 0.0;
    s->dt = 0.0;
    s->fctr = 0.0;
    s->lo_freq = 0.0;
    s->hi_freq = 0.0;
    s->orig_df = 0.0;
    s->chan_dm = 0.0;
    s->df = 0.0;
    s->BW = 0.0;
    s->ra2000 = 0.0;
    s->dec2000 = 0.0;
    s->azimuth = 0.0;
    s->zenith_ang = 0.0;
    s->beam_FWHM = 0.0;
    s->time_per_subint = 0.0;
    s->scan_number = 0;
    s->tracking = 1;
    s->orig_num_chan = 0;
    s->num_channels = 0;
    s->num_polns = 0;
    s->num_beams = 1;
    s->beamnum = 0;
    s->summed_polns = 1;
    s->FITS_typecode = 0;
    s->bits_per_sample = 0;
    s->bytes_per_spectra = 0;
    s->samples_per_spectra = 0;
    s->bytes_per_subint = 0;
    s->spectra_per_subint = 0;
    s->samples_per_subint = 0;
    s->min_spect_per_read = 0;
    s->num_files = 0;
    s->offs_sub_col = 0;
    s->dat_wts_col = 0;
    s->dat_offs_col = 0;
    s->dat_scl_col = 0;
    s->data_col = 0;
    s->apply_scale = 0;
    s->apply_offset = 0;
    s->apply_weight = 0;
    s->apply_flipband = 0;
    s->remote_zeroDM = 0;
    s->use_poln = 0;
    s->flip_bytes = 0;
    s->zero_offset = 0.0;
    s->clip_sigma = 0.0;
    s->start_MJD = NULL;
    s->files = NULL;
    s->fitsfiles = NULL;
    s->header_offset = NULL;
    s->start_subint = NULL;
    s->num_subint = NULL;
    s->start_spec = NULL;
    s->num_spec = NULL;
    s->num_pad = NULL;
};


int read_filterbank_rawblocks(FILE * infiles[], int numfiles,
                              float rawdata[], int numblocks, int *padding)
// This routine reads numblocks filterbank records from the input
// files *infiles.  The floating-point filterbank data is returned in
// rawdata which must have a size of numblocks * s->samples_per_subint.  The
// number of blocks read is returned.  If padding is returned as 1,
// then padding was added and statistics should not be calculated
{
    int ii, retval = 0, pad, numpad = 0;

    *padding = 0;
    for (ii = 0; ii < numblocks; ii++) {
        retval += read_filterbank_rawblock(infiles, numfiles,
                                           rawdata + ii * s->samples_per_subint, &pad);
        if (pad) numpad++;
    }
    /* Return padding 'true' if any block was padding */
    if (numpad) *padding = 1;
    return retval;
}


int read_filterbank(FILE * infiles[], int numfiles, float *data,
                    int numspect, int *delays, int *padding,
                    int *maskchans, int *nummasked, mask * obsmask)
// This routine reads numspect from the filterbank raw files *infiles.
// These files contain raw data in filterbank format.  Time delays and
// a mask are applied to each channel.  It returns the # of points
// read if successful, 0 otherwise.  If padding is returned as 1, then
// padding was added and statistics should not be calculated.
// maskchans is an array of length numchans contains a list of the
// number of channels that were masked.  The # of channels masked is
// returned in nummasked.  obsmask is the mask structure to use for
// masking.
{
   int ii, jj, numread = 0, offset;
   double starttime = 0.0;
   static float *tempzz, *rawdata1, *rawdata2;
   static float *currentdata, *lastdata;
   static int firsttime = 1, numblocks = 1, allocd = 0, mask = 0;
   static double duration = 0.0;

   *nummasked = 0;
   if (firsttime) {
      if (numspect % s->spectra_per_subint) {
         printf("numspect must be a multiple of %d in read_filterbank()!\n",
                s->spectra_per_subint);
         exit(1);
      } else
         numblocks = numspect / s->spectra_per_subint;
      if (obsmask->numchan)
         mask = 1;
      rawdata1 = gen_fvect(numblocks * s->samples_per_subint);
      rawdata2 = gen_fvect(numblocks * s->samples_per_subint);
      allocd = 1;
      duration = numblocks * s->time_per_subint;
      currentdata = rawdata1;
      lastdata = rawdata2;
   }

   /* Read, convert and de-disperse */
   if (allocd) {
      while (1) {
         numread = read_filterbank_rawblocks(infiles, numfiles, currentdata,
                                             numblocks, padding);
         if (mask) {
            starttime = currentblock * s->time_per_subint;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
         }

         /* Clip nasty RFI if requested and we're not masking all the channels */
         if ((clip_sigma_st > 0.0) && !(mask && (*nummasked == -1)))
            clip_times(currentdata, numspect, s->num_channels, clip_sigma_st, padvals);

         if (mask) {
            if (*nummasked == -1) {     /* If all channels are masked */
               for (ii = 0; ii < numspect; ii++)
                   memcpy(currentdata + ii * s->num_channels, 
                          padvals, s->num_channels * sizeof(float));
            } else if (*nummasked > 0) {        /* Only some of the channels are masked */
               int channum;
               for (ii = 0; ii < numspect; ii++) {
                  offset = ii * s->num_channels;
                  for (jj = 0; jj < *nummasked; jj++) {
                     channum = maskchans[jj];
                     currentdata[offset + channum] = padvals[channum];
                  }
               }
            }
         }
         if (!firsttime)
             float_dedisp(currentdata, lastdata, numspect, s->num_channels, delays, 0.0, data);

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
      return numblocks * s->spectra_per_subint;
   } else {
      return 0;
   }
}


void get_filterbank_channel(int channum, float chandat[],
                            float rawdata[], int numblocks)
// Return the values for channel 'channum' of a block of 'numblocks'
// filterbank floating-point data stored in 'rawdata' in 'chandat'.
// 'rawdata' should have been initialized using
// read_filterbank_rawblocks(), and 'chandat' must have at least
// 'numblocks' * 's->spectra_per_subint' spaces.  Channel 0 is assumed to be
// the lowest freq channel.
{
   int ii, jj;

   if (channum > s->num_channels || channum < 0) {
      printf("\nchannum = %d is out of range in get_GMR_channel()!\n\n", channum);
      exit(1);
   }
   /* Select the correct channel */
   for (ii = 0, jj = channum; ii < numblocks * s->spectra_per_subint; ii++, jj += s->num_channels)
      chandat[ii] = rawdata[jj];
}


int prep_filterbank_subbands(float *rawdata, float *data,
                             int *delays, int numsubbands,
                             int transpose, int *maskchans,
                             int *nummasked, mask * obsmask)
// This routine preps a block from from the filterbank system.  It
// uses dispersion delays in 'delays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'data' of length 'numsubbands' * 's->spectra_per_subint'.  The low freq
// subband is stored first, then the next highest subband etc, with
// 's->spectra_per_subint' floating points per subband. It returns the # of
// points read if succesful, 0 otherwise.  'maskchans' is an array of
// length numchans which contains a list of the number of channels
// that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
   int ii, jj, trtn, offset;
   double starttime = 0.0;
   static float *tempzz;
   static float rawdata1[MAXNUMCHAN * BLOCKLEN],
       rawdata2[MAXNUMCHAN * BLOCKLEN];
   static float *currentdata, *lastdata;
   static unsigned char *move;
   static int firsttime = 1, move_size = 0, mask = 0;

   *nummasked = 0;
   if (firsttime) {
      if (obsmask->numchan)
         mask = 1;
      move_size = (s->spectra_per_subint + numsubbands) / 2;
      move = gen_bvect(move_size);
      currentdata = rawdata1;
      lastdata = rawdata2;
   }

   /* Read and de-disperse */
   memcpy(currentdata, rawdata, s->samples_per_subint * sizeof(float));
   if (mask) {
      starttime = currentblock * s->time_per_subint;
      *nummasked = check_mask(starttime, s->time_per_subint, 
                              obsmask, maskchans);
   }

   /* Clip nasty RFI if requested and we're not masking all the channels*/
   if ((clip_sigma_st > 0.0) && !(mask && (*nummasked == -1)))
      clip_times(currentdata, s->spectra_per_subint, s->num_channels, clip_sigma_st, padvals);

   if (mask) {
      if (*nummasked == -1) {   /* If all channels are masked */
         for (ii = 0; ii < s->spectra_per_subint; ii++)
             memcpy(currentdata + ii * s->num_channels, 
                    padvals, s->num_channels * sizeof(float));
      } else if (*nummasked > 0) {      /* Only some of the channels are masked */
         int channum;
         for (ii = 0; ii < s->spectra_per_subint; ii++) {
            offset = ii * s->num_channels;
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
      dedisp_subbands(currentdata, lastdata, s->spectra_per_subint, s->num_channels,
                      delays, numsubbands, data);
      SWAP(currentdata, lastdata);
      /* Transpose the data into vectors in the result array */
      if (transpose) {
         if ((trtn = transpose_float(data, s->spectra_per_subint, numsubbands,
                                     move, move_size)) < 0)
            printf("Error %d in transpose_float().\n", trtn);
      }
      return s->spectra_per_subint;
   }
}


int read_filterbank_subbands(FILE *infiles[], int numfiles, float *data,
                             int *delays, int numsubbands,
                             int transpose, int *padding,
                             int *maskchans, int *nummasked, mask * obsmask)
// This routine reads a record from the input files *infiles[] in
// SIGPROC filterbank format.  The routine uses dispersion delays in
// 'delays' to de-disperse the data into 'numsubbands' subbands.
// It stores the resulting data in vector 'data' of length
// 'numsubbands' * 's->spectra_per_subint'.  The low freq subband is stored
// first, then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// succesful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
   static int firsttime = 1;
   static float rawdata[MAXNUMCHAN * BLOCKLEN];

   if (firsttime) {
      if (!read_filterbank_rawblock(infiles, numfiles, rawdata, padding)) {
         printf("Problem reading the raw filterbank data file.\n\n");
         return 0;
      }
      if (0 != prep_filterbank_subbands(rawdata, data, delays, numsubbands,
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
   return prep_filterbank_subbands(rawdata, data, delays, numsubbands,
                                   transpose, maskchans, nummasked, obsmask);
}
