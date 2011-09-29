#include "presto.h"
#include "mask.h"
#include "gmrt.h"

#define MAXNUMCHAN 1024
#define BLOCKLEN   1024

/* All of the following have an _st to indicate static */
static long long numpts_st[MAXPATCHFILES], padpts_st[MAXPATCHFILES], N_st;
static int numblks_st[MAXPATCHFILES], need_byteswap_st = 0, sampperblk_st;
static int numchan_st, ptsperblk_st, bytesperpt_st, bytesperblk_st;
static double times_st[MAXPATCHFILES], mjds_st[MAXPATCHFILES];
static double elapsed_st[MAXPATCHFILES], T_st, dt_st;
static double startblk_st[MAXPATCHFILES], endblk_st[MAXPATCHFILES];
static infodata idata_st[MAXPATCHFILES];
static unsigned char padvals[MAXNUMCHAN], newpadvals[MAXNUMCHAN], padval = 128;
static short sdatabuffer[MAXNUMCHAN * BLOCKLEN];
static unsigned char databuffer[MAXNUMCHAN * BLOCKLEN];
static int currentfile, currentblock;
static int bufferpts = 0, padnum = 0, shiftbuffer = 1;
static float clip_sigma_st = 0.0;
static unsigned char cdatabuffer[MAXNUMCHAN * BLOCKLEN];
static int sb_flag;
static int bitshift = 2;
static int bits;
static char badch[100];
static int using_MPI = 0;


void get_GMRT_static(int *bytesperpt, int *bytesperblk, float *clip_sigma)
{
   *bytesperpt = bytesperpt_st;
   *bytesperblk = bytesperblk_st;
   *clip_sigma = clip_sigma_st;
}


void set_GMRT_static(int ptsperblk, int bytesperpt, int bytesperblk,
                     int numchan, float clip_sigma, double dt)
{
   using_MPI = 1;
   currentblock = 0;
   ptsperblk_st = ptsperblk;
   bytesperpt_st = bytesperpt;
   bytesperblk_st = bytesperblk;
   numchan_st = numchan;
   sampperblk_st = ptsperblk_st * numchan_st;
   clip_sigma_st = clip_sigma;
   dt_st = dt;
}


void set_GMRT_padvals(float *fpadvals, int good_padvals)
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


void get_GMRT_file_info(FILE * files[], char *datfilenms[], int numfiles,
                        float clipsig, long long *N, int *ptsperblock, int *numchan,
                        double *dt, double *T, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of GMRT rawfiles that you want to patch     */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
   int ii;
  /** GMRT_hdr header; */

   if (numfiles > MAXPATCHFILES) {
      printf("\nThe number of input files (%d) is greater than \n", numfiles);
      printf("   MAXPATCHFILES=%d.  Exiting.\n\n", MAXPATCHFILES);
      exit(0);
   }
   GMRT_hdr_to_inf(datfilenms[0], &idata_st[0]);
   /* Are we going to clip the data? */
   if (clipsig > 0.0)
      clip_sigma_st = clipsig;
   *numchan = numchan_st = idata_st[0].num_chan;
   *ptsperblock = ptsperblk_st = BLOCKLEN;
   sampperblk_st = ptsperblk_st * numchan_st;
   bytesperblk_st = bytesperpt_st * sampperblk_st;
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
   for (ii = 1; ii < numfiles; ii++) {
      GMRT_hdr_to_inf(datfilenms[ii], &idata_st[ii]);
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
              ii + 1, startblk_st[ii], endblk_st[ii], numpts_st[ii], elapsed_st[ii],
              times_st[ii], mjds_st[ii], padpts_st[ii]);
      printf("\n");
   }
}

void GMRT_update_infodata(int numfiles, infodata * idata)
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


int skip_to_GMRT_rec(FILE * infiles[], int numfiles, int rec)
/* This routine skips to the record 'rec' in the input files   */
/* *infiles.  *infiles contains 1 bit digitized data from the  */
/* GMRT backend at Parkes.  Returns the record skipped to.     */
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
                     bytesperblk_st, SEEK_SET);
         bufferpts = (int) ((startblk_st[filenum] - floor_blk) * ptsperblk_st + 0.5);
         padnum = 0;
         /*
            printf("Data:  currentfile = %d  bufferpts = %d  padnum = %d\n", 
            currentfile, bufferpts, padnum);
          */
      }

   } else {
      printf("\n rec = %d out of range in skip_to_GMRT_rec()\n", rec);
      exit(1);
   }
   return rec;
}


int read_GMRT_rawblock(FILE * infiles[], int numfiles,
                       unsigned char *data, int *padding)
/* This routine reads a single record from the         */
/* input files *infiles which contain 2 byte digitized */
/* data from the GMRT pulsar backend.                  */
/* If padding is returned as 1, then padding was       */
/* added and statistics should not be calculated       */
{
   int offset, numtopad = 0, ii;
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
   if (bits == 16) {
      if (chkfread(sdatabuffer, bytesperblk_st, 1, infiles[currentfile])) {     /* Got Data */
         /* See if we need to byte-swap and if so, doit */
         if (need_byteswap_st) {
            short *sptr = sdatabuffer;
            for (ii = 0; ii < sampperblk_st; ii++, sptr++)
               *sptr = swap_short(*sptr);
         }
         convert_GMRT_block(sdatabuffer, dataptr);

         /* Clip nasty RFI if requested */
         if (clip_sigma_st > 0.0)
            clip_times(dataptr, ptsperblk_st, numchan_st, clip_sigma_st, newpadvals);
         *padding = 0;

         /* Put the new data into the databuffer if needed */
         if (bufferpts)
            memcpy(data, databuffer, sampperblk_st);
         currentblock++;
         return 1;
      } else {                  /* Didn't get data */
         if (feof(infiles[currentfile])) {      /* End of file? */
            numtopad = padpts_st[currentfile] - padnum;
            if (numtopad) {     /* Pad the data? */
               *padding = 1;
               if (numtopad >= ptsperblk_st - bufferpts) {      /* Lots of padding */
                  if (bufferpts) {      /* Buffer the padding? */
                     /* Add the amount of padding we need to */
                     /* make our buffer offset = 0           */
                     numtopad = ptsperblk_st - bufferpts;
                     for (ii = 0; ii < numtopad; ii++)
                        memcpy(dataptr + ii * numchan_st, newpadvals, numchan_st);
                     /* Copy the new data/padding into the output array */
                     memcpy(data, databuffer, sampperblk_st);
                     bufferpts = 0;
                  } else {      /* Add a full record of padding */
                     numtopad = ptsperblk_st;
                     for (ii = 0; ii < numtopad; ii++)
                        memcpy(data + ii * numchan_st, newpadvals, numchan_st);
                  }
                  padnum += numtopad;
                  currentblock++;
                  /* If done with padding reset padding variables */
                  if (padnum == padpts_st[currentfile]) {
                     padnum = 0;
                     currentfile++;
                  }
                  return 1;
               } else {         /* Need < 1 block (or remaining block) of padding */
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
                  return read_GMRT_rawblock(infiles, numfiles, data, &pad);
               }
            } else {            /* No padding needed.  Try reading the next file */
               currentfile++;
               shiftbuffer = 0;
               return read_GMRT_rawblock(infiles, numfiles, data, padding);
            }
         } else {
            printf("\nProblem reading record from GMRT data file:\n");
            printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                   currentfile, currentblock);
            exit(1);
         }
      }
   } else {                     /* here if handling 8 bit data files */

      /* First, attempt to read data from the current file */

      if (chkfread(cdatabuffer, bytesperblk_st, 1, infiles[currentfile])) {     /* Got Data */
         /* See if we need to byte-swap and if so, doit */
/*
      if (need_byteswap_st){
          short *sptr = cdatabuffer;
          for (ii=0; ii<sampperblk_st; ii++, sptr++)
  	  *sptr = swap_short(*sptr);
      }
*/
         convert_GMRT_block_8bit(cdatabuffer, dataptr);

         /* Clip nasty RFI if requested */
         if (clip_sigma_st > 0.0)
            clip_times(dataptr, ptsperblk_st, numchan_st, clip_sigma_st, newpadvals);
         *padding = 0;

         /* Put the new data into the databuffer if needed */
         if (bufferpts)
            memcpy(data, databuffer, sampperblk_st);
         currentblock++;
         return 1;
      } else {                  /* Didn't get data */
         if (feof(infiles[currentfile])) {      /* End of file? */
            numtopad = padpts_st[currentfile] - padnum;
            if (numtopad) {     /* Pad the data? */
               *padding = 1;
               if (numtopad >= ptsperblk_st - bufferpts) {      /* Lots of padding */
                  if (bufferpts) {      /* Buffer the padding? */
                     /* Add the amount of padding we need to */
                     /* make our buffer offset = 0           */
                     numtopad = ptsperblk_st - bufferpts;
                     for (ii = 0; ii < numtopad; ii++)
                        memcpy(dataptr + ii * numchan_st, newpadvals, numchan_st);
                     /* Copy the new data/padding into the output array */
                     memcpy(data, databuffer, sampperblk_st);
                     bufferpts = 0;
                  } else {      /* Add a full record of padding */
                     numtopad = ptsperblk_st;
                     for (ii = 0; ii < numtopad; ii++)
                        memcpy(data + ii * numchan_st, newpadvals, numchan_st);
                  }
                  padnum += numtopad;
                  currentblock++;
                  /* If done with padding reset padding variables */
                  if (padnum == padpts_st[currentfile]) {
                     padnum = 0;
                     currentfile++;
                  }
                  return 1;
               } else {         /* Need < 1 block (or remaining block) of padding */
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
                  return read_GMRT_rawblock(infiles, numfiles, data, &pad);
               }
            } else {            /* No padding needed.  Try reading the next file */
               currentfile++;
               shiftbuffer = 0;
               return read_GMRT_rawblock(infiles, numfiles, data, padding);
            }
         } else {
            printf("\nProblem reading record from GMRT data file:\n");
            printf("   currentfile = %d, currentblock = %d.  Exiting.\n",
                   currentfile, currentblock);
            exit(1);
         }
      }
   }
}


int read_GMRT_rawblocks(FILE * infiles[], int numfiles,
                        unsigned char rawdata[], int numblocks, int *padding)
/* This routine reads numblocks GMRT records from the input */
/* files *infiles.  The 8-bit filterbank data is returned   */
/* in rawdata which must have a size of numblocks*          */
/* sampperblk_st.  The number  of blocks read is returned.  */
/* If padding is returned as 1, then padding was added      */
/* and statistics should not be calculated                  */
{
   int ii, retval = 0, pad, numpad = 0;

   *padding = 0;
   for (ii = 0; ii < numblocks; ii++) {
      retval += read_GMRT_rawblock(infiles, numfiles,
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


int read_GMRT(FILE * infiles[], int numfiles, float *data,
              int numpts, double *dispdelays, int *padding,
              int *maskchans, int *nummasked, mask * obsmask)
/* This routine reads numpts from the GMRT raw input   */
/* files *infiles.  These files contain 16 bit data    */
/* from the GMRT backend.  Time delays and             */
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
         printf("numpts must be a multiple of %d in read_GMRT()!\n", ptsperblk_st);
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
         numread = read_GMRT_rawblocks(infiles, numfiles, currentdata,
                                       numblocks, padding);
         if (mask) {
            starttime = currentblock * timeperblk;
            *nummasked = check_mask(starttime, duration, obsmask, maskchans);
         }

         /* Only use the recently measured padding if all the channels aren't masked */
         if ((clip_sigma_st > 0.0) && 
             !(mask && (*nummasked == -1)) &&
             (padvals != newpadvals))
             memcpy(padvals, newpadvals, MAXNUMCHAN);

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


void get_GMRT_channel(int channum, float chandat[],
                      unsigned char rawdata[], int numblocks)
/* Return the values for channel 'channum' of a block of       */
/* 'numblocks' raw GMRT data stored in 'rawdata' in 'chandat'. */
/* 'rawdata' should have been initialized using                */
/* read_GMRT_rawblocks(), and 'chandat' must have at least     */
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


int prep_GMRT_subbands(unsigned char *rawdata, float *data,
                       double *dispdelays, int numsubbands,
                       int transpose, int *maskchans, int *nummasked, mask * obsmask)
/* This routine preps a block from from the GMRT system.  It uses         */
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

   /* Only use the recently measured padding if all the channels aren't masked */
   if ((clip_sigma_st > 0.0) && 
       !(mask && (*nummasked == -1)) &&
       (padvals != newpadvals))
       memcpy(padvals, newpadvals, MAXNUMCHAN);

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


int read_GMRT_subbands(FILE * infiles[], int numfiles, float *data,
                       double *dispdelays, int numsubbands,
                       int transpose, int *padding,
                       int *maskchans, int *nummasked, mask * obsmask)
/* This routine reads a record from the input files *infiles[]   */
/* which contain data from the GMRT system.  The routine uses    */
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
      if (!read_GMRT_rawblock(infiles, numfiles, rawdata, padding)) {
         printf("Problem reading the raw GMRT data file.\n\n");
         return 0;
      }
      if (0 != prep_GMRT_subbands(rawdata, data, dispdelays, numsubbands,
                                  transpose, maskchans, nummasked, obsmask)) {
         printf("Problem initializing prep_GMRT_subbands()\n\n");
         return 0;
      }
      firsttime = 0;
   }
   if (!read_GMRT_rawblock(infiles, numfiles, rawdata, padding)) {
      printf("Problem reading the raw GMRT data file.\n\n");
      return 0;
   }
   return prep_GMRT_subbands(rawdata, data, dispdelays, numsubbands,
                             transpose, maskchans, nummasked, obsmask);
}


void GMRT_hdr_to_inf(char *datfilenm, infodata * idata)
/* Convert GMRT header into an infodata structure */
{
   FILE *hdrfile;
   double ss;
   char line[200], ctmp[100], project[20], date[20];
   char *hdrfilenm;
   int numantennas, hh, mm, cliplen;

   hdrfilenm = calloc(strlen(datfilenm) + 1, 1);
   cliplen = strlen(datfilenm) - 3;
   strncpy(hdrfilenm, datfilenm, cliplen);
   strcpy(hdrfilenm + cliplen, "hdr");
   hdrfile = chkfopen(hdrfilenm, "r");
   while (fgets(line, 200, hdrfile)) {
      if (line[0] == '#') {
         continue;
      } else if (strncmp(line, "Site            ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", idata->telescope);
      } else if (strncmp(line, "Observer        ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", idata->observer);
      } else if (strncmp(line, "Proposal        ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", project);
      } else if (strncmp(line, "Array Mode      ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", idata->instrument);
      } else if (strncmp(line, "Observing Mode  ", 16) == 0) {
         continue;
      } else if (strncmp(line, "Date            ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", date);
      } else if (strncmp(line, "Num Antennas    ", 16) == 0) {
         sscanf(line, "%*[^:]: %d\n", &numantennas);
      } else if (strncmp(line, "Antenna List    ", 16) == 0) {
         continue;
      } else if (strncmp(line, "Num Channels    ", 16) == 0) {
         sscanf(line, "%*[^:]: %d\n", &idata->num_chan);
      } else if (strncmp(line, "Channel width   ", 16) == 0) {
         sscanf(line, "%*[^:]: %lf\n", &idata->chan_wid);
         /* Flag is made to handle LSB and USB data.... S.Sarala, 25 Apr 2005 */
         if (idata->chan_wid < 0.0)
            sb_flag = -1;
         else
            sb_flag = 1;
         if (idata->chan_wid < 0.0)
            idata->chan_wid = -(idata->chan_wid);
      } else if (strncmp(line, "Frequency Ch.1  ", 16) == 0) {
         sscanf(line, "%*[^:]: %lf\n", &idata->freq);
      } else if (strncmp(line, "Sampling Time   ", 16) == 0) {
         sscanf(line, "%*[^:]: %lf\n", &idata->dt);
         idata->dt /= 1000000.0;        /* Convert from us to s */
      } else if (strncmp(line, "Num bits/sample ", 16) == 0) {
         sscanf(line, "%*[^:]: %d\n", &bits);
         /* to handle 16 and 8 bit data.... S.Sarala, 28 May 2005 */
         if (bits == 16)
            bytesperpt_st = 2;
         else
            bytesperpt_st = 1;
      } else if (strncmp(line, "Data Format     ", 16) == 0) {
         {
            /* The following is from question 20.9 of the comp.lang.c FAQ */
            int x = 1, machine_is_little_endian = 0;
            if(*(char *)&x == 1) {
               machine_is_little_endian = 1;
            }

            if (strstr(line, "little")) {
               if (machine_is_little_endian)
                  need_byteswap_st = 0;
               else
                  need_byteswap_st = 1;
            } else {
               if (machine_is_little_endian)
                  need_byteswap_st = 1;
               else
                  need_byteswap_st = 0;
            }
         }
      } else if (strncmp(line, "Polarizations   ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", ctmp);
         if (strcmp(ctmp, "Total I")) {
            printf("\nWarning:  Cannot handle data with other than "
                   "'Total I' polarization.\n   Data is '%s'\n", ctmp);
         }
      } else if (strncmp(line, "MJD             ", 16) == 0) {
         sscanf(line, "%*[^:]: %d\n", &idata->mjd_i);
      } else if (strncmp(line, "UTC             ", 16) == 0) {
         sscanf(line, "%*[^:]: %d:%d:%lf\n", &hh, &mm, &ss);
         idata->mjd_f = (hh + (mm + ss / 60.0) / 60.0) / 24.0;
      } else if (strncmp(line, "Source          ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", idata->object);
      } else if (strncmp(line, "Coordinates     ", 16) == 0) {
         sscanf(line, "%*[^:]: %d:%d:%lf, %d:%d:%lf\n",
                &idata->ra_h, &idata->ra_m, &idata->ra_s,
                &idata->dec_d, &idata->dec_m, &idata->dec_s);
      } else if (strncmp(line, "Coordinate Sys  ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", ctmp);
         if (strcmp(ctmp, "J2000")) {
            printf("\nWarning:  Cannot non-J2000 coordinates!\n"
                   "   Data is '%s'\n", ctmp);
         }
      } else if (strncmp(line, "Drift Rate      ", 16) == 0) {
         continue;
      } else if (strncmp(line, "Obs. Length     ", 16) == 0) {
         continue;
      } else if (strncmp(line, "Bad Channels    ", 16) == 0) {
         sscanf(line, "%*[^:]: %[^\n]\n", ctmp);
      } else if (strncmp(line, "Bit shift value ", 16) == 0) {
         sscanf(line, "%*[^:]: %d\n", &bitshift);
      } else {
         continue;
      }
   }
   idata->freqband = idata->chan_wid * idata->num_chan;
   idata->dm = 0.0;
   idata->N = 0.0;
   idata->numonoff = 0;
   idata->bary = 0;
   idata->fov = 1.2 * SOL * 3600.0 / (1000000.0 * idata->freq * 45 * DEGTORAD);
   strcpy(idata->band, "Radio");
   strcpy(idata->analyzer, "Scott Ransom");
   strcpy(badch, ctmp);
   printf("bad channels are %s\n", badch);
   printf("bit shift value set to %d\n", bitshift);
   printf("No of bits/sample is %d\n", bits);

   sprintf(idata->notes, "%d antenna observation for Project %s\n"
           "    UT Date at file start = %s\n"
           "    Bad channels: %s\n", numantennas, project, date, ctmp);
   fclose(hdrfile);
   free(hdrfilenm);
}


void convert_GMRT_block(short *indata, unsigned char *outdata)
/* This routine converts 16 bit digitized data into bytes      */
/* Modified to handle both LSB and USB data...                 */
/* Also included featuring for zeroing bad channels given      */
/* in hdr file ...                                             */
/*               S.Sarala  and  Y. Gupta  29 Apr 2005          */
/* Included the clip option for GMRT data files...             */
/*               S.Sarala  and  Y. Gupta  9 May 2005           */
{
   int ii, jj;
   short inval;
   char str[100];
   char *ptr;
   int j, k, jump = 0;
   static int firsttime = 1, b[20][20], no;

   /* Note:  The LSB is the GPS bit.  The other bit that we throw */
   /*        away is actually significant.  But we need a byte... */

   if (firsttime) {
      printf("\n");
      no = 0;                   /* make sure that it does not crash evenif the bad channel string is empty... */
      printf("Inside convert_gmrt_block for FIRSTTIME \n");
      printf("BAD ch. string is %s\n", badch);
      sscanf(badch, "%d:%s", &no, str);

/*  parsing the string "str" for bands of bad channels   */
      if (no != 0) {            /* if the number of bands of bad channels are zero no need to go through this loop */
         j = 0;
         ptr = strtok(str, ",");
         sscanf(ptr, "%d-%d", &b[j][0], &b[j][1]);
         b[j][1] = (b[j][1] - b[j][0]) + 1;
         printf("band # %d : %d %d\n", j + 1, b[j][0], b[j][1]);

         while ((ptr = strtok(NULL, ",")) != NULL) {
            j++;
            sscanf(ptr, "%d-%d", &b[j][0], &b[j][1]);
            b[j][1] = (b[j][1] - b[j][0]) + 1;
            printf("band # %d : %d %d\n", j + 1, b[j][0], b[j][1]);
         }
      }

      firsttime = 0;
   }

   if (sb_flag < 0) {           /*   This Loop handles LSB data  */
      for (ii = 0; ii < numchan_st; ii++) {
         for (j = 0; j < no; j++) {
            if (ii == b[j][0] - 1) {    /* masking the band of bad channels  */
               for (k = ii; k < ii + b[j][1]; k++) {
                  for (jj = 0; jj < ptsperblk_st; jj++)
                     outdata[(jj * numchan_st) + (numchan_st - k - 1)] = 0.0;
               }
               ii = ii + b[j][1] - 1;
               jump = 1;        /* jump is introduced to make sure that ii does not exceed the value 255 */
            }
         }
         if (!jump) {
            for (jj = 0; jj < ptsperblk_st; jj++) {
               inval = indata[jj * numchan_st + ii] >> bitshift;
               outdata[(jj * numchan_st) + (numchan_st - ii - 1)] =
                   (inval > 255) ? 255 : inval;
            }
         } else {
            jump = 0;
         }
      }
   } else {                     /*   This Loop handles USB data   */
      for (ii = 0; ii < numchan_st; ii++) {
         for (j = 0; j < no; j++) {
            if (ii == b[j][0] - 1) {    /* masking the band of bad channels  */
               for (k = ii; k < ii + b[j][1]; k++) {
                  for (jj = 0; jj < ptsperblk_st; jj++)
                     outdata[jj * numchan_st + k] = 0.0;
               }
               ii = ii + b[j][1] - 1;
               jump = 1;
            }
         }
         if (!jump) {
            for (jj = 0; jj < ptsperblk_st; jj++) {
               inval = indata[jj * numchan_st + ii] >> bitshift;
               outdata[jj * numchan_st + ii] = (inval > 255) ? 255 : inval;
            }
         } else {
            jump = 0;
         }
      }                         /* end of the USB data loop     */
   }                            /* end of the loop covering all chans, all data points */
}


void convert_GMRT_block_8bit(unsigned char *indata, unsigned char *outdata)
/* This routine is modified to handle both LSB and USB data... */
/* Modified to handle both LSB and USB data...                 */
/* Also included featuring for zeroing bad channels given      */
/* in hdr file...                                              */
/*               S.Sarala  and  Y. Gupta  29 Apr 2005          */
/* Included the clip option for GMRT data files...             */
/*               S.Sarala  and  Y. Gupta  9 May 2005           */
/*                                                             */
/* Further modified convert_GMRT_block to work with 8-bit data */
/* input files...                                              */
/*               S.Sarala  and  Y. Gupta  20 May 2005	       */
{
   int ii, jj;
   unsigned char inval;
   char str[100];
   char *ptr;
   int j, k, jump = 0;
   static int firsttime = 1, b[20][20], no;

   /* Note:  The LSB is the GPS bit.  The other bit that we throw */
   /*        away is actually significant.  But we need a byte... */

   if (firsttime) {
      printf("\n");
      no = 0;                   /* make sure that it does not crash evenif the bad channel string is empty... */
      printf("Inside convert_gmrt_block for FIRSTTIME \n");
      printf("BAD ch. string is %s\n", badch);
      sscanf(badch, "%d:%s", &no, str);

/*  parsing the string "str" for bands of bad channels   */
      if (no != 0) {            /* if the number of bands of bad channels are zero no need to go through this loop */
         j = 0;
         ptr = strtok(str, ",");
         sscanf(ptr, "%d-%d", &b[j][0], &b[j][1]);
         b[j][1] = (b[j][1] - b[j][0]) + 1;
         printf("band # %d : %d %d\n", j + 1, b[j][0], b[j][1]);

         while ((ptr = strtok(NULL, ",")) != NULL) {
            j++;
            sscanf(ptr, "%d-%d", &b[j][0], &b[j][1]);
            b[j][1] = (b[j][1] - b[j][0]) + 1;
            printf("band # %d : %d %d\n", j + 1, b[j][0], b[j][1]);
         }
      }

      firsttime = 0;
   }

   if (sb_flag < 0) {           /*   This Loop handles LSB data  */
      for (ii = 0; ii < numchan_st; ii++) {
         for (j = 0; j < no; j++) {
            if (ii == b[j][0] - 1) {    /* masking the band of bad channels  */
               for (k = ii; k < ii + b[j][1]; k++) {
                  for (jj = 0; jj < ptsperblk_st; jj++)
                     outdata[(jj * numchan_st) + (numchan_st - k - 1)] = 0.0;
               }
               ii = ii + b[j][1] - 1;
               jump = 1;        /* jump is introduced to make sure that ii does not exceed the value 255 */
            }
         }
         if (!jump) {
            for (jj = 0; jj < ptsperblk_st; jj++) {
               inval = indata[jj * numchan_st + ii];
               outdata[(jj * numchan_st) + (numchan_st - ii - 1)] = inval;
            }
         } else {
            jump = 0;
         }
      }
   } else {                     /*   This Loop handles USB data   */
      for (ii = 0; ii < numchan_st; ii++) {
         for (j = 0; j < no; j++) {
            if (ii == b[j][0] - 1) {    /* masking the band of bad channels  */
               for (k = ii; k < ii + b[j][1]; k++) {
                  for (jj = 0; jj < ptsperblk_st; jj++)
                     outdata[jj * numchan_st + k] = 0.0;
               }
               ii = ii + b[j][1] - 1;
               jump = 1;
            }
         }
         if (!jump) {
            for (jj = 0; jj < ptsperblk_st; jj++) {
               inval = indata[jj * numchan_st + ii];
               outdata[jj * numchan_st + ii] = inval;
            }
         } else {
            jump = 0;
         }
      }                         /* end of the USB data loop     */
   }                            /* end of the loop covering all chans, all data points */
}
