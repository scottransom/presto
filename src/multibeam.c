#include "presto.h"
#include "multibeam.h"

/* All of the following have an _st to indicate static */
static long long numpts_st[MAXPATCHFILES], padpts_st[MAXPATCHFILES], N_st;
static int numblks_st[MAXPATCHFILES];
static int currentfile, currentblock, numchan_st, ptsperblk_st, bytesperpt_st;
static double times_st[MAXPATCHFILES], mjds_st[MAXPATCHFILES]g;
static double elapsed_st[MAXPATCHFILES], T_st, dt_st;
static double startblk_st[MAXPATCHFILES], endblk_st[MAXPATCHFILES];
static infodata idata_st[MAXPATCHFILES];
static unsigned char mask1[512], mask2[512], maskblock[DATLEN];

#define SWAPMASK(maskptr)(maskptr = (maskptr==mask1) ? mask2 : mask1)

void get_pkmb_file_info(FILE *files[], int numfiles, long long *N, 
			int *ptsperblock, int *numchan, double *dt, 
			double *T, int output)
/* Read basic information into static variables and make padding      */
/* calculations for a set of PKMB rawfiles that you want to patch     */
/* together.  N, numchan, dt, and T are return values and include all */
/* the files with the required padding.  If output is true, prints    */
/* a table showing a summary of the values.                           */
{
  int ii;
  multibeam_tapehdr header;

  if (numfiles > MAXPATCHFILES){
    printf("\nThe number of input files (%d) is greater than \n", numfiles);
    printf("   MAXPATCHFILES=%d.  Exiting.\n\n", MAXPATCHFILES);
    exit(0);
  }
  memset(mask1, 0xAA, 512); /* 10101010... */
  memset(mask2, 0x55, 512); /* 01010101... */
  chkfread(&header, 1, HDRLEN, files[0]);
  rewind(files[0]);
  multibeam_hdr_to_inf(&header, &idata_st[0]);
  numchan_st = *numchan = idata_st[0].num_chan;
  ptsperblk_st = *ptsperblock = DATLEN * 8 / numchan_st;
  bytesperpt_st = DATLEN / ptsperblk_st;
  for (ii=0; ii<ptsperblk_st; ii+=2){
    memset(maskblock + ii     * bytesperpt_st, 0xAA, bytesperpt_st);
    memset(maskblock + (ii+1) * bytesperpt_st, 0x55, bytesperpt_st);
  }
  numblks_st[0] = chkfilelen(files[0], RECLEN);
  numpts_st[0] = numblks_st[0] * ptsperblk_st;
  N_st = numpts_st[0];
  dt_st = *dt = idata_st[0].dt;
  times_st[0] = numpts_st[0] * dt_st;
  mjds_st[0] = idata_st[0].mjd_i + idata_st[0].mjd_f;
  elapsed_st[0] = 0.0;
  startblk_st[0] = 0.0;
  endblk_st[0] = (double) numpts_st[0] / ptsperblk_st;
  padpts_st[0] = padpts_st[numfiles-1] = 0;
  for (ii=1; ii<numfiles; ii++){
    chkfread(&header, 1, HDRLEN, files[ii]);
    rewind(files[ii]);
    multibeam_hdr_to_inf(&header, &idata_st[ii]);
    if (idata_st[ii].num_chan != numchan_st){
      printf("Number of channels (file %d) is not the same!\n\n", ii+1);
    }
    if (idata_st[ii].dt != dt_st){
      printf("Sample time (file %d) is not the same!\n\n", ii+1);
    }
    numblks_st[ii] = chkfilelen(files[ii], RECLEN);
    numpts_st[ii] = numblks_st[ii] * ptsperblk_st;
    times_st[ii] = numpts_st[ii] * dt_st;
    mjds_st[ii] = idata_st[ii].mjd_i + idata_st[ii].mjd_f;
    elapsed_st[ii] = mjd_sec_diff(idata_st[ii].mjd_i, idata_st[ii].mjd_f,
				  idata_st[ii-1].mjd_i, idata_st[ii-1].mjd_f);
    padpts_st[ii-1] = (long long)((elapsed_st[ii]-times_st[ii-1]) / 
				  dt_st + 0.5);
    elapsed_st[ii] += elapsed_st[ii-1];
    endblk_st[ii] = (double) (N_st + numpts_st[ii]) / ptsperblk_st;
    N_st += numpts_st[ii] + padpts_st[ii-1];
    startblk_st[ii] = (double) (numpts_st[ii-1] + padpts_st[ii-1]) / 
      ptsperblk_st;
  }
  padpts_st[numfiles-1] = (long long) ceiling(endblk_st[numfiles-1]) * 
    ptsperblk_st - N_st;
  N_st += padpts_st[numfiles-1];
  *N = N_st;
  *T = T_st = N_st * dt_st;
  currentfile = currentblock = 0;
  if (output){
    printf("  Number of files = %d\n", numfiles);
    printf("     Points/block = %d\n", ptsperblk_st);
    printf("  Num of channels = %d\n", numchan_st);
    printf(" Total points (N) = %lld\n", N_st);
    printf(" Sample time (dt) = %-14.14g\n", dt_st);
    printf("   Total time (s) = %-14.14g\n\n", T_st);
    printf("File  Start Block    End Block      Points      Elapsed (s)      Time (s)            MJD           Padding\n");
    printf("----  ------------  ------------  ----------  --------------  --------------  ------------------  ----------\n");
    for (ii=0; ii<numfiles; ii++)
      printf("%2d  %12.12g  %12.12g  %10lld  %14.14g  %14.14g  %17.12f  %10lld\n", 
	     ii+1, startblk_st[ii], endblk_st[ii], numpts_st[ii], 
	     elapsed_st[ii], times_st[ii], mjds_st[ii], padpts_st[ii]);
    printf("\n");
  }
}


int skip_to_multibeam_rec(FILE * infile, int rec)
/* This routine skips to the record 'rec' in the input file */
/* *infile.  *infile contains 1 bit digitized data from the */
/* multibeam receiver at Parkes                             */
/* Returns the record that was skipped to.                  */
{
  chkfileseek(infile, rec, RECLEN, SEEK_SET);
  return rec;
}


int read_multibeam_rawblock(FILE *infiles, int numfiles, 
			    multibeam_tapehdr *hdr, \
			    unsigned char *data)
/* This routine reads a single multibeam record from          */
/* the input files *infiles which contain 1 bit digitized     */
/* data from the multibeam pulsar backend at Parkes.          */
/* Length of a multibeam record is 640 bytes for the header   */
/* plus 48k of data = 49792 bytes.                            */
/* The header of the record read is placed in hdr.            */
/* *data must be pre-allocated with size 48k * blocks_to_read */
{
  int ii, blocksread=0, offset, numtopad=0;
  unsigned char record[RECLEN];
  static unsigned char databuffer[DATLEN*2];
  static int bufferpts=0, padnum=0;
  unsigned char *maskptr=mask1;

  /* If our buffer array is offset from last time */
  /* copy the second part into the first.         */

  if (bufferpts)
    for (ii=0; ii<bufferpts; ii++)
      *(databuffer+ii) = *(databuffer+DATLEN+ii);

  /* make sure our current file number is valid */

  if (currentfile >= numfiles)
    return 0;

  /* First try and read the current file */
  
  if (fread(record, 1, RECLEN, infiles[currentfile]) != 1){
    /* Are we at the end of the current file? */
    if (feof(infiles[currentfile])){
      /* Does this file need padding? */
      numtopad = padpts_st[currentfile] - padnum;
      if (numtopad){
	/* Pad the data */
	if (numtopad >= ptsperblk_st){
	  /* Add a full record of padding */
	  offset = bufferpts * bytesperpt_st;
	  for (ii=0; ii<DATLEN; ii++)
	    *(databuffer+offset+ii) = *(maskblock+ii);
	  padnum += ptsperblk_st;
	} else if (numtopad >= ptsperblk_st - bufferpts){
	  /* Add the amount of padding we need to */
	  /* make our buffer_offset = 0           */
	  numtopad = ptsperblk_st - bufferpts;
	  for (ii=0; ii<numtopad; ii++){
	    offset = (bufferpts + ii) * bytesperpt_st;
	    for (jj=0; jj<bytesperpt_st; jj++){
	      *(databuffer+offset+jj) = *(maskptr+jj);
	      SWAPMASK(maskptr);
	    }
	  }
	  padnum += ptsperblk_st;
	  bufferpts = 0;
	} else {
	  /* Add the remainder of the padding and */
	  /* then get a block from the next file. */
	  for (ii=0; ii<numtopad; ii++){
	    offset = (bufferpts + ii) * bytesperpt_st;
	    for (jj=0; jj<bytesperpt_st; jj++){
	      *(databuffer+offset+jj) = *(maskptr+jj);
	      SWAPMASK(maskptr);
	    }
	  }
	  padnum = 0;
	  bufferpts += numtopad;
	  currentfile++;
	  return read_multibeam_rawblock(infiles, numfiles, hdr, data);
	}
	/* If done with padding reset padding variables */
	if (padnum == padpts_st[currentfile]){
	  padnum = 0;
	  currentfile++;
	}
	/* Copy the new data into the output array */
	for (ii=0; ii<DATLEN; ii++)
	  *(data+ii) = *(databuffer+ii);
      } else {
	/* Try reading the next file */
	currentfile++;
	return read_multibeam_rawblock(infiles, numfiles, hdr, data);
      }
    } else {
      printf("\nProblem reading record from multibeam data file.\n");
      printf("Exiting\n");
      exit(1);
    }
  } else {
    /* Put the new header into the header array */
    for (ii=0; ii<HDRLEN; ii++)
      *(hdr+ii) = *(record+ii);
    /* Put the new data into the databuffer or directly */
    /* into the return array if the bufferoffset is 0.  */
    if (bufferpts){
      offset = bufferpts * bytesperpt_st;
      for (ii=0; ii<DATLEN; ii++){
	*(databuffer+offset+ii) = *(record+HDRLEN+ii);
	*(data+ii) = *(databuffer+ii);
      }
    } else {
      for (ii=0; ii<DATLEN; ii++)
	*(data+ii) = *(record+HDRLEN+ii);
    }
  }
  return 1;
}


int read_multibeam(FILE * infile, float *data, int numpts,
		   double *dispdelays, int numchan)
/* This routine reads a numpts record with numchan each from   */
/* the input file *infile which contains 1 bit digitized data  */
/* from the multibeam correlator at Parkes.                    */
/* It returns the number of points read.                       */
{
  static unsigned char *raw, *ptr, *rawdata1, *rawdata2;
  static unsigned char *currentdata, *lastdata, *tmp;
  static int firsttime = 1, recsize = 0, blocklen = 0, worklen = 0;
  static int decreasing_f = 0, blocks_per_read = 0;
  static multibeam_tapehdr hdr;
  int ii, numread;

  if (firsttime) {
    recsize = numchan / 8;    /* bytes per time interval           */
    blocklen = DATLEN * 8;    /* total samples / block (1 bit per) */
    blocks_per_read = numpts * numchan / blocklen;
    worklen = blocks_per_read * blocklen; /* 1 byte per data point */

    /* Create our data storage for the raw byte data */
 
    raw  = gen_bvect(blocks_per_read * DATLEN);
    rawdata1 = gen_bvect(worklen);
    rawdata2 = gen_bvect(worklen);
    firsttime = 0;

    /* Read the first multibeam records to cope with end effects */

    numread = read_multibeam_recs(infile, &hdr, raw, blocks_per_read);
    if (numread != blocks_per_read) {
      printf("Problem reading the raw data file.\n\n");
      free(rawdata1);
      free(rawdata2);
      return 0;
    }

    /* If decreasing_f is true, then the data has been  */
    /* recorded with high frequencies first.            */

    decreasing_f = (strtod(hdr.chanbw[1], NULL) > 0.0) ? 0 : 1;

    /* Initialize our data pointers */

    currentdata = rawdata1;
    lastdata = rawdata2;

    /* Convert the data to 1 byte points */

    ptr = raw;
    for (ii = 0; ii < numpts; ii++, ptr += recsize)
      convert_multibeam_point(ptr, currentdata + ii * numchan, numchan, 
			      decreasing_f);

    /* Swap our data pointers */

    tmp = currentdata;
    currentdata = lastdata;
    lastdata = tmp;
  }

  /* Read the multibeam records */

  numread = read_multibeam_recs(infile, &hdr, raw, blocks_per_read);

  /* Convert the data to 1 byte points */

  ptr = raw;
  for (ii = 0; ii < numpts; ii++, ptr += recsize)
    convert_multibeam_point(ptr, currentdata + ii * numchan, numchan, 
			    decreasing_f);

  /* De-disperse the data */

  dedisp(currentdata, lastdata, numpts, numchan, dispdelays, data);

  /* Exit if we have reached the EOF */

  if (numread != blocks_per_read){
    free(rawdata1);
    free(rawdata2);
    return 0;
  }

  /* Swap our data pointers */

  tmp = currentdata;
  currentdata = lastdata;
  lastdata = tmp;

  /* Return the number of points we actually read */

  return numpts;
}


int read_rawmultibeam(FILE *infile, float *data, int numchan, 
		      int numblocks)
/* This routine reads numblocks PKMB records with numchan */
/* channels each from the input file *infile.  The number */
/* of blocks read is returned.                            */
{
  static unsigned char *rawblocks, *sample;
  static int firsttime=1, decreasing_f=-1, recsize, ptsperblock;
  static multibeam_tapehdr hdr;
  int ii, jj, numread, dataindex;

  if (firsttime){
    recsize = numchan / 8; /* bytes per time interval */
    ptsperblock = DATLEN / recsize;
    rawblocks = gen_bvect(DATLEN * numblocks);
    sample = gen_bvect(numchan);
    firsttime = 0;
  }
   
  /* Read raw multibeam records */

  numread = read_multibeam_recs(infile, &hdr, rawblocks, numblocks);

  /* If decreasing_f is true, then the data has been  */
  /* recorded with high frequencies first.            */
  
  if (decreasing_f < 0)
    decreasing_f = (strtod(hdr.chanbw[1], NULL) > 0.0) ? 0 : 1;

  /* Convert the data to floats */

  for (ii=0; ii<ptsperblock*numread; ii++){
    convert_multibeam_point(rawblocks+ii*recsize, 
			    sample, numchan, decreasing_f);
    dataindex = ii * numchan;
    for (jj=0; jj<numchan; jj++)
      data[dataindex+jj] = (float) sample[jj];
  }
  return numread;
}


int read_multibeam_subbands(FILE * infile, float *data, int numpts,
			    double *dispdelays, int numsubbands, 
			    int numchan)
/* This routine reads a numpts record with numchan each from     */
/* the input file *infile which contains 1 bit digitized data    */
/* from the multibeam correlator at Parkes.  The routine uses    */
/* dispersion delays in 'dispdelays' to de-disperse the data     */
/* into 'numsubbands' subbands.  It stores the resulting         */
/* data in vector 'data' of length numsubbands * numpts.  The    */
/* low frequency subband is stored first, then the next highest  */
/* subband etc, with 'numpts' floating points per subband.       */
/* It returns the number of points (pts_per_read) read if        */
/* succesful or 0 if unsuccessful.                               */
{
  static unsigned char *raw, *ptr, *rawdata1, *rawdata2, *move;
  static unsigned char *currentdata, *lastdata, *tmp;
  static int firsttime = 1, recsize = 0, blocklen = 0, worklen = 0;
  static int decreasing_f = 0, blocks_per_read = 0, move_size = 0;
  static multibeam_tapehdr hdr;
  short trtn;
  int ii, numread;

  if (firsttime) {
    recsize = numchan / 8;    /* bytes per time interval           */
    blocklen = DATLEN * 8;    /* total samples / block (1 bit per) */
    blocks_per_read = numpts * numchan / blocklen;
    worklen = blocks_per_read * blocklen; /* 1 byte per data point */
    move_size = (numpts + numsubbands) / 2;

    /* Create our data storage for the raw byte data */

    raw  = gen_bvect(blocks_per_read * DATLEN);
    rawdata1 = gen_bvect(worklen);
    rawdata2 = gen_bvect(worklen);
    move = gen_bvect(move_size);
    firsttime = 0;

    /* Read the first multibeam records to cope with end effects */

    numread = read_multibeam_recs(infile, &hdr, raw, blocks_per_read);
    if (numread != blocks_per_read) {
      printf("Problem reading the raw data file.\n\n");
      free(rawdata1);
      free(rawdata2);
      return 0;
    }

    /* If decreasing_f is true, then the data has been  */
    /* recorded with high frequencies first.            */

    decreasing_f = (strtod(hdr.chanbw[1], NULL) > 0.0) ? 0 : 1;

    /* Initialize our data pointers */

    currentdata = rawdata1;
    lastdata = rawdata2;

    /* Convert the data to floats */

    ptr = raw;
    for (ii = 0; ii < numpts; ii++, ptr += recsize)
      convert_multibeam_point(ptr, currentdata + ii * numchan, numchan, 
			      decreasing_f);

    /* Swap our data pointers */

    tmp = currentdata;
    currentdata = lastdata;
    lastdata = tmp;

  }

  /* Read the multibeam records */

  numread = read_multibeam_recs(infile, &hdr, raw, blocks_per_read);

  /* Convert the data to floats */

  ptr = raw;
  for (ii = 0; ii < numpts; ii++, ptr += recsize)
    convert_multibeam_point(ptr, currentdata + ii * numchan, numchan, 
			    decreasing_f);

  /* De-disperse the data into subbands */

  dedisp_subbands(currentdata, lastdata, numpts, numchan, 
		  dispdelays, numsubbands, data);

  /* Transpose the data into vectors in the result array */

  if ((trtn = transpose_float(data, numpts, numsubbands,
			      move, move_size))<0){
    printf("Error %d in transpose_float().\n",trtn);
  }

  /* Exit if we have reached the EOF */

  if (numread != blocks_per_read){
    free(rawdata1);
    free(rawdata2);
    free(move);
    return 0;
  }

  /* Swap our data pointers */

  tmp = currentdata;
  currentdata = lastdata;
  lastdata = tmp;

  /* Return the number of points we actually read */

  return numpts;
}


void multibeam_hdr_to_inf(multibeam_tapehdr * hdr, infodata * idata)
/* Convert appropriate Multibeam header portions to make */
/* a basic '.inf' file                                   */
{
  double tmp1, tmp2;
  char ctmp[100];
  int ptsperrec = 0;
  int itmp1, itmp2;

  sprintf(idata->object, "%.16s", hdr->pname);
  sscanf(hdr->ra_start, \
	 "%3d:%2d:%9lf", &idata->ra_h, &idata->ra_m, &idata->ra_s);
  sscanf(hdr->dec_start, \
	 "%3d:%2d:%9lf", &idata->dec_d, &idata->dec_m, &idata->dec_s);
  strcpy(idata->telescope, "Parkes");
  itmp1 = strtol(hdr->ibeam, NULL, 10);  
  itmp2 = strtol(hdr->nbeam, NULL, 10);  
  sprintf(idata->instrument, "Multibeam (Beam %d of %d)", itmp1, itmp2);
  sscanf(hdr->nchan[0], "%4d", &idata->num_chan);
  ptsperrec = DATLEN / (idata->num_chan / 8);
  sscanf(hdr->samp_int[0], "%12lf", &idata->dt);
  idata->dt *= 0.001;
  sscanf(hdr->nblk_read, "%8lf", &tmp1);
  sprintf(ctmp, " %.8s ", hdr->obs_time);
  tmp2 = strtod(ctmp, NULL);
  if (tmp1) {
    idata->N = ptsperrec * tmp1;
  } else {
    idata->N = tmp2 / idata->dt;
  }
  tmp1 = strtod(hdr->chanbw[0], NULL);
  idata->chan_wid = fabs(tmp1);
  idata->freqband = idata->num_chan * idata->chan_wid;
  idata->freq = strtod(hdr->freq_ch1[0], NULL);
  if (tmp1 < 0.0) {
    idata->freq -= (idata->freqband - idata->chan_wid);
  }
  idata->fov = 1.2 * SOL * 3600.0 / (1000000.0 * idata->freq * 64 * DEGTORAD);
  idata->mjd_i = strtol(hdr->mjd, NULL, 10);
  sscanf(hdr->ut_start, "%lf:%lf:%lf", &idata->mjd_f, &tmp1, &tmp2);
  idata->mjd_f = (idata->mjd_f + (tmp1 + tmp2 / 60.0) / 60.0) / 24.0;
  sprintf(ctmp, " %.8s ", hdr->blk_cntr);
  tmp1 = strtod(ctmp, NULL) - 1.0;
  idata->mjd_f += (tmp1 * ptsperrec * idata->dt) / SECPERDAY;
  if (idata->mjd_f >= 1.0){
    idata->mjd_f -= 1.0;
    idata->mjd_i += 1;
  }
  idata->bary = 0;
  idata->numonoff = 1;
  idata->onoff = (double *)malloc(2*sizeof(double));
  idata->onoff[0] = 0;
  idata->onoff[1] = idata->N - 1;
  strcpy(idata->band, "Radio");
  strcpy(idata->analyzer, "--");
  strcpy(idata->observer, "--");
  sprintf(idata->notes, "Topo UT Date & Time at file start = %.8s, %.16s\n    From tape %.6s  file #%.4s  block #%.8s\n    Comment: %.64s\n", hdr->date, hdr->ut_start, hdr->tape_lbl, hdr->file_cntr, hdr->blk_cntr, hdr->comment);

}



void print_multibeam_hdr(multibeam_tapehdr * hdr)
/* Output in human readable form a multibeam header. */
{
  printf("\n");
  printf(" Program name                                      = %.6s\n", \
	 hdr->program);
  printf(" Version number 'vv.v'                             = %.4s\n", \
	 hdr->version);
  printf(" Tape label                                        = %.6s\n", \
	 hdr->tape_lbl);
  printf(" File number on tape                               = %.4s\n", \
	 hdr->file_cntr);
  printf(" Block counter per beam in file                    = %.8s\n", \
	 hdr->blk_cntr);
  printf(" UT Date 'yymmdd  '                                = %.8s\n", \
	 hdr->date);
  printf(" MJD at 00h UT                                     = %.8s\n", \
	 hdr->mjd);
  printf(" UT at file start  ' hh:mm:ss.ssssss'              = %.16s\n", \
	 hdr->ut_start);
  printf(" LST at file start ' hh:mm:ss.ss'                  = %.12s\n", \
	 hdr->lst_start);
  printf(" Coord system ' 4' = Galactic, ' 5' = J2000        = %.2s\n", \
	 hdr->coord_type);
  printf(" RA at file start  ' hh:mm:ss.ssss'                = %.16s\n", \
	 hdr->ra_start);
  printf(" Dec at file start '-dd:mm:ss.sss'                 = %.16s\n", \
	 hdr->dec_start);
  printf(" Galactic long at file start (deg) 'ddd.dddd'      = %.8s\n", \
	 hdr->gl_start);
  printf(" Galactic lat. at file start (deg) '-dd.dddd'      = %.8s\n", \
	 hdr->gb_start);
  printf(" FA, PA or GPA (deg) '-ddd.ddd'                    = %.8s\n", \
	 hdr->fangle);
  printf(" Obs length 'ssss.sss'                             = %.8s\n", \
	 hdr->obs_time);
  printf(" Comment                                           = %.64s\n", \
	 hdr->comment);
  printf(" Number of filter systems, ' 1' or ' 2'            = %.2s\n", \
	 hdr->nfilter);
  printf(" Channel incrt (MHz, -ve if inverted) '-b.bbbbb'   = %.8s\n", \
	 hdr->chanbw[0]);
  printf(" Nr of channels in each filter system              = %.4s\n", \
	 hdr->nchan[0]);
  printf(" RF of first channel centre 'ffffff.fffff'         = %.12s\n", \
	 hdr->freq_ch1[0]);
  printf(" Samp intval in ms (1 = 2 always) 'mm.mmmmmmmmm'   = %.12s\n", \
	 hdr->samp_int[0]);
  printf(" Samp per group (For backwards compat., now = 1)   = %.4s\n", \
	 hdr->samp_grp[0]);
  printf(" Samples per block (Was groups per block)          = %.8s\n", \
	 hdr->samp_blk);
  printf(" Seconds per tape block per beam 'ss.sssss'        = %.8s\n", \
	 hdr->blk_sec);
  printf(" 0 = none, 1 = fixed FA, 2 = fixed PA or GPA       = %.2s\n", \
	 hdr->fdctrl);
  printf(" Data type code ('2'=pol, '3'=normal, '5'=dedisp)  = %.1s\n", \
	 hdr->data[0]);
  printf(" UT at blk start ' hh:mm:ss.ssssss'                = %.16s\n", \
	 hdr->ut_blk);
  printf(" LST at blk start ' hh:mm:ss.ss'                   = %.12s\n", \
	 hdr->lst_blk);
  printf(" RA at blk start ' hh:mm:ss.ssss'                  = %.16s\n", \
	 hdr->ra_blk);
  printf(" Dec at blk start '-dd:mm:ss.sss'                  = %.16s\n", \
	 hdr->dec_blk);
  printf(" Gal long. at blk start 'ddd.dddd'                 = %.8s\n", \
	 hdr->gl_blk);
  printf(" Gal lat at blk start '-dd.dddd'                   = %.8s\n", \
	 hdr->gb_blk);
  printf(" Zenith angle at blk start (deg) 'ddd.dddd'        = %.8s\n", \
	 hdr->zen_blk);
  printf(" Azimuth angle at blk start (deg) 'ddd.dddd'       = %.8s\n", \
	 hdr->az_blk);
  printf(" Attenuator settings (A,B) for Sys 1,2 'dd.d'      = %.4s\n", \
	 hdr->atten[0]);
  printf(" Total powers for 32-channel groups 'iiii'         = %.4s\n", \
	 hdr->tpower[0]);
  printf(" Nr of tape blocks in disk file                    = %.8s\n", \
	 hdr->nblk_read);
  printf(" Scan rate in `long' (deg/min) '-r.rrrrr'          = %.8s\n", \
	 hdr->scan_ratel);
  printf(" Scan rate in `lat'  (deg/min) '-r.rrrrr'          = %.8s\n", \
	 hdr->scan_rateb);
  printf(" Total number of beams                             = %.4s\n", \
	 hdr->nbeam);
  printf(" Beam number                                       = %.4s\n", \
	 hdr->ibeam);
  printf(" Psrname (B or J) or pointing name                 = %.16s\n", \
	 hdr->pname);
  printf(" Config file name                                  = %.8s\n", \
	 hdr->cfg_file);
  printf(" Schedule file name                                = %.8s\n", \
	 hdr->schd_file);
  printf(" Nr of bits per sample                             = %.2s\n", \
	 hdr->nbits[0]);
  printf(" DM for dedispersion 'ddddd.dd'                    = %.8s\n", \
	 hdr->ddm);
  printf(" Nr of channels per dedispersed band               = %.4s\n", \
	 hdr->nddch[0]);
  printf(" Tel status 'k0': On posn k(0-9)                   = %.2s\n", \
	 hdr->move);
  printf("            '01': On grid posn\n");
  printf("            '02': Off position\n");
  printf("            '03': GC scan\n");
  printf("            '04': SC scan\n");
  printf(" Pointing error (arcmin) 'mmm.mm'                  = %.6s\n", \
	 hdr->pnterr);
  printf(" Dedisp type: '  ': Norm, ' T': Tree, ' D': PDM    = %.2s\n", \
	 hdr->tree);
  printf(" Filter system for disk file                       = %.2s\n", \
	 hdr->ifilter);
  printf(" Telescope ID (as in tel.def)                      = %.10s\n", \
	 hdr->telid);
  printf(" Parallactic angle (deg) '-ddd.ddd'                = %.8s\n", \
	 hdr->pangle);
  printf(" Nr of samples per beam switch                     = %.8s\n", \
	 hdr->bsw);
  printf(" Cal cycle period in samples                       = %.4s\n", \
	 hdr->cal_cycle);
  printf(" Spare                                             = %.22s\n", \
	 hdr->spare);
  printf("\n");
}


void convert_multibeam_point(unsigned char *rec, unsigned char *data, \
			     int numchan, int decreasing_f)
/* This routine converts 1 bit digitized data with 'numchan' */
/* channels to an array of 'numchan' floats.                 */
{
  int ii, jj;

  if (decreasing_f){
    for(ii = numchan/8-1, jj = 0; ii >= 0; ii--, jj+=8){
      data[jj] = (rec[ii] >> 7) & 1;
      data[jj+1] = (rec[ii] >> 6) & 1;
      data[jj+2] = (rec[ii] >> 5) & 1;
      data[jj+3] = (rec[ii] >> 4) & 1;
      data[jj+4] = (rec[ii] >> 3) & 1;
      data[jj+5] = (rec[ii] >> 2) & 1;
      data[jj+6] = (rec[ii] >> 1) & 1;
      data[jj+7] = rec[ii] & 1;
    }
  } else {
    for(ii = 0, jj = 0; ii < numchan/8; ii++, jj+=8){
      data[jj] = rec[ii] & 1;
      data[jj+1] = (rec[ii] >> 1) & 1;
      data[jj+2] = (rec[ii] >> 2) & 1;
      data[jj+3] = (rec[ii] >> 3) & 1;
      data[jj+4] = (rec[ii] >> 4) & 1;
      data[jj+5] = (rec[ii] >> 5) & 1;
      data[jj+6] = (rec[ii] >> 6) & 1;
      data[jj+7] = (rec[ii] >> 7) & 1;
    }
  }
}
