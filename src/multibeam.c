#include "presto.h"
#include "multibeam.h"

int skip_to_multibeam_rec(FILE * infile, long rec)
/* This routine skips to the record 'rec' in the input file */
/* *infile.  *infile contains 1 bit digitized data from the */
/* multibeam receiver at Parkes                             */
/* Returns the record that was skipped to.                  */
{
  chkfileseek(infile, rec, RECLEN, SEEK_SET);
  return rec;
}


int read_multibeam_rec(FILE * infile, multibeam_tapehdr * hdr, \
		       unsigned char *data)
/* This routine reads a single multibeam record from          */
/* the input file *infile which contains 1 bit digitized data */
/* from the multibeam correlator at Parkes.                   */
/* Length of a multibeam record is 640 bytes for the header   */
/* plus 48k of data = 49792 bytes.                            */
{
  static long lastctr = 0;
  long currentctr;
  char tmp[50];

  /* Read the header */

  if (fread(hdr, sizeof(unsigned char), HDRLEN, infile) != HDRLEN) {
    if (feof(infile)) {
      return 0;
    } else {
      printf("\nProblem reading header from multibeam data file.\n");
      printf("Exiting\n");
    }
  }
  /* Verify that the records are in order and all present */

  sprintf(tmp, " %.8s ", hdr->blk_cntr);
  currentctr = strtol(tmp, NULL, 10);
  /*
  if (currentctr - lastctr != 1) {
    printf("\nError:  Raw data file is missing a record between\n");
    printf("        blocks %ld and %ld.  Exiting.\n", lastctr, currentctr);
    exit(1);
  }
  */
  lastctr = currentctr;

  /* Read the raw data */

  if (fread(data, sizeof(unsigned char), DATLEN, infile) != DATLEN) {
    if (feof(infile)) {
      return 0;
    } else {
      printf("\nProblem reading raw data from multibeam data file.\n");
      printf("Exiting\n");
    }
  }
  return 1;
}


int read_multibeam(FILE * infile, float *data, long numpts,
		   double *dispdelays, long numchan)
/* This routine reads a numpts record with numchan each from   */
/* the input file *infile which contains 1 bit digitized data  */
/* from the multibeam correlator at Parkes.                    */
/* It returns the number of points read.                       */
{
  static unsigned char raw[DATLEN], *ptr;
  static int firsttime = 1, recsize = 0, order = 0;
  static float *rawdata1, *rawdata2, *currentdata, *lastdata, *tmp;
  static long worklen = 0;
  static multibeam_tapehdr hdr;
  long i;

  if (firsttime) {
    recsize = numchan / 8;
    worklen = DATLEN * 8;
    /*    raw = gen_bvect(recsize * numpts); */

    /* Create our data storage for the raw floating point data */

    rawdata1 = gen_fvect(worklen);
    rawdata2 = gen_fvect(worklen);
    for (i = 0; i < worklen; i++) {
      rawdata1[i] = 0.0;
      rawdata2[i] = 0.0;
    }

    firsttime = 0;

    /* Read the first multibeam record to cope with end effects */

    if (!read_multibeam_rec(infile, &hdr, raw)) {
      printf("Problem reading the raw data file.\n\n");
      free(rawdata1);
      free(rawdata2);
      return 0;
    }
    /* If order is negative, then the data has been recorded with */
    /* high frequencies first.                                    */

    order = (int) (strtod(hdr.chanbw[1], NULL) * 10000.0);
    order /= abs(order);

    /* Initialize our data pointers */

    currentdata = rawdata1;
    lastdata = rawdata2;

    /* Convert the data to floats */

    ptr = raw;
    for (i = 0; i < numpts; i++, ptr += recsize) {
      convert_multibeam_point(ptr, currentdata + i * numchan);
    }

    /* Swap our data pointers */

    tmp = currentdata;
    currentdata = lastdata;
    lastdata = tmp;

  }
  /* Read a multibeam record */

  if (!read_multibeam_rec(infile, &hdr, raw)) {
    for (i = 0; i < worklen; i++)
      currentdata[i] = 0.0;
    dedisp(rawdata1, lastdata, numpts, dispdelays, numchan, \
	   data, order);
    free(rawdata1);
    free(rawdata2);
    return 0;
  }
  /* Convert the data to floats */

  ptr = raw;
  for (i = 0; i < numpts; i++, ptr += recsize) {
    convert_multibeam_point(ptr, currentdata + i * numchan);
  }

  /* De-disperse the data */

  dedisp(currentdata, lastdata, numpts, dispdelays, numchan, \
	 data, order);

  /* Swap our data pointers */

  tmp = currentdata;
  currentdata = lastdata;
  lastdata = tmp;

  /* Return the number of points we actually read */

  return numpts;
}


int read_mb_chan_to_vecs(FILE * infile, float *data,
			 long numpts, long numchan)
/* This routine reads a numpts record with numchan each from   */
/* the input file *infile which contains 1 bit digitized data  */
/* from the multibeam correlator at Parkes.  It stores the     */
/* data in a vector *data of length numchan*numpts.  The low   */
/* frequency channel is stored first, then the next highest    */
/* freq etc, with numpts floating points per frequency.        */
/* It returns the number of points (numpts) read if successful */
/* or 0 if unsuccessful.                                       */
{
  static unsigned char raw[DATLEN], *ptr;
  static int firsttime = 1, recsize = 0, order = 0;
  static float *rawdata = NULL, *tmp = NULL;
  static long worklen = 0;
  static multibeam_tapehdr hdr;
  long i, j;

  if (firsttime) {
    recsize = numchan / 8;
    worklen = numchan * numpts;
    /*    raw = gen_bvect(recsize * numpts); */

    /* Create our data storage for the raw floating point data */

    rawdata = gen_fvect(worklen);
    tmp = gen_fvect(worklen);

    firsttime = 0;

    /* Read the first multibeam record to cope with end effects */

    if (!read_multibeam_rec(infile, &hdr, raw)) {
      printf("Problem reading the raw data file.\n\n");
      free(rawdata);
      return 0;
    }

    /* If order is negative, then the data has been recorded with */
    /* high frequencies first.                                    */

    order = (int) (strtod(hdr.chanbw[1], NULL) * 10000.0);
    order /= abs(order);

  } else {

    /* Read a multibeam record */
    
    if (!read_multibeam_rec(infile, &hdr, raw)) {
      printf("Problem reading the raw data file.\n\n");
      free(rawdata);
      return 0;
    }

  }    

  /* Convert the data to floats */
  
  ptr = raw;
  for (i = 0; i < numpts; i++, ptr += recsize) {
    convert_multibeam_point(ptr, rawdata + i * numchan);
  }

  /* Change ordering so that all of the data for an */
  /* individual frequency is grouped together.      */

  if (order > 0) {
    
    for (i = 0 ; i < numchan ; i++) {
      for (j = 0 ; j < numpts ; j++) {
	data[i * numpts + j] = rawdata[j * numchan + i];
      }
    }

  } else {

    for (i = 0 ; i < numchan ; i++) {
      for (j = 0 ; j < numpts ; j++) {
	data[i * numpts + j] = rawdata[j * numchan + numchan - i - 1];
      }
    }

  }
  
  /* Return the number of points we actually read */
  
  return numpts;
}


void multibeam_hdr_to_inf(multibeam_tapehdr * hdr, infodata * idata)
/* Convert appropriate Multibeam header portions to make */
/* a basic '.inf' file                                   */
{
  double tmp1, tmp2;
  char ctmp[100];
  long ptsperrec = 0;

  sprintf(idata->object, "%.16s", hdr->pname);
  sscanf(hdr->ra_start, \
	 "%3d:%2d:%9lf", &idata->ra_h, &idata->ra_m, &idata->ra_s);
  sscanf(hdr->dec_start, \
	 "%3d:%2d:%9lf", &idata->dec_d, &idata->dec_m, &idata->dec_s);
  strcpy(idata->telescope, "Parkes");
  strcpy(idata->instrument, "Multibeam (center beam)");
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
    idata->freq -= idata->freqband - idata->chan_wid;
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
  sprintf(idata->notes, "Topo UT Date & Time at file start = %.8s, %.16s  From tape %.6s  file #%.4s  block #%.8s  Comment: %.64s\n", hdr->date, hdr->ut_start, hdr->tape_lbl, hdr->file_cntr, hdr->blk_cntr, hdr->comment);

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


void convert_multibeam_point(unsigned char *rec, float *data)
/* This routine converts 1 bit digitized data with */
/* 256 channels to a 256 array of floats.          */
{
    data[0] = rec[0] & 0x80 ? 1.0 : 0.0;
    data[1] = rec[0] & 0x40 ? 1.0 : 0.0;
    data[2] = rec[0] & 0x20 ? 1.0 : 0.0;
    data[3] = rec[0] & 0x10 ? 1.0 : 0.0;
    data[4] = rec[0] & 0x08 ? 1.0 : 0.0;
    data[5] = rec[0] & 0x04 ? 1.0 : 0.0;
    data[6] = rec[0] & 0x02 ? 1.0 : 0.0;
    data[7] = rec[0] & 0x01 ? 1.0 : 0.0;
    data[8] = rec[1] & 0x80 ? 1.0 : 0.0;
    data[9] = rec[1] & 0x40 ? 1.0 : 0.0;
    data[10] = rec[1] & 0x20 ? 1.0 : 0.0;
    data[11] = rec[1] & 0x10 ? 1.0 : 0.0;
    data[12] = rec[1] & 0x08 ? 1.0 : 0.0;
    data[13] = rec[1] & 0x04 ? 1.0 : 0.0;
    data[14] = rec[1] & 0x02 ? 1.0 : 0.0;
    data[15] = rec[1] & 0x01 ? 1.0 : 0.0;
    data[16] = rec[2] & 0x80 ? 1.0 : 0.0;
    data[17] = rec[2] & 0x40 ? 1.0 : 0.0;
    data[18] = rec[2] & 0x20 ? 1.0 : 0.0;
    data[19] = rec[2] & 0x10 ? 1.0 : 0.0;
    data[20] = rec[2] & 0x08 ? 1.0 : 0.0;
    data[21] = rec[2] & 0x04 ? 1.0 : 0.0;
    data[22] = rec[2] & 0x02 ? 1.0 : 0.0;
    data[23] = rec[2] & 0x01 ? 1.0 : 0.0;
    data[24] = rec[3] & 0x80 ? 1.0 : 0.0;
    data[25] = rec[3] & 0x40 ? 1.0 : 0.0;
    data[26] = rec[3] & 0x20 ? 1.0 : 0.0;
    data[27] = rec[3] & 0x10 ? 1.0 : 0.0;
    data[28] = rec[3] & 0x08 ? 1.0 : 0.0;
    data[29] = rec[3] & 0x04 ? 1.0 : 0.0;
    data[30] = rec[3] & 0x02 ? 1.0 : 0.0;
    data[31] = rec[3] & 0x01 ? 1.0 : 0.0;
    data[32] = rec[4] & 0x80 ? 1.0 : 0.0;
    data[33] = rec[4] & 0x40 ? 1.0 : 0.0;
    data[34] = rec[4] & 0x20 ? 1.0 : 0.0;
    data[35] = rec[4] & 0x10 ? 1.0 : 0.0;
    data[36] = rec[4] & 0x08 ? 1.0 : 0.0;
    data[37] = rec[4] & 0x04 ? 1.0 : 0.0;
    data[38] = rec[4] & 0x02 ? 1.0 : 0.0;
    data[39] = rec[4] & 0x01 ? 1.0 : 0.0;
    data[40] = rec[5] & 0x80 ? 1.0 : 0.0;
    data[41] = rec[5] & 0x40 ? 1.0 : 0.0;
    data[42] = rec[5] & 0x20 ? 1.0 : 0.0;
    data[43] = rec[5] & 0x10 ? 1.0 : 0.0;
    data[44] = rec[5] & 0x08 ? 1.0 : 0.0;
    data[45] = rec[5] & 0x04 ? 1.0 : 0.0;
    data[46] = rec[5] & 0x02 ? 1.0 : 0.0;
    data[47] = rec[5] & 0x01 ? 1.0 : 0.0;
    data[48] = rec[6] & 0x80 ? 1.0 : 0.0;
    data[49] = rec[6] & 0x40 ? 1.0 : 0.0;
    data[50] = rec[6] & 0x20 ? 1.0 : 0.0;
    data[51] = rec[6] & 0x10 ? 1.0 : 0.0;
    data[52] = rec[6] & 0x08 ? 1.0 : 0.0;
    data[53] = rec[6] & 0x04 ? 1.0 : 0.0;
    data[54] = rec[6] & 0x02 ? 1.0 : 0.0;
    data[55] = rec[6] & 0x01 ? 1.0 : 0.0;
    data[56] = rec[7] & 0x80 ? 1.0 : 0.0;
    data[57] = rec[7] & 0x40 ? 1.0 : 0.0;
    data[58] = rec[7] & 0x20 ? 1.0 : 0.0;
    data[59] = rec[7] & 0x10 ? 1.0 : 0.0;
    data[60] = rec[7] & 0x08 ? 1.0 : 0.0;
    data[61] = rec[7] & 0x04 ? 1.0 : 0.0;
    data[62] = rec[7] & 0x02 ? 1.0 : 0.0;
    data[63] = rec[7] & 0x01 ? 1.0 : 0.0;
    data[64] = rec[8] & 0x80 ? 1.0 : 0.0;
    data[65] = rec[8] & 0x40 ? 1.0 : 0.0;
    data[66] = rec[8] & 0x20 ? 1.0 : 0.0;
    data[67] = rec[8] & 0x10 ? 1.0 : 0.0;
    data[68] = rec[8] & 0x08 ? 1.0 : 0.0;
    data[69] = rec[8] & 0x04 ? 1.0 : 0.0;
    data[70] = rec[8] & 0x02 ? 1.0 : 0.0;
    data[71] = rec[8] & 0x01 ? 1.0 : 0.0;
    data[72] = rec[9] & 0x80 ? 1.0 : 0.0;
    data[73] = rec[9] & 0x40 ? 1.0 : 0.0;
    data[74] = rec[9] & 0x20 ? 1.0 : 0.0;
    data[75] = rec[9] & 0x10 ? 1.0 : 0.0;
    data[76] = rec[9] & 0x08 ? 1.0 : 0.0;
    data[77] = rec[9] & 0x04 ? 1.0 : 0.0;
    data[78] = rec[9] & 0x02 ? 1.0 : 0.0;
    data[79] = rec[9] & 0x01 ? 1.0 : 0.0;
    data[80] = rec[10] & 0x80 ? 1.0 : 0.0;
    data[81] = rec[10] & 0x40 ? 1.0 : 0.0;
    data[82] = rec[10] & 0x20 ? 1.0 : 0.0;
    data[83] = rec[10] & 0x10 ? 1.0 : 0.0;
    data[84] = rec[10] & 0x08 ? 1.0 : 0.0;
    data[85] = rec[10] & 0x04 ? 1.0 : 0.0;
    data[86] = rec[10] & 0x02 ? 1.0 : 0.0;
    data[87] = rec[10] & 0x01 ? 1.0 : 0.0;
    data[88] = rec[11] & 0x80 ? 1.0 : 0.0;
    data[89] = rec[11] & 0x40 ? 1.0 : 0.0;
    data[90] = rec[11] & 0x20 ? 1.0 : 0.0;
    data[91] = rec[11] & 0x10 ? 1.0 : 0.0;
    data[92] = rec[11] & 0x08 ? 1.0 : 0.0;
    data[93] = rec[11] & 0x04 ? 1.0 : 0.0;
    data[94] = rec[11] & 0x02 ? 1.0 : 0.0;
    data[95] = rec[11] & 0x01 ? 1.0 : 0.0;
    data[96] = rec[12] & 0x80 ? 1.0 : 0.0;
    data[97] = rec[12] & 0x40 ? 1.0 : 0.0;
    data[98] = rec[12] & 0x20 ? 1.0 : 0.0;
    data[99] = rec[12] & 0x10 ? 1.0 : 0.0;
    data[100] = rec[12] & 0x08 ? 1.0 : 0.0;
    data[101] = rec[12] & 0x04 ? 1.0 : 0.0;
    data[102] = rec[12] & 0x02 ? 1.0 : 0.0;
    data[103] = rec[12] & 0x01 ? 1.0 : 0.0;
    data[104] = rec[13] & 0x80 ? 1.0 : 0.0;
    data[105] = rec[13] & 0x40 ? 1.0 : 0.0;
    data[106] = rec[13] & 0x20 ? 1.0 : 0.0;
    data[107] = rec[13] & 0x10 ? 1.0 : 0.0;
    data[108] = rec[13] & 0x08 ? 1.0 : 0.0;
    data[109] = rec[13] & 0x04 ? 1.0 : 0.0;
    data[110] = rec[13] & 0x02 ? 1.0 : 0.0;
    data[111] = rec[13] & 0x01 ? 1.0 : 0.0;
    data[112] = rec[14] & 0x80 ? 1.0 : 0.0;
    data[113] = rec[14] & 0x40 ? 1.0 : 0.0;
    data[114] = rec[14] & 0x20 ? 1.0 : 0.0;
    data[115] = rec[14] & 0x10 ? 1.0 : 0.0;
    data[116] = rec[14] & 0x08 ? 1.0 : 0.0;
    data[117] = rec[14] & 0x04 ? 1.0 : 0.0;
    data[118] = rec[14] & 0x02 ? 1.0 : 0.0;
    data[119] = rec[14] & 0x01 ? 1.0 : 0.0;
    data[120] = rec[15] & 0x80 ? 1.0 : 0.0;
    data[121] = rec[15] & 0x40 ? 1.0 : 0.0;
    data[122] = rec[15] & 0x20 ? 1.0 : 0.0;
    data[123] = rec[15] & 0x10 ? 1.0 : 0.0;
    data[124] = rec[15] & 0x08 ? 1.0 : 0.0;
    data[125] = rec[15] & 0x04 ? 1.0 : 0.0;
    data[126] = rec[15] & 0x02 ? 1.0 : 0.0;
    data[127] = rec[15] & 0x01 ? 1.0 : 0.0;
    data[128] = rec[16] & 0x80 ? 1.0 : 0.0;
    data[129] = rec[16] & 0x40 ? 1.0 : 0.0;
    data[130] = rec[16] & 0x20 ? 1.0 : 0.0;
    data[131] = rec[16] & 0x10 ? 1.0 : 0.0;
    data[132] = rec[16] & 0x08 ? 1.0 : 0.0;
    data[133] = rec[16] & 0x04 ? 1.0 : 0.0;
    data[134] = rec[16] & 0x02 ? 1.0 : 0.0;
    data[135] = rec[16] & 0x01 ? 1.0 : 0.0;
    data[136] = rec[17] & 0x80 ? 1.0 : 0.0;
    data[137] = rec[17] & 0x40 ? 1.0 : 0.0;
    data[138] = rec[17] & 0x20 ? 1.0 : 0.0;
    data[139] = rec[17] & 0x10 ? 1.0 : 0.0;
    data[140] = rec[17] & 0x08 ? 1.0 : 0.0;
    data[141] = rec[17] & 0x04 ? 1.0 : 0.0;
    data[142] = rec[17] & 0x02 ? 1.0 : 0.0;
    data[143] = rec[17] & 0x01 ? 1.0 : 0.0;
    data[144] = rec[18] & 0x80 ? 1.0 : 0.0;
    data[145] = rec[18] & 0x40 ? 1.0 : 0.0;
    data[146] = rec[18] & 0x20 ? 1.0 : 0.0;
    data[147] = rec[18] & 0x10 ? 1.0 : 0.0;
    data[148] = rec[18] & 0x08 ? 1.0 : 0.0;
    data[149] = rec[18] & 0x04 ? 1.0 : 0.0;
    data[150] = rec[18] & 0x02 ? 1.0 : 0.0;
    data[151] = rec[18] & 0x01 ? 1.0 : 0.0;
    data[152] = rec[19] & 0x80 ? 1.0 : 0.0;
    data[153] = rec[19] & 0x40 ? 1.0 : 0.0;
    data[154] = rec[19] & 0x20 ? 1.0 : 0.0;
    data[155] = rec[19] & 0x10 ? 1.0 : 0.0;
    data[156] = rec[19] & 0x08 ? 1.0 : 0.0;
    data[157] = rec[19] & 0x04 ? 1.0 : 0.0;
    data[158] = rec[19] & 0x02 ? 1.0 : 0.0;
    data[159] = rec[19] & 0x01 ? 1.0 : 0.0;
    data[160] = rec[20] & 0x80 ? 1.0 : 0.0;
    data[161] = rec[20] & 0x40 ? 1.0 : 0.0;
    data[162] = rec[20] & 0x20 ? 1.0 : 0.0;
    data[163] = rec[20] & 0x10 ? 1.0 : 0.0;
    data[164] = rec[20] & 0x08 ? 1.0 : 0.0;
    data[165] = rec[20] & 0x04 ? 1.0 : 0.0;
    data[166] = rec[20] & 0x02 ? 1.0 : 0.0;
    data[167] = rec[20] & 0x01 ? 1.0 : 0.0;
    data[168] = rec[21] & 0x80 ? 1.0 : 0.0;
    data[169] = rec[21] & 0x40 ? 1.0 : 0.0;
    data[170] = rec[21] & 0x20 ? 1.0 : 0.0;
    data[171] = rec[21] & 0x10 ? 1.0 : 0.0;
    data[172] = rec[21] & 0x08 ? 1.0 : 0.0;
    data[173] = rec[21] & 0x04 ? 1.0 : 0.0;
    data[174] = rec[21] & 0x02 ? 1.0 : 0.0;
    data[175] = rec[21] & 0x01 ? 1.0 : 0.0;
    data[176] = rec[22] & 0x80 ? 1.0 : 0.0;
    data[177] = rec[22] & 0x40 ? 1.0 : 0.0;
    data[178] = rec[22] & 0x20 ? 1.0 : 0.0;
    data[179] = rec[22] & 0x10 ? 1.0 : 0.0;
    data[180] = rec[22] & 0x08 ? 1.0 : 0.0;
    data[181] = rec[22] & 0x04 ? 1.0 : 0.0;
    data[182] = rec[22] & 0x02 ? 1.0 : 0.0;
    data[183] = rec[22] & 0x01 ? 1.0 : 0.0;
    data[184] = rec[23] & 0x80 ? 1.0 : 0.0;
    data[185] = rec[23] & 0x40 ? 1.0 : 0.0;
    data[186] = rec[23] & 0x20 ? 1.0 : 0.0;
    data[187] = rec[23] & 0x10 ? 1.0 : 0.0;
    data[188] = rec[23] & 0x08 ? 1.0 : 0.0;
    data[189] = rec[23] & 0x04 ? 1.0 : 0.0;
    data[190] = rec[23] & 0x02 ? 1.0 : 0.0;
    data[191] = rec[23] & 0x01 ? 1.0 : 0.0;
    data[192] = rec[24] & 0x80 ? 1.0 : 0.0;
    data[193] = rec[24] & 0x40 ? 1.0 : 0.0;
    data[194] = rec[24] & 0x20 ? 1.0 : 0.0;
    data[195] = rec[24] & 0x10 ? 1.0 : 0.0;
    data[196] = rec[24] & 0x08 ? 1.0 : 0.0;
    data[197] = rec[24] & 0x04 ? 1.0 : 0.0;
    data[198] = rec[24] & 0x02 ? 1.0 : 0.0;
    data[199] = rec[24] & 0x01 ? 1.0 : 0.0;
    data[200] = rec[25] & 0x80 ? 1.0 : 0.0;
    data[201] = rec[25] & 0x40 ? 1.0 : 0.0;
    data[202] = rec[25] & 0x20 ? 1.0 : 0.0;
    data[203] = rec[25] & 0x10 ? 1.0 : 0.0;
    data[204] = rec[25] & 0x08 ? 1.0 : 0.0;
    data[205] = rec[25] & 0x04 ? 1.0 : 0.0;
    data[206] = rec[25] & 0x02 ? 1.0 : 0.0;
    data[207] = rec[25] & 0x01 ? 1.0 : 0.0;
    data[208] = rec[26] & 0x80 ? 1.0 : 0.0;
    data[209] = rec[26] & 0x40 ? 1.0 : 0.0;
    data[210] = rec[26] & 0x20 ? 1.0 : 0.0;
    data[211] = rec[26] & 0x10 ? 1.0 : 0.0;
    data[212] = rec[26] & 0x08 ? 1.0 : 0.0;
    data[213] = rec[26] & 0x04 ? 1.0 : 0.0;
    data[214] = rec[26] & 0x02 ? 1.0 : 0.0;
    data[215] = rec[26] & 0x01 ? 1.0 : 0.0;
    data[216] = rec[27] & 0x80 ? 1.0 : 0.0;
    data[217] = rec[27] & 0x40 ? 1.0 : 0.0;
    data[218] = rec[27] & 0x20 ? 1.0 : 0.0;
    data[219] = rec[27] & 0x10 ? 1.0 : 0.0;
    data[220] = rec[27] & 0x08 ? 1.0 : 0.0;
    data[221] = rec[27] & 0x04 ? 1.0 : 0.0;
    data[222] = rec[27] & 0x02 ? 1.0 : 0.0;
    data[223] = rec[27] & 0x01 ? 1.0 : 0.0;
    data[224] = rec[28] & 0x80 ? 1.0 : 0.0;
    data[225] = rec[28] & 0x40 ? 1.0 : 0.0;
    data[226] = rec[28] & 0x20 ? 1.0 : 0.0;
    data[227] = rec[28] & 0x10 ? 1.0 : 0.0;
    data[228] = rec[28] & 0x08 ? 1.0 : 0.0;
    data[229] = rec[28] & 0x04 ? 1.0 : 0.0;
    data[230] = rec[28] & 0x02 ? 1.0 : 0.0;
    data[231] = rec[28] & 0x01 ? 1.0 : 0.0;
    data[232] = rec[29] & 0x80 ? 1.0 : 0.0;
    data[233] = rec[29] & 0x40 ? 1.0 : 0.0;
    data[234] = rec[29] & 0x20 ? 1.0 : 0.0;
    data[235] = rec[29] & 0x10 ? 1.0 : 0.0;
    data[236] = rec[29] & 0x08 ? 1.0 : 0.0;
    data[237] = rec[29] & 0x04 ? 1.0 : 0.0;
    data[238] = rec[29] & 0x02 ? 1.0 : 0.0;
    data[239] = rec[29] & 0x01 ? 1.0 : 0.0;
    data[240] = rec[30] & 0x80 ? 1.0 : 0.0;
    data[241] = rec[30] & 0x40 ? 1.0 : 0.0;
    data[242] = rec[30] & 0x20 ? 1.0 : 0.0;
    data[243] = rec[30] & 0x10 ? 1.0 : 0.0;
    data[244] = rec[30] & 0x08 ? 1.0 : 0.0;
    data[245] = rec[30] & 0x04 ? 1.0 : 0.0;
    data[246] = rec[30] & 0x02 ? 1.0 : 0.0;
    data[247] = rec[30] & 0x01 ? 1.0 : 0.0;
    data[248] = rec[31] & 0x80 ? 1.0 : 0.0;
    data[249] = rec[31] & 0x40 ? 1.0 : 0.0;
    data[250] = rec[31] & 0x20 ? 1.0 : 0.0;
    data[251] = rec[31] & 0x10 ? 1.0 : 0.0;
    data[252] = rec[31] & 0x08 ? 1.0 : 0.0;
    data[253] = rec[31] & 0x04 ? 1.0 : 0.0;
    data[254] = rec[31] & 0x02 ? 1.0 : 0.0;
    data[255] = rec[31] & 0x01 ? 1.0 : 0.0;
}
