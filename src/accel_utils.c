#include "accel.h"

#if defined (__GNUC__)
#  define inline __inline__
#else
#  undef inline
#endif

#define NEAREST_INT(x) (int) (x < 0 ? x - 0.5 : x + 0.5)

static inline int calc_required_z(int numharm, int harmnum, double zfull)
/* Calculate the 'z' you need for subharmonic     */
/* 'harmnum' out of 'numharm' subharmonics if the */
/* 'z' at the fundamental harmonic is 'zfull'.    */
{
  double zz;
  
  zz = ACCEL_RDZ * (zfull / numharm) * harmnum;
  return NEAREST_INT(zz) * ACCEL_DZ;
}

static inline double calc_required_r(int numharm, int harmnum, double rfull)
/* Calculate the 'r' you need for subharmonic     */
/* 'harmnum' out of 'numharm' subharmonics if the */
/* 'r' at the fundamental harmonic is 'rfull'.    */
{
  return (int) ((ACCEL_RDR * (rfull / numharm) * harmnum) + 0.5) * ACCEL_DR;
}


static inline int index_from_r(double r, double lor)
/* Return an index for a Fourier Freq given an array that */
/* has stepsize ACCEL_DR and low freq 'lor'.              */
{
  return (int) ((r - lor) * ACCEL_RDR + DBLCORRECT);
}


static inline int index_from_z(double z, double loz)
/* Return an index for a Fourier Fdot given an array that */
/* has stepsize ACCEL_DZ and low freq 'lor'.              */
{
  return (int) ((z - loz) * ACCEL_RDZ + DBLCORRECT);
}


static void compare_rzw_cands(fourierprops * list, int nlist, char *notes)
{
  int ii, jj, kk;
  char tmp[30];

  for (ii = 0; ii < nlist; ii++) {
    for (jj = 0; jj < nlist; jj++) {
      if (ii == jj)
	continue;
      if (fabs(list[ii].r - list[jj].r) < 15.0 && \
	  fabs(list[ii].z - list[jj].z) > 1.0 && \
	  list[ii].pow > list[jj].pow) {
	if (strncmp(notes+jj*20, "                      ", 20)==0){
	  sprintf(tmp, "SL? of Cand %d", ii+1);
	  strncpy(notes+jj*20, tmp, 20);
	}
	continue;
      }
      for (kk = 1; kk < 61; kk++) {
	if ((fabs(list[ii].r - list[jj].r / kk) < list[jj].rerr * 3) && \
	    (fabs(list[ii].z - list[jj].z / kk) < list[jj].zerr * 2)) {
	  if (strncmp(notes+jj*20, "                      ", 20)==0){
	    sprintf(tmp, "H %d of Cand %d", kk, ii+1);
	    strncpy(notes+jj*20, tmp, 20);
	    break;
	  }
	}
      }
    }
  }
}


static int calc_fftlen(int numharm, int harmnum, int max_zfull)
/* The fft length needed to properly process a subharmonic */
{
  int bins_needed, end_effects;

  bins_needed = (ACCEL_USELEN * harmnum) / numharm + 2;
  end_effects = z_resp_halfwidth(calc_required_z(numharm, harmnum, max_zfull), 
				 LOWACC) * 2 * ACCEL_NUMBETWEEN;
  return next2_to_n(bins_needed + end_effects);
}


static void init_kernel(int z, int fftlen, kernel *kern)
{
  int numkern;
  fcomplex *tempkern;

  kern->z = z;
  kern->fftlen = fftlen;
  kern->numbetween = ACCEL_NUMBETWEEN;
  kern->kern_half_width = z_resp_halfwidth((double) z, LOWACC);
  numkern = 2 * kern->numbetween * kern->kern_half_width;
  kern->numgoodbins = kern->fftlen - numkern;
  kern->data = gen_cvect(kern->fftlen);
  tempkern = gen_z_response(0.0, kern->numbetween, kern->z, numkern);
  place_complex_kernel(tempkern, numkern, kern->data, kern->fftlen);
  free(tempkern);
  COMPLEXFFT(kern->data, kern->fftlen, -1);
}


static void free_kernel(kernel *kern)
{
  free(kern->data);
}


static void init_subharminfo(int numharm, int harmnum, 
			     int zmax, subharminfo *shi)
/* Note:  'zmax' is the overall maximum 'z' in the search */
{
  int ii, fftlen;

  shi->numharm = numharm;
  shi->harmnum = harmnum;
  shi->zmax = calc_required_z(numharm, harmnum, zmax);
  shi->numkern = (shi->zmax / ACCEL_DZ) * 2 + 1;
  shi->kern = (kernel *)malloc(shi->numkern * sizeof(kernel));
  fftlen = calc_fftlen(numharm, harmnum, zmax); 
  for (ii=0; ii<shi->numkern; ii++)
    init_kernel(-shi->zmax+ii*ACCEL_DZ, fftlen, &shi->kern[ii]);
}


subharminfo **create_subharminfos(int numharm, int zmax)
{
  int ii, jj;
  subharminfo **shis;

  shis = (subharminfo **)malloc((numharm + 1) * sizeof(subharminfo *));
  /* Prep the fundamental */
  shis[1] = (subharminfo *)malloc(2 * sizeof(subharminfo));
  init_subharminfo(1, 1, zmax, &shis[1][1]);
  printf("  Fundamental  has %3d kernels from z = %4d to %4d\n", 
	 shis[1][1].numkern, -shis[1][1].zmax, shis[1][1].zmax);
  /* Prep the sub-harmonics */
  for (ii=2; ii<=numharm; ii++){
    shis[ii] = (subharminfo *)malloc(ii * sizeof(subharminfo));
    for (jj=1; jj<ii; jj++){
      if (jj==1 || ii % jj){
	init_subharminfo(ii, jj, zmax, &shis[ii][jj]);
	printf("  Harmonic %d/%d has %3d kernels from z = %4d to %4d\n", 
	       jj, ii, shis[ii][jj].numkern, 
	       -shis[ii][jj].zmax, shis[ii][jj].zmax);
      }
    }
  }
  return shis;
}


static void free_subharminfo(subharminfo *shi)
{
  int ii;

  for (ii=0; ii<shi->numkern; ii++)
    free_kernel(&shi->kern[ii]);
  free(shi->kern);
}


void free_subharminfos(int numharm, subharminfo **shis)
{
  int ii, jj;

  /* Free the sub-harmonics */
  for (ii=2; ii<=numharm; ii++){
    for (jj=1; jj<ii; jj++){
      if (jj==1 || ii % jj)
	free_subharminfo(&shis[ii][jj]);
    }
    free(shis[ii]);
  }
  /* Free the fundamental */
  free_subharminfo(&shis[1][1]);
  free(shis[1]);
  /* Free the container */
  free(shis);
}


static accelcand *create_accelcand(float power, float sigma, 
				   int numharm, double r, double z)
{
  accelcand *obj;

  obj = (accelcand *)malloc(sizeof(accelcand));
  obj->power = power;
  obj->sigma = sigma;
  obj->numharm = numharm;
  obj->r = r;
  obj->z = z;
  obj->pows = NULL;
  obj->hirs = NULL;
  obj->hizs = NULL;
  obj->derivs = NULL;
  return obj;
}

void free_accelcand(gpointer data, gpointer user_data)
{
  user_data = NULL;
  if (((accelcand *)data)->pows){
    free(((accelcand *)data)->pows);
    free(((accelcand *)data)->hirs);
    free(((accelcand *)data)->hizs);
    free(((accelcand *)data)->derivs);
  }
  free((accelcand *)data);
}


static int compare_accelcand_sigma(gconstpointer ca, gconstpointer cb)
/* Sorts from high to low sigma (ties are sorted by increasing r) */
{
  int result;
  accelcand *a, *b;

  a = (accelcand *) ca;
  b = (accelcand *) cb;
  result = (a->sigma < b->sigma) - (a->sigma > b->sigma);
  if (result) 
    return result;
  else
    return (a->r > b->r) - (a->r < b->r);
}


GSList *sort_accelcands(GSList *list)
/* Sort the candidate list by decreasing sigma */
{
  return g_slist_sort(list, compare_accelcand_sigma);
}


static GSList *insert_new_accelcand(GSList *list, float power, float sigma, 
				    int numharm, double rr, double zz)
/* Checks the current list to see if there is already */
/* a candidate within ACCEL_CLOSEST_R bins.  If not,  */
/* it adds it to the list in increasing freq order.   */
{
  GSList *tmp_list=list, *prev_list=NULL, *new_list;
  
  if (!list){
    new_list = g_slist_alloc();
    new_list->data = (gpointer *)create_accelcand(power, sigma, 
						  numharm, rr, zz);
    return new_list;
  }

  while ((tmp_list->next) && 
	 (((accelcand *)(tmp_list->data))->r < (rr - ACCEL_CLOSEST_R))){
    prev_list = tmp_list;
    tmp_list = tmp_list->next;
  }

  /* Similar candidate is present */

  if (fabs(rr - ((accelcand *)(tmp_list->data))->r) < ACCEL_CLOSEST_R){
    if (((accelcand *)(tmp_list->data))->sigma < sigma){
      free_accelcand(tmp_list->data, NULL);  /* Overwrite the old candidate */
      tmp_list->data = (gpointer *)create_accelcand(power, sigma, 
						    numharm, rr, zz);
    }
  } else {  /* This is a new candidate */
    new_list = g_slist_alloc();
    new_list->data = (gpointer *)create_accelcand(power, sigma, 
						  numharm, rr, zz);

    if (!tmp_list->next && 
	(((accelcand *)(tmp_list->data))->r < (rr - ACCEL_CLOSEST_R))){
      tmp_list->next = new_list;
      return list;
    }
    if (prev_list){
      prev_list->next = new_list;
      new_list->next = tmp_list;
    } else {
      new_list->next = list;
      return new_list;
    }
  }
  return list;
}


void optimize_accelcand(accelcand *cand, accelobs *obs)
{
  int ii;
  
  cand->pows = gen_dvect(cand->numharm);
  cand->hirs = gen_dvect(cand->numharm);
  cand->hizs = gen_dvect(cand->numharm);
  cand->derivs = (rderivs *)malloc(sizeof(rderivs)*cand->numharm);
  for (ii=0; ii<cand->numharm; ii++){
    cand->pows[ii] = max_rz_file(obs->fftfile, 
				 cand->r * (ii + 1) - obs->lobin,
				 cand->z * (ii + 1), 
				 &(cand->hirs[ii]), 
				 &(cand->hizs[ii]), 
				 &(cand->derivs[ii]));
    cand->hirs[ii] += obs->lobin;
  }
  cand->sigma = candidate_sigma(cand->power, cand->numharm, 
				obs->numindep[cand->numharm-1]);
}


static void center_string(char *outstring, char *instring, int width)
{
  int len;
  char *tmp;

  len = strlen(instring);
  if (width < len){
    printf("\nwidth < len (%d) in center_string(outstring, '%s', width=%d)\n",
	   len, instring, width);
  }
  tmp = memset(outstring, ' ', width);
  outstring[width] = '\0';
  if (len >= width){
    strncpy(outstring, instring, width);
  } else {
    strncpy(outstring+(width-len)/2, instring, len);
  }
}


static void write_val_with_err(FILE *outfile, double val, double err, 
			       int numerr, int width)
{
  int retval;
  char tmpstr[30], ctrstr[30];

  if (numerr==1)
    retval = nice_output_1(tmpstr, val, err, 0);
  else if (numerr==2)
    retval = nice_output_2(tmpstr, val, err, 0);
  else
    printf("\numerr = %d is out-of-range (1-2) in write_val_with_err()\n",
	   numerr);
  center_string(ctrstr, tmpstr, width);
  fprintf(outfile, "%s  ", ctrstr);
}


void output_fundamentals(fourierprops *props, GSList *list, 
			 accelobs *obs, infodata *idata)
{
  double accel=0.0, accelerr=0.0;
  int ii, numcols=11, numcands;
  int widths[11]={4, 5, 6, 4, 16, 15, 15, 15, 11, 15, 20};
  int errors[11]={0, 0, 0, 0,  1,  1,  2,  1,  2,  2,  0};
  char tmpstr[30], ctrstr[30], *notes;
  accelcand *cand;
  GSList *listptr;
  rzwerrs errs;
  static char *titles1[]={"", "", "Summed", "Num", "Period", 
			  "Frequency", "FFT 'r'", "Freq Deriv", "FFT 'z'", 
			  "Accel", ""};
  static char *titles2[]={"Cand", "Sigma", "Power", "Harm", "(ms)",
			  "(Hz)", "(bin)", "(Hz/s)", "(bins)", 
			  "(m/s)", "Notes"};

  numcands = g_slist_length(list);
  listptr = list;

  /* Close the old work file and open the cand file */
  
  fclose(obs->workfile);
  obs->workfile = chkfopen(obs->accelnm, "w");
  
  /* Set our candidate notes to all spaces */

  notes = (char *)malloc(numcands * widths[numcols-1]);
  memset(notes, ' ', numcands * widths[numcols-1]);
  
  /* Compare the candidates with the pulsar database */

  if (idata->ra_h && idata->dec_d){
    for (ii=0; ii<numcands; ii++){
      comp_psr_to_cand(props + ii, idata, notes + ii * 20, 0);
    }
  }

  /* Compare the candidates with themselves */

  compare_rzw_cands(props, numcands, notes);

  /* Print the header */
  
  for (ii=0; ii<numcols-1; ii++){
    center_string(ctrstr, titles1[ii], widths[ii]);
    fprintf(obs->workfile, "%s  ", ctrstr);
  }
  center_string(ctrstr, titles1[ii], widths[ii]);
  fprintf(obs->workfile, "%s\n", ctrstr);
  for (ii=0; ii<numcols-1; ii++){
    center_string(ctrstr, titles2[ii], widths[ii]);
    fprintf(obs->workfile, "%s  ", ctrstr);
  }
  center_string(ctrstr, titles2[ii], widths[ii]);
  fprintf(obs->workfile, "%s\n", ctrstr);
  for (ii=0; ii<numcols-1; ii++){
    memset(tmpstr, '-', widths[ii]);
    tmpstr[widths[ii]] = '\0';
    fprintf(obs->workfile, "%s--", tmpstr);
  }
  memset(tmpstr, '-', widths[ii]);
  tmpstr[widths[ii]] = '\0';
  fprintf(obs->workfile, "%s\n", tmpstr);
  
  /* Print the fundamentals */
  
  for (ii=0; ii<numcands; ii++){
    cand = (accelcand *)(listptr->data);
    calc_rzwerrs(props+ii, obs->T, &errs);
    sprintf(tmpstr, "%-4d", ii+1);
    center_string(ctrstr, tmpstr, widths[0]);
    fprintf(obs->workfile, "%s  ", ctrstr);
    sprintf(tmpstr, "%.2f", cand->sigma);
    center_string(ctrstr, tmpstr, widths[1]);
    fprintf(obs->workfile, "%s  ", ctrstr);
    sprintf(tmpstr, "%.2f", cand->power);
    center_string(ctrstr, tmpstr, widths[2]);
    fprintf(obs->workfile, "%s  ", ctrstr);
    sprintf(tmpstr, "%d", cand->numharm);
    center_string(ctrstr, tmpstr, widths[3]);
    fprintf(obs->workfile, "%s  ", ctrstr);
    write_val_with_err(obs->workfile, errs.p*1000.0, errs.perr*1000.0, 
		       errors[4], widths[4]);
    write_val_with_err(obs->workfile, errs.f, errs.ferr, 
		       errors[5], widths[5]);
    write_val_with_err(obs->workfile, props[ii].r, props[ii].rerr, 
		       errors[6], widths[6]);
    write_val_with_err(obs->workfile, errs.fd, errs.fderr, 
		       errors[7], widths[7]);
    write_val_with_err(obs->workfile, props[ii].z, props[ii].zerr, 
		       errors[8], widths[8]);
    accel = props[ii].z * SOL / (obs->T * obs->T * errs.f);
    accelerr = props[ii].zerr * SOL / (obs->T * obs->T * errs.f);
    write_val_with_err(obs->workfile, accel, accelerr, 
		       errors[9], widths[9]);
    fprintf(obs->workfile, "  %.20s\n", notes + ii * 20);
    fflush(obs->workfile);
    listptr = listptr->next;
  }
  fprintf(obs->workfile, "\n\n");
  free(notes);
}


void output_harmonics(GSList *list, accelobs *obs)
{
  int ii, jj, numcols=12, numcands;
  int widths[12]={4, 4, 5, 15, 15, 15, 12, 11, 8, 10, 10, 10};
  int errors[12]={0, 0, 0,  2,  2,  2,  0,  2, 0,  2,  2,  2};
  char tmpstr[30], ctrstr[30], command[200];
  accelcand *cand;
  GSList *listptr;
  fourierprops props;
  rzwerrs errs;
  static char *titles1[]={"", "", "", "Power /", "Raw", 
			  "FFT 'r'", "Pred 'r'", "FFT 'z'", "Pred 'z'", 
			  "Phase", "Centroid", "Purity"};
  static char *titles2[]={"Cand", "Harm", "Sigma", "Loc Pow", "Power", 
			  "(bin)", "(bin)", "(bins)", "(bins)", 
			  "(rad)", "(0-1)", "<p> = 1"};

  numcands = g_slist_length(list);
  listptr = list;

  /* Print the header */
  
  for (ii=0; ii<numcols-1; ii++){
    center_string(ctrstr, titles1[ii], widths[ii]);
    fprintf(obs->workfile, "%s  ", ctrstr);
  }
  center_string(ctrstr, titles1[ii], widths[ii]);
  fprintf(obs->workfile, "%s\n", ctrstr);
  for (ii=0; ii<numcols-1; ii++){
    center_string(ctrstr, titles2[ii], widths[ii]);
    fprintf(obs->workfile, "%s  ", ctrstr);
  }
  center_string(ctrstr, titles2[ii], widths[ii]);
  fprintf(obs->workfile, "%s\n", ctrstr);
  for (ii=0; ii<numcols-1; ii++){
    memset(tmpstr, '-', widths[ii]);
    tmpstr[widths[ii]] = '\0';
    fprintf(obs->workfile, "%s--", tmpstr);
  }
  memset(tmpstr, '-', widths[ii]);
  tmpstr[widths[ii]] = '\0';
  fprintf(obs->workfile, "%s\n", tmpstr);
  
  /* Print the fundamentals */
  
  for (ii=0; ii<numcands; ii++){
    cand = (accelcand *)(listptr->data);
    for (jj=0; jj<cand->numharm; jj++){
      calc_props(cand->derivs[jj], cand->hirs[jj], 
		 cand->hizs[jj], 0.0, &props);
      calc_rzwerrs(&props, obs->T, &errs);
      if (jj==0) sprintf(tmpstr, " %-3d", ii+1);
      else sprintf(tmpstr, "    ");
      center_string(ctrstr, tmpstr, widths[0]);
      fprintf(obs->workfile, "%s  ", ctrstr);
      sprintf(tmpstr, "%-4d", jj+1);
      center_string(ctrstr, tmpstr, widths[1]);
      fprintf(obs->workfile, "%s  ", ctrstr);
      sprintf(tmpstr, "%.2f", candidate_sigma(props.pow, 1, 1));
      center_string(ctrstr, tmpstr, widths[2]);
      fprintf(obs->workfile, "%s  ", ctrstr);
      write_val_with_err(obs->workfile, props.pow, props.powerr, 
			 errors[3], widths[3]);
      write_val_with_err(obs->workfile, props.rawpow, 
			 sqrt(2.0*props.rawpow), 
			 errors[4], widths[4]);
      write_val_with_err(obs->workfile, props.r, props.rerr, 
			 errors[5], widths[5]);
      sprintf(tmpstr, "%.2f", cand->r * (jj+1));
      center_string(ctrstr, tmpstr, widths[6]);
      fprintf(obs->workfile, "%s  ", ctrstr);
      write_val_with_err(obs->workfile, props.z, props.zerr, 
			 errors[7], widths[7]);
      sprintf(tmpstr, "%.2f", cand->z * (jj+1));
      center_string(ctrstr, tmpstr, widths[8]);
      fprintf(obs->workfile, "%s  ", ctrstr);
      write_val_with_err(obs->workfile, props.phs, props.phserr, 
			 errors[9], widths[9]);
      write_val_with_err(obs->workfile, props.cen, props.cenerr, 
			 errors[10], widths[10]);
      write_val_with_err(obs->workfile, props.pur, props.purerr, 
			 errors[11], widths[11]);
      fprintf(obs->workfile, "\n");
      fflush(obs->workfile);
    }     
    listptr = listptr->next;
  }
  fprintf(obs->workfile, "\n\n");
  fclose(obs->workfile);
  sprintf(command, "cat %s.inf >> %s", obs->rootfilenm, obs->accelnm);
  system(command);
}


void print_accelcand(gpointer data, gpointer user_data)
{
  accelcand *obj=(accelcand *)data;

  user_data = NULL;
  printf("sigma: %-7.4f  pow: %-7.2f  harm: %-2d  r: %-14.4f  z: %-10.4f\n", 
	 obj->sigma, obj->power, obj->numharm, obj->r, obj->z); 
}


ffdotpows *subharm_ffdot_plane(int numharm, int harmnum,
			       double fullrlo, double fullrhi, 
			       subharminfo *shi, accelobs *obs)
{
  int ii, lobin, hibin, numdata, nrs, fftlen, binoffset;
  float powargr, powargi;
  double drlo, drhi, norm;
  ffdotpows *ffdot;
  fcomplex *data, **result;
  presto_datainf datainf;

  ffdot = (ffdotpows *)malloc(sizeof(ffdotpows));

  /* Calculate and get the required amplitudes */

  drlo = calc_required_r(numharm, harmnum, fullrlo);
  drhi = calc_required_r(numharm, harmnum, fullrhi);
  ffdot->rlo = (int) floor(drlo);
  ffdot->zlo = calc_required_z(numharm, harmnum, obs->zlo);
  ffdot->numrs = (int) ((drhi - drlo) * ACCEL_RDR + DBLCORRECT) + 1;
  if (ffdot->numrs % ACCEL_RDR)
    ffdot->numrs = (ffdot->numrs / ACCEL_RDR + 1) * ACCEL_RDR;
  ffdot->numzs = shi->numkern;
  binoffset = shi->kern[0].kern_half_width;
  fftlen = shi->kern[0].fftlen;
  lobin = ffdot->rlo - binoffset;
  hibin = (int) ceil(drhi) + binoffset;
  numdata = hibin - lobin + 1;
  data = read_fcomplex_file(obs->fftfile, lobin - obs->lobin, numdata);

  /* Determine the mean local power level (via median) */

  {
    float *powers;

    powers = gen_fvect(numdata);
    for (ii=0; ii<numdata; ii++) 
      powers[ii] = POWER(data[ii].r, data[ii].i);
    norm = 1.0 / (1.442695 * median(powers, numdata));
    free(powers);
  }

  /* Perform the correlations */

  result = gen_cmatrix(ffdot->numzs, ffdot->numrs);
  datainf = RAW;
  for (ii=0; ii<ffdot->numzs; ii++){
    nrs = corr_complex(data, numdata, datainf, \
		       shi->kern[ii].data, fftlen, FFT, \
		       result[ii], ffdot->numrs, binoffset, \
		       ACCEL_NUMBETWEEN, binoffset, CORR);
    datainf = SAME;
  }
  free(data);

  /* Convert the amplitudes to normalized powers */

  ffdot->powers = gen_fmatrix(ffdot->numzs, ffdot->numrs);
  for (ii=0; ii<(ffdot->numzs*ffdot->numrs); ii++) 
    ffdot->powers[0][ii] = POWER(result[0][ii].r, result[0][ii].i) * norm;
  free(result[0]);
  free(result);
  return ffdot;
}


ffdotpows *copy_ffdotpows(ffdotpows *orig)
{
  int ii;
  ffdotpows *copy;

  copy = (ffdotpows *)malloc(sizeof(ffdotpows));
  copy->numrs = orig->numrs;
  copy->numzs = orig->numzs;
  copy->rlo = orig->rlo;
  copy->zlo = orig->zlo;
  copy->powers = gen_fmatrix(orig->numzs, orig->numrs);
  for (ii=0; ii<(orig->numzs * orig->numrs); ii++)
    copy->powers[0][ii] = orig->powers[0][ii];
  return copy;
}

void free_ffdotpows(ffdotpows *ffd)
{
  free(ffd->powers[0]);
  free(ffd->powers);
  free(ffd);
}


void add_ffdotpows(ffdotpows *fundamental,
		   ffdotpows *subharmonic, 
		   int numharm, int harmnum)
{
  int ii, jj, zz, rind, zind, subz, lastrind=-1;
  double rr, subr;
  float lastpow=0;
  
  for (ii=0; ii<fundamental->numzs; ii++){
    zz = fundamental->zlo + ii * ACCEL_DZ; 
    subz = calc_required_z(numharm, harmnum, zz);
    zind = index_from_z(subz, subharmonic->zlo);
    for (jj=0; jj<fundamental->numrs; jj++){
      rr = fundamental->rlo + jj * ACCEL_DR; 
      subr = calc_required_r(numharm, harmnum, rr);
      rind = index_from_r(subr, subharmonic->rlo);
      if (rind!=lastrind)
	lastpow = subharmonic->powers[zind][rind];
      fundamental->powers[ii][jj] += lastpow;
      lastrind = rind;
    }
  }
}


GSList *search_ffdotpows(ffdotpows *ffdot, int numharm, 
			 accelobs *obs, GSList *cands)
{
  int ii, jj;
  float powcut;
  long long numindep;
  
  powcut = obs->powcut[numharm-1];
  numindep = obs->numindep[numharm-1];

  for (ii=0; ii<ffdot->numzs; ii++){
    for (jj=0; jj<ffdot->numrs; jj++){
      if (ffdot->powers[ii][jj] > powcut){
	float pow, sig;
	double rr, zz;

	pow = ffdot->powers[ii][jj];
	sig = candidate_sigma(pow, numharm, numindep);
	rr = (ffdot->rlo + jj * ACCEL_DR) / numharm;
	zz = (ffdot->zlo + ii * ACCEL_DZ) / numharm;
	cands = insert_new_accelcand(cands, pow, sig, numharm, rr, zz);
	fprintf(obs->workfile,
		"%-7.2f  %-7.4f  %-2d  %-14.4f  %-14.9f  %-10.4f\n", 
		pow, sig, numharm, rr, rr / obs->T, zz);
      }
    }
  }
  return cands;
}


void create_accelobs(accelobs *obs, infodata *idata, Cmdline *cmd)
{
  int ii, rootlen;

  {
    int hassuffix=0;
    char *suffix;
    
    hassuffix = split_root_suffix(cmd->argv[0], 
				  &(obs->rootfilenm), &suffix);
    if (hassuffix){
      if (strcmp(suffix, "fft")!=0){
	printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
	       cmd->argv[0]);
	free(suffix);
	exit(0);
      }
      free(suffix);
    } else {
      printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
	     cmd->argv[0]);
      exit(0);
    }
  }
  
  /* Read the info file */

  readinf(idata, obs->rootfilenm);
  if (idata->object) {
    printf("Analyzing %s data from '%s'.\n\n", 
	   remove_whitespace(idata->object), cmd->argv[0]);
  } else {
    printf("Analyzing data from '%s'.\n\n", cmd->argv[0]);
  }

  /* Determine the output filenames */

  rootlen = strlen(obs->rootfilenm)+25;
  obs->candnm = (char *)calloc(rootlen, 1);
  obs->accelnm = (char *)calloc(rootlen, 1);
  obs->workfilenm = (char *)calloc(rootlen, 1);
  sprintf(obs->candnm, "%s_ACCEL_%d.cand", 
	  obs->rootfilenm, cmd->zmax);
  sprintf(obs->accelnm, "%s_ACCEL_%d", 
	  obs->rootfilenm, cmd->zmax);
  sprintf(obs->workfilenm, "%s_ACCEL_%d.txtcand", 
	  obs->rootfilenm, cmd->zmax);

  /* Determine the other parameters */

  if (cmd->zmax % ACCEL_DZ)
    cmd->zmax = (cmd->zmax / ACCEL_DZ + 1) * ACCEL_DZ;
  obs->fftfile = chkfopen(cmd->argv[0], "rb");
  obs->workfile = chkfopen(obs->workfilenm, "w");
  obs->N = (long long) idata->N;
  obs->nph = get_numphotons(obs->fftfile);
  obs->numbins = chkfilelen(obs->fftfile, sizeof(fcomplex));
  obs->lobin = cmd->lobin;
  if (obs->lobin > 0){
    obs->nph = 1.0;
    if (cmd->lobin > obs->numbins - 1) {
      printf("\n'lobin' is greater than the total number of\n");
      printf("   frequencies in the data set.  Exiting.\n\n");
      exit(1);
    }
  }
  obs->numharm = cmd->numharm;
  obs->dz = ACCEL_DZ;
  obs->numz = cmd->zmax * 2 + 1;
  obs->numbetween = ACCEL_NUMBETWEEN;
  obs->dt = idata->dt;
  obs->T = idata->dt * idata->N;
  if (cmd->floP){
    obs->rlo = floor(cmd->flo * obs->T);
    if (obs->rlo < obs->lobin) 
      obs->rlo = obs->lobin;
    if (obs->rlo > obs->numbins - 1) {
      printf("\nLow frequency to search 'flo' is greater than\n");
      printf("   the highest available frequency.  Exiting.\n\n");
      exit(1);
    }
  } else {
    obs->rlo = 1.0;
    if (obs->rlo < obs->lobin) 
      obs->rlo = obs->lobin;
    if (obs->rlo > obs->numbins - 1) {
      printf("\nLow frequency to search 'rlo' is greater than\n");
      printf("   the available number of points.  Exiting.\n\n");
      exit(1);
    }
  }
  obs->highestbin = obs->numbins - 1;
  if (cmd->fhiP){
    obs->highestbin = ceil(cmd->fhi * obs->T);
    if (obs->highestbin > obs->numbins - 1) 
      obs->highestbin = obs->numbins - 1;
    obs->rhi = obs->highestbin;
    if (obs->highestbin < obs->rlo){
      printf("\nHigh frequency to search 'fhi' is less than\n");
      printf("   the lowest frequency to search 'flo'.  Exiting.\n\n");
      exit(1);
    }
  } else if (cmd->rhiP){
    obs->highestbin = cmd->rhi;
    if (obs->highestbin > obs->numbins - 1) 
      obs->highestbin = obs->numbins - 1;
    obs->rhi = obs->highestbin;
    if (obs->highestbin < obs->rlo){
      printf("\nHigh frequency to search 'rhi' is less than\n");
      printf("   the lowest frequency to search 'rlo'.  Exiting.\n\n");
      exit(1);
    }
  }
  obs->dr = ACCEL_DR;
  obs->zhi = cmd->zmax;
  obs->zlo = -cmd->zmax;
  obs->sigma = cmd->sigma;
  obs->powcut = (float *)malloc(obs->numharm * sizeof(float));
  obs->numindep = (long long *)malloc(obs->numharm * sizeof(long long));
  for (ii=1; ii<=obs->numharm; ii++){
    if (obs->numz==1)
      obs->numindep[ii-1] = obs->rhi - obs->rlo;
    else
      /* The numz+1 takes care of the small amount of  */
      /* search we get above zmax and below zmin.      */
      obs->numindep[ii-1] = (obs->rhi - obs->rlo) * (obs->numz + 1) * 
	(obs->dz / 6.95) / ii;
    obs->powcut[ii-1] = power_for_sigma(obs->sigma, ii, 
					obs->numindep[ii-1]);
  }
  obs->numzap = 0;
  /*
  if (cmd->zapfileP)
    obs->numzap = get_birdies(cmd->zapfile, obs->T, obs->baryv, 
			      &(obs->lobins), &(obs->hibins));
  else
    obs->numzap = 0;
  */
}


void free_accelobs(accelobs *obs)
{
  fclose(obs->fftfile);
  free(obs->powcut);
  free(obs->numindep);
  free(obs->rootfilenm);
  free(obs->candnm);
  free(obs->accelnm);
  free(obs->workfilenm);
  if (obs->numzap){
    free(obs->lobins);
    free(obs->hibins);
  }
}
