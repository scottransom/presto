#include "accel.h"

#if defined (__GNUC__)
#  define inline __inline__
#else
#  undef inline
#endif

static inline int calc_required_z(int numharm, int harmnum, double zfull)
/* Calculate the 'z' you need for subharmonic     */
/* 'harmnum' out of 'numharm' subharmonics if the */
/* 'z' at the fundamental harmonic is 'zfull'.    */
{
  return (int) ((ACCEL_RDZ * (zfull / numharm) * harmnum) + 0.5) * ACCEL_DZ;
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
/* has stepsize ACCEL_DR and low freq 'lor'.              */
{
  return (int) ((z - loz) * ACCEL_RDZ + DBLCORRECT);
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
  printf("  Sub-harmonic %d has %d kernels from z = %d to %d.\n", 
	 harmnum, shi->numkern, -shi->zmax, shi->zmax);
  for (ii=0; ii<shi->numkern; ii++)
    init_kernel(-shi->zmax + ii * ACCEL_DZ, 
		fftlen, &(shi->kern[ii]));
}


static void free_subharminfo(subharminfo *shi)
{
  free_kernel(shi->kern);
  free(shi->kern);
}


subharminfo *create_subharminfo_vect(int numharm, int zmax)
{
  int ii;
  subharminfo *shi;

  if (numharm==1){  /* Fundamental */
    shi = (subharminfo *)malloc(sizeof(subharminfo));
    init_subharminfo(1, 1, zmax, shi);
  } else {          /* Sub-harmonics */
    shi = (subharminfo *)malloc((numharm-1) * sizeof(subharminfo));
    for (ii=1; ii<numharm; ii++)
      init_subharminfo(numharm, ii, zmax, &(shi[ii-1]));
  }
  return shi;
}


void free_subharminfo_vect(int numharm, subharminfo *shi)
{
  int ii;

  if (numharm==1){
    free_subharminfo(shi);
  } else {
    for (ii=1; ii<numharm; ii++)
      free_subharminfo(&(shi[ii-1]));
  }
  free(shi);
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
  return obj;
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
    new_list->data = create_accelcand(power, sigma, numharm, rr, zz);
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
      free(tmp_list->data);  /* Overwrite the old candidate */
      tmp_list->data = create_accelcand(power, sigma, numharm, rr, zz);
    }
  } else {  /* This is a new candidate */
    new_list = g_slist_alloc();
    new_list->data = create_accelcand(power, sigma, numharm, rr, zz);

    if (!tmp_list->next){
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
  for (ii=1; ii<=cand->numharm; ii++){
    cand->pows[ii] = max_rz_file(obs->fftfile, 
				 cand->r*ii, cand->z*ii, \
				 &(cand->hirs[ii]), 
				 &(cand->hizs[ii]), 
				 &(cand->derivs[ii]));
    cand->power
}


void print_accelcand(gpointer data, gpointer user_data)
{
  accelcand *obj=(accelcand *)data;

  user_data = NULL;
  printf("sigma: %-7.4f  pow: %-7.2f  harm: %-2d  r: %-14.4f  z: %-10.4f\n", 
	 obj->sigma, obj->power, obj->numharm, obj->r, obj->z); 
}


void free_accelcand(gpointer data, gpointer user_data)
{
  user_data = NULL;
  free(((accelcand *)data)->pows);
  free(((accelcand *)data)->hirs);
  free(((accelcand *)data)->hizs);
  free(((accelcand *)data)->rderivs);
  free((accelcand *)data);
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
  ffdot->zlo = calc_required_z(numharm, harmnum, -obs->zlo);
  ffdot->numrs = (int) ((drhi - drlo) * ACCEL_RDR + DBLCORRECT) + 1;
  ffdot->numzs = shi->numkern;
  binoffset = shi->kern[0].kern_half_width;
  fftlen = shi->kern[0].fftlen;
  lobin = ffdot->rlo - binoffset;
  hibin = (int) ceil(drhi) + binoffset;
  numdata = hibin - lobin + 1;
  data = read_fcomplex_file(obs->fftfile, lobin, numdata);

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
  for (ii=0; ii<(orig->numzs*orig->numrs); ii++)
    copy->powers[ii] = orig->powers[ii];
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
  int ii, jj, rr, zz, rind, zind, subz, lastrind=-1;
  double subr;
  float lastpow=0;
  
  for (ii=0; ii<fundamental->numzs; ii++){
    zz = fundamental->zlo; 
    subz = calc_required_z(numharm, harmnum, zz);
    zind = index_from_z(subz, subharmonic->zlo);
    for (jj=0; jj<fundamental->numrs; jj++){
      rr = fundamental->rlo; 
      subr = calc_required_r(numharm, harmnum, rr);
      rind = index_from_r(subr, subharmonic->rlo);
      if (rind!=lastrind)
	lastpow = subharmonic->powers[zind][rind];
      fundamental->powers[ii][jj] += lastpow;
      lastrind = rind;
      rr += ACCEL_DR;
    }
    zz += ACCEL_DZ;
  }
}


void search_ffdotpows(ffdotpows *ffdot, int numharm, 
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
}


void create_accelobs(FILE *infile, FILE *workfile, accelobs *obs, 
		     infodata *idata, Cmdline *cmd)
{
  int ii;

  if (cmd->zmax % ACCEL_DZ)
    cmd->zmax = (cmd->zmax / ACCEL_DZ + 1) * ACCEL_DZ;
  obs->fftfile = infile;
  obs->workfile = workfile;
  obs->N = (long long) idata->N;
  obs->nph = get_numphotons(infile);
  obs->numbins = chkfilelen(infile, sizeof(fcomplex));
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
    cmd->rlo = cmd->flo * obs->T;
    if (cmd->rlo < cmd->lobin) 
      cmd->rlo = cmd->lobin;
    if (cmd->rlo > obs->numbins - 1) {
      printf("\nLow frequency to search 'flo' is greater than\n");
      printf("   the highest available frequency.  Exiting.\n\n");
      exit(1);
    }
  } else {
    if (cmd->rlo < cmd->lobin) 
      cmd->rlo = cmd->lobin;
    if (cmd->rlo > obs->numbins - 1) {
      printf("\nLow frequency to search 'rlo' is greater than\n");
      printf("   the available number of points.  Exiting.\n\n");
      exit(1);
    }
  }
  obs->highestbin = obs->numbins - 1;
  if (cmd->fhiP){
    obs->highestbin = cmd->fhi * obs->T;
    if (obs->highestbin > obs->numbins - 1) 
      obs->highestbin = obs->numbins - 1;
    if (obs->highestbin < cmd->rlo){
      printf("\nHigh frequency to search 'fhi' is less than\n");
      printf("   the lowest frequency to search 'flo'.  Exiting.\n\n");
      exit(1);
    }
  } else if (cmd->rhiP){
    obs->highestbin = cmd->rhi;
    if (obs->highestbin > obs->numbins - 1) 
      obs->highestbin = obs->numbins - 1;
    if (obs->highestbin < cmd->rlo){
      printf("\nHigh frequency to search 'rhi' is less than\n");
      printf("   the lowest frequency to search 'rlo'.  Exiting.\n\n");
      exit(1);
    }
  }
  obs->rlo = cmd->rlo;
  obs->rhi = cmd->rhi;
  obs->dr = ACCEL_DR;
  obs->zhi = cmd->zmax;
  obs->zlo = -cmd->zmax;
  obs->sigma = cmd->sigma;
  obs->powcut = (float *)malloc(obs->numharm * sizeof(float));
  obs->numindep = (long long *)malloc(obs->numharm * sizeof(long long));
  for (ii=1; ii<=obs->numharm; ii++){
    obs->numindep[ii-1] = (obs->rhi - obs->rlo) * obs->numz * 
      (obs->dz / 6.95) / ii;
    obs->powcut[ii-1] = power_for_sigma(obs->sigma, ii, 
					obs->numindep[ii-1]);
  }
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
  fclose(obs->workfile);
  free(obs->numindep);
  free(obs->powcut);
  if (obs->numzap){
    free(obs->lobins);
    free(obs->hibins);
  }
}


