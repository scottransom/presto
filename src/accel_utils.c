#include <glib.h>
#include "accel.h"
#include "accelsearch_cmd.h"

static int calc_r_offset(int subharm)
/* Note:  This offset is in full bins towards lower freq */
{
  return subharm / 2;
}


static int calc_z_zero(int subharm)
/* The number of times the zero 'z' val is */
/* repeated for a subharmonic.             */
{
  return (subharm % 2) ? subharm : subharm - 1;
}


static int calc_z_for_subharm(int subharm, int zbase)
/* Calculate what 'z' you need for subharmonic 'subharm' */
/* if the 'z' at the fundamental harmonic is 'zbase'.    */
{
  int numnonzero, mod;

  mod = (zbase < 0) ? -ACCEL_DZ : ACCEL_DZ;
  if (zbase % 2)
    printf("zbase (%d) is not even in calc_z_for_subharm()\n", zbase);
  if (subharm==1 || zbase==0) 
    return zbase;
  else {
    zbase = abs(zbase);
    numnonzero = (zbase / 2) - ((calc_z_zero(subharm) - 1) / 2);
    if (numnonzero <= 0)
      return 0;
    else
      return (numnonzero % subharm) ? 
	mod * (numnonzero / subharm + 1) : 
	mod * (numnonzero / subharm);
  }
}


static int calc_fftlen(int subharm, int max_zbase)
/* The fft length needed to properly process a subharmonic */
{
  int bins_needed, end_effects;

  bins_needed = (ACCEL_USELEN % subharm) ?
    ACCEL_USELEN / subharm + 2 : ACCEL_USELEN / subharm + 1;
  end_effects = z_resp_halfwidth(calc_z_for_subharm(subharm, max_zbase), 
				 LOWACC) * 2 * ACCEL_NUMBETWEEN;
  return next2_to_n(bins_needed + end_effects);
}


static int calc_numkern(int subharm, int max_zbase)
/* The total number of kernels required for a subharmonic */
{
  return calc_z_for_subharm(subharm, abs(max_zbase)) + 1;
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


static void init_subharminfo(int subharm, int zmax, subharminfo *shi)
/* Note:  'zmax' is the overall maximum 'z' in the search */
{
  int ii, z, fftlen;

  shi->subharm = subharm;
  shi->subharm_zmax = calc_z_for_subharm(subharm, zmax);
  shi->num_z_zero = calc_z_zero(subharm);
  shi->num_r_offset = calc_r_offset(subharm);
  shi->numkern = calc_numkern(subharm, zmax);
  shi->kern = (kernel *)malloc(shi->numkern * sizeof(kernel));
  fftlen = calc_fftlen(subharm, zmax);
  for (ii=0, z=-zmax; ii<shi->numkern; ii++, z+=ACCEL_DZ)
    init_kernel(z, fftlen, &(shi->kern[ii]));
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

  shi = (subharminfo *)malloc(numharm * sizeof(subharminfo));
  for (ii=1; ii<=numharm; ii++)
    init_subharminfo(ii, zmax, &(shi[ii-1]));
  return shi;
}


void free_subharminfo_vect(int numharm, subharminfo *shi)
{
  int ii;

  for (ii=0; ii<numharm; ii++)
    free_subharminfo(&(shi[ii]));
  free(shi);
}


void create_accelobs(FILE *infile, accelobs *obs, 
		      infodata *idata, Cmdline *cmd)
{
  int ii;

  if (cmd->zmax % 2) cmd->zmax++;
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
  obs->baryv = cmd->baryv;
  obs->powcut = (float *)malloc(obs->numharm * sizeof(float));
  obs->numindep = (long long *)malloc(obs->numharm * sizeof(long long));
  for (ii=1; ii<=obs->numharm; ii++){
    obs->numindep[ii-1] = (obs->rhi - obs->rlo) * obs->numz * 
      (obs->dz / 6.95) / ii;
    obs->powcut[ii-1] = power_for_sigma(obs->sigma, ii, 
					obs->numindep[ii-1]);
  }
  if (cmd->zapfileP)
    obs->numzap = get_birdies(cmd->zapfile, obs->T, obs->baryv, 
			      &(obs->lobins), &(obs->hibins));
  else
    obs->numzap = 0;
}


void free_accelobs(accelobs *obs)
{
  free(obs->numindep);
  free(obs->powcut);
  if (obs->numzap){
    free(obs->lobins);
    free(obs->hibins);
  }
}


accelcand *create_accelcand(float power, float sigma, 
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


int compare_accelcand(gconstpointer ca, gconstpointer cb)
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
  free((accelcand *)data);
}


