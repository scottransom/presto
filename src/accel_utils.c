#include <glib.h>
#include "accel.h"
#include "accelsearch_cmd.h"

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
  return (int) ((ACCEL_RDZ * (fabs(zfull) / numharm) 
		 * harmnum) + 0.5) * ACCEL_DZ;
}


static inline double calc_required_r(int numharm, int harmnum, double rfull)
/* Calculate the 'r' you need for subharmonic     */
/* 'harmnum' out of 'numharm' subharmonics if the */
/* 'r' at the fundamental harmonic is 'rfull'.    */
{
  return (int) ((ACCEL_RDR * (rfull / numharm) 
		 * harmnum) + 0.5) * ACCEL_DR;
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
  int ii, z, fftlen;

  shi->numharm = numharm;
  shi->harmnum = harmnum;
  shi->zmax = calc_required_z(numharm, harmnum, zmax);
  shi->numkern = shi->zmax * 2 + 1;
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

  shi = (subharminfo *)malloc(numharm-1 * sizeof(subharminfo));
  for (ii=1; ii<numharm; ii++)
    init_subharminfo(ii, zmax, &(shi[ii-1]));
  return shi;
}


void free_subharminfo_vect(int numharm, subharminfo *shi)
{
  int ii;

  for (ii=1; ii<numharm; ii++)
    free_subharminfo(&(shi[ii-1]));
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


