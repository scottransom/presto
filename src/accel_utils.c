#include <glib.h>
#include "accel.h"
#include "accelsearch_cmd.h"

void create_accelobs(FILE *infile, accelobs *obs, 
		      infodata *idata, Cmdline *cmd)
{
  int ii;

  obs->N = (long long) idata->N;
  obs->numbins = chkfilelen(infile, sizeof(fcomplex));
  obs->lobin = cmd->lobin;
  obs->fftlen = cmd->fftlen;
  obs->numharm = cmd->numharm;
  obs->dz = cmd->dz;
  obs->numz = (int) ((cmd->zhi - cmd->zlo) / cmd->dz) + 1;
  obs->numbetween = cmd->numbetween;
  obs->dt = idata->dt;
  obs->T = idata->dt * idata->N;
  obs->rlo = cmd->rlo;
  obs->rhi = cmd->rhi;
  obs->dr = 1.0 / (double) cmd->numbetween;
  obs->zlo = cmd->zlo;
  obs->zhi = cmd->zhi;
  obs->sigma = cmd->sigma;
  obs->baryv = cmd->baryv;
  obs->powcut = (float *)malloc(obs->numharm * sizeof(float));
  obs->numindep = (long long *)malloc(obs->numharm * sizeof(long long));
  for (ii=1; ii<=obs->numharm; ii++){
    obs->numindep[ii-1] = (obs->rhi - obs->rlo) * numz * 
      (obs->dz / 6.95) / ii;
    obs->powcut[ii-1] = power_for_sigma(obs->sigma, ii, 
					obs->numindep[ii-1]);
  }
  if (cmd->zapfileP)
    obs->numzap = read_zapfile(cmd->zapfile, &(obs->zapfreqs), 
			       &(obs->zapwidths));
  else
    obs->numzap = 0;
}


void free_accelobs(accelobs *obs)
{
  free(obs->numindep);
  free(obs->powcut);
  if (obs->numzap){
    free(obs->zapfreqs);
    free(obs->zapwidths);
  }
}


accelcand *create_accelcand(float power, float sigma, 
			    int numharm, double r, double z);
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
  result = (a->val < b->val) - (a->val > b->val);
  if (result) 
    return result;
  else
    return (a->r > b->r) - (a->r < b->r);
}


void print_accelcand(gpointer data, gpointer user_data)
{
  accelcand *obj=(accelcand *)data;

  printf("sigma: %-7.4f  pow: %-7.2f  harm: %-2d  r: %-14.4f  z: %-10.4f\n", 
	 obj->sigma, obj->power, obj->numharm, obj->r, obj->z); 
}


void free_accelcand(gpointer data, gpointer user_data)
{
  free((accelcand *)data);
}


