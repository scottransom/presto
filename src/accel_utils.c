#include "accel.h"
#include "accel_cmd.h"

void create_accel_obs(FILE *infile, accel_obs *obs, 
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

void free_accel_obs(accel_obs *obs)
{
  free(obs->numindep);
  free(obs->powcut);
  if (obs->numzap){
    free(obs->zapfreqs);
    free(obs->zapwidths);
  }
}
