#include "prepfold.h"
#include "plot2d.h"
#include "float.h" 

#define TEST_EQUAL(a, b) (fabs(a) == 0.0 ? \
(fabs((a)-(b)) <= 2 * DBL_EPSILON ? 1 : 0) : \
(fabs((a)-(b))/fabs((a)) <= 2 * DBL_EPSILON ? 1 : 0))

/* This is a hack that allows the raw profile data */
/* to be printed on STDOUT so that they can be     */
/* easily imported into Python....                 */
/* #define PRINT_PROFILES 1 */

/********************************************/
/* The following is taken from ppgplot.c by */
/* Nick Patavalis (npat@ariadne.di.uoa.gr)  */
/********************************************/

void minmax (float *v, int nsz, float *min, float *max)
{
  register float *e;
  register float mn, mx;
  
  for (mn=mx=*v, e=v+nsz; v < e; v++)
    if (*v > mx) mx = *v;
    else if (*v < mn) mn = *v;
  *min = mn;
  *max = mx;
}


void scaleprof(double *in, float *out, int n)
/* Scales an input vector so that it goes from 0.0 to 1.0 */
{
  int ii;
  double min=1e100, max=-1e100;
  
  for (ii=0; ii<n; ii++){
    if (in[ii] > max) max = in[ii];
    else if (in[ii] < min) min = in[ii];
  }
  for (ii=0; ii<n; ii++)
    out[ii] = (float) ((in[ii] - min) / max);
}


void lininterp (float min, float max, int npts, float *v)
{
  register int i;
  register float step;
  register float lev;

  step = (max-min) / (npts-1);
  lev = min;
  for (i=0; i<npts; i++) {
    v[i] = lev;
    lev += step;
  }
}


static void autocal2d(float *a, int rn, int cn,
		      float *fg, float *bg, int nlevels, float *levels,
		      float *x1, float *x2, float *y1, float *y2, 
		      float *tr)
{
  /* int i; */
  float dx1, dx2, dy1, dy2;
  
  /* autocalibrate intensity-range. */
  if (*fg == *bg) {
    minmax(a,rn*cn,bg,fg);
    /* fprintf(stderr,"Intensity range:\n  fg=%f\n  bg=%f\n",*fg,*bg); */
  }
  
  if ((nlevels >= 2) && (levels))
    lininterp(*bg, *fg, nlevels, levels);
  
  /* autocalibrate x-y range. */
  if ((*x1 == *x2) || (*y1 == *y2)) cpgqwin(&dx1,&dx2,&dy1,&dy2);
  if (*x1 == *x2) {*x1=dx1; *x2=dx2;}
  if (*y1 == *y2) {*y1=dy1; *y2=dy2;}
  /* fprintf(stderr,"Xrange: [%f, %f]\nYrange[%f, %f]\n",*x1,*x2,*y1,*y2); */
  
  /* calculate transformation vector. */
  tr[2] = tr[4] = 0.0;    
  tr[1] = (*x2 - *x1) / cn;
  tr[0] = *x1 - (tr[1] / 2);
  tr[5] = (*y2 - *y1) / rn;
  tr[3] = *y1 - (tr[5] / 2);
  
  /* fprintf(stderr,"Tansformation vector:\n"); */
  /* for (i=0; i<6; fprintf(stderr,"  tr[%d]=%f\n",i,tr[i]),i++); */
}

/********************************************/


void prepfold_plot(prepfoldinfo *search, int xwin)
/* Make the beautiful 1 page prepfold output */
{
  int ii, jj, kk, ll, mm, profindex=0, loops=1, ct;
  int totpdelay=0, pdelay, pddelay;
  double profavg=0.0, profvar=0.0;
  double N=0.0, T, dphase, pofact, *currentprof, *lastprof;
  double parttime, *pdprofs, bestp, bestpd, bestpdd;
  double perr, pderr, pdderr, *dbestprof, fdot;
  double pfold, pdfold;
  float *ftmparr1;
  foldstats currentstats, beststats;
  /* Best Fold Plot */
  float *bestprof=NULL, *phasetwo=NULL;
  /* Profiles vs Time */
  float *timeprofs=NULL, *parttimes=NULL;
  /* RedChi vs Time */
  float *timechi=NULL;
  /* Profiles vs DM */
  float *dmprofs=NULL, *phaseone=NULL;
  /* DM vs RedChi */  
  float *dmchi=NULL;
  /* Period vs RedChi */  
  float *periodchi=NULL;
  /* P-dot vs RedChi */  
  float *pdotchi=NULL;
  /* Period P-dot 2D */  
  float *ppdot2d=NULL;

  if (xwin) loops = 2;

  if (TEST_EQUAL(search->fold.pow,1.0)){ /* Barycentric periods */
    bestp = search->bary.p1;
    bestpd = search->bary.p2;
    bestpdd = search->bary.p3;
  } else {                   /* Topocentric periods */
    bestp = search->topo.p1;
    bestpd = search->topo.p2;
    bestpdd = search->topo.p3;
  }
  pfold = 1.0 / search->fold.p1;
  pdfold = switch_pfdot(search->fold.p1, search->fold.p2);

  /* Time interval of 1 profile bin */

  dphase = 1.0 / (search->fold.p1 * search->proflen);

  /* Find out how many total points were folded */

  for (ii = 0; ii < search->npart; ii++)
    N += search->stats[ii * search->nsub].numdata;

  /* Calculate the time per part and the total observation time */

  parttime = search->stats[0].numdata * search->dt;
  T = N * search->dt;
  pofact = search->fold.p1 * search->fold.p1;

  /* Allocate the non-DM specific arrays */

  bestprof = gen_fvect(2 * search->proflen);
  phasetwo = gen_freqs(2 * search->proflen, 0.0, 
		       1.0 / search->proflen);
  timeprofs = gen_fvect(2 * search->proflen * search->npart);
  parttimes = gen_freqs(search->npart + 1, 0.0, parttime);
  timechi = gen_fvect(search->npart + 1);

  /* Since the reduced chi-sqrt of folded noise = 1 */

  timechi[0] = 1.0;
  periodchi = gen_fvect(search->numperiods);
  pdotchi = gen_fvect(search->numpdots);
  ppdot2d = gen_fvect(search->numperiods * search->numpdots);
  pdprofs = gen_dvect(search->npart * search->proflen);
  currentprof = gen_dvect(search->proflen);
  lastprof = gen_dvect(search->proflen);
  for (ii = 0; ii < search->proflen; ii++)
    lastprof[ii] = 0.0;

  /* Find the delays for the best periods and p-dots */
  
  for (ii = 0; ii < search->numperiods; ii++)
    if (TEST_EQUAL(search->periods[ii], bestp)){
      totpdelay = search->pstep * (ii - (search->numperiods - 1) / 2);
      break;
    }
  
  /* Correct profiles for best DM */

  if (search->nsub > 1){
    int *dmdelays;
    double *ddprofs, *subbanddelays, hif, hifdelay;
    foldstats *ddstats;

    /* Allocate DM specific arrays */

    dmprofs = gen_fvect(search->nsub * search->proflen);
    phaseone = gen_freqs(search->proflen + 1, 0.0, 
			 1.0 / search->proflen);
    dmchi = gen_fvect(search->numdms);
    
    /* Allocate local DM specific arrays*/

    ddprofs = gen_dvect(search->npart * search->proflen);
    ddstats = (foldstats *)malloc(search->npart * sizeof(foldstats));
    dmdelays = gen_ivect(search->nsub);

    /* Doppler corrected hi freq */

    hif = doppler(search->lofreq + (search->numchan - 1.0) * 
		  search->chan_wid, search->avgvoverc);

    /* De-disperse and combine the subbands */
    
    for (ii = 0; ii < search->numdms; ii++){  /* Loop over DMs */
      hifdelay = delay_from_dm(search->dms[ii], hif);
      subbanddelays = subband_delays(search->numchan, search->nsub, 
				     search->dms[ii], search->lofreq, 
				     search->chan_wid, search->avgvoverc);
      for (jj = 0; jj < search->nsub; jj++)
	dmdelays[jj] = NEAREST_INT((subbanddelays[jj] - hifdelay) / 
				   dphase) % search->proflen;
      free(subbanddelays);

      /* Make the DM vs subband plot */

      if (TEST_EQUAL(search->dms[ii], search->bestdm)){

	for (jj = 0; jj < search->nsub; jj++){

	  /* Copy the subband parts and stats into single arrays */
	  
	  for (kk = 0; kk < search->npart; kk++){
	    ll = kk * search->nsub + jj;
	    memcpy(ddprofs + kk * search->proflen, search->rawfolds + ll *
		   search->proflen, sizeof(double) * search->proflen);
	    ddstats[kk] = search->stats[ll];
	  }	  
	  
	  /* Correct each part for the best pdot and the DM delay */
	  
	  for (kk = 0; kk < search->npart; kk++){
	    profindex = kk * search->proflen;
	    fdot = switch_pfdot(pfold, bestpd) - search->fold.p2;
	    pddelay = NEAREST_INT(fdot2phasedelay(fdot, parttimes[kk]) * 
				  search->proflen);
	    shift_prof(ddprofs + profindex, search->proflen, 
		       pddelay - dmdelays[jj], pdprofs + profindex);
	  }
	  
	  /* Correct each part for the best period and sum */
	  
	  combine_profs(pdprofs, ddstats, search->npart, search->proflen, 
			totpdelay, currentprof, &currentstats);

	  /* Place the profile into the DM array */
	  
	  double2float(currentprof, dmprofs + jj * search->proflen, 
		       search->proflen);
	}
      }

      combine_subbands(search->rawfolds, search->stats, search->npart, 
		       search->nsub, search->proflen, dmdelays, 
		       ddprofs, ddstats);
      
      /* Perform the P-dot and Period searches */
      
      if (TEST_EQUAL(search->dms[ii], search->bestdm)){

	for (jj = 0; jj < search->numpdots; jj++){
	  
	  /* Correct each part for the current pdot */
	  
	  for (kk = 0; kk < search->npart; kk++){
	    profindex = kk * search->proflen;
	    fdot = switch_pfdot(pfold, search->pdots[jj]) - 
	      search->fold.p2;
	    pddelay = NEAREST_INT(fdot2phasedelay(fdot, parttimes[kk]) * 
				  search->proflen);
	    shift_prof(ddprofs + profindex, search->proflen, pddelay, 
		       pdprofs + profindex);
	  }
	
	  /* Search over the periods */
	  
	  for (kk = 0; kk < search->numperiods; kk++){
	    pdelay = search->pstep * (kk - (search->numperiods - 1) / 2);
	    combine_profs(pdprofs, ddstats, search->npart, search->proflen, 
			  pdelay, currentprof, &currentstats);

	    /* Add to the periodchi array */

	    if (TEST_EQUAL(search->pdots[jj], bestpd)) 
	      periodchi[kk] = currentstats.redchi;

	    /* Add to the pdotchi array */

	    if (TEST_EQUAL(search->periods[kk], bestp)) 
	      pdotchi[jj] = currentstats.redchi;

	    /* Add to the ppdot2d array */

	    ppdot2d[jj * search->numperiods + kk] = currentstats.redchi;

	    /* Generate the time based arrays */

	    if (TEST_EQUAL(search->periods[kk], bestp) && 
		TEST_EQUAL(search->pdots[jj], bestpd)){
	      int wrap;

	      /* The Best Prof */

	      double2float(currentprof, bestprof, search->proflen);
	      double2float(currentprof, bestprof + search->proflen, 
			   search->proflen);

	      /* Add this point to dmchi */

	      dmchi[ii] = currentstats.redchi;

	      /* Copy these statistics */

	      beststats = currentstats;

	      /* The profs at each of the npart times */

	      for (ll = 0; ll < search->npart; ll++){
		profindex = ll * search->proflen;
		wrap = (NEAREST_INT((double) (ll * pdelay) / 
				    ((double) search->npart)) % 
			search->proflen);
		shift_prof(pdprofs + profindex, search->proflen, wrap, 
			   currentprof);
		scaleprof(currentprof, timeprofs + 2 * profindex, 
			  search->proflen);
		memcpy(timeprofs + 2 * profindex + search->proflen,
		       timeprofs + 2 * profindex, 
		       search->proflen * sizeof(float));
		for (mm = 0; mm < search->proflen; mm++)
		  lastprof[mm] += currentprof[mm];
		profavg += ddstats[ll].prof_avg;
		profvar += ddstats[ll].prof_var;
		timechi[ll+1] = (chisqr(lastprof, search->proflen, 
					profavg, profvar) / 
				 (double) (search->proflen - 1.0));
	      }
	    }
	  }
	}
	
      /* Only check the best P and P-dot */

      }
      {

	/* Correct each part for the best pdot */
	  
	for (kk = 0; kk < search->npart; kk++){
	  profindex = kk * search->proflen;
	  fdot = switch_pfdot(pfold, bestpd) - search->fold.p2;
	  pddelay = NEAREST_INT(fdot2phasedelay(fdot, parttimes[kk]) * 
				search->proflen);
	  shift_prof(ddprofs + profindex, search->proflen, pddelay, 
		     pdprofs + profindex);
	}
	
	/* Correct each part for the current pdot */

	combine_profs(pdprofs, ddstats, search->npart, search->proflen, 
		      totpdelay, currentprof, &currentstats);
	dmchi[ii] = currentstats.redchi;
      }
    }
    free(ddprofs);
    free(ddstats);
    free(dmdelays);

  /* No DM corrections */

  } else {

    for (jj = 0; jj < search->numpdots; jj++){
      
      /* Correct each part for the current pdot */
      
      for (kk = 0; kk < search->npart; kk++){
	profindex = kk * search->proflen;
	fdot = switch_pfdot(pfold, search->pdots[jj]) - 
	  search->fold.p2;
	pddelay = NEAREST_INT(fdot2phasedelay(fdot, parttimes[kk]) * 
			      search->proflen);
	shift_prof(search->rawfolds + profindex, search->proflen, 
		   pddelay, pdprofs + profindex);
      }

      /* Search over the periods */
      
      for (kk = 0; kk < search->numperiods; kk++){
	pdelay = search->pstep * (kk - (search->numperiods - 1) / 2);
	combine_profs(pdprofs, search->stats, search->npart, 
		      search->proflen, pdelay, currentprof, 
		      &currentstats);
	
	/* Add to the periodchi array */
	
	if (TEST_EQUAL(search->pdots[jj], bestpd)) 
	  periodchi[kk] = currentstats.redchi;
	
	/* Add to the pdotchi array */
	
	if (TEST_EQUAL(search->periods[kk], bestp))
	  pdotchi[jj] = currentstats.redchi;
	
	/* Add to the ppdot2d array */
	
	ppdot2d[jj * search->numperiods + kk] = currentstats.redchi;
	
	/* Generate the time based arrays */
	
	if (TEST_EQUAL(search->periods[kk], bestp) && 
	    TEST_EQUAL(search->pdots[jj], bestpd)){
	  int wrap;
	  
	  /* The Best Prof */

	  double2float(currentprof, bestprof, search->proflen);
	  double2float(currentprof, bestprof + search->proflen, 
		       search->proflen);
	  
	  /* Copy these statistics */

	  beststats = currentstats;

	  /* The profs at each of the npart times */
	  
	  for (ll = 0; ll < search->npart; ll++){
	    profindex = ll * search->proflen;
	    wrap = (NEAREST_INT((double) (ll * pdelay) / 
				((double) search->npart)) % 
		    search->proflen);
	    shift_prof(pdprofs + profindex, search->proflen, wrap, 
		       currentprof);
	    double2float(currentprof, timeprofs + 2 * profindex, 
			 search->proflen);
	    double2float(currentprof, timeprofs + 2 * profindex + 
			 search->proflen, search->proflen);
	    for (mm = 0; mm < search->proflen; mm++)
	      lastprof[mm] += currentprof[mm];
	    profavg += search->stats[ll].prof_avg;
	    profvar += search->stats[ll].prof_var;
	    timechi[ll+1] = (chisqr(lastprof, search->proflen, 
				    profavg, profvar) / 
			     (double) (search->proflen - 1.0));
	  }
	}
      }
    }
  }

  /* Copy our best profile */

  dbestprof = gen_dvect(search->proflen);
  for (ii = 0; ii < search->proflen; ii++)
    dbestprof[ii] = bestprof[ii];
  
  /* Calculate the errors in the pulsation quantities */

  if (search->tepoch != 0.0 ||
      (search->tepoch == 0.0 && search->bary.p1 == 0.0))
    fold_errors(dbestprof, search->proflen, search->dt, N, 
		beststats.data_var, search->topo.p1, search->topo.p2, 
		search->topo.p3, &perr, &pderr, &pdderr);
  else
    fold_errors(dbestprof, search->proflen, search->dt, N, 
		beststats.data_var, search->bary.p1, search->bary.p2, 
		search->bary.p3, &perr, &pderr, &pdderr);
  free(dbestprof);

  /*
   *  Now plot the results
   */

  for (ct = 0; ct < loops; ct++){
    float min, max, over;

    /*Set the PGPLOT device to an X-Window */
    
    if (ct==1)
      strcpy(search->pgdev, "/XWIN");

    /* Open and prep our device */

    cpgopen(search->pgdev);
    cpgpap(10.25, 8.5/11.0);
    cpgpage();
    cpgiden();
    cpgsch(0.8);
    
    /* Time versus phase */

#ifdef PRINT_PROFILES
    {
      float min, max;
      FILE *proffile;
/*       printf("\n\ntemplate = Numeric.asarray(["); */
/*       for (ii = 0; ii < search->proflen; ii++){ */
/* 	if (ii == search->proflen - 1) */
/* 	  printf("%.2f])\n\n", bestprof[ii]); */
/* 	else */
/* 	  printf("%.2f, ", bestprof[ii]); */
/*       } */
      printf("\n\navgs = Numeric.asarray([");
      for (ii = 0; ii < search->npart; ii++){
	if (ii == search->npart - 1)
	  printf("%.2f])\n\n", search->stats[ii].prof_avg);
	else
	  printf("%.2f, ", search->stats[ii].prof_avg);
      }
      printf("\n\nstdevs = Numeric.asarray([");
      for (ii = 0; ii < search->npart; ii++){
	if (ii == search->npart - 1)
	  printf("%.2f])\n\n", sqrt(search->stats[ii].prof_var));
	else
	  printf("%.2f, ", sqrt(search->stats[ii].prof_var));
      }
      printf("\n\nmaxvals = Numeric.asarray([");
      for (ii = 0; ii < search->npart; ii++){
	minmax(timeprofs + ii * 2 * search->proflen, 
	       search->proflen, &min, &max);	
	if (ii == search->npart - 1)
	  printf("%.2f])\n\n", max);
	else
	  printf("%.2f, ", max);
      }
/*       proffile = chkfopen("NGC6544_profs.bin", "wb"); */
/*       for (ii = 0; ii < search->npart; ii++) */
/* 	chkfwrite(timeprofs + ii * 2 * search->proflen, sizeof(float), */
/* 		  search->proflen, proffile); */
/*       fclose(proffile); */
    }
#endif

    /* Scale each of the profiles in time so that the overall */
    /* scaling is not messed up.                              */

    if (ct==0) {
      int index;
      float min, max;
      for (ii = 0; ii < search->npart; ii++){
	minmax(timeprofs + ii * 2 * search->proflen, 
	       search->proflen, &min, &max);	
	for (jj = 0; jj < 2 * search->proflen; jj++){
	  index = ii * 2 * search->proflen + jj;
	  timeprofs[index] = (timeprofs[index] - min) / max;
	}
      }
    }

    cpgsvp (0.06, 0.27, 0.09, 0.68);
    cpgswin(0.0, 1.999, 0.0, T);
    {
      int mincol, maxcol, numcol, nr, nc;
      float l[2] = {0.0, 1.0};
      float r[2] = {1.0, 0.0};
      float g[2] = {1.0, 0.0};
      float b[2] = {1.0, 0.0};
      float fg = 0.0, bg = 0.0, tr[6], *levels;
      float x1 = 0.0, y1 = 0.0, x2 = 1.999, y2 = T;

      nr = search->npart;
      nc = 2 * search->proflen;
      cpgqcol(&mincol, &maxcol);
      mincol += 2;
      cpgscir(mincol, maxcol);
      numcol = maxcol - mincol + 1;
      levels = gen_fvect(numcol);
      cpgctab(l, r, g, b, numcol, 1.0, 0.5);
      autocal2d(timeprofs, nr, nc, &fg, &bg, numcol,
		levels, &x1, &x2, &y1, &y2, tr);
      cpgimag(timeprofs, nc, nr, 0+1, nc, 0+1, nr, bg, fg, tr);
      free(levels);
    }
    cpgbox ("BCNST", 0.0, 0, "BCNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Phase");
    cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");

    /*  Time versus Reduced chisqr */

    find_min_max_arr(search->npart+1, timechi, &min, &max);
    if (search->nsub > 1)
      cpgsvp (0.27, 0.36, 0.09, 0.68);
    else
      cpgsvp (0.27, 0.39, 0.09, 0.68);
    cpgswin(1.1 * max, 0.0, 0.0, T);
    cpgbox ("BCNST", 0.0, 0, "BST", 0.0, 0);
    cpgbox ("BCNST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Reduced \\gx\\u2\\d");
    cpgline(search->npart+1, timechi, parttimes);
    cpgswin(1.1 * max, 0.0, search->startT-0.0001, search->endT+0.0001);
    if (search->nsub > 1)
      cpgsch(0.7);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    cpgmtxt("R", 2.3, 0.5, 0.5, "Fraction of Observation");
    if (search->nsub > 1)
      cpgsch(0.8);

    /* Combined best profile */

    {
      float x[2] = {-0.2, 2.0}, avg[2];
      float errx = -0.1, erry = profavg, errlen;

      cpgsvp (0.06, 0.27, 0.68, 0.94);
      cpgswin(0.0, 1.999, 0.0, 1.0);
      cpgbox ("BST", 0.0, 0, "", 0.0, 0);
      cpgsvp (0.039, 0.27, 0.68, 0.94);
      find_min_max_arr(2 * search->proflen, bestprof, &min, &max);
      over = 0.1 * (max - min);
      cpgswin(-0.2, 2.0, min - over, max + over);
      cpgmtxt("T", 1.0, 0.5, 0.5, "2 Pulses of Best Profile");
      cpgline(2 * search->proflen, phasetwo, bestprof);
      cpgsls(4);
      avg[0] = avg[1] = profavg;
      cpgline(2, x, avg);
      cpgsls(1);
      errlen = sqrt(profvar);
      cpgerrb(6, 1, &errx, &erry, &errlen, 2);
      cpgpt(1, &errx, &erry, 5);
    }

    if (search->nsub > 1){

      /* DM vs reduced chisqr */

      cpgsvp (0.44, 0.66, 0.09, 0.22);
      find_min_max_arr(search->numdms, dmchi, &min, &max);
      cpgswin(search->dms[0], search->dms[search->numdms-1], 
	      0.0, 1.1 * max);
      cpgsch(0.7);
      cpgbox ("BCNST", 0.0, 0, "BCNST", 0.0, 0);
      cpgmtxt("L", 2.0, 0.5, 0.5, "Reduced \\gx\\u2\\d");
      cpgsch(0.8);
      cpgmtxt("B", 2.6, 0.5, 0.5, "DM");
      ftmparr1 = gen_fvect(search->numdms);
      double2float(search->dms, ftmparr1, search->numdms);
      cpgline(search->numdms, ftmparr1, dmchi);
      free(ftmparr1);

      /* Plots for each subband */

      {
	int chanpersb;
	double lofreq, hifreq, losubfreq, hisubfreq;
	float *tmpprof, dsubf, foffset, fnorm;

	tmpprof = gen_fvect(search->proflen + 1);
	chanpersb = search->numchan / search->nsub; 
	dsubf = chanpersb * search->chan_wid;
	lofreq = search->lofreq + dsubf - search->chan_wid;
	hifreq = search->lofreq + search->nsub * dsubf - search->chan_wid;
	losubfreq = doppler(lofreq, search->avgvoverc);
	hisubfreq = doppler(hifreq, search->avgvoverc);
	cpgsvp (0.44, 0.66, 0.3, 0.68);
	cpgswin(0.0-0.01, 1.0+0.01, 0.001, search->nsub + 0.999);
	cpgsch(0.7);
	cpgbox("BCNST", 0.2, 2, "BNST", 0.0, 0);
	cpgmtxt("L", 2.0, 0.5, 0.5, "Sub-band");
	cpgswin(0.0-0.01, 1.0+0.01, losubfreq - dsubf, hisubfreq + dsubf);
	cpgbox("", 0.2, 2, "CMST", 0.0, 0);
	cpgmtxt("R", 2.3, 0.5, 0.5, "Frequency (MHz)");
	cpgsch(0.8);
	cpgmtxt("B", 2.5, 0.5, 0.5, "Phase");
	for (ii = 0; ii < search->nsub; ii++){
	  find_min_max_arr(search->proflen, dmprofs + ii * 
			   search->proflen, &min, &max);
	  foffset = doppler(lofreq + (ii - 0.45) * dsubf, search->avgvoverc);
	  fnorm = 0.9 * dsubf / (max - min);
	  for (jj = 0; jj < search->proflen; jj++)
	    tmpprof[jj] = (dmprofs[ii * search->proflen + jj] - min) * 
	      fnorm + foffset;
	  tmpprof[search->proflen] = tmpprof[0];
	  cpgline(search->proflen + 1, phaseone, tmpprof);
	}
	free(tmpprof);
      }

    }

    {
      int mincol, maxcol, numcol, nr, nc;
      /* "Astro" (BG=Black, FG=White) (Better for X-win) */
      /*       float l[7] = {0.0, 0.167, 0.333, 0.5, 0.667, 0.833, 1.0}; */
      /*       float r[7] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; */
      /*       float g[7] = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0}; */
      /*       float b[7] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0}; */
      /* "Anti-Rainbow" (BG=White, FG=Red) (Better for printing) */
      float l[10] = {0.0, 0.035, 0.045, 0.225, 
		     0.4, 0.41, 0.6, 0.775, 0.985, 1.0};
      float r[10] = {1.0, 1.0, 0.947, 0.0, 0.0, 
		     0.0, 0.0, 1.0, 1.0, 1.0};
      float g[10] = {1.0, 0.844, 0.8, 0.0, 0.946, 
		     1.0, 1.0, 1.0, 0.0, 0.0};
      float b[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.95, 
		     0.0, 0.0, 0.0, 0.0};
      float fg = 0.0, bg = 0.0, tr[6], *levels, errlen;
      float x1l, x1h, y1l, y1h, x2l, x2h, y2l, y2h;
      char pout[100], pdout[100], fout[100], fdout[100];

      /* Plot Boundaries */ 

      /* Period / P-dot */
      x1l = (search->periods[0] - pfold) * 1000.0;
      x1h = (search->periods[search->numperiods-1] - pfold) * 1000.0;
      y1l = search->pdots[0] - pdfold;
      y1h = search->pdots[search->numpdots-1] - pdfold;
      /* Frequency / F-dot */
      x2l = 1.0 / search->periods[0] - search->fold.p1;
      x2h = 1.0 / search->periods[search->numperiods-1] - search->fold.p1;
      y2l = switch_pfdot(pfold, search->pdots[0]) - search->fold.p2;
      y2h = switch_pfdot(pfold, search->pdots[search->numperiods-1]) - 
	search->fold.p2;
      sprintf(pout, "Period - %-.8f (ms)", pfold * 1000.0);
      sprintf(fout, "Freq - %-.6f (Hz)", search->fold.p1);
      if (pdfold < 0.0)
	sprintf(pdout, "P-dot + %-.5g (s/s)", fabs(pdfold));
      else if (TEST_EQUAL(pdfold, 0.0))
	sprintf(pdout, "P-dot (s/s)");
      else
	sprintf(pdout, "P-dot - %-.5g (s/s)", pdfold);
      if (search->fold.p2 < 0.0)
	sprintf(fdout, "F-dot + %-.5g (Hz)", fabs(search->fold.p2));
      else if (TEST_EQUAL(search->fold.p2, 0.0))
	sprintf(fdout, "F-dot (Hz)");
      else
	sprintf(fdout, "F-dot - %-.5g (Hz)", search->fold.p2);

      /* Period vs reduced chisqr */

      cpgsch(0.8);
      ftmparr1 = gen_fvect(search->numperiods);
      for (ii = 0; ii < search->numperiods; ii++)
	ftmparr1[ii] = (search->periods[ii] - pfold) * 1000.0;
      find_min_max_arr(search->numperiods, periodchi, &min, &max);
      if (search->nsub > 1){
	cpgsvp (0.74, 0.94, 0.41, 0.51);
	cpgswin(x1l, x1h, 0.0, 1.1 * max);
	cpgline(search->numperiods, ftmparr1, periodchi);
	cpgsch(0.5);
	cpgbox ("BCNST", 0.0, 0, "BCMST", 0.0, 0);
	cpgsch(0.7);
	cpgmtxt("B", 2.2, 0.5, 0.5, pout);
	cpgmtxt("R", 2.4, 0.5, 0.5, "Reduced \\gx\\u2\\d");
      } else {
	cpgsvp (0.51, 0.82, 0.49, 0.63);
	cpgswin(x1l, x1h, 0.001, 1.1 * max);
	cpgline(search->numperiods, ftmparr1, periodchi);
	cpgsch(0.7);
	cpgbox ("BST", 0.0, 0, "BCMST", 0.0, 0);
	cpgswin(x2l, x2h, 0.001, 1.1 * max);
	cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
	cpgsch(0.8);
	cpgmtxt("T", 1.8, 0.5, 0.5, fout);
      }
      free(ftmparr1);

      /* P-dot vs reduced chisqr */

      ftmparr1 = gen_fvect(search->numpdots);
      for (ii = 0; ii < search->numpdots; ii++)
	ftmparr1[ii] = search->pdots[ii] - pdfold;
      find_min_max_arr(search->numpdots, pdotchi, &min, &max);
      if (search->nsub > 1){
	cpgsvp (0.74, 0.94, 0.58, 0.68);
	cpgswin(y1l, y1h, 0.0, 1.1 * max);
	cpgline(search->numpdots, ftmparr1, pdotchi);
	cpgsch(0.5);
	cpgbox ("BCNST", 0.0, 0, "BCMST", 0.0, 0);
	cpgsch(0.7);
	cpgmtxt("B", 2.2, 0.5, 0.5, pdout);
	cpgmtxt("R", 2.4, 0.5, 0.5, "Reduced \\gx\\u2\\d");
      } else {
	cpgsvp (0.82, 0.93, 0.09, 0.49);
	cpgswin(0.001, 1.1 * max, y1l, y1h);
	cpgline(search->numpdots, pdotchi, ftmparr1);
	cpgsch(0.7);
	cpgbox ("BCMST", 0.0, 0, "BST", 0.0, 0);
	cpgswin(0.001, 1.1 * max, y2l, y2h);
	cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
	cpgsch(0.8);
	cpgmtxt("T", 4.2, 0.5, 0.5, "Reduced");
	cpgmtxt("T", 2.8, 0.5, 0.5, "\\gx\\u2\\d");
	cpgmtxt("R", 2.4, 0.5, 0.5, fdout);
      }
      free(ftmparr1);

      /* P P-dot image */

      if (search->nsub > 1)
	cpgsvp (0.74, 0.94, 0.09, 0.29);
      else
	cpgsvp (0.51, 0.82, 0.09, 0.49);
      cpgswin(x1l, x1h, y1l, y1h);
      nr = search->numpdots;
      nc = search->numperiods;
      cpgqcol(&mincol, &maxcol);
      mincol += 2;
      cpgscir(mincol, maxcol);
      numcol = maxcol - mincol + 1;
      levels = gen_fvect(numcol);
      cpgctab(l, r, g, b, numcol, 1.0, 0.5);
      autocal2d(ppdot2d, nr, nc, &fg, &bg, numcol,
		levels, &x1l, &x1h, &y1l, &y1h, tr);
      cpgimag(ppdot2d, nc, nr, 0+1, nc, 0+1, nr, bg, fg, tr);
      x1l = (float) ((bestp - pfold) * 1000.0);
      y1l = (float) (bestpd - pdfold);
      /* Plot the error bars on the P-Pdot diagram */
      cpgpt(1, &x1l, &y1l, 5);
      errlen = (float) (perr * 1000.0);
      cpgerrb(5, 1, &x1l, &y1l, &errlen, 2);
      errlen = (float) (pderr);
      cpgerrb(6, 1, &x1l, &y1l, &errlen, 2);
      if (search->nsub > 1){
	cpgsch(0.5);
	cpgbox("BNST", 0.0, 0, "BNST", 0.0, 0);
	cpgsch(0.7);
	cpgmtxt("B", 2.4, 0.5, 0.5, pout);
	cpgmtxt("L", 2.0, 0.5, 0.5, pdout);
	cpgswin(x2l, x2h, y2l, y2h);
	cpgsch(0.5);
	cpgbox("CMST", 0.0, 0, "CMST", 0.0, 0);
	cpgsch(0.7);
	cpgmtxt("T", 1.8, 0.5, 0.5, fout);
	cpgmtxt("R", 2.3, 0.5, 0.5, fdout);
      } else {
	cpgsch(0.7);
	cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
	cpgsch(0.8);
	cpgmtxt("B", 2.6, 0.5, 0.5, pout);
	cpgmtxt("L", 2.1, 0.5, 0.5, pdout);
      }
      free(levels);
      cpgsch(0.8);
    }

    {
      char out[200], out2[100];

      /* Add the Data Info area */

      cpgsvp (0.27, 0.519, 0.68, 0.94);
      cpgswin(-0.1, 1.00, -0.1, 1.1);
      cpgsch(1.0);
      sprintf(out, "%-s", search->filenm);
      cpgmtxt("T", 1.0, 0.5, 0.5, out);
      cpgsch(0.7);
      sprintf(out, "Candidate:  %-s", search->candnm);
      cpgtext(0.0, 1.0, out);
      sprintf(out, "Telescope:  %-s", search->telescope);
      cpgtext(0.0, 0.9, out);
      if (TEST_EQUAL(search->tepoch, 0.0))
	sprintf(out, "Epoch\\dtopo\\u = N/A");
      else
	sprintf(out, "Epoch\\dtopo\\u = %-.11f", search->tepoch);
      cpgtext(0.0, 0.8, out);
      if (TEST_EQUAL(search->bepoch, 0.0))
	sprintf(out, "Epoch\\dbary\\u = N/A");
      else
	sprintf(out, "Epoch\\dbary\\u = %-.11f", search->bepoch);
      cpgtext(0.0, 0.7, out);
      sprintf(out, "T\\dsample\\u");
      cpgtext(0.0, 0.6, out);
      sprintf(out, "=  %f", search->dt);
      cpgtext(0.45, 0.6, out);
      sprintf(out, "Data Folded");
      cpgtext(0.0, 0.5, out);
      sprintf(out, "=  %-.0f", N);
      cpgtext(0.45, 0.5, out);
      sprintf(out, "Data Avg");
      cpgtext(0.0, 0.4, out);
      sprintf(out, "=  %.2f", beststats.data_avg);
      cpgtext(0.45, 0.4, out);
      sprintf(out, "Data StdDev"); 
      cpgtext(0.0, 0.3, out);
      sprintf(out, "=  %.2f", sqrt(beststats.data_var));
      cpgtext(0.45, 0.3, out);
      sprintf(out, "Profile Bins");
      cpgtext(0.0, 0.2, out);
      sprintf(out, "=  %d", search->proflen);
      cpgtext(0.45, 0.2, out);
      sprintf(out, "Profile Avg");
      cpgtext(0.0, 0.1, out);
      sprintf(out, "=  %.2f", beststats.prof_avg);
      cpgtext(0.45, 0.1, out);
      sprintf(out, "Profile StdDev");
      cpgtext(0.0, 0.0, out);
      sprintf(out, "=  %.2f", sqrt(beststats.prof_var));
      cpgtext(0.45, 0.0, out);

      /* Calculate the values of P and Q since we know X and DF */

      {
	int chiwhich=1, chistatus=0, goodsig=1;
	double chip=0.0, chiq=0.0, chixmeas=0.0, chidf=0.0, chitmp=0.0;
	double normz=0.0, normmean=0.0, normstdev=1.0;
	
	chidf = search->proflen - 1.0;
	chixmeas = beststats.redchi * chidf;
	cdfchi(&chiwhich, &chip, &chiq, &chixmeas, &chidf, &chistatus, &chitmp);
	if (chistatus != 0){
	  if (chistatus < 0)
	    printf("\nInput parameter %d to cdfchi() was out of range.\n", 
		   chistatus);
	  else if (chistatus == 3)
	    printf("\nP + Q do not equal 1.0 in cdfchi().\n");
	  else if (chistatus == 10)
	    printf("\nError in cdfgam().\n");
	  else printf("\nUnknown error in cdfchi().\n");
	}
	
	/* Calculate the equivalent sigma */
	
	chiwhich = 2;
	cdfnor(&chiwhich, &chip, &chiq, &normz, &normmean, &normstdev, 
	       &chistatus, &chitmp);
	if (chistatus != 0) goodsig=0;

	/* Add the Fold Info area */
      
	cpgsvp (0.519, 0.94, 0.68, 0.94);
	cpgswin(-0.05, 1.05, -0.1, 1.1);
	cpgsch(0.8);
	cpgmtxt("T", 1.0, 0.5, 0.5, "Search Information");
	cpgsch(0.7);
	cpgtext(0.0, 1.0, "   Best Fit Parameters");
	if (goodsig)
	  sprintf(out2, "(\\(0248)%.1f\\gs)", normz);
	else 
	  sprintf(out2, " ");
	sprintf(out, "Reduced \\gx\\u2\\d = %.3f   P(Noise) < %.3g   %s", 
		beststats.redchi, chiq, out2);
	cpgtext(0.0, 0.9, out);
	if (search->nsub > 1){
	  sprintf(out, "Dispersion Measure (DM) = %.3f", search->bestdm);
	  cpgtext(0.0, 0.8, out);
	}
	{
	  if (search->tepoch != 0.0){
	    cpgnice_output_2(out2, search->topo.p1*1000.0, perr*1000.0, 0);
	    sprintf(out, "P\\dtopo\\u (ms) = %s", out2);
	    cpgtext(0.0, 0.7, out);
	    cpgnice_output_2(out2, search->topo.p2, pderr, 0);
	    sprintf(out, "P'\\dtopo\\u (s/s) = %s", out2);
	    cpgtext(0.0, 0.6, out);
	    cpgnice_output_2(out2, search->topo.p3, pdderr, 0);
	    sprintf(out, "P''\\dtopo\\u (s/s\\u2\\d) = %s", out2);
	    cpgtext(0.0, 0.5, out);
	  } else {
	    cpgtext(0.0, 0.7, "P\\dtopo\\u (ms) = N/A");
	    cpgtext(0.0, 0.6, "P'\\dtopo\\u (s/s) = N/A");
	    cpgtext(0.0, 0.5, "P''\\dtopo\\u (s/s\\u2\\d) = N/A");
	  }

	  if (search->bepoch != 0.0){
	    cpgnice_output_2(out2, search->bary.p1*1000.0, perr*1000.0, 0);
	    sprintf(out, "P\\dbary\\u (ms) = %s", out2);
	    cpgtext(0.6, 0.7, out);
	    cpgnice_output_2(out2, search->bary.p2, pderr, 0);
	    sprintf(out, "P'\\dbary\\u (s/s) = %s", out2);
	    cpgtext(0.6, 0.6, out);
	    cpgnice_output_2(out2, search->bary.p3, pdderr, 0);
	    sprintf(out, "P''\\dbary\\u (s/s\\u2\\d) = %s", out2);
	    cpgtext(0.6, 0.5, out);
	  } else {
	    cpgtext(0.6, 0.7, "P\\dbary\\u (ms) = N/A");
	    cpgtext(0.6, 0.6, "P'\\dbary\\u (s/s) = N/A");
	    cpgtext(0.6, 0.5, "P''\\dbary\\u (s/s\\u2\\d) = N/A");
	  }
	}
	cpgtext(0.0, 0.3, "   Binary Parameters");
	if (TEST_EQUAL(search->orb.p, 0.0)){
	  cpgtext(0.0, 0.2, "P\\dorb\\u (s) = N/A");
	  cpgtext(0.0, 0.1, "a\\d1\\usin(i)/c (s) = N/A");
	  cpgtext(0.6, 0.2, "e = N/A");
	  cpgtext(0.6, 0.1, "\\gw (rad) = N/A");
	  cpgtext(0.0, 0.0, "T\\dperi\\u = N/A");
	} else {
	  sprintf(out, "P\\dorb\\u (s) = %f", search->orb.p);
	  cpgtext(0.0, 0.2, out);
	  sprintf(out, "a\\d1\\usin(i)/c (s) = %f", search->orb.x);
	  cpgtext(0.0, 0.1, out);
	  sprintf(out, "e = %f", search->orb.e);
	  cpgtext(0.6, 0.2, out);
	  sprintf(out, "\\gw (rad) = %f", search->orb.w);
	  cpgtext(0.6, 0.1, out);
	  sprintf(out, "T\\dperi\\u = %-.11f", search->orb.t);
	  cpgtext(0.0, 0.0, out);
	}
      }
    }
    cpgclos();
  }
  free(bestprof);
  free(phasetwo);
  free(timeprofs);
  free(parttimes);
  free(timechi);
  free(periodchi);
  free(pdotchi);
  free(ppdot2d);
  free(pdprofs);
  free(currentprof);
  free(lastprof);
  if (search->nsub > 1){
    free(dmprofs);
    free(phaseone);
    free(dmchi);
  }
}
