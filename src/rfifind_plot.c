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

double single_power_pdf(double power, int numpowers){
  double tmp;

  tmp = exp(-power);
  return tmp * numpowers * pow(1.0-tmp, numpowers-1);
}

/*
void clip_array(float *arr, int numarr, 
*/

void minmax(float *v, int nsz, float *min, float *max)
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


void rfifind_plot(int numchan, int numint, int ptsperint, float sigma, 
		  float **dataavg, float **datastd, float **datapow, 
		  unsigned char **datamask, infodata *idata, int xwin)
/* Make the beautiful 1 page rfifind output */
{
  int ii, jj, ct, loops=1;
  float *freqs, *chans, *times, *ints;
  float *avg_chan, *std_chan, *pow_chan;
  float *avg_int, *std_int, *pow_int;
  double inttim, T, lof, hif;
  double avgavg, stdavg, powavg;
  double dataavg_med, datastd_med, datapow_med;
  double avg_chan_std, std_chan_std, pow_chan_std;
  double avg_chan_med, std_chan_med, pow_chan_med;
  double avg_int_std, std_int_std, pow_int_std;
  double avg_int_med, std_int_med, pow_int_med;
  
  inttim = ptsperint * idata->dt;
  T = ptsperint * idata->dt * numint;
  lof = idata->freq - 0.5 * idata->chan_wid;
  hif = lof + idata->freqband;
  avg_chan = gen_fvect(numchan);
  avg_int = gen_fvect(numint);
  std_chan = gen_fvect(numchan);
  std_int = gen_fvect(numint);
  pow_chan = gen_fvect(numchan);
  pow_int = gen_fvect(numint);
  chans = gen_fvect(numchan);
  freqs = gen_fvect(numchan);
  for (ii=0; ii<numchan; ii++){
    chans[ii] = ii;
    freqs[ii] = idata->freq + ii * idata->chan_wid;
  }
  ints = gen_fvect(numint);
  times = gen_fvect(numint);
  for (ii=0; ii<numint; ii++){
    ints[ii] = ii;
    freqs[ii] = idata->freq + ii * idata->chan_wid;
  }

  /* Calculate the average vectors */

  for (ii=0; ii<numint; ii++){
    avgavg = stdavg = powavg = 0.0;
    for (jj=0; jj<numchan; jj++){
      avgavg += dataavg[ii][jj];
      stdavg += datastd[ii][jj];
      powavg += datapow[ii][jj];
    }
    avg_int[ii] = avgavg / numchan;
    std_int[ii] = stdavg / numchan;
    pow_int[ii] = powavg / numchan;
  }
  for (ii=0; ii<numchan; ii++){
    avgavg = stdavg = powavg = 0.0;
    for (jj=0; jj<numint; jj++){
      avgavg += dataavg[jj][ii];
      stdavg += datastd[jj][ii];
      powavg += datapow[jj][ii];
    }
    avg_chan[ii] = avgavg / numint;
    std_chan[ii] = stdavg / numint;
    pow_chan[ii] = powavg / numint;
  }
  
  /* Calculate the statistics of the average vectors */

  avg_var(avg_int,  numint,  &avg_int_med,  &avg_int_std);
  avg_var(std_int,  numint,  &std_int_med,  &std_int_std);
  avg_var(pow_int,  numint,  &pow_int_med,  &pow_int_std);
  avg_var(avg_chan, numchan, &avg_chan_med, &avg_chan_std);
  avg_var(std_chan, numchan, &std_chan_med, &std_chan_std);
  avg_var(pow_chan, numchan, &pow_chan_med, &pow_chan_std);
  /* Turn variances into standard deviations */
  avg_int_std  = sqrt(avg_int_std);
  std_int_std  = sqrt(std_int_std);
  pow_int_std  = sqrt(pow_int_std);
  avg_chan_std = sqrt(avg_chan_std);
  std_chan_std = sqrt(std_chan_std);
  pow_chan_std = sqrt(pow_chan_std);
  /* Use medians instead of averages */
  avg_int_med  = selectkth(numint/2,  numint,  avg_int-1);
  std_int_med  = selectkth(numint/2,  numint,  std_int-1);
  pow_int_med  = selectkth(numint/2,  numint,  pow_int-1);
  avg_chan_med = selectkth(numchan/2, numchan, avg_chan-1);
  std_chan_med = selectkth(numchan/2, numchan, std_chan-1);
  pow_chan_med = selectkth(numchan/2, numchan, pow_chan-1);
  ct = numchan * numint;
  dataavg_med = selectkth(ct/2, ct, dataavg[0]-1);
  datastd_med = selectkth(ct/2, ct, datastd[0]-1);
  datapow_med = selectkth(ct/2, ct, datapow[0]-1);

  /*
   *  Now plot the results
   */

  if (xwin) loops = 2;
  for (ct=0; ct<loops; ct++){
    float min, max;
    float tt, ft, th, fh;  /* thin and fat thicknesses and heights */
    float lm, rm, tm, bm;  /* LRTB margins */
    float xarr[2], yarr[2];
    char outdev[100];

    /*Set the PGPLOT device to an X-Window */
    
    if (ct==1)
      strcpy(outdev, "/XWIN");
    else
      sprintf(outdev, "%s.ps/VCPS", idata->name);

    /* Open and prep our device */

    cpgopen(outdev);
    cpgpap(10.25, 8.5/11.0);
    cpgpage();
    cpgiden();
    cpgsch(0.8);
    
    /* Set thicknesses and margins */

    lm = rm = 0.06;
    bm = tm = 0.09;
    ft = 3.0;  /* This sets fat thickness = 3 x thin thickness */
    tt = (1.0 - lm - rm) / (6.0 + 4.0 * ft);
    ft *= tt;
    fh = 0.55;
    th = tt * 8.5/11.0;

    /* Maximum Powers */

    cpgsvp (lm, lm+ft, bm, bm+fh);
    cpgswin(0.0, numchan, 0.0, T);
    {
      int mincol, maxcol, numcol;
      float l[2] = {0.0, 1.0};
      float r[2] = {1.0, 0.0};
      float g[2] = {1.0, 0.0};
      float b[2] = {1.0, 0.0};
      float fg = 0.0, bg = 0.0, tr[6], *levels;
      float x1 = 0.0, x2 = numchan, y1 = 0.0, y2 = T;
      cpgqcol(&mincol, &maxcol);
      mincol += 2;
      cpgscir(mincol, maxcol);
      numcol = maxcol - mincol + 1;
      levels = gen_fvect(numcol);
      cpgctab(l, r, g, b, numcol, 1.0, 0.5);
      autocal2d(datapow[0], numint, numchan, &fg, &bg, numcol,
		levels, &x1, &x2, &y1, &y2, tr);
      cpgimag(datapow[0], numchan, numint, 0+1, 
	      numchan, 0+1, numint, bg, fg, tr);
      free(levels);
    }
    cpgswin(0.0, numchan, 0.0, T);
    cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
    cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");
    cpgswin(lof, hif, 0.0, numint);
    cpgbox ("CST", 0.0, 0, "CST", 0.0, 0);

    /*  Max Power versus Max Power */

    find_min_max_arr(numint, pow_int, &min, &max);
    cpgsvp (lm+ft, lm+ft+tt, bm, bm+fh);
    cpgswin(1.1*max, 0.9*min, 0.0, T);
    cpgbox ("BCNST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Power");
    cpgline(numint, pow_int, times);
    cpgscr(2, 1.0, 0.0, 0.0); cpgsci(2); /* Red */
    yarr[0] = 0.0; yarr[1] = T;
    xarr[0] = xarr[1] = pow_int_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = pow_int_med + sigma * pow_int_std;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = pow_int_med - sigma * pow_int_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgswin(1.1*max, 0.9*min, 0.0, numint);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number");

    /*  Max Power versus Channel */

    find_min_max_arr(numchan, pow_chan, &min, &max);
    cpgsvp (lm, lm+ft, bm+fh, bm+fh+th);
    cpgswin(0.0, numchan, 0.9*min, 1.1*max);
    cpgbox ("BST", 0.0, 0, "BCNST", 0.0, 0);
    cpgmtxt("L", 2.1, 0.5, 0.5, "Power");
    cpgline(numchan, chans, pow_chan);
    cpgscr(2, 1.0, 0.0, 0.0); cpgsci(2); /* Red */
    xarr[0] = 0.0; xarr[1] = numchan;
    yarr[0] = yarr[1] = pow_chan_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = pow_chan_med + sigma * pow_chan_std;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = pow_chan_med - sigma * pow_chan_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgswin(lof, hif, 0.9*min, 1.1*max);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Powers Histogram */
    {
      float *theo, *hist, *hpows, *tpows, maxhist=0.0, maxtheo=0.0;
      int numhist=30, numtheo=100, bin, numpows;
      double dtheo, dhist;

      /* Calculate the predicted distribution of max powers */

      numpows = numint * numchan;
      find_min_max_arr(numpows, datapow[0], &min, &max);
      min *= 0.9;
      max *= 1.1;
      dhist = (max-min)/numhist;
      theo = gen_fvect(numtheo);
      tpows = gen_fvect(numtheo);
      hist = gen_fvect(numhist);
      hpows = gen_fvect(numhist);
      for (ii=0; ii<numhist; ii++){
	hist[ii] = 0.0;
	hpows[ii] = min + ii * dhist;
      }
      for (ii=0; ii<numpows; ii++){
	bin = (*(datapow[0]+ii)-min)/(max-min)*numhist;
	hist[bin] += 1.0;
      }
      for (ii=0; ii<numhist; ii++)
	if (hist[ii] > maxhist) maxhist = hist[ii];
      maxhist *= 1.1;
      dtheo = (max - min) / (double)(numtheo - 1);
      for (ii=0; ii<numtheo; ii++){
	tpows[ii] = min + ii * dtheo;
	theo[ii] = single_power_pdf(tpows[ii], 
				    ptsperint / 2) * numpows;
	if (theo[ii] > maxtheo) maxtheo = theo[ii];
      }
      maxtheo *= 1.1;
      if (maxtheo > maxhist) maxhist = maxtheo;
      cpgsvp (lm, lm+ft+tt, 1.0-tm-(0.66*(tm+bm+th+fh)), 1.0-tm);
      cpgswin(min, max, 0.0, maxhist);
      cpgbox ("BCNST", 0.0, 0, "BCNST", 0.0, 0);
      cpgmtxt("L", 2.1, 0.5, 0.5, "Number");
      cpgmtxt("B", 2.6, 0.5, 0.5, "Power");
      cpgbin(numhist, hpows, hist, 0);
      cpgscr(2, 1.0, 0.0, 0.0); cpgsci(2); /* Red */
      cpgline(numtheo, tpows, theo);
      cpgsci(1); /* Default color */
      free(hist);
      free(theo);
      free(tpows);
      free(hpows);
    }

    /* Standard Deviations */

    cpgsvp (lm+ft+2*tt, lm+2*ft+2*tt, bm, bm+fh);
    cpgswin(0.0, numchan, 0.0, T);
    {
      int mincol, maxcol, numcol;
      float l[2] = {0.0, 1.0};
      float r[2] = {1.0, 0.0};
      float g[2] = {1.0, 0.0};
      float b[2] = {1.0, 0.0};
      float fg = 0.0, bg = 0.0, tr[6], *levels;
      float x1 = 0.0, x2 = numchan, y1 = 0.0, y2 = T;

      cpgqcol(&mincol, &maxcol);
      mincol += 2;
      cpgscir(mincol, maxcol);
      numcol = maxcol - mincol + 1;
      levels = gen_fvect(numcol);
      cpgctab(l, r, g, b, numcol, 1.0, 0.5);
      autocal2d(datapow[0], numint, numchan, &fg, &bg, numcol,
		levels, &x1, &x2, &y1, &y2, tr);
      cpgimag(datapow[0], numchan, numint, 0+1, 
	      numchan, 0+1, numint, bg, fg, tr);
      free(levels);
    }
    cpgswin(0.0, numchan, 0.0, T);
    cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
    cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");
    cpgswin(lof, hif, 0.0, numint);
    cpgbox ("CST", 0.0, 0, "CST", 0.0, 0);

    /*  Max Power versus Max Power */

    find_min_max_arr(numint, pow_int, &min, &max);
    cpgsvp (lm+ft, lm+ft+tt, bm, bm+fh);
    cpgswin(1.1*max, 0.9*min, 0.0, T);
    cpgbox ("BCNST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Power");
    cpgline(numint, pow_int, times);
    cpgscr(2, 1.0, 0.0, 0.0); cpgsci(2); /* Red */
    yarr[0] = 0.0; yarr[1] = T;
    xarr[0] = xarr[1] = pow_int_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = pow_int_med + sigma * pow_int_std;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = pow_int_med - sigma * pow_int_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgswin(1.1*max, 0.9*min, 0.0, numint);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number");

    /*  Max Power versus Channel */

    find_min_max_arr(numchan, pow_chan, &min, &max);
    cpgsvp (lm, lm+ft, bm+fh, bm+fh+th);
    cpgswin(0.0, numchan, 0.9*min, 1.1*max);
    cpgbox ("BST", 0.0, 0, "BCNST", 0.0, 0);
    cpgmtxt("L", 2.1, 0.5, 0.5, "Power");
    cpgline(numchan, chans, pow_chan);
    cpgscr(2, 1.0, 0.0, 0.0); cpgsci(2); /* Red */
    xarr[0] = 0.0; xarr[1] = numchan;
    yarr[0] = yarr[1] = pow_chan_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = pow_chan_med + sigma * pow_chan_std;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = pow_chan_med - sigma * pow_chan_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgswin(lof, hif, 0.9*min, 1.1*max);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Add the Data Info area */

    {
      char out[200], out2[100];

      cpgsvp (lm+ft+2*tt, 1.0-rm, bm+fh+th, 1.0-tm);
      cpgswin(-0.1, 1.00, -0.1, 1.1);
      cpgsch(1.0);
      sprintf(out, "%-s", idata->name);
      cpgmtxt("T", 1.0, 0.5, 0.5, out);
      cpgsch(0.7);
      sprintf(out, "Object:  %-s", idata->object);
      cpgtext(0.0, 1.0, out);
      sprintf(out, "Telescope:  %-s", idata->telescope);
      cpgtext(0.0, 0.9, out);
      sprintf(out, "Instrument:  %-s", idata->instrument);
      cpgtext(0.0, 0.8, out);
      ra_dec_to_string(out2, idata->ra_h, idata->ra_m, idata->ra_s);
      sprintf(out, "RA:");
      cpgtext(0.0, 0.7, out);
      sprintf(out, "=  %-s", out2);
      cpgtext(0.15, 0.7, out);
      ra_dec_to_string(out2, idata->dec_d, idata->dec_m, idata->dec_s);
      sprintf(out, "DEC:");
      cpgtext(0.0, 0.6, out);
      sprintf(out, "=  %-s", out2);
      cpgtext(0.15, 0.6, out);
      sprintf(out, "Epoch\\dtopo\\u = %-.11f", 
	      idata->mjd_i + idata->mjd_f);
      cpgtext(0.0, 0.5, out);
      sprintf(out, "T\\dsample\\u");
      cpgtext(0.0, 0.4, out);
      sprintf(out, "=  %f", idata->dt);
      cpgtext(0.45, 0.4, out);
    }
    cpgclos();
  }

  /* Free our arrays */

  free(avg_chan);
  free(avg_int);
  free(std_chan);
  free(std_int);
  free(pow_chan);
  free(pow_int);
}
