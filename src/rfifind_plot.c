#include "presto.h"
#include "rfifind.h"
#include "plot2d.h"
#include "mask.h"

int compare_floats(const void *a, const void *b);
/* qsort comparison function for floats */

static double single_power_pdf(double power, int numpowers)
{
  double tmp;
  
  tmp = (power > 100.0) ? 0.0 : exp(-power);
  return tmp * numpowers * pow(1.0-tmp, numpowers-1);
}

static void plot_stats(float *median, float *average, float *stdev);

static void calc_avgmedstd(float *arr, int numarr, float fraction, 
			   int step, float *avg, float *med, float *std)
/* Calculates the median and middle-'fraction' std deviation  */
/* and average of the array 'arr'.  Values are returned in    */
/* 'avg', 'med' and 'std'.  The array is not modified.        */
{
  int ii, jj, len, start;
  float *tmparr;
  double davg, dstd;

  len = (int)(numarr * fraction + 0.5);
  if (len > numarr || len < 0){
    printf("fraction (%g) out-of-bounds in calc_avgmedstd()\n", 
	   fraction);
    exit(1);
  }
  start = (numarr - len) / 2;
  tmparr = gen_fvect(numarr);
  for (ii=0, jj=0; ii<numarr; ii++, jj+=step)
    tmparr[ii] = arr[jj];
  qsort(tmparr, numarr, sizeof(float), compare_floats);
  avg_var(tmparr+start, len, &davg, &dstd);
  *avg = (float) davg;
  *med = tmparr[numarr/2];
  *std = sqrt(dstd);
  free(tmparr);
}

void rfifind_plot(int numchan, int numint, int ptsperint, 
		  float timesigma, float freqsigma, 
		  float **dataavg, float **datastd, float **datapow,
		  int *userchan, int numuserchan, 
		  int *userints, int numuserints, 
		  infodata *idata, unsigned char **bytemask, 
		  mask *oldmask, mask *newmask, int xwin)
/* Make the beautiful multi-page rfifind plots */
{
  int ii, jj, ct, loops=1;
  float *freqs, *chans, *times, *ints;
  float *avg_chan_avg, *std_chan_avg, *pow_chan_avg;
  float *avg_chan_med, *std_chan_med, *pow_chan_med;
  float *avg_chan_std, *std_chan_std, *pow_chan_std;
  float *avg_int_avg,  *std_int_avg,  *pow_int_avg;
  float *avg_int_med,  *std_int_med,  *pow_int_med;
  float *avg_int_std,  *std_int_std,  *pow_int_std;
  float dataavg_avg, datastd_avg, datapow_avg;
  float dataavg_med, datastd_med, datapow_med;
  float dataavg_std, datastd_std, datapow_std;
  float avg_reject, std_reject, pow_reject;
  double inttim, T, lof, hif;
  
  inttim = ptsperint * idata->dt;
  T = inttim * numint;
  lof = idata->freq - 0.5 * idata->chan_wid;
  hif = lof + idata->freqband;
  avg_chan_avg = gen_fvect(numchan);
  std_chan_avg = gen_fvect(numchan);
  pow_chan_avg = gen_fvect(numchan);
  avg_int_avg = gen_fvect(numint);
  std_int_avg = gen_fvect(numint);
  pow_int_avg = gen_fvect(numint);
  avg_chan_med = gen_fvect(numchan);
  std_chan_med = gen_fvect(numchan);
  pow_chan_med = gen_fvect(numchan);
  avg_int_med = gen_fvect(numint);
  std_int_med = gen_fvect(numint);
  pow_int_med = gen_fvect(numint);
  avg_chan_std = gen_fvect(numchan);
  std_chan_std = gen_fvect(numchan);
  pow_chan_std = gen_fvect(numchan);
  avg_int_std = gen_fvect(numint);
  std_int_std = gen_fvect(numint);
  pow_int_std = gen_fvect(numint);
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
    times[ii] = 0.0 + ii * inttim;
  }

  /* Calculate the statistics of the full set */

  ct = numchan * numint;
  calc_avgmedstd(dataavg[0], ct, 0.8, 1, &dataavg_avg, 
		 &dataavg_med, &dataavg_std);
  calc_avgmedstd(datastd[0], ct, 0.8, 1, &datastd_avg, 
		 &datastd_med, &datastd_std);
  calc_avgmedstd(datapow[0], ct, 0.5, 1, &datapow_avg, 
		 &datapow_med, &datapow_std);
  avg_reject = timesigma * dataavg_std;
  std_reject = timesigma * datastd_std;
  pow_reject = power_for_sigma(freqsigma, 1, ptsperint/2);

  /* Calculate the channel/integration statistics vectors */

  for (ii=0; ii<numint; ii++){
    calc_avgmedstd(dataavg[0]+ii*numchan, numchan, 0.8, 1, 
		   avg_int_avg+ii, avg_int_med+ii, avg_int_std+ii);
    calc_avgmedstd(datastd[0]+ii*numchan, numchan, 0.8, 1, 
		   std_int_avg+ii, std_int_med+ii, std_int_std+ii);
    calc_avgmedstd(datapow[0]+ii*numchan, numchan, 0.5, 1, 
		   pow_int_avg+ii, pow_int_med+ii, pow_int_std+ii);
  }
  for (ii=0; ii<numchan; ii++){
    calc_avgmedstd(dataavg[0]+ii, numint, 0.8, numchan, 
		   avg_chan_avg+ii, avg_chan_med+ii, avg_chan_std+ii);
    calc_avgmedstd(datastd[0]+ii, numint, 0.8, numchan, 
		   std_chan_avg+ii, std_chan_med+ii, std_chan_std+ii);
    calc_avgmedstd(datapow[0]+ii, numint, 0.5, numchan, 
		   pow_chan_avg+ii, pow_chan_med+ii, pow_chan_std+ii);
    /*
    fprintf(stderr, "%12.7g  %12.7g  %12.7g    %12.7g  %12.7g  %12.7g    %12.7g  %12.7g  %12.7g    \n", 
	    avg_chan_avg[ii], avg_chan_med[ii], avg_chan_std[ii],
	    std_chan_avg[ii], std_chan_med[ii], std_chan_std[ii],
	    pow_chan_avg[ii], pow_chan_med[ii], pow_chan_std[ii]);
    */
  }

  /* Generate the byte mask */

  /* Set the channels/intervals picked by the user */
  if (numuserints)
    for (ii=0; ii<numuserints; ii++)
      if (userints[ii] >= 0 && userints[ii] < numint)
	for (jj=0; jj<numchan; jj++)
	  bytemask[userints[ii]][jj] |= USERINTS;
  if (numuserchan)
    for (ii=0; ii<numuserchan; ii++)
      if (userchan[ii] >= 0 && userchan[ii] < numchan)
	for (jj=0; jj<numint; jj++)
	  bytemask[jj][userchan[ii]] |= USERCHAN;

  /* Compare each point in an interval (or channel) with   */
  /* the interval's (or channel's) median and the overall  */
  /* standard deviation.  If the channel/integration       */
  /* medians are more than sigma different than the global */
  /* value, set them to the global.                        */
  {
    float int_med, chan_med;

    for (ii=0; ii<numint; ii++){
      for (jj=0; jj<numchan; jj++){
	{ /* Powers */
	  if (fabs(pow_int_med[ii] - datapow_med) > 
	      freqsigma * datapow_std)
	    int_med = datapow_med;
	  else
	    int_med = pow_int_med[ii];
	  if (fabs(pow_chan_med[ii] - datapow_med) > 
	      freqsigma * datapow_std)
	    chan_med = datapow_med;
	  else
	    chan_med = pow_chan_med[ii];
	  if (fabs(datapow[ii][jj] - int_med) > pow_reject ||
	      fabs(datapow[ii][jj] - chan_med) > pow_reject)
	    if (!(bytemask[ii][jj] & PADDING)) 
	      bytemask[ii][jj] |= BAD_POW;
	}
	{ /* Averages */
	  if (fabs(avg_int_med[ii] - dataavg_med) > 
	      timesigma * dataavg_std)
	    int_med = dataavg_med;
	  else
	    int_med = avg_int_med[ii];
	  if (fabs(avg_chan_med[ii] - dataavg_med) > 
	      timesigma * dataavg_std)
	    chan_med = dataavg_med;
	  else
	    chan_med = avg_chan_med[ii];
	  if (fabs(dataavg[ii][jj] - int_med) > avg_reject ||
	      fabs(dataavg[ii][jj] - chan_med) > avg_reject)
	    if (!(bytemask[ii][jj] & PADDING)) 
	      bytemask[ii][jj] |= BAD_AVG;
	}
	{ /* Standard Deviations */
	  if (fabs(std_int_med[ii] - datastd_med) > 
	      timesigma * datastd_std)
	    int_med = datastd_med;
	  else
	    int_med = std_int_med[ii];
	  if (fabs(std_chan_med[ii] - datastd_med) > 
	      timesigma * datastd_std)
	    chan_med = datastd_med;
	  else
	    chan_med = std_chan_med[ii];
	  if (fabs(datastd[ii][jj] - int_med) > std_reject ||
	      fabs(datastd[ii][jj] - chan_med) > std_reject)
	    if (!(bytemask[ii][jj] & PADDING)) 
	      bytemask[ii][jj] |= BAD_STD;
	}
      }
    }
  }

  /* Generate the New Mask */

  fill_mask(timesigma, freqsigma, idata->mjd_i+idata->mjd_f, 
	    ptsperint * idata->dt, idata->freq, idata->chan_wid, 
	    numchan, numint, ptsperint, numuserchan, userchan, 
	    numuserints, userints, bytemask, newmask);

  /* Place the oldmask over the newmask for plotting purposes */

  if (oldmask->numchan)
    set_oldmask_bits(oldmask, bytemask);

  /*
   *  Now plot the results
   */

  if (xwin) loops = 2;
  for (ct=0; ct<loops; ct++){
    float min, max, tr[6], locut, hicut;
    float left, right, top, bottom;
    float xl, xh, yl, yh;
    float tt, ft, th, fh;  /* thin and fat thicknesses and heights */
    float lm, rm, tm, bm;  /* LRTB margins */
    float xarr[2], yarr[2];
    char outdev[100];
    int ii, mincol, maxcol, numcol;

    /*Set the PGPLOT device to an X-Window */
    
    if (ct==1)
      strcpy(outdev, "/XWIN");
    else
      sprintf(outdev, "%s.ps/CPS", idata->name);

    /* Open and prep our device */

    cpgopen(outdev);
    cpgpap(10.25, 8.5/11.0);
    cpgpage();
    cpgiden();
    cpgsch(0.7);
    cpgqcir(&mincol, &maxcol);
    numcol = maxcol - mincol + 1;
    for (ii=mincol; ii<=maxcol; ii++){
      float color;
      color = (float) (maxcol - ii) / (float) numcol;
      cpgscr(ii, color, color, color);
    }
    
    /* Set thicknesses and margins */

    lm = 0.04;
    rm = 0.04;
    bm = 0.08;
    tm = 0.05;
    ft = 3.0;  /* This sets fat thickness = 3 x thin thickness */
    tt = 0.92 / (6.0 + 4.0 * ft);
    ft *= tt;
    fh = 0.55;
    th = tt * 11.0/8.5;

    /* Powers Histogram */
    {
      float *theo, *hist, *hpows, *tpows, maxhist=0.0, maxtheo=0.0;
      int numhist=40, numtheo=200, bin, numpows;
      double dtheo, dhist, spacing;

      /* Calculate the predicted distribution of max powers */


      numpows = numint * numchan;
      find_min_max_arr(numpows, datapow[0], &min, &max);
      min = (min==0.0) ? log10(5.0 * 0.95) : log10(min * 0.95);
      max = log10(max * 1.05);
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
	bin = (*(datapow[0]+ii)==0.0) ? 0 : \
	  (log10(*(datapow[0]+ii))-min)/dhist;
	if (bin < 0)
	  bin = 0;
	if (bin >= numhist)
	  bin = numhist;
	hist[bin] += 1.0;
      }
      for (ii=0; ii<numhist; ii++)
	if (hist[ii] > maxhist) maxhist = hist[ii];
      maxhist *= 1.1;
      dtheo = (max-min)/(double)(numtheo - 1);
      for (ii=0; ii<numtheo; ii++){
	tpows[ii] = min + ii * dtheo;
	theo[ii] = single_power_pdf(pow(10.0,tpows[ii]), 
				    ptsperint/2)*numpows;
	spacing = (pow(10.0, tpows[ii]+dhist) -
		   pow(10.0, tpows[ii]));
	theo[ii] *= spacing;
	if (theo[ii] > maxtheo) maxtheo = theo[ii];
      }
      maxtheo *= 1.1;
      if (maxtheo > maxhist) maxhist = maxtheo;
      left   = lm;
      right  = lm + ft + tt;
      bottom = 0.80;
      top    = 0.96;
      cpgsvp (left, right, bottom, top);
      xl = min;
      xh = max;
      yl = 0.0;
      yh = maxhist;
      cpgswin(xl, xh, yl, yh);
      cpgmtxt("L", 1.1, 0.5, 0.5, "Number");
      cpgmtxt("B", 2.1, 0.5, 0.5, "Max Power");
      cpgbin(numhist, hpows, hist, 0);
      cpgscr(maxcol, 0.5, 0.5, 0.5); cpgsci(maxcol); /* Grey */
      cpgline(numtheo, tpows, theo);
      xarr[0] = log10(power_for_sigma(freqsigma, 1, ptsperint/2));
      xarr[1] = xarr[0];
      yarr[0] = yl; yarr[1] = yh;
      cpgsls(4); /* Dotted line */
      cpgscr(maxcol, 1.0, 0.0, 0.0); cpgsci(maxcol); /* Red */
      cpgline(2, xarr, yarr);
      cpgsls(1); /* Solid line */
      cpgsci(1); /* Default color */
      cpgbox ("BCLNST", 0.0, 0, "BC", 0.0, 0);
      free(hist);
      free(theo);
      free(tpows);
      free(hpows);
    }

    /* Maximum Powers */

    left   = lm;
    right  = lm + ft;
    bottom = bm;
    top    = bm + fh;
    xl = 0.0;
    xh = numchan;
    yl = 0.0;
    yh = T;
    cpgsvp (left, right, bottom, top);
    cpgswin(xl, xh, yl, yh);
    cpgscr(maxcol, 1.0, 0.0, 0.0); /* Red */
    locut = 0.0;
    hicut = pow_reject;
    tr[2] = tr[4] = 0.0;    
    tr[1] = (xh - xl) / numchan;
    tr[0] = xl - (tr[1] / 2);
    tr[5] = (yh - yl) / numint;
    tr[3] = yl - (tr[5] / 2);
    cpgimag(datapow[0], numchan, numint, 1, 
	    numchan, 1, numint, locut, hicut, tr);
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
    cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");
    xl = lof;
    xh = hif;
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CST", 0.0, 0, "CST", 0.0, 0);

    /*  Max Power versus Time */

    left   = lm + ft;
    right  = lm + ft + tt;
    bottom = bm;
    top    = bm + fh;
    cpgsvp (left, right, bottom, top);
    find_min_max_arr(numint, pow_int_med, &min, &max);
    xl = 0.0;
    xh = 1.5 * pow_reject;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgscr(maxcol, 1.0, 0.0, 0.0); cpgsci(maxcol); /* Red */
    cpgmtxt("B", 1.2, 0.5, 0.5, "Max Power");
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = datapow_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = pow_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, pow_int_med, times);
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    /* cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number"); */

    /*  Max Power versus Channel */

    left   = lm;
    right  = lm + ft;
    bottom = bm + fh;
    top    = bm + fh + th;
    cpgsvp (left, right, bottom, top);
    find_min_max_arr(numint, pow_chan_med, &min, &max);
    xl = 0.0;
    xh = numchan;
    yl = 0.0;
    yh = 1.5 * pow_reject;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgscr(maxcol, 1.0, 0.0, 0.0); cpgsci(maxcol); /* Red */
    cpgmtxt("L", 0.8, 0.5, 0.5, "Max Power");
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = datapow_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = pow_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, pow_chan_med);
    xl = lof;
    xh = hif;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Standard Deviations */

    left   = lm + ft + 2.0 * tt;
    right  = lm + 2.0 * ft + 2.0 * tt;
    bottom = bm;
    top    = bm + fh;
    xl = 0.0;
    xh = numchan;
    yl = 0.0;
    yh = T;
    cpgsvp (left, right, bottom, top);
    cpgswin(xl, xh, yl, yh);
    cpgscr(mincol, 0.7, 1.0, 0.7); /* Light Green */
    cpgscr(maxcol, 0.3, 1.0, 0.3); /* Dark Green */
    locut = datastd_med - timesigma * datastd_std;
    hicut = datastd_med + timesigma * datastd_std;
    tr[2] = tr[4] = 0.0;    
    tr[1] = (xh - xl) / numchan;
    tr[0] = xl - (tr[1] / 2);
    tr[5] = (yh - yl) / numint;
    tr[3] = yl - (tr[5] / 2);
    cpgimag(datastd[0], numchan, numint, 1, 
	    numchan, 1, numint, locut, hicut, tr);
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
    xl = lof;
    xh = hif;
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CST", 0.0, 0, "CST", 0.0, 0);

    /*  Max sigma versus Time */

    left   = lm + 2.0 * ft + 2.0 * tt;
    right  = lm + 2.0 * ft + 3.0 * tt;
    bottom = bm;
    top    = bm + fh;
    cpgsvp (left, right, bottom, top);
    xl = datastd_med - 2.0 * std_reject;
    xh = datastd_med + 2.0 * std_reject;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgscr(maxcol, 0.0, 1.0, 0.0); cpgsci(maxcol); /* Green */
    cpgmtxt("B", 1.2, 0.5, 0.5, "Data \\gs");
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = datastd_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = datastd_med + std_reject;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = datastd_med - std_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, std_int_med, times);
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    /* cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number"); */

    /*  Max Sigma versus Channel */

    left   = lm + ft + 2.0 * tt;
    right  = lm + 2.0 * ft + 2.0 * tt;
    bottom = bm + fh;
    top    = bm + fh + th;
    cpgsvp (left, right, bottom, top);
    xl = 0.0;
    xh = numchan;
    yl = datastd_med - 2.0 * std_reject;
    yh = datastd_med + 2.0 * std_reject;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgscr(maxcol, 0.0, 1.0, 0.0); cpgsci(maxcol); /* Green */
    cpgmtxt("L", 0.8, 0.5, 0.5, "Data \\gs");
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = datastd_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = datastd_med + std_reject;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = datastd_med - std_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, std_chan_med);
    xl = lof;
    xh = hif;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Averages */

    left   = lm + 2.0 * ft + 4.0 * tt;
    right  = lm + 3.0 * ft + 4.0 * tt;
    bottom = bm;
    top    = bm + fh;
    xl = 0.0;
    xh = numchan;
    yl = 0.0;
    yh = T;
    cpgsvp (left, right, bottom, top);
    cpgswin(xl, xh, yl, yh);
    cpgscr(mincol, 0.7, 0.7, 1.0); /* Light Blue */
    cpgscr(maxcol, 0.3, 0.3, 1.0); /* Dark Blue */
    locut = dataavg_med - timesigma * dataavg_std;
    hicut = dataavg_med + timesigma * dataavg_std;
    tr[2] = tr[4] = 0.0;    
    tr[1] = (xh - xl) / numchan;
    tr[0] = xl - (tr[1] / 2);
    tr[5] = (yh - yl) / numint;
    tr[3] = yl - (tr[5] / 2);
    cpgimag(dataavg[0], numchan, numint, 1, 
	    numchan, 1, numint, locut, hicut, tr);
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
    cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
    xl = lof;
    xh = hif;
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CST", 0.0, 0, "CST", 0.0, 0);

    /*  Max Average versus Time */

    left   = lm + 3.0 * ft + 4.0 * tt;
    right  = lm + 3.0 * ft + 5.0 * tt;
    bottom = bm;
    top    = bm + fh;
    cpgsvp (left, right, bottom, top);
    xl = dataavg_med - 2.0 * avg_reject;
    xh = dataavg_med + 2.0 * avg_reject;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgscr(maxcol, 0.0, 0.0, 1.0); cpgsci(maxcol); /* Blue */
    cpgmtxt("B", 1.2, 0.5, 0.5, "Data Mean");
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = dataavg_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = dataavg_med + avg_reject;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = dataavg_med - avg_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, avg_int_med, times);
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);

    /*  Max Average versus Channel */

    left   = lm + 2.0 * ft + 4.0 * tt;
    right  = lm + 3.0 * ft + 4.0 * tt;
    bottom = bm + fh;
    top    = bm + fh + th;
    cpgsvp (left, right, bottom, top);
    xl = 0.0;
    xh = numchan;
    yl = dataavg_med - 2.0 * avg_reject;
    yh = dataavg_med + 2.0 * avg_reject;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgscr(maxcol, 0.0, 0.0, 1.0); cpgsci(maxcol); /* Blue */
    cpgmtxt("L", 0.8, 0.5, 0.5, "Data Mean");
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = dataavg_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = dataavg_med + avg_reject;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = dataavg_med - avg_reject;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, avg_chan_med);
    xl = lof;
    xh = hif;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Plot the Mask */

    {
      unsigned char byte;
      float **plotmask, rr, gg, bb;

      plotmask = gen_fmatrix(numint, numchan);
      for (ii=0; ii<numint; ii++){
	for (jj=0; jj<numchan; jj++){
	  byte = bytemask[ii][jj];
	  plotmask[ii][jj] = 0.0;
	  if (byte & PADDING) plotmask[ii][jj] = 1.0;
	  if (byte & OLDMASK) plotmask[ii][jj] = 2.0;
	  if (byte & USERZAP) plotmask[ii][jj] = 3.0;
	  if (byte & BAD_POW) plotmask[ii][jj] = 4.0;
	  else if (byte & BAD_AVG) plotmask[ii][jj] = 5.0;
	  else if (byte & BAD_STD) plotmask[ii][jj] = 6.0;
	}
      }
      /* Set the colors */
      numcol = 7;
      maxcol = mincol + numcol - 1;
      cpgscir(mincol, maxcol);
      cpgqcr(0, &rr, &gg, &bb);
      cpgscr(mincol+0, rr,  gg,  bb ); /* GOODDATA = background */
      cpgscr(mincol+1, 0.7, 0.7, 0.7); /* PADDING  = light grey */
      cpgscr(mincol+2, 0.3, 0.3, 0.3); /* OLDMASK  = dark grey */
      cpgqcr(1, &rr, &gg, &bb);
      cpgscr(mincol+3, rr,  gg,  bb ); /* USERZAP  = foreground */
      cpgscr(mincol+4, 1.0, 0.0, 0.0); /* BAD+POW  = red */
      cpgscr(mincol+5, 0.0, 0.0, 1.0); /* BAD+AVG  = blue */
      cpgscr(mincol+6, 0.0, 1.0, 0.0); /* BAD+STD  = green */
      /* Prep the image */
      left   = lm + 3.0 * ft + 6.0 * tt;
      right  = lm + 4.0 * ft + 6.0 * tt;
      bottom = bm;
      top    = bm + fh;
      xl = 0.0;
      xh = numchan;
      yl = 0.0;
      yh = T;
      cpgsvp (left, right, bottom, top);
      cpgswin(xl, xh, yl, yh);
      locut = 0.0;
      hicut = 6.0;
      tr[2] = tr[4] = 0.0;    
      tr[1] = (xh - xl) / numchan;
      tr[0] = xl - (tr[1] / 2);
      tr[5] = (yh - yl) / numint;
      tr[3] = yl - (tr[5] / 2);
      cpgimag(plotmask[0], numchan, numint, 1, 
	      numchan, 1, numint, locut, hicut, tr);
      cpgswin(xl, xh, yl, yh);
      cpgbox ("BNST", 0.0, 0, "BNST", 0.0, 0);
      cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
      xl = lof;
      xh = hif;
      yl = 0.0;
      yh = numint;
      cpgswin(xl, xh, yl, yh);
      cpgbox ("CMST", 0.0, 0, "CMST", 0.0, 0);
      cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");
      cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number");
      /* Add the Labels */
      cpgsvp (0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      cpgsch(0.7);
      cpgsci(mincol+1);
      cpgptxt(left, top+0.1, 0.0, 0.0, "Padding");
      cpgsci(mincol+2);
      cpgptxt(left, top+0.08, 0.0, 0.0, "Old Mask");
      cpgsci(mincol+3);
      cpgptxt(left, top+0.06, 0.0, 0.0, "User Zap");
      cpgsci(mincol+4);
      cpgptxt(right, top+0.1, 0.0, 1.0, "Max Power");
      cpgsci(mincol+6);
      cpgptxt(right, top+0.08, 0.0, 1.0, "Data \\gs");
      cpgsci(mincol+5);
      cpgptxt(right, top+0.06, 0.0, 1.0, "Data Mean");
      cpgsci(1);
      free(plotmask[0]); free(plotmask);
    }

    /* Add the Data Info area */

    {
      char out[200], out2[100];
      float dy = 0.025;
      
      cpgsvp (0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      left = lm + ft + 1.5 * tt;
      top = 1.0 - tm;
      cpgsch(1.0);
      sprintf(out, "%-s", idata->name);
      cpgtext(0.5, 1.0 - 0.5*tm, out);
      cpgsch(0.8);

      sprintf(out, "Object:");
      cpgtext(left+0.0,  top-0*dy, out);
      sprintf(out, "%-s", idata->object);
      cpgtext(left+0.1, top-0*dy, out);
      sprintf(out, "Telescope:");
      cpgtext(left+0.0,  top-1*dy, out);
      sprintf(out, "%-s", idata->telescope);
      cpgtext(left+0.1, top-1*dy, out);
      sprintf(out, "Instrument:");
      cpgtext(left+0.0,  top-2*dy, out);
      sprintf(out, "%-s", idata->instrument);
      cpgtext(left+0.1, top-2*dy, out);
      ra_dec_to_string(out2, idata->ra_h, idata->ra_m, idata->ra_s);
      sprintf(out, "RA\\dJ2000\\u");
      cpgtext(left+0.0,   top-3*dy, out);
      sprintf(out, "= %-s", out2);
      cpgtext(left+0.08, top-3*dy, out);
      ra_dec_to_string(out2, idata->dec_d, idata->dec_m, idata->dec_s);
      sprintf(out, "DEC\\dJ2000\\u");
      cpgtext(left+0.0,   top-4*dy, out);
      sprintf(out, "= %-s", out2);
      cpgtext(left+0.08, top-4*dy, out);
      sprintf(out, "Epoch\\dtopo\\u");
      cpgtext(left+0.0,   top-5*dy, out);
      sprintf(out, "= %-.11f", idata->mjd_i + idata->mjd_f);
      cpgtext(left+0.08, top-5*dy, out);
      sprintf(out, "T\\dsample\\u (s)");
      cpgtext(left+0.0,   top-6*dy, out);
      sprintf(out, "= %g", idata->dt);
      cpgtext(left+0.08, top-6*dy, out);
      sprintf(out, "T\\dtotal\\u (s)");
      cpgtext(left+0.0,   top-7*dy, out);
      sprintf(out, "= %g", T);
      cpgtext(left+0.08, top-7*dy, out);

      left = lm + ft + 7.8 * tt;
      sprintf(out, "Num channels");
      cpgtext(left+0.0,  top-0*dy, out);
      sprintf(out, "= %-d", numchan);
      cpgtext(left+0.12, top-0*dy, out);
      sprintf(out, "Pts per int");
      cpgtext(left+0.19, top-0*dy, out);
      sprintf(out, "= %-d", ptsperint);
      cpgtext(left+0.29, top-0*dy, out);
      sprintf(out, "Num intervals");
      cpgtext(left+0.0,  top-1*dy, out);
      sprintf(out, "= %-d", numint);
      cpgtext(left+0.12, top-1*dy, out);
      sprintf(out, "Time per int");
      cpgtext(left+0.19,  top-1*dy, out);
      sprintf(out, "= %-g", inttim);
      cpgtext(left+0.29, top-1*dy, out);
      sprintf(out, "Power:");
      cpgtext(left+0.0,  top-2*dy, out);
      sprintf(out, "median");
      cpgtext(left+0.06,  top-2*dy, out);
      sprintf(out, "= %-.3f", datapow_med);
      cpgtext(left+0.12,  top-2*dy, out);
      sprintf(out, "\\gs");
      cpgtext(left+0.21,   top-2*dy, out);
      sprintf(out, "= %-.3g", datapow_std);
      cpgtext(left+0.245,   top-2*dy, out);
      find_min_max_arr(numint*numchan, datapow[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-3*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-3*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.21,   top-3*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.245,   top-3*dy, out);
      sprintf(out, "Data \\gs:");
      cpgtext(left+0.0,  top-4*dy, out);
      sprintf(out, "median");
      cpgtext(left+0.06,  top-4*dy, out);
      sprintf(out, "= %-.3f", datastd_med);
      cpgtext(left+0.12,  top-4*dy, out);
      sprintf(out, "\\gs");
      cpgtext(left+0.21,   top-4*dy, out);
      sprintf(out, "= %-.3g", datastd_std);
      cpgtext(left+0.245,   top-4*dy, out);
      find_min_max_arr(numint*numchan, datastd[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-5*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-5*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.21,   top-5*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.245,   top-5*dy, out);
      sprintf(out, "Mean:");
      cpgtext(left+0.0,  top-6*dy, out);
      sprintf(out, "median");
      cpgtext(left+0.06,  top-6*dy, out);
      sprintf(out, "= %-.3f", dataavg_med);
      cpgtext(left+0.12,  top-6*dy, out);
      sprintf(out, "\\gs");
      cpgtext(left+0.21,   top-6*dy, out);
      sprintf(out, "= %-.3g", dataavg_std);
      cpgtext(left+0.245,   top-6*dy, out);
      find_min_max_arr(numint*numchan, dataavg[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-7*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-7*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.21,   top-7*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.245,   top-7*dy, out);
    }
    cpgclos();
  }

  /* Free our arrays */
    
  free(freqs); free(chans); free(times); free(ints);
  free(avg_chan_avg); free(std_chan_avg); free(pow_chan_avg);
  free(avg_int_avg); free(std_int_avg); free(pow_int_avg);
  free(avg_chan_med); free(std_chan_med); free(pow_chan_med);
  free(avg_int_med); free(std_int_med); free(pow_int_med);
  free(avg_chan_std); free(std_chan_std); free(pow_chan_std);
  free(avg_int_std); free(std_int_std); free(pow_int_std);
}
