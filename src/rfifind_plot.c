#include "presto.h"
#include "rfifind.h"
#include "plot2d.h"

double single_power_pdf(double power, int numpowers)
{
  double tmp;
  
  tmp = exp(-power);
  return tmp * numpowers * pow(1.0-tmp, numpowers-1);
}

void rfifind_plot(int numchan, int numint, int ptsperint, float sigma, 
		  float **dataavg, float **datastd, float **datapow,
		  int *userchan, int numuserchan, unsigned char **mask, 
		  infodata *idata, int xwin)
/* Make the beautiful 1 page rfifind output */
{
  int ii, jj, ct, loops=1;
  float *freqs, *chans, *times, *ints;
  float *avg_chan, *std_chan, *pow_chan;
  float *avg_int, *std_int, *pow_int;
  double inttim, T, lof, hif;
  double avgavg, stdavg, powavg, maxpow;
  double dataavg_med, datastd_med, datapow_med;
  double dataavg_std, datastd_std, datapow_std;
  double avg_chan_std, std_chan_std, pow_chan_std;
  double avg_chan_med, std_chan_med, pow_chan_med;
  double avg_int_std, std_int_std, pow_int_std;
  double avg_int_med, std_int_med, pow_int_med;
  
  inttim = ptsperint * idata->dt;
  T = inttim * numint;
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
    times[ii] = 0.0 + ii * inttim;
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
  avg_var(dataavg[0], numchan*numint, &dataavg_med, &dataavg_std);
  avg_var(datastd[0], numchan*numint, &datastd_med, &datastd_std);
  avg_var(datapow[0], numchan*numint, &datapow_med, &datapow_std);
  /* Turn variances into standard deviations */
  avg_int_std  = sqrt(avg_int_std);
  std_int_std  = sqrt(std_int_std);
  pow_int_std  = sqrt(pow_int_std);
  avg_chan_std = sqrt(avg_chan_std);
  std_chan_std = sqrt(std_chan_std);
  pow_chan_std = sqrt(pow_chan_std);
  dataavg_std = sqrt(dataavg_std);
  datastd_std = sqrt(datastd_std);
  datapow_std = sqrt(datapow_std);
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
  maxpow = power_for_sigma(sigma, 1, ptsperint/2);

  /* Generate the Mask */

  /*
  {
    int numbytes, ii;

    numbytes = (numchan % 8)  ? numchan  / 8 + 1 : numchan  / 8;
    *mask = gen_bvect(numbytes);
    for (ii=0; ii<numbytes; ii++)
      *mask[ii]
  }
  */

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
      color = (float) (ii - mincol) / (float) numcol;
      cpgscr(ii, color, color, color);
    }
    
    /* Set thicknesses and margins */

    lm = rm = 0.04;
    bm = 0.08;
    tm = 0.05;
    ft = 3.0;  /* This sets fat thickness = 3 x thin thickness */
    tt = 0.95 / (6.0 + 4.0 * ft);
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
      min = log10(min * 0.95);
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
	bin = (log10(*(datapow[0]+ii))-min)/dhist;
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
      cpgmtxt("B", 2.1, 0.5, 0.5, "Power");
      cpgbin(numhist, hpows, hist, 0);
      cpgscr(maxcol, 0.5, 0.5, 0.5); cpgsci(maxcol); /* Grey */
      cpgline(numtheo, tpows, theo);
      xarr[0] = log10(power_for_sigma(sigma, 1, ptsperint/2));
      xarr[1] = xarr[0];
      yarr[0] = yl; yarr[1] = yh;
      cpgsls(4); /* Dotted line */
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
    hicut = maxpow;
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
    find_min_max_arr(numint, pow_int, &min, &max);
    xl = 0.9 * min;
    xh = pow_int_med + 2.0 * sigma * pow_int_std;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 1.2, 0.5, 0.5, "Power");
    cpgscr(maxcol, 1.0, 0.0, 0.0); cpgsci(maxcol); /* Red */
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = pow_int_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = pow_int_med + sigma * pow_int_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, pow_int, times);
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
    find_min_max_arr(numint, pow_chan, &min, &max);
    xl = 0.0;
    xh = numchan;
    yl = 0.9 * min;
    yh = pow_chan_med + 2.0 * sigma * pow_chan_std;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgmtxt("L", 0.8, 0.5, 0.5, "Power");
    cpgscr(maxcol, 1.0, 0.0, 0.0); cpgsci(maxcol); /* Red */
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = pow_chan_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = pow_chan_med + sigma * pow_chan_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, pow_chan);
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
    locut = datastd_med - sigma * datastd_std;
    hicut = datastd_med + sigma * datastd_std;
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
    xl = std_int_med - 2.0 * sigma * std_int_std;
    xh = std_int_med + 2.0 * sigma * std_int_std;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 1.2, 0.5, 0.5, "Data \\gs");
    cpgscr(maxcol, 0.0, 1.0, 0.0); cpgsci(maxcol); /* Green */
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = std_int_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = std_int_med + sigma * std_int_std;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = std_int_med - sigma * std_int_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, std_int, times);
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
    yl = std_chan_med - 2.0 * sigma * std_chan_std;
    yh = std_chan_med + 2.0 * sigma * std_chan_std;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgmtxt("L", 0.8, 0.5, 0.5, "Data \\gs");
    cpgscr(maxcol, 0.0, 1.0, 0.0); cpgsci(maxcol); /* Green */
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = std_chan_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = std_chan_med + sigma * std_chan_std;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = std_chan_med - sigma * std_chan_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, std_chan);
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
    locut = dataavg_med - sigma * dataavg_std;
    hicut = dataavg_med + sigma * dataavg_std;
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
    xl = avg_int_med - 2.0 * sigma * avg_int_std;
    xh = avg_int_med + 2.0 * sigma * avg_int_std;
    yl = 0.0;
    yh = T;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BCST", 0.0, 0, "BST", 0.0, 0);
    cpgmtxt("B", 1.2, 0.5, 0.5, "Data Mean");
    cpgscr(maxcol, 0.0, 0.0, 1.0); cpgsci(maxcol); /* Blue */
    yarr[0] = yl; yarr[1] = yh;
    xarr[0] = xarr[1] = avg_int_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    xarr[0] = xarr[1] = avg_int_med + sigma * avg_int_std;
    cpgline(2, xarr, yarr);
    xarr[0] = xarr[1] = avg_int_med - sigma * avg_int_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numint, avg_int, times);
    yl = 0.0;
    yh = numint;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("", 0.0, 0, "CMST", 0.0, 0);
    /* cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number"); */

    /*  Max Average versus Channel */

    left   = lm + 2.0 * ft + 4.0 * tt;
    right  = lm + 3.0 * ft + 4.0 * tt;
    bottom = bm + fh;
    top    = bm + fh + th;
    cpgsvp (left, right, bottom, top);
    xl = 0.0;
    xh = numchan;
    yl = avg_chan_med - 2.0 * sigma * avg_chan_std;
    yh = avg_chan_med + 2.0 * sigma * avg_chan_std;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("BST", 0.0, 0, "BCST", 0.0, 0);
    cpgmtxt("L", 0.8, 0.5, 0.5, "Data Mean");
    cpgscr(maxcol, 0.0, 0.0, 1.0); cpgsci(maxcol); /* Blue */
    xarr[0] = xl; xarr[1] = xh;
    yarr[0] = yarr[1] = avg_chan_med;
    cpgline(2, xarr, yarr);
    cpgsls(4); /* Dotted line */
    yarr[0] = yarr[1] = avg_chan_med + sigma * avg_chan_std;
    cpgline(2, xarr, yarr);
    yarr[0] = yarr[1] = avg_chan_med - sigma * avg_chan_std;
    cpgline(2, xarr, yarr);
    cpgsls(1); /* Solid line */
    cpgsci(1); /* Default color */
    cpgline(numchan, chans, avg_chan);
    xl = lof;
    xh = hif;
    cpgswin(xl, xh, yl, yh);
    cpgbox ("CMST", 0.0, 0, "", 0.0, 0);
    cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

    /* Add the Data Info area */

    {
      char out[200], out2[100];
      float dy = 0.025;
      
      cpgsvp (0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      left = lm + ft + 2.0 * tt;
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
      cpgtext(left+0.20,   top-2*dy, out);
      sprintf(out, "= %-.3g", datapow_std);
      cpgtext(left+0.23,   top-2*dy, out);
      find_min_max_arr(numint*numchan, datapow[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-3*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-3*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.20,   top-3*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.23,   top-3*dy, out);
      sprintf(out, "Data \\gs:");
      cpgtext(left+0.0,  top-4*dy, out);
      sprintf(out, "median");
      cpgtext(left+0.06,  top-4*dy, out);
      sprintf(out, "= %-.3f", datastd_med);
      cpgtext(left+0.12,  top-4*dy, out);
      sprintf(out, "\\gs");
      cpgtext(left+0.20,   top-4*dy, out);
      sprintf(out, "= %-.3g", datastd_std);
      cpgtext(left+0.23,   top-4*dy, out);
      find_min_max_arr(numint*numchan, datastd[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-5*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-5*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.20,   top-5*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.23,   top-5*dy, out);
      sprintf(out, "Mean:");
      cpgtext(left+0.0,  top-6*dy, out);
      sprintf(out, "median");
      cpgtext(left+0.06,  top-6*dy, out);
      sprintf(out, "= %-.3f", dataavg_med);
      cpgtext(left+0.12,  top-6*dy, out);
      sprintf(out, "\\gs");
      cpgtext(left+0.20,   top-6*dy, out);
      sprintf(out, "= %-.3g", dataavg_std);
      cpgtext(left+0.23,   top-6*dy, out);
      find_min_max_arr(numint*numchan, dataavg[0], &min, &max);
      sprintf(out, "min");
      cpgtext(left+0.06,  top-7*dy, out);
      sprintf(out, "= %-.3f", min);
      cpgtext(left+0.12,  top-7*dy, out);
      sprintf(out, "max");
      cpgtext(left+0.20,   top-7*dy, out);
      sprintf(out, "= %-.3f", max);
      cpgtext(left+0.23,   top-7*dy, out);
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
