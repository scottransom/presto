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
    return tmp * numpowers * pow(1.0 - tmp, numpowers - 1);
}

static void plot_rfi(rfi * plotrfi, float top, int numint, int numchan,
                     float T, float lof, float hif)
{
    int ii;
    float period, perioderr, dy = 0.035, *temparr;
    float tr[6] = { -0.5, 1.0, 0.0, -0.5, 0.0, 1.0 };
    char temp[40];

    if (plotrfi->freq_avg == 0.0)
        period = 0.0;
    else
        period = 1000.0 / plotrfi->freq_avg;
    if (plotrfi->freq_var == 0.0)       /* Why are these zero? */
        perioderr = 0.0;
    else
        perioderr = 1000.0 * sqrt(plotrfi->freq_var) /
            (plotrfi->freq_avg * plotrfi->freq_avg);
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.0);
    cpgnice_output_2(temp, plotrfi->freq_avg, sqrt(plotrfi->freq_var), 0);
    cpgptxt(0.03, top - 0.6 * dy, 0.0, 0.0, temp);
    cpgnice_output_2(temp, period, perioderr, 0);
    cpgptxt(0.12, top - 0.6 * dy, 0.0, 0.0, temp);
    sprintf(temp, "%-5.2f", plotrfi->sigma_avg);
    cpgptxt(0.21, top - 0.6 * dy, 0.0, 0.0, temp);
    sprintf(temp, "%d", plotrfi->numobs);
    cpgptxt(0.27, top - 0.6 * dy, 0.0, 0.0, temp);
    ii = (numint > numchan) ? numint : numchan;
    temparr = gen_fvect(ii);
    for (ii = 0; ii < numchan; ii++)
        temparr[ii] = GET_BIT(plotrfi->chans, ii);
    cpgsvp(0.33, 0.64, top - dy, top);
    cpgswin(0.0, numchan, 0.0, 1.0);
    cpgimag(temparr, numchan, 1, 1, numchan, 1, 1, 0.0, 1.0, tr);
    cpgswin(0.0, numchan, 0.0, 1.0);
    cpgbox("BST", 0.0, 0, "BC", 0.0, 0);
    cpgswin(lof, hif, 0.0, 1.0);
    cpgbox("CST", 0.0, 0, "", 0.0, 0);
    for (ii = 0; ii < numint; ii++)
        temparr[ii] = GET_BIT(plotrfi->times, ii);
    cpgsvp(0.65, 0.96, top - dy, top);
    cpgswin(0.0, numint, 0.0, 1.0);
    cpgimag(temparr, numint, 1, 1, numint, 1, 1, 0.0, 1.0, tr);
    cpgswin(0.0, numint, 0.0, 1.0);
    cpgbox("BST", 0.0, 0, "BC", 0.0, 0);
    cpgswin(0.0, T, 0.0, 1.0);
    cpgbox("CST", 0.0, 0, "", 0.0, 0);
    vect_free(temparr);
}


void rfifind_plot(int numchan, int numint, int ptsperint,
                  float timesigma, float freqsigma,
                  float inttrigfrac, float chantrigfrac,
                  float **dataavg, float **datastd, float **datapow,
                  int *userchan, int numuserchan,
                  int *userints, int numuserints,
                  infodata * idata, unsigned char **bytemask,
                  mask * oldmask, mask * newmask,
                  rfi * rfivect, int numrfi, int rfixwin, int rfips, int xwin)
/* Make the beautiful multi-page rfifind plots */
{
    int ii, jj, ct, loops = 1;
    float *freqs, *chans, *times, *ints;
    float *avg_chan_avg, *std_chan_avg, *pow_chan_avg;
    float *avg_chan_med, *std_chan_med, *pow_chan_med;
    float *avg_chan_std, *std_chan_std, *pow_chan_std;
    float *avg_int_avg, *std_int_avg, *pow_int_avg;
    float *avg_int_med, *std_int_med, *pow_int_med;
    float *avg_int_std, *std_int_std, *pow_int_std;
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
    for (ii = 0; ii < numchan; ii++) {
        chans[ii] = ii;
        freqs[ii] = idata->freq + ii * idata->chan_wid;
    }
    ints = gen_fvect(numint);
    times = gen_fvect(numint);
    for (ii = 0; ii < numint; ii++) {
        ints[ii] = ii;
        times[ii] = 0.0 + ii * inttim;
    }

    /* Calculate the statistics of the full set */

    ct = numchan * numint;
    calc_avgmedstd(dataavg[0], ct, 0.8, 1, &dataavg_avg, &dataavg_med, &dataavg_std);
    calc_avgmedstd(datastd[0], ct, 0.8, 1, &datastd_avg, &datastd_med, &datastd_std);
    calc_avgmedstd(datapow[0], ct, 0.5, 1, &datapow_avg, &datapow_med, &datapow_std);
    avg_reject = timesigma * dataavg_std;
    std_reject = timesigma * datastd_std;
    pow_reject = power_for_sigma(freqsigma, 1, ptsperint / 2);

    /* Calculate the channel/integration statistics vectors */

    for (ii = 0; ii < numint; ii++) {
        calc_avgmedstd(dataavg[0] + ii * numchan, numchan, 0.8, 1,
                       avg_int_avg + ii, avg_int_med + ii, avg_int_std + ii);
        calc_avgmedstd(datastd[0] + ii * numchan, numchan, 0.8, 1,
                       std_int_avg + ii, std_int_med + ii, std_int_std + ii);
        calc_avgmedstd(datapow[0] + ii * numchan, numchan, 0.5, 1,
                       pow_int_avg + ii, pow_int_med + ii, pow_int_std + ii);
    }
    for (ii = 0; ii < numchan; ii++) {
        calc_avgmedstd(dataavg[0] + ii, numint, 0.8, numchan,
                       avg_chan_avg + ii, avg_chan_med + ii, avg_chan_std + ii);
        calc_avgmedstd(datastd[0] + ii, numint, 0.8, numchan,
                       std_chan_avg + ii, std_chan_med + ii, std_chan_std + ii);
        calc_avgmedstd(datapow[0] + ii, numint, 0.5, numchan,
                       pow_chan_avg + ii, pow_chan_med + ii, pow_chan_std + ii);
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
        for (ii = 0; ii < numuserints; ii++)
            if (userints[ii] >= 0 && userints[ii] < numint)
                for (jj = 0; jj < numchan; jj++)
                    bytemask[userints[ii]][jj] |= USERINTS;
    if (numuserchan)
        for (ii = 0; ii < numuserchan; ii++)
            if (userchan[ii] >= 0 && userchan[ii] < numchan)
                for (jj = 0; jj < numint; jj++)
                    bytemask[jj][userchan[ii]] |= USERCHAN;

    /* Compare each point in an interval (or channel) with   */
    /* the interval's (or channel's) median and the overall  */
    /* standard deviation.  If the channel/integration       */
    /* medians are more than sigma different than the global */
    /* value, set them to the global.                        */
    {
        float int_med, chan_med;

        for (ii = 0; ii < numint; ii++) {
            for (jj = 0; jj < numchan; jj++) {
                {               /* Powers */
                    if (datapow[ii][jj] > pow_reject)
                        if (!(bytemask[ii][jj] & PADDING))
                            bytemask[ii][jj] |= BAD_POW;
                }
                {               /* Averages */
                    if (fabs(avg_int_med[ii] - dataavg_med) >
                        timesigma * dataavg_std)
                        int_med = dataavg_med;
                    else
                        int_med = avg_int_med[ii];
                    if (fabs(avg_chan_med[jj] - dataavg_med) >
                        timesigma * dataavg_std)
                        chan_med = dataavg_med;
                    else
                        chan_med = avg_chan_med[jj];
                    if (fabs(dataavg[ii][jj] - int_med) > avg_reject ||
                        fabs(dataavg[ii][jj] - chan_med) > avg_reject)
                        if (!(bytemask[ii][jj] & PADDING))
                            bytemask[ii][jj] |= BAD_AVG;
                }
                {               /* Standard Deviations */
                    if (fabs(std_int_med[ii] - datastd_med) >
                        timesigma * datastd_std)
                        int_med = datastd_med;
                    else
                        int_med = std_int_med[ii];
                    if (fabs(std_chan_med[jj] - datastd_med) >
                        timesigma * datastd_std)
                        chan_med = datastd_med;
                    else
                        chan_med = std_chan_med[jj];
                    if (fabs(datastd[ii][jj] - int_med) > std_reject ||
                        fabs(datastd[ii][jj] - chan_med) > std_reject)
                        if (!(bytemask[ii][jj] & PADDING))
                            bytemask[ii][jj] |= BAD_STD;
                }
            }
        }
    }

    /* Step over the intervals and channels and count how many are set "bad". */
    /* For a given interval, if the number of bad channels is greater than    */
    /* chantrigfrac*numchan then reject the whole interval.                   */
    /* For a given channel, if the number of bad intervals is greater than    */
    /* inttrigfrac*numint then reject the whole channel.                      */
    {
        int badnum, trignum;

        /* Loop over the intervals */
        trignum = (int) (numchan * chantrigfrac);
        for (ii = 0; ii < numint; ii++) {
            if (!(bytemask[ii][0] & USERINTS)) {
                badnum = 0;
                for (jj = 0; jj < numchan; jj++)
                    if (bytemask[ii][jj] & BADDATA)
                        badnum++;
                if (badnum > trignum) {
                    userints[numuserints++] = ii;
                    for (jj = 0; jj < numchan; jj++)
                        bytemask[ii][jj] |= USERINTS;
                }
            }
        }

        /* Loop over the channels */
        trignum = (int) (numint * inttrigfrac);
        for (ii = 0; ii < numchan; ii++) {
            if (!(bytemask[0][ii] & USERCHAN)) {
                badnum = 0;
                for (jj = 0; jj < numint; jj++)
                    if (bytemask[jj][ii] & BADDATA)
                        badnum++;
                if (badnum > trignum) {
                    userchan[numuserchan++] = ii;
                    for (jj = 0; jj < numint; jj++)
                        bytemask[jj][ii] |= USERCHAN;
                }
            }
        }
    }

    /* Generate the New Mask */

    fill_mask(timesigma, freqsigma, idata->mjd_i + idata->mjd_f,
              ptsperint * idata->dt, idata->freq, idata->chan_wid,
              numchan, numint, ptsperint, numuserchan, userchan,
              numuserints, userints, bytemask, newmask);

    /* Place the oldmask over the newmask for plotting purposes */

    if (oldmask->numchan)
        set_oldmask_bits(oldmask, bytemask);

    /*
     *  Now plot the results
     */

    if (xwin)
        loops = 2;
    for (ct = 0; ct < loops; ct++) {    /* PS/XWIN Plot Loop */
        float min, max, tr[6], locut, hicut;
        float left, right, top, bottom;
        float xl, xh, yl, yh;
        float tt, ft, th, fh;   /* thin and fat thicknesses and heights */
        float lm, rm, tm, bm;   /* LRTB margins */
        float xarr[2], yarr[2];
        char outdev[300];
        int ii, mincol, maxcol, numcol;

        /*Set the PGPLOT device to an X-Window */

        if (ct == 1)
            strcpy(outdev, "/XWIN");
        else
            sprintf(outdev, "%s.ps/CPS", idata->name);

        /* Open and prep our device */

        cpgopen(outdev);
        cpgpap(10.25, 8.5 / 11.0);
        cpgpage();
        cpgiden();
        cpgsch(0.7);
        cpgqcir(&mincol, &maxcol);
        numcol = maxcol - mincol + 1;
        for (ii = mincol; ii <= maxcol; ii++) {
            float color;
            color = (float) (maxcol - ii) / (float) numcol;
            cpgscr(ii, color, color, color);
        }

        /* Set thicknesses and margins */

        lm = 0.04;
        rm = 0.04;
        bm = 0.08;
        tm = 0.05;
        ft = 3.0;               /* This sets fat thickness = 3 x thin thickness */
        tt = 0.92 / (6.0 + 4.0 * ft);
        ft *= tt;
        fh = 0.55;
        th = tt * 11.0 / 8.5;

        {                       /* Powers Histogram */
            float *theo, *hist, *hpows, *tpows, maxhist = 0.0, maxtheo = 0.0;
            int numhist = 40, numtheo = 200, bin, numpows;
            double dtheo, dhist, spacing;

            /* Calculate the predicted distribution of max powers */

            numpows = numint * numchan;
            {
                int badpows = 0;
                for (ii = 0; ii < numpows; ii++) {
                    if (!isnormal(*(datapow[0]+ii))) {
                        if (*(datapow[0]+ii) != 0.0) {
                            printf("WARNING:  bad power (%f) at int=%d, chan=%d\n",
                                   *(datapow[0]+ii), ii/numchan, ii%numchan);
                            *(datapow[0]+ii)=0.0;
                            badpows++;
                        }
                    }
                }
                if (badpows)
                    printf("WARNING:  Found %d bad powers in the datapow array.  "
                           "Zeroing them out.\n", badpows);
            }
            find_min_max_arr(numpows, datapow[0], &min, &max);
            min = (min < 5.0) ? log10(5.0 * 0.95) : log10(min * 0.95);
            max = log10(max * 1.05);
            dhist = (max - min) / numhist;
            theo = gen_fvect(numtheo);
            tpows = gen_fvect(numtheo);
            hist = gen_fvect(numhist);
            hpows = gen_fvect(numhist);
            for (ii = 0; ii < numhist; ii++) {
                hist[ii] = 0.0;
                hpows[ii] = min + ii * dhist;
            }
            for (ii = 0; ii < numpows; ii++) {
                bin = (*(datapow[0] + ii) == 0.0) ? 0 :
                    (log10(*(datapow[0] + ii)) - min) / dhist;
                if (bin < 0)
                    bin = 0;
                if (bin >= numhist)
                    bin = numhist - 1;
                hist[bin] += 1.0;
            }
            for (ii = 0; ii < numhist; ii++)
                if (hist[ii] > maxhist)
                    maxhist = hist[ii];
            maxhist *= 1.1;
            dtheo = (max - min) / (double) (numtheo - 1);
            for (ii = 0; ii < numtheo; ii++) {
                tpows[ii] = min + ii * dtheo;
                theo[ii] = single_power_pdf(pow(10.0, tpows[ii]),
                                            ptsperint / 2) * numpows;
                spacing = (pow(10.0, tpows[ii] + dhist) - pow(10.0, tpows[ii]));
                theo[ii] *= spacing;
                if (theo[ii] > maxtheo)
                    maxtheo = theo[ii];
            }
            maxtheo *= 1.1;
            if (maxtheo > maxhist)
                maxhist = maxtheo;
            left = lm;
            right = lm + ft + tt;
            bottom = 0.80;
            top = 0.96;
            cpgsvp(left, right, bottom, top);
            xl = min;
            xh = max;
            yl = 0.0;
            yh = maxhist;
            cpgswin(xl, xh, yl, yh);
            cpgmtxt("L", 1.1, 0.5, 0.5, "Number");
            cpgmtxt("B", 2.1, 0.5, 0.5, "Max Power");
            cpgbin(numhist, hpows, hist, 0);
            cpgscr(maxcol, 0.5, 0.5, 0.5);
            cpgsci(maxcol);     /* Grey */
            cpgline(numtheo, tpows, theo);
            xarr[0] = log10(power_for_sigma(freqsigma, 1, ptsperint / 2));
            xarr[1] = xarr[0];
            yarr[0] = yl;
            yarr[1] = yh;
            cpgsls(4);          /* Dotted line */
            cpgscr(maxcol, 1.0, 0.0, 0.0);
            cpgsci(maxcol);     /* Red */
            cpgline(2, xarr, yarr);
            cpgsls(1);          /* Solid line */
            cpgsci(1);          /* Default color */
            cpgbox("BCLNST", 0.0, 0, "BC", 0.0, 0);
            vect_free(hist);
            vect_free(theo);
            vect_free(tpows);
            vect_free(hpows);
        }

        /* Maximum Powers */

        left = lm;
        right = lm + ft;
        bottom = bm;
        top = bm + fh;
        xl = 0.0;
        xh = numchan;
        yl = 0.0;
        yh = T;
        cpgsvp(left, right, bottom, top);
        cpgswin(xl, xh, yl, yh);
        cpgscr(maxcol, 1.0, 0.0, 0.0);  /* Red */
        locut = 0.0;
        hicut = pow_reject;
        tr[2] = tr[4] = 0.0;
        tr[1] = (xh - xl) / numchan;
        tr[0] = xl - (tr[1] / 2);
        tr[5] = (yh - yl) / numint;
        tr[3] = yl - (tr[5] / 2);
        cpgimag(datapow[0], numchan, numint, 1, numchan, 1, numint, locut, hicut,
                tr);
        cpgswin(xl, xh, yl, yh);
        cpgbox("BNST", 0.0, 0, "BNST", 0.0, 0);
        cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
        cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");
        xl = lof;
        xh = hif;
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CST", 0.0, 0, "CST", 0.0, 0);

        /* Max Power Label */

        left = lm + ft;
        right = lm + ft + tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);         /* Red */
        cpgptxt(0.5, 0.7, 0.0, 0.5, "Max");
        cpgptxt(0.5, 0.3, 0.0, 0.5, "Power");
        cpgsci(1);              /* Default color */

        /*  Max Power versus Time */

        left = lm + ft;
        right = lm + ft + tt;
        bottom = bm;
        top = bm + fh;
        cpgsvp(left, right, bottom, top);
        find_min_max_arr(numint, pow_int_med, &min, &max);
        xl = 0.0;
        xh = 1.5 * pow_reject;
        yl = 0.0;
        yh = T;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCST", 0.0, 0, "BST", 0.0, 0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);         /* Red */
        yarr[0] = yl;
        yarr[1] = yh;
        xarr[0] = xarr[1] = datapow_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        xarr[0] = xarr[1] = pow_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numint, pow_int_med, times);
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("", 0.0, 0, "CMST", 0.0, 0);
        /* cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number"); */

        /*  Max Power versus Channel */

        left = lm;
        right = lm + ft;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        find_min_max_arr(numchan, pow_chan_med, &min, &max);
        xl = 0.0;
        xh = numchan;
        yl = 0.0;
        yh = 1.5 * pow_reject;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BST", 0.0, 0, "BCST", 0.0, 0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);         /* Red */
        xarr[0] = xl;
        xarr[1] = xh;
        yarr[0] = yarr[1] = datapow_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        yarr[0] = yarr[1] = pow_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numchan, chans, pow_chan_med);
        xl = lof;
        xh = hif;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CMST", 0.0, 0, "", 0.0, 0);
        cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

        /* Standard Deviations */

        left = lm + ft + 2.0 * tt;
        right = lm + 2.0 * ft + 2.0 * tt;
        bottom = bm;
        top = bm + fh;
        xl = 0.0;
        xh = numchan;
        yl = 0.0;
        yh = T;
        cpgsvp(left, right, bottom, top);
        cpgswin(xl, xh, yl, yh);
        cpgscr(mincol, 0.7, 1.0, 0.7);  /* Light Green */
        cpgscr(maxcol, 0.3, 1.0, 0.3);  /* Dark Green */
        locut = datastd_med - timesigma * datastd_std;
        hicut = datastd_med + timesigma * datastd_std;
        tr[2] = tr[4] = 0.0;
        tr[1] = (xh - xl) / numchan;
        tr[0] = xl - (tr[1] / 2);
        tr[5] = (yh - yl) / numint;
        tr[3] = yl - (tr[5] / 2);
        cpgimag(datastd[0], numchan, numint, 1, numchan, 1, numint, locut, hicut,
                tr);
        cpgswin(xl, xh, yl, yh);
        cpgbox("BNST", 0.0, 0, "BNST", 0.0, 0);
        cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
        xl = lof;
        xh = hif;
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CST", 0.0, 0, "CST", 0.0, 0);

        /* Data Sigma Label */

        left = lm + 2.0 * ft + 2.0 * tt;
        right = lm + 2.0 * ft + 3.0 * tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        cpgscr(maxcol, 0.0, 1.0, 0.0);
        cpgsci(maxcol);         /* Green */
        cpgptxt(0.5, 0.7, 0.0, 0.5, "Data");
        cpgptxt(0.5, 0.3, 0.0, 0.5, "Sigma");
        cpgsci(1);              /* Default color */

        /*  Data Sigma versus Time */

        left = lm + 2.0 * ft + 2.0 * tt;
        right = lm + 2.0 * ft + 3.0 * tt;
        bottom = bm;
        top = bm + fh;
        cpgsvp(left, right, bottom, top);
        xl = datastd_med - 2.0 * std_reject;
        xh = datastd_med + 2.0 * std_reject;
        yl = 0.0;
        yh = T;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCST", 0.0, 0, "BST", 0.0, 0);
        cpgscr(maxcol, 0.0, 1.0, 0.0);
        cpgsci(maxcol);         /* Green */
        yarr[0] = yl;
        yarr[1] = yh;
        xarr[0] = xarr[1] = datastd_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        xarr[0] = xarr[1] = datastd_med + std_reject;
        cpgline(2, xarr, yarr);
        xarr[0] = xarr[1] = datastd_med - std_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numint, std_int_med, times);
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("", 0.0, 0, "CMST", 0.0, 0);
        /* cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number"); */

        /*  Data Sigma versus Channel */

        left = lm + ft + 2.0 * tt;
        right = lm + 2.0 * ft + 2.0 * tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        xl = 0.0;
        xh = numchan;
        yl = datastd_med - 2.0 * std_reject;
        yh = datastd_med + 2.0 * std_reject;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BST", 0.0, 0, "BCST", 0.0, 0);
        cpgscr(maxcol, 0.0, 1.0, 0.0);
        cpgsci(maxcol);         /* Green */
        xarr[0] = xl;
        xarr[1] = xh;
        yarr[0] = yarr[1] = datastd_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        yarr[0] = yarr[1] = datastd_med + std_reject;
        cpgline(2, xarr, yarr);
        yarr[0] = yarr[1] = datastd_med - std_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numchan, chans, std_chan_med);
        xl = lof;
        xh = hif;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CMST", 0.0, 0, "", 0.0, 0);
        cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

        /* Data Mean */

        left = lm + 2.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 4.0 * tt;
        bottom = bm;
        top = bm + fh;
        xl = 0.0;
        xh = numchan;
        yl = 0.0;
        yh = T;
        cpgsvp(left, right, bottom, top);
        cpgswin(xl, xh, yl, yh);
        cpgscr(mincol, 0.7, 0.7, 1.0);  /* Light Blue */
        cpgscr(maxcol, 0.3, 0.3, 1.0);  /* Dark Blue */
        locut = dataavg_med - timesigma * dataavg_std;
        hicut = dataavg_med + timesigma * dataavg_std;
        tr[2] = tr[4] = 0.0;
        tr[1] = (xh - xl) / numchan;
        tr[0] = xl - (tr[1] / 2);
        tr[5] = (yh - yl) / numint;
        tr[3] = yl - (tr[5] / 2);
        cpgimag(dataavg[0], numchan, numint, 1, numchan, 1, numint, locut, hicut,
                tr);
        cpgswin(xl, xh, yl, yh);
        cpgbox("BNST", 0.0, 0, "BNST", 0.0, 0);
        cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
        xl = lof;
        xh = hif;
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CST", 0.0, 0, "CST", 0.0, 0);

        /* Data Mean Label */

        left = lm + 3.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 5.0 * tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        cpgscr(maxcol, 0.0, 0.0, 1.0);
        cpgsci(maxcol);         /* Blue */
        cpgptxt(0.5, 0.7, 0.0, 0.5, "Data");
        cpgptxt(0.5, 0.3, 0.0, 0.5, "Mean");
        cpgsci(1);              /* Default color */

        /*  Data Mean versus Time */

        left = lm + 3.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 5.0 * tt;
        bottom = bm;
        top = bm + fh;
        cpgsvp(left, right, bottom, top);
        xl = dataavg_med - 2.0 * avg_reject;
        xh = dataavg_med + 2.0 * avg_reject;
        yl = 0.0;
        yh = T;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCST", 0.0, 0, "BST", 0.0, 0);
        cpgscr(maxcol, 0.0, 0.0, 1.0);
        cpgsci(maxcol);         /* Blue */
        yarr[0] = yl;
        yarr[1] = yh;
        xarr[0] = xarr[1] = dataavg_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        xarr[0] = xarr[1] = dataavg_med + avg_reject;
        cpgline(2, xarr, yarr);
        xarr[0] = xarr[1] = dataavg_med - avg_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numint, avg_int_med, times);
        yl = 0.0;
        yh = numint;
        cpgswin(xl, xh, yl, yh);
        cpgbox("", 0.0, 0, "CMST", 0.0, 0);

        /*  Data Mean versus Channel */

        left = lm + 2.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 4.0 * tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        xl = 0.0;
        xh = numchan;
        yl = dataavg_med - 2.0 * avg_reject;
        yh = dataavg_med + 2.0 * avg_reject;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BST", 0.0, 0, "BCST", 0.0, 0);
        cpgscr(maxcol, 0.0, 0.0, 1.0);
        cpgsci(maxcol);         /* Blue */
        xarr[0] = xl;
        xarr[1] = xh;
        yarr[0] = yarr[1] = dataavg_med;
        cpgline(2, xarr, yarr);
        cpgsls(4);              /* Dotted line */
        yarr[0] = yarr[1] = dataavg_med + avg_reject;
        cpgline(2, xarr, yarr);
        yarr[0] = yarr[1] = dataavg_med - avg_reject;
        cpgline(2, xarr, yarr);
        cpgsls(1);              /* Solid line */
        cpgsci(1);              /* Default color */
        cpgline(numchan, chans, avg_chan_med);
        xl = lof;
        xh = hif;
        cpgswin(xl, xh, yl, yh);
        cpgbox("CMST", 0.0, 0, "", 0.0, 0);
        cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");

        {                       /* Add the Data Info area */
            char out[310], out2[300];
            float dy = 0.025;

            cpgsvp(0.0, 1.0, 0.0, 1.0);
            cpgswin(0.0, 1.0, 0.0, 1.0);
            left = lm + ft + 1.5 * tt;
            top = 1.0 - tm;
            cpgsch(1.0);
            sprintf(out, "%-s", idata->name);
            cpgptxt(0.5, 1.0 - 0.5 * tm, 0.0, 0.5, out);
            cpgsch(0.8);

            sprintf(out, "Object:");
            cpgtext(left + 0.0, top - 0 * dy, out);
            sprintf(out, "%-s", idata->object);
            cpgtext(left + 0.1, top - 0 * dy, out);
            sprintf(out, "Telescope:");
            cpgtext(left + 0.0, top - 1 * dy, out);
            sprintf(out, "%-s", idata->telescope);
            cpgtext(left + 0.1, top - 1 * dy, out);
            sprintf(out, "Instrument:");
            cpgtext(left + 0.0, top - 2 * dy, out);
            sprintf(out, "%-s", idata->instrument);
            cpgtext(left + 0.1, top - 2 * dy, out);
            ra_dec_to_string(out2, idata->ra_h, idata->ra_m, idata->ra_s);
            sprintf(out, "RA\\dJ2000\\u");
            cpgtext(left + 0.0, top - 3 * dy, out);
            sprintf(out, "= %-s", out2);
            cpgtext(left + 0.08, top - 3 * dy, out);
            ra_dec_to_string(out2, idata->dec_d, idata->dec_m, idata->dec_s);
            sprintf(out, "DEC\\dJ2000\\u");
            cpgtext(left + 0.0, top - 4 * dy, out);
            sprintf(out, "= %-s", out2);
            cpgtext(left + 0.08, top - 4 * dy, out);
            sprintf(out, "Epoch\\dtopo\\u");
            cpgtext(left + 0.0, top - 5 * dy, out);
            sprintf(out, "= %-.11f", idata->mjd_i + idata->mjd_f);
            cpgtext(left + 0.08, top - 5 * dy, out);
            sprintf(out, "T\\dsample\\u (s)");
            cpgtext(left + 0.0, top - 6 * dy, out);
            sprintf(out, "= %g", idata->dt);
            cpgtext(left + 0.08, top - 6 * dy, out);
            sprintf(out, "T\\dtotal\\u (s)");
            cpgtext(left + 0.0, top - 7 * dy, out);
            sprintf(out, "= %g", T);
            cpgtext(left + 0.08, top - 7 * dy, out);

            left = lm + ft + 7.8 * tt;
            sprintf(out, "Num channels");
            cpgtext(left + 0.0, top - 0 * dy, out);
            sprintf(out, "= %-d", numchan);
            cpgtext(left + 0.12, top - 0 * dy, out);
            sprintf(out, "Pts per int");
            cpgtext(left + 0.19, top - 0 * dy, out);
            sprintf(out, "= %-d", ptsperint);
            cpgtext(left + 0.29, top - 0 * dy, out);
            sprintf(out, "Num intervals");
            cpgtext(left + 0.0, top - 1 * dy, out);
            sprintf(out, "= %-d", numint);
            cpgtext(left + 0.12, top - 1 * dy, out);
            sprintf(out, "Time per int");
            cpgtext(left + 0.19, top - 1 * dy, out);
            sprintf(out, "= %-g", inttim);
            cpgtext(left + 0.29, top - 1 * dy, out);
            sprintf(out, "Power:");
            cpgtext(left + 0.0, top - 2 * dy, out);
            sprintf(out, "median");
            cpgtext(left + 0.06, top - 2 * dy, out);
            sprintf(out, "= %-.3f", datapow_med);
            cpgtext(left + 0.12, top - 2 * dy, out);
            sprintf(out, "\\gs");
            cpgtext(left + 0.21, top - 2 * dy, out);
            sprintf(out, "= %-.3g", datapow_std);
            cpgtext(left + 0.245, top - 2 * dy, out);
            find_min_max_arr(numint * numchan, datapow[0], &min, &max);
            sprintf(out, "min");
            cpgtext(left + 0.06, top - 3 * dy, out);
            sprintf(out, "= %-.3f", min);
            cpgtext(left + 0.12, top - 3 * dy, out);
            sprintf(out, "max");
            cpgtext(left + 0.21, top - 3 * dy, out);
            sprintf(out, "= %-.3f", max);
            cpgtext(left + 0.245, top - 3 * dy, out);
            sprintf(out, "Sigma:");
            cpgtext(left + 0.0, top - 4 * dy, out);
            sprintf(out, "median");
            cpgtext(left + 0.06, top - 4 * dy, out);
            sprintf(out, "= %-.3f", datastd_med);
            cpgtext(left + 0.12, top - 4 * dy, out);
            sprintf(out, "\\gs");
            cpgtext(left + 0.21, top - 4 * dy, out);
            sprintf(out, "= %-.3g", datastd_std);
            cpgtext(left + 0.245, top - 4 * dy, out);
            find_min_max_arr(numint * numchan, datastd[0], &min, &max);
            sprintf(out, "min");
            cpgtext(left + 0.06, top - 5 * dy, out);
            sprintf(out, "= %-.3f", min);
            cpgtext(left + 0.12, top - 5 * dy, out);
            sprintf(out, "max");
            cpgtext(left + 0.21, top - 5 * dy, out);
            sprintf(out, "= %-.3f", max);
            cpgtext(left + 0.245, top - 5 * dy, out);
            sprintf(out, "Mean:");
            cpgtext(left + 0.0, top - 6 * dy, out);
            sprintf(out, "median");
            cpgtext(left + 0.06, top - 6 * dy, out);
            sprintf(out, "= %-.3f", dataavg_med);
            cpgtext(left + 0.12, top - 6 * dy, out);
            sprintf(out, "\\gs");
            cpgtext(left + 0.21, top - 6 * dy, out);
            sprintf(out, "= %-.3g", dataavg_std);
            cpgtext(left + 0.245, top - 6 * dy, out);
            find_min_max_arr(numint * numchan, dataavg[0], &min, &max);
            sprintf(out, "min");
            cpgtext(left + 0.06, top - 7 * dy, out);
            sprintf(out, "= %-.3f", min);
            cpgtext(left + 0.12, top - 7 * dy, out);
            sprintf(out, "max");
            cpgtext(left + 0.21, top - 7 * dy, out);
            sprintf(out, "= %-.3f", max);
            cpgtext(left + 0.245, top - 7 * dy, out);
        }

        {                       /* Plot the Mask */
            unsigned char byte;
            char temp[300];
            float **plotmask, rr, gg, bb, page;

            plotmask = gen_fmatrix(numint, numchan);
            for (ii = 0; ii < numint; ii++) {
                for (jj = 0; jj < numchan; jj++) {
                    byte = bytemask[ii][jj];
                    plotmask[ii][jj] = 0.0;
                    if (byte & PADDING)
                        plotmask[ii][jj] = 1.0;
                    if (byte & OLDMASK)
                        plotmask[ii][jj] = 2.0;
                    if (byte & USERZAP)
                        plotmask[ii][jj] = 3.0;
                    if (byte & BAD_POW)
                        plotmask[ii][jj] = 4.0;
                    else if (byte & BAD_AVG)
                        plotmask[ii][jj] = 5.0;
                    else if (byte & BAD_STD)
                        plotmask[ii][jj] = 6.0;
                }
            }
            /* Set the colors */
            numcol = 7;
            maxcol = mincol + numcol - 1;
            cpgscir(mincol, maxcol);
            cpgqcr(0, &rr, &gg, &bb);
            cpgscr(mincol + 0, rr, gg, bb);     /* GOODDATA = background */
            cpgscr(mincol + 1, 0.7, 0.7, 0.7);  /* PADDING  = light grey */
            cpgscr(mincol + 2, 0.3, 0.3, 0.3);  /* OLDMASK  = dark grey */
            cpgqcr(1, &rr, &gg, &bb);
            cpgscr(mincol + 3, rr, gg, bb);     /* USERZAP  = foreground */
            cpgscr(mincol + 4, 1.0, 0.0, 0.0);  /* BAD+POW  = red */
            cpgscr(mincol + 5, 0.0, 0.0, 1.0);  /* BAD+AVG  = blue */
            cpgscr(mincol + 6, 0.0, 1.0, 0.0);  /* BAD+STD  = green */
            /* Prep the image */
            for (page = 0; page <= 1; page++) {
                xl = 0.0;
                xh = numchan;
                yl = 0.0;
                yh = T;
                locut = 0.0;
                hicut = 6.0;
                tr[2] = tr[4] = 0.0;
                tr[1] = (xh - xl) / numchan;
                tr[0] = xl - (tr[1] / 2);
                tr[5] = (yh - yl) / numint;
                tr[3] = yl - (tr[5] / 2);
                if (page == 0) {
                    left = lm + 3.0 * ft + 6.0 * tt;
                    right = lm + 4.0 * ft + 6.0 * tt;
                    bottom = bm;
                    top = bm + fh;
                } else {
                    cpgpage();
                    cpgiden();
                    left = 0.06;
                    right = 0.94;
                    bottom = 0.06;
                    top = 0.88;
                }
                cpgsvp(left, right, bottom, top);
                cpgswin(xl, xh, yl, yh);
                cpgimag(plotmask[0], numchan, numint, 1,
                        numchan, 1, numint, locut, hicut, tr);
                cpgswin(xl, xh, yl, yh);
                cpgbox("BNST", 0.0, 0, "BNST", 0.0, 0);
                cpgmtxt("B", 2.6, 0.5, 0.5, "Channel");
                if (page)
                    cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");
                xl = lof;
                xh = hif;
                yl = 0.0;
                yh = numint;
                cpgswin(xl, xh, yl, yh);
                cpgbox("CMST", 0.0, 0, "CMST", 0.0, 0);
                cpgmtxt("T", 1.8, 0.5, 0.5, "Frequency (MHz)");
                cpgmtxt("R", 2.3, 0.5, 0.5, "Interval Number");
                /* Add the Labels */
                cpgsvp(0.0, 1.0, 0.0, 1.0);
                cpgswin(0.0, 1.0, 0.0, 1.0);
                cpgsch(0.8);
                if (page == 0) {
                    cpgsci(mincol + 1);
                    cpgptxt(left, top + 0.1, 0.0, 0.0, "Padding");
                    cpgsci(mincol + 2);
                    cpgptxt(left, top + 0.08, 0.0, 0.0, "Old Mask");
                    cpgsci(mincol + 3);
                    cpgptxt(left, top + 0.06, 0.0, 0.0, "User Zap");
                    cpgsci(mincol + 4);
                    cpgptxt(right, top + 0.1, 0.0, 1.0, "Power");
                    cpgsci(mincol + 6);
                    cpgptxt(right, top + 0.08, 0.0, 1.0, "Sigma");
                    cpgsci(mincol + 5);
                    cpgptxt(right, top + 0.06, 0.0, 1.0, "Mean");
                    cpgsci(1);
                } else {
                    cpgsci(mincol + 1);
                    cpgptxt(1.0 / 12.0, 0.955, 0.0, 0.5, "Padding");
                    cpgsci(mincol + 2);
                    cpgptxt(3.0 / 12.0, 0.955, 0.0, 0.5, "Old Mask");
                    cpgsci(mincol + 3);
                    cpgptxt(5.0 / 12.0, 0.955, 0.0, 0.5, "User Zap");
                    cpgsci(mincol + 4);
                    cpgptxt(7.0 / 12.0, 0.955, 0.0, 0.5, "Max Power");
                    cpgsci(mincol + 6);
                    cpgptxt(9.0 / 12.0, 0.955, 0.0, 0.5, "Data Sigma");
                    cpgsci(mincol + 5);
                    cpgptxt(11.0 / 12.0, 0.955, 0.0, 0.5, "Data Mean");
                    cpgsci(1);
                    cpgsch(0.9);
                    sprintf(temp, "Recommended Mask for '%-s'", idata->name);
                    cpgptxt(0.5, 0.985, 0.0, 0.5, temp);
                }
            }
            vect_free(plotmask[0]);
            vect_free(plotmask);
        }

        if (ct == 0)
            printf("There are %d RFI instances.\n\n", numrfi);

        if ((ct == 0 && rfips) || (ct == 1 && rfixwin)) {       /* Plot the RFI instances */
            int maxcol, mincol, numperpage = 25, numtoplot;
            float dy = 0.035, top = 0.95, rr, gg, bb;
            char temp[300];

            qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_freq);
            /* qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_sigma); */
            for (ii = 0; ii <= (numrfi - 1) / numperpage; ii++) {
                cpgpage();
                cpgiden();
                cpgsvp(0.0, 1.0, 0.0, 1.0);
                cpgswin(0.0, 1.0, 0.0, 1.0);
                cpgsch(0.8);
                sprintf(temp, "%-s", idata->name);
                cpgtext(0.05, 0.985, temp);
                cpgsch(0.6);
                sprintf(temp, "Freq (Hz)");
                cpgptxt(0.03, 0.96, 0.0, 0.0, temp);
                sprintf(temp, "Period (ms)");
                cpgptxt(0.12, 0.96, 0.0, 0.0, temp);
                sprintf(temp, "Sigma");
                cpgptxt(0.21, 0.96, 0.0, 0.0, temp);
                sprintf(temp, "Number");
                cpgptxt(0.27, 0.96, 0.0, 0.0, temp);
                cpgsvp(0.33, 0.64, top - dy, top);
                cpgswin(lof, hif, 0.0, 1.0);
                cpgbox("CIMST", 0.0, 0, "", 0.0, 0);
                cpgmtxt("T", 2.5, 0.5, 0.5, "Frequency (MHz)");
                cpgsvp(0.65, 0.96, top - dy, top);
                cpgswin(0.0, T, 0.0, 1.0);
                cpgbox("CIMST", 0.0, 0, "", 0.0, 0);
                cpgmtxt("T", 2.5, 0.5, 0.5, "Time (s)");
                cpgqcir(&mincol, &maxcol);
                maxcol = mincol + 1;
                cpgscir(mincol, maxcol);
                cpgqcr(0, &rr, &gg, &bb);
                cpgscr(mincol, rr, gg, bb);     /* background */
                cpgqcr(1, &rr, &gg, &bb);
                /* cpgscr(maxcol, rr, gg, bb);  foreground */
                cpgscr(maxcol, 0.5, 0.5, 0.5);  /* grey */
                if (ii == (numrfi - 1) / numperpage)
                    numtoplot = numrfi % numperpage;
                else
                    numtoplot = numperpage;
                for (jj = 0; jj < numtoplot; jj++)
                    plot_rfi(rfivect + ii * numperpage + jj,
                             top - jj * dy, numint, numchan, T, lof, hif);
                cpgsvp(0.33, 0.64, top - jj * dy, top - (jj - 1) * dy);
                cpgswin(0.0, numchan, 0.0, 1.0);
                cpgbox("BINST", 0.0, 0, "", 0.0, 0);
                cpgmtxt("B", 2.5, 0.5, 0.5, "Channel");
                cpgsvp(0.65, 0.96, top - jj * dy, top - (jj - 1) * dy);
                cpgswin(0.0, numint, 0.0, 1.0);
                cpgbox("BINST", 0.0, 0, "", 0.0, 0);
                cpgmtxt("B", 2.5, 0.5, 0.5, "Interval");
            }
        }
        cpgclos();
    }                           /* Plot for loop */

    /* Free our arrays */

    vect_free(freqs);
    vect_free(chans);
    vect_free(times);
    vect_free(ints);
    vect_free(avg_chan_avg);
    vect_free(std_chan_avg);
    vect_free(pow_chan_avg);
    vect_free(avg_int_avg);
    vect_free(std_int_avg);
    vect_free(pow_int_avg);
    vect_free(avg_chan_med);
    vect_free(std_chan_med);
    vect_free(pow_chan_med);
    vect_free(avg_int_med);
    vect_free(std_int_med);
    vect_free(pow_int_med);
    vect_free(avg_chan_std);
    vect_free(std_chan_std);
    vect_free(pow_chan_std);
    vect_free(avg_int_std);
    vect_free(std_int_std);
    vect_free(pow_int_std);
}
