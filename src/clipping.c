#include <presto.h>

#define BLOCKSTOSKIP 100
int read_onoff_paris(FILE * onofffile, long long **onbins, long long **offbins)
// Read a text file containing rows of ASCII pairs of integers
// corresponding to where clipping should be on or off 
// The format for a line is: "CLIPON:CLIPOFF" with 0 indexing
// Lines beginning with '#' are ignored.
{
    char line[200], *sptr = NULL;
    int ii, numpairs;

    /* Read the input file once to count pairs */
    numpairs = 0;
    while (!feof(onofffile)) {
        sptr = fgets(line, 200, onofffile);
        if (sptr != NULL && sptr[0] != '\n') {
            if (line[0] == '#')
                continue;

            else
                numpairs++;
        }
    }
    *onbins = (long long *) malloc(sizeof(long long) * numpairs);
    *offbins = (long long *) malloc(sizeof(long long) * numpairs);

    /* Rewind and read the pairs for real */
    rewind(onofffile);
    ii = 0;
    while (ii < numpairs) {
        sptr = fgets(line, 200, onofffile);
        if (sptr != NULL && sptr[0] != '\n') {
            if (line[0] == '#')
                continue;

            else {
                sscanf(line, "%lld:%lld\n", &(*onbins)[ii], &(*offbins)[ii]);
                ii++;
            }
        }
    }
    return numpairs;
}


/* NEW Clipping Routine (uses channel running averages) */
int clip_times(float *rawdata, int ptsperblk, int numchan, float clip_sigma,
               float *good_chan_levels)
// Perform time-domain clipping of rawdata.  This is a 2D array with
// ptsperblk*numchan points, each of which is a float.  The clipping
// is done at clip_sigma sigma above/below the running mean.  The
// up-to-date running averages of the channels are returned in
// good_chan_levels (which must be pre-allocated).
{
    static float *chan_running_avg;
    static float running_avg = 0.0, running_std = 0.0;
    static int blocksread = 0, firsttime = 1;
    static long long current_point = 0;
    static int numonoff = 0, onoffindex = 0;
    static long long *onbins = NULL, *offbins = NULL;
    float *zero_dm_block, *ftmp, *powptr;
    double *chan_avg_temp;
    float current_med, trigger;
    double current_avg = 0.0, current_std = 0.0;
    int ii, jj, clipit = 0, clipped = 0;
    if (firsttime) {
        chan_running_avg = gen_fvect(numchan);
        firsttime = 0;
        {                       // This is experimental code to zap radar-filled data
            char *envval = getenv("CLIPBINSFILE");
            if (envval != NULL) {
                FILE *onofffile = chkfopen(envval, "r");
                numonoff = read_onoff_paris(onofffile, &onbins, &offbins);
                fclose(onofffile);
                printf("\nRead %d bin clipping pairs from '%s'.\n", numonoff,
                       envval);

                //for (ii=0;ii<numonoff;ii++) printf("%lld %lld\n", onbins[ii], offbins[ii]);
            }
        }
    }
    chan_avg_temp = gen_dvect(numchan);
    zero_dm_block = gen_fvect(ptsperblk);
    ftmp = gen_fvect(ptsperblk);

    /* Calculate the zero DM time series */
    for (ii = 0; ii < ptsperblk; ii++) {
        zero_dm_block[ii] = 0.0;
        powptr = rawdata + ii * numchan;
        for (jj = 0; jj < numchan; jj++)
            zero_dm_block[ii] += *powptr++;
        ftmp[ii] = zero_dm_block[ii];
    }
    avg_var(ftmp, ptsperblk, &current_avg, &current_std);
    current_std = sqrt(current_std);
    current_med = median(ftmp, ptsperblk);

    /* Calculate the current standard deviation and mean  */
    /* but only for data points that are within a certain */
    /* fraction of the median value.  This removes the    */
    /* really strong RFI from the calculation.            */
    {
        float lo_cutoff, hi_cutoff;
        int numgoodpts = 0;
        lo_cutoff = current_med - 3.0 * current_std;
        hi_cutoff = current_med + 3.0 * current_std;;
        for (jj = 0; jj < numchan; jj++)
            chan_avg_temp[jj] = 0.0;

        /* Find the "good" points */
        for (ii = 0; ii < ptsperblk; ii++) {
            if (zero_dm_block[ii] > lo_cutoff && zero_dm_block[ii] < hi_cutoff) {
                ftmp[numgoodpts] = zero_dm_block[ii];
                powptr = rawdata + ii * numchan;
                for (jj = 0; jj < numchan; jj++)
                    chan_avg_temp[jj] += *powptr++;
                numgoodpts++;
            }
        }

        //printf("avg = %f  med = %f  std = %f  numgoodpts = %d\n", 
        //       current_avg, current_med, current_std, numgoodpts);

        /* Calculate the current average and stddev */
        if (numgoodpts < 1) {
            current_avg = running_avg;
            current_std = running_std;
            for (jj = 0; jj < numchan; jj++)
                chan_avg_temp[jj] = chan_running_avg[jj];
        } else {
            avg_var(ftmp, numgoodpts, &current_avg, &current_std);
            current_std = sqrt(current_std);
            for (jj = 0; jj < numchan; jj++)
                chan_avg_temp[jj] /= numgoodpts;
        }
    }

    /* Update a pseudo running average and stdev */
    if (blocksread) {
        running_avg = 0.9 * running_avg + 0.1 * current_avg;
        running_std = 0.9 * running_std + 0.1 * current_std;
        for (ii = 0; ii < numchan; ii++)
            chan_running_avg[ii] =
                0.9 * chan_running_avg[ii] + 0.1 * chan_avg_temp[ii];
    } else {
        running_avg = current_avg;
        running_std = current_std;
        for (ii = 0; ii < numchan; ii++)
            chan_running_avg[ii] = chan_avg_temp[ii];
        if (current_avg == 0.0)
            printf("Warning:  problem with clipping in first block!!!\n\n");
    }

    /* See if any points need clipping */
    trigger = clip_sigma * running_std;
    for (ii = 0; ii < ptsperblk; ii++) {
        if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
            clipit = 1;
            break;
        }
    }

    /* or alternatively from the CLIPBINSFILE */
    if (numonoff && ((current_point > onbins[onoffindex]
                      && current_point <= offbins[onoffindex])
                     || (current_point + ptsperblk > onbins[onoffindex]
                         && current_point + ptsperblk <= offbins[onoffindex])
                     || (current_point < onbins[onoffindex]
                         && current_point + ptsperblk > offbins[onoffindex])))
        clipit = 1;

    /* Update the good channel levels */
    for (ii = 0; ii < numchan; ii++)
        good_chan_levels[ii] = chan_running_avg[ii];

    /* Replace the bad channel data with channel median values */
    /* that are scaled to equal the running_avg.               */
    if (clipit) {
        for (ii = 0; ii < ptsperblk; ii++) {
            if ((fabs(zero_dm_block[ii] - running_avg) > trigger) ||
                (numonoff && (current_point > onbins[onoffindex]
                              && current_point <= offbins[onoffindex]))) {
                powptr = rawdata + ii * numchan;
                for (jj = 0; jj < numchan; jj++)
                    *powptr++ = good_chan_levels[jj];
                clipped++;

                //fprintf(stderr, "zapping %lld\n", current_point);
            }
            current_point++;
            if (numonoff && current_point > offbins[onoffindex]
                && onoffindex < numonoff - 1) {
                while (current_point > offbins[onoffindex]
                       && onoffindex < numonoff - 1)
                    onoffindex++;

                //printf("updating index to %d\n", onoffindex);
            }
        }
    } else {
        current_point += ptsperblk;
        if (numonoff && current_point > offbins[onoffindex]
            && onoffindex < numonoff - 1) {
            while (current_point > offbins[onoffindex]
                   && onoffindex < numonoff - 1)
                onoffindex++;

            //printf("updating index to %d\n", onoffindex);
        }
    }
    blocksread++;
    vect_free(chan_avg_temp);
    vect_free(zero_dm_block);
    vect_free(ftmp);
    return clipped;
}


/* OLD Clipping Routine (uses channel medians) */
int old_clip_times(float *rawdata, int ptsperblk, int numchan, float clip_sigma,
                   float *good_chan_levels)
// Perform time-domain clipping of rawdata.  This is a 2D array with
// ptsperblk*numchan points, each of which is a float.  The clipping
// is done at clip_sigma sigma above/below the running mean.  The
// up-to-date running averages of the channels are returned in
// good_chan_levels (which must be pre-allocated).
{
    static float *median_chan_levels;
    static float running_avg = 0.0, running_std = 0.0, median_sum = 0.0;
    static int blocksread = 0, firsttime = 1;
    float *zero_dm_block, *median_temp;
    float current_med, trigger, running_wgt = 0.05;
    double current_avg = 0.0, current_std = 0.0, scaling;
    float *powptr;
    int ii, jj, clipit = 0, clipped = 0;
    if (firsttime) {
        median_chan_levels = gen_fvect(numchan);
        firsttime = 0;
    }
    zero_dm_block = gen_fvect(ptsperblk);
    median_temp = gen_fvect(ptsperblk);

    /* Calculate the zero DM time series */
    for (ii = 0; ii < ptsperblk; ii++) {
        zero_dm_block[ii] = 0.0;
        powptr = rawdata + ii * numchan;
        for (jj = 0; jj < numchan; jj++)
            zero_dm_block[ii] += *powptr++;
        median_temp[ii] = zero_dm_block[ii];
    }
    current_med = median(median_temp, ptsperblk);

    /* Calculate the current standard deviation and mean  */
    /* but only for data points that are within a certain */
    /* fraction of the median value.  This removes the    */
    /* really strong RFI from the calculation.            */
    {
        float lo_cutoff, hi_cutoff;
        int numgoodpts = 0;
        lo_cutoff = 0.7 * current_med;
        hi_cutoff = 1.3 * current_med;

        /* Find the "good" points */
        for (ii = 0; ii < ptsperblk; ii++) {
            if (zero_dm_block[ii] > lo_cutoff && zero_dm_block[ii] < hi_cutoff) {
                median_temp[numgoodpts] = zero_dm_block[ii];
                numgoodpts++;
            }
        }

        /* Calculate the current average and stddev */
        if (numgoodpts < 1) {
            current_avg = running_avg;
            current_std = running_std;
        } else {
            avg_var(median_temp, numgoodpts, &current_avg, &current_std);
            current_std = sqrt(current_std);
        }
    }

    /* Update a pseudo running average and stdev */
    if (blocksread) {
        running_avg =
            (running_avg * (1.0 - running_wgt) + running_wgt * current_avg);
        running_std =
            (running_std * (1.0 - running_wgt) + running_wgt * current_std);
    } else {
        running_avg = current_avg;
        running_std = current_std;
        if (running_avg == 0.0 || current_avg == 0.0)
            printf("BAD RFI IN BLOCK#1!!!\n\n");
    }

    /* See if any points need clipping */
    trigger = clip_sigma * running_std;
    for (ii = 0; ii < ptsperblk; ii++) {
        if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
            clipit = 1;
            break;
        }
    }

    /* Calculate the channel medians if required */
    if ((blocksread % BLOCKSTOSKIP == 0 && clipit == 0) || blocksread == 0) {
        median_sum = 0.0;
        for (ii = 0; ii < numchan; ii++) {
            powptr = rawdata + ii;
            for (jj = 0; jj < ptsperblk; jj++)
                median_temp[jj] = *(powptr + jj * numchan);
            median_chan_levels[ii] = median(median_temp, ptsperblk);
            median_sum += median_chan_levels[ii];
        }
    }

    /* Update the good channel levels */
    scaling = running_avg / median_sum;
    for (ii = 0; ii < numchan; ii++)
        good_chan_levels[ii] = median_chan_levels[ii] * scaling;

    /* Replace the bad channel data with channel median values */
    /* that are scaled to equal the running_avg.               */
    if (clipit) {
        for (ii = 0; ii < ptsperblk; ii++) {
            if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
                powptr = rawdata + ii * numchan;
                for (jj = 0; jj < numchan; jj++)
                    *powptr++ = good_chan_levels[jj];
                clipped++;
            }
        }
    }
    blocksread++;
    vect_free(zero_dm_block);
    vect_free(median_temp);
    return clipped;
}
