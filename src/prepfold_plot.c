#include "prepfold.h"
#include "prepfold_cmd.h"
#include "plot2d.h"
#include "float.h"

/*
#define OUTPUT_DMS 1
*/

#define TEST_EQUAL(a, b) (fabs(a) == 0.0 ? \
(fabs((a)-(b)) <= 2 * DBL_EPSILON ? 1 : 0) : \
(fabs((a)-(b))/fabs((a)) <= 2 * DBL_EPSILON ? 1 : 0))

#define TEST_CLOSE(a, b, c) (fabs(a) == 0.0 ? \
(fabs((a)-(b)) <= (c) ? 1 : 0) : \
(fabs((a)-(b))/fabs((a)) <= (c) ? 1 : 0))

/* This is a hack that allows the raw profile data */
/* to be printed on STDOUT so that they can be     */
/* easily imported into Python....                 */
/* #define PRINT_PROFILES 1 */

/********************************************/
/* The following is taken from ppgplot.c by */
/* Nick Patavalis (npat@ariadne.di.uoa.gr)  */
/********************************************/

void minmax(float *v, int nsz, float *min, float *max)
{
    register float *e;
    register float mn, mx;

    for (mn = mx = *v, e = v + nsz; v < e; v++)
        if (*v > mx)
            mx = *v;
        else if (*v < mn)
            mn = *v;
    *min = mn;
    *max = mx;
}

void dminmax(double *v, int nsz, double *min, double *max)
{
    register double *e;
    register double mn, mx;

    for (mn = mx = *v, e = v + nsz; v < e; v++)
        if (*v > mx)
            mx = *v;
        else if (*v < mn)
            mn = *v;
    *min = mn;
    *max = mx;
}


void scaleprof(double *in, float *out, int n, int scalesep, double gmax)
/* Scales an input vector so that it goes from 0.0 to 1.0 */
{
    int ii;
    double norm, min = 1e100, max = -1e100;

    dminmax(in, n, &min, &max);
    if (scalesep) {
        /* Each plot is normalized independently */
        norm = 1.0 / (max - min);
    } else {
        /* All plots are normalized together */
        norm = 1.0 / fabs(gmax);
    }
    if (TEST_CLOSE(min, max, 1.0e-7)) {
        for (ii = 0; ii < n; ii++)
            out[ii] = 0.0;
    } else {
        for (ii = 0; ii < n; ii++)
            out[ii] = (float) ((in[ii] - min) * norm);
    }
}


void fscaleprof(float *in, float *out, int n, int scalesep, double gmax)
/* Scales an input vector so that it goes from 0.0 to 1.0 */
{
    int ii;
    float norm, min = 1e100, max = -1e100;

    minmax(in, n, &min, &max);
    if (scalesep) {
        /* Each plot is normalized independently */
        norm = 1.0 / (max - min);
    } else {
        /* All plots are normalized together */
        norm = 1.0 / gmax;
    }
    if (TEST_CLOSE(min, max, 1.0e-7)) {
        for (ii = 0; ii < n; ii++)
            out[ii] = 0.0;
    } else {
        for (ii = 0; ii < n; ii++)
            out[ii] = (float) ((in[ii] - min) * norm);
    }
}


void lininterp(float min, float max, int npts, float *v)
{
    int ii;
    float step = (max - min) / (npts - 1);

    for (ii = 0; ii < npts; ii++)
        v[ii] = ii * step + min;
}


static void autocal2d(float *a, int rn, int cn,
                      float *fg, float *bg, int nlevels, float *levels,
                      float *x1, float *x2, float *y1, float *y2, float *tr)
{
    //int i;
    float dx1, dx2, dy1, dy2;

    /* autocalibrate intensity-range. */
    if (*fg == *bg) {
        minmax(a, rn * cn, bg, fg);
        //fprintf(stderr,"Intensity range:\n  fg=%f\n  bg=%f\n",*fg,*bg);
    }

    if ((nlevels >= 2) && (levels))
        lininterp(*bg, *fg, nlevels, levels);

    /* autocalibrate x-y range. */
    if ((*x1 == *x2) || (*y1 == *y2))
        cpgqwin(&dx1, &dx2, &dy1, &dy2);
    if (*x1 == *x2) {
        *x1 = dx1;
        *x2 = dx2;
    }
    if (*y1 == *y2) {
        *y1 = dy1;
        *y2 = dy2;
    }
    //fprintf(stderr,"Xrange: [%f, %f]\nYrange[%f, %f]\n",*x1,*x2,*y1,*y2);

    /* calculate transformation vector. */
    tr[2] = tr[4] = 0.0;
    tr[1] = (*x2 - *x1) / cn;
    tr[0] = *x1 - 0.5 * tr[1];
    tr[5] = (*y2 - *y1) / rn;
    tr[3] = *y1 - 0.5 * tr[5];

    //fprintf(stderr,"Tansformation vector:\n");
    //for (i=0; i<6; fprintf(stderr,"  tr[%d]=%f\n",i,tr[i]),i++);
}

/********************************************/

void write_bestprof(prepfoldinfo * search, foldstats * beststats,
                    float *bestprof, double N, double perr,
                    double pderr, double pdderr, double dofeff)
{
    FILE *outfile;
    char *outfilenm;
    int ii;

    outfilenm = (char *) malloc(strlen(search->pgdev) + 10);
    sprintf(outfilenm, "%.*s.bestprof",
            (int) strlen(search->pgdev) - 7, search->pgdev);
    outfile = chkfopen(outfilenm, "w");

    fprintf(outfile, "# Input file       =  %-s\n", search->filenm);
    fprintf(outfile, "# Candidate        =  %-s\n", search->candnm);
    fprintf(outfile, "# Telescope        =  %-s\n", search->telescope);
    if (TEST_EQUAL(search->tepoch, 0.0) || TEST_EQUAL(search->tepoch, -1))
        // -1.0 is for fake data made with makedata
        fprintf(outfile, "# Epoch_topo       =  N/A\n");
    else
        fprintf(outfile, "# Epoch_topo       =  %-.12f\n", search->tepoch);
    if (TEST_EQUAL(search->bepoch, 0.0))
        fprintf(outfile, "# Epoch_bary       =  N/A\n");
    else
        fprintf(outfile, "# Epoch_bary (MJD) =  %-.12f\n", search->bepoch);
    fprintf(outfile, "# T_sample         =  %.6g\n", search->dt);
    fprintf(outfile, "# Data Folded      =  %-.0f\n", N);
    fprintf(outfile, "# Data Avg         =  %-17.15g\n", beststats->data_avg);
    fprintf(outfile, "# Data StdDev      =  %-17.15g\n", sqrt(beststats->data_var));
    fprintf(outfile, "# Profile Bins     =  %d\n", search->proflen);
    fprintf(outfile, "# Profile Avg      =  %-17.15g\n", beststats->prof_avg);
    fprintf(outfile, "# Profile StdDev   =  %-17.15g\n", sqrt(beststats->prof_var));

    {
        int chidf;
        double chip, chi_lnp, chi_sig;
        char out2[80];

        chidf = (int)(dofeff + 1e-12);
        chi_lnp = chi2_logp(beststats->redchi * chidf, chidf);
        chip = (chi_lnp < -700) ? 0.0 : exp(chi_lnp);
        chi_sig = equivalent_gaussian_sigma(chi_lnp);
        sprintf(out2, "(~%.1f sigma)", chi_sig);
        fprintf(outfile, "# Reduced chi-sqr  =  %.3f\n", beststats->redchi);
        fprintf(outfile, "# Prob(Noise)      <  %.3g   %s\n", chip, out2);
        if (search->nsub > 1)
            fprintf(outfile, "# Best DM          =  %.3f\n", search->bestdm);
        {
            if (search->tepoch != 0.0) {
                fprintf(outfile, "# P_topo (ms)      =  %-17.15g +/- %-.3g\n",
                        search->topo.p1 * 1000.0, perr * 1000.0);
                fprintf(outfile, "# P'_topo (s/s)    =  %-17.15g +/- %-.3g\n",
                        search->topo.p2, pderr);
                fprintf(outfile, "# P''_topo (s/s^2) =  %-17.15g +/- %-.3g\n",
                        search->topo.p3, pdderr);
            } else {
                fprintf(outfile, "# P_topo (ms)      =  N/A\n");
                fprintf(outfile, "# P'_topo (s/s)    =  N/A\n");
                fprintf(outfile, "# P''_topo (s/s^2) =  N/A\n");
            }
            if (search->bepoch != 0.0) {
                fprintf(outfile, "# P_bary (ms)      =  %-17.15g +/- %-.3g\n",
                        search->bary.p1 * 1000.0, perr * 1000.0);
                fprintf(outfile, "# P'_bary (s/s)    =  %-17.15g +/- %-.3g\n",
                        search->bary.p2, pderr);
                fprintf(outfile, "# P''_bary (s/s^2) =  %-17.15g +/- %-.3g\n",
                        search->bary.p3, pdderr);
            } else {
                fprintf(outfile, "# P_bary (ms)      =  N/A\n");
                fprintf(outfile, "# P'_bary (s/s)    =  N/A\n");
                fprintf(outfile, "# P''_bary (s/s^2) =  N/A\n");
            }
        }
        if (TEST_EQUAL(search->orb.p, 0.0)) {
            fprintf(outfile, "# P_orb (s)        =  N/A\n");
            fprintf(outfile, "# asin(i)/c (s)    =  N/A\n");
            fprintf(outfile, "# eccentricity     =  N/A\n");
            fprintf(outfile, "# w (rad)          =  N/A\n");
            fprintf(outfile, "# T_peri           =  N/A\n");
        } else {
            fprintf(outfile, "# P_orb (s)        =  %-17.15g\n", search->orb.p);
            fprintf(outfile, "# asin(i)/c (s)    =  %-17.15g\n", search->orb.x);
            fprintf(outfile, "# eccentricity     =  %-17.15g\n", search->orb.e);
            fprintf(outfile, "# w (deg)          =  %-17.15g\n", search->orb.w);
            fprintf(outfile, "# T_peri           =  %-.12f\n", search->orb.t);
        }
    }
    fprintf(outfile, "######################################################\n");
    for (ii = 0; ii < search->proflen; ii++)
        fprintf(outfile, "%4d  %.7g\n", ii, bestprof[ii]);
    fclose(outfile);
    free(outfilenm);
}


void CSS_profs(double *inprofs, double *outprofs,
               foldstats * instats, int numprofs, int proflen,
               double *delays, double *sumprof, foldstats * sumstats,
               float *timechi, float chifact)
/* Combine, Scale and Shift 'numprofs' profiles, of length 'proflen',   */
/* into a single profile of length 'proflen'.  The profiles are         */
/* summed after the appropriate 'delays' are added to each profile.     */
{
    int ii, jj, index = proflen, offset;
    double rdof, redchi = 0.0, *local_delays;

    /* Initiate the output statistics */
    initialize_foldstats(sumstats);
    sumstats->numprof = proflen;
    rdof = 1.0 / (proflen - 1.0);
    local_delays = gen_dvect(numprofs);

    /* Convert all the delays to positive offsets from   */
    /* the phase=0 profile bin, in units of profile bins */
    for (ii = 0; ii < numprofs; ii++) {
        local_delays[ii] = fmod(delays[ii], proflen);
        if (local_delays[ii] < 0.0)
            local_delays[ii] += proflen;
    }

    /* Set the output arrays to all zeros */
    for (ii = 0; ii < proflen; ii++)
        sumprof[ii] = 0.0;

    /* Start defining the time vs chi-squared array */
    timechi[0] = 1.0;

    /* Loop over the profiles */
    for (ii = 0; ii < numprofs; ii++) {

        /* Calculate the appropriate offset into the profile array */
        offset = (int) (local_delays[ii] + 0.5);

        /* Shift and copy the profiles */
        index = ii * proflen;
        shift_prof(inprofs + index, proflen, offset, outprofs + index);

        /* Now add the current profile to the summed profile */
        for (jj = 0; jj < proflen; jj++, index++)
            sumprof[jj] += outprofs[index];

        /* Update the output statistics structure */
        sumstats->numdata += instats[ii].numdata;
        sumstats->data_avg += instats[ii].data_avg;
        sumstats->data_var += instats[ii].data_var;
        sumstats->prof_avg += instats[ii].prof_avg;
        sumstats->prof_var += instats[ii].prof_var;

        /* Calculate the current chi-squared */
        redchi = chisqr(sumprof, proflen, sumstats->prof_avg,
                        sumstats->prof_var) * rdof * chifact;
        timechi[ii + 1] = (float) redchi;
    }

    /* Profile information gets added together, but */
    /* data set info gets averaged together.        */
    sumstats->data_avg /= numprofs;
    sumstats->data_var /= numprofs;
    sumstats->redchi = redchi;
    vect_free(local_delays);
}


void prepfold_plot(prepfoldinfo * search, plotflags * flags, int xwin, float *ppdot)
/* Make the beautiful 1 page prepfold output */
{
    int ii, jj, profindex = 0, loops = 1, ct, bestidm = 0, bestip = 0, bestipd = 0;
    double N = 0.0, T, dofeff;
    double parttime, bestp, bestpd, bestpdd;
    double perr, pderr, pdderr;
    double pfold, pdfold, pddfold = 0.0;
    float ftmp, *ftmparr1, chifact, dt_per_bin;
    foldstats currentstats, beststats;
    /* Best Fold Plot */
    double *dbestprof = NULL;
    float *bestprof = NULL;
    /* Profiles vs Time */
    float *timeprofs = NULL, *parttimes = NULL;
    /* RedChi vs Time */
    float *timechi = NULL;
    /* Profiles vs DM */
    float *dmprofs = NULL;
    /* DM vs RedChi */
    float *dmchi = NULL;
    /* Period vs RedChi */
    float *periodchi = NULL;
    /* P-dot vs RedChi */
    float *pdotchi = NULL;
    /* Period P-dot 2D */
    float *ppdot2d = ppdot;

    if (xwin)
        loops = 2;

    if (flags->showfold) {
        switch_f_and_p(search->fold.p1, search->fold.p2, search->fold.p3,
                       &bestp, &bestpd, &bestpdd);
    } else if (TEST_EQUAL(search->fold.pow, 1.0)) {     /* Barycentric periods */
        bestp = search->bary.p1;
        bestpd = search->bary.p2;
        bestpdd = search->bary.p3;
    } else {                    /* Topocentric periods */
        bestp = search->topo.p1;
        bestpd = search->topo.p2;
        bestpdd = search->topo.p3;
    }
    // The number of samples that can fit across one profile bin
    dt_per_bin = bestp / search->proflen / search->dt;

    // Determine the chi^2 correction factor due to the correlations
    // in the profile bins caused by fold()
    dofeff = (search->proflen - 1.0) * DOF_corr(dt_per_bin);
    chifact = 1.0 / DOF_corr(dt_per_bin);
    if (flags->events || flags->samples) {
        chifact = 1.0;
        dofeff = search->proflen - 1.0;
    }
    // Get off-pulse reduced-chi^2
    ftmp = estimate_offpulse_redchi2(search->rawfolds, search->stats,
                                     search->npart, search->nsub,
                                     search->proflen, 50, dofeff);
    printf("Effective number of DOF = %.2f\n", dofeff);
    printf("Off-pulse Reduced chi^2 correction factor = %.2f\n", ftmp);
    if (flags->fixchi) {
        chifact *= ftmp;
    }

    /* Find the indices for the best periods, p-dots, and DMs */
    for (ii = 0; ii < search->numperiods; ii++) {
        if (TEST_EQUAL(search->periods[ii], bestp))
            bestip = ii;
        if (TEST_EQUAL(search->pdots[ii], bestpd))
            bestipd = ii;
    }
    if (search->nsub > 1) {
        for (ii = 0; ii < search->numdms; ii++) {
            if (TEST_EQUAL(search->dms[ii], search->bestdm))
                bestidm = ii;
        }
    } else {
        bestidm = 0;
    }

/* This zaps subbands and or parts in the folds */
#if (0)
    int part, sub, phase;

    for (part = 0; part < search->npart; part++) {
        for (sub = 0; sub < search->nsub; sub++) {
            for (phase = 0; phase < search->proflen; phase++) {
                if (part == 30 || sub > 31)
                    search->rawfolds[part * (search->proflen * search->nsub) +
                                     sub * search->proflen + phase] = 0.0;
            }
        }
    }
#endif

    /* Attempt to make show_pfd work on polyco-folded files */
    if (search->fold.pow == 0.0 && ppdot == NULL &&
        bestip == (search->numperiods - 1) / 2 &&
        bestipd == (search->numperiods - 1) / 2 && !flags->showfold) {
        printf("Assuming this was folded with polycos...\n");
        search->fold.p1 = 1.0 / search->topo.p1;
    }

    switch_f_and_p(search->fold.p1, search->fold.p2, search->fold.p3,
                   &pfold, &pdfold, &pddfold);

    /* Find out how many total points were folded */
    for (ii = 0; ii < search->npart; ii++)
        N += search->stats[ii * search->nsub].numdata;

    /* Calculate the time per part and the total observation time */
    T = N * search->dt;
    parttime = search->stats[0].numdata * search->dt;
    parttimes = gen_freqs(search->npart + 1, 0.0, parttime);

    {                           /* Generate the data we need for the plots */
        double df, dfd, dfdd = 0.0;
        double *currentprof, *ddprofs = search->rawfolds;
        double *delays, *pd_delays, *pdd_delays;
        foldstats *ddstats = search->stats;

        delays = gen_dvect(search->npart);
        pd_delays = gen_dvect(search->npart);
        pdd_delays = gen_dvect(search->npart);
        currentprof = gen_dvect(search->proflen);
        if (search->nsub > 1) {
            ddprofs = gen_dvect(search->npart * search->proflen);
            ddstats = (foldstats *) malloc(search->npart * sizeof(foldstats));
        }

        /* Calculate the delays for the pdotdot */
        if (bestpdd != 0.0)     /* bestpdd=0.0 only if there was no searching over pdd */
            dfdd = switch_pfdotdot(pfold, pdfold, bestpdd) - search->fold.p3;
        for (ii = 0; ii < search->npart; ii++)
            pdd_delays[ii] = fdotdot2phasedelay(dfdd, parttimes[ii]);

        {                       /* Correct for and fold the best profile */
            double *tmp_profs, gmin = 1e100, gmax = -1e100;

            bestprof = gen_fvect(2 * search->proflen);
            dbestprof = gen_dvect(search->proflen);

            /* De-disperse if needed */
            if (search->nsub > 1)
                correct_subbands_for_DM(search->dms[bestidm], search, ddprofs,
                                        ddstats);

            /* Determine the global min and max in the profiles */
            dminmax(ddprofs, search->npart * search->proflen, &gmin, &gmax);
            // Take care of a possible rare corner case
            if (gmax == 0.0)
                gmax = 1.0;

            /* Compute the errors in fdot, and f to correct */
            df = 1.0 / bestp - search->fold.p1;
            dfd = switch_pfdot(pfold, bestpd) - search->fold.p2;

            /* Compute the delays for the best profile */
            for (ii = 0; ii < search->npart; ii++)
                delays[ii] = (pdd_delays[ii] + fdot2phasedelay(dfd, parttimes[ii]) +
                              df * parttimes[ii]) * search->proflen;

            /* Create the best profile */
            timechi = gen_fvect(search->npart + 1);
            tmp_profs = gen_dvect(search->npart * search->proflen);
            CSS_profs(ddprofs, tmp_profs, ddstats, search->npart,
                      search->proflen, delays, dbestprof, &beststats, timechi,
                      chifact);
            double2float(dbestprof, bestprof, search->proflen);
            double2float(dbestprof, bestprof + search->proflen, search->proflen);

            /* Create the time vs phase plot data */
            timeprofs = gen_fvect(2 * search->proflen * search->npart);
            for (ii = 0; ii < search->npart; ii++) {
                profindex = ii * search->proflen;
                scaleprof(tmp_profs + profindex, timeprofs + 2 * profindex,
                          search->proflen, flags->scaleparts, gmax);
                memcpy(timeprofs + 2 * profindex + search->proflen,
                       timeprofs + 2 * profindex, search->proflen * sizeof(float));
            }
            vect_free(tmp_profs);
        }

        if (ppdot == NULL) {    /* Generate the p-pdot plane */
            int ip, ipd;

            ppdot2d = gen_fvect(search->numperiods * search->numpdots);
            for (ipd = 0; ipd < search->numpdots; ipd++) {      /* Loop over the pds */

                /* Compute the error in fdot, and its delays */
                dfd = switch_pfdot(pfold, search->pdots[ipd]) - search->fold.p2;
                for (ii = 0; ii < search->npart; ii++)
                    pd_delays[ii] =
                        pdd_delays[ii] + fdot2phasedelay(dfd, parttimes[ii]);

                for (ip = 0; ip < search->numperiods; ip++) {   /* Loop over the ps */
                    df = 1.0 / search->periods[ip] - search->fold.p1;

                    /* Compute the phase offsets for each subintegration */
                    for (ii = 0; ii < search->npart; ii++)
                        delays[ii] = (pd_delays[ii] +
                                      df * parttimes[ii]) * search->proflen;

                    /* Combine the profiles usingthe above computed delays */
                    combine_profs(ddprofs, ddstats, search->npart, search->proflen,
                                  delays, currentprof, &currentstats);
                    ppdot2d[ipd * search->numperiods + ip] =
                        currentstats.redchi * chifact;
                }
            }
        }

        {                       /* Create the p vs chi and pd vs chi plots */

            periodchi = gen_fvect(search->numperiods);
            for (ii = 0; ii < search->numperiods; ii++)
                periodchi[ii] = ppdot2d[bestipd * search->numperiods + ii];
            pdotchi = gen_fvect(search->numpdots);
            for (ii = 0; ii < search->numperiods; ii++)
                pdotchi[ii] = ppdot2d[ii * search->numperiods + bestip];
        }

        if (search->nsub > 1) { /* For data with subbands */

            /* Generate the DM vs chi plot */
            dmchi = gen_fvect(search->numdms);

            /* Compute the errors in fdotdot, fdot, and f to correct */
            df = 1.0 / bestp - search->fold.p1;
            dfd = switch_pfdot(pfold, bestpd) - search->fold.p2;

            /* Compute the delays for the best profile */
            for (ii = 0; ii < search->npart; ii++)
                delays[ii] = (pdd_delays[ii] + fdot2phasedelay(dfd, parttimes[ii]) +
                              df * parttimes[ii]) * search->proflen;

            /* De-disperse and fold */
            for (ii = 0; ii < search->numdms; ii++) {
                correct_subbands_for_DM(search->dms[ii], search, ddprofs, ddstats);
                combine_profs(ddprofs, ddstats, search->npart,
                              search->proflen, delays, currentprof, &currentstats);
                dmchi[ii] = currentstats.redchi * chifact;
            }

            {                   /* Generate the Subband vs phase plot */
                double *dmdelays, hif, dopplerhif, hifdelay, rdphase;
                double *tmpprofs, *totdelays, *subbanddelays;

                dmprofs = gen_fvect(search->nsub * search->proflen);
                tmpprofs = gen_dvect(search->npart * search->proflen);
                totdelays = gen_dvect(search->npart);

                /* Compute the DM-delays for each subband */
                rdphase = search->fold.p1 * search->proflen;
                hif = search->lofreq + (search->numchan - 1.0) * search->chan_wid;
                dopplerhif = doppler(hif, search->avgvoverc);
                hifdelay = delay_from_dm(search->dms[bestidm], dopplerhif);
                subbanddelays =
                    subband_delays(search->numchan, search->nsub,
                                   search->dms[bestidm], search->lofreq,
                                   search->chan_wid, search->avgvoverc);
                dmdelays = gen_dvect(search->nsub);
                for (ii = 0; ii < search->nsub; ii++)
                    dmdelays[ii] = -(subbanddelays[ii] - hifdelay) * rdphase;
                vect_free(subbanddelays);

                /* Fold each subband */
                for (ii = 0; ii < search->nsub; ii++) {

                    /* Create a temporary array filled with just the current subbands profiles */
                    for (jj = 0; jj < search->npart; jj++) {
                        memcpy(tmpprofs + jj * search->proflen,
                               search->rawfolds +
                               search->proflen * (jj * search->nsub + ii),
                               sizeof(double) * search->proflen);
                        totdelays[jj] = delays[jj] + dmdelays[ii];
                    }

                    /* Create the profile */
                    /* Note that the stats will be incorrect in the following... */
                    combine_profs(tmpprofs, ddstats, search->npart,
                                  search->proflen, totdelays, currentprof,
                                  &currentstats);
                    double2float(currentprof, dmprofs + ii * search->proflen,
                                 search->proflen);
                }
                vect_free(totdelays);
                vect_free(tmpprofs);
                vect_free(dmdelays);
            }
        }

        vect_free(delays);
        vect_free(pd_delays);
        vect_free(pdd_delays);
        vect_free(currentprof);
        if (search->nsub > 1) {
            vect_free(ddprofs);
            free(ddstats);
        }
    }

    /* Calculate the errors in the pulsation quantities */

    if (search->tepoch != 0.0 || (search->tepoch == 0.0 && search->bary.p1 == 0.0))
        fold_errors(dbestprof, search->proflen, search->dt, N,
                    beststats.data_var, search->topo.p1, search->topo.p2,
                    search->topo.p3, &perr, &pderr, &pdderr);
    else
        fold_errors(dbestprof, search->proflen, search->dt, N,
                    beststats.data_var, search->bary.p1, search->bary.p2,
                    search->bary.p3, &perr, &pderr, &pdderr);
    vect_free(dbestprof);

    write_bestprof(search, &beststats, bestprof, N, perr, pderr, pdderr, dofeff);

    /*
     *  Now plot the results
     */

    for (ct = 0; ct < loops; ct++) {
        float min, max, over;

        /*Set the PGPLOT device to an X-Window */

        if (ct == 1)
            strcpy(search->pgdev, "/XWIN");

        /* Open and prep our device */

        cpgopen(search->pgdev);
        if (!flags->justprofs) {
            cpgpap(10.25, 8.5 / 11.0);
            cpgpage();
            cpgiden();
        }                       /* else {
                                   cpgpap(10.0, 1.0);
                                   } */
        cpgslw(2);
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
            for (ii = 0; ii < search->npart; ii++) {
                if (ii == search->npart - 1)
                    printf("%.2f])\n\n", search->stats[ii].prof_avg);
                else
                    printf("%.2f, ", search->stats[ii].prof_avg);
            }
            printf("\n\nstdevs = Numeric.asarray([");
            for (ii = 0; ii < search->npart; ii++) {
                if (ii == search->npart - 1)
                    printf("%.2f])\n\n", sqrt(search->stats[ii].prof_var));
                else
                    printf("%.2f, ", sqrt(search->stats[ii].prof_var));
            }
            printf("\n\nmaxvals = Numeric.asarray([");
            for (ii = 0; ii < search->npart; ii++) {
                minmax(timeprofs + ii * 2 * search->proflen,
                       search->proflen, &min, &max);
                if (ii == search->npart - 1)
                    printf("%.2f])\n\n", max);
                else
                    printf("%.2f, ", max);
            }
        }
#endif

        cpgsvp(0.06, 0.27, 0.09, 0.68);
        cpgswin(0.0, 1.999, 0.0, T);
        {
            int mincol, maxcol, numcol, nr, nc;
            float l[2] = { 0.0, 1.0 };
            float r[2] = { 1.0, 0.0 };
            float g[2] = { 1.0, 0.0 };
            float b[2] = { 1.0, 0.0 };
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
            cpgimag(timeprofs, nc, nr, 0 + 1, nc, 0 + 1, nr, bg, fg, tr);
            vect_free(levels);
        }
        if (1) {                /* if 0 skip the chi-squared vs time plot */
            cpgbox("BCNST", 0.0, 0, "BNST", 0.0, 0);
            // Rescale window and provide ticks for each subint
            cpgswin(0.0, 1.999, 0.0, search->npart);
            cpgbox("", 0.0, 0, "CTSI", 5.0, 5);
        } else
            cpgbox("BCNST", 0.0, 0, "BNST", 0.0, 0);
        cpgmtxt("B", 2.6, 0.5, 0.5, "Phase");
        cpgmtxt("L", 2.1, 0.5, 0.5, "Time (s)");

        /*  Time versus Reduced chisqr */

        find_min_max_arr(search->npart + 1, timechi, &min, &max);
        if (!flags->justprofs) {        /* if 0 skip the chi-squared vs time plot */
            if (search->nsub > 1)
                cpgsvp(0.27, 0.36, 0.09, 0.68);
            else
                cpgsvp(0.27, 0.39, 0.09, 0.68);
            cpgswin(1.1 * max, 0.0, 0.0, T);
            //         cpgbox("BCNST", 0.0, 0, "BST", 0.0, 0);
            cpgbox("BCNST", 0.0, 0, "B", 0.0, 0);
            cpgmtxt("B", 2.6, 0.5, 0.5, "Reduced \\gx\\u2\\d");
            cpgline(search->npart + 1, timechi, parttimes);
            cpgswin(1.1 * max, 0.0, search->startT - 0.0001, search->endT + 0.0001);
            if (search->nsub > 1)
                cpgsch(0.7);
            cpgbox("", 0.0, 0, "CMST", 0.0, 0);
            cpgmtxt("R", 2.3, 0.5, 0.5, "Fraction of Observation");
            if (search->nsub > 1)
                cpgsch(0.8);
        }

        /* Combined best profile */

        {
            float x[2] = { -0.2, 2.0 }, avg[2];
            float errx = -0.1, erry = beststats.prof_avg, errlen;
            float *phasetwo = NULL;

            cpgsvp(0.06, 0.27, 0.68, 0.94);
            cpgswin(0.0, 1.999, 0.0, 1.0);
            cpgbox("BST", 0.0, 0, "", 0.0, 0);
            cpgsvp(0.039, 0.27, 0.68, 0.94);
            find_min_max_arr(2 * search->proflen, bestprof, &min, &max);
            over = 0.1 * (max - min);
            cpgswin(-0.2, 2.0, min - 2.0 * over, max + over);
            if (!flags->justprofs)
                cpgmtxt("T", 0.0, 0.5, 0.5, "2 Pulses of Best Profile");
            phasetwo = gen_freqs(2 * search->proflen, 0.0, 1.0 / search->proflen);
            cpgline(2 * search->proflen, phasetwo, bestprof);
            vect_free(phasetwo);
            if (!flags->justprofs) {
                cpgsls(4);
                avg[0] = avg[1] = beststats.prof_avg;
                cpgline(2, x, avg);
                cpgsls(1);
                errlen = sqrt(beststats.prof_var);
                cpgerr1(6, errx, erry, errlen, 2);
                cpgpt(1, &errx, &erry, 5);
                // Only do the following for radio data
                if (search->lofreq > 0.0 && search->chan_wid > 0.0) {
                    float fmid, tdms, tdt, tcbw, ttot = 0.0;
                    // Middle observing freq
                    fmid =
                        search->lofreq + 0.5 * (search->numchan -
                                                1.0) * search->chan_wid;
                    // Dispersion smearing from a single channel
                    tdms = smearing_from_bw(search->bestdm, fmid, search->chan_wid);
                    // Sample time
                    tdt = search->dt;
                    // 1/channel BW (i.e. best possible time res)
                    tcbw = 1.0 / (search->chan_wid * 1e6);
                    // total time resolution
                    ttot = sqrt(tdms * tdms + tdt * tdt + tcbw * tcbw);
                    // The above in units of a fraction of a period
                    if (search->bepoch != 0.0)
                        ttot /= search->bary.p1;
                    else
                        ttot /= search->topo.p1;
                    // Following is 0.5 * ttot since cpgerr1 gives +/- err in plot
                    cpgerr1(5, 1.0, min - over, 0.5 * ttot, 2);
                }
            }
        }

        if (!flags->justprofs) {

            if (search->nsub > 1) {

                /* DM vs reduced chisqr */

                cpgsvp(0.44, 0.66, 0.09, 0.22);
                find_min_max_arr(search->numdms, dmchi, &min, &max);
#ifdef OUTPUT_DMS
                {
                    int dmnum;

                    printf
                        ("\nThe raw data from the DM vs. Reduced chi-sqr plot:\n\n");
                    printf("   DM      Reduced chi-sqr\n");
                    for (dmnum = 0; dmnum < search->numdms; dmnum++)
                        printf("  %.3f     %.3f\n", search->dms[dmnum],
                               dmchi[dmnum]);
                    printf("\n");
                }
#endif
                cpgswin(search->dms[0], search->dms[search->numdms - 1], 0.0,
                        1.1 * max);
                cpgsch(0.7);
                cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
                cpgmtxt("L", 2.0, 0.5, 0.5, "Reduced \\gx\\u2\\d");
                cpgsch(0.8);
                cpgmtxt("B", 2.6, 0.5, 0.5, "DM (pc/cm\\u3\\d)");
                ftmparr1 = gen_fvect(search->numdms);
                double2float(search->dms, ftmparr1, search->numdms);
                cpgline(search->numdms, ftmparr1, dmchi);
                vect_free(ftmparr1);

                /* Plots for each subband */

                if (0) {
                    int chanpersb;
                    double lofreq, hifreq, losubfreq, hisubfreq;
                    float *tmpprof, dsubf, foffset, fnorm, *phaseone;

                    tmpprof = gen_fvect(search->proflen + 1);
                    chanpersb = search->numchan / search->nsub;
                    dsubf = chanpersb * search->chan_wid;
                    lofreq = search->lofreq + dsubf - search->chan_wid;
                    hifreq =
                        search->lofreq + search->nsub * dsubf - search->chan_wid;
                    losubfreq = doppler(lofreq, search->avgvoverc);
                    hisubfreq = doppler(hifreq, search->avgvoverc);
                    cpgsvp(0.44, 0.66, 0.3, 0.68);
                    cpgswin(0.0 - 0.01, 1.0 + 0.01, 0.001, search->nsub + 0.999);
                    cpgsch(0.7);
                    cpgbox("BCNST", 0.2, 2, "BNST", 0.0, 0);
                    cpgmtxt("L", 2.0, 0.5, 0.5, "Sub-band");
                    cpgswin(0.0 - 0.01, 1.0 + 0.01, losubfreq - dsubf,
                            hisubfreq + dsubf);
                    cpgbox("", 0.2, 2, "CMST", 0.0, 0);
                    cpgmtxt("R", 2.3, 0.5, 0.5, "Frequency (MHz)");
                    cpgsch(0.8);
                    cpgmtxt("B", 2.5, 0.5, 0.5, "Phase");
                    phaseone =
                        gen_freqs(search->proflen + 1, 0.0, 1.0 / search->proflen);
                    for (ii = 0; ii < search->nsub; ii++) {
                        find_min_max_arr(search->proflen,
                                         dmprofs + ii * search->proflen, &min, &max);
                        foffset =
                            doppler(lofreq + (ii - 0.45) * dsubf, search->avgvoverc);
                        if (min == max) {
                            for (jj = 0; jj < search->proflen; jj++)
                                tmpprof[jj] = 0.45 * dsubf + foffset;
                        } else {
                            fnorm = 0.9 * dsubf / (max - min);
                            for (jj = 0; jj < search->proflen; jj++)
                                tmpprof[jj] = (dmprofs[ii * search->proflen + jj] -
                                               min) * fnorm + foffset;
                        }
                        tmpprof[search->proflen] = tmpprof[0];
                        cpgline(search->proflen + 1, phaseone, tmpprof);
                    }
                    vect_free(phaseone);
                    vect_free(tmpprof);
                } else {
                    double lofreq, hifreq, losubfreq, hisubfreq;

                    lofreq = search->lofreq - 0.5 * search->chan_wid;
                    hifreq = lofreq + search->numchan * search->chan_wid;
                    losubfreq = doppler(lofreq, search->avgvoverc);
                    hisubfreq = doppler(hifreq, search->avgvoverc);
                    cpgsvp(0.44, 0.66, 0.3, 0.68);
                    cpgswin(0.0, 2.0, 0.0, search->nsub);
                    {
                        int mincol, maxcol, numcol, nr, nc;
                        float min, max, gmax = -1.0e100;
                        float l[2] = { 0.0, 1.0 };
                        float r[2] = { 1.0, 0.0 };
                        float g[2] = { 1.0, 0.0 };
                        float b[2] = { 1.0, 0.0 };
                        float fg = 0.0, bg = 0.0, tr[6], *levels;
                        float x1 = 0.0, y1 = 0.0, x2 = 1.0, y2 = search->nsub;
                        float *subprofs;

                        nr = search->nsub;
                        nc = search->proflen;
                        subprofs = gen_fvect(nr * nc);
                        cpgqcol(&mincol, &maxcol);
                        mincol += 2;
                        cpgscir(mincol, maxcol);
                        numcol = maxcol - mincol + 1;
                        levels = gen_fvect(numcol);
                        cpgctab(l, r, g, b, numcol, 1.0, 0.5);
                        for (ii = 0; ii < search->nsub; ii++) {
                            profindex = ii * search->proflen;
                            find_min_max_arr(search->proflen,
                                             dmprofs + profindex, &min, &max);
                            if (max > gmax)
                                gmax = max;
                        }
                        for (ii = 0; ii < search->nsub; ii++) {
                            profindex = ii * search->proflen;
                            fscaleprof(dmprofs + profindex, subprofs + profindex,
                                       search->proflen, flags->scaleparts, gmax);
                        }
                        autocal2d(subprofs, nr, nc, &fg, &bg, numcol,
                                  levels, &x1, &x2, &y1, &y2, tr);
                        cpgimag(subprofs, nc, nr, 0 + 1, nc, 0 + 1, nr, bg, fg, tr);
                        tr[0] += 1.0;
                        cpgimag(subprofs, nc, nr, 0 + 1, nc, 0 + 1, nr, bg, fg, tr);
                        vect_free(subprofs);
                        vect_free(levels);
                    }
                    cpgswin(0.0 - 0.01, 2.0 + 0.01, 0.0, search->nsub);
                    cpgsch(0.7);
                    cpgbox("CST", 0.4, 2, "", 0.0, 0);
                    cpgbox("BNSTI", 0.4, 2, "BNSTI", 0.0, 0);
                    cpgmtxt("L", 2.0, 0.5, 0.5, "Sub-band");
                    cpgswin(0.0 - 0.01, 2.0 + 0.01, losubfreq, hisubfreq);
                    cpgbox("", 0.4, 2, "CMSTI", 0.0, 0);
                    cpgmtxt("R", 2.5, 0.5, 0.5, "Frequency (MHz)");
                    cpgsch(0.8);
                    cpgmtxt("B", 2.5, 0.5, 0.5, "Phase");
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
                float l[10] = { 0.0, 0.035, 0.045, 0.225,
                    0.4, 0.41, 0.6, 0.775, 0.985, 1.0
                };
                float r[10] = { 1.0, 1.0, 0.947, 0.0, 0.0,
                    0.0, 0.0, 1.0, 1.0, 1.0
                };
                float g[10] = { 1.0, 0.844, 0.8, 0.0, 0.946,
                    1.0, 1.0, 1.0, 0.0, 0.0
                };
                float b[10] = { 1.0, 1.0, 1.0, 1.0, 1.0, 0.95,
                    0.0, 0.0, 0.0, 0.0
                };
                float fg = 0.0, bg = 0.0, tr[6], *levels, errlen;
                float x1l, x1h, y1l, y1h, x2l, x2h, y2l, y2h;
                char pout[100], pdout[100], fout[100], fdout[100];

                if (flags->allgrey) {
                    /* ApJ Grey  White to dark grey... */
                    l[0] = 0.0;
                    l[1] = 1.0;
                    r[0] = 1.0;
                    r[1] = 0.25;
                    g[0] = 1.0;
                    g[1] = 0.25;
                    b[0] = 1.0;
                    b[1] = 0.25;
                }

                /* Plot Boundaries */

                /* Period / P-dot */
                x1l = (search->periods[0] - pfold) * 1000.0;
                x1h = (search->periods[search->numperiods - 1] - pfold) * 1000.0;
                y1l = search->pdots[0] - pdfold;
                y1h = search->pdots[search->numpdots - 1] - pdfold;
                /* Frequency / F-dot */
                x2l = 1.0 / search->periods[0] - search->fold.p1;
                x2h =
                    1.0 / search->periods[search->numperiods - 1] - search->fold.p1;
                y2l = switch_pfdot(pfold, search->pdots[0]) - search->fold.p2;
                y2h = switch_pfdot(pfold, search->pdots[search->numperiods - 1]) -
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
                    sprintf(fdout, "F-dot + %-.5g (Hz/s)", fabs(search->fold.p2));
                else if (TEST_EQUAL(search->fold.p2, 0.0))
                    sprintf(fdout, "F-dot (Hz/s)");
                else
                    sprintf(fdout, "F-dot - %-.5g (Hz/s)", search->fold.p2);

                /* Period vs reduced chisqr */

                cpgsch(0.8);
                ftmparr1 = gen_fvect(search->numperiods);
                for (ii = 0; ii < search->numperiods; ii++)
                    ftmparr1[ii] = (search->periods[ii] - pfold) * 1000.0;
                find_min_max_arr(search->numperiods, periodchi, &min, &max);
                if (search->nsub > 1) {
                    cpgsvp(0.74, 0.94, 0.41, 0.51);
                    cpgswin(x1l, x1h, 0.0, 1.1 * max);
                    cpgline(search->numperiods, ftmparr1, periodchi);
                    cpgsch(0.5);
                    cpgbox("BCNST", 0.0, 0, "BCMST", 0.0, 0);
                    cpgsch(0.7);
                    cpgmtxt("B", 2.2, 0.5, 0.5, pout);
                    cpgmtxt("R", 2.4, 0.5, 0.5, "Reduced \\gx\\u2\\d");
                } else {
                    cpgsvp(0.51, 0.82, 0.49, 0.63);
                    cpgswin(x1l, x1h, 0.001, 1.1 * max);
                    cpgline(search->numperiods, ftmparr1, periodchi);
                    cpgsch(0.7);
                    cpgbox("BST", 0.0, 0, "BCMST", 0.0, 0);
                    cpgswin(x2l, x2h, 0.001, 1.1 * max);
                    cpgbox("CMST", 0.0, 0, "", 0.0, 0);
                    cpgsch(0.8);
                    cpgmtxt("T", 1.8, 0.5, 0.5, fout);
                }
                vect_free(ftmparr1);

                /* P-dot vs reduced chisqr */

                ftmparr1 = gen_fvect(search->numpdots);
                for (ii = 0; ii < search->numpdots; ii++)
                    ftmparr1[ii] = search->pdots[ii] - pdfold;
                find_min_max_arr(search->numpdots, pdotchi, &min, &max);
                if (search->nsub > 1) {
                    cpgsvp(0.74, 0.94, 0.58, 0.68);
                    cpgswin(y1l, y1h, 0.0, 1.1 * max);
                    cpgline(search->numpdots, ftmparr1, pdotchi);
                    cpgsch(0.5);
                    cpgbox("BCNST", 0.0, 0, "BCMST", 0.0, 0);
                    cpgsch(0.7);
                    cpgmtxt("B", 2.2, 0.5, 0.5, pdout);
                    cpgmtxt("R", 2.4, 0.5, 0.5, "Reduced \\gx\\u2\\d");
                } else {
                    cpgsvp(0.82, 0.93, 0.09, 0.49);
                    cpgswin(0.001, 1.1 * max, y1l, y1h);
                    cpgline(search->numpdots, pdotchi, ftmparr1);
                    cpgsch(0.7);
                    cpgbox("BCMST", 0.0, 0, "BST", 0.0, 0);
                    cpgswin(0.001, 1.1 * max, y2l, y2h);
                    cpgbox("", 0.0, 0, "CMST", 0.0, 0);
                    cpgsch(0.8);
                    cpgmtxt("T", 4.2, 0.5, 0.5, "Reduced");
                    cpgmtxt("T", 2.8, 0.5, 0.5, "\\gx\\u2\\d");
                    cpgmtxt("R", 2.4, 0.5, 0.5, fdout);
                }
                vect_free(ftmparr1);

                /* P P-dot image */

                if (search->nsub > 1)
                    cpgsvp(0.74, 0.94, 0.09, 0.29);
                else
                    cpgsvp(0.51, 0.82, 0.09, 0.49);
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
                cpgimag(ppdot2d, nc, nr, 0 + 1, nc, 0 + 1, nr, bg, fg, tr);
                x1l = (float) ((bestp - pfold) * 1000.0);
                y1l = (float) (bestpd - pdfold);
                /* Plot the error bars on the P-Pdot diagram */
                cpgpt(1, &x1l, &y1l, 5);
                errlen = (float) (perr * 1000.0);
                cpgerrb(5, 1, &x1l, &y1l, &errlen, 2);
                errlen = (float) (pderr);
                cpgerrb(6, 1, &x1l, &y1l, &errlen, 2);
                if (search->nsub > 1) {
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
                vect_free(levels);
                cpgsch(0.8);
            }

            {
                char out[200], out2[100];

                /* Add the Data Info area */

                cpgsvp(0.06, 0.94, 0.09, 0.94);
                cpgsch(0.8);
                sprintf(out, "%-s", search->filenm);
                cpgmtxt("B", 4.0, 0.0, 0.0, out);
                cpgsch(1.0);
                cpgsvp(0.27, 0.519, 0.68, 0.94);
                cpgswin(-0.1, 1.00, -0.1, 1.1);
                cpgsch(0.7);
                sprintf(out, "Candidate:  %-s", search->candnm);
                cpgtext(0.0, 1.0, out);
                sprintf(out, "Telescope:  %-s", search->telescope);
                cpgtext(0.0, 0.9, out);
                if (TEST_EQUAL(search->tepoch, 0.0)
                    || TEST_EQUAL(search->tepoch, -1))
                    // -1.0 is for fake data made with makedata
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
                if (flags->events)
                    sprintf(out, "=  N/A (Events)");
                else
                    sprintf(out, "=  %.5g", search->dt);
                cpgtext(0.45, 0.6, out);
                if (flags->events)
                    sprintf(out, "Events Folded");
                else
                    sprintf(out, "Data Folded");
                cpgtext(0.0, 0.5, out);
                if (flags->events)
                    sprintf(out, "=  %-.0f", beststats.prof_avg * search->proflen);
                else
                    sprintf(out, "=  %-.0f", N);
                cpgtext(0.45, 0.5, out);
                sprintf(out, "Data Avg");
                cpgtext(0.0, 0.4, out);
                sprintf(out, "=  %.4g", beststats.data_avg);
                cpgtext(0.45, 0.4, out);
                sprintf(out, "Data StdDev");
                cpgtext(0.0, 0.3, out);
                sprintf(out, "=  %.4g", sqrt(beststats.data_var));
                cpgtext(0.45, 0.3, out);
                sprintf(out, "Profile Bins");
                cpgtext(0.0, 0.2, out);
                sprintf(out, "=  %d", search->proflen);
                cpgtext(0.45, 0.2, out);
                sprintf(out, "Profile Avg");
                cpgtext(0.0, 0.1, out);
                sprintf(out, "=  %.4g", beststats.prof_avg);
                cpgtext(0.45, 0.1, out);
                sprintf(out, "Profile StdDev");
                cpgtext(0.0, 0.0, out);
                sprintf(out, "=  %.4g", sqrt(beststats.prof_var));
                cpgtext(0.45, 0.0, out);

                {
                    int chidf;
                    double chip, chi_lnp, chi_sig;

                    chidf = (int)(dofeff + 1e-12);
                    chi_lnp = chi2_logp(beststats.redchi * chidf, chidf);
                    chip = (chi_lnp < -700) ? 0.0 : exp(chi_lnp);
                    chi_sig = equivalent_gaussian_sigma(chi_lnp);

                    /* Add the Fold Info area */

                    cpgsvp(0.519, 0.94, 0.68, 0.94);
                    cpgswin(-0.05, 1.05, -0.1, 1.1);
                    cpgsch(0.8);
                    cpgmtxt("T", 0.0, 0.5, 0.5, "Search Information");
                    cpgsch(0.7);
                    sprintf(out, "RA\\dJ2000\\u = %s", search->rastr);
                    cpgtext(0.0, 1.0, out);
                    sprintf(out, "DEC\\dJ2000\\u = %s", search->decstr);
                    cpgtext(0.6, 1.0, out);
                    if (flags->nosearch)
                        cpgtext(0.0, 0.9, "        Folding Parameters");
                    else
                        cpgtext(0.0, 0.9, "        Best Fit Parameters");
                    sprintf(out2, "(%.1f\\gs)", chi_sig);
                    if (chip == 0.0)
                        sprintf(out,
                                "DOF\\deff\\u = %.2f  \\gx\\u2\\d\\dred\\u = %.3f  P(Noise) ~ 0   %s",
                                dofeff, beststats.redchi, out2);
                    else
                        sprintf(out,
                                "DOF\\deff\\u = %.2f  \\gx\\u2\\d\\dred\\u = %.3f  P(Noise) < %.3g  %s",
                                dofeff, beststats.redchi, chip, out2);
                    cpgtext(0.0, 0.8, out);
                    if (search->nsub > 1) {
                        sprintf(out, "Dispersion Measure (DM; pc/cm\\u3\\d) = %.3f",
                                search->bestdm);
                        cpgtext(0.0, 0.7, out);
                    } else {
                        sprintf(out, "Dispersion Measure (DM) = N/A");
                        cpgtext(0.0, 0.7, out);
                    }
                    {
                        if (search->tepoch != 0.0) {
                            cpgnice_output_2(out2, search->topo.p1 * 1000.0,
                                             perr * 1000.0, 0);
                            sprintf(out, "P\\dtopo\\u (ms) = %s", out2);
                            cpgtext(0.0, 0.6, out);
                            cpgnice_output_2(out2, search->topo.p2, pderr, 0);
                            sprintf(out, "P'\\dtopo\\u (s/s) = %s", out2);
                            cpgtext(0.0, 0.5, out);
                            cpgnice_output_2(out2, search->topo.p3, pdderr, 0);
                            sprintf(out, "P''\\dtopo\\u (s/s\\u2\\d) = %s", out2);
                            cpgtext(0.0, 0.4, out);
                        } else {
                            cpgtext(0.0, 0.6, "P\\dtopo\\u (ms) = N/A");
                            cpgtext(0.0, 0.5, "P'\\dtopo\\u (s/s) = N/A");
                            cpgtext(0.0, 0.4, "P''\\dtopo\\u (s/s\\u2\\d) = N/A");
                        }
                        if (search->bepoch != 0.0) {
                            cpgnice_output_2(out2, search->bary.p1 * 1000.0,
                                             perr * 1000.0, 0);
                            sprintf(out, "P\\dbary\\u (ms) = %s", out2);
                            cpgtext(0.6, 0.6, out);
                            cpgnice_output_2(out2, search->bary.p2, pderr, 0);
                            sprintf(out, "P'\\dbary\\u (s/s) = %s", out2);
                            cpgtext(0.6, 0.5, out);
                            cpgnice_output_2(out2, search->bary.p3, pdderr, 0);
                            sprintf(out, "P''\\dbary\\u (s/s\\u2\\d) = %s", out2);
                            cpgtext(0.6, 0.4, out);
                        } else {
                            cpgtext(0.6, 0.6, "P\\dbary\\u (ms) = N/A");
                            cpgtext(0.6, 0.5, "P'\\dbary\\u (s/s) = N/A");
                            cpgtext(0.6, 0.4, "P''\\dbary\\u (s/s\\u2\\d) = N/A");
                        }
                    }
                    cpgtext(0.0, 0.3, "        Binary Parameters");
                    if (TEST_EQUAL(search->orb.p, 0.0)) {
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
                        sprintf(out, "\\gw (deg) = %f", search->orb.w);
                        cpgtext(0.6, 0.1, out);
                        sprintf(out, "T\\dperi\\u = %-.11f", search->orb.t);
                        cpgtext(0.0, 0.0, out);
                    }
                }
            }
        }
        cpgclos();
        if (ct == 0) {
            // Attempt to change the .ps into a nice .png using latex2html...
            int retval = 0;
            char *command = (char *) malloc(2 * strlen(search->pgdev) + 60);
            sprintf(command, "pstoimg -density 200 -antialias -flip cw "
                    "-quiet -type png -out %.*s.png %.*s",
                    (int) strlen(search->pgdev) - 7, search->pgdev,
                    (int) strlen(search->pgdev) - 4, search->pgdev);
            if ((retval = system(command))) {
                perror("Error running pstoimg in prepfold_plot()");
                printf("\n");
            }
            free(command);
        }
    }
    vect_free(bestprof);
    vect_free(timeprofs);
    vect_free(parttimes);
    vect_free(timechi);
    vect_free(periodchi);
    vect_free(pdotchi);
    if (ppdot == NULL)
        vect_free(ppdot2d);
    if (search->nsub > 1) {
        vect_free(dmprofs);
        vect_free(dmchi);
    }
}
