#include "presto.h"
#include "ctype.h"

static psrparams pulsardata[NP];
static int np = 0, have_database = 0;

static char num[41][5] = { "0th", "1st", "2nd", "3rd", "4th", "5th", "6th",
    "7th", "8th", "9th", "10th", "11th", "12th",
    "13th", "14th", "15th", "16th", "17th", "18th",
    "19th", "20th", "21st", "22nd", "23rd", "24th",
    "25th", "26th", "27th", "28th", "29th", "30th",
    "31st", "32nd", "33rd", "34th", "35th", "36th",
    "37th", "38th", "39th", "40th"
};


int read_database(void)
/* Reads the full pulsar database into the static array psrdata */
{
    FILE *database;
    char databasenm[200];
    psrdata pdata;
    binpsrdata bpdata;

    /* Open the binary data file */
    sprintf(databasenm, "%s/lib/pulsars.cat", getenv("PRESTO"));
    database = chkfopen(databasenm, "rb");

    while (chkfread(&pdata, sizeof(psrdata), 1, database)) {
        if (np >= NP) {
            printf("NP value set to small (%d) in $PRESTO/include/database.h\n", NP);
            printf("Please increase it and recompile\n");
            exit(-1);
        }
        strncpy(pulsardata[np].jname, pdata.jname, 13);
        strncpy(pulsardata[np].bname, pdata.bname, 9);
        strncpy(pulsardata[np].alias, pdata.alias, 10);
        pulsardata[np].ra2000 = pdata.ra2000;
        pulsardata[np].dec2000 = pdata.dec2000;
        pulsardata[np].dm = pdata.dm;
        pulsardata[np].timepoch = pdata.timepoch;
        pulsardata[np].p = pdata.p;
        pulsardata[np].pd = pdata.pd;
        pulsardata[np].pdd = 0.0;
        pulsardata[np].f = 1.0 / pdata.p;
        pulsardata[np].fd = -pdata.pd / (pdata.p * pdata.p);
        pulsardata[np].fdd = 0.0;
        if (pdata.binary == 1.0) {
            chkfread(&bpdata, sizeof(binpsrdata), 1, database);
            pulsardata[np].orb.p = bpdata.pb;
            pulsardata[np].orb.e = bpdata.e;
            pulsardata[np].orb.x = bpdata.x;
            pulsardata[np].orb.w = bpdata.w;
            pulsardata[np].orb.t = bpdata.To;
        } else {
            pulsardata[np].orb.p = 0.0;
            pulsardata[np].orb.e = 0.0;
            pulsardata[np].orb.x = 0.0;
            pulsardata[np].orb.w = 0.0;
            pulsardata[np].orb.t = 0.0;
        }
        pulsardata[np].orb.pd = 0.0;
        pulsardata[np].orb.wd = 0.0;
        np++;
    };

    /* Close the database */
    fclose(database);

    have_database = 1;
    return np;
}


void get_psr(int psrnumber, psrparams * psr)
/* Returns a full database entry for the pulsar #psrnumber in */
/* the database psrdata.  Returns *psr completed.             */
{
    int ii = psrnumber;

    if (ii < 0) {
        printf("psrnumber < 0 in get_psr().  Exiting.\n\n");
        exit(1);
    }
    strncpy(psr->jname, pulsardata[ii].jname, 13);
    strncpy(psr->bname, pulsardata[ii].bname, 9);
    strncpy(psr->alias, pulsardata[ii].alias, 10);
    psr->ra2000 = pulsardata[ii].ra2000;
    psr->dec2000 = pulsardata[ii].dec2000;
    psr->dm = pulsardata[ii].dm;
    psr->timepoch = pulsardata[ii].timepoch;
    psr->p = pulsardata[ii].p;
    psr->pd = pulsardata[ii].pd;
    psr->pdd = pulsardata[ii].pdd;
    psr->f = pulsardata[ii].f;
    psr->fd = pulsardata[ii].fd;
    psr->fdd = pulsardata[ii].fdd;
    psr->orb.p = pulsardata[ii].orb.p;
    psr->orb.e = pulsardata[ii].orb.e;
    psr->orb.x = pulsardata[ii].orb.x;
    psr->orb.w = pulsardata[ii].orb.w;
    psr->orb.t = pulsardata[ii].orb.t;
    psr->orb.pd = pulsardata[ii].orb.pd;
    psr->orb.wd = pulsardata[ii].orb.wd;
}


void get_psrparams(psrparams * psr, char *psrname)
/* Read a full pulsar database entry for pulsar psrname. */
/* Return the data in a psrparams structure.             */
{
    int pnum = 0;

    /* Read the database if needed */
    if (!have_database)
        np = read_database();

    /* Find the pulsar of interest */
    pnum = psr_number_from_name(psrname);

    /* Fill the structure with data */
    if (pnum >= 0)
        get_psr(pnum, psr);
    else {
        printf("Could not find the PSR in the database in get_psrparams().\n");
        exit(2);
    }
}


int psr_number_from_name(char *psrname)
/* Returns the pulsar number of psrname from the database */
/* This number can be from zero to the total number       */
/* of pulsars minus 1.  This way you can use this number  */
/* as an index from the result of collect_psrparams().    */
/* Return -1 if no pulsar is found.                       */
{
    int ii, psrnumber = -1;
    char *matchname, jname[13], bname[9];

    matchname = strlower(psrname);
    if (matchname[0] == 'j' || matchname[0] == 'b')
        matchname++;

    /* Read the database if needed */
    if (!have_database)
        np = read_database();

    /* Search for the J-name, the B-name, or the alias */
    for (ii = 0; ii < np; ii++) {
        strncpy(jname, pulsardata[ii].jname, 13);
        strncpy(bname, pulsardata[ii].bname, 9);
        if (!strcmp(strlower(jname), matchname) ||
            !strcmp(strlower(bname), matchname) ||
            !strcmp(pulsardata[ii].alias, matchname)) {
            psrnumber = ii;
            break;
        }
    }

    /* Return the pulsar number */

    return psrnumber;
}


int get_psr_at_epoch(char *psrname, double epoch, psrparams * psr)
/* Converts info from the pulsar database to "current" epoch.       */
/* Returned values go in *psr.                                      */
/* psrname is the pulsar name we are looking for (no J or B prefix) */
/* epoch is the time in question in MJD.                            */
/* The int returned is the number of the pulsar in the database.    */
/* If the int = 0, then no match was found.                         */
{
    int ii;
    double difft, f, fd;

    /* Get the pulsar's number from the database */

    ii = psr_number_from_name(psrname);

    /* Pulsar not in the database. */

    if (ii < 0)
        return 0;

    /* Pulsar is in the database. */

    else {
        get_psr(ii, psr);
        difft = SECPERDAY * (epoch - psr->timepoch);
        psr->timepoch = epoch;
        f = psr->f;
        fd = psr->fd;
        psr->f = f + fd * difft + 0.5 * psr->fdd * difft * difft;
        psr->fd = fd + psr->fdd * difft;
        psr->p = 1.0 / psr->f;
        psr->pd = -psr->fd * psr->p * psr->p;
        psr->pdd = (2.0 * (fd * fd) / f - psr->fdd) / (f * f);

        /* Selected pulsar is binary... */

        if (psr->orb.p != 0.0) {
            difft = SECPERDAY * (epoch - psr->orb.t);
            psr->orb.p = psr->orb.p * SECPERDAY + psr->orb.pd * difft;
            /* psr->orb.t is in seconds, _not_ MJD.  It represents the time    */
            /*      in sec _since_ the last periastron passage, _not_ when the */
            /*      next periastron will occur....                             */
            psr->orb.t = fmod(difft, psr->orb.p);
            if (psr->orb.t < 0.0)
                psr->orb.t += psr->orb.p;
            psr->orb.w =
                fmod((psr->orb.w + difft * psr->orb.wd / SECPERJULYR), 360.0);
        }
    }

    /* Return the number of the pulsar in the database. */

    return ii;
}


int comp_psr_to_cand(fourierprops * cand, infodata * idata, char *output, int full)
/* Compares a pulsar candidate defined by its properties found in   */
/*   *cand, and *idata with all of the pulsars in the pulsar        */
/*   database.  It returns a string (verbose if full==1) describing */
/*   the results of the search in *output.                          */
{
    int ii, jj;
    static infodata *old_idata;
    double theor, theoz, sidedr = 20.0;
    double r_criteria, z_criteria, rdiff, zdiff, difft = 0;
    static double T, beam2, ra, dec, epoch;
    char tmp1[80], tmp2[80], tmp3[80], shortout[30], psrname[20];
    rzwerrs rzws;

    /* Read the database if needed */

    if (!have_database)
        np = read_database();

    /* If calling for the first time for a certain data set, */
    /* initialize some values.                               */

    if (idata != old_idata) {
        /* Convert the beam width to radians */

        beam2 = 2.0 * ARCSEC2RAD * idata->fov;

        /* Convert RA and DEC to radians  (Use J2000) */

        ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
        dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);
        T = idata->N * idata->dt;
        epoch = (double) idata->mjd_i + idata->mjd_f + T / (2.0 * SECPERDAY);

        /* Set up old_idata for next time */

        old_idata = idata;
    }
    /* Calculate the measured r, z, w's and their derivatives */

    calc_rzwerrs(cand, T, &rzws);

    /* Run through RAs in database looking for things close  */
    /* If find one, check the DEC as well (the angle between */
    /* the sources < 2*beam diam).  If find one, check its   */
    /* period.  If this matches within 2*perr, return the    */
    /* number of the pulsar.  If no matches, return 0.       */

    for (ii = 0; ii < np; ii++) {

        /* See if we're close in RA */

        if (fabs(pulsardata[ii].ra2000 - ra) < 5 * beam2) {

            /* See if we're close in RA and DEC */

            if (sphere_ang_diff(pulsardata[ii].ra2000, pulsardata[ii].dec2000,
                                ra, dec) < 5 * beam2) {

                /* Predict the period of the pulsar at the observation MJD */

                difft = SECPERDAY * (epoch - pulsardata[ii].timepoch);
                theor = T / (pulsardata[ii].p + pulsardata[ii].pd * difft);
                theoz = -pulsardata[ii].pd * theor * theor;

                /* Check the predicted period and its harmonics against the */
                /* measured period.                                         */

                for (jj = 1; jj < 41; jj++) {

                    /* If the psr from the database is in a             */
                    /* binary orbit, loosen the match criteria.         */
                    /* This accounts for Doppler variations in period.  */

                    if (pulsardata[ii].orb.p != 0.0) {
                        r_criteria = 0.001 * theor * jj;        /* 0.1% fractional error   */
                        z_criteria = 9999999999.0;      /* Always match for binary */
                        strcpy(tmp1, "?");
                        if (full) {
                            strcpy(tmp3, "Possibly (large error) ");
                        }
                    } else {
                        r_criteria = 5.0;       /* 5 bin error matching... */
                        z_criteria = 9999999999.0;      /* Always match for binary */
                        /* z_criteria = 2.5 * cand->zerr; */
                        strcpy(tmp1, "");
                        if (full) {
                            strcpy(tmp3, "Looks like ");
                        }
                    }

                    if (theor * jj > 1.5 * cand->r)
                        break;

                    rdiff = fabs(theor * jj - cand->r);
                    zdiff = fabs(theoz * jj - cand->z);

                    if (rdiff < r_criteria && zdiff < z_criteria) {
                        if (strlen(pulsardata[ii].bname) == 0)
                            sprintf(psrname, "J%s", pulsardata[ii].jname);
                        else
                            sprintf(psrname, "B%s", pulsardata[ii].bname);
                        if (jj == 1) {
                            if (full) {
                                sprintf(tmp1, "the fundamental of ");
                                sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n",
                                        psrname, T / theor);
                                sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
                            } else {
                                sprintf(shortout, "PSR %s%s", psrname, tmp1);
                                strncpy(output, shortout, 20);
                            }
                        } else {
                            if (full) {
                                sprintf(tmp1, "the %s harmonic of ", num[jj]);
                                sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n",
                                        psrname, T / theor);
                                sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
                            } else {
                                sprintf(shortout, "%s H %s%s", num[jj], psrname,
                                        tmp1);
                                strncpy(output, shortout, 20);
                            }
                        }
                        return ii + 1;
                    } else if (rdiff < sidedr) {
                        if (strlen(pulsardata[ii].bname) == 0)
                            sprintf(psrname, "J%s", pulsardata[ii].jname);
                        else
                            sprintf(psrname, "B%s", pulsardata[ii].bname);
                        if (full) {
                            sprintf(tmp1, "a sidelobe of the %s harmonic of ",
                                    num[jj]);
                            sprintf(tmp2, "PSR %s. (predicted p = %11.7f s).\n",
                                    psrname, T / theor);
                            sprintf(output, "%s%s\n     %s", tmp3, tmp1, tmp2);
                        } else {
                            sprintf(shortout, "SL H%d %s", jj, psrname);
                            strncpy(output, shortout, 20);
                        }
                        return ii + 1;
                    }
                }
            }
        }
    }

    /* Didn't find a match */

    if (full) {
        sprintf(output,
                "I don't recognize this candidate in the pulsar database.\n");
    } else {
        strncpy(output, "                       ", 20);
    }
    return 0;
}


int comp_bin_to_cand(binaryprops * cand, infodata * idata, char *output, int full)
/* Compares a binary PSR candidate defined by its props found in    */
/*   *cand, and *idata with all of the pulsars in the pulsar        */
/*   database.  It returns a string (verbose if full==1) describing */
/*   the results of the search in *output.                          */
{
    int ii, jj, kk;
    double theop, ra, dec, beam2, difft = 0.0, epoch;
    double bmod, pmod;
    char tmp1[80], tmp2[80], tmp3[80], psrname[20], shortout[30];

    /* Read the database if needed */

    if (!have_database)
        np = read_database();

    /* Convert the beam width to radians */

    beam2 = 2.0 * ARCSEC2RAD * idata->fov;

    /* Convert RA and DEC to radians  (Use J2000) */

    ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
    dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);

    /* Calculate the time related variables  */

    epoch = (double) idata->mjd_i + idata->mjd_f;

    /* Run through RAs in database looking for things close  */
    /* If find one, check the DEC as well (the angle between */
    /* the sources < 2*beam diam).  If find one, check its   */
    /* period.  If this matches within 2*perr, return the    */
    /* number of the pulsar.  If no matches, return 0.       */

    for (ii = 0; ii < np; ii++) {

        /* See if we're close in RA */

        if (fabs(pulsardata[ii].ra2000 - ra) < beam2) {

            /* See if we're close in RA and DEC */

            if (sphere_ang_diff(pulsardata[ii].ra2000, pulsardata[ii].dec2000,
                                ra, dec) < beam2) {

                /* Check that the psr in the database is in a binary   */

                if (pulsardata[ii].orb.p != 0.0) {

                    /* Predict the period of the pulsar at the observation MJD */

                    difft = SECPERDAY * (epoch - pulsardata[ii].timepoch);
                    theop = pulsardata[ii].p + pulsardata[ii].pd * difft;

                    /* Check the predicted period and its harmonics against the */
                    /* measured period.  Use both pulsar and binary periods.    */

                    for (jj = 1, pmod = 1.0; jj < 41; jj++, pmod = 1.0 / (double) jj) {
                        if (fabs(theop * pmod - cand->ppsr) < (4 * cand->ppsrerr)) {
                            for (kk = 1, bmod = 1.0; kk < 10;
                                 kk++, bmod = 1.0 / (double) kk) {
                                if (fabs
                                    (pulsardata[ii].orb.p * bmod -
                                     cand->pbin / SECPERDAY) <
                                    (4 * cand->pbinerr / SECPERDAY)) {
                                    if (strlen(pulsardata[ii].bname) == 0)
                                        sprintf(psrname, "J%s",
                                                pulsardata[ii].jname);
                                    else
                                        sprintf(psrname, "B%s",
                                                pulsardata[ii].bname);
                                    if (jj > 1) {
                                        sprintf(tmp1,
                                                "Possibly the %s phasemod harmonic ",
                                                num[kk]);
                                        if (full) {
                                            sprintf(tmp2,
                                                    "of the %s harmonic of PSR ",
                                                    num[jj]);
                                            sprintf(tmp3,
                                                    "%s (p = %11.7f s, pbin = %9.4f d).\n",
                                                    psrname, theop,
                                                    pulsardata[ii].orb.p);
                                            sprintf(output, "%s%s%s", tmp1, tmp2,
                                                    tmp3);
                                        } else {
                                            sprintf(shortout, "%s H %s", num[kk],
                                                    psrname);
                                            strncpy(output, shortout, 20);
                                        }
                                    } else {
                                        if (full) {
                                            sprintf(tmp2,
                                                    "of PSR %s (p = %11.7f s, pbin = %9.4f d).\n",
                                                    psrname, theop,
                                                    pulsardata[ii].orb.p);
                                            sprintf(output, "%s%s", tmp1, tmp2);
                                        } else {
                                            sprintf(shortout, "PSR %s", psrname);
                                            strncpy(output, shortout, 20);
                                        }
                                    }
                                }
                                return ii + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Didn't find a match */

    if (full) {
        sprintf(output,
                "I don't recognize this candidate in the pulsar database.\n");
    } else {
        strncpy(output, "                       ", 20);
    }
    return 0;
}


int comp_rawbin_to_cand(rawbincand * cand, infodata * idata, char *output, int full)
/* Compares a binary PSR candidate defined by its props found in    */
/*   *cand, and *idata with all of the pulsars in the pulsar        */
/*   database.  It returns a string (verbose if full==1) describing */
/*   the results of the search in *output.                          */
{
    int ii, jj, kk;
    static int np;
    double theop, ra, dec, beam2, difft = 0.0, epoch;
    double bmod, pmod, orbperr, psrperr;
    char tmp1[80], tmp2[80], tmp3[80];

    /* Read the database if needed */

    if (!have_database)
        np = read_database();

    /* Convert the beam width to radians */

    beam2 = 2.0 * ARCSEC2RAD * idata->fov;

    /* Convert RA and DEC to radians  (Use J2000) */

    ra = hms2rad(idata->ra_h, idata->ra_m, idata->ra_s);
    dec = dms2rad(idata->dec_d, idata->dec_m, idata->dec_s);

    /* Calculate the time related variables  */

    epoch = (double) idata->mjd_i + idata->mjd_f;

    /* Calculate the approximate error in our value of orbital period */

    orbperr = 0.5 * cand->full_T / cand->mini_N;

    /* Calculate the approximate error in our value of spin period */

    if (cand->full_lo_r == 0.0)
        psrperr = cand->psr_p;
    else
        psrperr = fabs(cand->full_T /
                       (cand->full_lo_r + 0.5 * cand->mini_N) -
                       cand->full_T / cand->full_lo_r);

    /* Run through RAs in database looking for things close  */
    /* If find one, check the DEC as well (the angle between */
    /* the sources < 2*beam diam).  If find one, check its   */
    /* period.  If this matches within 2*perr, return the    */
    /* number of the pulsar.  If no matches, return 0.       */

    for (ii = 0; ii < np; ii++) {

        /* See if we're close in RA */

        if (fabs(pulsardata[ii].ra2000 - ra) < 5.0 * beam2) {

            /* See if we're close in RA and DEC */

            if (sphere_ang_diff(pulsardata[ii].ra2000, pulsardata[ii].dec2000,
                                ra, dec) < 5.0 * beam2) {

                /* Check that the psr in the database is in a binary   */

                if (pulsardata[ii].orb.p != 0.0) {

                    /* Predict the period of the pulsar at the observation MJD */

                    difft = SECPERDAY * (epoch - pulsardata[ii].timepoch);
                    theop = pulsardata[ii].p + pulsardata[ii].pd * difft;

                    /* Check the predicted period and its harmonics against the */
                    /* measured period.  Use both pulsar and binary periods.    */

                    for (jj = 1; jj < 41; jj++) {
                        pmod = 1.0 / (double) jj;
                        if (fabs(theop * pmod - cand->psr_p) < psrperr) {
                            for (kk = 1; kk < 10; kk++) {
                                bmod = (double) kk;
                                if (fabs
                                    (pulsardata[ii].orb.p * bmod -
                                     cand->orb_p / SECPERDAY) < orbperr) {
                                    if (strlen(pulsardata[ii].bname) == 0) {
                                        if (jj > 1) {
                                            if (full) {
                                                sprintf(tmp1,
                                                        "Possibly the %s phasemod harmonic ",
                                                        num[kk]);
                                                sprintf(tmp2,
                                                        "of the %s harmonic of PSR ",
                                                        num[jj]);
                                                sprintf(tmp3,
                                                        "J%s (p = %11.7f s, pbin = %9.4f d).\n",
                                                        pulsardata[ii].jname, theop,
                                                        pulsardata[ii].orb.p);
                                                sprintf(output, "%s%s%s", tmp1, tmp2,
                                                        tmp3);
                                            } else {
                                                sprintf(output, "%s H J%.12s",
                                                        num[kk],
                                                        pulsardata[ii].jname);
                                            }
                                        } else {
                                            if (full) {
                                                sprintf(tmp1,
                                                        "Possibly the %s phasemod harmonic ",
                                                        num[kk]);
                                                sprintf(tmp2,
                                                        "of PSR J%s (p = %11.7f s, pbin = %9.4f d).\n",
                                                        pulsardata[ii].jname, theop,
                                                        pulsardata[ii].orb.p);
                                                sprintf(output, "%s%s", tmp1, tmp2);
                                            } else {
                                                sprintf(output, "PSR J%.12s",
                                                        pulsardata[ii].jname);
                                            }
                                        }
                                    } else {
                                        if (jj > 1) {
                                            if (full) {
                                                sprintf(tmp1,
                                                        "Possibly the %s modulation harmonic ",
                                                        num[kk]);
                                                sprintf(tmp2,
                                                        "of the %s harmonic of PSR ",
                                                        num[jj]);
                                                sprintf(tmp3,
                                                        "B%s (p = %11.7f s, pbin = %9.4f d).\n",
                                                        pulsardata[ii].bname, theop,
                                                        pulsardata[ii].orb.p);
                                                sprintf(output, "%s%s%s", tmp1, tmp2,
                                                        tmp3);
                                            } else {
                                                sprintf(output, "%s H B%s", num[kk],
                                                        pulsardata[ii].bname);
                                            }
                                        } else {
                                            if (full) {
                                                sprintf(tmp1,
                                                        "Possibly the %s phasemod harmonic ",
                                                        num[kk]);
                                                sprintf(tmp2,
                                                        "of PSR B%s (p = %11.7f s, pbin = %9.4f d).\n",
                                                        pulsardata[ii].bname, theop,
                                                        pulsardata[ii].orb.p);
                                                sprintf(output, "%s%s", tmp1, tmp2);
                                            } else {
                                                sprintf(output, "PSR B%s",
                                                        pulsardata[ii].bname);
                                            }
                                        }
                                    }
                                }
                                return ii + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Didn't find a match */

    if (full) {
        sprintf(output,
                "I don't recognize this candidate in the pulsar database.\n");
    } else {
        sprintf(output, "                  ");
    }
    return 0;
}
