#include "prepfold.h"
#include "prepfold_cmd.h"
#include "f77.h"

int compare_doubles(const void *a, const void *b);

static void slaStringExport(const char *source_c, char *dest_f, int dest_len)
{
    int i;

    /* Check the supplied pointers. */
    if (!source_c || !dest_f)
        return;

    /* Copy the characters of the input C string to the output FORTRAN
       string, taking care not to go beyond the end of the FORTRAN
       string. */
    for (i = 0; source_c[i] && (i < dest_len); i++) {
        dest_f[i] = source_c[i];
    }

    /* Fill the rest of the output FORTRAN string with blanks. */
    for (; i < dest_len; i++)
        dest_f[i] = ' ';
}


F77_SUBROUTINE(dgels) (CHARACTER(TRANS),
                       INTEGER(M),
                       INTEGER(N),
                       INTEGER(NRHS),
                       DOUBLE_ARRAY(A),
                       INTEGER(LDA),
                       DOUBLE_ARRAY(B),
                       INTEGER(LDB),
                       DOUBLE_ARRAY(WORK), INTEGER(LWORK), INTEGER(INFO)
                       TRAIL(TRANS));

void call_dgels(char *trans, int m, int n, int nrhs,
                double *a, int lda, double *b, int ldb,
                double *work, int lwork, int *info)
{
    DECLARE_CHARACTER(TRANS, 1);
    DECLARE_INTEGER(M);
    DECLARE_INTEGER(N);
    DECLARE_INTEGER(NRHS);
    DECLARE_INTEGER(LDA);
    DECLARE_INTEGER(LDB);
    DECLARE_INTEGER(LWORK);
    DECLARE_INTEGER(INFO);

    slaStringExport(trans, TRANS, 1);
    M = m;
    N = n;
    NRHS = nrhs;
    LDA = lda;
    LDB = ldb;
    LWORK = lwork;
    F77_CALL(dgels) (CHARACTER_ARG(TRANS),
                     INTEGER_ARG(&M),
                     INTEGER_ARG(&N),
                     INTEGER_ARG(&NRHS),
                     DOUBLE_ARRAY_ARG(a),
                     INTEGER_ARG(&LDA),
                     DOUBLE_ARRAY_ARG(b),
                     INTEGER_ARG(&LDB),
                     DOUBLE_ARRAY_ARG(work), INTEGER_ARG(&LWORK), INTEGER_ARG(&INFO)
                     TRAIL_ARG(TRANS));
    *info = INFO;
}

double switch_pfdot(double pf, double pfdot)
{
    double retval;

    if (pf == 0.0)
        return 0.0;
    else {
        retval = -pfdot / (pf * pf);
        if (retval == -0)
            return 0.0;
        else
            return retval;
    }
}

double switch_pfdotdot(double pf, double pfdot, double pfdotdot)
{
    double retval;

    if (pf == 0.0 || pfdotdot == 0.0)
        return 0.0;
    else {
        retval = 2.0 * pfdot * pfdot / (pf * pf * pf) - pfdotdot / (pf * pf);
        if (retval == -0)
            return 0.0;
        else
            return retval;
    }
}

double fdot2phasedelay(double fdot, double time)
{
    double retval;

    retval = fdot * time * time / 2.0;
    if (retval == -0)
        return 0.0;
    else
        return retval;
}

double fdotdot2phasedelay(double fdotdot, double time)
{
    double retval;

    retval = fdotdot * time * time * time / 6.0;
    if (retval == -0)
        return 0.0;
    else
        return retval;
}

double phasedelay2fdot(double phasedelay, double time)
{
    if (time == 0.0)
        return 0.0;
    else
        return 2.0 * phasedelay / (time * time);
}

double phasedelay2fdotdot(double phasedelay, double time)
{
    if (time == 0.0)
        return 0.0;
    else
        return 6.0 * phasedelay / (time * time * time);
}

void set_posn(prepfoldinfo * in, infodata * idata)
{
    if (strcmp(idata->telescope, "None (Artificial Data Set)") != 0) {
        ra_dec_to_string(in->rastr, idata->ra_h, idata->ra_m, idata->ra_s);
        ra_dec_to_string(in->decstr, idata->dec_d, idata->dec_m, idata->dec_s);
    }
}

int read_floats(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
    return chkfread(data, sizeof(float),
                    (unsigned long) (numpts * numchan), file) / numchan;
}


int read_shorts(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains short integer data.  */
/* The equivalent floats are placed in *data.               */
/* It returns the number of points read.                    */
{
    short *sdata;
    int ii, numread;

    sdata = (short *) malloc((size_t) (sizeof(short) * (numpts * numchan)));
    if (!sdata) {
        perror("\nError allocating short array in read_shorts()");
        printf("\n");
        exit(-1);
    }
    numread = chkfread(sdata, sizeof(short),
                       (unsigned long) (numpts * numchan), file) / numchan;
    for (ii = 0; ii < numread; ii++)
        data[ii] = (float) sdata[ii];
    free(sdata);
    return numread;
}


int read_PRESTO_subbands(FILE * infiles[], int numfiles, float *subbanddata,
                         double timeperblk, int *maskchans,
                         int *nummasked, mask * obsmask, float *padvals)
/* Read short int subband data written by prepsubband     */
/* Note:  This version returns a transpose of the version */
/*        listed in prepsubband.c                         */
{
    int ii, jj, index, numread = 0;
    short subsdata[SUBSBLOCKLEN];
    static int currentblock = 0;

    /* Read the data */
    for (ii = 0; ii < numfiles; ii++) {
        index = ii * SUBSBLOCKLEN;
        numread = chkfread(subsdata, sizeof(short), SUBSBLOCKLEN, infiles[ii]);
        for (jj = 0; jj < numread; jj++, index++)
            subbanddata[index] = (float) subsdata[jj];
        for (jj = numread; jj < SUBSBLOCKLEN; jj++, index++)
            subbanddata[index] = 0.0;
    }

    /* Mask it if required */
    if (obsmask->numchan && numread) {
        double starttime;
        starttime = currentblock * timeperblk;
        *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
        if (*nummasked == -1) { /* If all channels are masked */
            for (ii = 0; ii < numfiles; ii++) {
                index = ii * SUBSBLOCKLEN;
                for (jj = 0; jj < SUBSBLOCKLEN; jj++, index++)
                    subbanddata[index] = padvals[ii];
            }
        } else if (*nummasked > 0) {    /* Only some of the channels are masked */
            int channum;
            for (ii = 0; ii < *nummasked; ii++) {
                channum = maskchans[ii];
                index = channum * SUBSBLOCKLEN;
                for (jj = 0; jj < SUBSBLOCKLEN; jj++, index++)
                    subbanddata[index] = padvals[channum];
            }
        }
    }
    currentblock += 1;
    return numread;
}


double *read_events(FILE * infile, int bin, int days, int *numevents,
                    double MJD0, double Ttot, double startfrac, double endfrac,
                    double offset)
/* This routine reads a set of events from the open file 'infile'.     */
/* It returns a double precision vector of events in seconds from the  */
/* first event.  If 'bin' is true the routine treats the data as       */
/* binary double precision (otherwise text).  If 'days' is 1 then the  */
/* data is assumed to be in days since the 'inf' EPOCH (0 is sec from  */
/* EPOCH in 'inf').  If 'days' is 2, the data are assumed to be MJDs.  */
/* The number of events read is placed in 'numevents', and the raw     */
/* event is placed in 'firstevent'.  MJD0 is the time to use for the   */
/* reference time.  Ttot is the duration of the observation. 'start'   */
/* and 'end' are define the fraction of the observation that we are    */
/* interested in.  'offset' is a time offset to apply to the events.   */
{
    int N = 0, nn = 0, goodN = 0;
    double *ts, *goodts, dtmp, lotime, hitime;
    char line[80], *sptr = NULL;

    if (bin) {
        N = chkfilelen(infile, sizeof(double));
    } else {
        /* Read the input file once to count events */
        while (1) {
            sptr = fgets(line, 80, infile);
            if (!feof(infile) && sptr != NULL && sptr[0] != '\n') {
                if (line[0] != '#' && sscanf(line, "%lf", &dtmp) == 1)
                    N++;
            } else {
                break;
            }
        }
    }

    /* Allocate the event arrays */

    ts = (double *) malloc(N * sizeof(double));

    /* Rewind and read the events for real */

    rewind(infile);
    if (bin) {
        chkfread(ts, sizeof(double), N, infile);
    } else {
        while (1) {
            sptr = fgets(line, 80, infile);
            if (!feof(infile) && sptr != NULL && sptr[0] != '\n') {
                if (line[0] != '#' && sscanf(line, "%lf", &ts[nn]) == 1)
                    nn++;
            } else {
                break;
            }
        }
    }

    /* Sort the events  */

    qsort(ts, N, sizeof(double), compare_doubles);

    /* If there is no offset specified and the data are non-MJD */
    /* days or seconds, then set the offset to be the first event */

    if (offset == 0.0 && days < 2)
        offset = -ts[0];

    /* Convert all the events to MJD */

    if (days == 0) {            /* Events are in seconds since MJD0 */
        for (nn = 0; nn < N; nn++)
            ts[nn] = MJD0 + (ts[nn] + offset) / SECPERDAY;
    } else if (days == 1) {     /* Events are in days since MJD0 */
        for (nn = 0; nn < N; nn++)
            ts[nn] = MJD0 + (ts[nn] + offset);
    } else if (days == 2 && offset != 0.0) {    /* Events are in MJD with an offset */
        for (nn = 0; nn < N; nn++)
            ts[nn] += offset;
    }

    /* Count how many events are within our range and only keep them */

    lotime = MJD0 + startfrac * Ttot / SECPERDAY;
    hitime = MJD0 + endfrac * Ttot / SECPERDAY;
    for (nn = 0; nn < N; nn++)
        if (ts[nn] >= lotime && ts[nn] < hitime)
            goodN++;
    if (goodN != N) {
        goodts = (double *) malloc(goodN * sizeof(double));
        goodN = 0;
        for (nn = 0; nn < N; nn++) {
            if (ts[nn] >= lotime && ts[nn] < hitime) {
                goodts[goodN] = ts[nn];
                goodN++;
            }
        }
        free(ts);
        ts = goodts;
        N = goodN;
    } else {
        goodts = ts;
    }
    *numevents = N;

    /* Convert the events to seconds from MJD0 */

    for (nn = 0; nn < N; nn++)
        goodts[nn] = (goodts[nn] - MJD0) * SECPERDAY;

    return goodts;
}



int bary2topo(double *topotimes, double *barytimes, int numtimes,
              double fb, double fbd, double fbdd,
              double *ft, double *ftd, double *ftdd)
/* Convert a set of barycentric pulsar spin parameters (fb, fbd, fbdd) */
/* into topocentric spin parameters (ft, ftd, ftdd) by performing      */
/* a linear least-squares fit (using LAPACK routine DGELS).  The       */
/* routine equates the pulse phase using topcentric parameters and     */
/* times to the pulse phase using barycentric parameters and times.    */
{
    double *work, *aa, *bb, dtmp;
    int ii, mm = 3, nn, nrhs = 1, lwork, info, index;
    char trans = 'T';

    if (numtimes < 4) {
        printf("\n'numtimes' < 4 in bary2topo():  Cannot solve.\n\n");
        exit(0);
    }
    nn = numtimes;
    lwork = mm + nn * 9;
    aa = gen_dvect(mm * nn);
    bb = gen_dvect(nn);
    work = gen_dvect(lwork);
    for (ii = 0; ii < nn; ii++) {
        index = ii * 3;
        dtmp = (topotimes[ii] - topotimes[0]) * SECPERDAY;
        aa[index] = dtmp;
        aa[index + 1] = 0.5 * dtmp * dtmp;
        aa[index + 2] = dtmp * dtmp * dtmp / 6.0;
        dtmp = (barytimes[ii] - barytimes[0]) * SECPERDAY;
        bb[ii] = dtmp * (fb + dtmp * (0.5 * fbd + fbdd * dtmp / 6.0));
    }
    // dgels_(&trans, &mm, &nn, &nrhs, aa, &mm, bb, &nn, work, &lwork, &info);
    call_dgels(&trans, mm, nn, nrhs, aa, mm, bb, nn, work, lwork, &info);
    *ft = bb[0];
    *ftd = bb[1];
    *ftdd = bb[2];
    vect_free(aa);
    vect_free(bb);
    vect_free(work);
    return info;
}


void init_prepfoldinfo(prepfoldinfo * in)
/* Set all values to 0 or NULL */
{
    in->rawfolds = NULL;
    in->dms = NULL;
    in->periods = NULL;
    in->pdots = NULL;
    in->stats = NULL;
    in->numdms = 0;
    in->numperiods = 0;
    in->numpdots = 0;
    in->nsub = 0;
    in->npart = 0;
    in->proflen = 0;
    in->numchan = 1;
    in->ndmfact = 2;
    in->npfact = 1;
    in->pstep = 1;
    in->pdstep = 1;
    in->dmstep = 1;
    in->filenm = NULL;
    in->candnm = NULL;
    in->telescope = NULL;
    in->pgdev = NULL;
    {
        int ii;
        for (ii = 0; ii < 16; ii++) {
            in->rastr[ii] = '\0';
            in->decstr[ii] = '\0';
        }
        strcpy(in->rastr, "Unknown");
        strcpy(in->decstr, "Unknown");
    }
    in->dt = 0.0;
    in->startT = 0.0;
    in->endT = 0.0;
    in->tepoch = 0.0;
    in->bepoch = 0.0;
    in->avgvoverc = 0.0;
    in->lofreq = 0.0;
    in->chan_wid = 0.0;
    in->bestdm = 0.0;
    in->topo.pow = 0.0;
    in->topo.p1 = 0.0;
    in->topo.p2 = 0.0;
    in->topo.p3 = 0.0;
    in->bary.pow = 0.0;
    in->bary.p1 = 0.0;
    in->bary.p2 = 0.0;
    in->bary.p3 = 0.0;
    in->fold.pow = 0.0;
    in->fold.p1 = 0.0;
    in->fold.p2 = 0.0;
    in->fold.p3 = 0.0;
    in->orb.p = 0.0;
    in->orb.e = 0.0;
    in->orb.x = 0.0;
    in->orb.w = 0.0;
    in->orb.t = 0.0;
    in->orb.pd = 0.0;
    in->orb.wd = 0.0;
}

void print_prepfoldinfo(prepfoldinfo * in)
/* Print a prepfoldinfo data structure to STDOUT */
{
    printf("\n\n");
    printf("numdms      =  %d\n", in->numdms);
    printf("numperiods  =  %d\n", in->numperiods);
    printf("numpdots    =  %d\n", in->numpdots);
    printf("nsub        =  %d\n", in->nsub);
    printf("npart       =  %d\n", in->npart);
    printf("proflen     =  %d\n", in->proflen);
    printf("numchan     =  %d\n", in->numchan);
    printf("pstep       =  %d\n", in->pstep);
    printf("pdstep      =  %d\n", in->pdstep);
    printf("dmstep      =  %d\n", in->dmstep);
    printf("ndmfact     =  %d\n", in->ndmfact);
    printf("npfact      =  %d\n", in->npfact);
    printf("filenm      =  '%s'\n", in->filenm);
    printf("candnm      =  '%s'\n", in->candnm);
    printf("telescope   =  '%s'\n", in->telescope);
    printf("pgdev       =  '%s'\n", in->pgdev);
    printf("rastr       =  '%s'\n", in->rastr);
    printf("decstr      =  '%s'\n", in->decstr);
    printf("dt          =  %.14g\n", in->dt);
    printf("startT      =  %.14g\n", in->startT);
    printf("endT        =  %.14g\n", in->endT);
    printf("tepoch      =  %.14g\n", in->tepoch);
    printf("bepoch      =  %.14g\n", in->bepoch);
    printf("avgvoverc   =  %.14g\n", in->avgvoverc);
    printf("lofreq      =  %.14g\n", in->lofreq);
    printf("chan_wid    =  %.14g\n", in->chan_wid);
    printf("bestdm      =  %.14g\n", in->bestdm);
    printf("topo.pow    =  %.14g\n", in->topo.pow);
    printf("topo.p1     =  %.14g\n", in->topo.p1);
    printf("topo.p2     =  %.14g\n", in->topo.p2);
    printf("topo.p3     =  %.14g\n", in->topo.p3);
    printf("bary.pow    =  %.14g\n", in->bary.pow);
    printf("bary.p1     =  %.14g\n", in->bary.p1);
    printf("bary.p2     =  %.14g\n", in->bary.p2);
    printf("bary.p3     =  %.14g\n", in->bary.p3);
    printf("fold.pow    =  %.14g\n", in->fold.pow);
    printf("fold.p1     =  %.14g\n", in->fold.p1);
    printf("fold.p2     =  %.14g\n", in->fold.p2);
    printf("fold.p3     =  %.14g\n", in->fold.p3);
    printf("orb.p       =  %.14g\n", in->orb.p);
    printf("orb.e       =  %.14g\n", in->orb.e);
    printf("orb.x       =  %.14g\n", in->orb.x);
    printf("orb.w       =  %.14g\n", in->orb.w);
    printf("orb.t       =  %.14g\n", in->orb.t);
    printf("orb.pd      =  %.14g\n", in->orb.pd);
    printf("orb.wd      =  %.14g\n", in->orb.wd);
    printf("\n\n");
}

void write_prepfoldinfo(prepfoldinfo * in, char *filename)
/* Write a prepfoldinfo data structure to a binary file */
{
    FILE *outfile;
    int itmp;

    outfile = chkfopen(filename, "wb");
    chkfwrite(&in->numdms, sizeof(int), 1, outfile);
    chkfwrite(&in->numperiods, sizeof(int), 1, outfile);
    chkfwrite(&in->numpdots, sizeof(int), 1, outfile);
    chkfwrite(&in->nsub, sizeof(int), 1, outfile);
    chkfwrite(&in->npart, sizeof(int), 1, outfile);
    chkfwrite(&in->proflen, sizeof(int), 1, outfile);
    chkfwrite(&in->numchan, sizeof(int), 1, outfile);
    chkfwrite(&in->pstep, sizeof(int), 1, outfile);
    chkfwrite(&in->pdstep, sizeof(int), 1, outfile);
    chkfwrite(&in->dmstep, sizeof(int), 1, outfile);
    chkfwrite(&in->ndmfact, sizeof(int), 1, outfile);
    chkfwrite(&in->npfact, sizeof(int), 1, outfile);
    itmp = strlen(in->filenm);
    chkfwrite(&itmp, sizeof(int), 1, outfile);
    chkfwrite(in->filenm, sizeof(char), itmp, outfile);
    itmp = strlen(in->candnm);
    chkfwrite(&itmp, sizeof(int), 1, outfile);
    chkfwrite(in->candnm, sizeof(char), itmp, outfile);
    itmp = strlen(in->telescope);
    chkfwrite(&itmp, sizeof(int), 1, outfile);
    chkfwrite(in->telescope, sizeof(char), itmp, outfile);
    itmp = strlen(in->pgdev);
    chkfwrite(&itmp, sizeof(int), 1, outfile);
    chkfwrite(in->pgdev, sizeof(char), itmp, outfile);
    chkfwrite(in->rastr, sizeof(char), 16, outfile);
    chkfwrite(in->decstr, sizeof(char), 16, outfile);
    chkfwrite(&in->dt, sizeof(double), 1, outfile);
    chkfwrite(&in->startT, sizeof(double), 1, outfile);
    chkfwrite(&in->endT, sizeof(double), 1, outfile);
    chkfwrite(&in->tepoch, sizeof(double), 1, outfile);
    chkfwrite(&in->bepoch, sizeof(double), 1, outfile);
    chkfwrite(&in->avgvoverc, sizeof(double), 1, outfile);
    chkfwrite(&in->lofreq, sizeof(double), 1, outfile);
    chkfwrite(&in->chan_wid, sizeof(double), 1, outfile);
    chkfwrite(&in->bestdm, sizeof(double), 1, outfile);
    chkfwrite(&(in->topo.pow), sizeof(double), 1, outfile);
    chkfwrite(&(in->topo.p1), sizeof(double), 1, outfile);
    chkfwrite(&(in->topo.p2), sizeof(double), 1, outfile);
    chkfwrite(&(in->topo.p3), sizeof(double), 1, outfile);
    chkfwrite(&(in->bary.pow), sizeof(double), 1, outfile);
    chkfwrite(&(in->bary.p1), sizeof(double), 1, outfile);
    chkfwrite(&(in->bary.p2), sizeof(double), 1, outfile);
    chkfwrite(&(in->bary.p3), sizeof(double), 1, outfile);
    chkfwrite(&(in->fold.pow), sizeof(double), 1, outfile);
    chkfwrite(&(in->fold.p1), sizeof(double), 1, outfile);
    chkfwrite(&(in->fold.p2), sizeof(double), 1, outfile);
    chkfwrite(&(in->fold.p3), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.p), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.e), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.x), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.w), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.t), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.pd), sizeof(double), 1, outfile);
    chkfwrite(&(in->orb.wd), sizeof(double), 1, outfile);
    chkfwrite(in->dms, sizeof(double), in->numdms, outfile);
    chkfwrite(in->periods, sizeof(double), in->numperiods, outfile);
    chkfwrite(in->pdots, sizeof(double), in->numpdots, outfile);
    chkfwrite(in->rawfolds, sizeof(double), in->nsub *
              in->npart * in->proflen, outfile);
    chkfwrite(in->stats, sizeof(foldstats), in->nsub * in->npart, outfile);
    fclose(outfile);
}

void read_prepfoldinfo(prepfoldinfo * in, char *filename)
/* Read a prepfoldinfo data structure from a binary file */
{
    FILE *infile;
    int itmp, byteswap = 0;
    char temp[16];

    infile = chkfopen(filename, "rb");
    in->numdms = read_int(infile, byteswap);
    in->numperiods = read_int(infile, byteswap);
    in->numpdots = read_int(infile, byteswap);
    in->nsub = read_int(infile, byteswap);
    in->npart = read_int(infile, byteswap);
    /* The following is not exactly the most robust, but it should work... */
    if (in->npart < 1 || in->npart > 10000) {
        byteswap = 1;
        in->numdms = swap_int(in->numdms);
        in->numperiods = swap_int(in->numperiods);
        in->numpdots = swap_int(in->numpdots);
        in->nsub = swap_int(in->nsub);
        in->npart = swap_int(in->npart);
    }
    in->proflen = read_int(infile, byteswap);
    in->numchan = read_int(infile, byteswap);
    in->pstep = read_int(infile, byteswap);
    in->pdstep = read_int(infile, byteswap);
    in->dmstep = read_int(infile, byteswap);
    in->ndmfact = read_int(infile, byteswap);
    in->npfact = read_int(infile, byteswap);
    itmp = read_int(infile, byteswap);
    in->filenm = calloc(itmp + 1, sizeof(char));
    chkfread(in->filenm, sizeof(char), itmp, infile);
    itmp = read_int(infile, byteswap);
    in->candnm = calloc(itmp + 1, sizeof(char));
    chkfread(in->candnm, sizeof(char), itmp, infile);
    itmp = read_int(infile, byteswap);
    in->telescope = calloc(itmp + 1, sizeof(char));
    chkfread(in->telescope, sizeof(char), itmp, infile);
    itmp = read_int(infile, byteswap);
    in->pgdev = calloc(itmp + 1, sizeof(char));
    chkfread(in->pgdev, sizeof(char), itmp, infile);
    //chkfread(in->rastr, sizeof(char), 16, infile);
    {
        int has_posn = 1, ii;
        chkfread(temp, sizeof(char), 16, infile);
        /* Check to see if a position string was written */
        for (ii = 0; ii < 16; ii++) {
            if (!isdigit(temp[ii]) &&
                temp[ii] != ':' &&
                temp[ii] != '.' && temp[ii] != '-' && temp[ii] != '\0') {
                has_posn = 0;
                break;
            }
        }
        if (has_posn) {
            strcpy(in->rastr, temp);
            chkfread(in->decstr, sizeof(char), 16, infile);
            in->dt = read_double(infile, byteswap);
            in->startT = read_double(infile, byteswap);
        } else {
            strcpy(in->rastr, "Unknown");
            strcpy(in->decstr, "Unknown");
            in->dt = *(double *) (temp + 0);
            if (byteswap)
                in->dt = swap_double(in->dt);
            in->startT = *(double *) (temp + sizeof(double));
            if (byteswap)
                in->startT = swap_double(in->startT);
        }
    }
    in->endT = read_double(infile, byteswap);
    in->tepoch = read_double(infile, byteswap);
    in->bepoch = read_double(infile, byteswap);
    in->avgvoverc = read_double(infile, byteswap);
    in->lofreq = read_double(infile, byteswap);
    in->chan_wid = read_double(infile, byteswap);
    in->bestdm = read_double(infile, byteswap);
    /* The .pow elements were written as doubles (Why??) */
    in->topo.pow = read_float(infile, byteswap);
    read_float(infile, byteswap);
    in->topo.p1 = read_double(infile, byteswap);
    in->topo.p2 = read_double(infile, byteswap);
    in->topo.p3 = read_double(infile, byteswap);
    /* The .pow elements were written as doubles (Why??) */
    in->bary.pow = read_float(infile, byteswap);
    read_float(infile, byteswap);
    in->bary.p1 = read_double(infile, byteswap);
    in->bary.p2 = read_double(infile, byteswap);
    in->bary.p3 = read_double(infile, byteswap);
    /* The .pow elements were written as doubles (Why??) */
    in->fold.pow = read_float(infile, byteswap);
    read_float(infile, byteswap);
    in->fold.p1 = read_double(infile, byteswap);
    in->fold.p2 = read_double(infile, byteswap);
    in->fold.p3 = read_double(infile, byteswap);
    in->orb.p = read_double(infile, byteswap);
    in->orb.e = read_double(infile, byteswap);
    in->orb.x = read_double(infile, byteswap);
    in->orb.w = read_double(infile, byteswap);
    in->orb.t = read_double(infile, byteswap);
    in->orb.pd = read_double(infile, byteswap);
    in->orb.wd = read_double(infile, byteswap);
    in->dms = gen_dvect(in->numdms);
    chkfread(in->dms, sizeof(double), in->numdms, infile);
    in->periods = gen_dvect(in->numperiods);
    chkfread(in->periods, sizeof(double), in->numperiods, infile);
    in->pdots = gen_dvect(in->numpdots);
    chkfread(in->pdots, sizeof(double), in->numpdots, infile);
    in->rawfolds = gen_dvect(in->nsub * in->npart * in->proflen);
    chkfread(in->rawfolds, sizeof(double), in->nsub *
             in->npart * in->proflen, infile);
    in->stats = (foldstats *) malloc(sizeof(foldstats) * in->nsub * in->npart);
    chkfread(in->stats, sizeof(foldstats), in->nsub * in->npart, infile);
    fclose(infile);
    if (byteswap) {
        int ii;
        for (ii = 0; ii < in->numdms; ii++)
            in->dms[ii] = swap_double(in->dms[ii]);
        for (ii = 0; ii < in->numperiods; ii++)
            in->periods[ii] = swap_double(in->periods[ii]);
        for (ii = 0; ii < in->numpdots; ii++)
            in->pdots[ii] = swap_double(in->pdots[ii]);
        for (ii = 0; ii < in->nsub * in->npart * in->proflen; ii++)
            in->rawfolds[ii] = swap_double(in->rawfolds[ii]);
        for (ii = 0; ii < in->nsub * in->npart; ii++) {
            in->stats[ii].numdata = swap_double(in->stats[ii].numdata);
            in->stats[ii].data_avg = swap_double(in->stats[ii].data_avg);
            in->stats[ii].data_var = swap_double(in->stats[ii].data_var);
            in->stats[ii].numprof = swap_double(in->stats[ii].numprof);
            in->stats[ii].prof_avg = swap_double(in->stats[ii].prof_avg);
            in->stats[ii].prof_var = swap_double(in->stats[ii].prof_var);
            in->stats[ii].redchi = swap_double(in->stats[ii].redchi);
        }
    }
}

void delete_prepfoldinfo(prepfoldinfo * in)
/* Free all dynamic arrays in the prepfold array */
{
    vect_free(in->rawfolds);
    if (in->nsub > 1)
        vect_free(in->dms);
    vect_free(in->periods);
    vect_free(in->pdots);
    free(in->stats);
    free(in->filenm);
    free(in->candnm);
    free(in->telescope);
    free(in->pgdev);
}

void double2float(double *in, float *out, int numpts)
/* Copy a double vector into a float vector */
{
    int ii;

    for (ii = 0; ii < numpts; ii++)
        out[ii] = (float) in[ii];
}


void correct_subbands_for_DM(double dm, prepfoldinfo * search,
                             double *ddprofs, foldstats * ddstats)
/* Calculate the DM delays and apply them to the subbands */
/* to create de-dispersed profiles.                      */
{
    int ii, *dmdelays;
    double *subbanddelays, hif, dopplerhif, hifdelay, rdphase;

    rdphase = search->fold.p1 * search->proflen;
    hif = search->lofreq + (search->numchan - 1.0) * search->chan_wid;
    dopplerhif = doppler(hif, search->avgvoverc);
    hifdelay = delay_from_dm(dm, dopplerhif);
    subbanddelays = subband_delays(search->numchan, search->nsub, dm,
                                   search->lofreq, search->chan_wid,
                                   search->avgvoverc);
    dmdelays = gen_ivect(search->nsub);
    for (ii = 0; ii < search->nsub; ii++)
        dmdelays[ii] =
            NEAREST_INT((subbanddelays[ii] - hifdelay) * rdphase) % search->proflen;
    vect_free(subbanddelays);
    combine_subbands(search->rawfolds, search->stats, search->npart,
                     search->nsub, search->proflen, dmdelays, ddprofs, ddstats);
    vect_free(dmdelays);
}


double DOF_corr(double dt_per_bin)
// Return a multiplicative correction for the effective number of
// degrees of freedom in the chi^2 measurement resulting from a pulse
// profile folded by PRESTO's fold() function (i.e. prepfold).  This
// is required because there are correlations between the bins caused
// by the way that prepfold folds data (i.e. treating a sample as
// finite duration and smearing it over potenitally several bins in
// the profile as opposed to instantaneous and going into just one
// profile bin).  The correction is semi-analytic (thanks to Paul
// Demorest and Walter Brisken) but the values for 'power' and
// 'factor' have been determined from Monte Carlos.  The correction is
// good to a fractional error of less than a few percent as long as
// dt_per_bin is > 0.5 or so (which it usually is for pulsar
// candidates).  There is a very minimal number-of-bins dependence,
// which is apparent when dt_per_bin < 0.7 or so.  dt_per_bin is the
// width of a profile bin in samples (a float), and so for prepfold is
// pulse period / nbins / sample time.  Note that the sqrt of this
// factor can be used to 'inflate' the RMS of the profile as well, for
// radiometer eqn flux density estimates, for instance.
{
    double power = 1.806;       // From Monte Carlos
    double factor = 0.96;       // From Monte Carlos
    return dt_per_bin * factor * pow(1.0 + pow(dt_per_bin, power), -1.0 / power);
}


float estimate_offpulse_redchi2(double *inprofs, foldstats * stats,
                                int numparts, int numsubbands,
                                int proflen, int numtrials, double dofeff)
// Randomly offset each pulse profile in a .pfd data square or cube
// and combine them to estimate a "true" off-pulse level.  Do this
// numtrials times in order to improve the statistics.  Return the
// inverse of the average of the off-pulse reduced-chi^2 (i.e. the
// correction factor).  dofeff is the effective number of DOF as
// returned by DOF_corr().
{
    int ii, jj, kk, offset, trialnum, phsindex, statindex;
    float *chis;
    double chi_avg, chi_var, redchi;
    double prof_avg, prof_var, *prof_ptr, *sumprof;

    sumprof = gen_dvect(proflen);
    chis = gen_fvect(numtrials);

    for (trialnum = 0; trialnum < numtrials; trialnum++) {
        // Initialize the summed profile
        for (ii = 0; ii < proflen; ii++)
            sumprof[ii] = 0.0;
        prof_avg = 0.0;
        prof_var = 0.0;
        prof_ptr = inprofs;
        for (ii = 0; ii < numparts; ii++) {     // parts
            for (jj = 0; jj < numsubbands; jj++) {      // subbands
                statindex = ii * numsubbands + jj;
                offset = random() % proflen;
                phsindex = 0;
                for (kk = offset; kk < proflen; kk++, phsindex++)       // phases
                    sumprof[phsindex] += prof_ptr[kk];
                for (kk = 0; kk < offset; kk++, phsindex++)     // phases
                    sumprof[phsindex] += prof_ptr[kk];
                prof_ptr += proflen;
                prof_avg += stats[statindex].prof_avg;
                prof_var += stats[statindex].prof_var;
            }
        }
        /* Calculate the current chi-squared */
        redchi = chisqr(sumprof, proflen, prof_avg, prof_var) / dofeff;
        chis[trialnum] = (float) redchi;
    }
    avg_var(chis, numtrials, &chi_avg, &chi_var);
    vect_free(chis);
    vect_free(sumprof);
    return 1.0 / chi_avg;
}
