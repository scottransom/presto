#include "presto.h"
#include "cpgplot.h"
#ifdef USEMMAP
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define DEBUGOUT 0

/* Note:  zoomlevel is simply (LOGMAXDISPNUM-Log_2(numsamps)) */
#define LOGMAXDISPNUM      13   /* 8192: Maximum number of points to display at once */
#define LOGMINDISPNUM      3    /* 8: Minimum number of points to display at once */
#define LOGMINCHUNKLEN     3    /* 8: The minimum number of real points in a stats chunk */
#define LOGMAXPTS          28   /* 256M points */
#define LOGINITIALNUMPTS   17   /* 131072: The initial number of samples to plot */
#define MAXDISPNUM      (1<<LOGMAXDISPNUM)
#define MINDISPNUM      (1<<LOGMINDISPNUM)
#define MINCHUNKLEN     (1<<LOGMINCHUNKLEN)
#define MAXPTS          (1<<LOGMAXPTS)
#define INITIALNUMPTS   (1<<LOGINITIALNUMPTS)
#define AVGMED_COLOR 11
#define STDDEV_COLOR 10
#define LARGENUM 1.0e32
#define SMALLNUM -1.0e32

static long long Ndat;
static infodata idata;
#ifdef USEMMAP
static int mmap_file;
#else
static FILE *datfile;
#endif
static int plotstats = 0, usemedian = 0;
/* plotstats: 0 = both, 1 = stats only, 2 = data only */
/* usemedian: 0 = average, 1 = median */

typedef struct datapart {
    double tlo;                 /* Elapsed time (s) of the first point */
    long nlo;                   /* The sample number of the first point */
    long nn;                    /* The total number samples in *data */
    float *data;                 /* Raw data  */
} datapart;

typedef struct dataview {
    double vdt;                 /* Data view time step (2.0**(-zoomlevel))*dt */
    float minval;                /* The minimum sample value in this view */
    float maxval;                /* The maximum sample value in this view */
    long centern;               /* The center sample to plot */
    long lon;                   /* The lowest sample to plot */
    long zoomlevel;             /* Positive = zoomed in, Negative = zoomed out */
    long numsamps;              /* The number of samples covered by the display */
    long dispnum;               /* The number of points actually plotted */
    long chunklen;              /* The length of the chunk of samples used to calculate stats */
    long numchunks;             /* The number of chunks that are being displayed */
    float avgmeds[MAXDISPNUM];   /* The average or median samples for each chunk */
    float stds[MAXDISPNUM];      /* The standard deviation of the samples for each chunk */
    float maxs[MAXDISPNUM];      /* The maximum samples for each chunk */
    float mins[MAXDISPNUM];      /* The minimum samples for each chunk */
    float vals[MAXDISPNUM];      /* The raw data values when zoomlevel > 0 */
} dataview;


typedef struct basicstats {
    double average;
    double stdev;
    double skewness;
    double kurtosis;
    double median;
    double min;
    double max;
} basicstats;


static basicstats *calc_stats(dataview * dv, datapart * dp)
{
    int ii, jj;
    float *tmpdata;
    basicstats *tmpstats;

    tmpstats = (basicstats *) malloc(sizeof(basicstats));
    tmpstats->max = SMALLNUM;
    tmpstats->min = LARGENUM;
    tmpdata = gen_fvect(dv->numsamps);
    for (ii = 0, jj = dv->lon - dp->nlo; ii < dv->numsamps; ii++, jj++) {
        tmpdata[ii] = dp->data[jj];
        if (tmpdata[ii] > tmpstats->max)
            tmpstats->max = tmpdata[ii];
        if (tmpdata[ii] < tmpstats->min)
            tmpstats->min = tmpdata[ii];
    }
    stats(tmpdata, dv->numsamps, &tmpstats->average, &tmpstats->stdev,
          &tmpstats->skewness, &tmpstats->kurtosis);
    tmpstats->stdev = sqrt(tmpstats->stdev);
    tmpstats->median = median(tmpdata, dv->numsamps);
    vect_free(tmpdata);
    return tmpstats;
}

static int plot_dataview(dataview * dv, float minval, float maxval, float charhgt)
/* The return value is offsetn */
{
    long ii, lon, hin, offsetn = 0, tmpn;
    double lot, hit, offsett = 0.0;
    float ns[MAXDISPNUM], hiavg[MAXDISPNUM], loavg[MAXDISPNUM];
    float scalemin = 0.0, scalemax = 0.0, dscale;

    cpgsave();
    cpgbbuf();

    /* Set the "Normal" plotting attributes */

    cpgsls(1);
    cpgslw(1);
    cpgsch(charhgt);
    cpgsci(1);
    cpgvstd();

    /* Autoscale for the maximum value */
    if (maxval > 0.5 * LARGENUM)
        scalemax = dv->maxval;
    else
        scalemax = maxval;
    /* Autoscale for the minimum value */
    if (minval < 0.5 * SMALLNUM)
        scalemin = dv->minval;
    else
        scalemin = minval;
    dscale = 0.1 * (scalemax - scalemin);
    if (maxval > 0.5 * LARGENUM)
        maxval = scalemax + dscale;
    if (minval < 0.5 * SMALLNUM)
        minval = scalemin - dscale;

    lon = dv->lon;
    lot = lon * idata.dt;
    hin = lon + dv->numsamps;
    hit = hin * idata.dt;

    /* Time Labels (top of box) */

    if ((hit - lot) / hit < 0.001) {
        int numchar;
        char label[50];

        offsett = 0.5 * (hit + lot);
        numchar = snprintf(label, 50, "Time - %.15g (s)", offsett);
        cpgmtxt("T", 2.5, 0.5, 0.5, label);
    } else {
        cpgmtxt("T", 2.5, 0.5, 0.5, "Time (s)");
    }
    cpgswin(lot - offsett, hit - offsett, minval, maxval);
    cpgbox("CMST", 0.0, 0, "", 0.0, 0);

    /* Sample number labels */

    if (lon > 10000000 || (double) (hin - lon) / (double) hin < 0.001) {
        int numchar;
        char label[50];

        offsetn = (lon / 10000) * 10000;
        numchar = snprintf(label, 50, "Sample - %ld", offsetn);
        cpgmtxt("B", 2.8, 0.5, 0.5, label);
    } else {
        cpgmtxt("B", 2.8, 0.5, 0.5, "Sample");
    }
    cpgswin(lon - offsetn, hin - offsetn, minval, maxval);
    cpgbox("BNST", 0.0, 0, "BCNST", 0.0, 0);

    /* Plot the rawdata if required */

    tmpn = lon - offsetn;
    if (plotstats == 0 || plotstats == 2) {
        if (dv->zoomlevel > 0) {
            for (ii = 0; ii < dv->dispnum; ii++)
                ns[ii] = tmpn + ii;
            cpgbin(dv->dispnum, ns, dv->vals, 0);
        } else {                /* Plot the min/max values */
            for (ii = 0; ii < dv->numchunks; ii++, tmpn += dv->chunklen) {
                cpgmove((float) tmpn, dv->mins[ii]);
                cpgdraw((float) tmpn, dv->maxs[ii]);
            }
        }
    }

    /* Plot the other statistics if requested */

    if (plotstats == 0 || plotstats == 1) {
        tmpn = lon - offsetn;
        for (ii = 0; ii < dv->numchunks; ii++, tmpn += dv->chunklen) {
            ns[ii] = tmpn;
            hiavg[ii] = dv->avgmeds[ii] + dv->stds[ii];
            loavg[ii] = dv->avgmeds[ii] - dv->stds[ii];
        }
        if (dv->numchunks > 512) {
            if (plotstats == 1) {
                cpgline(dv->numchunks, ns, dv->mins);
                cpgline(dv->numchunks, ns, dv->maxs);
            }
            cpgsci(AVGMED_COLOR);
            cpgline(dv->numchunks, ns, dv->avgmeds);
            if (usemedian)
                cpgmtxt("T", -1.4, 0.02, 0.0, "Median");
            else
                cpgmtxt("T", -1.4, 0.02, 0.0, "Average");
            cpgsci(STDDEV_COLOR);
            cpgline(dv->numchunks, ns, hiavg);
            cpgline(dv->numchunks, ns, loavg);
            cpgmtxt("T", -2.6, 0.02, 0.0, "+/- 1 Std Dev");
        } else {
            if (plotstats == 1) {
                cpgbin(dv->numchunks, ns, dv->mins, 0);
                cpgbin(dv->numchunks, ns, dv->maxs, 0);
            }
            cpgsci(AVGMED_COLOR);
            cpgbin(dv->numchunks, ns, dv->avgmeds, 0);
            if (usemedian)
                cpgmtxt("T", -1.4, 0.02, 0.0, "Median");
            else
                cpgmtxt("T", -1.4, 0.02, 0.0, "Average");
            cpgsci(STDDEV_COLOR);
            cpgbin(dv->numchunks, ns, hiavg, 0);
            cpgbin(dv->numchunks, ns, loavg, 0);
            cpgmtxt("T", -2.6, 0.02, 0.0, "+/- 1 Std Dev");
        }
    }
    cpgsci(1);
    cpgmtxt("L", 2.5, 0.5, 0.5, "Sample Value");
    cpgebuf();
    cpgunsa();
    return offsetn;
}


static dataview *get_dataview(long centern, long zoomlevel, datapart * dp)
{
    long ii, jj, offset;
    double tmpavg, tmpvar;
    float *tmpchunk;
    dataview *dv;

    dv = (dataview *) malloc(sizeof(dataview));
    dv->zoomlevel = zoomlevel;
    dv->numsamps = (1 << (LOGMAXDISPNUM - zoomlevel));
    if (dv->numsamps > dp->nn)
        dv->numsamps = next2_to_n(dp->nn) / 2;
    dv->chunklen = (zoomlevel < -LOGMINCHUNKLEN) ?
        (1 << abs(zoomlevel)) : (1 << LOGMINCHUNKLEN);
    dv->dispnum = (dv->numsamps > MAXDISPNUM) ? MAXDISPNUM : dv->numsamps;
    if (DEBUGOUT)
        printf("zoomlevel = %ld  numsamps = %ld  chunklen = %ld  dispnum %ld  nn = %ld\n",
               dv->zoomlevel, dv->numsamps, dv->chunklen, dv->dispnum, dp->nn);
    dv->numchunks = dv->numsamps / dv->chunklen;
    dv->centern = centern;
    dv->lon = centern - dv->numsamps / 2;
    dv->vdt = dv->chunklen * idata.dt;
    dv->maxval = SMALLNUM;
    dv->minval = LARGENUM;
    if (dv->lon < 0) {
        dv->lon = 0;
        dv->centern = dv->lon + dv->numsamps / 2;
    }
    if (dv->lon + dv->numsamps >= dp->nn) {
        dv->lon = dp->nn - dv->numsamps;
        dv->centern = dv->lon + dv->numsamps / 2;
    }
    tmpchunk = gen_fvect(dv->chunklen);
    for (ii = 0; ii < dv->numchunks; ii++) {
        float tmpmin = LARGENUM, tmpmax = SMALLNUM, tmpval;
        offset = dv->lon + ii * dv->chunklen;
        memcpy(tmpchunk, dp->data + offset, sizeof(float) * dv->chunklen);
        avg_var(dp->data + offset, dv->chunklen, &tmpavg, &tmpvar);
        if (usemedian)
            dv->avgmeds[ii] = median(tmpchunk, dv->chunklen);
        else
            dv->avgmeds[ii] = tmpavg;
        dv->stds[ii] = sqrt(tmpvar);
        for (jj = 0; jj < dv->chunklen; jj++, offset++) {
            tmpval = dp->data[offset];
            if (tmpval > tmpmax)
                tmpmax = tmpval;
            if (tmpval < tmpmin)
                tmpmin = tmpval;
        }
        dv->maxs[ii] = tmpmax;
        if (tmpmax > dv->maxval)
            dv->maxval = tmpmax;
        dv->mins[ii] = tmpmin;
        if (tmpmin < dv->minval)
            dv->minval = tmpmin;
    }
    vect_free(tmpchunk);
    offset = dv->lon;
    if (zoomlevel > 0) {
        for (ii = 0, offset = dv->lon; ii < dv->numsamps; ii++, offset++)
            *(dv->vals + ii) = *(dp->data + offset);
    }
    return dv;
}


static datapart *get_datapart(long nlo, long numn)
{
    datapart *dp;

    if (nlo + numn > Ndat)
        return NULL;
    else {
        dp = (datapart *) malloc(sizeof(datapart));
        dp->nn = numn;
        dp->nlo = nlo;
        dp->tlo = idata.dt * nlo;
#ifdef USEMMAP
        dp->data = (float *) mmap(0, sizeof(float) * numn, PROT_READ,
                                  MAP_SHARED, mmap_file, 0);
#else
        dp->data = read_float_file(datfile, nlo, numn);
#endif
        return dp;
    }
}


static void free_datapart(datapart * dp)
{
#ifdef USEMMAP
    munmap(dp->data, sizeof(float) * dp->nn);
#else
    vect_free(dp->data);
#endif
    free(dp);
}


static void print_help(void)
{
    printf("\n"
           " Button or Key            Effect\n"
           " -------------            ------\n"
           " Left Mouse or I or A     Zoom in  by a factor of 2\n"
           " Right Mouse or O or X    Zoom out by a factor of 2\n"
           " <                        Shift left  by a full screen width\n"
           " >                        Shift right by a full screen width\n"
           " ,                        Shift left  by 1/8 of the screen width\n"
           " .                        Shift right by 1/8 of the screen width\n"
           " +/_                      Increase/Decrease the top edge\n"
           " =/-                      Increase/Decrease the bottom edge\n"
           " SPACE                    Toggle statistics and sample plotting on/off\n"
           " M                        Toggle between median and average\n"
           " S                        Scale the powers automatically\n"
           " V                        Print the statistics for the current view\n"
           " P                        Print the current plot to a file\n"
           " G                        Go to a specified time\n"
           " ?                        Show this help screen\n"
           " Q                        Quit\n" "\n");
}


int main(int argc, char *argv[])
{
    float minval = SMALLNUM, maxval = LARGENUM, inx = 0, iny = 0;
    long centern, offsetn;
    long zoomlevel, maxzoom = 0, minzoom, xid, psid;
    char *rootfilenm, inchar;
    datapart *lodp;
    dataview *dv;
    basicstats *statvals;

    if (argc == 1) {
        printf("\nusage:  exploredat datafilename\n\n");
        exit(0);
    }

    printf("\n\n");
    printf("      Interactive Data Explorer\n");
    printf("         by Scott M. Ransom\n");
    printf("            November, 2001\n");
    print_help();

    {
        int hassuffix = 0;
        char *suffix;

        hassuffix = split_root_suffix(argv[1], &rootfilenm, &suffix);
        if (hassuffix) {
            if (strcmp(suffix, "dat") != 0) {
                printf
                    ("\nInput file ('%s') must be a single PRESTO data file ('.dat')!\n\n",
                     argv[1]);
                free(suffix);
                exit(0);
            }
            free(suffix);
        } else {
            printf("\nInput file ('%s') must be a PRESTO data file ('.dat')!\n\n",
                   argv[1]);
            exit(0);
        }
    }

    /* Read the info file */

    readinf(&idata, rootfilenm);
    if (strlen(remove_whitespace(idata.object)) > 0) {
        printf("Examining %s data from '%s'.\n\n",
               remove_whitespace(idata.object), argv[1]);
    } else {
        printf("Examining data from '%s'.\n\n", argv[1]);
    }
#ifdef USEMMAP
    mmap_file = open(argv[1], O_RDONLY);
    {
        int rt;
        struct stat buf;

        rt = fstat(mmap_file, &buf);
        if (rt == -1) {
            perror("\nError in fstat() in exploredat.c");
            printf("\n");
            exit(-1);
        }
        Ndat = buf.st_size / sizeof(float);
    }
    lodp = get_datapart(0, Ndat);
#else
    {
        long numsamp;

        datfile = chkfopen(argv[1], "rb");
        Ndat = chkfilelen(datfile, sizeof(float));
        numsamp = (Ndat > MAXPTS) ? (long) MAXPTS : (long) Ndat;
        lodp = get_datapart(0, numsamp);
    }
#endif

    /* Plot the initial data */

    centern = 0.5 * INITIALNUMPTS;
    if (centern > lodp->nn)
        centern = lodp->nn / 2;
    zoomlevel = LOGMAXDISPNUM - LOGINITIALNUMPTS;
    minzoom = LOGMAXDISPNUM - LOGMAXPTS;
    maxzoom = LOGMAXDISPNUM - LOGMINDISPNUM;
    dv = get_dataview(centern, zoomlevel, lodp);

    /* Prep the XWIN device for PGPLOT */

    xid = cpgopen("/XWIN");
    if (xid <= 0) {
        free_datapart(lodp);
#ifdef USEMMAP
        close(mmap_file);
#else
        fclose(datfile);
#endif
        free(dv);
        exit(EXIT_FAILURE);
    }
    cpgask(0);
    cpgpage();
    offsetn = plot_dataview(dv, minval, maxval, 1.0);

    do {
        cpgcurs(&inx, &iny, &inchar);
        if (DEBUGOUT)
            printf("You pressed '%c'\n", inchar);

        switch (inchar) {
        case ' ':              /* Toggle stats and sample plotting on/off */
            /* 0 = both, 1 = stats only, 2 = data only */
            plotstats++;
            plotstats = plotstats % 3;
            cpgpage();
            offsetn = plot_dataview(dv, minval, maxval, 1.0);
            break;
        case 'M':              /* Toggle between median and average */
            /* FALLTHRU */
        case 'm':
            usemedian = (usemedian) ? 0 : 1;
            free(dv);
            dv = get_dataview(centern, zoomlevel, lodp);
            cpgpage();
            offsetn = plot_dataview(dv, minval, maxval, 1.0);
            break;
        case 'A':              /* Zoom in */
            /* FALLTHRU */
        case 'a':
            centern = inx + offsetn;
            /* FALLTHRU */
        case 'I':
            /* FALLTHRU */
        case 'i':
            if (DEBUGOUT)
                printf("  Zooming in  (zoomlevel = %ld)...\n", zoomlevel);
            if (zoomlevel < maxzoom) {
                zoomlevel++;
                free(dv);
                dv = get_dataview(centern, zoomlevel, lodp);
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
            } else
                printf("  Already at maximum zoom level (%ld).\n", zoomlevel);
            break;
        case 'X':              /* Zoom out */
            /* FALLTHRU */
        case 'x':
            /* FALLTHRU */
        case 'O':
            /* FALLTHRU */
        case 'o':
            if (DEBUGOUT)
                printf("  Zooming out  (zoomlevel = %ld)...\n", zoomlevel);
            if (zoomlevel > minzoom) {
                zoomlevel--;
                free(dv);
                dv = get_dataview(centern, zoomlevel, lodp);
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
            } else
                printf("  Already at minimum zoom level (%ld).\n", zoomlevel);
            break;
        case '<':              /* Shift left 1 full screen */
            centern -= dv->numsamps + dv->numsamps / 8;
            /* FALLTHRU */
        case ',':              /* Shift left 1/8 screen */
            if (DEBUGOUT)
                printf("  Shifting left...\n");
            centern -= dv->numsamps / 8;
            {                   /* Should probably get the previous chunk from the datfile... */
                double lowestr;

                lowestr = 0.5 * dv->numsamps;
                if (centern < lowestr)
                    centern = lowestr;
            }
            free(dv);
            dv = get_dataview(centern, zoomlevel, lodp);
            cpgpage();
            offsetn = plot_dataview(dv, minval, maxval, 1.0);
            break;
        case '>':              /* Shift right 1 full screen */
            centern += dv->numsamps - dv->numsamps / 8;
            /* FALLTHRU */
        case '.':              /* Shift right 1/8 screen */
            centern += dv->numsamps / 8;
            if (DEBUGOUT)
                printf("  Shifting right...\n");
            {                   /* Should probably get the next chunk from the datfile... */
                double highestr;

                highestr = lodp->nlo + lodp->nn - 0.5 * dv->numsamps;
                if (centern > highestr)
                    centern = highestr;
            }
            free(dv);
            dv = get_dataview(centern, zoomlevel, lodp);
            cpgpage();
            offsetn = plot_dataview(dv, minval, maxval, 1.0);
            break;
        case '+':              /* Increase height of top edge */
            {
                float dy;

                if (maxval > 0.5 * LARGENUM) {
                    printf("  Auto-scaling of top edge is off.\n");
                    if (minval < 0.5 * SMALLNUM)
                        dy = dv->maxval - dv->minval;
                    else
                        dy = dv->maxval - minval;
                    maxval = dv->maxval + 0.1 * dy;
                } else {
                    if (minval < 0.5 * SMALLNUM)
                        dy = maxval - dv->minval;
                    else
                        dy = maxval - minval;
                    maxval += 0.1 * dy;
                }
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
                break;
            }
        case '_':              /* Decrease height of top edge */
            {
                float dy;

                if (maxval > 0.5 * LARGENUM) {
                    printf("  Auto-scaling of top edge is off.\n");
                    if (minval < 0.5 * SMALLNUM)
                        dy = dv->maxval - dv->minval;
                    else
                        dy = dv->maxval - minval;
                    maxval = dv->maxval - 0.1 * dy;
                } else {
                    if (minval < 0.5 * SMALLNUM)
                        dy = maxval - dv->minval;
                    else
                        dy = maxval - minval;
                    maxval -= 0.1 * dy;
                }
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
                break;
            }
        case '=':              /* Increase height of bottom edge */
            {
                float dy;

                if (minval < 0.5 * SMALLNUM) {
                    printf("  Auto-scaling of bottom edge is off.\n");
                    if (maxval > 0.5 * LARGENUM)
                        dy = dv->maxval - dv->minval;
                    else
                        dy = maxval - dv->minval;
                    minval = dv->minval + 0.1 * dy;
                } else {
                    if (maxval > 0.5 * LARGENUM)
                        dy = dv->maxval - minval;
                    else
                        dy = maxval - minval;
                    minval += 0.1 * dy;
                }
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
                break;
            }
        case '-':              /* Decrease height of bottom edge */
            {
                float dy;

                if (minval < 0.5 * SMALLNUM) {
                    printf("  Auto-scaling of bottom edge is off.\n");
                    if (maxval > 0.5 * LARGENUM)
                        dy = dv->maxval - dv->minval;
                    else
                        dy = maxval - dv->minval;
                    minval = dv->minval - 0.1 * dy;
                } else {
                    if (maxval > 0.5 * LARGENUM)
                        dy = dv->maxval - minval;
                    else
                        dy = maxval - minval;
                    minval -= 0.1 * dy;
                }
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
                break;
            }
        case 'S':              /* Auto-scale */
            /* FALLTHRU */
        case 's':
            printf("  Auto-scaling is on.\n");
            minval = SMALLNUM;
            maxval = LARGENUM;
            cpgpage();
            offsetn = plot_dataview(dv, minval, maxval, 1.0);
            break;
        case 'G':              /* Goto a time */
            /* FALLTHRU */
        case 'g':
            {
                char timestr[50];
                double time = -1.0;

                while (time < 0.0) {
                    printf
                        ("  Enter the time (s) from the beginning of the file to go to:\n");
                    fgets(timestr, 50, stdin);
                    timestr[strlen(timestr) - 1] = '\0';
                    time = atof(timestr);
                }
                offsetn = 0.0;
                centern = (long) (time / idata.dt + 0.5);
                printf("  Moving to time %.15g (data point %ld).\n", time, centern);
                free(dv);
                dv = get_dataview(centern, zoomlevel, lodp);
                cpgpage();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
            }
            break;
        case '?':              /* Print help screen */
            print_help();
            break;
        case 'P':              /* Print the current plot */
            /* FALLTHRU */
        case 'p':
            {
                int len;
                char filename[200];

                printf("  Enter the filename to save the plot as:\n");
                fgets(filename, 195, stdin);
                len = strlen(filename) - 1;
                filename[len + 0] = '/';
                filename[len + 1] = 'C';
                filename[len + 2] = 'P';
                filename[len + 3] = 'S';
                filename[len + 4] = '\0';
                psid = cpgopen(filename);
                cpgslct(psid);
                cpgpap(10.25, 8.5 / 11.0);
                cpgiden();
                offsetn = plot_dataview(dv, minval, maxval, 1.0);
                cpgclos();
                cpgslct(xid);
                filename[len] = '\0';
                printf("  Wrote the plot to the file '%s'.\n", filename);
            }
            break;
        case 'V':              /* Show the basic statistics for the current dataview */
            /* FALLTHRU */
        case 'v':
            statvals = calc_stats(dv, lodp);
            printf("\n  Statistics:\n"
                   "    Low sample               %ld\n"
                   "    Number of samples        %ld\n"
                   "    Low time (s)             %.7g\n"
                   "    Duration of samples (s)  %.7g\n"
                   "    Maximum value            %.7g\n"
                   "    Minimum value            %.7g\n"
                   "    Average value            %.7g\n"
                   "    Median value             %.7g\n"
                   "    Standard Deviation       %.7g\n"
                   "    Skewness                 %.7g\n"
                   "    Kurtosis                 %.7g\n\n",
                   dv->lon, dv->numsamps, dv->lon * idata.dt,
                   dv->numsamps * idata.dt, statvals->max, statvals->min,
                   statvals->average, statvals->median, statvals->stdev,
                   statvals->skewness, statvals->kurtosis);
            free(statvals);
            break;
        case 'Q':              /* Quit */
            /* FALLTHRU */
        case 'q':
            printf("  Quitting...\n");
            free(dv);
            cpgclos();
            break;
        default:
            printf("  Unrecognized option '%c'.\n", inchar);
            break;
        }
    } while (inchar != 'Q' && inchar != 'q');

    free_datapart(lodp);
#ifdef USEMMAP
    close(mmap_file);
#else
    fclose(datfile);
#endif
    printf("Done\n\n");
    return 0;
}
