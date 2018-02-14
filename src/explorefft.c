#include "presto.h"
#include "cpgplot.h"

/*#undef USEMMAP*/

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

int compare_birds(const void *ca, const void *cb);

/* Note:  zoomlevel is simply (LOGDISPLAYNUM-Log_2(numbins)) */
#define LOGDISPLAYNUM     13    /* 1024: Maximum number of points to display at once */
#define LOGLOCALCHUNK     4     /* 16: Chunk size for polynomial fit */
#define LOGMINBINS        5     /* 32 points */
#define LOGMAXBINS        27    /* 128M points */
#define LOGINITIALNUMBINS 17    /* 131072: The initial number of bins to plot */
#define DISPLAYNUM     (1<<LOGDISPLAYNUM)
#define LOCALCHUNK     (1<<LOGLOCALCHUNK)
#define MINBINS        (1<<LOGMINBINS)
#define MAXBINS        (1<<LOGMAXBINS)
#define INITIALNUMBINS (1<<LOGINITIALNUMBINS)

static long long N;             /* Number of points in the original time series */
static long long Nfft;          /* Number of bins in the FFT */
static double T;                /* The time duration of FFT */
static float r0;                /* The value of the zeroth Fourier freq */
static infodata idata;
#ifdef USEMMAP
static int mmap_file;
#else
static FILE *fftfile;
#endif
static double norm_const = 0.0; /* Used if the user specifies a power normalization. */
static int numzaplist = 0;      /* The number of actual lobin/hibin pairs in zaplist */
static int lenzaplist = 0;      /* The number of possible lobin/hibin pairs in zaplist */
static bird *zaplist = NULL;

typedef struct fftpart {
    long rlo;                  /* Lowest Fourier freq displayed */
    long numamps;              /* Number of raw amplitudes */
    float maxrawpow;            /* The highest raw power present */
    float *rawpowers;           /* The raw powers */
    float *medians;             /* The local median values (chunks of size LOCALCHUNK bins) */
    float *normvals;            /* The values to use for normalization (default is median/-log(0.5)) */
    fcomplex *amps;            /* Raw FFT amplitudes    */
} fftpart;

typedef struct fftview {
    double dr;                  /* Fourier frequency stepsize (2.0**(-zoomlevel)) */
    double centerr;             /* The center Fourier freq to plot */
    long lor;                   /* The lowest Fourier freq to plot */
    int zoomlevel;              /* Positive = zoomed in, Negative = zoomed out */
    int numbins;                /* The number of full bins from low to high to display */
    float maxpow;               /* The maximum normalized power in the view */
    float powers[DISPLAYNUM];   /* The normalized powers to display */
    double rs[DISPLAYNUM];      /* The corresponding Fourier freqs */
} fftview;


/*
static int floor_log2(int nn)
{
  int ii=0;
  
  if (nn <= 0)
    return -1;
  while (nn >> 1){
    ii++;
    nn >>= 1;
  }
  return ii;
}
*/

static double plot_fftview(fftview * fv, float maxpow, float charhgt,
                           float vertline, int vertline_color)
/* The return value is offsetf */
{
    int ii;
    double lor, lof, hir, hif, offsetf = 0.0;
    float *freqs;

    cpgsave();
    cpgbbuf();

    /* Set the "Normal" plotting attributes */

    cpgsls(1);
    cpgslw(1);
    cpgsch(charhgt);
    cpgsci(1);
    cpgvstd();

    if (maxpow == 0.0)          /* Autoscale for the maximum value */
        maxpow = 1.1 * fv->maxpow;

    lor = fv->lor;
    lof = lor / T;
    hir = lor + fv->dr * DISPLAYNUM;
    hif = hir / T;
    offsetf = 0.0;

    /* Period Labels */

    if (fv->zoomlevel >= 0 && lof > 1.0) {
        double lop, hip, offsetp = 0.0;
        lop = 1.0 / lof;
        hip = 1.0 / hif;
        offsetp = 0.0;

        if ((lop - hip) / hip < 0.001) {
            int numchar;
            char label[50];

            offsetp = 0.5 * (hip + lop);
            numchar = snprintf(label, 50, "Period - %.15g (s)", offsetp);
            cpgmtxt("T", 2.5, 0.5, 0.5, label);
        } else {
            cpgmtxt("T", 2.5, 0.5, 0.5, "Period (s)");
        }
        cpgswin(lop - offsetp, hip - offsetp, 0.0, maxpow);
        cpgbox("CIMST", 0.0, 0, "", 0.0, 0);
    }

    /* Frequency Labels */

    if ((hif - lof) / hif < 0.001) {
        int numchar;
        char label[50];

        offsetf = 0.5 * (hif + lof);
        numchar = snprintf(label, 50, "Frequency - %.15g (Hz)", offsetf);
        cpgmtxt("B", 2.8, 0.5, 0.5, label);
    } else {
        cpgmtxt("B", 2.8, 0.5, 0.5, "Frequency (Hz)");
    }
    cpgswin(lof - offsetf, hif - offsetf, 0.0, maxpow);

    /* Add zapboxes if required */

    if (numzaplist) {
        double zaplo, zaphi;

        cpgsave();
        cpgsci(15);
        cpgsfs(1);
        for (ii = 0; ii < numzaplist; ii++) {
            zaplo = zaplist[ii].lobin;
            zaphi = zaplist[ii].hibin;
            if ((zaplo < hir && zaplo > lor) || (zaphi < hir && zaphi > lor)) {
                cpgrect(zaplo / T - offsetf, zaphi / T - offsetf, 0.0,
                        0.95 * maxpow);
            }
        }
        cpgunsa();
    }

    /* Add a background vertical line if requested */

    if (vertline != 0.0 && vertline_color != 0) {
        cpgsave();
        cpgsci(vertline_color);
        cpgmove(vertline / T - offsetf, 0.0);
        cpgdraw(vertline / T - offsetf, maxpow);
        cpgunsa();
    }

    if (fv->zoomlevel >= 0 && lof > 1.0)
        cpgbox("BINST", 0.0, 0, "BCNST", 0.0, 0);
    else
        cpgbox("BCINST", 0.0, 0, "BCNST", 0.0, 0);

    /* Plot the spectrum */

    freqs = gen_fvect(DISPLAYNUM);
    for (ii = 0; ii < DISPLAYNUM; ii++)
        freqs[ii] = fv->rs[ii] / T - offsetf;
    if (fv->zoomlevel > 0) {    /* Magnified power spectrum */
        cpgline(DISPLAYNUM, freqs, fv->powers);
    } else {                    /* Down-sampled power spectrum */
        for (ii = 0; ii < DISPLAYNUM; ii++) {
            cpgmove(freqs[ii], 0.0);
            cpgdraw(freqs[ii], fv->powers[ii]);
        }
    }
    vect_free(freqs);
    cpgmtxt("L", 2.5, 0.5, 0.5, "Normalized Power");
    cpgebuf();
    cpgunsa();
    return offsetf;
}


static fftview *get_fftview(double centerr, int zoomlevel, fftpart * fp)
{
    int ii;
    double split = 1.0 / (double) LOCALCHUNK;
    float powargr, powargi;
    fftview *fv;

    fv = (fftview *) malloc(sizeof(fftview));
    fv->zoomlevel = zoomlevel;
    fv->centerr = centerr;
    if (zoomlevel > 0) {        /* Magnified power spectrum */
        int numbetween, nextbin, index;
        fcomplex *interp;

        numbetween = (1 << zoomlevel);
        fv->numbins = DISPLAYNUM / numbetween;
        fv->dr = 1.0 / (double) numbetween;
        fv->lor = centerr - fv->numbins / 2;
        if (fv->lor + fv->numbins > Nfft) {
            fv->lor = Nfft - fv->numbins;
            fv->centerr = fv->lor + fv->numbins / 2;
        }
        if (fv->lor < 0)
            fv->lor = 0;
        interp = corr_rz_interp(fp->amps, fp->numamps, numbetween,
                                fv->lor - fp->rlo, 0.0, DISPLAYNUM * 2,
                                LOWACC, &nextbin);
        if (norm_const == 0.0) {
            for (ii = 0; ii < DISPLAYNUM; ii++) {
                fv->rs[ii] = fv->lor + ii * fv->dr;
                index = (long) ((fv->rs[ii] - fp->rlo) * split + 0.5);
                fv->powers[ii] =
                    POWER(interp[ii].r, interp[ii].i) * fp->normvals[index];
            }
        } else {
            for (ii = 0; ii < DISPLAYNUM; ii++) {
                fv->rs[ii] = fv->lor + ii * fv->dr;
                fv->powers[ii] = POWER(interp[ii].r, interp[ii].i) * norm_const;
            }
        }
        vect_free(interp);
    } else {                    /* Down-sampled power spectrum */
        long jj, powindex, normindex, binstocombine;
        float *tmprawpwrs, maxpow;

        binstocombine = (1 << abs(zoomlevel));
        fv->numbins = DISPLAYNUM * binstocombine;
        fv->dr = (double) binstocombine;
        fv->lor = (long) floor(centerr - 0.5 * fv->numbins);
        if (fv->lor + fv->numbins > Nfft) {
            fv->lor = Nfft - fv->numbins;
            fv->centerr = fv->lor + fv->numbins / 2;
        }
        if (fv->lor < 0)
            fv->lor = 0;
        tmprawpwrs = gen_fvect(fv->numbins);
        if (norm_const == 0.0) {
            for (ii = 0; ii < fv->numbins; ii++) {
                powindex = (long) (fv->lor - fp->rlo + ii + 0.5);
                normindex = (long) (powindex * split);
                tmprawpwrs[ii] = fp->rawpowers[powindex] * fp->normvals[normindex];
            }
        } else {
            for (ii = 0; ii < fv->numbins; ii++) {
                powindex = (long) (fv->lor - fp->rlo + ii + 0.5);
                tmprawpwrs[ii] = fp->rawpowers[powindex] * norm_const;
            }
        }
        powindex = 0;
        for (ii = 0; ii < DISPLAYNUM; ii++) {
            maxpow = 0.0;
            for (jj = 0; jj < binstocombine; jj++, powindex++)
                if (tmprawpwrs[powindex] > maxpow)
                    maxpow = tmprawpwrs[powindex];
            fv->rs[ii] = fv->lor + ii * fv->dr;
            fv->powers[ii] = maxpow;
        }
    }
    fv->maxpow = 0.0;
    for (ii = 0; ii < DISPLAYNUM; ii++)
        if (fv->powers[ii] > fv->maxpow)
            fv->maxpow = fv->powers[ii];
    return fv;
}


static fftpart *get_fftpart(long rlo, long numr)
{
    int ii, jj, index;
    float powargr, powargi, tmppwr, chunk[LOCALCHUNK];
    fftpart *fp;

    if (rlo + numr > Nfft) {
        return NULL;
    } else {
        fp = (fftpart *) malloc(sizeof(fftpart));
        fp->rlo = rlo;
        fp->numamps = numr;
#ifdef USEMMAP
        fp->amps = (fcomplex *) mmap(0, sizeof(fcomplex) * numr, PROT_READ,
                                     MAP_SHARED, mmap_file, 0);
#else
        fp->amps = read_fcomplex_file(fftfile, rlo, numr);
#endif
        if (rlo == 0)
            r0 = fp->amps[0].r;
        fp->rawpowers = gen_fvect(fp->numamps);
        fp->medians = gen_fvect(fp->numamps / LOCALCHUNK);
        fp->normvals = gen_fvect(fp->numamps / LOCALCHUNK);
        fp->maxrawpow = 0.0;
        for (ii = 0; ii < fp->numamps / LOCALCHUNK; ii++) {
            index = ii * LOCALCHUNK;
            for (jj = 0; jj < LOCALCHUNK; jj++, index++) {
                tmppwr = POWER(fp->amps[index].r, fp->amps[index].i);
                if (tmppwr > fp->maxrawpow)
                    if (rlo || (rlo == 0 && index > 0))
                        fp->maxrawpow = tmppwr;
                fp->rawpowers[index] = tmppwr;
                chunk[jj] = tmppwr;
            }
            if (rlo == 0 && ii == 0) {
                chunk[0] = 1.0;
                fp->rawpowers[0] = 1.0;
            }
            fp->medians[ii] = median(chunk, LOCALCHUNK);
            fp->normvals[ii] = 1.0 / (1.4426950408889634 * fp->medians[ii]);
        }
        return fp;
    }
}

static void free_fftpart(fftpart * fp)
{
    vect_free(fp->normvals);
    vect_free(fp->medians);
    vect_free(fp->rawpowers);
#ifdef USEMMAP
    munmap(fp->amps, sizeof(fcomplex) * fp->numamps);
#else
    vect_free(fp->amps);
#endif
    free(fp);
}


static double find_peak(float inf, fftview * fv, fftpart * fp)
{
    long ii, lobin, hibin, maxbin = 0;
    float maxpow = 0.0;
    double inr, viewfrac = 0.05, newmaxr, newmaxz;
    rderivs derivs;
    fourierprops props;

    inr = inf * T;
    lobin = inr - (fv->numbins * 0.5 * viewfrac);
    hibin = inr + (fv->numbins * 0.5 * viewfrac);
    for (ii = lobin - fp->rlo; ii < hibin - fp->rlo + 1; ii++) {
        if (fp->rawpowers[ii] > maxpow) {
            maxpow = fp->rawpowers[ii];
            maxbin = ii + fp->rlo;
        }
    }
    maxpow = max_rz_arr(fp->amps, fp->numamps, maxbin, 0.0,
                        &newmaxr, &newmaxz, &derivs);
    newmaxr += fp->rlo;
    calc_props(derivs, newmaxr, newmaxz, 0.0, &props);
    print_candidate(&props, idata.dt, N, r0, 2);
    return newmaxr;
}


static void print_help(void)
{
    printf("\n"
           " Button or Key       Effect\n"
           " -------------       ------\n"
           " Left Mouse or I     Zoom in  by a factor of 2\n"
           " Right Mouse or O    Zoom out by a factor of 2\n"
           " Middle Mouse or D   Show details about a selected frequency\n"
           " <                   Shift left  by a full screen width\n"
           " >                   Shift right by a full screen width\n"
           " ,                   Shift left  by 1/8 of the screen width\n"
           " .                   Shift right by 1/8 of the screen width\n"
           " + or =              Increase the power scale (make them taller)\n"
           " - or _              Decrease the power scale (make them shorter)\n"
           " S                   Scale the powers automatically\n"
           " N                   Re-normalize the nowers by one of several methods\n"
           " P                   Print the current plot to a file\n"
           " G                   Go to a specified frequency\n"
           " L                   Load a zaplist showing potential RFI locations\n"
           " Z                   Add a frequency chunk to the RFI zaplist\n"
           " H                   Show the harmonics of the center frequency\n"
           " ?                   Show this help screen\n"
           " Q                   Quit\n" "\n");
}


static fftview *get_harmonic(double rr, int zoomlevel, fftpart * fp)
{
    int numharmbins;

#ifdef USEMMAP
    numharmbins = (1 << (LOGDISPLAYNUM - zoomlevel));
    if (rr + numharmbins > Nfft)
        return NULL;
    else
        return get_fftview(rr, zoomlevel, fp);
#else
    fftpart *harmpart;
    fftview *harmview;

    numharmbins = (1 << (LOGDISPLAYNUM - zoomlevel));
    harmpart = get_fftpart((long) (rr - numharmbins), 2 * numharmbins);
    if (harmpart != NULL) {
        harmview = get_fftview(rr, zoomlevel, harmpart);
        free_fftpart(harmpart);
        return harmview;
    } else {
        return NULL;
    }
#endif
}

static void plot_harmonics(double rr, int zoomlevel, fftpart * fp)
{
    int ii, hh;
    double offsetf;
    char label[20];
    fftview *harmview;

    cpgsubp(4, 4);
    for (ii = 0, hh = 2; ii < 8; ii++, hh++) {
        cpgpanl(ii % 4 + 1, ii / 4 + 1);
        harmview = get_harmonic(hh * rr, zoomlevel, fp);
        if (harmview != NULL) {
            offsetf = plot_fftview(harmview, 0.0, 2.0, hh * rr, 2);
            snprintf(label, 20, "Harmonic %d", hh);
            cpgsave();
            cpgsch(2.0);
            cpgmtxt("T", -1.8, 0.05, 0.0, label);
            cpgunsa();
            free(harmview);
        }
    }
    for (ii = 8, hh = 2; ii < 16; ii++, hh++) {
        cpgpanl(ii % 4 + 1, ii / 4 + 1);
        harmview = get_harmonic(rr / (double) hh, zoomlevel, fp);
        if (harmview != NULL) {
            offsetf = plot_fftview(harmview, 0.0, 2.0, rr / (double) hh, 2);
            snprintf(label, 20, "Harmonic 1/%d", hh);
            cpgsave();
            cpgsch(2.0);
            cpgmtxt("T", -1.8, 0.05, 0.0, label);
            cpgunsa();
            free(harmview);
        }
    }
    cpgsubp(1, 1);
    cpgpanl(1, 1);
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(2.0, 6.0, -2.0, 2.0);
    cpgmove(2.0, 0.0);
    cpgslw(3);
    cpgdraw(6.0, 0.0);
    cpgslw(1);
}

static double harmonic_loop(int xid, double rr, int zoomlevel, fftpart * fp)
{
    float inx = 0.0, iny = 0.0;
    double retval = 0.0;
    int xid2, psid, badchoice = 1;
    char choice;

    xid2 = cpgopen("/XWIN");
    cpgpap(10.25, 8.5 / 11.0);
    cpgask(0);
    cpgslct(xid2);
    plot_harmonics(rr, zoomlevel, fp);
    printf("  Click on the harmonic to go it,\n"
           "    press 'P' to print, or press 'Q' to close.\n");
    while (badchoice) {
        cpgcurs(&inx, &iny, &choice);
        if (choice == 'Q' || choice == 'q') {
            badchoice = 0;
        } else if (choice == 'P' || choice == 'p') {
            int len, numharmbins;
            double offsetf;
            char filename[200];
            fftpart *harmpart;
            fftview *harmview;

            printf("  Enter the filename to save the plot as:\n");
            fgets(filename, 195, stdin);
            len = strlen(filename) - 1;
            strcpy(filename + len, "/CPS");
            psid = cpgopen(filename);
            cpgslct(psid);
            cpgpap(10.25, 8.5 / 11.0);
            cpgiden();
            cpgscr(15, 0.8, 0.8, 0.8);
            numharmbins = (1 << (LOGDISPLAYNUM - zoomlevel));
            harmpart = get_fftpart((long) (rr - numharmbins), 2 * numharmbins);
            harmview = get_fftview(rr, zoomlevel, harmpart);
            free_fftpart(harmpart);
            offsetf = plot_fftview(harmview, 0.0, 1.0, rr, 2);
            cpgpage();
            plot_harmonics(rr, zoomlevel, fp);
            cpgclos();
            cpgslct(xid2);
            cpgscr(15, 0.4, 0.4, 0.4);
            filename[len] = '\0';
            printf("  Wrote the plot to the file '%s'.\n", filename);
        } else if (choice == 'A' || choice == 'a') {
            if (iny > 1.0)
                retval = rr * (long) (inx);
            else if (iny > 0.0)
                retval = rr * ((long) (inx) + 4.0);
            else if (iny > -1.0)
                retval = rr / (long) (inx);
            else
                retval = rr / ((long) (inx) + 4.0);
            badchoice = 0;
        } else {
            printf("  Option not recognized.\n");
        }
    };
    cpgclos();
    cpgslct(xid);
    return retval;
}


int main(int argc, char *argv[])
{
    float maxpow = 0.0, inx = 0.0, iny = 0.0;
    double centerr, offsetf;
    int zoomlevel, maxzoom, minzoom, xid, psid;
    char *rootfilenm, inchar;
    fftpart *lofp;
    fftview *fv;

    if (argc == 1) {
        printf("\nusage:  explorefft fftfilename\n\n");
        exit(0);
    }

    printf("\n\n");
    printf("      Interactive FFT Explorer\n");
    printf("         by Scott M. Ransom\n");
    printf("            October, 2001\n");
    print_help();

    {
        int hassuffix = 0;
        char *suffix;

        hassuffix = split_root_suffix(argv[1], &rootfilenm, &suffix);
        if (hassuffix) {
            if (strcmp(suffix, "fft") != 0) {
                printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
                       argv[1]);
                free(suffix);
                exit(0);
            }
            free(suffix);
        } else {
            printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n", argv[1]);
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
    N = idata.N;
    T = idata.dt * idata.N;
#ifdef USEMMAP
    printf("Memory mapping the input FFT.  This may take a while...\n");
    mmap_file = open(argv[1], O_RDONLY);
    {
        int rt;
        struct stat buf;

        rt = fstat(mmap_file, &buf);
        if (rt == -1) {
            perror("\nError in fstat() in explorefft.c");
            printf("\n");
            exit(-1);
        }
        Nfft = buf.st_size / sizeof(fcomplex);
    }
    lofp = get_fftpart(0, Nfft);
#else
    {
        long numamps;

        fftfile = chkfopen(argv[1], "rb");
        Nfft = chkfilelen(fftfile, sizeof(fcomplex));
        numamps = (Nfft > MAXBINS) ? (long) MAXBINS : (long) Nfft;
        lofp = get_fftpart(0, numamps);
    }
#endif

    /* Plot the initial data */

    {
        long initnumbins = INITIALNUMBINS;

        if (initnumbins > Nfft) {
            initnumbins = next2_to_n(Nfft) / 2;
            zoomlevel = LOGDISPLAYNUM - (int) (log(initnumbins) / log(2.0));
            minzoom = zoomlevel;
        } else {
            zoomlevel = LOGDISPLAYNUM - LOGINITIALNUMBINS;
            minzoom = LOGDISPLAYNUM - LOGMAXBINS;
        }
        maxzoom = LOGDISPLAYNUM - LOGMINBINS;
        centerr = initnumbins / 2;
    }
    fv = get_fftview(centerr, zoomlevel, lofp);

    /* Prep the XWIN device for PGPLOT */

    xid = cpgopen("/XWIN");
    if (xid <= 0) {
        free(fv);
#ifdef USEMMAP
        close(mmap_file);
#else
        fclose(fftfile);
#endif
        free_fftpart(lofp);
        exit(EXIT_FAILURE);
    }
    cpgscr(15, 0.4, 0.4, 0.4);
    cpgask(0);
    cpgpage();
    offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);

    do {
        cpgcurs(&inx, &iny, &inchar);
        if (DEBUGOUT)
            printf("You pressed '%c'\n", inchar);

        switch (inchar) {
        case 'A':              /* Zoom in */
            /* FALLTHRU */
        case 'a':
            centerr = (inx + offsetf) * T;
            /* FALLTHRU */
        case 'I':
            /* FALLTHRU */
        case 'i':
            if (DEBUGOUT)
                printf("  Zooming in  (zoomlevel = %d)...\n", zoomlevel);
            if (zoomlevel < maxzoom) {
                zoomlevel++;
                free(fv);
                fv = get_fftview(centerr, zoomlevel, lofp);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            } else
                printf("  Already at maximum zoom level (%d).\n", zoomlevel);
            break;
        case 'X':              /* Zoom out */
            /* FALLTHRU */
        case 'x':
            /* FALLTHRU */
        case 'O':
            /* FALLTHRU */
        case 'o':
            if (DEBUGOUT)
                printf("  Zooming out  (zoomlevel = %d)...\n", zoomlevel);
            if (zoomlevel > minzoom) {
                zoomlevel--;
                free(fv);
                fv = get_fftview(centerr, zoomlevel, lofp);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            } else
                printf("  Already at minimum zoom level (%d).\n", zoomlevel);
            break;
        case '<':              /* Shift left 1 full screen */
            centerr -= fv->numbins + fv->numbins / 8;
            /* FALLTHRU */
        case ',':              /* Shift left 1/8 screen */
            if (DEBUGOUT)
                printf("  Shifting left...\n");
            centerr -= fv->numbins / 8;
            {                   /* Should probably get the previous chunk from the fftfile... */
                double lowestr;

                lowestr = 0.5 * fv->numbins;
                if (centerr < lowestr)
                    centerr = lowestr;
            }
            free(fv);
            fv = get_fftview(centerr, zoomlevel, lofp);
            cpgpage();
            offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            break;
        case '>':              /* Shift right 1 full screen */
            centerr += fv->numbins - fv->numbins / 8;
            /* FALLTHRU */
        case '.':              /* Shift right 1/8 screen */
            if (DEBUGOUT)
                printf("  Shifting right...\n");
            centerr += fv->numbins / 8;
            {                   /* Should probably get the next chunk from the fftfile... */
                double highestr;

                highestr = lofp->rlo + lofp->numamps - 0.5 * fv->numbins;
                if (centerr > highestr)
                    centerr = highestr;
            }
            free(fv);
            fv = get_fftview(centerr, zoomlevel, lofp);
            cpgpage();
            offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            break;
        case '+':              /* Increase height of powers */
            /* FALLTHRU */
        case '=':
            if (maxpow == 0.0) {
                printf("  Auto-scaling is off.\n");
                maxpow = 1.1 * fv->maxpow;
            }
            maxpow = 3.0 / 4.0 * maxpow;
            cpgpage();
            offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            break;
        case '-':              /* Decrease height of powers */
            /* FALLTHRU */
        case '_':
            if (maxpow == 0.0) {
                printf("  Auto-scaling is off.\n");
                maxpow = 1.1 * fv->maxpow;
            }
            maxpow = 4.0 / 3.0 * maxpow;
            cpgpage();
            offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            break;
        case 'S':              /* Auto-scale */
            /* FALLTHRU */
        case 's':
            if (maxpow == 0.0)
                break;
            else {
                printf("  Auto-scaling is on.\n");
                maxpow = 0.0;
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
                break;
            }
        case 'G':              /* Goto a frequency */
            /* FALLTHRU */
        case 'g':
            {
                char freqstr[50];
                double freq = -1.0;

                while (freq < 0.0) {
                    printf("  Enter the frequency (Hz) to go to:\n");
                    fgets(freqstr, 50, stdin);
                    freqstr[strlen(freqstr) - 1] = '\0';
                    freq = atof(freqstr);
                }
                offsetf = 0.0;
                centerr = freq * T;
                printf("  Moving to frequency %.15g.\n", freq);
                free(fv);
                fv = get_fftview(centerr, zoomlevel, lofp);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, centerr, 2);
            }
            break;
        case 'H':              /* Show harmonics */
            /* FALLTHRU */
        case 'h':
            {
                double retval;

                retval = harmonic_loop(xid, centerr, zoomlevel, lofp);
                if (retval > 0.0) {
                    offsetf = 0.0;
                    centerr = retval;
                    free(fv);
                    fv = get_fftview(centerr, zoomlevel, lofp);
                    cpgpage();
                    offsetf = plot_fftview(fv, maxpow, 1.0, centerr, 2);
                }
            }
            break;
        case '?':              /* Print help screen */
            print_help();
            break;
        case 'D':              /* Show details about a selected point  */
            /* FALLTHRU */
        case 'd':
            {
                double newr;

                printf("  Searching for peak near freq = %.7g Hz...\n",
                       (inx + offsetf));
                newr = find_peak(inx + offsetf, fv, lofp);
                centerr = newr;
                free(fv);
                fv = get_fftview(centerr, zoomlevel, lofp);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, centerr, 2);
            }
            break;
        case 'L':              /* Load a zaplist */
            /* FALLTHRU */
        case 'l':
            {
                int ii, len;
                char filename[200];
                double *lobins, *hibins;

                printf("  Enter the filename containing the zaplist to load:\n");
                fgets(filename, 199, stdin);
                len = strlen(filename) - 1;
                filename[len] = '\0';
                numzaplist = get_birdies(filename, T, 0.0, &lobins, &hibins);
                lenzaplist = numzaplist + 20;   /* Allow some room to add more */
                if (lenzaplist)
                    free(zaplist);
                zaplist = (bird *) malloc(sizeof(bird) * lenzaplist);
                for (ii = 0; ii < numzaplist; ii++) {
                    zaplist[ii].lobin = lobins[ii];
                    zaplist[ii].hibin = hibins[ii];
                }
                vect_free(lobins);
                vect_free(hibins);
                printf("\n");
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            }
            break;
        case 'Z':              /* Add a birdie to a zaplist */
            /* FALLTHRU */
        case 'z':
            {
                int badchoice = 2;
                float lox, hix, loy, hiy;
                double rs[2];
                char choice;

                if (numzaplist + 1 > lenzaplist) {
                    lenzaplist += 10;
                    zaplist = (bird *) realloc(zaplist, sizeof(bird) * lenzaplist);
                }
                cpgqwin(&lox, &hix, &loy, &hiy);
                printf
                    ("  Click the left mouse button on the first frequency limit.\n");
                while (badchoice) {
                    cpgcurs(&inx, &iny, &choice);
                    if (choice == 'A' || choice == 'a') {
                        rs[2 - badchoice] = ((double) inx + offsetf) * T;
                        cpgsave();
                        cpgsci(7);
                        cpgmove(inx, 0.0);
                        cpgdraw(inx, hiy);
                        cpgunsa();
                        badchoice--;
                        if (badchoice == 1)
                            printf
                                ("  Click the left mouse button on the second frequency limit.\n");
                    } else {
                        printf("  Option not recognized.\n");
                    }
                };
                if (rs[1] > rs[0]) {
                    zaplist[numzaplist].lobin = rs[0];
                    zaplist[numzaplist].hibin = rs[1];
                } else {
                    zaplist[numzaplist].lobin = rs[1];
                    zaplist[numzaplist].hibin = rs[0];
                }
                printf("    The new birdie has:  f_avg = %.15g  f_width = %.15g\n\n",
                       0.5 * (zaplist[numzaplist].hibin +
                              zaplist[numzaplist].lobin) / T,
                       (zaplist[numzaplist].hibin - zaplist[numzaplist].lobin) / T);
                numzaplist++;
                qsort(zaplist, numzaplist, sizeof(bird), compare_birds);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            }
            break;
        case 'P':              /* Print the current plot */
            /* FALLTHRU */
        case 'p':
            {
                int len;
                char filename[200];

                printf("  Enter the filename to save the plot as:\n");
                fgets(filename, 196, stdin);
                len = strlen(filename) - 1;
                filename[len + 0] = '/';
                filename[len + 1] = 'P';
                filename[len + 2] = 'S';
                filename[len + 3] = '\0';
                psid = cpgopen(filename);
                cpgslct(psid);
                cpgpap(10.25, 8.5 / 11.0);
                cpgiden();
                cpgscr(15, 0.8, 0.8, 0.8);
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
                cpgclos();
                cpgslct(xid);
                cpgscr(15, 0.4, 0.4, 0.4);
                filename[len] = '\0';
                printf("  Wrote the plot to the file '%s'.\n", filename);
            }
            break;
        case 'N':              /* Changing power normalization */
            /* FALLTHRU */
        case 'n':
            {
                float inx2 = 0.0, iny2 = 0.0;
                char choice;
                unsigned char badchoice = 1;

                printf("  Specify the type of power normalization:\n"
                       "       m,M  :  Median values determined locally\n"
                       "       d,D  :  DC frequency amplitude\n"
                       "       r,R  :  Raw powers (i.e. no normalization)\n"
                       "       u,U  :  User specified interval (the average powers)\n");
                while (badchoice) {
                    cpgcurs(&inx2, &iny2, &choice);
                    switch (choice) {
                    case 'M':
                        /* FALLTHRU */
                    case 'm':
                        norm_const = 0.0;
                        maxpow = 0.0;
                        badchoice = 0;
                        printf
                            ("  Using local median normalization.  Autoscaling is on.\n");
                        break;
                    case 'D':
                        /* FALLTHRU */
                    case 'd':
                        norm_const = 1.0 / r0;
                        maxpow = 0.0;
                        badchoice = 0;
                        printf
                            ("  Using DC frequency (%f) normalization.  Autoscaling is on.\n",
                             r0);
                        break;
                    case 'R':
                        /* FALLTHRU */
                    case 'r':
                        norm_const = 1.0;
                        maxpow = 0.0;
                        badchoice = 0;
                        printf
                            ("  Using raw powers (i.e. no normalization).  Autoscaling is on.\n");
                        break;
                    case 'U':
                        /* FALLTHRU */
                    case 'u':
                        {
                            char choice2;
                            float xx = inx, yy = iny;
                            long lor, hir, numr;
                            double avg, var;

                            printf
                                ("  Use the left mouse button to select a left and right boundary\n"
                                 "  of a region to calculate the average power.\n");
                            do {
                                cpgcurs(&xx, &yy, &choice2);
                            } while (choice2 != 'A' && choice2 != 'a');
                            lor = (long) ((xx + offsetf) * T);
                            cpgsci(7);
                            cpgmove(xx, 0.0);
                            cpgdraw(xx, 10.0 * fv->maxpow);
                            do {
                                cpgcurs(&xx, &yy, &choice2);
                            } while (choice2 != 'A' && choice2 != 'a');
                            hir = (long) ((xx + offsetf) * T);
                            cpgmove(xx, 0.0);
                            cpgdraw(xx, 10.0 * fv->maxpow);
                            cpgsci(1);
                            if (lor > hir) {
                                long tempr;
                                tempr = hir;
                                hir = lor;
                                lor = tempr;
                            }
                            numr = hir - lor + 1;
                            avg_var(lofp->rawpowers + lor - lofp->rlo, numr, &avg,
                                    &var);
                            printf("  Selection has:  average = %.5g\n"
                                   "                  std dev = %.5g\n", avg,
                                   sqrt(var));
                            norm_const = 1.0 / avg;
                            maxpow = 0.0;
                            badchoice = 0;
                            printf
                                ("  Using %.5g as the normalization constant.  Autoscaling is on.\n",
                                 avg);
                            break;
                        }
                    default:
                        printf("  Unrecognized choice '%c'.\n", choice);
                        break;
                    }
                }
                free(fv);
                fv = get_fftview(centerr, zoomlevel, lofp);
                cpgpage();
                offsetf = plot_fftview(fv, maxpow, 1.0, 0.0, 0);
            }
            break;
        case 'Q':              /* Quit */
            /* FALLTHRU */
        case 'q':
            printf("  Quitting...\n");
            free(fv);
            cpgclos();
            break;
        default:
            printf("  Unrecognized option '%c'.\n", inchar);
            break;
        }
    } while (inchar != 'Q' && inchar != 'q');

    free_fftpart(lofp);
#ifdef USEMMAP
    close(mmap_file);
#else
    fclose(fftfile);
#endif
    if (lenzaplist)
        free(zaplist);
    printf("Done\n\n");
    return 0;
}
