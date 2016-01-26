#include <glib.h>
#include "plot2d.h"
#include "presto.h"
#include "zapbirds_cmd.h"

#define NUMBETWEEN     4
#define FFTLEN         262144
#define BINSTOGET      FFTLEN/NUMBETWEEN
#define MAXPTSTOSHOW   64000
#define MAXBINSTOSHOW  MAXPTSTOSHOW/NUMBETWEEN

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

typedef struct birdie {
    double freq;
    double width;
} birdie;

static FILE *fftfile;
static int khw;
static fcomplex *kernel;
static double T, dr;

extern float calc_median_powers(fcomplex * amplitudes, int numamps);
extern fcomplex *get_rawbins(FILE * fftfile, double bin,
                             int numtoget, float *med, int *lobin);
extern void zapbirds(double lobin, double hibin, FILE * fftfile, fcomplex * fft);

static birdie *birdie_create(double lofreq, double hifreq, double baryv)
/* If baryv corrects barycentric freqs to topocentric */
{
    birdie *obj;

    lofreq /= (1.0 + baryv);
    hifreq /= (1.0 + baryv);
    obj = (birdie *) malloc(sizeof(birdie));
    obj->freq = 0.5 * (hifreq + lofreq);
    obj->width = hifreq - lofreq;
    return obj;
}

static void birdie_free(gpointer data, gpointer user_data)
{
    user_data = NULL;
    free((birdie *) data);
}


static int birdie_compare(gconstpointer ca, gconstpointer cb)
/* Sorts from high to low */
{
    birdie *a, *b;

    a = (birdie *) ca;
    b = (birdie *) cb;
    return (a->freq > b->freq) - (a->freq < b->freq);
}


static void birdie_print(gpointer data, gpointer user_data)
{
    birdie *obj = (birdie *) data;
    FILE *file = (FILE *) user_data;

    fprintf(file, "%17.14g  %17.14g\n", obj->freq, obj->width);
}


static void process_bird(double basebin, int harm, double *lofreq, double *hifreq)
{
    int ii, plotnumpts = 1000, not_done_yet = 1, plotoffset;
    int lodatabin, firstcorrbin, numgoodpts, replot = 1;
    char inchar;
    float med, xx[2], yy[2], inx, iny;
    float powargr, powargi, pwr, maxpow = 0.0, maxbin = 0.0;
    double truebin, pred_freq, average;
    double firstbin = 0.0, lastbin = 0.0, numbins = 0.0;
    fcomplex *data, *result;

    /* 'bin' means normal resolution FFT amplitude */
    /* 'pt'  means an interpolated FFT amplitude   */
    /* 'pts'=='bins' only if NUMBETWEEN==1         */

    *lofreq = *hifreq = 0.0;
    truebin = basebin * harm;
    pred_freq = truebin / T;
    xx[0] = xx[1] = pred_freq;
    data = get_rawbins(fftfile, truebin, BINSTOGET, &med, &lodatabin);
    if (lodatabin <= 0) {
        data[abs(lodatabin)].r = 1.0;
        data[abs(lodatabin)].i = 1.0;
    }
    firstcorrbin = (int) truebin - MAXBINSTOSHOW / 2;
    average = med / -log(0.5);
    result = gen_cvect(FFTLEN);
    numgoodpts = corr_complex(data, BINSTOGET, RAW,
                              kernel, FFTLEN, FFT,
                              result, MAXPTSTOSHOW,
                              firstcorrbin - lodatabin, NUMBETWEEN, khw, CORR);
    for (ii = 0; ii < numgoodpts; ii++) {
        pwr = POWER(result[ii].r, result[ii].i) / average;
        if (pwr > maxpow) {
            maxpow = pwr;
            maxbin = firstcorrbin + dr * ii;
        }
    }
    printf("\nHarmonic %d of %.15g Hz (%.15g Hz, bin = %.15g)\n",
           harm, basebin / T, pred_freq, truebin);
    printf("  Max power = %.2f at %.15g Hz (bin = %.15g)\n",
           maxpow, maxbin / T, maxbin);
    do {
        cpgsci(1);
        if (replot) {
            plotoffset = MAXPTSTOSHOW / 2 - plotnumpts / 2;
            firstbin = firstcorrbin + dr * plotoffset;
            numbins = dr * plotnumpts;
            lastbin = firstbin + numbins;
            plot_spectrum(result + plotoffset, plotnumpts, firstbin, dr, T, average);
            cpgswin(0.0, 1.0, 0.0, 1.0);
            xx[0] = xx[1] = (truebin - firstbin) / numbins;
            yy[0] = 0.0;
            yy[1] = 1.0;
            cpgsci(2);          /* red */
            cpgline(2, xx, yy); /* Predicted freq */
            cpgsci(7);          /* yellow */
            if (*lofreq) {
                xx[0] = xx[1] = ((*lofreq * T) - firstbin) / numbins;
                cpgline(2, xx, yy);     /* Boundary */
            }
            if (*hifreq) {
                xx[0] = xx[1] = ((*hifreq * T) - firstbin) / numbins;
                cpgline(2, xx, yy);     /* Boundary */
            }
        }
        replot = 1;
        cpgsci(7);              /* yellow */
        cpgcurs(&inx, &iny, &inchar);
        switch (inchar) {
        case ' ':
        case 'A':
        case 'a':
            xx[0] = xx[1] = inx;
            cpgline(2, xx, yy); /* New boundary */
            if (*lofreq == 0.0) {
                *lofreq = (inx * numbins + firstbin) / T;
                printf("  Added 1st boundary at %.12g Hz\n", *lofreq);
            } else {
                *hifreq = (inx * numbins + firstbin) / T;
                printf("  Added 2nd boundary at %.12g Hz\n", *hifreq);
            }
            replot = 0;
            break;
        case 'I':              /* Zoom in */
        case 'i':
            plotnumpts /= 2;
            if (plotnumpts <= 8)
                plotnumpts = 8;
            printf("  Zooming in...\n");
            break;
        case 'O':              /* Zoom out */
        case 'o':
            plotnumpts *= 2;
            if (plotnumpts > MAXPTSTOSHOW)
                plotnumpts = MAXPTSTOSHOW;
            printf("  Zooming out...\n");
            break;
        case 'C':              /* Clear/Delete the points */
        case 'c':
        case 'D':
        case 'd':
        case 'X':
        case 'x':
            *lofreq = *hifreq = 0.0;
            printf("  Clearing boundaries.\n");
            break;
        case 'Q':              /* Quit/Next birdie */
        case 'q':
        case 'N':
        case 'n':
            *lofreq = *hifreq = 0.0;
            free(data);
            vect_free(result);
            printf("  Skipping to next harmonic.\n");
            return;
        default:
            printf("  Unrecognized option '%c'.\n", inchar);
            break;
        }
        if (*lofreq && *hifreq)
            not_done_yet = 0;
    } while (not_done_yet);
    if (*hifreq < *lofreq) {
        double tmpfreq;
        tmpfreq = *lofreq;
        *lofreq = *hifreq;
        *hifreq = tmpfreq;
    }
    free(data);
    vect_free(result);
}


int main(int argc, char *argv[])
{
    int ii, jj, numbirds;
    double lofreq, hifreq;
    char *rootfilenm;
    birdie *newbird;
    GSList *zapped = NULL;
    infodata idata;
    Cmdline *cmd;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("     Interactive/Automatic Birdie Zapping Program\n");
    printf("              by Scott M. Ransom\n");
    printf("                 January, 2001\n\n");

    if (!cmd->zapP && !cmd->inzapfileP && !cmd->outzapfileP) {
        printf("You must specify '-in' and '-out' if you are not\n");
        printf("automatically zapping a file (with '-zap').\n\n");
        exit(0);
    }

    {
        int hassuffix = 0;
        char *suffix;

        hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);
        if (hassuffix) {
            if (strcmp(suffix, "fft") != 0) {
                printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
                       cmd->argv[0]);
                free(suffix);
                exit(0);
            }
            free(suffix);
        } else {
            printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
                   cmd->argv[0]);
            exit(0);
        }
    }

    /* Read the info file */

    readinf(&idata, rootfilenm);
    if (strlen(remove_whitespace(idata.object)) > 0) {
        printf("Examining %s data from '%s'.\n\n",
               remove_whitespace(idata.object), cmd->argv[0]);
    } else {
        printf("Examining data from '%s'.\n\n", cmd->argv[0]);
    }
    T = idata.dt * idata.N;
    dr = 1.0 / NUMBETWEEN;

    if (cmd->zapP) {            /* Automatic  */
        double *bird_lobins, *bird_hibins, hibin;

        if (!cmd->zapfileP) {
            printf("You must specify a 'zapfile' containing freqs\n");
            printf("and widths if you want to write to the FFT file.\n\n");
            free(rootfilenm);
            exit(0);
        }
        hibin = idata.N / 2;

        /* Read the Standard bird list */

        numbirds = get_birdies(cmd->zapfile, T, cmd->baryv,
                               &bird_lobins, &bird_hibins);

        /* Zap the birdies */

        fftfile = chkfopen(cmd->argv[0], "rb+");
        for (ii = 0; ii < numbirds; ii++) {
            if (bird_lobins[ii] >= hibin - 1)
                break;
            if (bird_hibins[ii] >= hibin)
                bird_hibins[ii] = hibin - 1;
            zapbirds(bird_lobins[ii], bird_hibins[ii], fftfile, NULL);
        }

        vect_free(bird_lobins);
        vect_free(bird_hibins);

    } else {                    /* Interactive */

        int *bird_numharms;
        double *bird_basebins;

        /* Read the Standard bird list */

        numbirds = get_std_birds(cmd->inzapfile, T, cmd->baryv,
                                 &bird_basebins, &bird_numharms);

        /* Create our correlation kernel */

        {
            int numkern;
            fcomplex *resp;

            khw = r_resp_halfwidth(LOWACC);
            numkern = 2 * NUMBETWEEN * khw;
            resp = gen_r_response(0.0, NUMBETWEEN, numkern);
            kernel = gen_cvect(FFTLEN);
            place_complex_kernel(resp, numkern, kernel, FFTLEN);
            COMPLEXFFT(kernel, FFTLEN, -1);
            vect_free(resp);
        }

        /* Loop over the birdies */

        fftfile = chkfopen(cmd->argv[0], "rb");
        cpgstart_x("landscape");
        cpgask(0);
        for (ii = 0; ii < numbirds; ii++) {
            for (jj = 0; jj < bird_numharms[ii]; jj++) {
                process_bird(bird_basebins[ii], jj + 1, &lofreq, &hifreq);
                if (lofreq && hifreq) {
                    newbird = birdie_create(lofreq, hifreq, cmd->baryv);
                    zapped = g_slist_insert_sorted(zapped, newbird, birdie_compare);
                }
            }
        }
        cpgclos();

        /* Output the birdies */

        {
            FILE *outfile;

            outfile = chkfopen(cmd->outzapfile, "w");
            fprintf(outfile, "#\n");
            fprintf(outfile,
                    "# Topocentric birdies found using 'zapbirds' for '%s'\n",
                    cmd->argv[0]);
            fprintf(outfile, "#\n");
            fprintf(outfile, "# Frequency (Hz)       Width (Hz)\n");
            fprintf(outfile, "#\n");
            g_slist_foreach(zapped, birdie_print, outfile);
            fclose(outfile);
        }

        printf("\nOutput birdie file is '%s'.\n\n", cmd->outzapfile);

        /* Free the memory */

        g_slist_foreach(zapped, birdie_free, NULL);
        g_slist_free(zapped);
        vect_free(kernel);
        vect_free(bird_numharms);
        vect_free(bird_basebins);
    }
    fclose(fftfile);
    free(rootfilenm);
    printf("Done\n\n");
    return 0;
}
