#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "rfifind_cmd.h"
#include "mask.h"
#include "backend_common.h"
#include "rfifind.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

/* Some function definitions */

void rfifind_plot(int numchan, int numint, int ptsperint,
                  float timesigma, float freqsigma,
                  float inttrigfrac, float chantrigfrac,
                  float **dataavg, float **datastd, float **datapow,
                  int *userchan, int numuserchan,
                  int *userints, int numuserints,
                  infodata * idata, unsigned char **bytemask,
                  mask * oldmask, mask * newmask,
                  rfi * rfivect, int numrfi, int rfixwin, int rfips, int xwin);
static void write_rfifile(char *rfifilenm, rfi * rfivect, int numrfi,
                          int numchan, int numint, int ptsperint,
                          int lobin, int numbetween, int harmsum,
                          float fracterror, float freqsigma);
static void write_statsfile(char *statsfilenm, float *datapow,
                            float *dataavg, float *datastd,
                            int numchan, int numint, int ptsperint,
                            int lobin, int numbetween);
static void read_rfifile(char *rfifilenm, rfi ** rfivect, int *numrfi,
                         int *numchan, int *numint, int *ptsperint,
                         int *lobin, int *numbetween, int *harmsum,
                         float *fracterror, float *freqsigma);
static void read_statsfile(char *statsfilenm, float ***datapow,
                           float ***dataavg, float ***datastd,
                           int *numchan, int *numint, int *ptsperint,
                           int *lobin, int *numbetween);
static int *get_prime_factors(int nn, int *numfactors);
int compare_rfi_sigma(const void *ca, const void *cb);
int compare_rfi_numobs(const void *ca, const void *cb);
int read_subband_rawblocks(FILE * infiles[], int numfiles, short *subbanddata,
                           int numsamples, int *padding);
void get_subband(int subbandnum, float chandat[], short srawdata[], int numsamples);
extern int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals);

/* The main program */

int main(int argc, char *argv[])
{
    FILE *bytemaskfile;
    float **dataavg = NULL, **datastd = NULL, **datapow = NULL;
    float *chandata = NULL, powavg, powstd, powmax;
    float inttime, norm = 0.0, fracterror = RFI_FRACTERROR;
    float *rawdata = NULL;
    unsigned char **bytemask = NULL;
    short *srawdata = NULL;
    char *outfilenm, *statsfilenm, *maskfilenm;
    char *bytemaskfilenm, *rfifilenm;
    int numchan = 0, numint = 0, newper = 0, oldper = 0;
    int blocksperint, ptsperint = 0, ptsperblock = 0, padding = 0;
    int numcands, candnum, numrfi = 0, numrfivect = NUM_RFI_VECT;
    int ii, jj, kk, slen, insubs = 0;
    int harmsum = RFI_NUMHARMSUM, lobin = RFI_LOBIN, numbetween = RFI_NUMBETWEEN;
    double davg, dvar, freq;
    struct spectra_info s;
    presto_interptype interptype;
    rfi *rfivect = NULL;
    mask oldmask, newmask;
    fftcand *cands = NULL;
    infodata idata;
    Cmdline *cmd;

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        printf("\n");
        usage();
        exit(0);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);
    spectra_info_set_defaults(&s);
    s.filenames = cmd->argv;
    s.num_files = cmd->argc;
    // If we are zeroDMing, make sure that clipping is off.
    if (cmd->zerodmP)
        cmd->noclipP = 1;
    s.clip_sigma = cmd->clip;
    // -1 causes the data to determine if we use weights, scales, &
    // offsets for PSRFITS or flip the band for any data type where
    // we can figure that out with the data
    s.apply_flipband = (cmd->invertP) ? 1 : -1;
    s.apply_weight = (cmd->noweightsP) ? 0 : -1;
    s.apply_scale = (cmd->noscalesP) ? 0 : -1;
    s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
    s.remove_zerodm = (cmd->zerodmP) ? 1 : 0;
    if (cmd->noclipP) {
        cmd->clip = 0.0;
        s.clip_sigma = 0.0;
    }
    if (cmd->ifsP) {
        // 0 = default or summed, 1-4 are possible also
        s.use_poln = cmd->ifs + 1;
    }
    slen = strlen(cmd->outfile) + 20;

    if (cmd->ncpus > 1) {
#ifdef _OPENMP
        int maxcpus = omp_get_num_procs();
        int openmp_numthreads = (cmd->ncpus <= maxcpus) ? cmd->ncpus : maxcpus;
        // Make sure we are not dynamically setting the number of threads
        omp_set_dynamic(0);
        omp_set_num_threads(openmp_numthreads);
        printf("Using %d threads with OpenMP\n\n", openmp_numthreads);
#endif
    } else {
#ifdef _OPENMP
        omp_set_num_threads(1); // Explicitly turn off OpenMP
#endif
    }

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n\n");
    printf("               Pulsar Data RFI Finder\n");
    printf("                 by Scott M. Ransom\n\n");

    /* The following is the root of all the output files */

    outfilenm = (char *) calloc(slen, sizeof(char));
    sprintf(outfilenm, "%s_rfifind", cmd->outfile);

    /* And here are the output file names */

    maskfilenm = (char *) calloc(slen, sizeof(char));
    sprintf(maskfilenm, "%s.mask", outfilenm);
    bytemaskfilenm = (char *) calloc(slen, sizeof(char));
    sprintf(bytemaskfilenm, "%s.bytemask", outfilenm);
    rfifilenm = (char *) calloc(slen, sizeof(char));
    sprintf(rfifilenm, "%s.rfi", outfilenm);
    statsfilenm = (char *) calloc(slen, sizeof(char));
    sprintf(statsfilenm, "%s.stats", outfilenm);
    sprintf(idata.name, "%s", outfilenm);

    if (RAWDATA) {
        if (cmd->filterbankP)
            s.datatype = SIGPROCFB;
        else if (cmd->psrfitsP)
            s.datatype = PSRFITS;
    } else {                    // Attempt to auto-identify the data
        identify_psrdatatype(&s, 1);
        if (s.datatype == SIGPROCFB)
            cmd->filterbankP = 1;
        else if (s.datatype == PSRFITS)
            cmd->psrfitsP = 1;
        else if (s.datatype == SUBBAND)
            insubs = 1;
        else {
            printf
                ("Error:  Unable to identify input data files.  Please specify type.\n\n");
            exit(1);
        }
    }

    /* Read an input mask if wanted */
    if (cmd->maskfileP) {
        read_mask(cmd->maskfile, &oldmask);
        printf("Read old mask information from '%s'\n\n", cmd->maskfile);
    } else {
        oldmask.numchan = oldmask.numint = 0;
    }

    if (!cmd->nocomputeP) {

        if (RAWDATA || insubs) {
            char description[40];
            psrdatatype_description(description, s.datatype);
            if (s.num_files > 1)
                printf("Reading %s data from %d files:\n", description, s.num_files);
            else
                printf("Reading %s data from 1 file:\n", description);
            if (insubs)
                s.files = (FILE **) malloc(sizeof(FILE *) * s.num_files);
            for (ii = 0; ii < s.num_files; ii++) {
                printf("  '%s'\n", cmd->argv[ii]);
                if (insubs)
                    s.files[ii] = chkfopen(cmd->argv[ii], "rb");
            }
            printf("\n");
        }

        if (RAWDATA) {
            read_rawdata_files(&s);
            if (cmd->ignorechanstrP) {
                s.ignorechans = get_ignorechans(cmd->ignorechanstr, 0, s.num_channels-1,
                                                &s.num_ignorechans, &s.ignorechans_str);
                if (s.ignorechans_str==NULL) {
                    s.ignorechans_str = (char *)malloc(strlen(cmd->ignorechanstr)+1);
                    strcpy(s.ignorechans_str, cmd->ignorechanstr);
                }
            }
            print_spectra_info_summary(&s);
            spectra_info_to_inf(&s, &idata);
            ptsperblock = s.spectra_per_subint;
            numchan = s.num_channels;
            idata.dm = 0.0;
            writeinf(&idata);
        }

        if (insubs) {
            /* Set-up values if we are using subbands */
            char *tmpname, *root, *suffix;
            if (split_root_suffix(s.filenames[0], &root, &suffix) == 0) {
                printf("Error:  The input filename (%s) must have a suffix!\n\n",
                       s.filenames[0]);
                exit(1);
            }
            if (strncmp(suffix, "sub", 3) == 0) {
                tmpname = calloc(strlen(root) + 6, 1);
                sprintf(tmpname, "%s.sub", root);
                readinf(&idata, tmpname);
                free(tmpname);
            } else {
                printf
                    ("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                     s.filenames[0]);
                exit(1);
            }
            free(root);
            free(suffix);
            ptsperblock = 1;
            /* Compensate for the fact that we have subbands and not channels */
            idata.freq = idata.freq - 0.5 * idata.chan_wid +
                0.5 * idata.chan_wid * (idata.num_chan / s.num_files);
            idata.chan_wid = idata.num_chan / s.num_files * idata.chan_wid;
            idata.num_chan = numchan = s.num_files;
            idata.dm = 0.0;
            sprintf(idata.name, "%s", outfilenm);
            writeinf(&idata);
            s.padvals = gen_fvect(s.num_files);
            for (ii = 0; ii < s.num_files; ii++)
                s.padvals[ii] = 0.0;
        }

        if (cmd->maskfileP)
            determine_padvals(cmd->maskfile, &oldmask, s.padvals);

        /* The number of data points and blocks to work with at a time */

        if (cmd->blocksP) {
            blocksperint = cmd->blocks;
            cmd->time = blocksperint * ptsperblock * idata.dt;
        } else {
            blocksperint = (int) (cmd->time / (ptsperblock * idata.dt) + 0.5);
            // Must process at least 1 block at a time
            if (blocksperint == 0)
                blocksperint = 1;
        }
        ptsperint = blocksperint * ptsperblock;
        numint = (long long) idata.N / ptsperint;
        if ((long long) idata.N % ptsperint)
            numint++;
        inttime = ptsperint * idata.dt;
        printf("Analyzing data sections of length %d points (%.6g sec).\n",
               ptsperint, inttime);
        {
            int *factors, numfactors;

            factors = get_prime_factors(ptsperint, &numfactors);
            printf("  Prime factors are:  ");
            for (ii = 0; ii < numfactors; ii++)
                printf("%d ", factors[ii]);
            printf("\n");
            if (factors[numfactors - 1] > 13) {
                printf
                    ("  WARNING:  The largest prime factor is pretty big!  This will\n"
                     "            cause the FFTs to take a long time to compute.  I\n"
                     "            recommend choosing a different -time value.\n");
            }
            printf("\n");
            free(factors);
        }

        /* Allocate our workarrays */

        if (RAWDATA)
            rawdata = gen_fvect(idata.num_chan * ptsperblock * blocksperint);
        else if (insubs)
            srawdata = gen_svect(idata.num_chan * ptsperblock * blocksperint);
        dataavg = gen_fmatrix(numint, numchan);
        datastd = gen_fmatrix(numint, numchan);
        datapow = gen_fmatrix(numint, numchan);
        chandata = gen_fvect(ptsperint);
        bytemask = gen_bmatrix(numint, numchan);
        for (ii = 0; ii < numint; ii++)
            for (jj = 0; jj < numchan; jj++)
                bytemask[ii][jj] = GOODDATA;
        rfivect = rfi_vector(rfivect, numchan, numint, 0, numrfivect);
        if (numbetween == 2)
            interptype = INTERBIN;
        else
            interptype = INTERPOLATE;

        /* Main loop */

        printf("Writing mask data  to '%s'.\n", maskfilenm);
        printf("Writing  RFI data  to '%s'.\n", rfifilenm);
        printf("Writing statistics to '%s'.\n\n", statsfilenm);
        printf("Massaging the data ...\n\n");
        printf("Amount Complete = %3d%%", oldper);
        fflush(stdout);

        for (ii = 0; ii < numint; ii++) {       /* Loop over the intervals */
            newper = (int) ((float) ii / numint * 100.0 + 0.5);
            if (newper > oldper) {
                printf("\rAmount Complete = %3d%%", newper);
                fflush(stdout);
                oldper = newper;
            }

            /* Read a chunk of data */

            if (RAWDATA) {
                read_rawblocks(rawdata, blocksperint, &s, &padding);
                // Clip nasty RFI if requested (we are not masking)
                if (s.clip_sigma > 0.0)
                    clip_times(rawdata, ptsperint, s.num_channels, s.clip_sigma,
                               s.padvals);
            } else if (insubs) {
                read_subband_rawblocks(s.files, s.num_files,
                                       srawdata, blocksperint, &padding);
                // TODO: should implement clipping for subbands
            }

            if (padding)
                for (jj = 0; jj < numchan; jj++)
                    bytemask[ii][jj] |= PADDING;

            for (jj = 0; jj < numchan; jj++) {  /* Loop over the channels */

                if (RAWDATA)
                    get_channel(chandata, jj, blocksperint, rawdata, &s);
                else if (insubs)
                    get_subband(jj, chandata, srawdata, blocksperint);

                /* Calculate the averages and standard deviations */
                /* for each point in time.                        */

                if (padding) {
                    dataavg[ii][jj] = 0.0;
                    datastd[ii][jj] = 0.0;
                    datapow[ii][jj] = 1.0;
                } else {
                    avg_var(chandata, ptsperint, &davg, &dvar);
                    dataavg[ii][jj] = davg;
                    datastd[ii][jj] = sqrt(dvar);
                    numcands = 0;
                    powmax = 0.0;
                    // Don't search the power spectrum if there is little to no variance
                    if (datastd[ii][jj] > 1e-4) {
                        realfft(chandata, ptsperint, -1);
                        norm = datastd[ii][jj] * datastd[ii][jj] * ptsperint;
                        cands = search_fft((fcomplex *) chandata, ptsperint / 2,
                                           lobin, ptsperint / 2, harmsum,
                                           numbetween, interptype, norm, cmd->freqsigma,
                                           &numcands, &powavg, &powstd, &powmax);
                        // Make sure that nothing bad happened in the FFT search
                        if (!isnormal(powmax)) {
                            printf("WARNING:  FFT search returned bad powmax (%f) in"
                                   "int=%d and chan=%d.  Fixing.\n",
                                   powmax, ii, jj);
                            powmax = 0.0;
                            numcands = 0;
                            free(cands);
                        }
                    }
                    datapow[ii][jj] = powmax;

                    /* Record the birdies */

                    if (numcands) {
                        for (kk = 0; kk < numcands; kk++) {
                            freq = cands[kk].r / inttime;
                            candnum =
                                find_rfi(rfivect, numrfi, freq, RFI_FRACTERROR);
                            if (candnum >= 0) {
                                update_rfi(rfivect + candnum, freq, cands[kk].sig,
                                           jj, ii);
                            } else {
                                update_rfi(rfivect + numrfi, freq, cands[kk].sig, jj,
                                           ii);
                                numrfi++;
                                if (numrfi == numrfivect) {
                                    numrfivect *= 2;
                                    rfivect = rfi_vector(rfivect, numchan, numint,
                                                         numrfivect / 2, numrfivect);
                                }
                            }
                        }
                        free(cands);
                    }
                }
            }
        }
        printf("\rAmount Complete = 100%%\n");

        /* Write the data to the output files */

        write_rfifile(rfifilenm, rfivect, numrfi, numchan, numint,
                      ptsperint, lobin, numbetween, harmsum,
                      fracterror, cmd->freqsigma);
        write_statsfile(statsfilenm, datapow[0], dataavg[0], datastd[0],
                        numchan, numint, ptsperint, lobin, numbetween);

    } else {                    /* If "-nocompute" */
        float freqsigma;

        /* Read the data from the output files */

        printf("Reading  RFI data  from '%s'.\n", rfifilenm);
        printf("Reading statistics from '%s'.\n", statsfilenm);
        readinf(&idata, outfilenm);
        read_rfifile(rfifilenm, &rfivect, &numrfi, &numchan, &numint,
                     &ptsperint, &lobin, &numbetween, &harmsum,
                     &fracterror, &freqsigma);
        numrfivect = numrfi;
        read_statsfile(statsfilenm, &datapow, &dataavg, &datastd,
                       &numchan, &numint, &ptsperint, &lobin, &numbetween);
        bytemask = gen_bmatrix(numint, numchan);
        printf("Reading  bytemask  from '%s'.\n\n", bytemaskfilenm);
        bytemaskfile = chkfopen(bytemaskfilenm, "rb");
        chkfread(bytemask[0], numint * numchan, 1, bytemaskfile);
        fclose(bytemaskfile);
        for (ii = 0; ii < numint; ii++)
            for (jj = 0; jj < numchan; jj++)
                bytemask[ii][jj] &= PADDING;    /* Clear all but the PADDING bits */
        inttime = ptsperint * idata.dt;
    }

    /* Make the plots and set the mask */

    {
        int *zapints, *zapchan;
        int numzapints = 0, numzapchan = 0;

        if (cmd->zapintsstrP) {
            zapints = ranges_to_ivect(cmd->zapintsstr, 0, numint - 1, &numzapints);
            zapints = (int *) realloc(zapints, (size_t) (sizeof(int) * numint));
        } else {
            zapints = gen_ivect(numint);
        }
        if (cmd->zapchanstrP) {
            zapchan = ranges_to_ivect(cmd->zapchanstr, 0, numchan - 1, &numzapchan);
            zapchan = (int *) realloc(zapchan, (size_t) (sizeof(int) * numchan));
        } else {
            zapchan = gen_ivect(numchan);
        }
        rfifind_plot(numchan, numint, ptsperint, cmd->timesigma, cmd->freqsigma,
                     cmd->inttrigfrac, cmd->chantrigfrac,
                     dataavg, datastd, datapow, zapchan, numzapchan,
                     zapints, numzapints, &idata, bytemask,
                     &oldmask, &newmask, rfivect, numrfi,
                     cmd->rfixwinP, cmd->rfipsP, cmd->xwinP);

        vect_free(zapints);
        vect_free(zapchan);
    }

    /* Write the new mask and bytemask to the file */

    write_mask(maskfilenm, &newmask);
    bytemaskfile = chkfopen(bytemaskfilenm, "wb");
    chkfwrite(bytemask[0], numint * numchan, 1, bytemaskfile);
    fclose(bytemaskfile);

    /* Determine the percent of good and bad data */

    {
        int numpad = 0, numbad = 0, numgood = 0;

        for (ii = 0; ii < numint; ii++) {
            for (jj = 0; jj < numchan; jj++) {
                if (bytemask[ii][jj] == GOODDATA) {
                    numgood++;
                } else {
                    if (bytemask[ii][jj] & PADDING)
                        numpad++;
                    else
                        numbad++;
                }
            }
        }
        printf("\nTotal number of intervals in the data:  %d\n\n", numint * numchan);
        printf("  Number of padded intervals:  %7d  (%6.3f%%)\n",
               numpad, (float) numpad / (float) (numint * numchan) * 100.0);
        printf("  Number of  good  intervals:  %7d  (%6.3f%%)\n",
               numgood, (float) numgood / (float) (numint * numchan) * 100.0);
        printf("  Number of  bad   intervals:  %7d  (%6.3f%%)\n\n",
               numbad, (float) numbad / (float) (numint * numchan) * 100.0);
        qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_sigma);
        printf("  Ten most significant birdies:\n");
        printf("#  Sigma     Period(ms)      Freq(Hz)       Number \n");
        printf("----------------------------------------------------\n");
        for (ii = 0; ii < 10; ii++) {
            double pperr;
            char temp1[40], temp2[40];

            if (rfivect[ii].freq_var == 0.0) {
                pperr = 0.0;
                sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
                sprintf(temp2, " %-14g", 1000.0 / rfivect[ii].freq_avg);
            } else {
                pperr = 1000.0 * sqrt(rfivect[ii].freq_var) /
                    (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
                nice_output_2(temp1, rfivect[ii].freq_avg,
                              sqrt(rfivect[ii].freq_var), -15);
                nice_output_2(temp2, 1000.0 / rfivect[ii].freq_avg, pperr, -15);
            }
            printf("%-2d %-8.2f %13s %13s %-8d\n", ii + 1, rfivect[ii].sigma_avg,
                   temp2, temp1, rfivect[ii].numobs);
        }
        qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_numobs);
        printf("\n  Ten most numerous birdies:\n");
        printf("#  Number    Period(ms)      Freq(Hz)       Sigma \n");
        printf("----------------------------------------------------\n");
        for (ii = 0; ii < 10; ii++) {
            double pperr;
            char temp1[40], temp2[40];

            if (rfivect[ii].freq_var == 0.0) {
                pperr = 0.0;
                sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
                sprintf(temp2, " %-14g", 1000.0 / rfivect[ii].freq_avg);
            } else {
                pperr = 1000.0 * sqrt(rfivect[ii].freq_var) /
                    (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
                nice_output_2(temp1, rfivect[ii].freq_avg,
                              sqrt(rfivect[ii].freq_var), -15);
                nice_output_2(temp2, 1000.0 / rfivect[ii].freq_avg, pperr, -15);
            }
            printf("%-2d %-8d %13s %13s %-8.2f\n", ii + 1, rfivect[ii].numobs,
                   temp2, temp1, rfivect[ii].sigma_avg);
        }
        printf("\nDone.\n\n");
    }

    /* Close the files and cleanup */

    free_rfi_vector(rfivect, numrfivect);
    free_mask(newmask);
    if (cmd->maskfileP)
        free_mask(oldmask);
    free(outfilenm);
    free(statsfilenm);
    free(bytemaskfilenm);
    free(maskfilenm);
    free(rfifilenm);
    vect_free(dataavg[0]);
    vect_free(dataavg);
    vect_free(datastd[0]);
    vect_free(datastd);
    vect_free(datapow[0]);
    vect_free(datapow);
    vect_free(bytemask[0]);
    vect_free(bytemask);
    if (!cmd->nocomputeP) {
        //  Close all the raw files and free their vectors
        close_rawfiles(&s);
        vect_free(chandata);
        if (insubs)
            vect_free(srawdata);
        else
            vect_free(rawdata);
    }
    return (0);
}

static void write_rfifile(char *rfifilenm, rfi * rfivect, int numrfi,
                          int numchan, int numint, int ptsperint,
                          int lobin, int numbetween, int harmsum,
                          float fracterror, float freqsigma)
{
    FILE *outfile;
    int ii;

    outfile = chkfopen(rfifilenm, "wb");
    chkfwrite(&numchan, sizeof(int), 1, outfile);
    chkfwrite(&numint, sizeof(int), 1, outfile);
    chkfwrite(&ptsperint, sizeof(int), 1, outfile);
    chkfwrite(&lobin, sizeof(int), 1, outfile);
    chkfwrite(&numbetween, sizeof(int), 1, outfile);
    chkfwrite(&harmsum, sizeof(int), 1, outfile);
    chkfwrite(&numrfi, sizeof(int), 1, outfile);
    chkfwrite(&fracterror, sizeof(float), 1, outfile);
    chkfwrite(&freqsigma, sizeof(float), 1, outfile);
    for (ii = 0; ii < numrfi; ii++)
        write_rfi(outfile, rfivect + ii, numchan, numint);
    fclose(outfile);
}

static void write_statsfile(char *statsfilenm, float *datapow,
                            float *dataavg, float *datastd,
                            int numchan, int numint, int ptsperint,
                            int lobin, int numbetween)
{
    FILE *outfile;

    outfile = chkfopen(statsfilenm, "wb");
    chkfwrite(&numchan, sizeof(int), 1, outfile);
    chkfwrite(&numint, sizeof(int), 1, outfile);
    chkfwrite(&ptsperint, sizeof(int), 1, outfile);
    chkfwrite(&lobin, sizeof(int), 1, outfile);
    chkfwrite(&numbetween, sizeof(int), 1, outfile);
    chkfwrite(datapow, sizeof(float), numchan * numint, outfile);
    chkfwrite(dataavg, sizeof(float), numchan * numint, outfile);
    chkfwrite(datastd, sizeof(float), numchan * numint, outfile);
    fclose(outfile);
}

static void read_rfifile(char *rfifilenm, rfi ** rfivect, int *numrfi,
                         int *numchan, int *numint, int *ptsperint,
                         int *lobin, int *numbetween, int *harmsum,
                         float *fracterror, float *freqsigma)
{
    FILE *outfile;
    int ii;

    outfile = chkfopen(rfifilenm, "rb");
    chkfread(numchan, sizeof(int), 1, outfile);
    chkfread(numint, sizeof(int), 1, outfile);
    chkfread(ptsperint, sizeof(int), 1, outfile);
    chkfread(lobin, sizeof(int), 1, outfile);
    chkfread(numbetween, sizeof(int), 1, outfile);
    chkfread(harmsum, sizeof(int), 1, outfile);
    chkfread(numrfi, sizeof(int), 1, outfile);
    chkfread(fracterror, sizeof(float), 1, outfile);
    chkfread(freqsigma, sizeof(float), 1, outfile);
    *rfivect = rfi_vector(*rfivect, *numchan, *numint, 0, *numrfi);
    for (ii = 0; ii < *numrfi; ii++)
        read_rfi(outfile, *rfivect + ii, *numchan, *numint);
    fclose(outfile);
}

static void read_statsfile(char *statsfilenm, float ***datapow,
                           float ***dataavg, float ***datastd,
                           int *numchan, int *numint, int *ptsperint,
                           int *lobin, int *numbetween)
{
    FILE *outfile;

    outfile = chkfopen(statsfilenm, "rb");
    chkfread(numchan, sizeof(int), 1, outfile);
    chkfread(numint, sizeof(int), 1, outfile);
    chkfread(ptsperint, sizeof(int), 1, outfile);
    chkfread(lobin, sizeof(int), 1, outfile);
    chkfread(numbetween, sizeof(int), 1, outfile);
    *dataavg = gen_fmatrix(*numint, *numchan);
    *datastd = gen_fmatrix(*numint, *numchan);
    *datapow = gen_fmatrix(*numint, *numchan);
    chkfread(*(datapow[0]), sizeof(float), *numchan * *numint, outfile);
    chkfread(*(dataavg[0]), sizeof(float), *numchan * *numint, outfile);
    chkfread(*(datastd[0]), sizeof(float), *numchan * *numint, outfile);
    fclose(outfile);
}

int read_subband_rawblocks(FILE * infiles[], int numfiles, short *subbanddata,
                           int numsamples, int *padding)
{
    int ii, jj, index, numread = 0;

    for (ii = 0; ii < numfiles; ii++) {
        index = ii * numsamples;
        numread =
            chkfread(subbanddata + index, sizeof(short), numsamples, infiles[ii]);
        for (jj = numread; jj < numsamples; jj++)
            subbanddata[index + jj] = 0.0;
    }
    /* Subband data cannot currently be padded */
    *padding = 0;
    return numread;
}

void get_subband(int subbandnum, float chandat[], short srawdata[], int numsamples)
{
    int ii, index;

    index = subbandnum * numsamples;
    for (ii = 0; ii < numsamples; ii++, index++)
        chandat[ii] = (float) srawdata[index];
}

static int *get_prime_factors(int nn, int *numfactors)
{
    int ii, factors[30], *cfactors;

    *numfactors = 0;
    while ((nn & 1) == 0) {
        nn >>= 1;
        factors[(*numfactors)++] = 2;
    }
    if (nn == 1)
        goto finish;

    for (ii = 3; ii * ii <= nn; ii += 2) {
        while (!(nn % ii)) {
            nn /= ii;
            factors[(*numfactors)++] = ii;
        }
    }
    if (nn == 1)
        goto finish;

    factors[(*numfactors)++] = nn;

  finish:
    cfactors = (int *) malloc(*numfactors * sizeof(int));
    for (ii = 0; ii < *numfactors; ii++)
        cfactors[ii] = factors[ii];
    return cfactors;
}
