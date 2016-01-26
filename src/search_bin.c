#include "presto.h"
#include "search_bin_cmd.h"

/* The number of candidates to return from the search of each miniFFT */
#define MININCANDS 6

/* Minimum binary period (s) to accept as 'real' */
#define MINORBP 300.0

/* Function definitions */
int not_already_there_rawbin(rawbincand newcand, rawbincand * list, int nlist);
int comp_rawbin_to_cand(rawbincand * cand, infodata * idata, char *output, int full);
void compare_rawbin_cands(rawbincand * list, int nlist, char *notes);
void file_rawbin_candidates(rawbincand * cand, char *notes,
                            int numcands, int numharm, char name[]);
float percolate_rawbincands(rawbincand * cands, int numcands);

/* Main routine */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *fftfile, *candfile;
    float powargr, powargi, *powers = NULL, *minifft;
    float norm, numchunks, *powers_pos;
    int nbins, newncand, nfftsizes, fftlen, halffftlen, binsleft;
    int numtoread, filepos = 0, loopct = 0, powers_offset, ncand2;
    int ii, ct, newper = 0, oldper = 0, numsumpow = 1;
    double T, totnumsearched = 0.0, minsig = 0.0, min_orb_p, max_orb_p;
    char *rootfilenm, *notes;
    fcomplex *data = NULL;
    rawbincand tmplist[MININCANDS], *list;
    infodata idata;
    struct tms runtimes;
    double ttim, utim, stim, tott;
    Cmdline *cmd;
    fftwf_plan fftplan;

    /* Prep the timer */

    tott = times(&runtimes) / (double) CLK_TCK;

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
    printf("     Phase Modulation Pulsar Search Routine\n");
    printf("              by Scott M. Ransom\n\n");

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
    T = idata.N * idata.dt;
    if (strlen(remove_whitespace(idata.object)) > 0) {
        printf("Analyzing '%s' data from '%s'.\n\n",
               remove_whitespace(idata.object), cmd->argv[0]);
    } else {
        printf("Analyzing data from '%s'.\n\n", cmd->argv[0]);
    }
    min_orb_p = MINORBP;
    if (cmd->noaliasP)
        max_orb_p = T / 2.0;
    else
        max_orb_p = T / 1.2;

    /* open the FFT file and get its length */

    fftfile = chkfopen(cmd->argv[0], "rb");
    nbins = chkfilelen(fftfile, sizeof(fcomplex));

    /* Check that cmd->maxfft is an acceptable power of 2 */

    ct = 4;
    ii = 1;
    while (ct < MAXREALFFT || ii) {
        if (ct == cmd->maxfft)
            ii = 0;
        ct <<= 1;
    }
    if (ii) {
        printf("\n'maxfft' is out of range or not a power-of-2.\n\n");
        exit(1);
    }

    /* Check that cmd->minfft is an acceptable power of 2 */

    ct = 4;
    ii = 1;
    while (ct < MAXREALFFT || ii) {
        if (ct == cmd->minfft)
            ii = 0;
        ct <<= 1;
    }
    if (ii) {
        printf("\n'minfft' is out of range or not a power-of-2.\n\n");
        exit(1);
    }

    /* Low and high Fourier freqs to check */

    if (cmd->floP) {
        cmd->rlo = floor(cmd->flo * T);
        if (cmd->rlo < cmd->lobin)
            cmd->rlo = cmd->lobin;
        if (cmd->rlo > cmd->lobin + nbins - 1) {
            printf("\nLow frequency to search 'flo' is greater than\n");
            printf("   the highest available frequency.  Exiting.\n\n");
            exit(1);
        }
    } else {
        cmd->rlo = 1.0;
        if (cmd->rlo < cmd->lobin)
            cmd->rlo = cmd->lobin;
        if (cmd->rlo > cmd->lobin + nbins - 1) {
            printf("\nLow frequency to search 'rlo' is greater than\n");
            printf("   the available number of points.  Exiting.\n\n");
            exit(1);
        }
    }
    if (cmd->fhiP) {
        cmd->rhi = ceil(cmd->fhi * T);
        if (cmd->rhi > cmd->lobin + nbins - 1)
            cmd->rhi = cmd->lobin + nbins - 1;
        if (cmd->rhi < cmd->rlo) {
            printf("\nHigh frequency to search 'fhi' is less than\n");
            printf("   the lowest frequency to search 'flo'.  Exiting.\n\n");
            exit(1);
        }
    } else if (cmd->rhiP) {
        if (cmd->rhi > cmd->lobin + nbins - 1)
            cmd->rhi = cmd->lobin + nbins - 1;
        if (cmd->rhi < cmd->rlo) {
            printf("\nHigh frequency to search 'rhi' is less than\n");
            printf("   the lowest frequency to search 'rlo'.  Exiting.\n\n");
            exit(1);
        }
    }

    /* Determine how many different mini-fft sizes we will use */

    nfftsizes = 1;
    ii = cmd->maxfft;
    while (ii > cmd->minfft) {
        ii >>= 1;
        nfftsizes++;
    }

    /* Allocate some memory and prep some variables.             */
    /* For numtoread, the 6 just lets us read extra data at once */

    numtoread = 6 * cmd->maxfft;
    if (cmd->stack == 0)
        powers = gen_fvect(numtoread);
    minifft = (float *) fftwf_malloc(sizeof(float) *
                                     (cmd->maxfft * cmd->numbetween + 2));
    ncand2 = 2 * cmd->ncand;
    list = (rawbincand *) malloc(sizeof(rawbincand) * ncand2);
    for (ii = 0; ii < ncand2; ii++)
        list[ii].mini_sigma = 0.0;
    for (ii = 0; ii < MININCANDS; ii++)
        tmplist[ii].mini_sigma = 0.0;
    filepos = cmd->rlo - cmd->lobin;
    numchunks = (float) (cmd->rhi - cmd->rlo) / numtoread;
    printf("Searching...\n");
    printf("   Amount complete = %3d%%", 0);
    fflush(stdout);

    /* Prep FFTW */
    read_wisdom();

    /* Loop through fftfile */

    while ((filepos + cmd->lobin) < cmd->rhi) {

        /* Calculate percentage complete */

        newper = (int) (loopct / numchunks * 100.0);

        if (newper > oldper) {
            newper = (newper > 99) ? 100 : newper;
            printf("\r   Amount complete = %3d%%", newper);
            oldper = newper;
            fflush(stdout);
        }

        /* Adjust our search parameters if close to end of zone to search */

        binsleft = cmd->rhi - (filepos + cmd->lobin);
        if (binsleft < cmd->minfft)
            break;
        if (binsleft < numtoread) {     /* Change numtoread */
            numtoread = cmd->maxfft;
            while (binsleft < numtoread) {
                cmd->maxfft /= 2;
                numtoread = cmd->maxfft;
            }
        }
        fftlen = cmd->maxfft;

        /* Read from fftfile */

        if (cmd->stack == 0) {
            data = read_fcomplex_file(fftfile, filepos, numtoread);
            for (ii = 0; ii < numtoread; ii++)
                powers[ii] = POWER(data[ii].r, data[ii].i);
            numsumpow = 1;
        } else {
            powers = read_float_file(fftfile, filepos, numtoread);
            numsumpow = cmd->stack;
        }
        if (filepos == 0)
            powers[0] = 1.0;

        /* Chop the powers that are way above the median level */

        prune_powers(powers, numtoread, numsumpow);

        /* Loop through the different small FFT sizes */

        while (fftlen >= cmd->minfft) {

            halffftlen = fftlen / 2;
            powers_pos = powers;
            powers_offset = 0;

            /* Create the appropriate FFT plan */

            fftplan = fftwf_plan_dft_r2c_1d(cmd->interbinP ? fftlen : 2 * fftlen,
                                            minifft, (fftwf_complex *) minifft,
                                            FFTW_PATIENT);

            /* Perform miniffts at each section of the powers array */

            while ((numtoread - powers_offset) >
                   (int) ((1.0 - cmd->overlap) * cmd->maxfft + DBLCORRECT)) {

                /* Copy the proper amount and portion of powers into minifft */

                memcpy(minifft, powers_pos, fftlen * sizeof(float));
                /* For Fourier interpolation use a zeropadded FFT */
                if (cmd->numbetween > 1 && !cmd->interbinP) {
                    for (ii = fftlen; ii < cmd->numbetween * fftlen; ii++)
                        minifft[ii] = 0.0;
                }

                /* Perform the minifft */

                fftwf_execute(fftplan);

                /* Normalize and search the miniFFT */

                norm = sqrt(fftlen * numsumpow) / minifft[0];
                for (ii = 0; ii < (cmd->interbinP ? fftlen + 1 : 2 * fftlen + 1);
                     ii++)
                    minifft[ii] *= norm;
                search_minifft((fcomplex *) minifft, halffftlen, min_orb_p,
                               max_orb_p, tmplist, MININCANDS, cmd->harmsum,
                               cmd->numbetween, idata.N, T,
                               (double) (powers_offset + filepos + cmd->lobin),
                               cmd->interbinP ? INTERBIN : INTERPOLATE,
                               cmd->noaliasP ? NO_CHECK_ALIASED : CHECK_ALIASED);

                /* Check if the new cands should go into the master cand list */

                for (ii = 0; ii < MININCANDS; ii++) {
                    if (tmplist[ii].mini_sigma > minsig) {

                        /* Check to see if another candidate with these properties */
                        /* is already in the list.                                 */

                        if (not_already_there_rawbin(tmplist[ii], list, ncand2)) {
                            list[ncand2 - 1] = tmplist[ii];
                            minsig = percolate_rawbincands(list, ncand2);
                        }
                    } else {
                        break;
                    }
                    /* Mini-fft search for loop */
                }

                totnumsearched += fftlen;
                powers_pos += (int) (cmd->overlap * fftlen);
                powers_offset = powers_pos - powers;

                /* Position of mini-fft in data set while loop */
            }

            fftwf_destroy_plan(fftplan);
            fftlen >>= 1;

            /* Size of mini-fft while loop */
        }

        if (cmd->stack == 0)
            vect_free(data);
        else
            vect_free(powers);
        filepos += (numtoread - (int) ((1.0 - cmd->overlap) * cmd->maxfft));
        loopct++;

        /* File position while loop */
    }

    /* Print the final percentage update */

    printf("\r   Amount complete = %3d%%\n\n", 100);

    /* Print the number of frequencies searched */

    printf("Searched %.0f pts (including interbins).\n\n", totnumsearched);

    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n",
           ttim, utim, stim);
    printf("  Total time: %.3f sec\n\n", tott);

    printf("Writing result files and cleaning up.\n");

    /* Count how many candidates we actually have */

    ii = 0;
    while (ii < ncand2 && list[ii].mini_sigma != 0)
        ii++;
    newncand = (ii > cmd->ncand) ? cmd->ncand : ii;

    /* Set our candidate notes to all spaces */

    notes = malloc(sizeof(char) * newncand * 18 + 1);
    for (ii = 0; ii < newncand; ii++)
        strncpy(notes + ii * 18, "                     ", 18);

    /* Check the database for possible known PSR detections */

    if (idata.ra_h && idata.dec_d) {
        for (ii = 0; ii < newncand; ii++) {
            comp_rawbin_to_cand(&list[ii], &idata, notes + ii * 18, 0);
        }
    }

    /* Compare the candidates with each other */

    compare_rawbin_cands(list, newncand, notes);

    /* Send the candidates to the text file */

    file_rawbin_candidates(list, notes, newncand, cmd->harmsum, rootfilenm);

    /* Write the binary candidate file */
    {
        char *candnm;

        candnm = (char *) calloc(strlen(rootfilenm) + 15, sizeof(char));
        sprintf(candnm, "%s_bin%d.cand", rootfilenm, cmd->harmsum);
        candfile = chkfopen(candnm, "wb");
        chkfwrite(list, sizeof(rawbincand), (unsigned long) newncand, candfile);
        fclose(candfile);
        free(candnm);
    }

    /* Free our arrays and close our files */

    if (cmd->stack == 0)
        vect_free(powers);
    free(list);
    fftwf_free(minifft);
    free(notes);
    free(rootfilenm);
    fclose(fftfile);
    printf("Done.\n\n");
    return (0);
}
