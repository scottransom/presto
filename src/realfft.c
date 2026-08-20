/*     Real-Valued Data FFT Program        */
/*          by Scott Ransom                */
/*            Version 3.0                  */

#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include "misc_utils.h"
#include "chkio.h"
#include "makeinf.h"
#include "ransomfft.h"
#include "vectors.h"
#include "realfft_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/*  This program calculates the FFT of a file containing    */
/*  a number of single-precision floats representing        */
/*  real numbers.  (i.e. a normal time series)              */
/*  The data is assumed to be located in the directory,     */
/*  specified in the input filename(s).                     */
/*  Input filename(s) must include '.dat' or '.fft'         */
/*        suffixes.  The output file(s) will have the       */
/*        appropriate other suffix.                         */
/*  Do not end paths in '/'.                                */
/*  Scratch file(s) are the same size as the input file(s). */
/*  If '-inv' is specified, the file to be transformed      */
/*        should end in '.fft'.  Otherwise, it should end   */
/*        in '.dat'.                                        */

int main(int argc, char *argv[])
{
    FILE *datfile, *tmpfile = NULL, *outfile;
    char *datdir = NULL, *datfilenm = NULL, *tmpfilenm = NULL, *outfilenm = NULL;
    char *inffilenm = NULL, *outroot = NULL;
    float *data;
    int isign = -1, numfiles;
    long long numdata = 0, filelen, maxfilelen = 0;
    double T = 0.0;
    infodata idata;
    struct tms runtimes;
    double ttim, stim, utim, tott;
    double clk_tck = (double) sysconf(_SC_CLK_TCK);
    Cmdline *cmd;
    char datsuffix[] = "dat";
    char outsuffix[] = "fft";
    char tmpsuffix[] = "tmp";

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);

#ifdef DEBUG
    showOptionValues();
#endif

    tott = times(&runtimes) / clk_tck;
    printf("\n");
    printf("   Real-Valued Data FFT Program v3.0\n");
    printf("        by Scott M. Ransom\n\n");

    /* Get our file information */

    numfiles = cmd->argc;
    for (int fi = 0 ; fi < numfiles ; fi++) {
        printf("%4d:  Processing data in '%s'\n", fi + 1, cmd->argv[fi]);
        {
            int hassuffix = 0, filenmlen;
            char *filenm, *root, *suffix;

            /* Release the names built for the previous file, if any */
            free(datfilenm);
            free(tmpfilenm);
            free(outfilenm);
            free(inffilenm);
            free(outroot);
            split_path_file(cmd->argv[fi], &datdir, &filenm);
            hassuffix = split_root_suffix(filenm, &root, &suffix);
            if (hassuffix) {
                if (strcmp(suffix, "fft") == 0) {
                    isign = 1;
                    strcpy(datsuffix, "fft");
                    strcpy(outsuffix, "dat");
                }
                free(suffix);
            }
            free(filenm);
            filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
            datfilenm = (char *) calloc(filenmlen, 1);
            sprintf(datfilenm, "%s/%s.%s", datdir, root, datsuffix);
            filenmlen = strlen(datdir) + 1 + strlen(root) + 1;
            inffilenm = (char *) calloc(filenmlen, 1);
            sprintf(inffilenm, "%s/%s", datdir, root);
            if (cmd->tmpdirP) {
                filenmlen = strlen(cmd->tmpdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", cmd->tmpdir, root, tmpsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", datdir, root, tmpsuffix);
            }
            {
                /* Red noise removal appends '_red' to the output root, */
                /* exactly as rednoise(1) does.                         */
                char *redtag = cmd->rednoiseP ? "_red" : "";

                if (cmd->outdirP) {
                    filenmlen = strlen(cmd->outdir) + 1 + strlen(root) +
                        strlen(redtag) + 1;
                    outroot = (char *) calloc(filenmlen, 1);
                    sprintf(outroot, "%s/%s%s", cmd->outdir, root, redtag);
                } else {
                    /* Keep the path the user gave us, rather than the */
                    /* absolute one, so that the names stay short.     */
                    char *argroot, *argsuffix;

                    if (split_root_suffix(cmd->argv[fi], &argroot, &argsuffix))
                        free(argsuffix);
                    outroot = (char *) calloc(strlen(argroot) + strlen(redtag) + 1, 1);
                    sprintf(outroot, "%s%s", argroot, redtag);
                    free(argroot);
                }
                outfilenm = (char *) calloc(strlen(outroot) + 5, 1);
                sprintf(outfilenm, "%s.%s", outroot, outsuffix);
            }
            free(root);
            free(datdir);
        }

        /* Force a forward or inverse transform.   */
        /* Note that the suffixes do _not_ change! */

        if (cmd->forwardP)
            isign = -1;
        if (cmd->inverseP)
            isign = 1;
        if (cmd->diskfftP && cmd->memfftP) {
            printf("\nYou cannot take both an in- and out-of-core FFT!\n\n");
            exit(1);
        }

        /* Red noise removal only makes sense for a forward transform, and */
        /* it needs the '.inf' file in order to know how long the          */
        /* observation was.                                                */

        if (cmd->rednoiseP) {
            char *tmpnm;
            size_t tmpnmlen;

            if (isign == 1) {
                printf("\nYou cannot remove red noise from an inverse FFT!\n\n");
                exit(1);
            }
            tmpnmlen = strlen(inffilenm);
            tmpnm = (char *) calloc(tmpnmlen + 5, 1);
            memcpy(tmpnm, inffilenm, tmpnmlen);
            memcpy(tmpnm + tmpnmlen, ".inf", 4);
            if (access(tmpnm, R_OK)) {
                printf("\nCannot read '%s', which is required by '-rednoise'!\n\n",
                       tmpnm);
                exit(1);
            }
            free(tmpnm);
            readinf(&idata, inffilenm);
            T = idata.N * idata.dt;
            if (strlen(outroot) >= sizeof(idata.name)) {
                printf("\nOutput file name '%s' is too long for a '.inf' file!\n\n",
                       outroot);
                exit(1);
            }
            strcpy(idata.name, outroot);
        }

        /* Open and check data files */

        datfile = chkfopen(datfilenm, "r");
        filelen = chkfilelen(datfile, 1);
        if (filelen > maxfilelen)
            maxfilelen = filelen;
        numdata = filelen / sizeof(float);
        if (isign == -1) {
            if (filelen % sizeof(float)) {
                printf("\nInput file does not contain the correct number of\n");
                printf("   bytes for it to be floating point data!  Skipping.\n\n");
                continue;
            }
            printf("\nData OK.  There are %lld floats.\n\n", numdata);
        } else {
            if (filelen % sizeof(fcomplex)) {
                printf("\nInput file does not contain the correct number of\n");
                printf("   bytes for it to be single precision complex data!  Skipping.\n\n");
                continue;
            }
            printf("\nData OK.  There are %lld complex points.\n\n", numdata / 2);
        }
        printf("%4d:   Result will be in '%s'\n", fi + 1, outfilenm);

        /*  Start the transform sequence  */

        if ((numdata > MAXREALFFT || cmd->diskfftP) && !cmd->memfftP) {

            /*  Perform Two-Pass, Out-of-Core, FFT  */

            if (isign == -1) {
                printf("\nPerforming out-of-core two-pass forward FFT on data.\n");
            } else {
                printf("\nPerforming out-of-core two-pass inverse FFT on data.\n");
            }
            if (cmd->outdirP)
                printf("\nWarning:  '-outdir' is not used for out-of-core FFTs.\n"
                       "          The results will be next to the input file.\n");

            /* Copy the input files if we want to keep them */

            if (!cmd->deleteP) {
                int status, slen, suf;
                char *cmd, *root, *suffix;

                slen = strlen(datfilenm);
                cmd = (char *) calloc(2 * slen + 10, 1);
                suf = split_root_suffix(datfilenm, &root, &suffix);
                sprintf(cmd, "cp %s %s.bak", datfilenm, root);
                if ((status = (system(cmd))) == -1 || status == 127) {
                    perror("\nSystem call (cp) failed");
                    printf("\n");
                    exit(1);
                }
                if (suf)
                    free(suffix);
                free(root);
                free(cmd);
            }

            /* Close the input files and re-open them in write mode */

            fclose(datfile);
            datfile = chkfopen(datfilenm, "r+");
            tmpfile = chkfopen(tmpfilenm, "w+");
            if (isign == 1) {
                realfft_scratch_inv(datfile, tmpfile, numdata);
            } else {
                realfft_scratch_fwd(datfile, tmpfile, numdata);
            }

            /* Remove the scratch files.  Note that datfile holds the    */
            /* transformed data at this point, so close it (and flush it) */
            /* before anything else touches the file.                     */

            fclose(datfile);
            datfile = NULL;
            fclose(tmpfile);
            remove(tmpfilenm);

            /* Change the output filename to the correct suffix and   */
            /* rename the back-up data files if needed.               */

            {
                int slen, suf;
                char *file1, *file2, *root, *suffix;

                slen = strlen(datfilenm);
                file1 = (char *) calloc(slen + 5, 1);
                file2 = (char *) calloc(slen + 5, 1);
                suf = split_root_suffix(datfilenm, &root, &suffix);
                sprintf(file1, "%s.%s", root, datsuffix);
                sprintf(file2, "%s.%s", root, outsuffix);
                rename(file1, file2);
                if (!cmd->deleteP) {
                    sprintf(file1, "%s.bak", root);
                    rename(file1, datfilenm);
                }

                /* Deredden the transform in place, then rename it '_red' */

                if (cmd->rednoiseP) {
                    FILE *redfile;
                    long long numbins;
                    char *redfilenm;

                    printf("\nRemoving red noise.\n");
                    redfilenm = (char *) calloc(slen + 10, 1);
                    sprintf(redfilenm, "%s_red.%s", root, outsuffix);
                    redfile = chkfopen(file2, "rb+");
                    numbins = chkfilelen(redfile, sizeof(fcomplex));
                    if (deredden_file(redfile, redfile, numbins, cmd->startwidth,
                                      cmd->endwidth, cmd->endfreq, T) < 0) {
                        printf("\nRed noise removal failed.\n\n");
                        exit(1);
                    }
                    fclose(redfile);
                    if (rename(file2, redfilenm)) {
                        perror("\nError renaming the dereddened file");
                        printf("\n");
                        exit(1);
                    }
                    if (snprintf(idata.name, sizeof(idata.name), "%s_red", root) >=
                        (int) sizeof(idata.name)) {
                        printf("\nOutput file name is too long for a '.inf' file!\n\n");
                        exit(1);
                    }
                    writeinf(&idata);
                    printf("\nResult is in '%s'\n", redfilenm);
                    free(redfilenm);
                }
                if (suf)
                    free(suffix);
                free(root);
                free(file1);
                free(file2);
            }

        } else {

            /* Perform standard FFT for real functions  */

            outfile = chkfopen(outfilenm, "w");
            if (isign == -1) {
                printf("\nPerforming in-core forward FFT on data:\n");
            } else {
                printf("\nPerforming in-core inverse FFT on data:\n");
            }
            printf("   Reading.\n");
            data = gen_fvect(numdata);
            chkfread(data, sizeof(float), numdata, datfile);
            printf("   Transforming.\n");
            realfft(data, numdata, isign);
            /* fftwcall((fcomplex *)data, numdata/2, isign); */
            /* tablesixstepfft((fcomplex *)data, numdata/2, isign); */
            if (cmd->rednoiseP) {
                printf("   Removing red noise.\n");
                if (deredden((fcomplex *) data, numdata / 2, cmd->startwidth,
                             cmd->endwidth, cmd->endfreq, T) < 0) {
                    printf("\nRed noise removal failed.\n\n");
                    exit(1);
                }
            }
            printf("   Writing.\n");
            chkfwrite(data, sizeof(float), numdata, outfile);
            fclose(outfile);
            vect_free(data);
            if (cmd->rednoiseP)
                writeinf(&idata);

            /* Delete the input files if requested */

            if (cmd->deleteP)
                remove(datfilenm);
        }

        /* Close our input files */

        if (datfile != NULL)
            fclose(datfile);
    }

    /* Output the timing information */

    printf("Finished.\n\n");
    printf("Timing summary:\n");
    tott = times(&runtimes) / clk_tck - tott;
    utim = runtimes.tms_utime / clk_tck;
    stim = runtimes.tms_stime / clk_tck;
    ttim = utim + stim;
    printf("  CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n",
           ttim, utim, stim);
    printf("  Total time elapsed:  %.3f sec (%.3f sec/file)\n\n", tott, tott / numfiles);

    /*
       fftw_print_max_memory_usage();
       fftw_check_memory_leaks();
     */

    free(datfilenm);
    free(tmpfilenm);
    free(outfilenm);
    free(inffilenm);
    free(outroot);
    exit(0);
}
