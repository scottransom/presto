/*     Real-Valued Data FFT Program        */
/*          by Scott Ransom                */
/*            Version 3.0                  */

#include <time.h>
#include <sys/times.h>
#include "clk_tck.h"
#include "misc_utils.h"
#include "chkio.h"
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
    float *data;
    int isign = -1, numfiles;
    long long numdata = 0, filelen, maxfilelen = 0;
    struct tms runtimes;
    double ttim, stim, utim, tott;
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

    tott = times(&runtimes) / (double) CLK_TCK;
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
            if (cmd->tmpdirP) {
                filenmlen = strlen(cmd->tmpdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", cmd->tmpdir, root, tmpsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", datdir, root, tmpsuffix);
            }
            if (cmd->outdirP) {
                filenmlen = strlen(cmd->outdir) + 1 + strlen(root) + 5;
                outfilenm = (char *) calloc(filenmlen, 1);
                sprintf(outfilenm, "%s/%s.%s", cmd->outdir, root, outsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                outfilenm = (char *) calloc(filenmlen, 1);
                sprintf(outfilenm, "%s/%s.%s", datdir, root, outsuffix);
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

            /* Remove the scratch files */

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
            printf("   Writing.\n");
            chkfwrite(data, sizeof(float), numdata, outfile);
            fclose(outfile);

            /* Delete the input files if requested */

            if (cmd->deleteP)
                remove(datfilenm);
        }

        /* Close our input files */

        fclose(datfile);
    }

    /* Output the timing information */

    printf("Finished.\n\n");
    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
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
    exit(0);
}
