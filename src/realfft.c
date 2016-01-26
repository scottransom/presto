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
    multifile *datfile, *tmpfile = NULL, *outfile;
    char *datdir, **datfilenms, **tmpfilenms, **outfilenms;
    float *data;
    int ii, isign = -1, numfiles;
    long long numdata = 0, maxfilelen = 0;
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
    printf("        by Scott M. Ransom\n");
    printf("            3 Oct, 2000\n\n");

    /* Get our file information */

    numfiles = cmd->argc;
    printf("Checking data in %d input file(s):\n", numfiles);
    datfilenms = (char **) malloc(numfiles * sizeof(char *));
    tmpfilenms = (char **) malloc(numfiles * sizeof(char *));
    outfilenms = (char **) malloc(numfiles * sizeof(char *));
    {
        int hassuffix = 0, new_hassuffix = 0, filenmlen;
        char *dir, *filenm, *root, *suffix;

        for (ii = 0; ii < numfiles; ii++) {
            if (ii == 0) {
                split_path_file(cmd->argv[0], &datdir, &filenm);
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
            } else {
                split_path_file(cmd->argv[ii], &dir, &filenm);
                new_hassuffix = split_root_suffix(filenm, &root, &suffix);
                if (new_hassuffix && hassuffix) {
                    if (strcmp(datsuffix, suffix)) {
                        printf("\nAll input files must have the same suffix!\n\n");
                        exit(1);
                    }
                }
                if (strcmp(datdir, dir)) {
                    printf("\nAll input files must be in the same directory!\n\n");
                    exit(1);
                }
                if (new_hassuffix)
                    free(suffix);
                free(dir);
                free(filenm);
            }
            filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
            datfilenms[ii] = (char *) calloc(filenmlen, 1);
            sprintf(datfilenms[ii], "%s/%s.%s", datdir, root, datsuffix);
            if (cmd->tmpdirP) {
                filenmlen = strlen(cmd->tmpdir) + 1 + strlen(root) + 5;
                tmpfilenms[ii] = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenms[ii], "%s/%s.%s", cmd->tmpdir, root, tmpsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                tmpfilenms[ii] = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenms[ii], "%s/%s.%s", datdir, root, tmpsuffix);
            }
            if (cmd->outdirP) {
                filenmlen = strlen(cmd->outdir) + 1 + strlen(root) + 5;
                outfilenms[ii] = (char *) calloc(filenmlen, 1);
                sprintf(outfilenms[ii], "%s/%s.%s", cmd->outdir, root, outsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                outfilenms[ii] = (char *) calloc(filenmlen, 1);
                sprintf(outfilenms[ii], "%s/%s.%s", datdir, root, outsuffix);
            }
            free(root);
        }
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

    datfile = fopen_multifile(numfiles, datfilenms, "r", 0);
    for (ii = 0; ii < numfiles; ii++) {
        if (datfile->filelens[ii] > maxfilelen)
            maxfilelen = datfile->filelens[ii];
    }
    for (ii = 0; ii < numfiles; ii++)
        printf("   %d:  '%s'\n", ii + 1, datfile->filenames[ii]);
    numdata = datfile->length / sizeof(float);
    if (isign == -1) {
        if (datfile->length % sizeof(float)) {
            printf("\nInput file does not contain the correct number of\n");
            printf("   bytes for it to be floating point data!\n\n");
            exit(1);
        }
        printf("\nData OK.  There are %lld floats.\n\n", numdata);
    } else {
        if (datfile->length % sizeof(fcomplex)) {
            printf("\nInput file does not contain the correct number of\n");
            printf("   bytes for it to be single precision complex data!\n\n");
            exit(1);
        }
        printf("\nData OK.  There are %lld complex points.\n\n", numdata / 2);
    }
    printf("Result will be written to %d output file(s):\n", numfiles);
    for (ii = 0; ii < numfiles; ii++)
        printf("   %d:  '%s'\n", ii + 1, outfilenms[ii]);

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
            int status, maxslen = 0, slen, suf;
            char *cmd, *root, *suffix;

            for (ii = 0; ii < datfile->numfiles; ii++) {
                slen = strlen(datfile->filenames[ii]);
                if (slen > maxslen)
                    maxslen = slen;
            }
            maxslen = 2 * maxslen + 10;
            cmd = (char *) calloc(maxslen, 1);
            for (ii = 0; ii < datfile->numfiles; ii++) {
                suf = split_root_suffix(datfile->filenames[ii], &root, &suffix);
                sprintf(cmd, "cp %s %s.bak", datfile->filenames[ii], root);
                if ((status = (system(cmd))) == -1 || status == 127) {
                    perror("\nSystem call (cp) failed");
                    printf("\n");
                    exit(1);
                }
                if (suf)
                    free(suffix);
                free(root);
            }
            free(cmd);
        }

        /* Close the input files and re-open them in write mode */

        fclose_multifile(datfile);
        datfile = fopen_multifile(numfiles, datfilenms, "r+", maxfilelen);
        tmpfile = fopen_multifile(numfiles, tmpfilenms, "w+", maxfilelen);
        if (isign == 1) {
            realfft_scratch_inv(datfile, tmpfile, numdata);
        } else {
            realfft_scratch_fwd(datfile, tmpfile, numdata);
        }
        /* twopassfft_scratch(datfile, tmpfile, numdata, isign); */

        /* Remove the scratch files */

        fclose_multifile(tmpfile);
        for (ii = 0; ii < numfiles; ii++)
            remove(tmpfilenms[ii]);

        /* Change the output filename to the correct suffix and   */
        /* rename the back-up data files if needed.               */

        {
            int maxslen = 0, slen, suf;
            char *file1, *file2, *root, *suffix;

            for (ii = 0; ii < datfile->numfiles; ii++) {
                slen = strlen(datfile->filenames[ii]);
                if (slen > maxslen)
                    maxslen = slen;
            }
            maxslen = maxslen + 5;
            file1 = (char *) calloc(maxslen, 1);
            file2 = (char *) calloc(maxslen, 1);
            for (ii = 0; ii < datfile->numfiles; ii++) {
                suf = split_root_suffix(datfile->filenames[ii], &root, &suffix);
                sprintf(file1, "%s.%s", root, datsuffix);
                sprintf(file2, "%s.%s", root, outsuffix);
                rename(file1, file2);
                if (!cmd->deleteP) {
                    sprintf(file1, "%s.bak", root);
                    rename(file1, datfile->filenames[ii]);
                }
                if (suf)
                    free(suffix);
                free(root);
            }
            free(file1);
            free(file2);
        }

    } else {

        /* Perform standard FFT for real functions  */

        outfile = fopen_multifile(numfiles, outfilenms, "w", maxfilelen);
        if (isign == -1) {
            printf("\nPerforming in-core forward FFT on data:\n");
        } else {
            printf("\nPerforming in-core inverse FFT on data:\n");
        }
        printf("   Reading.\n");
        data = gen_fvect(numdata);
        fread_multifile(data, sizeof(float), numdata, datfile);
        printf("   Transforming.\n");
        realfft(data, numdata, isign);
        /* fftwcall((fcomplex *)data, numdata/2, isign); */
        /* tablesixstepfft((fcomplex *)data, numdata/2, isign); */
        printf("   Writing.\n");
        fwrite_multifile(data, sizeof(float), numdata, outfile);
        fclose_multifile(outfile);

        /* Delete the input files if requested */

        if (cmd->deleteP) {
            for (ii = 0; ii < numfiles; ii++)
                remove(datfilenms[ii]);
        }
    }

    /* Close our input files */

    fclose_multifile(datfile);

    /* Output the timing information */

    printf("Finished.\n\n");
    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("  CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n",
           ttim, utim, stim);
    printf("  Total time elapsed:  %.3f sec\n\n", tott);

    /*
       fftw_print_max_memory_usage();
       fftw_check_memory_leaks();
     */

    for (ii = 0; ii < numfiles; ii++) {
        free(datfilenms[ii]);
        free(tmpfilenms[ii]);
        free(outfilenms[ii]);
    }
    free(datfilenms);
    free(tmpfilenms);
    free(outfilenms);
    exit(0);
}
