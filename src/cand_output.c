#define COPYSIGN(a,b) ((b) < 0.0 ? (-fabs(a)) : (fabs(a)))

/*  NOTE:  See presto.h for function definitions. */

/*  NEED TO ADD A METHOD TO TELL WHERE a2x RESIDES  */

#include "presto.h"

void file_reg_candidates(fourierprops cand[], char *notes, int numcands,
                         double dt, long N, double nph, char name[], char longname[])
{
    FILE *fname;
    double T, pownph, pownpherr;
    int i, j, k = 0;
    int nlines = 77, pages, extralines, nwrit, linestoprint;
    char output[40], infonm[100], command[200], tmp[] = "";
    rzwerrs rzws;

    sprintf(infonm, "%s.inf", name);
    fname = chkfopen(longname, "w");

    if (numcands <= 0) {
        printf(" Must have at least 1 candidate in ");
        printf("file_reg_candidates().\n\n");
        exit(1);
    }
    pages = numcands / nlines + 1;
    extralines = numcands % nlines;

    T = dt * N;

    for (i = 1; i <= pages; i++) {

        /*              1234xx1234567xx123456789012xx123456789012345678xx */
        fprintf(fname, " Cand   Sigma     Power/      Candidate Period   ");
        /*              123456789012345678xx1234567xx123456789012xx */
        fprintf(fname, "Period Derivative    Phase     Centroid   ");
        /*              123456789012xx123456789012345678xx123456789012345xx */
        fprintf(fname, "  Purity      FFT Frequency     FFT f-dot     ");
        /*              12345678901234xx12345678901xx1234568901234567890 */
        fprintf(fname, " FFT f-dotdot       Power/            \n");

        /*              1234xx1234567xx123456789012xx123456789012345678xx */
        fprintf(fname, "Number  Local  Local Power       (seconds)       ");
        /*              123456789012345678xx1234567xx123456789012xx */
        fprintf(fname, "    (sec/sec)        (rad)     (0.0-1.0)  ");
        /*              123456789012xx123456789012345678xx123456789012345xx */
        fprintf(fname, " <p> = 1.0     (bin-number)     (bins/FFT)    ");
        /*              12345678901234xx12345678901xx1234568901234567890 */
        fprintf(fname, " (bins/FFT^2)      Bin Zero         Notes:\n");

        /*              1234xx1234567xx123456789012xx123456789012345678xx */
        fprintf(fname, "-------------------------------------------------");
        /*              123456789012345678xx1234567xx123456789012xx */
        fprintf(fname, "------------------------------------------");
        /*              123456789012xx123456789012345678xx123456789012345xx */
        fprintf(fname, "--------------------------------------------------");
        /*              12345678901234xx12345678901xx1234568901234567890 */
        fprintf(fname, "-------------------------------------------------\n");

        if (i == pages) {
            linestoprint = extralines;
        } else {
            linestoprint = nlines;
        }

        for (j = 0; j < linestoprint; j++, k++) {

            pownph = cand[k].rawpow / nph;
            pownpherr = sqrt(2.0 * pownph);
            calc_rzwerrs(&cand[k], T, &rzws);

            /*  Now output it... */

            fprintf(fname, " %-5d", k + 1);
            fprintf(fname, "%7.2f  ", cand[k].sig);
            nwrit = nice_output_1(output, cand[k].pow, cand[k].powerr, 0);
            nwrit += (12 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 12 - nwrit, tmp);
            nwrit = nice_output_1(output, rzws.p, rzws.perr, 0);
            nwrit += (18 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 18 - nwrit, tmp);
            nwrit = nice_output_1(output, rzws.pd, rzws.pderr, 0);
            nwrit += (18 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 18 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].phs, cand[k].phserr, 0);
            nwrit += (7 - nwrit) / 2;
            fprintf(fname, "%*s%*s ", nwrit, output, 7 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].cen, cand[k].cenerr, 0);
            nwrit += (12 - nwrit) / 2;
            fprintf(fname, "  %*s%*s ", nwrit, output, 11 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].pur, cand[k].purerr, 0);
            nwrit += (11 - nwrit) / 2;
            fprintf(fname, "%*s%*s", nwrit, output, 12 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].r, cand[k].rerr, 0);
            nwrit += (18 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 18 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].z, cand[k].zerr, 0);
            nwrit += (12 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 12 - nwrit, tmp);
            nwrit = nice_output_1(output, cand[k].w, cand[k].werr, 0);
            nwrit += (14 - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, 14 - nwrit, tmp);
            if (pownph > 1.0e7) {
                fprintf(fname, "%12.5g    ", pownph);
            } else {
                nwrit = nice_output_1(output, pownph, pownpherr, 0);
                nwrit += (16 - nwrit) / 2;
                fprintf(fname, "%*s%*s", nwrit, output, 16 - nwrit, tmp);
            }
            fprintf(fname, "  %.20s\n", notes + k * 20);
            fflush(fname);
        }
    }
    fprintf(fname, "\n Notes:  H = Harmonic.  SL = Possible sidelobe.\n\n");
    fclose(fname);
    sprintf(command, "cat %s >> %s", infonm, longname);
    system(command);
    sprintf(command,
            "$PRESTO/bin/a2x -l -c1 -n80 -title -date -num %s > %s.ps",
            longname, longname);
    system(command);
}


void file_bin_candidates(binaryprops cand[], char *notes, int numcands, char name[])
/* Outputs a .ps file describing all the binary candidates from a    */
/*   binary search. */
{
    FILE *fname;
    int i, j, k = 0, wide;
    int nlines = 77, pages, extralines, nwrit, linestoprint;
    char output[40], filenm[100], infonm[100], command[200], tmp[] = "";

    sprintf(filenm, "%s_bin", name);
    sprintf(infonm, "%s.inf", name);
    fname = chkfopen(filenm, "w");

    if (numcands <= 0) {
        printf(" Must have at least 1 candidate in ");
        printf("file_bin_candidates().\n\n");
        exit(1);
    }
    pages = numcands / nlines + 1;
    extralines = numcands % nlines;

    for (i = 1; i <= pages; i++) {

        fprintf(fname, " Cand   MiniFFT        Normalized    Cand PSR Period  ");
        fprintf(fname, "  Bin Period     a*sin(i)/c       Phasemod     ");
        fprintf(fname, "Centroid    Purity     PSR FFT Freq     Bin FFT Freq ");
        fprintf(fname, " Detect FFT Freq\n");
        fprintf(fname, "Number   bins   Sigma    Power         (seconds)     ");
        fprintf(fname, "     (sec)         (lt-sec)      Amplitude 'z'  ");
        fprintf(fname, "(0.0-1.0)  <p> = 1.0   (bin-number)      (bin-number)");
        fprintf(fname, "  (bin-number)        Notes\n");
        fprintf(fname, "-------------------------------------------------------");
        fprintf(fname, "-----------------------------------------------");
        fprintf(fname, "--------------------------------------------------------");
        fprintf(fname, "---------------------------------\n");

        if (i == pages) {
            linestoprint = extralines;
        } else {
            linestoprint = nlines;
        }

        for (j = 0; j < linestoprint; j++, k++) {

            /*  Now output it... */

            fprintf(fname, " %-4d %7ld ", k + 1, cand[k].nfftbins);
            fprintf(fname, "%7.2f  ", cand[k].sig);

            nwrit = nice_output_1(output, cand[k].pow, cand[k].powerr, 0);
            wide = 12;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].ppsr, cand[k].ppsrerr, 0);
            wide = 15;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].pbin, cand[k].pbinerr, 0);
            wide = 14;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].asinic, cand[k].asinicerr, 0);
            wide = 12;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].z, cand[k].zerr, 0);
            wide = 14;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "  %*s%*s ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].cen, cand[k].cenerr, 0);
            wide = 11;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].pur, cand[k].purerr, 0);
            wide = 11;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].rpsr, cand[k].rpsrerr, 0);
            wide = 17;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].rbin, cand[k].rbinerr, 0);
            wide = 13;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s  ", nwrit, output, wide - nwrit, tmp);

            nwrit = nice_output_1(output, cand[k].rdetect, cand[k].rdetecterr, 0);
            wide = 14;
            nwrit += (wide - nwrit) / 2;
            fprintf(fname, "%*s%*s", nwrit, output, wide - nwrit, tmp);

            fprintf(fname, " %.18s\n", notes + k * 18);
            fflush(fname);
        }
    }
    fprintf(fname, "\n Notes:  MH = Modulation harmonic.  ");
    fprintf(fname, "H = Pulsar harmonic.  # indicates the candidate number.\n\n");
    fclose(fname);
    sprintf(command, "cat %s >> %s", infonm, filenm);
    system(command);
    sprintf(command,
            "$PRESTO/bin/a2x -l -c1 -n80 -title -date -num %s > %s.ps",
            filenm, filenm);
    system(command);
}

#undef COPYSIGN
