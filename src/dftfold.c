#include "dftfold.h"
#include "presto.h"
/* #include "plot2d.h" */
#include "dftfold_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
/* dftfold:  Does complex plane vector addition of a DFT freq */
/* Written by Scott Ransom on 31 Aug 00 based on Ransom and   */
/* Eikenberry paper I (to be completed sometime...).          */
{
    FILE *infile;
    char infilenm[200], outfilenm[200];
    int dataperread;
    unsigned long N;
    double T, rr = 0.0, norm = 1.0;
    dftvector dftvec;
    infodata idata;
    Cmdline *cmd;

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

    printf("\n\n");
    printf("        DFT Vector Folding Routine\n");
    printf("            by Scott M. Ransom\n");
    printf("              31 August, 2000\n\n");

    /* Open the datafile and read the info file */

    sprintf(infilenm, "%s.dat", cmd->argv[0]);
    infile = chkfopen(infilenm, "rb");
    readinf(&idata, cmd->argv[0]);

    /* The number of points in datafile */

    N = chkfilelen(infile, sizeof(float));
    dataperread = N / cmd->numvect;
/*   N = cmd->numvect * dataperread; */
    T = N * idata.dt;

    /* Calculate the Fourier frequency */

    if (!cmd->rrP) {
        if (cmd->ffP)
            rr = cmd->ff;
        else if (cmd->ppP)
            rr = T / cmd->pp;
        else {
            printf("\n  You must specify a frequency to fold!  Exiting.\n\n");
        }
    } else
        rr = cmd->rr;

    /* Calculate the amplitude normalization if required */

    if (cmd->normP)
        norm = 1.0 / sqrt(cmd->norm);
    else if (cmd->fftnormP) {
        FILE *fftfile;
        int kern_half_width, fftdatalen, startbin;
        double rrfrac, rrint;
        char fftfilenm[200];
        fcomplex *fftdata;

        sprintf(fftfilenm, "%s.fft", cmd->argv[0]);
        fftfile = chkfopen(fftfilenm, "rb");
        kern_half_width = r_resp_halfwidth(HIGHACC);
        fftdatalen = 2 * kern_half_width + 10;
        rrfrac = modf(rr, &rrint);
        startbin = (int) rrint - fftdatalen / 2;
        fftdata = read_fcomplex_file(fftfile, startbin, fftdatalen);
        norm = 1.0 / sqrt(get_localpower3d(fftdata, fftdatalen,
                                           rrfrac + fftdatalen / 2, 0.0, 0.0));
        vect_free(fftdata);
        fclose(fftfile);
    }

    /* Initialize the dftvector */

    init_dftvector(&dftvec, dataperread, cmd->numvect, idata.dt, rr, norm, T);

    /* Show our folding values */

    printf("\nFolding data from '%s':\n", infilenm);
    printf("   Folding Fourier Freq = %.5f\n", rr);
    printf("      Folding Freq (Hz) = %-.11f\n", rr / T);
    printf("     Folding Period (s) = %-.14f\n", T / rr);
    printf("  Points per sub-vector = %d\n", dftvec.n);
    printf("  Number of sub-vectors = %d\n", dftvec.numvect);
    printf(" Normalization constant = %g\n", norm * norm);

    /* Perform the actual vector addition */

    {
        int ii, jj;
        float *data;
        double real, imag, sumreal = 0.0, sumimag = 0.0;
        double theta, aa, bb, cc, ss, dtmp;
        double powargr, powargi, phsargr, phsargi, phstmp;

        data = gen_fvect(dftvec.n);
        theta = -TWOPI * rr / (double) N;
        dtmp = sin(0.5 * theta);
        aa = -2.0 * dtmp * dtmp;
        bb = sin(theta);
        cc = 1.0;
        ss = 0.0;
        for (ii = 0; ii < dftvec.numvect; ii++) {
            chkfread(data, sizeof(float), dftvec.n, infile);
            real = 0.0;
            imag = 0.0;
            for (jj = 0; jj < dftvec.n; jj++) {
                real += data[jj] * cc;
                imag += data[jj] * ss;
                cc = aa * (dtmp = cc) - bb * ss + cc;
                ss = aa * ss + bb * dtmp + ss;
            }
            dftvec.vector[ii].r = norm * real;
            dftvec.vector[ii].i = norm * imag;
            sumreal += dftvec.vector[ii].r;
            sumimag += dftvec.vector[ii].i;
        }
        vect_free(data);
        printf("\nDone:\n");
        printf("             Vector sum = %.3f + %.3fi\n", sumreal, sumimag);
        printf("      Total phase (deg) = %.2f\n", PHASE(sumreal, sumimag));
        printf("            Total power = %.2f\n", POWER(sumreal, sumimag));
        printf("\n");
    }
    fclose(infile);

    /* Write the output structure */

    sprintf(outfilenm, "%s_%.3f.dftvec", cmd->argv[0], rr);
    write_dftvector(&dftvec, outfilenm);

    /* Free our vector and return */

    free_dftvector(&dftvec);
    return (0);
}


void read_dftvector(dftvector * data, char *filename)
/* Read a dftvector data structure from a binary file */
{
    FILE *infile;
    int ii;
    double dtmp;

    infile = chkfopen(filename, "rb");
    chkfread(&dtmp, sizeof(double), 1, infile);
    data->n = (int) dtmp;
    chkfread(&dtmp, sizeof(double), 1, infile);
    data->numvect = (int) dtmp;
    chkfread(&data->dt, sizeof(double), 1, infile);
    chkfread(&data->r, sizeof(double), 1, infile);
    chkfread(&data->norm, sizeof(double), 1, infile);
    chkfread(&data->T, sizeof(double), 1, infile);
    data->vector = gen_cvect(data->n);
    for (ii = 0; ii < data->numvect; ii++) {
        chkfread(&dtmp, sizeof(double), 1, infile);
        data->vector[ii].r = (float) dtmp;
        chkfread(&dtmp, sizeof(double), 1, infile);
        data->vector[ii].i = (float) dtmp;
    }
    fclose(infile);
}

void write_dftvector(dftvector * data, char *filename)
/* Write a dftvector data structure to a binary file */
{
    FILE *outfile;
    int ii;
    double dtmp;

    outfile = chkfopen(filename, "wb");
    dtmp = (double) data->n;
    chkfwrite(&dtmp, sizeof(double), 1, outfile);
    dtmp = (double) data->numvect;
    chkfwrite(&dtmp, sizeof(double), 1, outfile);
    chkfwrite(&data->dt, sizeof(double), 1, outfile);
    chkfwrite(&data->r, sizeof(double), 1, outfile);
    chkfwrite(&data->norm, sizeof(double), 1, outfile);
    chkfwrite(&data->T, sizeof(double), 1, outfile);
    for (ii = 0; ii < data->numvect; ii++) {
        dtmp = (double) data->vector[ii].r;
        chkfwrite(&dtmp, sizeof(double), 1, outfile);
        dtmp = (double) data->vector[ii].i;
        chkfwrite(&dtmp, sizeof(double), 1, outfile);
    }
    fclose(outfile);
}

void init_dftvector(dftvector * data, int n, int numvect,
                    double dt, double r, double norm, double T)
/* Initialize a dftvector and allocate its vector part */
{
    data->n = n;
    data->numvect = numvect;
    data->dt = dt;
    data->r = r;
    data->norm = norm;
    data->T = T;
    data->vector = gen_cvect(n);
}

void free_dftvector(dftvector * data)
/* Free the dynamically allocated vector in data */
{
    vect_free(data->vector);
}
