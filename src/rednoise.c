#include "presto.h"
#include "rednoise_cmd.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *infile, *outfile = NULL;
    char *rootfilenm, *newroot, *outname;
    long long numbins;
    long numwrote;
    double T;
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
    printf("     Rednoise Removal Routine v3.0\n");
    printf("           August, 2026\n\n");

    /* The input file must be a '.fft' file */

    {
        int hassuffix = 0;
        char *suffix;

        hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);
        if (hassuffix && strcmp(suffix, "fft") == 0) {
            free(suffix);
        } else {
            printf("\nInput file ('%s') must be a fourier transform ('.fft')!\n\n",
                   cmd->argv[0]);
            if (hassuffix)
                free(suffix);
            exit(1);
        }
    }

    /* The output files are written next to the input file, with '_red' */
    /* appended to the root of its name.                                */

    newroot = (char *) calloc(strlen(rootfilenm) + strlen("_red") + 1, 1);
    sprintf(newroot, "%s_red", rootfilenm);
    outname = (char *) calloc(strlen(newroot) + strlen(".fft") + 1, 1);
    sprintf(outname, "%s.fft", newroot);

    /* Read the info file so that we know how long the observation was */

    readinf(&idata, rootfilenm);
    T = idata.N * idata.dt;
    if (strlen(newroot) >= sizeof(idata.name)) {
        printf("\nOutput file name '%s' is too long for a '.inf' file!\n\n", newroot);
        exit(1);
    }
    strcpy(idata.name, newroot);

    /* De-redden the amplitudes */

    if (cmd->inplaceP) {
        infile = chkfopen(cmd->argv[0], "rb+");
        numbins = chkfilelen(infile, sizeof(fcomplex));
        printf("Dereddening '%s' in place.\n", cmd->argv[0]);
        numwrote = deredden_file(infile, infile, numbins, cmd->startwidth,
                                 cmd->endwidth, cmd->endfreq, T);
        fclose(infile);
        if (numwrote > 0 && rename(cmd->argv[0], outname)) {
            perror("\nError renaming the dereddened file");
            printf("\n");
            exit(1);
        }
    } else {
        infile = chkfopen(cmd->argv[0], "rb");
        numbins = chkfilelen(infile, sizeof(fcomplex));
        outfile = chkfopen(outname, "wb");
        numwrote = deredden_file(infile, outfile, numbins, cmd->startwidth,
                                 cmd->endwidth, cmd->endfreq, T);
        fclose(infile);
        fclose(outfile);
    }
    if (numwrote < 0) {
        printf("\nRednoise removal failed.\n\n");
        exit(1);
    }
    writeinf(&idata);

    printf("\nDone.  Rednoise removed from %ld of %lld points.\n\n",
           numwrote, numbins);
    printf("   Dereddened fft file: %s\n", outname);
    printf("Corresponding inf file: %s.inf\n\n", newroot);

    free(rootfilenm);
    free(newroot);
    free(outname);
    exit(0);
}
