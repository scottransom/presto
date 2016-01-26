#include "presto.h"
#include "plot2d.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *infile;
    int numchan, numtodisplay = 0, numprofs = 0;
    long i, j, proflen, offset;
    double *profs, *sumprof, nc, pl, tt, bt, p, f, df, rotate = 0.0;
    double dm = 0.0, bindelay = 0.0;
    char device[200], output[200];

    if (argc <= 2) {
        printf
            ("usage:  showmulti_dm filename dm [rotate] [numtodisplay] [device]\n");
        printf("   'filename'     (required) is the multi-profile save file.\n");
        printf("   'dm'           (required) DM to base profile delays on.\n");
        printf("   'rotate'       (optional) is the number of bins to rotate\n");
        printf("                             each profile to the left.\n");
        printf("                             Can be fractional.  Default is 0.\n");
        printf("   'numtodisplay' (optional) is the number of profiles to\n");
        printf("                             display at once.  Defaults to\n");
        printf("                             the number of channels.\n");
        printf("   'device'       (optional) is the pgplot device to use ('x' or\n");
        printf("                             'ps').  Defaults to 'x'\n");
        exit(1);
    }

    infile = chkfopen(argv[1], "rb");
    sprintf(output, "%s.ps", argv[1]);
    chkfread(&nc, sizeof(double), 1, infile);
    chkfread(&pl, sizeof(double), 1, infile);
    chkfread(&p, sizeof(double), 1, infile);
    chkfread(&tt, sizeof(double), 1, infile);
    chkfread(&bt, sizeof(double), 1, infile);
    chkfread(&f, sizeof(double), 1, infile);
    chkfread(&df, sizeof(double), 1, infile);
    numchan = nc;
    proflen = pl;

    if (argc == 3) {
        dm = strtod(argv[2], NULL);
        rotate = 0.0;
        strcpy(device, "x");
    } else if (argc == 4) {
        dm = strtod(argv[2], NULL);
        rotate = strtod(argv[3], NULL);
        numtodisplay = numchan;
        strcpy(device, "x");
    } else if (argc == 5) {
        dm = strtod(argv[2], NULL);
        rotate = strtod(argv[3], NULL);
        numtodisplay = (int) strtol(argv[4], NULL, 10);
        strcpy(device, "x");
    } else if (argc == 6) {
        dm = strtod(argv[2], NULL);
        rotate = strtod(argv[3], NULL);
        numtodisplay = (int) strtol(argv[4], NULL, 10);
        strcpy(device, argv[5]);
    }

    printf("\n      Multi-Profile Display Program\n");
    printf("         With DM Delay Correction\n");
    printf("              Scott M. Ransom\n");
    printf("               20 July 1998\n");
    printf("\nProfile properties:\n");
    printf("Initial folding period   (s)  =  %-15.13f\n", p);
    printf("Topocentric time     (start)  =  %-15.10f\n", tt);
    printf("Barycentric time     (start)  =  %-15.10f\n", bt);
    printf("Profile length        (bins)  =  %-ld\n", proflen);
    printf("Number of channels            =  %-d\n", numchan);
    printf("Channel 1 frequency    (MHz)  =  %-10.5f\n", f);
    printf("Channel freq width     (MHz)  =  %-10.5f\n", df);
    printf("Dispersion Measure (cm-3 pc)  =  %-10.5f\n\n", dm);

    /* Read the profiles. */

    profs = gen_dvect(proflen * numchan);
    chkfread(profs, sizeof(double), (unsigned long) (numchan * proflen), infile);
    fclose(infile);

    /* Create a Summed-Profile vector */

    sumprof = gen_dvect(proflen);
    for (i = 0; i < proflen; i++) {
        sumprof[i] = 0.0;
    }

    /* Rotate the vectors and sum the profiles */

    for (i = 0; i < numchan; i++) {
        bindelay = delay_from_dm(dm, f + i * df) * (double) proflen / p;
        drotate(&profs[i * proflen], proflen, bindelay);
        if (rotate)
            drotate(&profs[i * proflen], proflen, rotate);
        for (j = 0; j < proflen; j++) {
            sumprof[j] += profs[i * proflen + j];
        }
    }

    /* Plot the profiles */

    if (0 == strcmp("x", device))
        cpgstart_x("portrait");
    else
        cpgstart_ps(output, "portrait");
    for (i = 0; i <= numchan / numtodisplay; i++) {
        offset = i * numtodisplay;
        numprofs =
            (numchan - offset) > numtodisplay ? numtodisplay : numchan - offset;
        if (numprofs > 0) {
            cpgpage();
            multi_prof_plot(proflen, numprofs, profs + offset * proflen,
                            sumprof, "Pulse Phase (Periods)",
                            (double) (1.0 + offset), 1.0, "Channel Number",
                            f + offset * df, df, "Frequency (MHz)");
        }
    }
    cpgend();

    /* Cleanup */

    vect_free(profs);
    vect_free(sumprof);

    return 0;
}
