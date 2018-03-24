#include "presto.h"

#define MAXTIMES 5000

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
/* Convert topocentric arrival times (from stdin) to */
/* barycentric times using TEMPO.                    */
{
    long i, N = 0;
    char obs[3], ra[30], dec[30], ephem[10], tmp[40];
    double topof = 0.0, dm = 0.0;
    double topotimes[MAXTIMES], barytimes[MAXTIMES];
    double voverc[MAXTIMES], dmdelays[MAXTIMES];

    if (argc < 4) {
        printf("\n   Usage:  bary observatory ra dec [freq (MHz)] [dm] [ephem]\n\n");
        printf
            ("   Convert topocentric times (from stdin) to barycentric times.\n\n");
        printf("   Arguments:\n");
        printf
            ("     observatory is a 2 letter ITOA observatory code in obsys.dat\n");
        printf("     ra is a string with J2000 RA in the format 'hh:mm:ss.ssss'\n");
        printf
            ("     dec is a string with J2000 DEC in the format 'dd:mm:ss.ssss'\n");
        printf("     freq is the (optional) observing frequency in MHz\n");
        printf("        freq defaults to infinite frequency (0.0).\n");
        printf
            ("     dm is the (optional) Dispersion Measure of the obs (cm^-3 pc)\n");
        printf("        dm defaults to 0.0 cm^-3 pc.\n");
        printf("     ephem is the (optional) ephemeris to use.  Must be supported by TEMPO.\n");
        printf("        Examples include 'DE200', 'DE421', or 'DE436'.  Defaults is 'DE421'.\n\n");
        printf("   Notes:  The topocentric times must be in UTC MJD format.\n");
        printf("     There is a maximum limit of 5000 input times.\n\n");
        exit(0);
    }

    strcpy(ephem, "DE421"); // Default choice
    if (argc >= 4) {
        strcpy(obs, argv[1]);
        strcpy(ra, argv[2]);
        strcpy(dec, argv[3]);
    }
    if (argc >= 5)
        topof = strtod(argv[4], NULL);
    if (argc >= 6)
        dm = strtod(argv[5], NULL);
    if (argc >= 7) {
        if (strlen(argv[6])==5 && argv[6][0]=='D' && argv[6][1]=='E')
            strcpy(ephem, argv[6]);
    }

    /* Read in the topocentric times */

    while (fgets(tmp, 40, stdin) != NULL) {
        if (strcmp(".\n", tmp) == 0 || strcmp("\n", tmp) == 0)
            break;
        topotimes[N] = strtod(tmp, NULL);
        N++;
    }

    fprintf(stderr, "\n\n      Barycentric Time Converter Hack\n");
    fprintf(stderr, "            by Scott M. Ransom\n");
    fprintf(stderr, "               20 July 1998\n\n");

    fprintf(stderr, "  Calling TEMPO with %ld topocentic time(s).  Using:\n", N);
    fprintf(stderr, "                  RA = '%s'\n", ra);
    fprintf(stderr, "                 DEC = '%s'\n", dec);
    fprintf(stderr, "      Obs freq (MHz) = %.6g\n", topof);
    fprintf(stderr, "        DM (pc/cm^3) = %.4g\n", dm);
    fprintf(stderr, "           Ephemeris = '%s'\n\n", ephem);
    /* Call TEMPO */

    barycenter(topotimes, barytimes, voverc, N, ra, dec, obs, ephem);

    /* Determine DM Delays */

    if (dm == 0.0 || topof == 0.0 || topof >= 10000.0) {
        for (i = 0; i < N; i++)
            dmdelays[i] = 0.0;
    } else {
        for (i = 0; i < N; i++)
            dmdelays[i] = delay_from_dm(dm, doppler(topof, voverc[i]));
    }

    fprintf(stderr,
            "   Obs Freq          Inf Freq          Radial Vel     Emitted    DM Delay \n");
    fprintf(stderr,
            " Topotimes (MJD)   Barytimes (MJD)       (v/c)       Freq (MHz)   (sec)   \n");
    fprintf(stderr,
            "--------------------------------------------------------------------------\n");

    for (i = 0; i < N; i++)
        fprintf(stdout, "%16.10f  %16.10f  %14.7g  %11.6f  %9.7f\n",
                topotimes[i], barytimes[i] - dmdelays[i]/86400.0, voverc[i],
                doppler(topof, voverc[i]), dmdelays[i]);

    fprintf(stderr, "\n\n");
    exit(0);
}
