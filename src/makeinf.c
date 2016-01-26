#include "makeinf.h"

extern char bands[NUMBANDS][40];
extern char scopes[NUMSCOPES][40];

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(void)
{
    static char filters[NUMFILTERS][7] = { "U", "B", "V", "R", "I", "J", "H",
        "K", "L", "M", "Other"
    };

    char temp1[100], temp2[100];
    int i, itemp, itemp2;
    infodata data;

    printf("\n\n");
    printf("   Info File  Generation Program\n");
    printf("        by Scott M. Ransom\n");
    printf("           23 Aug, 1997\n");

    printf("\n Enter data filename (no suffix) (100 char or less):\n");
    fgets(data.name, 100, stdin);
    data.name[strlen(data.name) - 1] = '\0';

    do {
        printf("\n Enter telescope used:\n");
        for (i = 1; i <= NUMSCOPES; i++) {
            printf("    %3d = %s\n", i, scopes[i - 1]);
        }
        fgets(temp1, 100, stdin);
        temp1[strlen(temp1) - 1] = '\0';
        itemp = atoi(temp1);
        strcpy(data.telescope, scopes[itemp - 1]);
    } while ((itemp > NUMSCOPES) || (itemp < 1));

    if (itemp > 1) {            /* If using makedata */

        printf("\n Enter instrument used:\n");
        fgets(data.instrument, 100, stdin);
        data.instrument[strlen(data.instrument) - 1] = '\0';

        printf("\n Enter object being observed (100 char or less):\n");
        fgets(data.object, 100, stdin);
        data.object[strlen(data.object) - 1] = '\0';

        do {
            printf("\n Enter J2000 Right Ascension (hh:mm:ss.ssss):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            sscanf(temp1, "%d:%d:%lf", &data.ra_h, &data.ra_m, &data.ra_s);
            if (data.ra_h >= 24.0 || data.ra_h < 0.0 ||
                data.ra_m >= 60.0 || data.ra_m < 0.0 ||
                data.ra_s >= 60.0 || data.ra_s < 0.0) {
                printf("\n   Not a valid Right ascension.\n");
                itemp2 = 1;
            } else {
                itemp2 = 0;
            }
        } while (itemp2);

        do {
            printf("\n Enter J2000 Declination (dd:mm:ss.ssss):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            sscanf(temp1, "%d:%d:%lf", &data.dec_d, &data.dec_m, &data.dec_s);
            if (data.dec_d > 90.0 || data.dec_d < -90.0 ||
                data.dec_m >= 60.0 || data.dec_m < 0.0 ||
                data.dec_s >= 60.0 || data.dec_s < 0.0) {
                printf("\n   Not a valid Declination.\n");
                itemp2 = 1;
            } else {
                itemp2 = 0;
            }
        } while (itemp2);

        printf("\n Who observed the data? (100 char or less):\n");
        fgets(data.observer, 100, stdin);
        data.observer[strlen(data.observer) - 1] = '\0';

        do {
            printf("\n Enter Epoch of observation (MJD) (include a decimal):\n");
            fgets(temp2, 100, stdin);
            temp2[strlen(temp2) - 1] = '\0';
            sscanf(temp2, "%d.%s", &data.mjd_i, temp1);
            sprintf(temp2, "0.%s", temp1);
            data.mjd_f = atof(temp2);
            if (data.mjd_i < 0) {
                printf("\n   Not a MJD.  Should be a positive number.\n");
                itemp2 = 1;
            } else {
                itemp2 = 0;
            }
        } while (itemp2);


        printf("\n Is the data barycentered?  1=yes, 0=no\n");
        fgets(temp1, 100, stdin);
        temp1[strlen(temp1) - 1] = '\0';
        data.bary = atoi(temp1);

    }
    printf("\n Enter number of bins in the time series :\n");
    fgets(temp1, 100, stdin);
    temp1[strlen(temp1) - 1] = '\0';
    data.N = atof(temp1);

    printf("\n Enter width of each time series bin (sec):\n");
    fgets(temp1, 100, stdin);
    temp1[strlen(temp1) - 1] = '\0';
    data.dt = atof(temp1);

    printf("\n Any breaks in the data?  1=yes, 0=no\n");
    fgets(temp1, 100, stdin);
    temp1[strlen(temp1) - 1] = '\0';
    data.numonoff = atoi(temp1);

    if (data.numonoff) {
        printf("\n   Type the on/off pairs starting from\n");
        printf("bin 0 and ending at bin N-1 (%11.0f).\n", data.N - 1);
        printf("Max of 20 on/off pairs. (The bins between\n");
        printf("paired numbers are where the data is \"on\"\n");
        printf("or un-padded.)\n");
        i = 0;
        do {
            scanf("%lf %lf", &data.onoff[i], &data.onoff[i + 1]);
            i += 2;
        } while (data.onoff[i - 1] != data.N - 1 && i < 40);
        data.numonoff = i / 2;
        fgets(temp1, 100, stdin);
        temp1[strlen(temp1) - 1] = '\0';
    } else {
        data.numonoff = 0;
        data.onoff[0] = 0;
        data.onoff[1] = data.N - 1;
    }
    if (itemp > 1) {            /* If using makedata */

        do {
            printf("\n Enter type of observation (EM band):\n");
            for (i = 1; i <= NUMBANDS; i++) {
                printf("    %3d = %s\n", i, bands[i - 1]);
            }
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            itemp = atoi(temp1);
            strcpy(data.band, bands[itemp - 1]);
        } while ((itemp > NUMBANDS) || (itemp < 1));

        if (itemp == 1) {

            printf("\n Enter approximate beam diameter (arcsec):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.fov = atof(temp1);

            printf("\n Enter dispersion measure (cm-3 pc):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.dm = atof(temp1);

            printf("\n Enter central frequency of the low channel (Mhz):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.freq = atof(temp1);

            printf("\n Enter total bandwidth (Mhz):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.freqband = atof(temp1);

            printf("\n Enter number of channels:\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.num_chan = atoi(temp1);

            printf("\n Enter channel bandwidth (Mhz):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.chan_wid = atof(temp1);

        } else if (itemp == 2 || itemp == 3 || itemp == 4) {

            printf("\n Enter approximate field-of-view diameter (arcsec):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.fov = atof(temp1);

            do {
                printf("\n Enter photometric filter used:\n");
                for (i = 1; i <= NUMFILTERS; i++) {
                    printf("    %3d = %s\n", i, filters[i - 1]);
                }
                fgets(temp1, 100, stdin);
                temp1[strlen(temp1) - 1] = '\0';
                itemp = atoi(temp1);
                strcpy(data.filt, filters[itemp - 1]);
            } while ((itemp > NUMFILTERS) || (itemp < 1));

            printf("\n Enter central wavelength (nm):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.wavelen = atof(temp1);

            printf("\n Enter bandpass (nm):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.waveband = atof(temp1);

        } else {

            printf("\n Enter approximate field-of-view diameter (arcsec):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.fov = atof(temp1);

            printf("\n Enter central energy (kev):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.energy = atof(temp1);

            printf("\n Enter energy bandpass (kev):\n");
            fgets(temp1, 100, stdin);
            temp1[strlen(temp1) - 1] = '\0';
            data.energyband = atof(temp1);

        }

    }
    printf("\n Who analyzed the data? (100 char or less):\n");
    fgets(data.analyzer, 100, stdin);
    data.analyzer[strlen(data.analyzer) - 1] = '\0';

    printf("\n Enter any additional notes: (200 char or less):\n");
    fgets(data.notes, 200, stdin);
    data.notes[strlen(data.notes) - 1] = '\0';

    writeinf(&data);
    printf("\nFinished.\n\n");

    exit(0);
}
