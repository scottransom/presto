#include "makedata.h"

void read_mak_input(makedata * mdata)
/* Get information for makedata from stdin.  */
{
    int i;
    double tmponoff[40];

    printf("\nEnter the filename to save the data (no ");
    printf("extensions):\n");
    fgets(mdata->basefilenm, 200, stdin);
    mdata->basefilenm[strlen(mdata->basefilenm) - 1] = '\0';

    printf("\nEnter a 200 character (or less) description ");
    printf("of the data set.\n");
    fgets(mdata->description, 200, stdin);
    mdata->description[strlen(mdata->description) - 1] = '\0';

    printf("\nEnter the number of data points and the bin ");
    printf("time dt (sec):\n");
    scanf("%ld %lf", &mdata->N, &mdata->dt);
    while (mdata->N <= 0 || mdata->dt <= 0) {
        printf("\nBoth N and dt must be >= 0 in read_mak_input()\n");
        printf("\nEnter the number of data points and the bin ");
        printf("time dt (sec):\n");
        scanf("%ld %lf", &mdata->N, &mdata->dt);
    }
    mdata->T = mdata->N * mdata->dt;
    mdata->next2_to_n = 1;
    while (mdata->next2_to_n < mdata->N) {
        mdata->next2_to_n <<= 1;
    }

    printf("\nWhat type of Pulse shape?\n");
    printf("  1=Sine, 2=Crab-like, 3=Spike, 4=Gauss:\n");
    scanf("%d", &mdata->pnum);
    while (mdata->pnum < 1 || mdata->pnum > 4) {
        printf("\nNot a proper pulse shape choice in read_mak_input()\n");
        printf("\nWhat type of Pulse shape?\n");
        printf("  1=Sine, 2=Crab-like, 3=Spike, 4=Gauss:\n");
        scanf("%d", &mdata->pnum);
    }
    switch (mdata->pnum) {
    case 1:
        strcpy(mdata->ptype, "Sine");
        break;
    case 2:
        strcpy(mdata->ptype, "Crab-like");
        break;
    case 3:
        strcpy(mdata->ptype, "Spike");
        printf("\nEnter the phase width (0.0-0.5) FWHM of the pulse.\n");
        scanf("%lf", &mdata->fwhm);
        while (mdata->fwhm < 0 || mdata->fwhm > 0.5) {
            printf("\nNot a proper FWHM in read_mak_input()\n");
            printf("\nEnter the phase width (0.0-0.5) FWHM of the pulse.\n");
            scanf("%lf", &mdata->fwhm);
        }
        break;
    case 4:
        strcpy(mdata->ptype, "Gauss");
        printf("\nEnter the phase width (0.0-0.5) FWHM of the pulse.\n");
        scanf("%lf", &mdata->fwhm);
        while (mdata->fwhm < 0 || mdata->fwhm > 0.5) {
            printf("\nNot a proper FWHM in read_mak_input()\n");
            printf("\nEnter the phase width (0.0-0.5) FWHM of the pulse.\n");
            scanf("%lf", &mdata->fwhm);
        }
        break;
    }

    printf("\nDo you want the binned data rounded to ");
    printf("whole numbers?\n");
    printf("  1=yes (Whole number (Photon) data), ");
    printf("0=no (Fractional data):\n");
    scanf("%d", &mdata->roundnum);
    if (mdata->roundnum)
        strcpy(mdata->round, "Whole Numbers");
    else
        strcpy(mdata->round, "Fractional");

    printf("\nEnter pulsation freq (hz), freq deriv (s-2),\n");
    printf("  and freq 2nd deriv (s-3):\n");
    scanf("%lf %lf %lf", &mdata->f, &mdata->fd, &mdata->fdd);
    while (mdata->f <= 0) {
        printf("\nNot a proper pulsation frequency in read_mak_input()\n");
        printf("\nEnter pulsation freq (hz), freq deriv (s-2),\n");
        printf("  and freq 2nd deriv (s-3):\n");
        scanf("%lf %lf %lf", &mdata->f, &mdata->fd, &mdata->fdd);
    }
    mdata->p = 1.0 / mdata->f;
    mdata->pd = -mdata->fd / (mdata->f * mdata->f);
    if (mdata->fdd == 0.0)
        mdata->pdd = 0.0;
    else
        mdata->pdd = (2 * mdata->fd * mdata->fd / mdata->f -
                      mdata->fdd) / mdata->f * mdata->f;
    mdata->r = mdata->f * mdata->T;
    mdata->z = mdata->fd * mdata->T * mdata->T;
    mdata->w = mdata->fdd * mdata->T * mdata->T * mdata->T;

    printf("\nEnter pulsation amplitude, phase (deg), and\n");
    printf("  average background level per bin (DC offset):\n");
    scanf("%lf %lf %lf", &mdata->amp, &mdata->phs, &mdata->dc);
    while (mdata->amp < 0 || mdata->phs < 0 || mdata->phs >= 360) {
        if (mdata->amp < 0) {
            printf("\nAmplitude must be >= 0.0 in read_mak_input()\n");
        } else {
            printf("\nPhase must be (0 <= phs < 360) in read_mak_input()\n");
        }
        printf("\nEnter pulsation amplitude, phase (deg), and\n");
        printf("  average background level per bin (DC offset):\n");
        scanf("%lf %lf %lf", &mdata->amp, &mdata->phs, &mdata->dc);
    }

    printf("\nPulsar in a binary orbit?\n");
    printf("  1=yes, 0=no:\n");
    scanf("%d", &mdata->binary);
    if (mdata->binary) {
        printf("Input binary period (s), asini/c (s), ");
        printf("e, omega (deg),\n");
        printf("  and time since last periastron (s):\n");
        scanf("%lf %lf %lf %lf %lf", &mdata->orb.p, &mdata->orb.x,
              &mdata->orb.e, &mdata->orb.w, &mdata->orb.t);
        while (mdata->orb.p < 0.0) {
            printf("\nPeriod must be > 0.  Input binary period (s):\n");
            scanf("%lf", &mdata->orb.p);
        }
        while (mdata->orb.x < 0.0) {
            printf("\na_sin(i)/c must be > 0.  Input a_sin(i)/c (lt-s):\n");
            scanf("%lf", &mdata->orb.x);
        }
        while (mdata->orb.e < 0.0 || mdata->orb.e >= 1.0) {
            printf("\nEccentricity must be (0 <= e < 1.0).  Input e:\n");
            scanf("%lf", &mdata->orb.e);
        }
        while (mdata->orb.w < 0.0 || mdata->orb.w >= 360.0) {
            printf("\nOmega must be (0 <= omega < 360).  Input omega (deg):\n");
            scanf("%lf", &mdata->orb.w);
        }
    }

    printf("\nAmplitude modulation of the signal?\n");
    printf("  1=yes, 0=no:\n");
    scanf("%d", &mdata->ampmod);
    if (mdata->ampmod) {
        printf("Input modulation amplitude, frequency (hz)");
        printf("and phase (deg):\n");
        scanf("%lf %lf %lf", &mdata->ampmoda, &mdata->ampmodf, &mdata->ampmodp);
        while (mdata->ampmodf < 0 || mdata->ampmodp < 0 || mdata->ampmodp >= 360) {
            if (mdata->ampmodf < 0) {
                printf("\nFrequency must be >= 0.0 in read_mak_input()\n");
            } else {
                printf("\nPhase must be (0 <= phs < 360) in read_mak_input()\n");
            }
            printf("Input modulation amplitude, frequency (hz)");
            printf(" and phase (deg):\n");
            scanf("%lf %lf %lf", &mdata->ampmoda, &mdata->ampmodf, &mdata->ampmodp);
        }
    }

    printf("\nWhat type of Noise?\n");
    printf("  1=Standard noise level (poissonian ");
    printf("with siqma=sqrt of bin level),\n");
    printf("  2=Other noise level:\n");
    scanf("%d", &mdata->noise);
    if (mdata->noise == 2) {
        strcpy(mdata->noisetype, "Other");
        printf("Enter standard deviation of normally distributed noise:\n");
        scanf("%lf", &mdata->noisesig);
        while (mdata->noisesig < 0) {
            printf("Noise stdev must be >= 0.0.  Enter noise stdev:\n");
            scanf("%lf", &mdata->noisesig);
        }
    } else {
        strcpy(mdata->noisetype, "Standard");
        mdata->noisesig = sqrt(mdata->dc);
    }

    printf("\nDo you want the signal to turn off and ");
    printf("on during the run?\n");
    printf("  1=yes, 0=no:\n");
    scanf("%d", &mdata->numonoff);
    if (mdata->numonoff != 0) {
        printf("Type the on/off pairs starting from ");
        printf("0.0 and ending at 1.0.\n");
        printf("  Max of 40 on/off pairs (Time between ");
        printf("paired numbers is on).\n");
        i = 0;
        do {
            scanf("%lf %lf", &tmponoff[i], &tmponoff[i + 1]);
            i += 2;
        } while (tmponoff[i - 1] < 1.0 && i < 40);
        if (tmponoff[i - 1] > 1.0)
            tmponoff[i - 1] = 1.0;
        mdata->numonoff = i / 2;
        mdata->onoff = (double *) malloc(mdata->numonoff * 2 * sizeof(double));
        for (i = 0; i < mdata->numonoff; i++) {
            mdata->onoff[2 * i] = tmponoff[2 * i];
            mdata->onoff[2 * i + 1] = tmponoff[2 * i + 1];
        }
    } else {
        mdata->numonoff = 1;
        mdata->onoff = (double *) malloc(2 * sizeof(double));
        mdata->onoff[0] = 0.0;
        mdata->onoff[1] = 1.0;
    }
    printf("\n");
}


void read_mak_file(char basefilenm[], makedata * mdata)
/* Read the data for makedata from the makefile. */
{
    FILE *makefile;
    char makefilenm[200], tmp[50];
    double tmponoff[40];
    int i;

    sprintf(makefilenm, "%s.mak", basefilenm);
    printf("Reading make information from \"%s\".\n", makefilenm);
    makefile = chkfopen(makefilenm, "r");

    strcpy(mdata->basefilenm, basefilenm);
    fscanf(makefile, "%[^\n]", mdata->description);

    fscanf(makefile, "%*[^=]= %ld", &mdata->N);
    if (mdata->N <= 0) {
        printf("\nmdata->N must be > 0 in read_mak_file()\n");
        exit(1);
    }

    mdata->next2_to_n = 1;
    while (mdata->next2_to_n < mdata->N) {
        mdata->next2_to_n <<= 1;
    }

    fscanf(makefile, "%*[^=]= %lf", &mdata->dt);
    if (mdata->dt <= 0) {
        printf("\nmdata->dt must be > 0 in read_mak_file()\n");
        exit(1);
    }

    mdata->T = mdata->N * mdata->dt;

    i = fscanf(makefile, "%*[^=]= %s%[^\n]", mdata->ptype, tmp);
    if (strcmp(mdata->ptype, "Sine") == 0)
        mdata->pnum = 1;
    else if (strcmp(mdata->ptype, "Crab-like") == 0)
        mdata->pnum = 2;
    else {
        if (strcmp(mdata->ptype, "Spike") == 0) {
            mdata->pnum = 3;
        } else {
            mdata->pnum = 4;
        }
        /* Default value for fwhm */
        if (i == 2)
            sscanf(tmp, "%*[^=]= %lf", &mdata->fwhm);
        else
            mdata->fwhm = 0.1;
        if (mdata->fwhm <= 0.0 || mdata->fwhm > 0.5) {
            printf("\nmdata->fwhm must be between 0 and 0.5 in ");
            printf("read_mak_file()\n");
            exit(1);
        }
    }

    fscanf(makefile, "%*[^=]= %s", mdata->round);
    if (strcmp(mdata->round, "Fractional") == 0)
        mdata->roundnum = 0;
    else
        mdata->roundnum = 1;

    fscanf(makefile, "%*[^=]= %lf", &mdata->f);
    if (mdata->f <= 0.0) {
        printf("\nmdata->f must be > 0.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->fd);
    fscanf(makefile, "%*[^=]= %lf", &mdata->fdd);
    mdata->p = 1.0 / mdata->f;
    mdata->pd = -mdata->fd / (mdata->f * mdata->f);
    if (mdata->fdd == 0.0)
        mdata->pdd = 0.0;
    else
        mdata->pdd = (2 * mdata->fd * mdata->fd / mdata->f -
                      mdata->fdd) / mdata->f * mdata->f;
    mdata->r = mdata->f * mdata->T;
    mdata->z = mdata->fd * mdata->T * mdata->T;
    mdata->w = mdata->fdd * mdata->T * mdata->T * mdata->T;

    fscanf(makefile, "%*[^=]= %lf", &mdata->amp);
    if (mdata->amp < 0.0) {
        printf("\nmdata->amp must be >= 0.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->phs);
    if (mdata->phs < 0.0 || mdata->phs >= 360) {
        printf("\nmdata->phs must be 0 <= phs < 360 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->dc);

    fscanf(makefile, "%*[^=]= %lf", &mdata->orb.p);
    if (mdata->orb.p < 0.0) {
        printf("\nmdata->orb.p must be >= 0.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->orb.x);
    if (mdata->orb.x < 0.0) {
        printf("\nmdata->orb.x must be >= 0.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->orb.e);
    if (mdata->orb.e < 0.0 || mdata->orb.e >= 1.0) {
        printf("\nmdata->e must be 0 <= e < 1.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->orb.w);
    if (mdata->orb.w < 0.0 || mdata->orb.w >= 360) {
        printf("\nmdata->orb.w must be 0 <= w < 360 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->orb.t);
    if (mdata->orb.p == 0.0 || mdata->orb.x == 0.0)
        mdata->binary = 0;
    else
        mdata->binary = 1;

    fscanf(makefile, "%*[^=]= %lf", &mdata->ampmoda);
    if (mdata->ampmoda < 0.0) {
        printf("\nmdata->ampmoda must be >= 0.0 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->ampmodp);
    if (mdata->ampmodp < 0.0 || mdata->ampmodp >= 360) {
        printf("\nmdata->ampmodp must be 0 <= phs < 360 in read_mak_file()\n");
        exit(1);
    }
    fscanf(makefile, "%*[^=]= %lf", &mdata->ampmodf);
    if (mdata->ampmodf < 0.0) {
        printf("\nmdata->ampmodf must be >= 0.0 in read_mak_file()\n");
        exit(1);
    }
    if (mdata->ampmoda == 0.0 || mdata->ampmodf == 0.0)
        mdata->ampmod = 0;
    else
        mdata->ampmod = 1;

    fscanf(makefile, "%*[^=]= %s", mdata->noisetype);
    if (strcmp(mdata->noisetype, "Standard") == 0)
        mdata->noise = 1;
    else
        mdata->noise = 2;
    fscanf(makefile, "%*[^=]= %lf", &mdata->noisesig);

    i = 0;
    do {
        fscanf(makefile, "%*[^=]= %lf %lf", &tmponoff[i], &tmponoff[i + 1]);
        i += 2;
    } while (tmponoff[i - 1] < 1.0 && i < 40);
    if (tmponoff[i - 1] > 1.0)
        tmponoff[i - 1] = 1.0;
    mdata->numonoff = i / 2;
    mdata->onoff = (double *) malloc(mdata->numonoff * 2 * sizeof(double));
    for (i = 0; i < mdata->numonoff; i++) {
        mdata->onoff[2 * i] = tmponoff[2 * i];
        mdata->onoff[2 * i + 1] = tmponoff[2 * i + 1];
    }
    fclose(makefile);
}


void write_mak_file(makedata * mdata)
/* Write the data for makedata to the makefile.  */
{
    FILE *makefile;
    char makefilenm[200];
    int i;

    sprintf(makefilenm, "%s.mak", mdata->basefilenm);
    makefile = chkfopen(makefilenm, "w");

    fprintf(makefile, "%s\n", mdata->description);
    fprintf(makefile, "Num data pts      = %ld\n", mdata->N);
    fprintf(makefile, "dt per bin (s)    = %25.15g\n", mdata->dt);
    fprintf(makefile, "Pulse shape       = %s", mdata->ptype);
    if (mdata->pnum > 2) {
        fprintf(makefile, " (FWHM = %g)", mdata->fwhm);
    }
    fprintf(makefile, "\nRounding format   = %s\n", mdata->round);
    fprintf(makefile, "Pulse freq (hz)   = %25.15g\n", mdata->f);
    fprintf(makefile, "fdot (s-2)        = %25.15g\n", mdata->fd);
    fprintf(makefile, "fdotdot (s-3)     = %25.15g\n", mdata->fdd);
    fprintf(makefile, "Pulse amp         = %25.15g\n", mdata->amp);
    fprintf(makefile, "Pulse phs (deg)   = %25.15g\n", mdata->phs);
    fprintf(makefile, "DC backgrnd level = %25.15g\n", mdata->dc);
    fprintf(makefile, "Binary period (s) = %25.15g\n",
            ((mdata->binary) ? mdata->orb.p : 0.0));
    fprintf(makefile, "Bin asini/c (s)   = %25.15g\n",
            ((mdata->binary) ? mdata->orb.x : 0.0));
    fprintf(makefile, "Bin eccentricity  = %25.15g\n",
            ((mdata->binary) ? mdata->orb.e : 0.0));
    fprintf(makefile, "Ang of Peri (deg) = %25.15g\n",
            ((mdata->binary) ? mdata->orb.w : 0.0));
    fprintf(makefile, "Tm since peri (s) = %25.15g\n",
            ((mdata->binary) ? mdata->orb.t : 0.0));
    fprintf(makefile, "Amp Mod amplitude = %25.15g\n",
            ((mdata->ampmod) ? mdata->ampmoda : 0.0));
    fprintf(makefile, "Amp Mod phs (deg) = %25.15g\n",
            ((mdata->ampmod) ? mdata->ampmodp : 0.0));
    fprintf(makefile, "Amp Mod freq (hz) = %25.15g\n",
            ((mdata->ampmod) ? mdata->ampmodf : 0.0));
    fprintf(makefile, "Noise type        = %s\n", mdata->noisetype);
    fprintf(makefile, "Noise sigma       = %25.15g\n", mdata->noisesig);
    i = 0;
    do {
        fprintf(makefile, "On/Off Pair %2d    = %g %g\n",
                i + 1, mdata->onoff[2 * i], mdata->onoff[2 * i + 1]);
        i++;
    }
    while (mdata->onoff[2 * i - 1] != 1.0);

    fclose(makefile);
}
