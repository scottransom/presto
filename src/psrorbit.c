#include "presto.h"
#include "plot2d.h"
#include "orbint.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    FILE *makefile;
    double p_psr = 0.0, pdot_psr, N, dt, epoch = 0.0, avgz;
    double *Ep, *tmpE, startE, endtime = 0.0, orbdt = 2.0;
    int i, ps = 0, pnum;
    long numpoints;
    float *trl, *Erl, *Erl2;
    double onoffpairs[30];
    char type = 'a', zz;
    char filenm[100], infofilenm[100], makefilenm[100], psfilenm[100];
    makedata mdata;
    orbitparams orb;
    infodata idata;
    psrparams psr;

    if ((argc < 2) || (argc > 5)) {
        printf("\nUsage:  'psrorbit [filename] [PSR name] [type] [ps]'\n\n");
        printf("  Optional arguments:\n");
        printf("         'filename' = the data's infofile name (without '.inf').\n");
        printf("         'PSR Name' = the pulsar's J2000 or B1950 name.\n");
        printf("                      Do not include 'PSR', 'J', or 'B'.\n");
        printf("                      If PSR Name is not present, porbit\n");
        printf("                      will look for the data in 'filename.mak'\n");
        printf
            ("             'type' = a character describing which graphic output:\n");
        printf("                      'v' = velocity in km/sec\n");
        printf("                      'd' = pulse time delay in sec\n");
        printf("                      'p' = pulsar period in ms\n");
        printf("                      'z' = instantaneous fourier f-dot in bins\n");
        printf("                      'a' = (the default) all of the above\n");
        printf("               'ps' = if 'ps' is present, the output graphs will\n");
        printf
            ("                         go to postscript file 'filename_orb.ps'\n\n");
        printf("If used with at least 1 argument, psorbit will either find\n");
        printf("the data it needs, or it will ask you for all required\n");
        printf("data to determine the pulsar orbit.\n");
        printf("                                        26 Sep 1997\n\n");
        exit(0);
    }
    printf("\n\n");
    printf("  Pulsar Binary Orbit Integration Routine\n");
    printf("           by Scott M. Ransom\n");
    printf("              26 Sep, 1997\n\n");

    sprintf(psfilenm, "psrorbit.ps");
    if (argc > 1 && strlen(argv[1]) > 2) {
        sprintf(filenm, "%s", argv[1]);
        sprintf(infofilenm, "%s.inf", argv[1]);
        sprintf(makefilenm, "%s.mak", argv[1]);
        sprintf(psfilenm, "%s_orb.ps", argv[1]);
        readinf(&idata, argv[1]);
        epoch = (double) idata.mjd_i + idata.mjd_f;
        dt = idata.dt;
        N = idata.N;
        endtime = dt * N;
    }
    i = 2;
    if (argc >= 3) {
        while (i < argc) {

            /* Read from the pulsar database */

            if (strlen(argv[i]) > 2) {

                pnum = get_psr_at_epoch(argv[i], epoch, &psr);
                if (!pnum) {
                    printf("\nThe pulsar is not in the database.  Exiting.\n\n");
                    exit(1);
                }
                if (psr.orb.p != 0.0) {
                    orb = psr.orb;
                } else {
                    printf("\nPulsar is not in a binary.  Exiting.\n\n");
                    exit(1);
                }
                p_psr = psr.p;
                pdot_psr = psr.pd;
            }

            /* Determine which plots to make.  Default is all. */

            if (strlen(argv[i]) == 1) {
                type = argv[i][0];
                if ((type != 'v') && (type != 'd') && (type != 'p')
                    && (type != 'z') && (type != 'a')) {
                    printf("The 'type' of plot selection is invalid.  Exiting.\n\n");
                    exit(1);
                }
            }
            /* Send output to the screen or to a file. */

            if ((strlen(argv[i]) == 2) && (!strcmp(argv[i], "ps"))) {
                ps = 1;
            }
            i++;
        }
    }
    /* Read from the makefile or ask for parameters */

    if (!p_psr) {
        if (!(makefile = fopen(makefilenm, "r"))) {
            printf("\nYou need to enter the required data:\n");
            printf
                ("\nEnter the number of data points and the time between bins (s):\n");
            scanf("%lf", &N);
            scanf("%lf", &dt);
            printf("\nEnter the pulsar period in seconds:\n");
            scanf("%lf", &p_psr);
            printf("\nEnter the binary info in the following order (no commas):\n");
            printf("  period (s), a_sin(i)/c (s), eccen, ang of peri");
            printf(" (deg), time of peri (s)\n");
            scanf("%lf", &orb.p);
            scanf("%lf", &orb.x);
            scanf("%lf", &orb.e);
            scanf("%lf", &orb.w);
            scanf("%lf", &orb.t);
            zz = getchar();
            orb.t = fmod(orb.t, orb.p);
            printf("\n");
            onoffpairs[0] = 0.0;
            onoffpairs[1] = 1.0;
            endtime = N * dt;
        } else {
            fclose(makefile);
            read_mak_file(filenm, &mdata);
            orb = mdata.orb;
            p_psr = mdata.p;
            N = mdata.N;
            dt = mdata.dt;
            endtime = mdata.T;
        }
    }
    if (ps) {
        cpgstart_ps(psfilenm, "landscape");
    } else {
        cpgstart_x("landscape");
    }

    /* Prep the orbit integrator */

    startE = keplers_eqn(orb.t, orb.p, orb.e, 1.0E-15);
    if (endtime > 2048)
        orbdt = 2.0;
    else
        orbdt = endtime / 1024.0;
    numpoints = (long) floor(endtime / orbdt + 0.5) + 1;
    Ep = dorbint(startE, numpoints, orbdt, &orb);
    tmpE = gen_dvect(numpoints);
    trl = gen_fvect(numpoints);
    Erl = gen_fvect(numpoints);
    Erl2 = gen_fvect(numpoints);

    for (i = 0; i < numpoints; i++) {
        trl[i] = orbdt * i;
        tmpE[i] = Ep[i];
    }

    /* Calculate the period change and velocity */

    if (type == 'v' || type == 'p' || type == 'a') {
        E_to_p(Ep, numpoints, p_psr, &orb);
        for (i = 0; i < numpoints; i++)
            Erl[i] = Ep[i] * 1000.0;
        for (i = 0; i < numpoints; i++)
            Ep[i] = tmpE[i];
        E_to_v(Ep, numpoints, &orb);
        for (i = 0; i < numpoints; i++)
            Erl2[i] = Ep[i];
        xyline2lab(numpoints, trl, Erl, Erl2, "Time (s)", "Pulsar Period (ms)",
                   "Velocity (km/s)", 1);
    }
    /* Calculate the phase delay */

    if (type == 'd' || type == 'a') {
        for (i = 0; i < numpoints; i++)
            Ep[i] = tmpE[i];
        E_to_phib(Ep, numpoints, &orb);
        for (i = 0; i < numpoints; i++)
            Erl[i] = Ep[i];
        xyline(numpoints, trl, Erl, "Time (s)", "TOA Delay (s)", 1);
    }
    /* Calculate the Fourier f-dot */

    for (i = 0; i < numpoints; i++)
        Ep[i] = tmpE[i];
    E_to_z(Ep, numpoints, p_psr, endtime, &orb);
    avgz = 0.0;
    for (i = 0; i < numpoints; i++) {
        Erl[i] = Ep[i];
        avgz += Ep[i];
    }
    avgz /= numpoints;
    if (type == 'z' || type == 'a') {
        xyline(numpoints, trl, Erl, "Time (s)", "Fourier F-dot (bins)", 1);
    }
    cpgend();

    printf("\n\nData used for orbit:\n");
    if (idata.N) {
        printf("  Starting Epoch (MJD)    =  %-15.9f\n", epoch);
    }
    printf("  Initial PSR period  (s) =  %-15.12f\n", p_psr);
    printf("  Initial Fourier bin     =  %-14.4f\n", endtime / p_psr);
    printf("  Orbital period      (s) =  %-13.4f\n", orb.p);
    printf("  eccentricity            =  %-8.6f\n", orb.e);
    printf("  a_sin(i)/c          (s) =  %-9.4f\n", orb.x);
    printf("  Angle of Peri     (deg) =  %-11.7f\n", orb.w);
    printf("  Time since Peri     (s) =  %-13.4f\n", orb.t);
    printf("  Length of Obs       (s) =  %-9.2f\n", endtime);
    printf("  Average 'z'      (bins) =  %-9.3f\n\n", avgz);

    vect_free(Erl);
    vect_free(Erl2);
    vect_free(trl);
    vect_free(tmpE);
    vect_free(Ep);
    exit(0);
}
