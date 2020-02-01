/*

Sample code:

 if(!getpoly(mjd, &dm, fppoly, pulsarname)) 
 {printf("Polycos not found!\n"); exit(1);}

 mjd is the MJD of the start of the observation (or middle)
 fppoly is the FILE *

 ....

 phcalc(mjd0, mjd1, &phase, &psrfreq);
 p = 1./psrfreq;
 fbin = dt*((float)nbins)/p;
 tbin = phase*((float)nbins);

 mjd0 is integer day
 mjd1 is the fractional day (topocentric UTC)
 tbin is current phase (of central frequency; you need to add dm offset)
 fbin is bin increment per sample.  
 
 Generally it's a good idea to
 recalculate phase and p every now and then, depending on how
 fast/relativistic the pulsar is.  I often do it every ~100 msec.

 from Ingrid Stairs

*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "presto.h"

#define MAX_POLYCOS 500

static double f0[MAX_POLYCOS], z4[MAX_POLYCOS], rphase[MAX_POLYCOS],
    mjdmid[MAX_POLYCOS], mjd1mid[MAX_POLYCOS], coeff[MAX_POLYCOS][15];
static int isets, nblk, ncoeff, icurr;

extern int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);

char *make_polycos(char *parfilenm, infodata * idata, char *polycofilenm, int debug_tempo)
{
    FILE *tmpfile;
    int tracklen;
    double T, fmid = 0.0, epoch;
    char *command, *psrname, scopechar;
    char *pcpathnm, *pcfilenm;
    psrparams psr;

    /* Get the path and name of the output polycofilenm */
    split_path_file(polycofilenm, &pcpathnm, &pcfilenm);

    /* Read the parfile */
    epoch = idata->mjd_i + idata->mjd_f;
    T = (idata->dt * idata->N) / SECPERDAY;
    if (!get_psr_from_parfile(parfilenm, epoch, &psr)) {
        fprintf(stderr,
                "\nError:  Cannot read parfile '%s' in make_polycos()\n\n",
                parfilenm);
        exit(-1);
    }

    /* Generate temp directory */
    char tmpdir[] = "/tmp/polycoXXXXXX";
    if (mkdtemp(tmpdir) == NULL) {
        fprintf(stderr,
                "\nError:  Cannot generate temp dir '%s' in make_polycos()\n\n",
                tmpdir);
        exit(-1);
    }
    if (debug_tempo) {
        fprintf(stderr, "Debugging TEMPO call:  Using temp directory '%s'\n",
                tmpdir);
    }

    
    /* Copy the parfile to the temp directory */
    command = (char *) calloc(strlen(parfilenm) + strlen(tmpdir) +
                              strlen(pcfilenm) + strlen(pcpathnm) + 200, 1);
    sprintf(command, "cp %s %s/pulsar.par", parfilenm, tmpdir);
    if (system(command) != 0) {
        fprintf(stderr,
                "\nError:  Cannot copy parfile '%s' to tmpdir '%s' in make_polycos()\n\n",
                parfilenm, tmpdir);
        exit(-1);
    }

    /* change to temp dir */
    char *origdir = getcwd(NULL, 0);
    chdir(tmpdir);

    /* Write tz.in */
    if (strcmp(idata->telescope, "GBT") == 0) {
        scopechar = '1';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "Arecibo") == 0) {
        scopechar = '3';
        tracklen = 3;
    } else if (strcmp(idata->telescope, "VLA") == 0) {
        scopechar = '6';
        tracklen = 6;
    } else if (strcmp(idata->telescope, "Parkes") == 0) {
        scopechar = '7';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "Jodrell") == 0) {
        scopechar = '8';
        tracklen = 12;
    } else if ((strcmp(idata->telescope, "GB43m") == 0) ||
               (strcmp(idata->telescope, "GB 140FT") == 0) ||
               (strcmp(idata->telescope, "NRAO20") == 0)) {
        scopechar = 'a';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "Nancay") == 0) {
        scopechar = 'f';
        tracklen = 4;
    } else if (strcmp(idata->telescope, "Effelsberg") == 0) {
        scopechar = 'g';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "ATA") == 0) {
        scopechar = 's';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "LOFAR") == 0) {
        scopechar = 't';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "WSRT") == 0) {
        scopechar = 'i';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "FAST") == 0) {
        scopechar = 'k';
        tracklen = 5;
    } else if (strcmp(idata->telescope, "GMRT") == 0) {
        scopechar = 'r';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "CHIME") == 0) {
        scopechar = 'y';
        tracklen = 1;
    } else if (strcmp(idata->telescope, "MWA") == 0) {
        scopechar = 'u';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "LWA") == 0) {
        scopechar = 'x';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "SRT") == 0) {
        scopechar = 'z';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "MeerKAT") == 0) {
        scopechar = 'm';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "KAT-7") == 0) {
        scopechar = 'k';
        tracklen = 12;
    } else if (strcmp(idata->telescope, "Geocenter") == 0) {
        scopechar = 'o';
        tracklen = 12;
    } else {                    /*  Barycenter */
        printf("Defaulting to barycenter for polyco generation...\n");
        scopechar = '@';
        tracklen = 12;
    }
    /* For optical, X-ray, or gamma-ray data */
    if (scopechar != '@' && scopechar != 'o') {
        fmid = idata->freq + (idata->num_chan / 2 - 0.5) * idata->chan_wid;
    } else {
        fmid = 0.0;
    }
    printf("Generating polycos for PSR %s.\n", psr.jname);
    tmpfile = chkfopen("tz.in", "w");
    fprintf(tmpfile, "%c %d 60 12 430\n\n\n%s 60 12 %d %.5f\n",
            scopechar, tracklen, psr.jname, tracklen, fmid);
    fclose(tmpfile);
    if (debug_tempo) {
        sprintf(command, "echo %d %d | tempo -z -f pulsar.par > tempo.out",
                idata->mjd_i - 1, (int) ceil(epoch + T));
        fprintf(stderr, "Debugging TEMPO call:  '%s'\n", command);
    } else {
        sprintf(command, "echo %d %d | tempo -z -f pulsar.par > /dev/null",
                idata->mjd_i - 1, (int) ceil(epoch + T));
    }
    if (system(command) != 0) {
        fprintf(stderr,
                "\nError:  Problem running TEMPO in '%s' for make_polycos()\n\n",
                tmpdir);
        exit(-1);
    } else {
        sprintf(command, "cp polyco.dat %s/%s", pcpathnm, pcfilenm);
        if (system(command) != 0) {
            fprintf(stderr,
                    "\nError:  Cannot copy polyco.dat in '%s' to '%s' make_polycos()\n\n",
                    tmpdir, polycofilenm);
            exit(-1);
        }
        if (!debug_tempo) {
            remove("polyco.dat");
            remove("pulsar.par");
            remove("tempo.lis");
            remove("tz.in");
            remove("tz.tmp");
        }
    }
    chdir(origdir);
    free(origdir);
    free(pcpathnm);
    free(pcfilenm);
    free(command);
    if (!debug_tempo) remove(tmpdir);
    psrname = (char *) calloc(strlen(psr.jname) + 1, sizeof(char));
    strcpy(psrname, psr.jname);
    return psrname;
}



int getpoly(double mjd, double duration, double *dm, FILE * fp, char *pname)
{

/*
  Given pname, mjd, and duration (days).  returns coeffs and dm
  Reads polynomial coeffs for current source and time from polyco.dat.
  DM and binary phase info added 23 Apr 90 (new version TZ)
  Puts out two sets of mjdmids for better precision. 16 Nov 1996.
  Recovers Earth-Doppler-shift factor and polyco frequency 13 Feb 1997.
*/

    char name0[15], testname[15], date0[15], binpha[16];
    char dummy[3][30];
    char buffer[160];
    double aphi0b[30], adphib[30];
    double dm0, z40, mjdend;
    float r, tfreq;
    long int mjddummy;
    int binary;

    int j, k, kk, len, jobs;
    int nblk0, ncoeff0;
    double mjdcheck;

    mjdend = mjd + duration;
    j = 0;
    while (fgets(buffer, 90, fp) != NULL) {
        sscanf(buffer, "%s", testname);
        if (strncmp(pname, testname, 4) == 0) {
            sscanf(buffer, "%s%s%f%ld%lf%lf%lf",
                   name0, date0, &r, &mjddummy, &mjd1mid[j], &dm0, &z40);
            fgets(buffer, 80, fp);
            sscanf(buffer, "%lf%lf%i%d%d", &rphase[j], &f0[j], &jobs, &nblk0,
                   &ncoeff0);
            fgets(buffer, 80, fp);
            sscanf(buffer, "%f%16c", &r, binpha);
            for (k = 0; k < ncoeff0 / 3; k++) {
                fgets(buffer, 80, fp);
                sscanf(buffer, "%s%s%s", dummy[0], dummy[1], dummy[2]);
                for (kk = 0; kk < 3; kk++) {
                    len = strlen(dummy[kk]);
                    if (dummy[kk][len - 4] == 'D')
                        dummy[kk][len - 4] = 'e';
                    sscanf(dummy[kk], "%lf", &coeff[j][3 * k + kk]);
                }
            }
            mjdmid[j] = mjddummy;

            /* just in case its a futmid we want an mjdmid */
            if (mjdmid[j] < 20000)
                mjdmid[j] += 39126.;
            z4[j] = z40;
            tfreq = r;
            *dm = dm0;
            binary = (binpha[0] != ' ');
            mjdcheck = mjdmid[j] + mjd1mid[j];
            // printf("mjd, mjdcheck, j: %lf %lf %d\n",mjd,mjdcheck,j);
            if (mjdcheck > mjd - 0.5 && mjdcheck < mjdend + 0.5) {
                nblk = nblk0;
                ncoeff = ncoeff0;
                if (binary)
                    sscanf(binpha, "%lf%lf", &aphi0b[j], &adphib[j]);
                if (ncoeff > 15) {
                    printf("ncoeff too big in polyco.dat.\n");
                    exit(1);
                }
                if (ncoeff < 15)
                    for (k = ncoeff; k < 15; k++)
                        coeff[j][k] = 0.;
                rphase[j] -= floor(rphase[j]);
                if ((rphase[j] < 0.) || (rphase[j] > 1.)) {
                    printf("rphase[%d] = %f\n", j, rphase[j]);
                    exit(1);
                }
                j++;
                if (j >= MAX_POLYCOS) {
                    fprintf(stderr,
                            "Too many polycos in polycos.c (MAX_POLYCO = %d)",
                            MAX_POLYCOS);
                    exit(32);
                }
            }
        }
    }

    isets = j;
    return (isets);

}


/*  Compute pulsar phase and frequency at time mjd0+mjd1. */

int phcalc(double mjd0, double mjd1, int last_index, double *phase, double *psrfreq)
{
    double dtmin;
    int i, j;
    int show;

    show = 0;

    icurr = -1;
    *psrfreq = f0[0];           /* Default psrfreq */
    for (j = last_index; j < isets; j++) {
        dtmin = (mjd0 - mjdmid[j]);
        dtmin = (dtmin + (mjd1 - mjd1mid[j])) * 1440.;  /* Time from center of this set */
        // printf("j, dtmin, comp = %d  %.16g  %.16g  %.16g\n", j, fabs(dtmin), nblk / 2.0 + 1e-10, fabs(dtmin)-(nblk / 2.0 + 1e-10));
        if (fabs(dtmin) < (nblk / 2.0 + 1e-7)) {        /* The extra bit avoids a subtle bug since
                                                           roundoff can cause fabs(dtmin)>(nblk/2) */
            // printf("j, mjds: %d %f %f %f %f\n",j,mjd0,mjd1,mjdmid[j],mjd1mid[j]);
            *psrfreq = 0.;      /* Compute psrfreq and phase from */
            *phase = coeff[j][ncoeff - 1];      /* the polynomial coeffs. */
            if (show)
                printf("phase = %21.15e   :%21.15e\n", *phase, coeff[j][ncoeff - 1]);
            for (i = ncoeff - 1; i > 0; --i) {
                *psrfreq = dtmin * (*psrfreq) + i * coeff[j][i];
                *phase = dtmin * (*phase) + coeff[j][i - 1];
                if (show)
                    printf("phase = %21.15e   :%21.15e\n", *phase, coeff[j][i - 1]);
            }

            *psrfreq = f0[j] + *psrfreq / 60.;  /* Add in the DC terms and scale */
            *phase += rphase[j] + dtmin * 60. * f0[j];
            if (show)
                printf("phase = %21.15e   f0: %21.15e\n", *phase, f0[j]);
            *phase -= floor(*phase);
            if ((*phase < 0.) || (*phase > 1.)) {
                printf("phase = %21.15f\n", *phase);
                exit(1);
            }
            icurr = j;
            break;
        }
    }
    if (icurr == -1) {
        printf("MJD %9.3f out of range (%9.3f to %9.3f)\n",
               (mjd0 + mjd1), mjdmid[0] - nblk / 2880.,
               mjdmid[isets - 1] + nblk / 2880.);
        *phase = -999.;
        printf("isets = %d\n", isets);
        exit(1);
    }
    return icurr;
}
