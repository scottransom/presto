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
#include <unistd.h>
#include "prepfold.h"


/* The polycos loaded by getpoly() and used by phcalc() live in the       */
/* default set below.  Programs folding several pulsars at once (e.g.     */
/* prepfold_multi) give each one its own polycoset and use the _set()     */
/* versions, which share no state, so different sets can be used from     */
/* different threads at the same time.                                    */
static polycoset default_set;

extern int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);

/* PRESTO telescope name (idata->telescope) to the tempo2 observatory */
/* name, the old TEMPO tz site code, and the track length (i.e. max   */
/* hour angle) used for polyco generation.  Note that the old TEMPO   */
/* codes for ATA ('s' = SHAO 65m!) and KAT-7 ('k' = FAST!) were       */
/* actually wrong -- the tempo2 names fix that.                       */
typedef struct polyco_scope {
    char name[12];              /* PRESTO telescope name */
    char t2name[12];            /* tempo2 observatory name */
    char t1code;                /* TEMPO tz.in site code */
    int maxha;                  /* Track length (hours) */
} polyco_scope;

static polyco_scope pscopes[] = {
    {"GBT", "gbt", '1', 12},
    {"Arecibo", "ao", '3', 3},
    {"VLA", "vla", '6', 6},
    {"Parkes", "pks", '7', 12},
    {"Jodrell", "jb", '8', 12},
    {"GB43m", "gb140", 'a', 12},
    {"GB 140FT", "gb140", 'a', 12},
    {"NRAO20", "gb140", 'a', 12},
    {"Nancay", "ncy", 'f', 4},
    {"Effelsberg", "eff", 'g', 12},
    {"ATA", "ata", 's', 12},
    {"LOFAR", "lofar", 't', 12},
    {"WSRT", "wsrt", 'i', 12},
    {"FAST", "fast", 'k', 5},
    {"GMRT", "gmrt", 'r', 12},
    {"CHIME", "chime", 'y', 1},
    {"MWA", "mwa", 'u', 12},
    {"LWA", "lwa1", 'x', 12},
    {"SRT", "srt", 'z', 12},
    {"MeerKAT", "meerkat", 'm', 12},
    {"KAT-7", "k7", 'k', 12},
    {"Geocenter", "coe", '0', 12},
    {"Barycenter", "@", '@', 12},   /* Also the default */
};

static polyco_scope *find_polyco_scope(char *telescope)
/* Return the site info for polyco generation, defaulting to the barycenter */
{
    int ii, numscopes = sizeof(pscopes) / sizeof(polyco_scope);

    for (ii = 0; ii < numscopes; ii++) {
        if (strcmp(telescope, pscopes[ii].name) == 0)
            return &pscopes[ii];
    }
    printf("Defaulting to barycenter for polyco generation...\n");
    return &pscopes[numscopes - 1];
}

char *make_polycos(char *parfilenm, infodata * idata, char *polycofilenm, int debug)
/* Generate polycos for the observation described by idata using the   */
/* parfile parfilenm, and copy them to polycofilenm.  Returns the      */
/* pulsar name from the parfile.  This runs "tempo2 -tempo1 -polyco",  */
/* which writes TEMPO1-format polycos, so the external 'tempo2'        */
/* executable (with its $TEMPO2 runtime directory, both available      */
/* from conda-forge) is required.  TEMPO itself is no longer used.     */
/* tempo2 runs in its own temporary directory (via "cd dir && ..." in  */
/* the shell, not chdir(), which would change the directory of every   */
/* thread in the process), so this is safe to call from threads.       */
{
    FILE *checkfile;
    double T, fmid = 0.0, epoch;
    char *command, *psrname, *tmpfilenm;
    psrparams psr;
    polyco_scope *site;
    size_t cmdlen;

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
    if (debug) {
        fprintf(stderr, "Debugging tempo2 call:  Using temp directory '%s'\n",
                tmpdir);
    }

    /* Copy the parfile to the temp directory */
    cmdlen = strlen(parfilenm) + strlen(polycofilenm) + 2 * strlen(tmpdir) + 300;
    command = (char *) calloc(cmdlen, 1);
    tmpfilenm = (char *) calloc(strlen(tmpdir) + 20, 1);
    snprintf(command, cmdlen, "cp '%s' %s/pulsar.par", parfilenm, tmpdir);
    if (system(command) != 0) {
        fprintf(stderr,
                "\nError:  Cannot copy parfile '%s' to tmpdir '%s' in make_polycos()\n\n",
                parfilenm, tmpdir);
        exit(-1);
    }

    /* Identify the observatory and observing band */
    site = find_polyco_scope(idata->telescope);
    if (strcmp(site->t2name, "@") != 0 && strcmp(site->t2name, "coe") != 0) {
        fmid = idata->freq + (idata->num_chan / 2 - 0.5) * idata->chan_wid;
    } else {
        fmid = 0.0;
    }
    printf("Generating polycos for PSR %s.\n", psr.jname);
    snprintf(command, cmdlen,
             "cd %s && tempo2 -tempo1 -f pulsar.par -polyco \"%d %d 60 12 %d %s %.5f\" %s",
             tmpdir, idata->mjd_i - 1, (int) ceil(epoch + T), site->maxha,
             site->t2name, fmid, debug ? "> tempo2.out 2>&1" : "> /dev/null 2>&1");
    if (debug)
        fprintf(stderr, "Debugging tempo2 call:  '%s'\n", command);
    if (system(command) != 0) {
        fprintf(stderr,
                "\nError:  Problem running tempo2 in '%s' for make_polycos().\n"
                "        Is 'tempo2' on your PATH, with $TEMPO2 pointing at its\n"
                "        runtime directory?  (Both are available on conda-forge.)\n\n",
                tmpdir);
        exit(-1);
    }
    /* tempo2 -polyco writes the TEMPO1-format polycos to polyco_new.dat */
    sprintf(tmpfilenm, "%s/polyco_new.dat", tmpdir);
    if ((checkfile = fopen(tmpfilenm, "r")) == NULL) {
        fprintf(stderr,
                "\nError:  tempo2 did not generate polycos in '%s' for make_polycos().\n"
                "        Re-run with -debug to see the tempo2 output.\n\n", tmpdir);
        exit(-1);
    }
    fclose(checkfile);
    snprintf(command, cmdlen, "cp %s '%s'", tmpfilenm, polycofilenm);
    if (system(command) != 0) {
        fprintf(stderr,
                "\nError:  Cannot copy polyco_new.dat in '%s' to '%s' in make_polycos()\n\n",
                tmpdir, polycofilenm);
        exit(-1);
    }
    if (!debug) {
        const char *tmpfiles[] = { "polyco_new.dat", "newpolyco.dat",
            "polyco.tim", "pulsar.par"
        };
        int ii;
        for (ii = 0; ii < 4; ii++) {
            sprintf(tmpfilenm, "%s/%s", tmpdir, tmpfiles[ii]);
            remove(tmpfilenm);
        }
        remove(tmpdir);
    }
    free(tmpfilenm);
    free(command);
    psrname = (char *) calloc(strlen(psr.jname) + 1, sizeof(char));
    strcpy(psrname, psr.jname);
    return psrname;
}


char *make_polycos_tempo1(char *parfilenm, infodata * idata, char *polycofilenm,
                          int debug_tempo)
/* The original TEMPO-based version of make_polycos(), retained only   */
/* so that the tempo2 implementation can be validated against it (see  */
/* check_polycos_t2.c).  Requires the external 'tempo' executable.     */
{
    FILE *tmpfile;
    int tracklen;
    double T, fmid = 0.0, epoch;
    char *command, *psrname, scopechar;
    char *pcpathnm, *pcfilenm;
    psrparams psr;
    polyco_scope *site;

    /* Get the path and name of the output polycofilenm */
    split_path_file(polycofilenm, &pcpathnm, &pcfilenm);

    /* Read the parfile */
    epoch = idata->mjd_i + idata->mjd_f;
    T = (idata->dt * idata->N) / SECPERDAY;
    if (!get_psr_from_parfile(parfilenm, epoch, &psr)) {
        fprintf(stderr,
                "\nError:  Cannot read parfile '%s' in make_polycos_tempo1()\n\n",
                parfilenm);
        exit(-1);
    }

    /* Generate temp directory */
    char tmpdir[] = "/tmp/polycoXXXXXX";
    if (mkdtemp(tmpdir) == NULL) {
        fprintf(stderr,
                "\nError:  Cannot generate temp dir '%s' in make_polycos_tempo1()\n\n",
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
                "\nError:  Cannot copy parfile '%s' to tmpdir '%s' in make_polycos_tempo1()\n\n",
                parfilenm, tmpdir);
        exit(-1);
    }

    /* change to temp dir */
    char *origdir = getcwd(NULL, 0);
    chdir(tmpdir);

    /* Write tz.in */
    site = find_polyco_scope(idata->telescope);
    scopechar = site->t1code;
    tracklen = site->maxha;
    /* For optical, X-ray, or gamma-ray data */
    if (scopechar != '@' && scopechar != '0') {
        fmid = idata->freq + (idata->num_chan / 2 - 0.5) * idata->chan_wid;
    } else {
        fmid = 0.0;
    }
    printf("Generating TEMPO polycos for PSR %s.\n", psr.jname);
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
                "\nError:  Problem running TEMPO in '%s' for make_polycos_tempo1()\n\n",
                tmpdir);
        exit(-1);
    } else {
        sprintf(command, "cp polyco.dat %s/%s", pcpathnm, pcfilenm);
        if (system(command) != 0) {
            fprintf(stderr,
                    "\nError:  Cannot copy polyco.dat in '%s' to '%s' make_polycos_tempo1()\n\n",
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
    if (!debug_tempo)
        remove(tmpdir);
    psrname = (char *) calloc(strlen(psr.jname) + 1, sizeof(char));
    strcpy(psrname, psr.jname);
    return psrname;
}


/* Make room for at least n polyco sets in set */
static void grow_polycoset(polycoset * set, int n)
{
    if (n <= set->nalloc)
        return;
    n = (set->nalloc == 0) ? 16 : 2 * set->nalloc;
    set->f0 = (double *) realloc(set->f0, n * sizeof(double));
    set->z4 = (double *) realloc(set->z4, n * sizeof(double));
    set->rphase = (double *) realloc(set->rphase, n * sizeof(double));
    set->mjdmid = (double *) realloc(set->mjdmid, n * sizeof(double));
    set->mjd1mid = (double *) realloc(set->mjd1mid, n * sizeof(double));
    set->coeff = (double (*)[15]) realloc(set->coeff, n * sizeof(*set->coeff));
    if (!set->f0 || !set->z4 || !set->rphase || !set->mjdmid ||
        !set->mjd1mid || !set->coeff) {
        fprintf(stderr, "\nError:  out of memory reading polycos.\n\n");
        exit(1);
    }
    set->nalloc = n;
}


void free_polycoset(polycoset * set)
/* Free the polycos held in set (but not set itself) */
{
    free(set->f0);
    free(set->z4);
    free(set->rphase);
    free(set->mjdmid);
    free(set->mjd1mid);
    free(set->coeff);
    memset(set, 0, sizeof(polycoset));
}


int getpoly(double mjd, double duration, double *dm, FILE * fp, char *pname)
/* getpoly_set() into the default polyco set used by phcalc() */
{
    return getpoly_set(&default_set, mjd, duration, dm, fp, pname);
}


int phcalc(double mjd0, double mjd1, int last_index, double *phase, double *psrfreq)
/* phcalc_set() with the default polyco set loaded by getpoly() */
{
    return phcalc_set(&default_set, mjd0, mjd1, last_index, phase, psrfreq);
}


int getpoly_set(polycoset * set, double mjd, double duration, double *dm,
                FILE * fp, char *pname)
/* A NULL set means the default set (as used by getpoly()) */
{

/*
  Given pname, mjd, and duration (days).  returns coeffs and dm
  Reads polynomial coeffs for current source and time from polyco.dat.
  DM and binary phase info added 23 Apr 90 (new version TZ)
  Puts out two sets of mjdmids for better precision. 16 Nov 1996.
  Recovers Earth-Doppler-shift factor and polyco frequency 13 Feb 1997.
*/

    char name0[15], testname[15], date0[15], binpha[17];
    char obsstr[16];
    char dummy[3][30];
    char buffer[160];
    double aphi0b, adphib;      /* binary phase info (unused) */
    double dm0, z40, mjdend;
    float r;
    long int mjddummy;
    int binary;

    int j, k, kk, len;
    int nblk0, ncoeff0;
    double mjdcheck;

    if (set == NULL)
        set = &default_set;
    mjdend = mjd + duration;
    j = 0;
    grow_polycoset(set, 1);
    while (fgets(buffer, 90, fp) != NULL) {
        sscanf(buffer, "%s", testname);
        if (strncmp(pname, testname, 4) == 0) {
            sscanf(buffer, "%s%s%f%ld%lf%lf%lf",
                   name0, date0, &r, &mjddummy, &set->mjd1mid[j], &dm0, &z40);
            fgets(buffer, 80, fp);
            /* The observatory field is a numeric TEMPO code from TEMPO */
            /* polycos but a site name (e.g. "gbt") from tempo2 ones,   */
            /* so read it as a string (it is unused regardless).        */
            sscanf(buffer, "%lf%lf%15s%d%d", &set->rphase[j], &set->f0[j], obsstr, &nblk0,
                   &ncoeff0);
            if (ncoeff0 > 15) { /* coeff[] holds at most 15 per set */
                printf("ncoeff too big in polyco.dat.\n");
                exit(1);
            }
            fgets(buffer, 80, fp);
            /* %16c writes no NUL, and nothing at all if the line is short */
            /* (as in tempo2 polycos), so start from a blank string.      */
            memset(binpha, ' ', 16);
            binpha[16] = '\0';
            sscanf(buffer, "%f%16c", &r, binpha);
            for (k = 0; k < ncoeff0 / 3; k++) {
                fgets(buffer, 80, fp);
                sscanf(buffer, "%s%s%s", dummy[0], dummy[1], dummy[2]);
                for (kk = 0; kk < 3; kk++) {
                    len = strlen(dummy[kk]);
                    if (dummy[kk][len - 4] == 'D')
                        dummy[kk][len - 4] = 'e';
                    sscanf(dummy[kk], "%lf", &set->coeff[j][3 * k + kk]);
                }
            }
            set->mjdmid[j] = mjddummy;

            /* just in case its a futmid we want an mjdmid */
            if (set->mjdmid[j] < 20000)
                set->mjdmid[j] += 39126.;
            set->z4[j] = z40;
            *dm = dm0;
            binary = (binpha[0] != ' ');
            mjdcheck = set->mjdmid[j] + set->mjd1mid[j];
            // printf("mjd, mjdcheck, j: %lf %lf %d\n",mjd,mjdcheck,j);
            if (mjdcheck > mjd - 0.5 && mjdcheck < mjdend + 0.5) {
                set->nblk = nblk0;
                set->ncoeff = ncoeff0;
                if (binary)
                    sscanf(binpha, "%lf%lf", &aphi0b, &adphib);
                if (set->ncoeff > 15) {
                    printf("ncoeff too big in polyco.dat.\n");
                    exit(1);
                }
                if (set->ncoeff < 15)
                    for (k = set->ncoeff; k < 15; k++)
                        set->coeff[j][k] = 0.;
                set->rphase[j] -= floor(set->rphase[j]);
                if ((set->rphase[j] < 0.) || (set->rphase[j] > 1.)) {
                    printf("rphase[%d] = %f\n", j, set->rphase[j]);
                    exit(1);
                }
                j++;
                grow_polycoset(set, j + 1);
            }
        }
    }

    set->isets = j;
    return (set->isets);

}


/*  Compute pulsar phase and frequency at time mjd0+mjd1. */

int phcalc_set(polycoset * set, double mjd0, double mjd1, int last_index,
               double *phase, double *psrfreq)
/* A NULL set means the default set (as used by phcalc()) */
{
    double dtmin;
    int i, j, icurr;
    int show;

    if (set == NULL)
        set = &default_set;

    show = 0;

    icurr = -1;
    *psrfreq = set->f0[0];           /* Default psrfreq */
    for (j = last_index; j < set->isets; j++) {
        dtmin = (mjd0 - set->mjdmid[j]);
        dtmin = (dtmin + (mjd1 - set->mjd1mid[j])) * 1440.;  /* Time from center of this set */
        // printf("j, dtmin, comp = %d  %.16g  %.16g  %.16g\n", j, fabs(dtmin), set->nblk / 2.0 + 1e-10, fabs(dtmin)-(set->nblk / 2.0 + 1e-10));
        if (fabs(dtmin) < (set->nblk / 2.0 + 1e-7)) {        /* The extra bit avoids a subtle bug since
                                                           roundoff can cause fabs(dtmin)>(set->nblk/2) */
            // printf("j, mjds: %d %f %f %f %f\n",j,mjd0,mjd1,set->mjdmid[j],set->mjd1mid[j]);
            *psrfreq = 0.;      /* Compute psrfreq and phase from */
            *phase = set->coeff[j][set->ncoeff - 1];      /* the polynomial coeffs. */
            if (show)
                printf("phase = %21.15e   :%21.15e\n", *phase, set->coeff[j][set->ncoeff - 1]);
            for (i = set->ncoeff - 1; i > 0; --i) {
                *psrfreq = dtmin * (*psrfreq) + i * set->coeff[j][i];
                *phase = dtmin * (*phase) + set->coeff[j][i - 1];
                if (show)
                    printf("phase = %21.15e   :%21.15e\n", *phase, set->coeff[j][i - 1]);
            }

            *psrfreq = set->f0[j] + *psrfreq / 60.;  /* Add in the DC terms and scale */
            *phase += set->rphase[j] + dtmin * 60. * set->f0[j];
            if (show)
                printf("phase = %21.15e   f0: %21.15e\n", *phase, set->f0[j]);
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
               (mjd0 + mjd1), set->mjdmid[0] - set->nblk / 2880.,
               set->mjdmid[set->isets - 1] + set->nblk / 2880.);
        *phase = -999.;
        printf("isets = %d\n", set->isets);
        exit(1);
    }
    return icurr;
}
