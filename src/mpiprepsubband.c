#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <errno.h>
#include "presto.h"
#include "mpiprepsubband_cmd.h"
#include "mask.h"
#include "backend_common.h"
#include "mpi.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAWDATA (cmd->filterbankP || cmd->psrfitsP)

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_LONG(x) (long) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

extern void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite);
extern void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite);
extern void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins, int downsamp);
extern void print_percent_complete(int current, int number);
extern void make_infodata_struct(void);
extern void make_maskbase_struct(void);
extern void make_spectra_info_struct(void);
extern void broadcast_mask(mask * obsmask, int myid);
extern void broadcast_spectra_info(struct spectra_info *s, int myid);
extern void print_dms(char *hostname, int myid, int numprocs,
                      int local_numdms, double *dms);
static int get_data(float **outdata, int blocksperread,
                    struct spectra_info *s,
                    mask * obsmask, int *idispdts, int **offsets, int *padding);
extern void set_using_MPI(void);

MPI_Datatype maskbase_type;
MPI_Datatype spectra_info_type;

static Cmdline *cmd;
static int blocksperread = 0, worklen = 0;
static int insubs = 0;
static int local_numdms = 1, myid = 0, numprocs = 1;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
    /* Any variable that begins with 't' means topocentric */
    /* Any variable that begins with 'b' means barycentric */
    FILE **outfiles = NULL;
    float **outdata;
    double dtmp, *dms, avgdm = 0.0, dsdt = 0, maxdm;
    double *dispdt, tlotoa = 0.0, blotoa = 0.0, BW_ddelay = 0.0;
    double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
    double *btoa = NULL, *ttoa = NULL, avgvoverc = 0.0;
    char obs[3], ephem[10], rastring[50], decstring[50];
    long totnumtowrite, totwrote = 0, padwrote = 0, datawrote = 0;
    int *idispdt, **offsets;
    int ii, jj, numadded = 0, numremoved = 0, padding = 0, good_inputs = 1;
    int numbarypts = 0, numread = 0, numtowrite = 0;
    int padtowrite = 0, statnum = 0;
    int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL, good_padvals = 0;
    double local_lodm;
    char *datafilenm, *outpath, *outfilenm, hostname[256];
    struct spectra_info s;
    infodata idata;
    mask obsmask;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#ifdef _OPENMP
    omp_set_num_threads(1);     // Explicitly turn off OpenMP
#endif
    set_using_MPI();

    /* Get hostname on Unix machine */
    {
        int retval = gethostname(hostname, 255);
        if (retval == -1) {
            printf("Warning:  error determining hostname: %s\n", strerror(errno));
            sprintf(hostname, "unknown");
        }
    }

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        if (myid == 0) {
            Program = argv[0];
            usage();
        }
        MPI_Finalize();
        exit(1);
    }

    make_maskbase_struct();
    make_spectra_info_struct();

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);
    spectra_info_set_defaults(&s);
    // If we are zeroDMing, make sure that clipping is off.
    if (cmd->zerodmP)
        cmd->noclipP = 1;
    s.clip_sigma = cmd->clip;
    if (cmd->noclipP) {
        cmd->clip = 0.0;
        s.clip_sigma = 0.0;
    }
    if (cmd->ifsP) {
        // 0 = default or summed, 1-4 are possible also
        s.use_poln = cmd->ifs + 1;
    }
    if (!cmd->numoutP)
        cmd->numout = LONG_MAX;

#ifdef DEBUG
    showOptionValues();
#endif

    if (myid == 0) {            /* Master node only */
        printf("\n\n");
        printf("      Parallel Pulsar Subband De-dispersion Routine\n");
        printf("                 by Scott M. Ransom\n\n");

        s.filenames = cmd->argv;
        s.num_files = cmd->argc;
        s.clip_sigma = cmd->clip;
        // -1 causes the data to determine if we use weights, scales, & 
        // offsets for PSRFITS or flip the band for any data type where
        // we can figure that out with the data
        s.apply_flipband = (cmd->invertP) ? 1 : -1;
        s.apply_weight = (cmd->noweightsP) ? 0 : -1;
        s.apply_scale = (cmd->noscalesP) ? 0 : -1;
        s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
        s.remove_zerodm = (cmd->zerodmP) ? 1 : 0;

        if (RAWDATA) {
            if (cmd->filterbankP)
                s.datatype = SIGPROCFB;
            else if (cmd->psrfitsP)
                s.datatype = PSRFITS;
        } else {                // Attempt to auto-identify the data
            identify_psrdatatype(&s, 1);
            if (s.datatype == SIGPROCFB)
                cmd->filterbankP = 1;
            else if (s.datatype == PSRFITS)
                cmd->psrfitsP = 1;
            else if (s.datatype == SUBBAND)
                insubs = 1;
            else {
                printf
                    ("\nError:  Unable to identify input data files.  Please specify type.\n\n");
                good_inputs = 0;
            }
        }
        // So far we can only handle PSRFITS, filterbank, and subbands
        if (s.datatype != PSRFITS &&
            s.datatype != SIGPROCFB && s.datatype != SUBBAND)
            good_inputs = 0;

        // For subbanded data
        if (!RAWDATA)
            s.files = (FILE **) malloc(sizeof(FILE *) * s.num_files);

        if (good_inputs && (RAWDATA || insubs)) {
            char description[40];
            psrdatatype_description(description, s.datatype);
            if (s.num_files > 1)
                printf("Reading %s data from %d files:\n", description, s.num_files);
            else
                printf("Reading %s data from 1 file:\n", description);
            for (ii = 0; ii < s.num_files; ii++) {
                printf("  '%s'\n", cmd->argv[ii]);
                if (insubs)
                    s.files[ii] = chkfopen(s.filenames[ii], "rb");
            }
            printf("\n");
            if (RAWDATA) {
                read_rawdata_files(&s);
                // Make sure that the requested number of subbands divides into the
                // the raw number of channels.
                if (s.num_channels % cmd->nsub) {
                    printf("Error:  The number of subbands (-nsub %d) must divide into the\n"
                           "        number of channels (%d)\n\n",
                           cmd->nsub, s.num_channels);
                    exit(1);
                }
                if (cmd->ignorechanstrP) {
                    s.ignorechans = get_ignorechans(cmd->ignorechanstr, 0, s.num_channels-1,
                                                    &s.num_ignorechans, &s.ignorechans_str);
                    if (s.ignorechans_str==NULL) {
                        s.ignorechans_str = (char *)malloc(strlen(cmd->ignorechanstr)+1);
                        strcpy(s.ignorechans_str, cmd->ignorechanstr);
                    }
                }
                print_spectra_info_summary(&s);
                spectra_info_to_inf(&s, &idata);
            } else {            // insubs
                char *root, *suffix;
                cmd->nsub = s.num_files;
                s.N = chkfilelen(s.files[0], sizeof(short));
                s.start_subint = gen_ivect(1);
                s.num_subint = gen_ivect(1);
                s.start_MJD = (long double *) malloc(sizeof(long double));
                s.start_spec = (long long *) malloc(sizeof(long long));
                s.num_spec = (long long *) malloc(sizeof(long long));
                s.num_pad = (long long *) malloc(sizeof(long long));
                s.start_spec[0] = 0L;
                s.start_subint[0] = 0;
                s.num_spec[0] = s.N;
                s.num_subint[0] = s.N / SUBSBLOCKLEN;
                s.num_pad[0] = 0L;
                s.padvals = gen_fvect(s.num_files);
                for (ii = 0; ii < s.num_files; ii++)
                    s.padvals[ii] = 0.0;
                if (split_root_suffix(s.filenames[0], &root, &suffix) == 0) {
                    printf
                        ("\nError:  The input filename (%s) must have a suffix!\n\n",
                         s.filenames[0]);
                    exit(1);
                }
                if (strncmp(suffix, "sub", 3) == 0) {
                    char *tmpname;
                    tmpname = calloc(strlen(root) + 10, 1);
                    sprintf(tmpname, "%s.sub", root);
                    readinf(&idata, tmpname);
                    free(tmpname);
                    // Both are length 40
                    strcpy(s.telescope, idata.telescope);
                    // Both are length 100
                    strcpy(s.backend, idata.instrument);
                    // Both are length 100
                    strcpy(s.observer, idata.observer);
                    // Both are length 100
                    strcpy(s.source, idata.object);
                    s.ra2000 = hms2rad(idata.ra_h, idata.ra_m,
                                       idata.ra_s) * RADTODEG;
                    s.dec2000 = dms2rad(idata.dec_d, idata.dec_m,
                                        idata.dec_s) * RADTODEG;
                    ra_dec_to_string(s.ra_str, idata.ra_h, idata.ra_m, idata.ra_s);
                    ra_dec_to_string(s.dec_str,
                                     idata.dec_d, idata.dec_m, idata.dec_s);
                    s.num_channels = idata.num_chan;
                    s.start_MJD[0] = idata.mjd_i + idata.mjd_f;
                    s.dt = idata.dt;
                    s.T = s.N * s.dt;
                    s.lo_freq = idata.freq;
                    s.df = idata.chan_wid;
                    s.hi_freq = s.lo_freq + (s.num_channels - 1.0) * s.df;
                    s.BW = s.num_channels * s.df;
                    s.fctr = s.lo_freq - 0.5 * s.df + 0.5 * s.BW;
                    s.beam_FWHM = idata.fov / 3600.0;
                    s.spectra_per_subint = SUBSBLOCKLEN;
                    print_spectra_info_summary(&s);
                } else {
                    printf
                        ("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                         cmd->argv[0]);
                    MPI_Finalize();
                    exit(1);
                }
                free(root);
                free(suffix);
            }
        }
    }
    //  If we don't have good input data, exit
    MPI_Bcast(&good_inputs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!good_inputs) {
        MPI_Finalize();
        exit(1);
    }

    MPI_Bcast(&insubs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (insubs)
        cmd->nsub = cmd->argc;

    /* Determine the output file names and open them */

    local_numdms = cmd->numdms / (numprocs - 1);
    dms = gen_dvect(local_numdms);
    if (cmd->numdms % (numprocs - 1)) {
        if (myid == 0)
            printf
                ("\nThe number of DMs (%d) must be divisible by (the number of processors (%d) - 1).\n\n", cmd->numdms, numprocs);
        MPI_Finalize();
        exit(1);
    }
    local_lodm = cmd->lodm + (myid - 1) * local_numdms * cmd->dmstep;

    split_path_file(cmd->outfile, &outpath, &outfilenm);
    datafilenm = (char *) calloc(strlen(outfilenm) + 20, 1);
    if (myid > 0) {
        if (chdir(outpath) == -1) {
            printf("\nProcess %d on %s cannot chdir() to '%s'.  Exiting.\n\n",
                   myid, hostname, outpath);
            MPI_Finalize();
            exit(1);
        }
        outfiles = (FILE **) malloc(local_numdms * sizeof(FILE *));
        for (ii = 0; ii < local_numdms; ii++) {
            dms[ii] = local_lodm + ii * cmd->dmstep;
            avgdm += dms[ii];
            sprintf(datafilenm, "%s_DM%.2f.dat", outfilenm, dms[ii]);
            outfiles[ii] = chkfopen(datafilenm, "wb");
        }
        avgdm /= local_numdms;
    }
    // Broadcast the raw data information

    broadcast_spectra_info(&s, myid);
    if (myid > 0) {
        spectra_info_to_inf(&s, &idata);
        if (s.datatype == SIGPROCFB)
            cmd->filterbankP = 1;
        else if (s.datatype == PSRFITS)
            cmd->psrfitsP = 1;
        else if (s.datatype == SUBBAND)
            insubs = 1;
    }
    s.filenames = cmd->argv;

    /* Read an input mask if wanted */

    if (myid > 0) {
        int numpad = s.num_channels;
        if (insubs)
            numpad = s.num_files;
        s.padvals = gen_fvect(numpad);
        for (ii = 0; ii < numpad; ii++)
            s.padvals[ii] = 0.0;
    }
    if (cmd->maskfileP) {
        if (myid == 0) {
            read_mask(cmd->maskfile, &obsmask);
            printf("Read mask information from '%s'\n\n", cmd->maskfile);
            good_padvals = determine_padvals(cmd->maskfile, &obsmask, s.padvals);
        }
        broadcast_mask(&obsmask, myid);
        MPI_Bcast(&good_padvals, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(s.padvals, obsmask.numchan, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else {
        obsmask.numchan = obsmask.numint = 0;
        MPI_Bcast(&good_padvals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // The number of topo to bary time points to generate with TEMPO
    numbarypts = (int) (s.T * 1.1 / TDT + 5.5) + 1;

    // Identify the TEMPO observatory code
    {
        char *outscope = (char *) calloc(40, sizeof(char));
        telescope_to_tempocode(idata.telescope, outscope, obs);
        free(outscope);
    }

    // Broadcast or calculate a few extra important values
    if (insubs)
        avgdm = idata.dm;
    idata.dm = avgdm;
    dsdt = cmd->downsamp * idata.dt;
    maxdm = cmd->lodm + cmd->numdms * cmd->dmstep;
    BW_ddelay = delay_from_dm(maxdm, idata.freq) -
        delay_from_dm(maxdm, idata.freq + (idata.num_chan - 1) * idata.chan_wid);
    blocksperread = ((int) (BW_ddelay / idata.dt) / s.spectra_per_subint + 1);
    worklen = s.spectra_per_subint * blocksperread;

    /* If we are offsetting into the file, change inf file start time */
    if (cmd->start > 0.0 || cmd->offset > 0) {
        if (cmd->start > 0.0) /* Offset in units of worklen */
            cmd->offset = (long) (cmd->start *
                                  idata.N / worklen) * worklen;
        add_to_inf_epoch(&idata, cmd->offset * idata.dt);
        if (myid==0) {
            printf("Offsetting into the input files by %ld spectra (%.6g sec)\n",
                   cmd->offset, cmd->offset * idata.dt);
            if (RAWDATA)
                offset_to_spectra(cmd->offset, &s);
            else { // subbands
                for (ii = 0; ii < s.num_files; ii++)
                    chkfileseek(s.files[ii], cmd->offset, sizeof(short), SEEK_SET);
                if (cmd->maskfileP)
                    printf("WARNING!:  masking does not work with old-style subbands and -start or -offset!\n");
            }
        } else { // This is needed for masking on the nodes
            if (RAWDATA)
                set_currentspectra(cmd->offset);
        }
    }

    if (cmd->nsub > s.num_channels) {
        printf
            ("Warning:  The number of requested subbands (%d) is larger than the number of channels (%d).\n",
             cmd->nsub, s.num_channels);
        printf("          Re-setting the number of subbands to %d.\n\n",
               s.num_channels);
        cmd->nsub = s.num_channels;
    }

    if (s.spectra_per_subint % cmd->downsamp) {
        if (myid == 0) {
            printf
                ("\nError:  The downsample factor (%d) must be a factor of the\n",
                 cmd->downsamp);
            printf("        blocklength (%d).  Exiting.\n\n", s.spectra_per_subint);
        }
        MPI_Finalize();
        exit(1);
    }

    tlotoa = idata.mjd_i + idata.mjd_f; /* Topocentric epoch */

    /* Set the output length to a good number if it wasn't requested */
    if (!cmd->numoutP) {
        cmd->numoutP = 1;
        cmd->numout = choose_good_N((long long)(idata.N/cmd->downsamp));
        printf("Setting a 'good' output length of %ld samples\n", cmd->numout);
    }
    totnumtowrite = cmd->numout;

    if (cmd->nobaryP) {         /* Main loop if we are not barycentering... */

        /* Dispersion delays (in bins).  The high freq gets no delay   */
        /* All other delays are positive fractions of bin length (dt)  */

        dispdt = subband_search_delays(s.num_channels, cmd->nsub, avgdm,
                                       idata.freq, idata.chan_wid, 0.0);
        idispdt = gen_ivect(s.num_channels);
        for (ii = 0; ii < s.num_channels; ii++)
            idispdt[ii] = NEAREST_LONG(dispdt[ii] / idata.dt);
        vect_free(dispdt);

        /* The subband dispersion delays (see note above) */

        offsets = gen_imatrix(local_numdms, cmd->nsub);
        for (ii = 0; ii < local_numdms; ii++) {
            double *subdispdt;

            subdispdt = subband_delays(s.num_channels, cmd->nsub, dms[ii],
                                       idata.freq, idata.chan_wid, 0.0);
            dtmp = subdispdt[cmd->nsub - 1];
            for (jj = 0; jj < cmd->nsub; jj++)
                offsets[ii][jj] = NEAREST_LONG((subdispdt[jj] - dtmp) / dsdt);
            vect_free(subdispdt);
        }

        /* Allocate our data array and start getting data */

        if (myid == 0) {
            printf("De-dispersing using %d subbands.\n", cmd->nsub);
            if (cmd->downsamp > 1)
                printf("Downsampling by a factor of %d (new dt = %.10g)\n",
                       cmd->downsamp, dsdt);
            printf("\n");
        }

        /* Print the nodes and the DMs they are handling */
        print_dms(hostname, myid, numprocs, local_numdms, dms);

        outdata = gen_fmatrix(local_numdms, worklen / cmd->downsamp);
        numread = get_data(outdata, blocksperread, &s,
                           &obsmask, idispdt, offsets, &padding);

        while (numread == worklen) {

            numread /= cmd->downsamp;
            if (myid == 0)
                print_percent_complete(totwrote, totnumtowrite);

            /* Write the latest chunk of data, but don't   */
            /* write more than cmd->numout points.         */

            numtowrite = numread;
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            if (myid > 0) {
                write_data(outfiles, local_numdms, outdata, 0, numtowrite);
                /* Update the statistics */
                if (!padding) {
                    for (ii = 0; ii < numtowrite; ii++)
                        update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg,
                                     &var);
                    statnum += numtowrite;
                }
            }
            totwrote += numtowrite;

            /* Stop if we have written out all the data we need to */

            if (totwrote == cmd->numout)
                break;

            numread = get_data(outdata, blocksperread, &s,
                               &obsmask, idispdt, offsets, &padding);
        }
        datawrote = totwrote;

    } else {                    /* Main loop if we are barycentering... */

        /* What ephemeris will we use?  (Default is DE405) */
        strcpy(ephem, "DE405");

        /* Define the RA and DEC of the observation */

        ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
        ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

        /* Allocate some arrays */

        btoa = gen_dvect(numbarypts);
        ttoa = gen_dvect(numbarypts);
        for (ii = 0; ii < numbarypts; ii++)
            ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

        /* Call TEMPO for the barycentering */

        if (myid == 0) {
            double maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;

            printf("\nGenerating barycentric corrections...\n");
            voverc = gen_dvect(numbarypts);
            barycenter(ttoa, btoa, voverc, numbarypts, rastring, decstring, obs,
                       ephem);
            for (ii = 0; ii < numbarypts; ii++) {
                if (voverc[ii] > maxvoverc)
                    maxvoverc = voverc[ii];
                if (voverc[ii] < minvoverc)
                    minvoverc = voverc[ii];
                avgvoverc += voverc[ii];
            }
            avgvoverc /= numbarypts;
            vect_free(voverc);

            printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
            printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
            printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
            printf("De-dispersing using %d subbands.\n", cmd->nsub);
            if (cmd->downsamp > 1) {
                printf("     Downsample = %d\n", cmd->downsamp);
                printf("  New sample dt = %.10g\n", dsdt);
            }
            printf("\n");
        }

        /* Print the nodes and the DMs they are handling */
        print_dms(hostname, myid, numprocs, local_numdms, dms);

        MPI_Bcast(btoa, numbarypts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&avgvoverc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        blotoa = btoa[0];

        /* Dispersion delays (in bins).  The high freq gets no delay   */
        /* All other delays are positive fractions of bin length (dt)  */

        dispdt = subband_search_delays(s.num_channels, cmd->nsub, avgdm,
                                       idata.freq, idata.chan_wid, avgvoverc);
        idispdt = gen_ivect(s.num_channels);
        for (ii = 0; ii < s.num_channels; ii++)
            idispdt[ii] = NEAREST_LONG(dispdt[ii] / idata.dt);
        vect_free(dispdt);

        /* The subband dispersion delays (see note above) */

        offsets = gen_imatrix(local_numdms, cmd->nsub);
        for (ii = 0; ii < local_numdms; ii++) {
            double *subdispdt;

            subdispdt = subband_delays(s.num_channels, cmd->nsub, dms[ii],
                                       idata.freq, idata.chan_wid, avgvoverc);
            dtmp = subdispdt[cmd->nsub - 1];
            for (jj = 0; jj < cmd->nsub; jj++)
                offsets[ii][jj] = NEAREST_LONG((subdispdt[jj] - dtmp) / dsdt);
            vect_free(subdispdt);
        }

        /* Convert the bary TOAs to differences from the topo TOAs in */
        /* units of bin length (dt) rounded to the nearest integer.   */

        dtmp = (btoa[0] - ttoa[0]);
        for (ii = 0; ii < numbarypts; ii++)
            btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

        /* Find the points where we need to add or remove bins */
        {
            int oldbin = 0, currentbin;
            double lobin, hibin, calcpt;

            numdiffbins = labs(NEAREST_LONG(btoa[numbarypts - 1])) + 1;
            diffbins = gen_ivect(numdiffbins);
            diffbinptr = diffbins;
            for (ii = 1; ii < numbarypts; ii++) {
                currentbin = NEAREST_LONG(btoa[ii]);
                if (currentbin != oldbin) {
                    if (currentbin > 0) {
                        calcpt = oldbin + 0.5;
                        lobin = (ii - 1) * TDT / dsdt;
                        hibin = ii * TDT / dsdt;
                    } else {
                        calcpt = oldbin - 0.5;
                        lobin = -((ii - 1) * TDT / dsdt);
                        hibin = -(ii * TDT / dsdt);
                    }
                    while (fabs(calcpt) < fabs(btoa[ii])) {
                        /* Negative bin number means remove that bin */
                        /* Positive bin number means add a bin there */
                        *diffbinptr =
                            NEAREST_LONG(LININTERP
                                         (calcpt, btoa[ii - 1], btoa[ii], lobin,
                                          hibin));
                        diffbinptr++;
                        calcpt = (currentbin > 0) ? calcpt + 1.0 : calcpt - 1.0;
                    }
                    oldbin = currentbin;
                }
            }
            *diffbinptr = cmd->numout;  /* Used as a marker */
        }
        diffbinptr = diffbins;

        /* Now perform the barycentering */

        outdata = gen_fmatrix(local_numdms, worklen / cmd->downsamp);
        numread = get_data(outdata, blocksperread, &s,
                           &obsmask, idispdt, offsets, &padding);

        while (numread == worklen) {    /* Loop to read and write the data */
            int numwritten = 0;
            double block_avg, block_var;

            numread /= cmd->downsamp;
            /* Determine the approximate local average */
            avg_var(outdata[0], numread, &block_avg, &block_var);
            if (myid == 0)
                print_percent_complete(totwrote, totnumtowrite);

            /* Simply write the data if we don't have to add or */
            /* remove any bins from this batch.                 */
            /* OR write the amount of data up to cmd->numout or */
            /* the next bin that will be added or removed.      */

            numtowrite = abs(*diffbinptr) - datawrote;
            if ((totwrote + numtowrite) > cmd->numout)
                numtowrite = cmd->numout - totwrote;
            if (numtowrite > numread)
                numtowrite = numread;
            if (myid > 0) {
                write_data(outfiles, local_numdms, outdata, 0, numtowrite);
                /* Update the statistics */
                if (!padding) {
                    for (ii = 0; ii < numtowrite; ii++)
                        update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg,
                                     &var);
                    statnum += numtowrite;
                }
            }
            datawrote += numtowrite;
            totwrote += numtowrite;
            numwritten += numtowrite;

            if ((datawrote == abs(*diffbinptr)) && (numwritten != numread) && (totwrote < cmd->numout)) {       /* Add/remove a bin */
                int skip, nextdiffbin;

                skip = numtowrite;

                /* Write the rest of the data after adding/removing a bin  */
                do {

                    if (*diffbinptr > 0) {
                        /* Add a bin */
                        if (myid > 0)
                            write_padding(outfiles, local_numdms, block_avg, 1);
                        numadded++;
                        totwrote++;
                    } else {
                        /* Remove a bin */
                        numremoved++;
                        datawrote++;
                        numwritten++;
                        skip++;
                    }
                    diffbinptr++;

                    /* Write the part after the diffbin */

                    numtowrite = numread - numwritten;
                    if ((totwrote + numtowrite) > cmd->numout)
                        numtowrite = cmd->numout - totwrote;
                    nextdiffbin = abs(*diffbinptr) - datawrote;
                    if (numtowrite > nextdiffbin)
                        numtowrite = nextdiffbin;
                    if (myid > 0) {
                        write_data(outfiles, local_numdms, outdata, skip,
                                   numtowrite);
                        /* Update the statistics and counters */
                        if (!padding) {
                            for (ii = 0; ii < numtowrite; ii++)
                                update_stats(statnum + ii,
                                             outdata[0][skip + ii], &min, &max, &avg,
                                             &var);
                            statnum += numtowrite;
                        }
                    }
                    numwritten += numtowrite;
                    datawrote += numtowrite;
                    totwrote += numtowrite;
                    skip += numtowrite;

                    /* Stop if we have written out all the data we need to */

                    if (totwrote == cmd->numout)
                        break;
                } while (numwritten < numread);
            }
            /* Stop if we have written out all the data we need to */

            if (totwrote == cmd->numout)
                break;

            numread = get_data(outdata, blocksperread, &s,
                               &obsmask, idispdt, offsets, &padding);
        }
    }

    if (myid > 0) {

        /* Calculate the amount of padding we need  */

        if (cmd->numout > totwrote)
            padwrote = padtowrite = cmd->numout - totwrote;

        /* Write the new info file for the output data */

        idata.dt = dsdt;
        update_infodata(&idata, totwrote, padtowrite, diffbins,
                        numdiffbins, cmd->downsamp);
        for (ii = 0; ii < local_numdms; ii++) {
            idata.dm = dms[ii];
            if (!cmd->nobaryP) {
                double baryepoch, barydispdt, baryhifreq;

                baryhifreq = idata.freq + (s.num_channels - 1) * idata.chan_wid;
                barydispdt = delay_from_dm(dms[ii], doppler(baryhifreq, avgvoverc));
                baryepoch = blotoa - (barydispdt / SECPERDAY);
                idata.bary = 1;
                idata.mjd_i = (int) floor(baryepoch);
                idata.mjd_f = baryepoch - idata.mjd_i;
            }
            sprintf(idata.name, "%s_DM%.2f", outfilenm, dms[ii]);
            writeinf(&idata);
        }

        /* Set the padded points equal to the average data point */

        if (idata.numonoff >= 1) {
            int index, startpad, endpad;

            for (ii = 0; ii < local_numdms; ii++) {
                fclose(outfiles[ii]);
                sprintf(datafilenm, "%s_DM%.2f.dat", outfilenm, dms[ii]);
                outfiles[ii] = chkfopen(datafilenm, "rb+");
            }
            for (ii = 0; ii < idata.numonoff; ii++) {
                index = 2 * ii;
                startpad = idata.onoff[index + 1];
                if (ii == idata.numonoff - 1)
                    endpad = idata.N - 1;
                else
                    endpad = idata.onoff[index + 2];
                for (jj = 0; jj < local_numdms; jj++)
                    chkfseek(outfiles[jj], (startpad + 1) * sizeof(float), SEEK_SET);
                padtowrite = endpad - startpad;
                write_padding(outfiles, local_numdms, avg, padtowrite);
            }
        }
    }

    /* Print simple stats and results */

    var /= (datawrote - 1);
    if (myid == 0)
        print_percent_complete(1, 1);
    if (myid == 1) {
        printf("\n\nDone.\n\nSimple statistics of the output data:\n");
        printf("             Data points written:  %ld\n", totwrote);
        if (padwrote)
            printf("          Padding points written:  %ld\n", padwrote);
        if (!cmd->nobaryP) {
            if (numadded)
                printf("    Bins added for barycentering:  %d\n", numadded);
            if (numremoved)
                printf("  Bins removed for barycentering:  %d\n", numremoved);
        }
        printf("           Maximum value of data:  %.2f\n", max);
        printf("           Minimum value of data:  %.2f\n", min);
        printf("              Data average value:  %.2f\n", avg);
        printf("         Data standard deviation:  %.2f\n", sqrt(var));
        printf("\n");
    }

    /* Close the files and cleanup */

    if (cmd->maskfileP)
        free_mask(obsmask);
    if (myid > 0) {
        for (ii = 0; ii < local_numdms; ii++)
            fclose(outfiles[ii]);
        free(outfiles);
    }
    vect_free(outdata[0]);
    vect_free(outdata);
    vect_free(dms);
    vect_free(idispdt);
    vect_free(offsets[0]);
    vect_free(offsets);
    free(datafilenm);
    free(outfilenm);
    free(outpath);
    if (!cmd->nobaryP) {
        vect_free(btoa);
        vect_free(ttoa);
        vect_free(diffbins);
    }
    MPI_Finalize();
    return (0);
}


static int simple_read_subbands(FILE * infiles[], int numfiles, short *subbanddata)
/* Read short int subband data written by prepsubband */
{
    int ii, numread = 0;

    /* Read the data */
    for (ii = 0; ii < numfiles; ii++) {
        numread = chkfread(subbanddata + ii * SUBSBLOCKLEN,
                           sizeof(short), SUBSBLOCKLEN, infiles[ii]);
    }
    return numread;
}

static void convert_subbands(int numfiles, short *shortdata,
                             float *subbanddata, double timeperblk,
                             int *maskchans, int *nummasked, mask * obsmask,
                             float clip_sigma, float *padvals)
/* Convert and transpose the subband data, then mask it*/
{
    int ii, jj, index, shortindex, mask = 0, offset;
    double starttime, run_avg;
    float subband_sum;
    static int currentblock = 0;

    *nummasked = 0;
    if (obsmask->numchan)
        mask = 1;

    for (ii = 0; ii < numfiles; ii++) {
        shortindex = ii * SUBSBLOCKLEN;
        run_avg = 0.0;
        if (cmd->runavgP == 1) {
            for (jj = 0; jj < SUBSBLOCKLEN; jj++)
                run_avg += (float) shortdata[shortindex + jj];
            run_avg /= SUBSBLOCKLEN;
        }
        for (jj = 0, index = ii; jj < SUBSBLOCKLEN; jj++, index += numfiles)
            subbanddata[index] = (float) shortdata[shortindex + jj] - run_avg;
    }

    if (mask) {
        starttime = currentblock * timeperblk;
        *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
    }

    /* Clip nasty RFI if requested and we're not masking all the channels */
    if ((clip_sigma > 0.0) && !(mask && (*nummasked == -1))) {
        clip_times(subbanddata, SUBSBLOCKLEN, numfiles, clip_sigma, padvals);
    }

    /* Mask it if required */
    if (mask) {
        if (*nummasked == -1) { /* If all channels are masked */
            for (ii = 0; ii < SUBSBLOCKLEN; ii++)
                memcpy(subbanddata + ii * numfiles, padvals,
                       sizeof(float) * numfiles);
        } else if (*nummasked > 0) {    /* Only some of the channels are masked */
            int channum;
            for (ii = 0; ii < SUBSBLOCKLEN; ii++) {
                offset = ii * numfiles;
                for (jj = 0; jj < *nummasked; jj++) {
                    channum = maskchans[jj];
                    subbanddata[offset + channum] = padvals[channum];
                }
            }
        }
    }

    /* Zero-DM removal if required */
    if (cmd->zerodmP == 1) {
        for (ii = 0; ii < SUBSBLOCKLEN; ii++) {
            offset = ii * numfiles;
            subband_sum = 0.0;
            for (jj = offset; jj < offset + numfiles; jj++) {
                subband_sum += subbanddata[jj];
            }
            subband_sum /= (float) numfiles;
            /* Remove the channel average */
            for (jj = offset; jj < offset + numfiles; jj++) {
                subbanddata[jj] -= subband_sum;
            }
        }
    }

    currentblock += 1;
}


static int get_data(float **outdata, int blocksperread,
                    struct spectra_info *s,
                    mask * obsmask, int *idispdts, int **offsets, int *padding)
{
    static int firsttime = 1, *maskchans = NULL, blocksize;
    static int worklen, dsworklen;
    static float *tempzz, *data1, *data2, *dsdata1 = NULL, *dsdata2 = NULL;
    static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
    static float *frawdata;
    static double blockdt;
    int totnumread = 0, gotblock = 0, numread = 0;
    int ii, jj, tmppad = 0, nummasked = 0;

    if (firsttime) {
        if (cmd->maskfileP)
            maskchans = gen_ivect(insubs ? s->num_files : s->num_channels);
        worklen = s->spectra_per_subint * blocksperread;
        dsworklen = worklen / cmd->downsamp;
        blocksize = s->spectra_per_subint * cmd->nsub;
        blockdt = s->spectra_per_subint * s->dt;
        if (RAWDATA) {
            frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);
            // To initialize the data reading and prep_subbands routines
            firsttime = 2;
        }
        data1 = gen_fvect(cmd->nsub * worklen);
        data2 = gen_fvect(cmd->nsub * worklen);
        currentdata = data1;
        lastdata = data2;
        if (cmd->downsamp > 1) {
            dsdata1 = gen_fvect(cmd->nsub * dsworklen);
            dsdata2 = gen_fvect(cmd->nsub * dsworklen);
            currentdsdata = dsdata1;
            lastdsdata = dsdata2;
        } else {
            currentdsdata = data1;
            lastdsdata = data2;
        }

        {                       // Make sure that our working blocks are long enough...
            int testlen = insubs ? s->num_files : s->num_channels;
            for (ii = 0; ii < testlen; ii++) {
                if (idispdts[ii] > worklen)
                    printf
                        ("WARNING! (myid = %d):  (idispdts[%d] = %d) > (worklen = %d)\n",
                         myid, ii, idispdts[ii], worklen);
            }
            for (ii = 0; ii < local_numdms; ii++) {
                for (jj = 0; jj < cmd->nsub; jj++) {
                    if (offsets[ii][jj] > dsworklen)
                        printf
                            ("WARNING! (myid = %d):  (offsets[%d][%d] = %d) > (dsworklen = %d)\n",
                             myid, ii, jj, offsets[ii][jj], dsworklen);
                }
            }
        }

    }
    while (1) {
        if (RAWDATA) {
            for (ii = 0; ii < blocksperread; ii++) {
                if (myid == 0) {
                    gotblock = s->get_rawblock(frawdata, s, &tmppad);
                    if (gotblock && !firsttime)
                        totnumread += s->spectra_per_subint;
                }
                MPI_Bcast(&gotblock, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&tmppad, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(frawdata, s->num_channels * s->spectra_per_subint,
                          MPI_FLOAT, 0, MPI_COMM_WORLD);

                if (myid > 0) {
                    if (gotblock) {
                        numread =
                            prep_subbands(currentdata + ii * blocksize,
                                          frawdata, idispdts, cmd->nsub, s,
                                          0, maskchans, &nummasked, obsmask);
                        if (!firsttime)
                            totnumread += numread;
                    } else {
                        *padding = 1;
                        for (jj = ii * blocksize; jj < (ii + 1) * blocksize; jj++)
                            currentdata[jj] = 0.0;
                    }
                    if (tmppad)
                        *padding = 1;
                }
            }
        } else if (insubs) {
            short *subsdata = NULL;

            subsdata = gen_svect(SUBSBLOCKLEN * s->num_files);
            for (ii = 0; ii < blocksperread; ii++) {
                if (myid == 0)
                    numread = simple_read_subbands(s->files, s->num_files, subsdata);
                MPI_Bcast(&numread, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(subsdata, SUBSBLOCKLEN * s->num_files, MPI_SHORT, 0,
                          MPI_COMM_WORLD);
                convert_subbands(s->num_files, subsdata,
                                 currentdata + ii * blocksize, blockdt,
                                 maskchans, &nummasked, obsmask,
                                 cmd->clip, s->padvals);
                if (!firsttime)
                    totnumread += numread;
            }
            vect_free(subsdata);
        }
        /* Downsample the subband data if needed */
        if (myid > 0) {
            if (cmd->downsamp > 1) {
                int kk, offset, dsoffset, index, dsindex;
                float ftmp;
                for (ii = 0; ii < dsworklen; ii++) {
                    dsoffset = ii * cmd->nsub;
                    offset = dsoffset * cmd->downsamp;
                    for (jj = 0; jj < cmd->nsub; jj++) {
                        dsindex = dsoffset + jj;
                        index = offset + jj;
                        currentdsdata[dsindex] = 0.0;
                        for (kk = 0, ftmp = 0.0; kk < cmd->downsamp; kk++) {
                            ftmp += currentdata[index];
                            index += cmd->nsub;
                        }
                        currentdsdata[dsindex] += ftmp / cmd->downsamp;
                    }
                }
            }
        }
        if (firsttime) {
            SWAP(currentdata, lastdata);
            SWAP(currentdsdata, lastdsdata);
            firsttime -= 1;
        } else
            break;
    }
    if (myid > 0) {
        for (ii = 0; ii < local_numdms; ii++)
            float_dedisp(currentdsdata, lastdsdata, dsworklen,
                         cmd->nsub, offsets[ii], 0.0, outdata[ii]);
    }
    SWAP(currentdata, lastdata);
    SWAP(currentdsdata, lastdsdata);
    if (totnumread != worklen) {
        if (cmd->maskfileP)
            vect_free(maskchans);
        vect_free(data1);
        vect_free(data2);
        if (RAWDATA)
            vect_free(frawdata);
        if (cmd->downsamp > 1) {
            vect_free(dsdata1);
            vect_free(dsdata2);
        }
    }
    return totnumread;
}
