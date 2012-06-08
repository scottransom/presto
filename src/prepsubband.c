#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "prepsubband_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "gmrt.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "psrfits.h"

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP || cmd->gmrtP || cmd->spigotP || cmd->filterbankP || cmd->psrfitsP)

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

static void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite);
static void write_subs(FILE * outfiles[], int numfiles, short **subsdata,
                       int startpoint, int numtowrite);
static void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite);
static int read_subbands(FILE * infiles[], int numfiles,
                         float *subbanddata, double timeperblk,
                         int *maskchans, int *nummasked, mask * obsmask,
                         float clip_sigma, float *padvals);
static int get_data(FILE * infiles[], int numfiles, float **outdata,
                    int numchan, int blocklen, int blocksperread,
                    mask * obsmask, float *padvals, double dt,
                    double *dispdts, int **offsets, int *padding, short **subsdata);
static void update_infodata(infodata * idata, int datawrote, int padwrote,
                            int *barybins, int numbarybins, int downsamp);
static void print_percent_complete(int current, int number);

/* From CLIG */
static int insubs = 0;
static Cmdline *cmd;
static IFs ifs = SUMIFS;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   /* Any variable that begins with 't' means topocentric */
   /* Any variable that begins with 'b' means barycentric */
   FILE **infiles, **outfiles;
   float **outdata = NULL, *padvals = NULL;
   short **subsdata = NULL;
   double dtmp, *dms = NULL, avgdm = 0.0, maxdm, dsdt = 0;
   double *dispdt, tlotoa = 0.0, blotoa = 0.0, BW_ddelay = 0.0;
   double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
   double *btoa = NULL, *ttoa = NULL, avgvoverc = 0.0;
   char obs[3], ephem[10], rastring[50], decstring[50];
   int numinfiles, numchan = 1, totnumtowrite, **offsets;
   int ii, jj, numadded = 0, numremoved = 0, padding = 0;
   int numbarypts = 0, blocklen = 0, blocksperread = 0, worklen = 0;
   int numread = 0, numtowrite = 0, totwrote = 0, datawrote = 0;
   int padwrote = 0, padtowrite = 0, statnum = 0, good_padvals = 0;
   int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL;
   char *datafilenm;
   PKMB_tapehdr hdr;
   infodata idata;
   mask obsmask;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(1);
   }

   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);
   numinfiles = cmd->argc;
   if (cmd->noclipP)
      cmd->clip = 0.0;
   /* Which IFs will we use? */
   if (cmd->ifsP) {
      if (cmd->ifs == 0)
         ifs = IF0;
      else if (cmd->ifs == 1)
         ifs = IF1;
      else
         ifs = SUMIFS;
   }
#ifdef DEBUG
   showOptionValues();
#endif

   printf("\n\n");
   printf("          Pulsar Subband De-dispersion Routine\n");
   printf("                 by Scott M. Ransom\n\n");

   if (!RAWDATA) {
      char *root, *suffix;
      /* Split the filename into a rootname and a suffix */
      if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
         printf("\nThe input filename (%s) must have a suffix!\n\n", cmd->argv[0]);
         exit(1);
      } else {
         if (strcmp(suffix, "bcpm1") == 0 || strcmp(suffix, "bcpm2") == 0) {
            printf("Assuming the data is from a GBT BCPM...\n");
            cmd->bcpmP = 1;
         } else if (strcmp(suffix, "fil") == 0 || strcmp(suffix, "fb") == 0) {
            printf("Assuming the data is in SIGPROC filterbank format...\n");
            cmd->filterbankP = 1;
         } else if (strcmp(suffix, "fits") == 0) {
             if (strstr(root, "spigot_5") != NULL) {
                 printf("Assuming the data is from the NRAO/Caltech Spigot card...\n");
                 cmd->spigotP = 1;
             } else if (is_PSRFITS(cmd->argv[0])) {
                 printf("Assuming the data is in PSRFITS format.\n");
                 cmd->psrfitsP = 1;
             } 
         } else if (strcmp(suffix, "pkmb") == 0) {
            printf
                ("Assuming the data is from the Parkes/Jodrell 1-bit filterbank system...\n");
            cmd->pkmbP = 1;
         } else if (strncmp(suffix, "gmrt", 4) == 0) {
            printf("Assuming the data is from the GMRT Phased Array system...\n");
            cmd->gmrtP = 1;
         } else if (isdigit(suffix[0]) && isdigit(suffix[1]) && isdigit(suffix[2])) {
            printf("Assuming the data is from the Arecibo WAPP system...\n");
            cmd->wappP = 1;
         } else if (strcmp(suffix, "sub0") == 0 ||
                    strcmp(suffix, "sub00") == 0 ||
                    strcmp(suffix, "sub000") == 0 ||
                    strcmp(suffix, "sub0000") == 0) {
            printf
                ("Assuming the data is subband data generated by prepsubband...\n");
            insubs = 1;
         } else {
            printf
                ("\nCannot determine the format of the input files '%s'...\n\n",
                 cmd->argv[0]);
            exit(1);
         }
         free(root);
         free(suffix);
      }
   }

   if (cmd->pkmbP) {
      if (numinfiles > 1)
         printf
             ("Reading 1-bit filterbank (Parkes/Jodrell) data from %d files:\n",
              numinfiles);
      else
         printf("Reading 1-bit filterbank (Parkes/Jodrell) data from 1 file:\n");
   } else if (cmd->bcpmP) {
      if (numinfiles > 1)
         printf("Reading Green Bank BCPM data from %d files:\n", numinfiles);
      else
         printf("Reading Green Bank BCPM data from 1 file:\n");
   } else if (cmd->filterbankP) {
      if (numinfiles > 1)
         printf("Reading SIGPROC filterbank data from %d files:\n", numinfiles);
      else
         printf("Reading SIGPROC filterbank data from 1 file:\n");
   } else if (cmd->psrfitsP) {
      if (numinfiles > 1)
         printf("Reading PSRFITS search-mode data from %d files:\n", numinfiles);
      else
         printf("Reading PSRFITS search-mode data from 1 file:\n");
   } else if (cmd->spigotP) {
      if (numinfiles > 1)
         printf("Reading Green Bank Spigot data from %d files:\n", numinfiles);
      else
         printf("Reading Green Bank Spigot data from 1 file:\n");
   } else if (cmd->gmrtP) {
      if (numinfiles > 1)
         printf("Reading GMRT Phased Array data from %d files:\n", numinfiles);
      else
         printf("Reading GMRT Phased Array data from 1 file:\n");
   } else if (cmd->wappP) {
      if (numinfiles > 1)
         printf("Reading Arecibo WAPP data from %d files:\n", numinfiles);
      else
         printf("Reading Arecibo WAPP data from 1 file:\n");
   } else if (insubs) {
      if (numinfiles > 1)
         printf("Reading subband data from %d files:\n", numinfiles);
      else
         printf("Reading subband data from 1 file:\n");
   }

   /* Open the raw data files */

   infiles = (FILE **) malloc(numinfiles * sizeof(FILE *));
   for (ii = 0; ii < numinfiles; ii++) {
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
   }
   if (!cmd->numoutP)
      cmd->numout = INT_MAX;

   /* Determine the output file names and open them */

   datafilenm = (char *) calloc(strlen(cmd->outfile) + 20, 1);
   if (!cmd->subP) {
      printf("\nWriting output data to:\n");
      outfiles = (FILE **) malloc(cmd->numdms * sizeof(FILE *));
      dms = gen_dvect(cmd->numdms);
      for (ii = 0; ii < cmd->numdms; ii++) {
         dms[ii] = cmd->lodm + ii * cmd->dmstep;
         avgdm += dms[ii];
         sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
         outfiles[ii] = chkfopen(datafilenm, "wb");
         printf("   '%s'\n", datafilenm);
      }
      avgdm /= cmd->numdms;
      maxdm = dms[cmd->numdms - 1];
   } else {
      char format_str[30];
      int num_places;

      if (!cmd->nobaryP) {
         printf("\nWarning:  You cannot (currently) barycenter subbands.\n"
                "          Setting the '-nobary' flag automatically.\n");
         cmd->nobaryP = 1;
      }
      printf("\nWriting subbands to:\n");
      cmd->numdms = 1;
      dms = gen_dvect(cmd->numdms);
      dms[0] = cmd->subdm;
      cmd->lodm = cmd->subdm;
      avgdm = cmd->subdm;
      maxdm = cmd->subdm;
      outfiles = (FILE **) malloc(cmd->nsub * sizeof(FILE *));
      num_places = (int) ceil(log10(cmd->nsub));
      sprintf(format_str, "%%s_DM%%.2f.sub%%0%dd", num_places);
      for (ii = 0; ii < cmd->nsub; ii++) {
         sprintf(datafilenm, format_str, cmd->outfile, avgdm, ii);
         outfiles[ii] = chkfopen(datafilenm, "wb");
         printf("   '%s'\n", datafilenm);
      }
   }

   /* Read an input mask if wanted */

   if (cmd->maskfileP) {
      read_mask(cmd->maskfile, &obsmask);
      printf("Read mask information from '%s'\n\n", cmd->maskfile);
      good_padvals = determine_padvals(cmd->maskfile, &obsmask, &padvals);
   } else {
      obsmask.numchan = obsmask.numint = 0;
   }

   /* Set-up values if we are reading subband data  */

   if (insubs) {
      char *root, *suffix, *tmpname;
      /* Split the filename into a rootname and a suffix */
      if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
         printf("\nThe input filename (%s) must have a suffix!\n\n", cmd->argv[0]);
         exit(1);
      } else {
         if (strncmp(suffix, "sub", 3) == 0) {
            tmpname = calloc(strlen(root) + 6, 1);
            sprintf(tmpname, "%s.sub", root);
            readinf(&idata, tmpname);
            free(tmpname);
            free(root);
            free(suffix);
         } else {
            printf
                ("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                 cmd->argv[0]);
            exit(1);
         }
      }
      cmd->nsub = cmd->argc;
      if (padvals==NULL) {
         padvals = gen_fvect(cmd->nsub);
         for (ii = 0 ; ii < cmd->nsub ; ii++)
            padvals[ii] = 0.0;
      }
      numchan = idata.num_chan;
      dsdt = cmd->downsamp * idata.dt;
      avgdm = idata.dm;
      blocklen = SUBSBLOCKLEN;
      BW_ddelay = delay_from_dm(maxdm, idata.freq) - 
         delay_from_dm(maxdm, idata.freq + (idata.num_chan-1) * idata.chan_wid);
      blocksperread = ((int) (BW_ddelay / idata.dt) / blocklen + 1);
      worklen = blocklen * blocksperread;

      /* The number of topo to bary time points to generate with TEMPO */
      numbarypts = (int) (idata.N * idata.dt * 1.1 / TDT + 5.5) + 1;
   }

   if (RAWDATA) {
      double dt, T;
      int ptsperblock;
      long long N;

      /* Set-up values if we are using the Parkes multibeam */
      if (cmd->pkmbP) {
         printf("\nFilterbank input file information:\n");
         get_PKMB_file_info(infiles, numinfiles, cmd->clip, &N, 
                            &ptsperblock, &numchan, &dt, &T, 1);
         chkfread(&hdr, 1, HDRLEN, infiles[0]);
         rewind(infiles[0]);
         PKMB_hdr_to_inf(&hdr, &idata);
         PKMB_update_infodata(numinfiles, &idata);
      }

      /* Set-up values if we are using the GMRT Phased Array system */
      if (cmd->gmrtP) {
         printf("\nGMRT input file information:\n");
         get_GMRT_file_info(infiles, argv + 1, numinfiles, cmd->clip,
                            &N, &ptsperblock, &numchan, &dt, &T, 1);
         /* Read the first header file and generate an infofile from it */
         GMRT_hdr_to_inf(argv[1], &idata);
         GMRT_update_infodata(numinfiles, &idata);
         set_GMRT_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using SIGPROC filterbank-style data */
      if (cmd->filterbankP) {
         int headerlen;
         sigprocfb fb;

         /* Read the first header file and generate an infofile from it */
         rewind(infiles[0]);
         headerlen = read_filterbank_header(&fb, infiles[0]);
         sigprocfb_to_inf(&fb, &idata);
         rewind(infiles[0]);
         printf("\nSIGPROC filterbank input file information:\n");
         get_filterbank_file_info(infiles, numinfiles, cmd->clip,
                                  &N, &ptsperblock, &numchan, &dt, &T, 1);
         filterbank_update_infodata(numinfiles, &idata);
         set_filterbank_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using the Berkeley-Caltech */
      /* Pulsar Machine (or BPP) format.                    */
      if (cmd->bcpmP) {
         printf("\nBCPM input file information:\n");
         get_BPP_file_info(infiles, numinfiles, cmd->clip, &N,
                           &ptsperblock, &numchan, &dt, &T, &idata, 1);
         BPP_update_infodata(numinfiles, &idata);
         set_BPP_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using the NRAO-Caltech Spigot card */
      if (cmd->spigotP) {
         SPIGOT_INFO *spigots;

         printf("\nSpigot card input file information:\n");
         spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numinfiles);
         for (ii = 0; ii < numinfiles; ii++)
            read_SPIGOT_header(cmd->argv[ii], spigots + ii);
         get_SPIGOT_file_info(infiles, spigots, numinfiles, cmd->windowP,
                              cmd->clip, &N, &ptsperblock, &numchan, &dt,
                              &T, &idata, 1);
         SPIGOT_update_infodata(numinfiles, &idata);
         set_SPIGOT_padvals(padvals, good_padvals);
         free(spigots);
      }

      /* Set-up values if we are using search-mode PSRFITS data */
      if (cmd->psrfitsP) {
         struct spectra_info s;
         
         printf("PSRFITS input file information:\n");
          // -1 causes the data to determine if we use weights, scales, & offsets
         s.apply_weight = (cmd->noweightsP) ? 0 : -1;
         s.apply_scale  = (cmd->noscalesP) ? 0 : -1;
         s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
         read_PSRFITS_files(cmd->argv, cmd->argc, &s);
         N = s.N;
         ptsperblock = s.spectra_per_subint;
         numchan = s.num_channels;
         dt = s.dt;
         T = s.T;
         get_PSRFITS_file_info(cmd->argv, cmd->argc, cmd->clip, 
                               &s, &idata, 1);
         PSRFITS_update_infodata(&idata);
         set_PSRFITS_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using the Arecobo WAPP */
      if (cmd->wappP) {
         printf("\nWAPP input file information:\n");
         get_WAPP_file_info(infiles, cmd->numwapps, numinfiles,
                            cmd->windowP, cmd->clip, &N, &ptsperblock,
                            &numchan, &dt, &T, &idata, 1);
         WAPP_update_infodata(numinfiles, &idata);
         set_WAPP_padvals(padvals, good_padvals);
      }

      /* Finish setting up stuff common to all raw formats */
      dsdt = cmd->downsamp * idata.dt;
      idata.dm = avgdm;
      blocklen = ptsperblock;

      /* Hack to shift to a later part in the raw data */
      {
         char *envval = getenv("BCPM_OBS_OFFSET");

         if (envval != NULL) {
            double startT = strtod(envval, NULL);
            if (startT) {
               double recdt, newstartT, startTday;
               long lorec, numrec;

               numrec = idata.N / blocklen;
               recdt = blocklen * idata.dt;
               lorec = (long) (startT * numrec + DBLCORRECT);
               newstartT = (double) (lorec * blocklen) / idata.N;
               startTday = lorec * recdt / SECPERDAY;
               idata.mjd_f += startTday;
               if (idata.mjd_f > 1.0) {
                  idata.mjd_f -= 1.0;
                  idata.mjd_i += 1;
               }
               printf("\n  Found startT = %.10g, newstartT = %.10g:\n",
                      startT, newstartT);
               printf
                   ("    numrec = %ld, lorec = %ld, timeoffset = %.10g:\n",
                    numrec, lorec, startTday);
               skip_to_BPP_rec(infiles, numinfiles, lorec + 1);
            }
         }
      }

      BW_ddelay = delay_from_dm(maxdm, idata.freq) - 
         delay_from_dm(maxdm, idata.freq + (idata.num_chan-1) * idata.chan_wid);
      blocksperread = ((int) (BW_ddelay / idata.dt) / blocklen + 1);
      worklen = blocklen * blocksperread;
      /* The number of topo to bary time points to generate with TEMPO */
      numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
   }

   // Identify the TEMPO observatory code
   {
       char *outscope = (char *) calloc(40, sizeof(char));
       telescope_to_tempocode(idata.telescope, outscope, obs);
       free(outscope);
   }

   if (cmd->nsub > numchan) {
      printf
          ("Warning:  The number of requested subbands (%d) is larger than the number of channels (%d).\n",
           cmd->nsub, numchan);
      printf("          Re-setting the number of subbands to %d.\n\n", numchan);
      cmd->nsub = numchan;
   }

   if (blocklen % cmd->downsamp) {
      printf
          ("Error:  The downsample factor (%d) must be a factor of the\n",
           cmd->downsamp);
      printf("        blocklength (%d).  Exiting.\n\n", blocklen);
      exit(1);
   }

   tlotoa = idata.mjd_i + idata.mjd_f;  /* Topocentric epoch */

   if (cmd->numoutP)
      totnumtowrite = cmd->numout;
   else
      totnumtowrite = (int) idata.N / cmd->downsamp;

   if (cmd->nobaryP) {          /* Main loop if we are not barycentering... */

      /* Dispersion delays (in bins).  The high freq gets no delay   */
      /* All other delays are positive fractions of bin length (dt)  */

      dispdt = subband_search_delays(numchan, cmd->nsub, avgdm,
                                     idata.freq, idata.chan_wid, 0.0);
      for (ii = 0; ii < numchan; ii++)
         dispdt[ii] /= idata.dt;

      /* The subband dispersion delays (see note above) */

      offsets = gen_imatrix(cmd->numdms, cmd->nsub);
      for (ii = 0; ii < cmd->numdms; ii++) {
         double *subdispdt;

         subdispdt = subband_delays(numchan, cmd->nsub, dms[ii],
                                    idata.freq, idata.chan_wid, 0.0);
         dtmp = subdispdt[cmd->nsub - 1];
         for (jj = 0; jj < cmd->nsub; jj++)
            offsets[ii][jj] = NEAREST_INT((subdispdt[jj] - dtmp) / dsdt);
         vect_free(subdispdt);
      }

      /* Allocate our data array and start getting data */

      printf("De-dispersing using:\n");
      printf("       Subbands = %d\n", cmd->nsub);
      printf("     Average DM = %.7g\n", avgdm);
      if (cmd->downsamp > 1) {
         printf("     Downsample = %d\n", cmd->downsamp);
         printf("  New sample dt = %.10g\n", dsdt);
      }
      printf("\n");

      if (cmd->subP)
         subsdata = gen_smatrix(cmd->nsub, worklen / cmd->downsamp);
      else
         outdata = gen_fmatrix(cmd->numdms, worklen / cmd->downsamp);
      numread = get_data(infiles, numinfiles, outdata,
                         numchan, blocklen, blocksperread,
                         &obsmask, padvals, idata.dt, dispdt,
                         offsets, &padding, subsdata);

      while (numread == worklen) {

         numread /= cmd->downsamp;
         print_percent_complete(totwrote, totnumtowrite);

         /* Write the latest chunk of data, but don't   */
         /* write more than cmd->numout points.         */

         numtowrite = numread;
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         if (cmd->subP)
            write_subs(outfiles, cmd->nsub, subsdata, 0, numtowrite);
         else
            write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
         totwrote += numtowrite;

         /* Update the statistics */

         if (!padding && !cmd->subP) {
            for (ii = 0; ii < numtowrite; ii++)
               update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg, &var);
            statnum += numtowrite;
         }

         /* Stop if we have written out all the data we need to */

         if (cmd->numoutP && (totwrote == cmd->numout))
            break;

         numread = get_data(infiles, numinfiles, outdata,
                            numchan, blocklen, blocksperread,
                            &obsmask, padvals, idata.dt, dispdt,
                            offsets, &padding, subsdata);
      }
      datawrote = totwrote;

   } else {                     /* Main loop if we are barycentering... */
      double maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;

      /* What ephemeris will we use?  (Default is DE200) */

      if (cmd->de405P)
         strcpy(ephem, "DE405");
      else
         strcpy(ephem, "DE200");

      /* Define the RA and DEC of the observation */

      ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
      ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

      /* Allocate some arrays */

      btoa = gen_dvect(numbarypts);
      ttoa = gen_dvect(numbarypts);
      voverc = gen_dvect(numbarypts);
      for (ii = 0; ii < numbarypts; ii++)
         ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

      /* Call TEMPO for the barycentering */

      printf("\nGenerating barycentric corrections...\n");
      barycenter(ttoa, btoa, voverc, numbarypts, rastring, decstring, obs, ephem);
      for (ii = 0; ii < numbarypts; ii++) {
         if (voverc[ii] > maxvoverc)
            maxvoverc = voverc[ii];
         if (voverc[ii] < minvoverc)
            minvoverc = voverc[ii];
         avgvoverc += voverc[ii];
      }
      avgvoverc /= numbarypts;
      vect_free(voverc);
      blotoa = btoa[0];

      printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
      printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
      printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
      printf("De-dispersing and barycentering using:\n");
      printf("       Subbands = %d\n", cmd->nsub);
      printf("     Average DM = %.7g\n", avgdm);
      if (cmd->downsamp > 1) {
         printf("     Downsample = %d\n", cmd->downsamp);
         printf("  New sample dt = %.10g\n", dsdt);
      }
      printf("\n");

      /* Dispersion delays (in bins).  The high freq gets no delay   */
      /* All other delays are positive fractions of bin length (dt)  */

      dispdt = subband_search_delays(numchan, cmd->nsub, avgdm,
                                     idata.freq, idata.chan_wid, avgvoverc);
      for (ii = 0; ii < numchan; ii++)
         dispdt[ii] /= idata.dt;

      /* The subband dispersion delays (see note above) */

      offsets = gen_imatrix(cmd->numdms, cmd->nsub);
      for (ii = 0; ii < cmd->numdms; ii++) {
         double *subdispdt;

         subdispdt = subband_delays(numchan, cmd->nsub, dms[ii],
                                    idata.freq, idata.chan_wid, avgvoverc);
         dtmp = subdispdt[cmd->nsub - 1];
         for (jj = 0; jj < cmd->nsub; jj++)
            offsets[ii][jj] = NEAREST_INT((subdispdt[jj] - dtmp) / dsdt);
         vect_free(subdispdt);
      }

      /* Convert the bary TOAs to differences from the topo TOAs in */
      /* units of bin length (dt) rounded to the nearest integer.   */

      dtmp = (btoa[0] - ttoa[0]);
      for (ii = 0; ii < numbarypts; ii++)
         btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

      {                         /* Find the points where we need to add or remove bins */

         int oldbin = 0, currentbin;
         double lobin, hibin, calcpt;

         numdiffbins = abs(NEAREST_INT(btoa[numbarypts - 1])) + 1;
         diffbins = gen_ivect(numdiffbins);
         diffbinptr = diffbins;
         for (ii = 1; ii < numbarypts; ii++) {
            currentbin = NEAREST_INT(btoa[ii]);
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
                  *diffbinptr = NEAREST_INT(LININTERP(calcpt, btoa[ii - 1],
                                                      btoa[ii], lobin, hibin));
                  diffbinptr++;
                  calcpt = (currentbin > 0) ? calcpt + 1.0 : calcpt - 1.0;
               }
               oldbin = currentbin;
            }
         }
         *diffbinptr = cmd->numout; /* Used as a marker */
      }
      diffbinptr = diffbins;

      /* Now perform the barycentering */

      if (cmd->subP)
         subsdata = gen_smatrix(cmd->nsub, worklen / cmd->downsamp);
      else
         outdata = gen_fmatrix(cmd->numdms, worklen / cmd->downsamp);
      numread = get_data(infiles, numinfiles, outdata,
                         numchan, blocklen, blocksperread,
                         &obsmask, padvals, idata.dt, dispdt,
                         offsets, &padding, subsdata);

      while (numread == worklen) {      /* Loop to read and write the data */
         int numwritten = 0;
         double block_avg, block_var;

         numread /= cmd->downsamp;
         /* Determine the approximate local average */
         avg_var(outdata[0], numread, &block_avg, &block_var);
         print_percent_complete(totwrote, totnumtowrite);

         /* Simply write the data if we don't have to add or */
         /* remove any bins from this batch.                 */
         /* OR write the amount of data up to cmd->numout or */
         /* the next bin that will be added or removed.      */

         numtowrite = abs(*diffbinptr) - datawrote;
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         if (numtowrite > numread)
            numtowrite = numread;
         if (cmd->subP)
            write_subs(outfiles, cmd->nsub, subsdata, 0, numtowrite);
         else
            write_data(outfiles, cmd->numdms, outdata, 0, numtowrite);
         datawrote += numtowrite;
         totwrote += numtowrite;
         numwritten += numtowrite;

         /* Update the statistics */

         if (!padding && !cmd->subP) {
            for (ii = 0; ii < numtowrite; ii++)
               update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg, &var);
            statnum += numtowrite;
         }

         if ((datawrote == abs(*diffbinptr)) && 
             (numwritten != numread) && 
             (totwrote < cmd->numout)) {  /* Add/remove a bin */
            int skip, nextdiffbin;

            skip = numtowrite;

            do {  /* Write the rest of the data after adding/removing a bin  */

               if (*diffbinptr > 0) {
                  /* Add a bin */
                  write_padding(outfiles, cmd->numdms, block_avg, 1);
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
               if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
                  numtowrite = cmd->numout - totwrote;
               nextdiffbin = abs(*diffbinptr) - datawrote;
               if (numtowrite > nextdiffbin)
                  numtowrite = nextdiffbin;
               if (cmd->subP)
                  write_subs(outfiles, cmd->nsub, subsdata, skip, numtowrite);
               else
                  write_data(outfiles, cmd->numdms, outdata, skip, numtowrite);
               numwritten += numtowrite;
               datawrote += numtowrite;
               totwrote += numtowrite;

               /* Update the statistics and counters */

               if (!padding && !cmd->subP) {
                  for (ii = 0; ii < numtowrite; ii++)
                     update_stats(statnum + ii, outdata[0][skip + ii],
                                  &min, &max, &avg, &var);
                  statnum += numtowrite;
               }
               skip += numtowrite;

               /* Stop if we have written out all the data we need to */

               if (cmd->numoutP && (totwrote == cmd->numout))
                  break;
            } while (numwritten < numread);
         }
         /* Stop if we have written out all the data we need to */

         if (cmd->numoutP && (totwrote == cmd->numout))
            break;

         numread = get_data(infiles, numinfiles, outdata,
                            numchan, blocklen, blocksperread,
                            &obsmask, padvals, idata.dt, dispdt,
                            offsets, &padding, subsdata);
      }
   }

   /* Calculate the amount of padding we need  */

   if (cmd->numoutP && (cmd->numout > totwrote))
      padwrote = padtowrite = cmd->numout - totwrote;

   /* Write the new info file for the output data */

   idata.dt = dsdt;
   update_infodata(&idata, totwrote, padtowrite, diffbins,
                   numdiffbins, cmd->downsamp);
   for (ii = 0; ii < cmd->numdms; ii++) {
      idata.dm = dms[ii];
      if (!cmd->nobaryP) {
         double baryepoch, barydispdt, baryhifreq;

         baryhifreq = idata.freq + (numchan - 1) * idata.chan_wid;
         barydispdt = delay_from_dm(dms[ii], doppler(baryhifreq, avgvoverc));
         baryepoch = blotoa - (barydispdt / SECPERDAY);
         idata.bary = 1;
         idata.mjd_i = (int) floor(baryepoch);
         idata.mjd_f = baryepoch - idata.mjd_i;
      }
      if (cmd->subP)
         sprintf(idata.name, "%s_DM%.2f.sub", cmd->outfile, dms[ii]);
      else
         sprintf(idata.name, "%s_DM%.2f", cmd->outfile, dms[ii]);
      writeinf(&idata);
   }

   /* Set the padded points equal to the average data point */

   if (idata.numonoff >= 1) {
      int index, startpad, endpad;

      for (ii = 0; ii < cmd->numdms; ii++) {
         fclose(outfiles[ii]);
         sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
         outfiles[ii] = chkfopen(datafilenm, "rb+");
      }
      for (ii = 0; ii < idata.numonoff; ii++) {
         index = 2 * ii;
         startpad = idata.onoff[index + 1];
         if (ii == idata.numonoff - 1)
            endpad = idata.N - 1;
         else
            endpad = idata.onoff[index + 2];
         for (jj = 0; jj < cmd->numdms; jj++)
            chkfseek(outfiles[jj], (startpad + 1) * sizeof(float), SEEK_SET);
         padtowrite = endpad - startpad;
         write_padding(outfiles, cmd->numdms, avg, padtowrite);
      }
   }

   /* Print simple stats and results */

   if (!cmd->subP) {
      var /= (datawrote - 1);
      print_percent_complete(1, 1);
      printf("\n\nDone.\n\nSimple statistics of the output data:\n");
      printf("             Data points written:  %d\n", totwrote);
      if (padwrote)
         printf("          Padding points written:  %d\n", padwrote);
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
   } else {
      printf("\n\nDone.\n");
      printf("             Data points written:  %d\n", totwrote);
      if (padwrote)
         printf("          Padding points written:  %d\n", padwrote);
      if (!cmd->nobaryP) {
         if (numadded)
            printf("    Bins added for barycentering:  %d\n", numadded);
         if (numremoved)
            printf("  Bins removed for barycentering:  %d\n", numremoved);
      }
      printf("\n");
   }

   /* Close the files and cleanup */

   if (cmd->maskfileP) {
      free_mask(obsmask);
      vect_free(padvals);
   }
   for (ii = 0; ii < numinfiles; ii++)
      fclose(infiles[ii]);
   free(infiles);
   for (ii = 0; ii < cmd->numdms; ii++)
      fclose(outfiles[ii]);
   if (cmd->subP) {
      vect_free(subsdata[0]);
      vect_free(subsdata);
   } else {
      vect_free(outdata[0]);
      vect_free(outdata);
   }
   free(outfiles);
   vect_free(dms);
   vect_free(dispdt);
   vect_free(offsets[0]);
   vect_free(offsets);
   free(datafilenm);
   if (!cmd->nobaryP) {
      vect_free(btoa);
      vect_free(ttoa);
      vect_free(diffbins);
   }
   return (0);
}

static void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite)
{
   int ii;

   for (ii = 0; ii < numfiles; ii++)
      chkfwrite(outdata[ii] + startpoint, sizeof(float), numtowrite, outfiles[ii]);
}


static void write_subs(FILE * outfiles[], int numfiles, short **subsdata,
                       int startpoint, int numtowrite)
{
   int ii;

   for (ii = 0; ii < numfiles; ii++)
      chkfwrite(subsdata[ii] + startpoint, sizeof(short), numtowrite, outfiles[ii]);
}


static void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite)
{
   int ii;

   if (numtowrite <= 0) {
      return;
   } else if (numtowrite == 1) {
      for (ii = 0; ii < numfiles; ii++)
         chkfwrite(&value, sizeof(float), 1, outfiles[ii]);
   } else {
      int maxatonce = 8192, veclen, jj;
      float *buffer;
      veclen = (numtowrite > maxatonce) ? maxatonce : numtowrite;
      buffer = gen_fvect(veclen);
      for (ii = 0; ii < veclen; ii++)
         buffer[ii] = value;
      if (veclen == numtowrite) {
         for (ii = 0; ii < numfiles; ii++)
            chkfwrite(buffer, sizeof(float), veclen, outfiles[ii]);
      } else {
         for (ii = 0; ii < numtowrite / veclen; ii++) {
            for (jj = 0; jj < numfiles; jj++)
               chkfwrite(buffer, sizeof(float), veclen, outfiles[jj]);
         }
         for (jj = 0; jj < numfiles; jj++)
            chkfwrite(buffer, sizeof(float), numtowrite % veclen, outfiles[jj]);
      }
      vect_free(buffer);
   }
}


static int read_subbands(FILE * infiles[], int numfiles,
                         float *subbanddata, double timeperblk,
                         int *maskchans, int *nummasked, mask * obsmask,
                         float clip_sigma, float *padvals)
/* Read short int subband data written by prepsubband */
{
   int ii, jj, index, numread = 0, mask = 0,offset;
   short subsdata[SUBSBLOCKLEN]; 
   double starttime, run_avg;
   float subband_sum;
   static int currentblock = 0;

   if (obsmask->numchan) mask = 1;
   
   /* Read the data */
   for (ii = 0; ii < numfiles; ii++) {
      numread = chkfread(subsdata, sizeof(short), SUBSBLOCKLEN, infiles[ii]);
      run_avg = 0.0;
      if (cmd->runavgP==1) {
          for (jj = 0; jj < numread ; jj++)
              run_avg += (float) subsdata[jj];
          run_avg /= numread;
      }
      for (jj = 0, index = ii; jj < numread; jj++, index += numfiles)
         subbanddata[index] = (float) subsdata[jj] - run_avg;
      for (jj = numread; jj < SUBSBLOCKLEN; jj++, index += numfiles)
         subbanddata[index] = 0.0;
      index += numread;
   }

   if (mask) {
      starttime = currentblock * timeperblk;
      *nummasked = check_mask(starttime, timeperblk, obsmask, maskchans);
   }

   /* Clip nasty RFI if requested and we're not masking all the channels*/
   if ((clip_sigma > 0.0) && !(mask && (*nummasked == -1))){
      subs_clip_times(subbanddata, SUBSBLOCKLEN, numfiles, clip_sigma, padvals);
   }

   /* Mask it if required */
   if (mask && numread) {
      if (*nummasked == -1) {   /* If all channels are masked */
         for (ii = 0; ii < SUBSBLOCKLEN; ii++)
            memcpy(subbanddata + ii * numfiles, padvals, sizeof(float) * numfiles);
      } else if (*nummasked > 0) {      /* Only some of the channels are masked */
         int offset, channum;
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
   if (cmd->zerodmP==1) {
       for (ii = 0; ii < SUBSBLOCKLEN; ii++) {
           offset = ii * numfiles;
           subband_sum = 0.0; 
           for (jj = offset; jj < offset+numfiles; jj++) {
               subband_sum += subbanddata[jj];
           }    
           subband_sum /= (float) numfiles;
           /* Remove the channel average */
           for (jj = offset; jj < offset+numfiles; jj++) {
               subbanddata[jj] -= subband_sum;
           }    
       }    
   }

   currentblock += 1;
   return numread;
}


static int get_data(FILE * infiles[], int numfiles, float **outdata,
                    int numchan, int blocklen, int blocksperread,
                    mask * obsmask, float *padvals, double dt,
                    double *dispdts, int **offsets, int *padding, short **subsdata)
{
   static int firsttime = 1, *maskchans = NULL, blocksize;
   static int worklen, dsworklen;
   static float *tempzz, *data1, *data2, *dsdata1 = NULL, *dsdata2 = NULL;
   static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
   static double blockdt;
   int totnumread = 0, numread = 0, ii, jj, tmppad = 0, nummasked = 0;

   if (firsttime) {
      if (cmd->maskfileP)
         maskchans = gen_ivect(numchan);
      worklen = blocklen * blocksperread;
      dsworklen = worklen / cmd->downsamp;
      { // Make sure that out working blocks are long enough...
         for (ii = 0; ii < numchan; ii++) {
            if (dispdts[ii] > worklen)
               printf("WARNING!:  (dispdts[%d] = %.0f) > (worklen = %d)\n", 
                      ii, dispdts[ii], worklen);
         }
         for (ii = 0; ii < cmd->numdms; ii++) {
            for (jj = 0; jj < cmd->nsub; jj++) {
               if (offsets[ii][jj] > dsworklen)
                  printf("WARNING!:  (offsets[%d][%d] = %d) > (dsworklen = %d)\n", 
                         ii, jj, offsets[ii][jj], dsworklen);
            }
         }
      }

      blocksize = blocklen * cmd->nsub;
      blockdt = blocklen * dt;
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
   }
   while (1) {
      if (RAWDATA || insubs) {
         for (ii = 0; ii < blocksperread; ii++) {
            if (cmd->pkmbP)
               numread = read_PKMB_subbands(infiles, numfiles,
                                            currentdata + ii * blocksize,
                                            dispdts, cmd->nsub, 0, &tmppad,
                                            maskchans, &nummasked, obsmask);
            else if (cmd->bcpmP)
               numread = read_BPP_subbands(infiles, numfiles,
                                           currentdata + ii * blocksize,
                                           dispdts, cmd->nsub, 0, &tmppad,
                                           maskchans, &nummasked, obsmask, ifs);
            else if (cmd->spigotP)
               numread = read_SPIGOT_subbands(infiles, numfiles,
                                              currentdata + ii * blocksize,
                                              dispdts, cmd->nsub, 0,
                                              &tmppad, maskchans,
                                              &nummasked, obsmask, ifs);
            else if (cmd->psrfitsP)
               numread = read_PSRFITS_subbands(currentdata + ii * blocksize,
                                               dispdts, cmd->nsub, 0,
                                               &tmppad, maskchans,
                                               &nummasked, obsmask);
            else if (cmd->wappP)
               numread = read_WAPP_subbands(infiles, numfiles,
                                            currentdata + ii * blocksize,
                                            dispdts, cmd->nsub, 0, &tmppad,
                                            maskchans, &nummasked, obsmask, ifs);
            else if (cmd->gmrtP)
               numread = read_GMRT_subbands(infiles, numfiles,
                                            currentdata + ii * blocksize,
                                            dispdts, cmd->nsub, 0, &tmppad,
                                            maskchans, &nummasked, obsmask);
            else if (cmd->filterbankP)
               numread = read_filterbank_subbands(infiles, numfiles,
                                                  currentdata +
                                                  ii * blocksize, dispdts,
                                                  cmd->nsub, 0, &tmppad,
                                                  maskchans, &nummasked, obsmask);
            else if (insubs)
               numread = read_subbands(infiles, numfiles,
                                       currentdata + ii * blocksize,
                                       blockdt, maskchans, &nummasked,
                                       obsmask, cmd->clip, padvals);
            if (!firsttime)
               totnumread += numread;
            if (numread != blocklen) {
               for (jj = ii * blocksize; jj < (ii + 1) * blocksize; jj++)
                  currentdata[jj] = 0.0;
            }
            if (tmppad)
               *padding = 1;
         }
      }
      /* Downsample the subband data if needed */
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
               /* Keep the short ints from overflowing */
               currentdsdata[dsindex] += ftmp / cmd->downsamp;
            }
         }
      }
      if (firsttime) {
         SWAP(currentdata, lastdata);
         SWAP(currentdsdata, lastdsdata);
         firsttime = 0;
      } else
         break;
   }
   if (!cmd->subP) {
      for (ii = 0; ii < cmd->numdms; ii++)
         float_dedisp(currentdsdata, lastdsdata, dsworklen,
                      cmd->nsub, offsets[ii], 0.0, outdata[ii]);
   } else {
      /* Input format is sub1[0], sub2[0], sub3[0], ..., sub1[1], sub2[1], sub3[1], ... */
      float infloat;
      for (ii = 0; ii < cmd->nsub; ii++) {
         for (jj = 0; jj < dsworklen; jj++) {
            infloat = lastdsdata[ii + (cmd->nsub * jj)];
            subsdata[ii][jj] = (short) (infloat + 0.5);
            //if ((float) subsdata[ii][jj] != infloat)
            //   printf
            //       ("Warning:  We are incorrectly converting subband data! float = %f  short = %d\n",
            //         infloat, subsdata[ii][jj]);
         }
      }
   }
   SWAP(currentdata, lastdata);
   SWAP(currentdsdata, lastdsdata);
   if (totnumread != worklen) {
      if (cmd->maskfileP)
         vect_free(maskchans);
      vect_free(data1);
      vect_free(data2);
      if (cmd->downsamp > 1) {
         vect_free(dsdata1);
         vect_free(dsdata2);
      }
   }
   return totnumread;
}


static void print_percent_complete(int current, int number)
{
   static int newper = 0, oldper = -1;

   newper = (int) (current / (float) (number) * 100.0);
   if (newper < 0)
      newper = 0;
   if (newper > 100)
      newper = 100;
   if (newper > oldper) {
      printf("\rAmount complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
   }
}


static void update_infodata(infodata * idata, int datawrote, int padwrote,
                            int *barybins, int numbarybins, int downsamp)
/* Update our infodata for barycentering and padding */
{
   int ii, jj, index;

   idata->N = datawrote + padwrote;
   if (idata->numonoff == 0) {
      if (padwrote) {
         idata->numonoff = 2;
         idata->onoff[0] = 0.0;
         idata->onoff[1] = datawrote - 1;
         idata->onoff[2] = idata->N - 1;
         idata->onoff[3] = idata->N - 1;
      }
      return;
   } else {
      for (ii = 0; ii < idata->numonoff; ii++) {
         idata->onoff[ii * 2] /= downsamp;
         idata->onoff[ii * 2 + 1] /= downsamp;
      }
   }

   /* Determine the barycentric onoff bins (approximate) */

   if (numbarybins) {
      int numadded = 0, numremoved = 0;

      ii = 1;                   /* onoff index    */
      jj = 0;                   /* barybins index */
      while (ii < idata->numonoff * 2) {
         while (abs(barybins[jj]) <= idata->onoff[ii] && jj < numbarybins) {
            if (barybins[jj] < 0)
               numremoved++;
            else
               numadded++;
            jj++;
         }
         idata->onoff[ii] += numadded - numremoved;
         ii++;
      }
   }

   /* Now cut off the extra onoff bins */

   for (ii = 1, index = 1; ii <= idata->numonoff; ii++, index += 2) {
      if (idata->onoff[index - 1] > idata->N - 1) {
         idata->onoff[index - 1] = idata->N - 1;
         idata->onoff[index] = idata->N - 1;
         break;
      }
      if (idata->onoff[index] > datawrote - 1) {
         idata->onoff[index] = datawrote - 1;
         idata->numonoff = ii;
         if (padwrote) {
            idata->numonoff++;
            idata->onoff[index + 1] = idata->N - 1;
            idata->onoff[index + 2] = idata->N - 1;
         }
         break;
      }
   }
}
