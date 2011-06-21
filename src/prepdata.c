#include <limits.h>
#include <ctype.h>
#include "prepdata_cmd.h"
#include "presto.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "gmrt.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "psrfits.h"

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP || cmd->gmrtP || cmd->spigotP || cmd->filterbankP || cmd->psrfitsP)

/* Some function definitions */

static int read_floats(FILE * file, float *data, int numpts, int numchan);
static int read_shorts(FILE * file, float *data, int numpts, int numchan);
static int downsample(float outdata[], int numread, int downsampfact);
static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   /* Any variable that begins with 't' means topocentric */
   /* Any variable that begins with 'b' means barycentric */
   FILE **infiles, *outfile;
   float *outdata = NULL, *padvals;
   double tdf = 0.0, dtmp = 0.0, barydispdt = 0.0, dsdt = 0.0;
   double *dispdt, *tobsf = NULL, tlotoa = 0.0, blotoa = 0.0;
   double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
   char obs[3], ephem[10], *datafilenm, *outinfonm, *root, *suffix;
   char rastring[50], decstring[50];
   int numfiles, numchan = 1, newper = 0, oldper = 0, nummasked = 0, useshorts = 0;
   int slen, numadded = 0, numremoved = 0, padding = 0, *maskchans = NULL, offset =
       0;
   long ii, numbarypts = 0, worklen = 65536;
   long numread = 0, numtowrite = 0, totwrote = 0, datawrote = 0;
   long padwrote = 0, padtowrite = 0, statnum = 0;
   int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL, good_padvals = 0;
   IFs ifs = SUMIFS;
   infodata idata;
   Cmdline *cmd;
   mask obsmask;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(1);
   }

   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);
   numfiles = cmd->argc;
   infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
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
   printf("           Pulsar Data Preparation Routine\n");
   printf("    Type conversion, de-dispersion, barycentering.\n");
   printf("                 by Scott M. Ransom\n\n");

   if (!RAWDATA) {
      /* Split the filename into a rootname and a suffix */
      if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
         printf("\nThe input filename (%s) must have a suffix!\n\n", cmd->argv[0]);
         exit(1);
      } else {
         if (strcmp(suffix, "sdat") == 0) {
            useshorts = 1;
         } else if (strcmp(suffix, "bcpm1") == 0 || strcmp(suffix, "bcpm2") == 0) {
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
         }
      }
      if (RAWDATA) {            /* Clean-up a bit */
         free(root);
         free(suffix);
      }
   }

   if (!RAWDATA) {
      printf("Reading input data from '%s'.\n", cmd->argv[0]);
      printf("Reading information from '%s.inf'.\n\n", root);
      /* Read the info file if available */
      readinf(&idata, root);
      free(root);
      free(suffix);
      infiles[0] = chkfopen(cmd->argv[0], "rb");
   } else if (cmd->pkmbP) {
      if (numfiles > 1)
         printf("Reading 1-bit filterbank (Parkes/Jodrell) data from %d files:\n",
                numfiles);
      else
         printf("Reading 1-bit filterbank (Parkes/Jodrell) data from 1 file:\n");
   } else if (cmd->bcpmP) {
      if (numfiles > 1)
         printf("Reading Green Bank BCPM data from %d files:\n", numfiles);
      else
         printf("Reading Green Bank BCPM data from 1 file:\n");
   } else if (cmd->filterbankP) {
      if (numfiles > 1)
         printf("Reading SIGPROC filterbank data from %d files:\n", numfiles);
      else
         printf("Reading SIGPROC filterbank data from 1 file:\n");
   } else if (cmd->psrfitsP) {
      if (numfiles > 1)
         printf("Reading PSRFITS search-mode data from %d files:\n", numfiles);
      else
         printf("Reading PSRFITS search-mode data from 1 file:\n");
   } else if (cmd->spigotP) {
      if (numfiles > 1)
         printf("Reading Green Bank Spigot data from %d files:\n", numfiles);
      else
         printf("Reading Green Bank Spigot data from 1 file:\n");
   } else if (cmd->gmrtP) {
      if (numfiles > 1)
         printf("Reading GMRT Phased Array data from %d files:\n", numfiles);
      else
         printf("Reading GMRT Phased Array data from 1 file:\n");
   } else if (cmd->wappP) {
      if (numfiles > 1)
         printf("Reading Arecibo WAPP data from %d files:\n", numfiles);
      else
         printf("Reading Arecibo WAPP data from 1 file:\n");
   }

   /* Open the raw data files */

   if (RAWDATA) {
      for (ii = 0; ii < numfiles; ii++) {
         printf("  '%s'\n", cmd->argv[ii]);
         infiles[ii] = chkfopen(cmd->argv[ii], "rb");
      }
      printf("\n");
   }

   /* Determine the other file names and open the output data file */
   slen = strlen(cmd->outfile) + 8;
   datafilenm = (char *) calloc(slen, 1);
   sprintf(datafilenm, "%s.dat", cmd->outfile);
   outfile = chkfopen(datafilenm, "wb");
   sprintf(idata.name, "%s", cmd->outfile);
   outinfonm = (char *) calloc(slen, 1);
   sprintf(outinfonm, "%s.inf", cmd->outfile);

   /* Read an input mask if wanted */

   if (cmd->maskfileP) {
      read_mask(cmd->maskfile, &obsmask);
      printf("Read mask information from '%s'\n\n", cmd->maskfile);
      good_padvals = determine_padvals(cmd->maskfile, &obsmask, &padvals);
   } else {
      obsmask.numchan = obsmask.numint = 0;
   }

   if (RAWDATA) {
      double dt, T;
      int ptsperblock;
      long long N;

      /* Set-up values if we are using the Parkes multibeam */
      if (cmd->pkmbP) {
         PKMB_tapehdr hdr;

         printf("Filterbank input file information:\n");
         get_PKMB_file_info(infiles, numfiles, cmd->clip, &N, &ptsperblock, 
                            &numchan, &dt, &T, 1);
         /* Read the first header file and generate an infofile from it */
         chkfread(&hdr, 1, HDRLEN, infiles[0]);
         rewind(infiles[0]);
         PKMB_hdr_to_inf(&hdr, &idata);
         PKMB_update_infodata(numfiles, &idata);
      }

      /* Set-up values if we are using the GMRT Phased Array system */
      if (cmd->gmrtP) {
         printf("GMRT input file information:\n");
         get_GMRT_file_info(infiles, argv + 1, numfiles, cmd->clip, &N, &ptsperblock,
                            &numchan, &dt, &T, 1);
         /* Read the first header file and generate an infofile from it */
         GMRT_hdr_to_inf(argv[1], &idata);
         GMRT_update_infodata(numfiles, &idata);
         set_GMRT_padvals(padvals, good_padvals);
         /* OBS code for TEMPO for the GMRT */
         strcpy(obs, "GM");
      }

      if (cmd->filterbankP) {
         int headerlen;
         sigprocfb fb;

         /* Read the first header file and generate an infofile from it */
         rewind(infiles[0]);
         headerlen = read_filterbank_header(&fb, infiles[0]);
         sigprocfb_to_inf(&fb, &idata);
         rewind(infiles[0]);
         printf("SIGPROC filterbank input file information:\n");
         get_filterbank_file_info(infiles, numfiles, cmd->clip,
                                  &N, &ptsperblock, &numchan, &dt, &T, 1);
         filterbank_update_infodata(numfiles, &idata);
         set_filterbank_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using the Berkeley-Caltech */
      /* Pulsar Machine (or BPP) format.                    */
      if (cmd->bcpmP) {
         printf("BCPM input file information:\n");
         get_BPP_file_info(infiles, numfiles, cmd->clip, &N, &ptsperblock, &numchan,
                           &dt, &T, &idata, 1);
         BPP_update_infodata(numfiles, &idata);
         set_BPP_padvals(padvals, good_padvals);
      }

      /* Set-up values if we are using the NRAO-Caltech Spigot card */
      if (cmd->spigotP) {
         SPIGOT_INFO *spigots;

         printf("Spigot card input file information:\n");
         spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
         for (ii = 0; ii < numfiles; ii++)
            read_SPIGOT_header(cmd->argv[ii], spigots + ii);
         get_SPIGOT_file_info(infiles, spigots, numfiles, cmd->windowP, cmd->clip,
                              &N, &ptsperblock, &numchan, &dt, &T, &idata, 1);
         SPIGOT_update_infodata(numfiles, &idata);
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

      /* Set-up values if we are using the Arecibo WAPP */
      if (cmd->wappP) {
         printf("WAPP input file information:\n");
         get_WAPP_file_info(infiles, cmd->numwapps, numfiles, cmd->windowP,
                            cmd->clip, &N, &ptsperblock, &numchan, &dt, &T, &idata,
                            1);
         WAPP_update_infodata(numfiles, &idata);
         set_WAPP_padvals(padvals, good_padvals);
      }

      /* Finish setting up stuff common to all raw formats */
      idata.dm = cmd->dm;
      worklen = ptsperblock;
      if (cmd->maskfileP)
         maskchans = gen_ivect(idata.num_chan);
      /* Compare the size of the data to the size of output we request */
      if (cmd->numoutP) {
         dtmp = idata.N;
         idata.N = cmd->numout;
         writeinf(&idata);
         idata.N = dtmp;
      } else {
         cmd->numout = INT_MAX;
         writeinf(&idata);
      }
      /* The number of topo to bary time points to generate with TEMPO */
      numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;

      // Identify the TEMPO observatory code
      {
          char *outscope = (char *) calloc(40, sizeof(char));
          telescope_to_tempocode(idata.telescope, outscope, obs);
          free(outscope);
      }
   }

   /* Determine our initialization data if we do _not_ have Parkes, */
   /* Green Bank BCPM, or Arecibo WAPP data sets.                   */
   if (!RAWDATA) {

      /* If we will be barycentering... */
      if (!cmd->nobaryP) {

         /* The number of topo to bary time points to generate with TEMPO */
         numbarypts = (long) (idata.dt * idata.N * 1.1 / TDT + 5.5) + 1;
         
         // Identify the TEMPO observatory code
         {
             char *outscope = (char *) calloc(40, sizeof(char));
             telescope_to_tempocode(idata.telescope, outscope, obs);
             free(outscope);
         }

      }

      /* The number of data points to work with at a time */
      if (worklen > idata.N)
          worklen = idata.N;
      worklen = (long) (worklen / 1024) * 1024;

      /* Compare the size of the data to the size of output we request */
      if (!cmd->numoutP)
         cmd->numout = INT_MAX;
   }
   
   /* Check if we are downsampling */
   dsdt = idata.dt * cmd->downsamp;
   if (cmd->downsamp > 1) {
      printf("Downsampling by a factor of %d\n", cmd->downsamp);
      printf("New sample dt = %.10g\n\n", dsdt);
      if (worklen % cmd->downsamp) {
         printf("Error:  The downsample factor (%d) must be a factor of the\n",
                cmd->downsamp);
         printf("        worklength (%ld).  Exiting.\n\n", worklen);
         exit(1);
      }
   }
   printf("Writing output data to '%s'.\n", datafilenm);
   printf("Writing information to '%s'.\n\n", outinfonm);

   /* The topocentric epoch of the start of the data */
   tlotoa = (double) idata.mjd_i + idata.mjd_f;

   if (!strcmp(idata.band, "Radio") && RAWDATA) {

      /* The topocentric spacing between channels */
      tdf = idata.chan_wid;
      numchan = idata.num_chan;

      /* The topocentric observation frequencies */
      tobsf = gen_dvect(numchan);
      tobsf[0] = idata.freq;
      for (ii = 0; ii < numchan; ii++)
         tobsf[ii] = tobsf[0] + ii * tdf;

      /* The dispersion delays (in time bins) */
      dispdt = gen_dvect(numchan);

      if (cmd->nobaryP) {

         /* Determine our dispersion time delays for each channel */
         for (ii = 0; ii < numchan; ii++)
            dispdt[ii] = delay_from_dm(cmd->dm, tobsf[ii]);

         /* The highest frequency channel gets no delay                 */
         /* All other delays are positive fractions of bin length (dt)  */
         dtmp = dispdt[numchan - 1];
         for (ii = 0; ii < numchan; ii++)
            dispdt[ii] = (dispdt[ii] - dtmp) / idata.dt;
         worklen *= ((int) (fabs(dispdt[0])) / worklen) + 1;
      }

   } else {                     /* For unknown radio raw data (Why is this here?) */
      tobsf = gen_dvect(numchan);
      dispdt = gen_dvect(numchan);
      dispdt[0] = 0.0;
      if (!strcmp(idata.band, "Radio")) {
         tobsf[0] = idata.freq + (idata.num_chan - 1) * idata.chan_wid;
         cmd->dm = idata.dm;
      } else {
         tobsf[0] = 0.0;
         cmd->dm = 0.0;
      }
   }

   if (cmd->nobaryP) {          /* Main loop if we are not barycentering... */

      /* Allocate our data array */
      outdata = gen_fvect(worklen);

      printf("Massaging the data ...\n\n");
      printf("Amount Complete = 0%%");

      do {

         if (cmd->pkmbP)
            numread = read_PKMB(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked, &obsmask);
         else if (cmd->bcpmP)
            numread = read_BPP(infiles, numfiles, outdata, worklen,
                               dispdt, &padding, maskchans, &nummasked,
                               &obsmask, ifs);
         else if (cmd->spigotP)
            numread = read_SPIGOT(infiles, numfiles, outdata, worklen,
                                  dispdt, &padding, maskchans, &nummasked,
                                  &obsmask, ifs);
         else if (cmd->psrfitsP)
            numread = read_PSRFITS(outdata, worklen, dispdt, &padding, 
                                   maskchans, &nummasked, &obsmask);
         else if (cmd->wappP)
            numread = read_WAPP(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked,
                                &obsmask, ifs);
         else if (cmd->gmrtP)
            numread = read_GMRT(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked, &obsmask);
         else if (cmd->filterbankP)
            numread = read_filterbank(infiles, numfiles, outdata, worklen,
                                      dispdt, &padding, maskchans, &nummasked,
                                      &obsmask);
         else
            numread = read_floats(infiles[0], outdata, worklen, numchan);
         if (numread == 0)
            break;

         /* Downsample if requested */
         if (cmd->downsamp > 1)
            numread = downsample(outdata, numread, cmd->downsamp);

         /* Print percent complete */
         if (cmd->numoutP)
            newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
         else
            newper = (int) ((float) totwrote * cmd->downsamp * 100.0 / idata.N) + 1;
         if (newper > oldper) {
            printf("\rAmount Complete = %3d%%", newper);
            fflush(stdout);
            oldper = newper;
         }

         /* Write the latest chunk of data, but don't   */
         /* write more than cmd->numout points.         */
         numtowrite = numread;
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         chkfwrite(outdata, sizeof(float), numtowrite, outfile);
         totwrote += numtowrite;

         /* Update the statistics */
         if (!padding) {
            for (ii = 0; ii < numtowrite; ii++)
               update_stats(statnum + ii, outdata[ii], &min, &max, &avg, &var);
            statnum += numtowrite;
         }

         /* Stop if we have written out all the data we need to */
         if (cmd->numoutP && (totwrote == cmd->numout))
            break;

      } while (numread);

      datawrote = totwrote;

   } else {                     /* Main loop if we are barycentering... */

      double avgvoverc = 0.0, maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;
      double *bobsf = NULL, *btoa = NULL, *ttoa = NULL;

      /* What ephemeris will we use?  (Default is DE200) */
      if (cmd->de405P)
         strcpy(ephem, "DE405");
      else
         strcpy(ephem, "DE200");

      /* Define the RA and DEC of the observation */
      ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
      ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

      /* Allocate some arrays */
      bobsf = gen_dvect(numchan);
      btoa = gen_dvect(numbarypts);
      ttoa = gen_dvect(numbarypts);
      voverc = gen_dvect(numbarypts);
      for (ii = 0; ii < numbarypts; ii++)
         ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

      /* Call TEMPO for the barycentering */
      printf("Generating barycentric corrections...\n");
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
      printf("Collecting and barycentering %s...\n\n", cmd->argv[0]);

      /* Determine the initial dispersion time delays for each channel */
      for (ii = 0; ii < numchan; ii++) {
         bobsf[ii] = doppler(tobsf[ii], avgvoverc);
         dispdt[ii] = delay_from_dm(cmd->dm, bobsf[ii]);
      }

      /* The highest frequency channel gets no delay                   */
      /* All other delays are positive fractions of bin length (dt)    */
      barydispdt = dispdt[numchan - 1];
      for (ii = 0; ii < numchan; ii++)
         dispdt[ii] = (dispdt[ii] - barydispdt) / idata.dt;
      if (RAWDATA)
         worklen *= ((int) (dispdt[0]) / worklen) + 1;

      /* If the data is de-dispersed radio data... */
      if (!strcmp(idata.band, "Radio")) {
         printf("The DM of %.2f at the barycentric observing freq of %.3f MHz\n",
                idata.dm, bobsf[numchan - 1]);
         printf("   causes a delay of %f seconds compared to infinite freq.\n",
                barydispdt);
         printf("   This delay is removed from the barycented times.\n\n");
      }
      printf("Topocentric epoch (at data start) is:\n");
      printf("   %17.11f\n\n", tlotoa);
      printf("Barycentric epoch (infinite obs freq at data start) is:\n");
      printf("   %17.11f\n\n", blotoa - (barydispdt / SECPERDAY));

      /* Convert the bary TOAs to differences from the topo TOAs in  */
      /* units of bin length (dsdt) rounded to the nearest integer.  */
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

      printf("Massaging the data...\n\n");
      printf("Amount Complete = 0%%");

      /* Allocate our data array */

      outdata = gen_fvect(worklen);

      do {                      /* Loop to read and write the data */
         int numwritten = 0;
         double block_avg, block_var;

         if (cmd->pkmbP)
            numread = read_PKMB(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked, &obsmask);
         else if (cmd->bcpmP)
            numread = read_BPP(infiles, numfiles, outdata, worklen,
                               dispdt, &padding, maskchans, &nummasked,
                               &obsmask, ifs);
         else if (cmd->spigotP)
            numread = read_SPIGOT(infiles, numfiles, outdata, worklen,
                                  dispdt, &padding, maskchans, &nummasked,
                                  &obsmask, ifs);
         else if (cmd->psrfitsP)
            numread = read_PSRFITS(outdata, worklen, dispdt, &padding, 
                                   maskchans, &nummasked, &obsmask);
         else if (cmd->wappP)
            numread = read_WAPP(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked,
                                &obsmask, ifs);
         else if (cmd->gmrtP)
            numread = read_GMRT(infiles, numfiles, outdata, worklen,
                                dispdt, &padding, maskchans, &nummasked, &obsmask);
         else if (cmd->filterbankP)
            numread = read_filterbank(infiles, numfiles, outdata, worklen,
                                      dispdt, &padding, maskchans, &nummasked,
                                      &obsmask);
         else if (useshorts)
            numread = read_shorts(infiles[0], outdata, worklen, numchan);
         else
            numread = read_floats(infiles[0], outdata, worklen, numchan);
         if (numread == 0)
            break;

         /* Downsample if requested */
         if (cmd->downsamp > 1)
            numread = downsample(outdata, numread, cmd->downsamp);

         /* Determine the approximate local average */
         avg_var(outdata, numread, &block_avg, &block_var);

         /* Print percent complete */

         if (cmd->numoutP)
            newper = (int) ((float) totwrote / cmd->numout * 100.0) + 1;
         else
            newper = (int) ((float) totwrote * cmd->downsamp * 100.0 / idata.N) + 1;
         if (newper > oldper) {
            printf("\rAmount Complete = %3d%%", newper);
            fflush(stdout);
            oldper = newper;
         }

         /* Simply write the data if we don't have to add or */
         /* remove any bins from this batch.                 */
         /* OR write the amount of data up to cmd->numout or */
         /* the next bin that will be added or removed.      */

         numtowrite = abs(*diffbinptr) - datawrote;
         /* FIXME: numtowrite+totwrote can wrap! */
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         if (numtowrite > numread)
            numtowrite = numread;
         chkfwrite(outdata, sizeof(float), numtowrite, outfile);
         datawrote += numtowrite;
         totwrote += numtowrite;
         numwritten += numtowrite;

         /* Update the statistics */

         if (!padding) {
            for (ii = 0; ii < numtowrite; ii++)
               update_stats(statnum + ii, outdata[ii], &min, &max, &avg, &var);
            statnum += numtowrite;
         }

         if ((datawrote == abs(*diffbinptr)) && (numwritten != numread) && (totwrote < cmd->numout)) {  /* Add/remove a bin */
            float favg;
            int skip, nextdiffbin;

            skip = numtowrite;

            do {                /* Write the rest of the data after adding/removing a bin  */

               if (*diffbinptr > 0) {

                  /* Add a bin */

                  favg = (float) block_avg;
                  chkfwrite(&favg, sizeof(float), 1, outfile);
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
               chkfwrite(outdata + skip, sizeof(float), numtowrite, outfile);
               numwritten += numtowrite;
               datawrote += numtowrite;
               totwrote += numtowrite;

               /* Update the statistics and counters */

               if (!padding) {
                  for (ii = 0; ii < numtowrite; ii++)
                     update_stats(statnum + ii, outdata[skip + ii],
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

      } while (numread);

      /* Free the arrays used in barycentering */

      vect_free(bobsf);
      vect_free(btoa);
      vect_free(ttoa);
   }

   /* Calculate what the amount of padding we need  */

   if (cmd->numoutP && (cmd->numout > totwrote))
      padwrote = padtowrite = cmd->numout - totwrote;


   /* Write the new info file for the output data */

   if (!cmd->nobaryP) {
      idata.bary = 1;
      idata.mjd_i = (int) floor(blotoa - (barydispdt / SECPERDAY));
      idata.mjd_f = blotoa - (barydispdt / SECPERDAY) - idata.mjd_i;
   }
   if (cmd->downsamp > 1)
      idata.dt = dsdt;
   update_infodata(&idata, totwrote, padtowrite, diffbins, numdiffbins);
   writeinf(&idata);

   /* Set the padded points equal to the average data point */

   if (idata.numonoff >= 1) {
      int jj, index, startpad, endpad;

      for (ii = 0; ii < worklen; ii++)
         outdata[ii] = avg;
      fclose(outfile);
      outfile = chkfopen(datafilenm, "rb+");
      for (ii = 0; ii < idata.numonoff; ii++) {
         index = 2 * ii;
         startpad = idata.onoff[index + 1];
         if (ii == idata.numonoff - 1)
            endpad = idata.N - 1;
         else
            endpad = idata.onoff[index + 2];
         chkfseek(outfile, (startpad + 1) * sizeof(float), SEEK_SET);
         padtowrite = endpad - startpad;
         for (jj = 0; jj < padtowrite / worklen; jj++)
            chkfwrite(outdata, sizeof(float), worklen, outfile);
         chkfwrite(outdata, sizeof(float), padtowrite % worklen, outfile);
      }
   }
   vect_free(outdata);

   /* Close the files */

   for (ii = 0; ii < numfiles; ii++)
      fclose(infiles[ii]);
   free(infiles);
   fclose(outfile);

   /* Print simple stats and results */

   var /= (datawrote - 1);

   /* Conver the '.dat' file to '.sdat' if requested */

   if (cmd->shortsP) {
      FILE *infile;
      int safe_convert = 1, bufflen = 65536;
      char *sdatafilenm;
      float *fbuffer;
      short *sbuffer;

      offset = (int) (floor(avg));
      if ((max - min) > (SHRT_MAX - SHRT_MIN)) {
         if ((max - min) < 1.5 * (SHRT_MAX - SHRT_MIN)) {
            printf("Warning:  There is more dynamic range in the data\n"
                   "          than can be handled perfectly:\n"
                   "               max - min = %.2f - %.2f = %.2f\n"
                   "          Clipping the low values...\n\n", max, min, max - min);
            offset = max - SHRT_MAX;
         } else {
            printf("Error:  There is way too much dynamic range in the data:\n"
                   "               max - min = %.2f - %.2f = %.2f\n"
                   "        Not converting to shorts.\n\n", max, min, max - min);
            safe_convert = 0;
         }
      }

      if (safe_convert) {
         fbuffer = gen_fvect(bufflen);
         sbuffer = gen_svect(bufflen);
         sdatafilenm = (char *) calloc(slen, 1);
         sprintf(sdatafilenm, "%s.sdat", cmd->outfile);
         printf("\n\nConverting floats in '%s' to shorts in '%s'.",
                datafilenm, sdatafilenm);
         fflush(NULL);

         infile = chkfopen(datafilenm, "rb");
         outfile = chkfopen(sdatafilenm, "wb");
         while ((numread = chkfread(fbuffer, sizeof(float), bufflen, infile))) {
            for (ii = 0; ii < numread; ii++)
               sbuffer[ii] = (short) (fbuffer[ii] + 1e-20 - offset);
            chkfwrite(sbuffer, sizeof(short), numread, outfile);
         }
         fclose(infile);
         fclose(outfile);
         remove(datafilenm);
         vect_free(fbuffer);
         vect_free(sbuffer);
      }
   }

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
   if (cmd->shortsP && offset != 0)
      printf("          Offset applied to data:  %d\n", -offset);
   printf("\n");

   /* Cleanup */

   if (cmd->maskfileP) {
      free_mask(obsmask);
      vect_free(maskchans);
   }
   vect_free(tobsf);
   vect_free(dispdt);
   free(outinfonm);
   free(datafilenm);
   if (!cmd->nobaryP)
      vect_free(diffbins);
   return (0);
}


static int read_floats(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains normal floating      */
/* point data.                                              */
/* It returns the number of points read.                    */
{
   /* Read the raw data and return numbar read */

   return chkfread(data, sizeof(float), (numpts * numchan), file) / numchan;
}


static int read_shorts(FILE * file, float *data, int numpts, int numchan)
/* This routine reads a numpts records of numchan each from */
/* the input file *file which contains short integer data.  */
/* The equivalent floats are placed in *data.               */
/* It returns the number of points read.                    */
{
   short *sdata;
   int ii, numread;

   sdata = (short *) malloc((size_t) (sizeof(short) * (numpts * numchan)));
   if (!sdata) {
      perror("\nError allocating short array in read_shorts()");
      printf("\n");
      exit(-1);
   }
   numread = chkfread(sdata, sizeof(short),
                      (unsigned long) (numpts * numchan), file) / numchan;
   for (ii = 0; ii < numread; ii++)
      data[ii] = (float) sdata[ii];
   free(sdata);
   return numread;
}

static int downsample(float outdata[], int numread, int downsampfact)
/* Downsample the floating point data by a factor downsampfact */
{
   float tmpout;
   int ii, jj, index;

   numread /= downsampfact;
   for (ii = 0, index = 0; ii < numread; ii++, index += downsampfact) {
      tmpout = 0.0;
      for (jj = 0; jj < downsampfact; jj++)
         tmpout += outdata[index + jj];
      outdata[ii] = tmpout;
   }
   return numread;
}


static void update_infodata(infodata * idata, long datawrote, long padwrote,
                            int *barybins, int numbarybins)
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
