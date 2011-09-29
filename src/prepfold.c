#include <ctype.h>
#include "prepfold.h"
#include "prepfold_cmd.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "gmrt.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "psrfits.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP || cmd->gmrtP || cmd->spigotP || cmd->filterbankP || cmd->psrfitsP)

extern int getpoly(double mjd, double duration, double *dm, FILE * fp, char *pname);
extern int phcalc(double mjd0, double mjd1, int last_index,
                  double *phase, double *psrfreq);
extern int get_psr_from_parfile(char *parfilenm, double epoch, psrparams * psr);
extern char *make_polycos(char *parfilenm, infodata * idata);
void set_posn(prepfoldinfo * in, infodata * idata);

/* 
 * The main program 
 */

int main(int argc, char *argv[])
{
   FILE **infiles, *filemarker;
   float *data = NULL, *padvals, *ppdot = NULL;
   double f = 0.0, fd = 0.0, fdd = 0.0, foldf = 0.0, foldfd = 0.0, foldfdd = 0.0;
   double recdt = 0.0, barydispdt = 0.0, N = 0.0, T = 0.0, proftime, startTday = 0.0;
   double polyco_phase = 0.0, polyco_phase0 = 0.0;
   double *obsf = NULL, *dispdts = NULL, *parttimes = NULL, *Ep = NULL, *tp = NULL;
   double *barytimes = NULL, *topotimes = NULL, *bestprof, dtmp;
   double *buffers, *phasesadded, *events = NULL, orig_foldf = 0.0;
   char *plotfilenm, *outfilenm, *rootnm, *root, *suffix;
   char obs[3], ephem[6], pname[30], rastring[50], decstring[50];
   int numfiles, numevents, numchan = 1, binary = 0, numdelays = 0, numbarypts = 0;
   int info, ptsperrec = 1, flags = 1, padding = 0, arrayoffset = 0, useshorts = 0;
   int *maskchans = NULL, nummasked = 0, polyco_index = 0, insubs = 0, good_padvals =
       0;
   long ii = 0, jj, kk, worklen = 0, numread = 0, reads_per_part = 0;
   long totnumfolded = 0, lorec = 0, hirec = 0, numbinpoints = 0;
   unsigned long numrec = 0;
   IFs ifs = SUMIFS;
   infodata idata;
   foldstats beststats;
   prepfoldinfo search;
   Cmdline *cmd;
   plotflags pflags;
   mask obsmask;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(0);
   }
   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);
   if (cmd->noclipP)
      cmd->clip = 0.0;
   obsmask.numchan = obsmask.numint = 0;
   if (cmd->timingP) {
      cmd->nosearchP = 1;
      cmd->nopsearchP = 1;
      cmd->nopdsearchP = 1;
      cmd->nodmsearchP = 1;
      if (cmd->npart == 64)     /* The default value */
         cmd->npart = 60;
      cmd->fineP = 1;
   }
   if (cmd->slowP) {
      cmd->fineP = 1;
      if (!cmd->proflen) {
         cmd->proflenP = 1;
         cmd->proflen = 100;
      }
      if (cmd->nsub == 32)      /* The default value */
         cmd->nsub = 16;
   }
   if (cmd->fineP) {
      cmd->ndmfact = 1;
      cmd->dmstep = 1;
      cmd->npfact = 1;
      cmd->pstep = 1;
      cmd->pdstep = 2;
   }
   if (cmd->coarseP) {
      cmd->npfact = 4;
      if (cmd->pstep == 1)
         cmd->pstep = 2;
      else
         cmd->pstep = 3;
      if (cmd->pdstep == 2)
         cmd->pdstep = 4;
      else
         cmd->pdstep = 6;
   }
   /* Which IFs will we use? */
   if (cmd->ifsP) {
      if (cmd->ifs == 0)
         ifs = IF0;
      else if (cmd->ifs == 1)
         ifs = IF1;
      else
         ifs = SUMIFS;
   }
   pflags.events = cmd->eventsP;
   pflags.nosearch = cmd->nosearchP;
   pflags.scaleparts = cmd->scalepartsP;
   pflags.justprofs = cmd->justprofsP;
   pflags.allgrey = cmd->allgreyP;

#ifdef DEBUG
   showOptionValues();
#endif

   printf("\n\n");
   printf("        Pulsar Raw-Data Folding Search Routine\n");
   printf(" Used for DM, Period, and P-dot tweaking of PSR candidates.\n");
   printf("                 by Scott M. Ransom\n\n");

   init_prepfoldinfo(&search);

   numfiles = cmd->argc;
   infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
   {
      char *path;

      split_path_file(cmd->argv[0], &path, &search.filenm);
      free(path);

      if (!cmd->outfileP) {
         char *tmprootnm;

         split_root_suffix(cmd->argv[0], &tmprootnm, &suffix);
         if ((cmd->startT != 0.0) || (cmd->endT != 1.0)) {
            rootnm = (char *) calloc(strlen(tmprootnm) + 11, sizeof(char));
            sprintf(rootnm, "%s_%4.2f-%4.2f", tmprootnm, cmd->startT, cmd->endT);
         } else {
            rootnm = tmprootnm;
         }
         free(suffix);
      } else {
         rootnm = cmd->outfile;
      }
   }

   if (!RAWDATA) {
      /* Split the filename into a rootname and a suffix */
      if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
         printf("\nThe input filename (%s) must have a suffix!\n\n", search.filenm);
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
         } else if (strcmp(suffix, "sub0") == 0 ||
                    strcmp(suffix, "sub00") == 0 ||
                    strcmp(suffix, "sub000") == 0 ||
                    strcmp(suffix, "sub0000") == 0) {
            printf
                ("Assuming the data is subband data generated by prepsubband...\n");
            insubs = 1;
            useshorts = 1;
         }
      }
      if (RAWDATA) {            /* Clean-up a bit */
         free(root);
         free(suffix);
      }
   }

   if (!RAWDATA) {
      if (insubs) {
         char *tmpname;
         if (strncmp(suffix, "sub", 3) == 0) {
            tmpname = calloc(strlen(root) + 6, 1);
            sprintf(tmpname, "%s.sub", root);
            readinf(&idata, tmpname);
            free(tmpname);
         } else {
            printf("\nThe input files (%s) must be subbands!  (i.e. *.sub##)\n\n",
                   cmd->argv[0]);
            exit(1);
         }
         if (numfiles > 1)
            printf("Reading subband data from %d files:\n", numfiles);
         else
            printf("Reading subband data from 1 file:\n");
      } else {
         printf("Reading input data from '%s'.\n", cmd->argv[0]);
         printf("Reading information from '%s.inf'.\n\n", root);
         /* Read the info file if available */
         readinf(&idata, root);
      }
      free(root);
      free(suffix);
      /* Use events instead of a time series */
      if (cmd->eventsP) {
         int eventtype = 0;     /* 0=sec since .inf, 1=days since .inf, 2=MJDs */

         /* The following allows using inf files from searches of a subset */
         /* of events from an event file.                                  */
         if (cmd->rzwcandP || cmd->accelcandP) {
            infodata rzwidata;
            char *cptr = NULL;

            if (cmd->rzwcandP) {
               if (!cmd->rzwfileP) {
                  printf("\nYou must enter a name for the rzw candidate ");
                  printf("file (-rzwfile filename)\n");
                  printf("Exiting.\n\n");
                  exit(1);
               } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))) {
                  ii = (long) (cptr - cmd->rzwfile);
               } else if (NULL != (cptr = strstr(cmd->rzwfile, "_ACCEL"))) {
                  ii = (long) (cptr - cmd->rzwfile);
               }
               cptr = (char *) calloc(ii + 1, sizeof(char));
               strncpy(cptr, cmd->rzwfile, ii);
            }
            if (cmd->accelcandP) {
               if (!cmd->accelfileP) {
                  printf("\nYou must enter a name for the ACCEL candidate ");
                  printf("file (-accelfile filename)\n");
                  printf("Exiting.\n\n");
                  exit(1);
               } else if (NULL != (cptr = strstr(cmd->accelfile, "_rzw"))) {
                  ii = (long) (cptr - cmd->accelfile);
               } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))) {
                  ii = (long) (cptr - cmd->accelfile);
               }
               cptr = (char *) calloc(ii + 1, sizeof(char));
               strncpy(cptr, cmd->accelfile, ii);
            }
            readinf(&rzwidata, cptr);
            free(cptr);
            idata.mjd_i = rzwidata.mjd_i;
            idata.mjd_f = rzwidata.mjd_f;
            idata.N = rzwidata.N;
            idata.dt = rzwidata.dt;
         }
         printf("Assuming the events are barycentered or geocentered.\n");
         if (!cmd->proflenP) {
            cmd->proflenP = 1;
            cmd->proflen = 20;
            printf("Using %d bins in the profile since not specified.\n",
                   cmd->proflen);
         }
         if (cmd->doubleP)
            infiles[0] = chkfopen(cmd->argv[0], "rb");
         else
            infiles[0] = chkfopen(cmd->argv[0], "r");
         if (cmd->daysP)
            eventtype = 1;
         else if (cmd->mjdsP)
            eventtype = 2;
         events = read_events(infiles[0], cmd->doubleP, eventtype, &numevents,
                              idata.mjd_i + idata.mjd_f, idata.N * idata.dt,
                              cmd->startT, cmd->endT, cmd->offset);
         if (cmd->rzwcandP || cmd->accelcandP)
            T = idata.N * idata.dt;
         else {
            /* The 1e-8 prevents floating point rounding issues from
               causing the last event to go into slice npart+1.
               Thanks to Paul Ray for finding this. */
            T = events[numevents - 1] + 1e-8;
         }
      } else {
         infiles[0] = chkfopen(cmd->argv[0], "rb");
      }
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
   } else if (cmd->wappP) {
      if (numfiles > 1)
         printf("Reading Arecibo WAPP data from %d files:\n", numfiles);
      else
         printf("Reading Arecibo WAPP data from 1 file:\n");
   } else if (cmd->gmrtP) {
      if (numfiles > 1)
         printf("Reading GMRT Phased Array data from %d files:\n", numfiles);
      else
         printf("Reading GMRT Phased Array data from 1 file:\n");
   }

   /* Open the raw data files */

   if (RAWDATA || insubs) {
      for (ii = 0; ii < numfiles; ii++) {
         printf("  '%s'\n", cmd->argv[ii]);
         infiles[ii] = chkfopen(cmd->argv[ii], "rb");
      }
      printf("\n");
      /* Read an input mask if wanted */
      if (cmd->maskfileP) {
         read_mask(cmd->maskfile, &obsmask);
         printf("Read mask information from '%s'\n\n", cmd->maskfile);
         good_padvals = determine_padvals(cmd->maskfile, &obsmask, &padvals);
      } else {
         obsmask.numchan = obsmask.numint = 0;
      }
   }

   /* Manipulate the file names we will use  */

   {
      int slen;

      /* Determine the candidate name */

      if (cmd->psrnameP) {
         search.candnm = (char *) calloc(strlen(cmd->psrname) + 5, sizeof(char));
         sprintf(search.candnm, "PSR_%s", cmd->psrname);
      } else if (cmd->parnameP || cmd->timingP) {
         int retval;
         psrparams psr;
         if (cmd->timingP) {
            cmd->parnameP = 1;
            cmd->parname = cmd->timing;
         }
         /* Read the par file just to get the PSR name */
         retval = get_psr_from_parfile(cmd->parname, 51000.0, &psr);
         search.candnm = (char *) calloc(strlen(psr.jname) + 5, sizeof(char));
         sprintf(search.candnm, "PSR_%s", psr.jname);
      } else if (cmd->rzwcandP) {
         slen = 20;
         search.candnm = (char *) calloc(slen, sizeof(char));
         sprintf(search.candnm, "RZW_Cand_%d", cmd->rzwcand);
      } else if (cmd->accelcandP) {
         slen = 22;
         search.candnm = (char *) calloc(slen, sizeof(char));
         sprintf(search.candnm, "ACCEL_Cand_%d", cmd->accelcand);
      } else {
         slen = 20;
         search.candnm = (char *) calloc(slen, sizeof(char));
         if (cmd->pP)
            sprintf(search.candnm, "%.2fms_Cand", cmd->p * 1000.0);
         else if (cmd->fP)
            sprintf(search.candnm, "%.2fHz_Cand", cmd->f);
         else {
            printf("\nYou must specify candidate parameters (i.e period).\n\n");
            exit(1);
         }
      }

      /* Determine the output and plot file names */

      slen = strlen(rootnm) + strlen(search.candnm) + 6;
      outfilenm = (char *) calloc(slen, sizeof(char));
      sprintf(outfilenm, "%s_%s.pfd", rootnm, search.candnm);
      plotfilenm = (char *) calloc(slen + 3, sizeof(char));
      sprintf(plotfilenm, "%s_%s.pfd.ps", rootnm, search.candnm);
      search.pgdev = (char *) calloc(slen + 7, sizeof(char));
      sprintf(search.pgdev, "%s/CPS", plotfilenm);
   }

   /* What ephemeris will we use?  (Default is DE200) */

   if (cmd->de405P)
      strcpy(ephem, "DE405");
   else
      strcpy(ephem, "DE200");

   /* Set-up values if we are using the Parkes Multibeam System, */
   /* a Berkeley-Caltech Pulsar Machine, or the Arecibo WAPP.    */

   if (RAWDATA || insubs) {
      double local_dt, local_T;
      long long local_N;

      if (cmd->pkmbP) {
         PKMB_tapehdr hdr;

         printf("Filterbank input file information:\n");
         get_PKMB_file_info(infiles, numfiles, cmd->clip, &local_N, 
                            &ptsperrec, &numchan, &local_dt, &local_T, 1);

         /* Read the first header file and generate an infofile from it */

         chkfread(&hdr, 1, HDRLEN, infiles[0]);
         rewind(infiles[0]);
         PKMB_hdr_to_inf(&hdr, &idata);
         PKMB_update_infodata(numfiles, &idata);

      } else if (cmd->bcpmP) {

         printf("BCPM input file information:\n");
         get_BPP_file_info(infiles, numfiles, cmd->clip, &local_N, &ptsperrec,
                           &numchan, &local_dt, &local_T, &idata, 1);
         BPP_update_infodata(numfiles, &idata);
         set_BPP_padvals(padvals, good_padvals);

      } else if (cmd->spigotP) {
         SPIGOT_INFO *spigots;

         printf("Spigot card input file information:\n");
         spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
         for (ii = 0; ii < numfiles; ii++)
            read_SPIGOT_header(cmd->argv[ii], spigots + ii);
         get_SPIGOT_file_info(infiles, spigots, numfiles,
                              cmd->windowP, cmd->clip,
                              &local_N, &ptsperrec, &numchan,
                              &local_dt, &local_T, &idata, 1);
         SPIGOT_update_infodata(numfiles, &idata);
         set_SPIGOT_padvals(padvals, good_padvals);
         free(spigots);

      } else if (cmd->psrfitsP) {
         struct spectra_info s;
         
         printf("PSRFITS input file information:\n");
          // -1 causes the data to determine if we use weights, scales, & offsets
         s.apply_weight = (cmd->noweightsP) ? 0 : -1;
         s.apply_scale  = (cmd->noscalesP) ? 0 : -1;
         s.apply_offset = (cmd->nooffsetsP) ? 0 : -1;
         read_PSRFITS_files(cmd->argv, cmd->argc, &s);
         local_N = s.N;
         ptsperrec = s.spectra_per_subint;
         numchan = s.num_channels;
         local_dt = s.dt;
         local_T = s.T;
         get_PSRFITS_file_info(cmd->argv, cmd->argc, cmd->clip, 
                               &s, &idata, 1);
         PSRFITS_update_infodata(&idata);
         set_PSRFITS_padvals(padvals, good_padvals);

      } else if (cmd->wappP) {

         printf("WAPP input file information:\n");
         get_WAPP_file_info(infiles, cmd->numwapps, numfiles,
                            cmd->windowP, cmd->clip,
                            &local_N, &ptsperrec, &numchan,
                            &local_dt, &local_T, &idata, 1);
         WAPP_update_infodata(numfiles, &idata);
         set_WAPP_padvals(padvals, good_padvals);

      } else if (cmd->gmrtP) {

         printf("GMRT input file information:\n");
         get_GMRT_file_info(infiles, argv + 1, numfiles, cmd->clip,
                            &local_N, &ptsperrec, &numchan, &local_dt, &local_T, 1);
         /* Read the first header file and generate an infofile from it */
         GMRT_hdr_to_inf(argv[1], &idata);
         GMRT_update_infodata(numfiles, &idata);
         set_GMRT_padvals(padvals, good_padvals);

      } else if (cmd->filterbankP) {
         int headerlen;
         sigprocfb fb;

         /* Read the first header file and generate an infofile from it */
         rewind(infiles[0]);
         headerlen = read_filterbank_header(&fb, infiles[0]);
         sigprocfb_to_inf(&fb, &idata);
         rewind(infiles[0]);
         printf("SIGPROC filterbank input file information:\n");
         get_filterbank_file_info(infiles, numfiles, cmd->clip,
                                  &local_N, &ptsperrec, &numchan,
                                  &local_dt, &local_T, 1);
         filterbank_update_infodata(numfiles, &idata);
         set_filterbank_padvals(padvals, good_padvals);

      } else if (insubs) {

         cmd->nsub = numfiles;
         local_N = chkfilelen(infiles[0], sizeof(short));
         local_dt = idata.dt;
         local_T = local_N * local_T;
         numchan = idata.num_chan;
         ptsperrec = SUBSBLOCKLEN;

      }

      // Identify the TEMPO observatory code
      search.telescope = (char *) calloc(40, sizeof(char));
      telescope_to_tempocode(idata.telescope, search.telescope, obs);

      idata.dm = cmd->dm;
      numrec = local_N / ptsperrec;
      if (cmd->maskfileP)
         maskchans = gen_ivect(obsmask.numchan);

      /* Define the RA and DEC of the observation */

      ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
      ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

      /* Define some variables */

      search.dt = idata.dt;
      recdt = search.dt * ptsperrec;

      /* Determine the number of records to use from the command line */

      search.startT = cmd->startT;
      search.endT = cmd->endT;
      lorec = (long) (cmd->startT * numrec + DBLCORRECT);
      hirec = (long) (cmd->endT * numrec + DBLCORRECT);
      startTday = lorec * recdt / SECPERDAY;
      numrec = hirec - lorec;

      /* The number of reads from the file we need for */
      /* each sub-integration.                         */

      reads_per_part = numrec / cmd->npart;

      /* Correct numrec so that each part will contain */
      /* the same number of records.                   */

      numrec = reads_per_part * cmd->npart;
      T = numrec * recdt;
      N = numrec * ptsperrec;

      /* Topocentric and barycentric times of folding epoch data */

      if (idata.mjd_i) {
         search.tepoch = idata.mjd_i + idata.mjd_f + startTday;

         if (!cmd->polycofileP && !cmd->timingP && !cmd->topoP && !cmd->parnameP) {
            barycenter(&search.tepoch, &search.bepoch, &dtmp, 1, rastring,
                       decstring, obs, ephem);

            /* Correct the barycentric time for the dispersion delay.     */
            /* This converts the barycentric time to infinite frequency.  */
            if (cmd->dm > 0.0) {
               barydispdt = delay_from_dm(cmd->dm, idata.freq +
                                          (idata.num_chan - 1) * idata.chan_wid);
               search.bepoch -= (barydispdt / SECPERDAY);
            }
         }
      }
      worklen = ptsperrec;

   } else {  /* Raw floating point or event data (already de-dispersed if radio data) */

      cmd->nsub = 1;
      search.startT = cmd->startT;
      search.endT = cmd->endT;

      if (!cmd->eventsP) {

         /* Some information about the size of the records */

         numchan = 1;
         worklen = SUBSBLOCKLEN;
         search.dt = idata.dt;
         if (useshorts)
            N = chkfilelen(infiles[0], sizeof(short));
         else
            N = chkfilelen(infiles[0], sizeof(float));

         /* Determine the number of records to use from the command line */

         lorec = (long) (cmd->startT * N + DBLCORRECT);
         hirec = (long) (cmd->endT * N + DBLCORRECT);
         startTday = lorec * search.dt / SECPERDAY;
         numrec = (hirec - lorec) / worklen;
         recdt = worklen * search.dt;

         /* The number of reads from the file we need for */
         /* each sub-integration.                         */

         reads_per_part = numrec / cmd->npart;

         /* Correct numrec so that each part will contain */
         /* the same number of records.                   */

         numrec = reads_per_part * cmd->npart;
         N = numrec * worklen;
         T = N * search.dt;
      }

      /* Until I figure out a better way to do this... */

      search.telescope = (char *) calloc(strlen(idata.telescope) + 1, sizeof(char));
      strcpy(search.telescope, idata.telescope);

      if (idata.mjd_i) {
         if (idata.bary)
            search.bepoch = idata.mjd_i + idata.mjd_f + startTday;
         else
            search.tepoch = idata.mjd_i + idata.mjd_f + startTday;
      }
   }

   /* Make sure that the number of subbands evenly divides the number of channels */
   if (numchan % cmd->nsub != 0) {
       printf("Error:  # of channels (%d) not divisible by # of subbands (%d)!\n", 
              numchan, cmd->nsub);
       exit(1);
   }

   set_posn(&search, &idata);
   printf("Folding a %s candidate.\n\n", search.candnm);
   printf("Output data file is '%s'.\n", outfilenm);
   printf("Output plot file is '%s'.\n", plotfilenm);
   printf("Best profile is in  '%s.bestprof'.\n", outfilenm);

   /* Generate polycos if required and set the pulsar name */
   if ((cmd->timingP || cmd->parnameP) && !idata.bary) {
      char *polycofilenm;
      cmd->psrnameP = 1;
      if (cmd->timingP)
         cmd->psrname = make_polycos(cmd->timing, &idata);
      else
         cmd->psrname = make_polycos(cmd->parname, &idata);
      polycofilenm = (char *) calloc(strlen(outfilenm) + 9, sizeof(char));
      sprintf(polycofilenm, "%s.polycos", outfilenm);
      printf("Polycos used are in '%s'.\n", polycofilenm);
      rename("polyco.dat", polycofilenm);
      cmd->polycofileP = 1;
      cmd->polycofile = (char *) calloc(strlen(polycofilenm) + 1, sizeof(char));
      strcpy(cmd->polycofile, polycofilenm);
      free(polycofilenm);
   }

   /* Read the pulsar database if needed */
   if (cmd->psrnameP) {
      if (cmd->polycofileP) {
         FILE *polycofileptr;
         int numsets;
         double polyco_dm;

         polycofileptr = chkfopen(cmd->polycofile, "r");
         numsets = getpoly(search.tepoch, T / SECPERDAY, &polyco_dm,
                           polycofileptr, cmd->psrname);
         fclose(polycofileptr);
         if (cmd->dm > 0.0) {
            printf("\nRead %d set(s) of polycos for PSR %s at %18.12f\n",
                   numsets, cmd->psrname, search.tepoch);
            printf("Overriding polyco DM = %f with %f\n", polyco_dm, cmd->dm);
         } else {
            printf("\nRead %d set(s) of polycos for PSR %s at %18.12f (DM = %.5g)\n",
                   numsets, cmd->psrname, search.tepoch, polyco_dm);
            cmd->dm = polyco_dm;
         }
         polyco_index = phcalc(idata.mjd_i, idata.mjd_f + startTday,
                               polyco_index, &polyco_phase0, &f);
         search.topo.p1 = 1.0 / f;
         search.topo.p2 = fd = 0.0;
         search.topo.p3 = fdd = 0.0;
         strcpy(pname, cmd->psrname);
      } else {                  /* Use the database */
         int pnum;
         psrparams psr;

         if (search.bepoch == 0.0) {
            printf("\nYou cannot fold topocentric data with the pulsar database.\n");
            printf("Use '-timing' or polycos instead.  Exiting.\n\n");
            exit(1);
         }
         pnum = get_psr_at_epoch(cmd->psrname, search.bepoch, &psr);
         if (!pnum) {
            printf("The pulsar is not in the database.  Exiting.\n\n");
            exit(1);
         }
         if (psr.orb.p != 0.0) {        /* Checks if the pulsar is in a binary */
            binary = 1;
            search.orb = psr.orb;
         }
         search.bary.p1 = psr.p;
         search.bary.p2 = psr.pd;
         search.bary.p3 = psr.pdd;
         if (cmd->dm == 0.0)
            cmd->dm = psr.dm;
         f = psr.f;
         fd = psr.fd;
         fdd = psr.fdd;
         strcpy(pname, psr.jname);
      }

   } else if (cmd->parnameP) {  /* Read ephemeris from a par file */
      int retval;
      psrparams psr;

      if (search.bepoch == 0.0) {
         printf("\nYou cannot fold topocentric data with a par file.\n");
         printf("Use '-timing' or polycos instead.  Exiting.\n\n");
         exit(1);
      }

      /* Read the par file */
      retval = get_psr_from_parfile(cmd->parname, search.bepoch, &psr);

      if (psr.orb.p != 0.0) {   /* Checks if the pulsar is in a binary */
         binary = 1;
         search.orb = psr.orb;
      }
      search.bary.p1 = psr.p;
      search.bary.p2 = psr.pd;
      search.bary.p3 = psr.pdd;
      f = psr.f;
      fd = psr.fd;
      fdd = psr.fdd;
      strcpy(pname, psr.jname);
      if (cmd->dm == 0.0)
         cmd->dm = psr.dm;

      /* If the user specifies all of the binaries parameters */

   } else if (cmd->binaryP) {

      /* Assume that the psr characteristics were measured at the time */
      /* of periastron passage (cmd->To)                               */

      if (search.bepoch == 0.0)
         dtmp = SECPERDAY * (search.tepoch - cmd->To);
      else
         dtmp = SECPERDAY * (search.bepoch - cmd->To);
      search.orb.p = cmd->pb;
      search.orb.x = cmd->asinic;
      search.orb.e = cmd->e;
      search.orb.t = fmod(dtmp, search.orb.p);
      if (search.orb.t < 0.0)
         search.orb.t += search.orb.p;
      search.orb.w = (cmd->w + dtmp * cmd->wdot / SECPERJULYR);
      binary = 1;

   } else if (cmd->rzwcandP || cmd->accelcandP) {
      fourierprops rzwcand;
      infodata rzwidata;
      char *cptr = NULL;

      if (cmd->rzwcandP) {
         if (!cmd->rzwfileP) {
            printf("\nYou must enter a name for the rzw candidate ");
            printf("file (-rzwfile filename)\n");
            printf("Exiting.\n\n");
            exit(1);
         } else if (NULL != (cptr = strstr(cmd->rzwfile, "_rzw"))) {
            ii = (long) (cptr - cmd->rzwfile);
         } else if (NULL != (cptr = strstr(cmd->rzwfile, "_ACCEL"))) {
            ii = (long) (cptr - cmd->rzwfile);
         }
         cptr = (char *) calloc(ii + 1, sizeof(char));
         strncpy(cptr, cmd->rzwfile, ii);
      }
      if (cmd->accelcandP) {
         if (!cmd->accelfileP) {
            printf("\nYou must enter a name for the ACCEL candidate ");
            printf("file (-accelfile filename)\n");
            printf("Exiting.\n\n");
            exit(1);
         } else if (NULL != (cptr = strstr(cmd->accelfile, "_rzw"))) {
            ii = (long) (cptr - cmd->accelfile);
         } else if (NULL != (cptr = strstr(cmd->accelfile, "_ACCEL"))) {
            ii = (long) (cptr - cmd->accelfile);
         }
         cptr = (char *) calloc(ii + 1, sizeof(char));
         strncpy(cptr, cmd->accelfile, ii);
      }
      printf("\nAttempting to read '%s.inf'.  ", cptr);
      readinf(&rzwidata, cptr);
      free(cptr);
      printf("Successful.\n");
      if (cmd->rzwfileP)
         get_rzw_cand(cmd->rzwfile, cmd->rzwcand, &rzwcand);
      if (cmd->accelfileP)
         get_rzw_cand(cmd->accelfile, cmd->accelcand, &rzwcand);
      f = (rzwcand.r - 0.5 * rzwcand.z) / (rzwidata.dt * rzwidata.N);
      fd = rzwcand.z / ((rzwidata.dt * rzwidata.N) * (rzwidata.dt * rzwidata.N));

      /* Now correct for the fact that we may not be starting */
      /* to fold at the same start time as the rzw search.    */

      if (RAWDATA || insubs)
         f += lorec * recdt * fd;
      else
         f += lorec * search.dt * fd;
      if (rzwidata.bary)
         switch_f_and_p(f, fd, fdd, &search.bary.p1,
                        &search.bary.p2, &search.bary.p3);
      else
         switch_f_and_p(f, fd, fdd, &search.topo.p1,
                        &search.topo.p2, &search.topo.p3);
   }

   /* Determine the pulsar parameters to fold if we are not getting   */
   /* the data from a .cand file, the pulsar database, or a makefile. */

   if (!cmd->rzwcandP && !cmd->accelcandP && !cmd->psrnameP && !cmd->parnameP) {
      double p = 0.0, pd = 0.0, pdd = 0.0;

      if (cmd->pP) {
         p = cmd->p;
         f = 1.0 / p;
      }
      if (cmd->fP) {
         f = cmd->f;
         p = 1.0 / f;
      }
      if (cmd->pd != 0.0) {
         pd = cmd->pd;
         fd = -pd / (p * p);
      }
      if (cmd->fd != 0.0) {
         fd = cmd->fd;
         pd = -fd / (f * f);
      }
      if (cmd->pdd != 0.0) {
         pdd = cmd->pdd;
         fdd = 2 * pd * pd / (p * p * p) - pdd / (p * p);
      }
      if (cmd->fdd != 0.0) {
         fdd = cmd->fdd;
         pdd = 2 * fd * fd / (f * f * f) - fdd / (f * f);
      }
      if (idata.bary) {
         search.bary.p1 = p;
         search.bary.p2 = pd;
         search.bary.p3 = pdd;
      } else {
         search.topo.p1 = p;
         search.topo.p2 = pd;
         search.topo.p3 = pdd;
      }
   }

   /* Modify the fold period or frequency */

   if (cmd->pfact != 1.0 || cmd->ffact != 1.0) {
      if (cmd->pfact == 0.0 || cmd->ffact == 0.0) {
         printf("\nFolding factors (-pfact or -ffact) cannot be 0!  Exiting\n");
         exit(1);
      }
      if (cmd->pfact != 1.0)
         cmd->ffact = 1.0 / cmd->pfact;
      else if (cmd->ffact != 1.0)
         cmd->pfact = 1.0 / cmd->ffact;
      search.bary.p1 *= cmd->pfact;
      search.bary.p2 /= cmd->pfact;
      search.topo.p1 *= cmd->pfact;
      search.topo.p2 /= cmd->pfact;
      f *= cmd->ffact;
      fd *= cmd->ffact;
      fdd *= cmd->ffact;
   }

   /* Determine the length of the profile */

   if (cmd->proflenP) {
      search.proflen = cmd->proflen;
   } else {
      if (search.topo.p1 == 0.0)
         search.proflen = (long) (search.bary.p1 / search.dt + 0.5);
      else
         search.proflen = (long) (search.topo.p1 / search.dt + 0.5);
      if (cmd->timingP)
         search.proflen = next2_to_n(search.proflen);
      if (search.proflen > 64)
         search.proflen = 64;
   }

   /* Determine the phase delays caused by the orbit if needed */

   if (binary && !cmd->eventsP) {
      double orbdt = 1.0, startE = 0.0;

      /* Save the orbital solution every half second               */
      /* The times in *tp are now calculated as barycentric times. */
      /* Later, we will change them to topocentric times after     */
      /* applying corrections to Ep using TEMPO.                   */

      startE = keplars_eqn(search.orb.t, search.orb.p, search.orb.e, 1.0E-15);
      if (T > 2048)
         orbdt = 0.5;
      else
         orbdt = T / 4096.0;
      numbinpoints = (long) floor(T / orbdt + 0.5) + 1;
      Ep = dorbint(startE, numbinpoints, orbdt, &search.orb);
      tp = gen_dvect(numbinpoints);
      for (ii = 0; ii < numbinpoints; ii++)
         tp[ii] = ii * orbdt;

      /* Convert Eccentric anomaly to time delays */

      search.orb.w *= DEGTORAD;
      E_to_phib(Ep, numbinpoints, &search.orb);
      numdelays = numbinpoints;
      if (search.bepoch == 0.0)
         search.orb.t = -search.orb.t / SECPERDAY + search.tepoch;
      else
         search.orb.t = -search.orb.t / SECPERDAY + search.bepoch;
   }

   if (((RAWDATA || insubs) && !cmd->topoP)
       && cmd->dm == 0.0 && !cmd->polycofileP) {
      /* Correct the barycentric time for the dispersion delay.     */
      /* This converts the barycentric time to infinite frequency.  */
      barydispdt = delay_from_dm(cmd->dm, idata.freq +
                                 (idata.num_chan - 1) * idata.chan_wid);
      search.bepoch -= (barydispdt / SECPERDAY);
   }

   if (cmd->eventsP) {
      search.dt = (search.bary.p1 + 0.5 * T * search.bary.p2) / search.proflen;
      N = ceil(T / search.dt);
   }

   /* Output some informational data on the screen and to the */
   /* output file.                                            */

   fprintf(stdout, "\n");
   filemarker = stdout;
   for (ii = 0; ii < 1; ii++) {
      double p, pd, pdd;

      if (cmd->psrnameP)
         fprintf(filemarker, "Pulsar                       =  %s\n", pname);
      if (search.tepoch != 0.0)
         fprintf(filemarker,
                 "Folding (topo) epoch  (MJD)  =  %-17.12f\n", search.tepoch);
      if (search.bepoch != 0.0)
         fprintf(filemarker,
                 "Folding (bary) epoch  (MJD)  =  %-17.12f\n", search.bepoch);
      fprintf(filemarker, "Data pt duration (dt)   (s)  =  %-.12g\n", search.dt);
      fprintf(filemarker, "Total number of data points  =  %-.0f\n", N);
      fprintf(filemarker, "Number of profile bins       =  %d\n", search.proflen);
      switch_f_and_p(f, fd, fdd, &p, &pd, &pdd);
      fprintf(filemarker, "Folding period          (s)  =  %-.15g\n", p);
      if (pd != 0.0)
         fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10g\n", pd);
      if (pdd != 0.0)
         fprintf(filemarker, "Folding p-dotdot    (s/s^2)  =  %-.10g\n", pdd);
      fprintf(filemarker, "Folding frequency      (hz)  =  %-.12g\n", f);
      if (fd != 0.0)
         fprintf(filemarker, "Folding f-dot        (hz/s)  =  %-.8g\n", fd);
      if (fdd != 0.0)
         fprintf(filemarker, "Folding f-dotdot   (hz/s^2)  =  %-.8g\n", fdd);
      if (binary) {
         double tmpTo;
         fprintf(filemarker,
                 "Orbital period          (s)  =  %-.10g\n", search.orb.p);
         fprintf(filemarker,
                 "a*sin(i)/c (x)     (lt-sec)  =  %-.10g\n", search.orb.x);
         fprintf(filemarker,
                 "Eccentricity                 =  %-.10g\n", search.orb.e);
         fprintf(filemarker,
                 "Longitude of peri (w) (deg)  =  %-.10g\n",
                 search.orb.w / DEGTORAD);
         tmpTo = search.orb.t;
         if (cmd->eventsP) {
            if (search.bepoch == 0.0)
               tmpTo = -search.orb.t / SECPERDAY + search.tepoch;
            else
               tmpTo = -search.orb.t / SECPERDAY + search.bepoch;
         }
         fprintf(filemarker, "Epoch of periapsis    (MJD)  =  %-17.11f\n", tmpTo);
      }
   }

   /* Check to see if the folding period is close to the time in a sub-interval */
   if (!cmd->eventsP) {
      if (T / cmd->npart < 5.0 / f)
         printf
             ("\nWARNING:  The pulse period is close to the duration of a folding\n"
              "  sub-interval.  This may cause artifacts in the plot and/or excessive\n"
              "  loss of data during the fold.  I recommend re-running with -npart set\n"
              "  to a significantly smaller value than the current value of %d.\n",
              cmd->npart);
   } else {
      if (T / cmd->npart < 5.0 / f) {
         cmd->npart = (int) (T * f / 5.0);
         printf("\nOverriding default number of sub-intervals due to low number\n"
                "  of pulses during the observation.  Current npart = %d\n",
                cmd->npart);
      }
   }

   /* Allocate and initialize some arrays and other information */

   search.nsub = cmd->nsub;
   search.npart = cmd->npart;
   search.rawfolds = gen_dvect(cmd->nsub * cmd->npart * search.proflen);
   search.stats = (foldstats *) malloc(sizeof(foldstats) * cmd->nsub * cmd->npart);
   for (ii = 0; ii < cmd->npart * cmd->nsub * search.proflen; ii++)
      search.rawfolds[ii] = 0.0;
   for (ii = 0; ii < cmd->npart * cmd->nsub; ii++) {
      search.stats[ii].numdata = 0.0;
      search.stats[ii].data_avg = 0.0;
      search.stats[ii].data_var = 0.0;
   }
   if (numdelays == 0)
      flags = 0;

   if (cmd->eventsP) {          /* Fold events instead of a time series */
      double event, dtmp, cts, phase, begphs, endphs, dphs, lphs, rphs;
      double tf, tfd, tfdd, totalphs, calctotalphs, numwraps;
      int partnum, binnum;

      foldf = f;
      foldfd = fd;
      foldfdd = fdd;
      search.fold.pow = 1.0;
      search.fold.p1 = f;
      search.fold.p2 = fd;
      search.fold.p3 = fdd;
      tf = f;
      tfd = fd / 2.0;
      tfdd = fdd / 6.0;
      dtmp = cmd->npart / T;
      parttimes = gen_dvect(cmd->npart);
      if (binary)
         search.orb.w *= DEGTORAD;
      for (ii = 0; ii < numevents; ii++) {
         event = events[ii];
         if (binary) {
            double tt, delay;
            tt = fmod(search.orb.t + event, search.orb.p);
            delay = keplars_eqn(tt, search.orb.p, search.orb.e, 1.0E-15);
            E_to_phib(&delay, 1, &search.orb);
            event -= delay;
         }
         partnum = (int) floor(event * dtmp);
         phase = event * (event * (event * tfdd + tfd) + tf);
         binnum = (int) ((phase - (long long) phase) * search.proflen);
         search.rawfolds[partnum * search.proflen + binnum] += 1.0;
      }
      if (binary) {
         if (search.bepoch == 0.0)
            search.orb.t = -search.orb.t / SECPERDAY + search.tepoch;
         else
            search.orb.t = -search.orb.t / SECPERDAY + search.bepoch;
      }
      for (ii = 0; ii < cmd->npart; ii++) {
         parttimes[ii] = (T * ii) / (double) (cmd->npart);
         /* Correct each part for the "exposure".  This gives us a count rate. */
         event = parttimes[ii];
         begphs = event * (event * (event * tfdd + tfd) + tf);
         event = (T * (ii + 1)) / (double) (cmd->npart);
         endphs = event * (event * (event * tfdd + tfd) + tf);
         totalphs = endphs - begphs;
         begphs = begphs < 0.0 ? fmod(begphs, 1.0) + 1.0 : fmod(begphs, 1.0);
         endphs = endphs < 0.0 ? fmod(endphs, 1.0) + 1.0 : fmod(endphs, 1.0);
         dphs = 1.0 / search.proflen;
         /* printf("%.4f %.4f %5.4f\n", begphs, endphs, totalphs); */
         calctotalphs = 0.0;
         for (jj = 0; jj < search.proflen; jj++) {
            numwraps = floor(totalphs);
            lphs = jj * dphs;
            rphs = (jj + 1) * dphs;
            if (begphs <= lphs) {       /* BLR */
               if (rphs <= endphs) {    /* LRE */
                  numwraps += 1.0;
               } else if (endphs <= lphs) {     /* ELR */
                  if (endphs <= begphs)
                     numwraps += 1.0;
               } else {         /* LER */
                  numwraps += (endphs - lphs) * search.proflen;
               }
            } else if (rphs <= begphs) {        /* LRB */
               if (rphs <= endphs) {    /* LRE */
                  if (endphs <= begphs)
                     numwraps += 1.0;
               } else if (lphs <= endphs && endphs <= rphs) {   /* LER */
                  numwraps += (endphs - lphs) * search.proflen;
               }
            } else {            /* LBR */
               numwraps += (rphs - begphs) * search.proflen;    /* All E's */
               if (lphs <= endphs && endphs <= rphs) {  /* LER */
                  numwraps += (endphs - lphs) * search.proflen;
                  if (begphs <= endphs)
                     numwraps -= 1.0;
               }
            }
            /* printf("%.2f ", numwraps); */
            calctotalphs += numwraps;
            if (numwraps > 0)
               search.rawfolds[ii * search.proflen + jj] *= (totalphs / numwraps);
         }
         /* printf("\n"); */
         calctotalphs /= search.proflen;
         if (fabs(totalphs - calctotalphs) > 0.00001)
            printf("\nThere seems to be a problem in the \"exposure\" calculation\n"
                   "  in prepfold (npart = %ld):  totalphs = %.6f but calctotalphs = %0.6f\n",
                   ii, totalphs, calctotalphs);
         cts = 0.0;
         for (jj = ii * search.proflen; jj < (ii + 1) * search.proflen; jj++)
            cts += search.rawfolds[jj];
         search.stats[ii].numdata = ceil((T / cmd->npart) / search.dt);
         search.stats[ii].numprof = search.proflen;
         search.stats[ii].prof_avg = search.stats[ii].prof_var =
             cts / search.proflen;
         search.stats[ii].data_avg = search.stats[ii].data_var =
             search.stats[ii].prof_avg / search.stats[ii].numdata;
         /* Compute the Chi-Squared probability that there is a signal */
         /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */
         search.stats[ii].redchi = 0.0;
         for (jj = ii * search.proflen; jj < (ii + 1) * search.proflen; jj++) {
            dtmp = search.rawfolds[jj] - search.stats[ii].prof_avg;
            search.stats[ii].redchi += dtmp * dtmp;
         }
         search.stats[ii].redchi /= (search.stats[ii].prof_var *
                                     (search.proflen - 1));
      }
      printf("\r  Folded %d events.", numevents);
      fflush(NULL);

   } else {                     /* Fold a time series */

      buffers = gen_dvect(cmd->nsub * search.proflen);
      phasesadded = gen_dvect(cmd->nsub);
      for (ii = 0; ii < cmd->nsub * search.proflen; ii++)
         buffers[ii] = 0.0;
      for (ii = 0; ii < cmd->nsub; ii++)
         phasesadded[ii] = 0.0;

      /* Move to the correct starting record */

      data = gen_fvect(cmd->nsub * worklen);
      if (cmd->pkmbP)
         skip_to_PKMB_rec(infiles, numfiles, lorec + 1);
      else if (cmd->bcpmP)
         skip_to_BPP_rec(infiles, numfiles, lorec + 1);
      else if (cmd->spigotP)
         skip_to_SPIGOT_rec(infiles, numfiles, lorec + 1);
      else if (cmd->wappP)
         skip_to_WAPP_rec(infiles, numfiles, lorec + 1);
      else if (cmd->gmrtP)
         skip_to_GMRT_rec(infiles, numfiles, lorec + 1);
      else if (cmd->filterbankP)
         skip_to_filterbank_rec(infiles, numfiles, lorec + 1);
      else if (cmd->psrfitsP)
          skip_to_PSRFITS_samp(lorec*worklen);
      else {
         if (useshorts) {
            int reclen = 1;

            if (insubs)
               reclen = SUBSBLOCKLEN;
            /* Use a loop to accommodate subband data */
            for (ii = 0; ii < numfiles; ii++)
               chkfileseek(infiles[ii], lorec * reclen, sizeof(short), SEEK_SET);
         } else {
            chkfileseek(infiles[0], lorec, sizeof(float), SEEK_SET);
         }
      }

      /* Data is already barycentered or 
         the polycos will take care of the barycentering or
         we are folding topocentrically.  */
      if (!(RAWDATA || insubs) || cmd->polycofileP || cmd->topoP) {
         foldf = f;
         foldfd = fd;
         foldfdd = fdd;
         orig_foldf = foldf;
      } else {                  /* Correct our fold parameters if we are barycentering */
         double *voverc;

         /* The number of topo to bary points to generate with TEMPO */

         numbarypts = T / TDT + 10;
         barytimes = gen_dvect(numbarypts);
         topotimes = gen_dvect(numbarypts);
         voverc = gen_dvect(numbarypts);

         /* topocentric times in days from data start */

         for (ii = 0; ii < numbarypts; ii++)
            topotimes[ii] = search.tepoch + (double) ii *TDT / SECPERDAY;

         /* Call TEMPO for the barycentering */

         printf("\nGenerating barycentric corrections...\n");
         barycenter(topotimes, barytimes, voverc, numbarypts,
                    rastring, decstring, obs, ephem);

         /* Determine the avg v/c of the Earth's motion during the obs */

         for (ii = 0; ii < numbarypts - 1; ii++)
            search.avgvoverc += voverc[ii];
         search.avgvoverc /= (numbarypts - 1.0);
         vect_free(voverc);
         printf("The average topocentric velocity is %.6g (units of c).\n\n",
                search.avgvoverc);
         printf("Barycentric folding frequency    (hz)  =  %-.12g\n", f);
         printf("Barycentric folding f-dot      (hz/s)  =  %-.8g\n", fd);
         printf("Barycentric folding f-dotdot (hz/s^2)  =  %-.8g\n", fdd);

         /* Convert the barycentric folding parameters into topocentric */

         if ((info = bary2topo(topotimes, barytimes, numbarypts,
                               f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
            printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
         printf("Topocentric folding frequency    (hz)  =  %-.12g\n", foldf);
         if (foldfd != 0.0)
            printf("Topocentric folding f-dot      (hz/s)  =  %-.8g\n", foldfd);
         if (foldfdd != 0.0)
            printf("Topocentric folding f-dotdot (hz/s^2)  =  %-.8g\n", foldfdd);
         printf("\n");

         /* Modify the binary delay times so they refer to */
         /* topocentric reference times.                   */

         if (binary) {
            for (ii = 0; ii < numbinpoints; ii++) {
               arrayoffset++;   /* Beware nasty NR zero-offset kludges! */
               dtmp = search.bepoch + tp[ii] / SECPERDAY;
               hunt(barytimes, numbarypts, dtmp, &arrayoffset);
               arrayoffset--;
               tp[ii] = LININTERP(dtmp, barytimes[arrayoffset],
                                  barytimes[arrayoffset + 1],
                                  topotimes[arrayoffset],
                                  topotimes[arrayoffset + 1]);
            }
            numdelays = numbinpoints;
            dtmp = tp[0];
            for (ii = 0; ii < numdelays; ii++)
               tp[ii] = (tp[ii] - dtmp) * SECPERDAY;
         }
      }

      /* Record the f, fd, and fdd we used to do the raw folding */

      if (idata.bary)
         search.fold.pow = 1.0;
      search.fold.p1 = foldf;
      search.fold.p2 = foldfd;
      search.fold.p3 = foldfdd;

      /* If this is 'raw' radio data, determine the dispersion delays */

      if (!strcmp(idata.band, "Radio")) {

         /* The observation frequencies */

         obsf = gen_dvect(numchan);
         obsf[0] = idata.freq;
         search.numchan = numchan;
         search.lofreq = idata.freq;
         search.bestdm = idata.dm;
         search.chan_wid = idata.chan_wid;
         for (ii = 1; ii < numchan; ii++)
            obsf[ii] = obsf[0] + ii * idata.chan_wid;
         if (RAWDATA || insubs) {
            for (ii = 0; ii < numchan; ii++)
               obsf[ii] = doppler(obsf[ii], search.avgvoverc);
         }
         dispdts = subband_search_delays(numchan, cmd->nsub, cmd->dm,
                                         idata.freq, idata.chan_wid,
                                         search.avgvoverc);

         /* Convert the delays in seconds to delays in bins */

         for (ii = 0; ii < numchan; ii++)
            dispdts[ii] /= search.dt;

         if (cmd->nsub > 1 && (RAWDATA || insubs)) {
            int numdmtrials;
            double dphase, lodm, hidm, ddm;

            dphase = 1 / (foldf * search.proflen);
            ddm = dm_from_delay(dphase * cmd->dmstep, obsf[0]);
            numdmtrials = 2 * cmd->ndmfact * search.proflen + 1;
            lodm = cmd->dm - (numdmtrials - 1) / 2 * ddm;
            if (lodm < 0.0)
               lodm = 0.0;
            hidm = lodm + numdmtrials * ddm;
            printf("Will search %d DMs from %.3f to %.3f (ddm = %.4f)\n",
                   numdmtrials, lodm, hidm, ddm);
         }
      }

      /* 
       *   Perform the actual folding of the data
       */

      printf("\nStarting work on '%s'...\n\n", search.filenm);
      proftime = worklen * search.dt;
      parttimes = gen_dvect(cmd->npart);
      printf("  Folded %ld points of %.0f", totnumfolded, N);

      /* sub-integrations in time  */

      dtmp = (double) cmd->npart;
      for (ii = 0; ii < cmd->npart; ii++) {
         parttimes[ii] = ii * reads_per_part * proftime;

         /* reads per sub-integration */

         for (jj = 0; jj < reads_per_part; jj++) {
            double fold_time0;

            if (cmd->pkmbP)
               numread = read_PKMB_subbands(infiles, numfiles, data,
                                            dispdts, cmd->nsub, 1, &padding,
                                            maskchans, &nummasked, &obsmask);
            else if (cmd->bcpmP)
               numread = read_BPP_subbands(infiles, numfiles, data,
                                           dispdts, cmd->nsub, 1, &padding,
                                           maskchans, &nummasked, &obsmask, ifs);
            else if (cmd->spigotP)
               numread = read_SPIGOT_subbands(infiles, numfiles, data,
                                              dispdts, cmd->nsub, 1, &padding,
                                              maskchans, &nummasked, &obsmask, ifs);
            else if (cmd->wappP)
               numread = read_WAPP_subbands(infiles, numfiles, data,
                                            dispdts, cmd->nsub, 1, &padding,
                                            maskchans, &nummasked, &obsmask, ifs);
            else if (cmd->gmrtP)
               numread = read_GMRT_subbands(infiles, numfiles, data,
                                            dispdts, cmd->nsub, 1, &padding,
                                            maskchans, &nummasked, &obsmask);
            else if (cmd->filterbankP)
               numread = read_filterbank_subbands(infiles, numfiles, data,
                                                  dispdts, cmd->nsub, 1, &padding,
                                                  maskchans, &nummasked, &obsmask);
            else if (cmd->psrfitsP)
               numread = read_PSRFITS_subbands(data, dispdts, cmd->nsub, 1, &padding,
                                               maskchans, &nummasked, &obsmask);
            else if (insubs)
               numread = read_subbands(infiles, numfiles, data, recdt,
                                       maskchans, &nummasked, &obsmask, padvals);
            else {
               int mm;
               float runavg = 0.0;
               static float oldrunavg = 0.0;
               static int firsttime = 1;

               if (useshorts)
                  numread = read_shorts(infiles[0], data, worklen, numchan);
               else
                  numread = read_floats(infiles[0], data, worklen, numchan);
               if (cmd->runavgP) {
                  for (mm = 0; mm < numread; mm++)
                     runavg += data[mm];
                  runavg /= numread;
                  if (firsttime) {
                      firsttime = 0;
                  } else {
                      // Use a running average of the block averages to subtract...
                      runavg = 0.95 * oldrunavg + 0.05 * runavg;
                  }
                  oldrunavg = runavg;
                  for (mm = 0; mm < numread; mm++)
                     data[mm] -= runavg;
               }
            }

            if (cmd->polycofileP) {     /* Update the period/phase */
               double mjdf, currentsec, currentday, offsetphase, orig_cmd_phs = 0.0;

               if (ii == 0 && jj == 0)
                  orig_cmd_phs = cmd->phs;
               currentsec = parttimes[ii] + jj * proftime;
               currentday = currentsec / SECPERDAY;
               mjdf = idata.mjd_f + startTday + currentday;
               /* Calculate the pulse phase at the start of the current block */
               polyco_index =
                   phcalc(idata.mjd_i, mjdf, polyco_index, &polyco_phase, &foldf);
               if (!cmd->absphaseP)
                  polyco_phase -= polyco_phase0;
               if (polyco_phase < 0.0)
                  polyco_phase += 1.0;
               /* Calculate the folding frequency at the middle of the current block */
               polyco_index = phcalc(idata.mjd_i, mjdf + 0.5 * proftime / SECPERDAY,
                                     polyco_index, &offsetphase, &foldf);
               cmd->phs = orig_cmd_phs + polyco_phase;
               fold_time0 = 0.0;
            } else {
               fold_time0 = parttimes[ii] + jj * proftime;
            }

            /* Fold the frequency sub-bands */

            for (kk = 0; kk < cmd->nsub; kk++) {
               /* This is a quick hack to see if it will remove power drifts */
               if (cmd->runavgP && (numread > 0)) {
                  int dataptr;
                  double avg, var;
                  avg_var(data + kk * worklen, numread, &avg, &var);
                  for (dataptr = 0; dataptr < worklen; dataptr++)
                     data[kk * worklen + dataptr] -= avg;
               }
               fold(data + kk * worklen, numread, search.dt,
                    fold_time0,
                    search.rawfolds + (ii * cmd->nsub + kk) * search.proflen,
                    search.proflen, cmd->phs, buffers + kk * search.proflen,
                    phasesadded + kk, foldf, foldfd, foldfdd, flags, Ep, tp,
                    numdelays, NULL, &(search.stats[ii * cmd->nsub + kk]));
            }
            totnumfolded += numread;
         }

         printf("\r  Folded %ld points of %.0f", totnumfolded, N);
         fflush(NULL);
      }
      vect_free(buffers);
      vect_free(phasesadded);
   }
   /* This resets foldf (which is used below) to the original value */
   if (cmd->polycofileP)
      foldf = orig_foldf;
   for (ii = 0; ii < numfiles; ii++)
      fclose(infiles[ii]);
   free(infiles);

   /*
    *   Perform the candidate optimization search
    */

   printf("\n\nOptimizing...\n\n");
   bestprof = gen_dvect(search.proflen);
   {
      int numtrials, totpdelay;
      int good_idm = 0, good_ip = 0, good_ipd = 0, good_ipdd = 0;
      double dphase, po, pdo, pddo;
      double *currentprof, *fdots, *fdotdots = NULL;
      foldstats currentstats;

      search.ndmfact = cmd->ndmfact;
      search.npfact = cmd->npfact;
      search.pstep = cmd->pstep;
      search.pdstep = cmd->pdstep;
      search.dmstep = cmd->dmstep;

      /* The number of trials for the P-dot and P searches */

      numtrials = 2 * search.npfact * search.proflen + 1;

      /* If we really don't need to search for the pulsations, */
      /* Don't make us search (and display) a hugh p/pd plane  */

      if (cmd->nosearchP && search.proflen > 100)
         numtrials = 201;

      /* Initialize a bunch of variables */

      search.numperiods = numtrials;
      search.periods = gen_dvect(numtrials);
      search.pdots = gen_dvect(numtrials);
      fdots = gen_dvect(numtrials);
      fdotdots = gen_dvect(numtrials);
      if (cmd->searchfddP)
         cmd->searchpddP = 1;
      if (cmd->nopsearchP && cmd->nopdsearchP) {
         if (cmd->nsub > 1 && cmd->nodmsearchP)
            cmd->nosearchP = 1;
         else if (cmd->nsub == 1)
            cmd->nosearchP = 1;
         pflags.nosearch = cmd->nosearchP;
      }
      if (cmd->nosearchP) {
         cmd->nopsearchP = cmd->nopdsearchP = 1;
         if (cmd->nsub > 1)
            cmd->nodmsearchP = 1;
      }
      search.numpdots = numtrials;
      currentprof = gen_dvect(search.proflen);
      initialize_foldstats(&beststats);
      if (cmd->nopsearchP)
         good_ip = (numtrials - 1) / 2;
      if (cmd->nopdsearchP)
         good_ipd = (numtrials - 1) / 2;
      if (cmd->nosearchP && cmd->searchpddP)
         good_ipdd = (numtrials - 1) / 2;

      /* Convert the folding freqs and derivs to periods */

      switch_f_and_p(foldf, foldfd, foldfdd, &po, &pdo, &pddo);

      /* Our P and P-dot steps are the changes that cause the pulse      */
      /* to be delayed a number of bins between the first and last time. */

      dphase = po / search.proflen;
      for (ii = 0; ii < numtrials; ii++) {
         totpdelay = ii - (numtrials - 1) / 2;
         dtmp = (double) (totpdelay * search.pstep) / search.proflen;
         search.periods[ii] = 1.0 / (foldf + dtmp / T);
         dtmp = (double) (totpdelay * search.pdstep) / search.proflen;
         fdots[ii] = phasedelay2fdot(dtmp, T);
         search.pdots[ii] = switch_pfdot(foldf, foldfd + fdots[ii]);
         fdotdots[ii] = phasedelay2fdotdot(dtmp, T);
      }

      {                         /* Do the optimization */
         int numdmtrials = 1, numpdds = 1, lodmnum = 0;
         int totnumtrials, currtrial = 0;
         int oldper = -1, newper = 0;
         int idm, ip, ipd, ipdd, bestidm = 0, bestip = 0, bestipd = 0, bestipdd = 0;
         double lodm = 0.0, ddm = 0.0;
         double *delays, *pd_delays, *pdd_delays, *ddprofs = search.rawfolds;
         foldstats *ddstats = search.stats;

         if (cmd->searchpddP)
            numpdds = numtrials;
         delays = gen_dvect(cmd->npart);
         pd_delays = gen_dvect(cmd->npart);
         pdd_delays = gen_dvect(cmd->npart);

         if (cmd->nsub > 1) {   /* This is only for doing DM searches */
            ddprofs = gen_dvect(cmd->npart * search.proflen);
            ddstats = (foldstats *) malloc(cmd->npart * sizeof(foldstats));
            numdmtrials = 2 * search.ndmfact * search.proflen + 1;
            search.numdms = numdmtrials;
            search.dms = gen_dvect(numdmtrials);

            /* Our DM step is the change in DM that would cause the pulse   */
            /* to be delayed a number of phasebins at the lowest frequency. */

            ddm = dm_from_delay(dphase * search.dmstep, obsf[0]);

            /* Insure that we don't try a dm < 0.0 */

            lodm = cmd->dm - (numdmtrials - 1) / 2 * ddm;
            if (cmd->nodmsearchP)
               good_idm = (numdmtrials - 1) / 2;
            if (lodm < 0.0) {
               lodm = 0.0;
               /* Find the closest DM to the requested DM */
               if (cmd->nodmsearchP) {
                  double mindmerr = 1000.0, dmerr, trialdm;
                  for (idm = 0; idm < numdmtrials; idm++) {     /* Loop over DMs */
                     trialdm = lodm + idm * ddm;
                     dmerr = fabs(trialdm - cmd->dm);
                     if (dmerr < mindmerr) {
                        good_idm = idm;
                        mindmerr = dmerr;
                     }
                  }
               }
            }
            for (idm = 0; idm < numdmtrials; idm++)
               search.dms[idm] = lodm + idm * ddm;

            if (cmd->searchpddP)
               printf
                   ("  Searching %d DMs, %d periods, %d p-dots, and %d p-dotdots...\n",
                    search.numdms, search.numperiods, search.numpdots,
                    search.numpdots);
            else
               printf("  Searching %d DMs, %d periods, and %d p-dots...\n",
                      search.numdms, search.numperiods, search.numpdots);
         } else {               /* No DM searches are to be done */
            /* Allocate our p-pdot plane to speed up plotting */
            if (!cmd->searchpddP)
               ppdot = gen_fvect(search.numpdots * search.numperiods);
            if (cmd->searchpddP)
               printf("  Searching %d periods, %d p-dots, and %d p-dotdots...\n",
                      search.numperiods, search.numpdots, search.numpdots);
            else
               printf("  Searching %d periods, and %d p-dots...\n",
                      search.numperiods, search.numpdots);
         }

         /* Skip searching over all the DMs if we know that we're */
         /* only going to look at one of them.                    */
         if (cmd->nodmsearchP) {
            lodmnum = good_idm;
            numdmtrials = lodmnum + 1;
         }
         /* The total number of search trials */
         totnumtrials=numdmtrials*numpdds*numtrials*numtrials;

         for (idm = lodmnum; idm < numdmtrials; idm++) {        /* Loop over DMs */
            if (cmd->nsub > 1) {        /* This is only for doing DM searches */
               if (!cmd->nodmsearchP)
                  good_idm = idm;
               correct_subbands_for_DM(search.dms[idm], &search, ddprofs, ddstats);
            }

            for (ipdd = 0; ipdd < numpdds; ipdd++) {    /* Loop over pdds */
               if (!cmd->nosearchP)
                  good_ipdd = ipdd;
               for (ii = 0; ii < cmd->npart; ii++)
                  pdd_delays[ii] = cmd->searchpddP ?
                      fdotdot2phasedelay(fdotdots[ipdd],
                                         parttimes[ii]) * search.proflen : 0.0;

               for (ipd = 0; ipd < numtrials; ipd++) {  /* Loop over the pds */
                  if (!cmd->nopdsearchP)
                     good_ipd = ipd;
                  else if (cmd->nopdsearchP && cmd->nsub > 1 && ipd != good_ipd)
                     continue;
                  for (ii = 0; ii < cmd->npart; ii++)
                     pd_delays[ii] = pdd_delays[ii] +
                         fdot2phasedelay(fdots[ipd], parttimes[ii]) * search.proflen;

                  for (ip = 0; ip < numtrials; ip++) {  /* Loop over the ps */
                     if (!cmd->nopsearchP)
                        good_ip = ip;
                     else if (cmd->nopsearchP && cmd->nsub > 1 && ip != good_ip)
                        continue;
                     totpdelay = search.pstep * (ip - (numtrials - 1) / 2);
                     for (ii = 0; ii < cmd->npart; ii++)
                        delays[ii] =
                            pd_delays[ii] + (double) (ii * totpdelay) / cmd->npart;

                     /* Combine the profiles usingthe above computed delays */
                     combine_profs(ddprofs, ddstats, cmd->npart, search.proflen,
                                   delays, currentprof, &currentstats);
                    
                     /* If this is a simple fold, create the chi-square p-pdot plane */
                     if (cmd->nsub == 1 && !cmd->searchpddP)
                        ppdot[ipd * search.numpdots + ip] = currentstats.redchi;

                     /* If this is the best profile or if it is the specific */
                     /* profile that we are looking for, save it.            */
                     if (idm == good_idm && ipdd == good_ipdd && ipd == good_ipd
                         && ip == good_ip) {
                        if (cmd->nosearchP
                            || (currentstats.redchi > beststats.redchi
                                && !cmd->nosearchP)) {
                           bestidm = idm;
                           bestipdd = ipdd;
                           bestipd = ipd;
                           bestip = ip;
                           beststats = currentstats;
                           memcpy(bestprof, currentprof,
                                  sizeof(double) * search.proflen);
                        }
                     }
                     currtrial+=1;
                     newper = (int) ((float) currtrial/totnumtrials*100.0 + 0.5);
                     if (newper > oldper) {
                        printf("\r  Amount Complete = %3d%%", newper);
                        fflush(stdout);
                        oldper = newper;
                     }
                  }
               }
            }
         }

         /* Convert the indices of the folds into p, pd, pdd and DM values */
         if (cmd->nsub > 1)
            search.bestdm = search.dms[bestidm];
         if (idata.bary) {
            search.bary.p1 = search.periods[bestip];
            search.bary.p2 = search.pdots[bestipd];
            if (cmd->searchpddP)
               search.bary.p3 = switch_pfdotdot(1.0 / search.periods[bestip],
                                                foldfd + fdots[bestipd],
                                                foldfdd + fdotdots[bestipdd]);
         } else {
            search.topo.p1 = search.periods[bestip];
            search.topo.p2 = search.pdots[bestipd];
            if (cmd->searchpddP)
               search.topo.p3 = switch_pfdotdot(1.0 / search.periods[bestip],
                                                foldfd + fdots[bestipd],
                                                foldfdd + fdotdots[bestipdd]);
         }

         vect_free(delays);
         vect_free(pd_delays);
         vect_free(pdd_delays);
         if (cmd->nsub > 1) {
            vect_free(ddprofs);
            free(ddstats);
         }
      }
      vect_free(currentprof);
      vect_free(fdots);
      vect_free(fdotdots);
   }
   printf("  Done searching.\n\n");

   /* Write and plot the results and cleanup */

   {
      double perr, pderr, pdderr;
      char out[100];

      printf("Maximum reduced chi-squared found  =  %-.5f\n", beststats.redchi);
      if (cmd->nsub > 1)
         printf("Best DM     (pc cm^-3)  =  %-.4f\n", search.bestdm);

      /* Convert best params from/to barycentric to/from topocentric */

      if (!(RAWDATA || insubs) || cmd->polycofileP || cmd->topoP) {

         /* Data was barycentered */

         if (idata.bary) {

            /* Calculate the errors in our new pulsation quantities */

            fold_errors(bestprof, search.proflen, search.dt, N,
                        beststats.data_var, search.bary.p1, search.bary.p2,
                        search.bary.p3, &perr, &pderr, &pdderr);

            nice_output_2(out, search.bary.p1, perr, 0);
            printf("Best period        (s)  =  %s\n", out);
            nice_output_2(out, search.bary.p2, pderr, 0);
            printf("Best p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.bary.p3, pdderr, 0);
            printf("Best p-dotdot  (s/s^2)  =  %s\n", out);

            /* Data was topocentric */

         } else {

            /* Calculate the errors in our new pulsation quantities */

            fold_errors(bestprof, search.proflen, search.dt, N,
                        beststats.data_var, search.topo.p1, search.topo.p2,
                        search.topo.p3, &perr, &pderr, &pdderr);

            nice_output_2(out, search.topo.p1, perr, 0);
            printf("Best period        (s)  =  %s\n", out);
            nice_output_2(out, search.topo.p2, pderr, 0);
            printf("Best p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.topo.p3, pdderr, 0);
            printf("Best p-dotdot  (s/s^2)  =  %s\n", out);

         }

      } else {

         if (idata.bary) {

            /* Calculate the errors in our new pulsation quantities */

            fold_errors(bestprof, search.proflen, search.dt, N,
                        beststats.data_var, search.bary.p1, search.bary.p2,
                        search.bary.p3, &perr, &pderr, &pdderr);

            nice_output_2(out, search.bary.p1, perr, 0);
            printf("Best barycentric period        (s)  =  %s\n", out);
            nice_output_2(out, search.bary.p2, pderr, 0);
            printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.bary.p3, pdderr, 0);
            printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);

            /* Convert the barycentric folding parameters into topocentric */

            switch_f_and_p(search.bary.p1, search.bary.p2, search.bary.p3,
                           &f, &fd, &fdd);
            if ((info = bary2topo(topotimes, barytimes, numbarypts,
                                  f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
               printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
            switch_f_and_p(foldf, foldfd, foldfdd,
                           &search.topo.p1, &search.topo.p2, &search.topo.p3);

            nice_output_2(out, search.topo.p1, perr, 0);
            printf("Best topocentric period        (s)  =  %s\n", out);
            nice_output_2(out, search.topo.p2, pderr, 0);
            printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.topo.p3, pdderr, 0);
            printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

         } else {

            /* Calculate the errors in our new pulsation quantities */

            fold_errors(bestprof, search.proflen, search.dt, N,
                        beststats.data_var, search.topo.p1, search.topo.p2,
                        search.topo.p3, &perr, &pderr, &pdderr);

            nice_output_2(out, search.topo.p1, perr, 0);
            printf("Best topocentric period        (s)  =  %s\n", out);
            nice_output_2(out, search.topo.p2, pderr, 0);
            printf("Best topocentric p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.topo.p3, pdderr, 0);
            printf("Best topocentric p-dotdot  (s/s^2)  =  %s\n", out);

            /* Convert the barycentric folding parameters into topocentric */

            switch_f_and_p(search.topo.p1, search.topo.p2, search.topo.p3,
                           &f, &fd, &fdd);
            if ((info = bary2topo(barytimes, topotimes, numbarypts,
                                  f, fd, fdd, &foldf, &foldfd, &foldfdd)) < 0)
               printf("\nError in bary2topo().  Argument %d was bad.\n\n", -info);
            switch_f_and_p(foldf, foldfd, foldfdd,
                           &search.bary.p1, &search.bary.p2, &search.bary.p3);

            nice_output_2(out, search.bary.p1, perr, 0);
            printf("Best barycentric period        (s)  =  %s\n", out);
            nice_output_2(out, search.bary.p2, pderr, 0);
            printf("Best barycentric p-dot       (s/s)  =  %s\n", out);
            nice_output_2(out, search.bary.p3, pdderr, 0);
            printf("Best barycentric p-dotdot  (s/s^2)  =  %s\n", out);
         }
      }
   }

   printf("\nMaking plots.\n\n");

   /*
    *   Write the raw prepfoldinfo structure
    */

   write_prepfoldinfo(&search, outfilenm);

   /*
    *   Plot our results
    */

   prepfold_plot(&search, &pflags, !cmd->noxwinP, ppdot);

   /* Free our memory  */

   if (cmd->nsub == 1 && !cmd->searchpddP)
      vect_free(ppdot);
   delete_prepfoldinfo(&search);
   vect_free(data);
   if (!cmd->outfileP)
      free(rootnm);
   free(outfilenm);
   free(plotfilenm);
   vect_free(parttimes);
   vect_free(bestprof);
   if (binary) {
      vect_free(Ep);
      vect_free(tp);
   }
   if (cmd->maskfileP) {
      free_mask(obsmask);
      vect_free(padvals);
   }
   if (RAWDATA || insubs) {
      vect_free(barytimes);
      vect_free(topotimes);
   }
   if (!strcmp(idata.band, "Radio")) {
      vect_free(obsf);
      vect_free(dispdts);
   }
   printf("Done.\n\n");
   return (0);
}
