#include <limits.h>
#include <ctype.h>
#include "presto.h"
#include "rfifind_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "gmrt.h"
#include "spigot.h"
#include "sigproc_fb.h"
#include "psrfits.h"
#include "rfifind.h"

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP || cmd->gmrtP || cmd->spigotP || cmd->filterbankP || cmd->psrfitsP)

/* Some function definitions */

void rfifind_plot(int numchan, int numint, int ptsperint,
                  float timesigma, float freqsigma,
                  float inttrigfrac, float chantrigfrac,
                  float **dataavg, float **datastd, float **datapow,
                  int *userchan, int numuserchan,
                  int *userints, int numuserints,
                  infodata * idata, unsigned char **bytemask,
                  mask * oldmask, mask * newmask,
                  rfi * rfivect, int numrfi, int rfixwin, int rfips, int xwin);
static void write_rfifile(char *rfifilenm, rfi * rfivect, int numrfi,
                          int numchan, int numint, int ptsperint,
                          int lobin, int numbetween, int harmsum,
                          float fracterror, float freqsigma);
static void write_statsfile(char *statsfilenm, float *datapow,
                            float *dataavg, float *datastd,
                            int numchan, int numint, int ptsperint,
                            int lobin, int numbetween);
static void read_rfifile(char *rfifilenm, rfi ** rfivect, int *numrfi,
                         int *numchan, int *numint, int *ptsperint,
                         int *lobin, int *numbetween, int *harmsum,
                         float *fracterror, float *freqsigma);
static void read_statsfile(char *statsfilenm, float ***datapow,
                           float ***dataavg, float ***datastd,
                           int *numchan, int *numint, int *ptsperint,
                           int *lobin, int *numbetween);
static int *get_prime_factors(int nn, int *numfactors);
int compare_rfi_sigma(const void *ca, const void *cb);
int compare_rfi_numobs(const void *ca, const void *cb);
int read_subband_rawblocks(FILE * infiles[], int numfiles, short *subbanddata,
                           int numsamples, int *padding);
void get_subband(int subbandnum, float chandat[], short srawdata[], int numsamples);
extern int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals);
extern void read_offsets(float **loptr, float **optr, int numpts, int numchan);

/* The main program */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   FILE **infiles = NULL, *bytemaskfile;
   float **dataavg = NULL, **datastd = NULL, **datapow = NULL, *padvals;
   float *chandata = NULL, powavg, powstd, powmax;
   float inttime, norm, fracterror = RFI_FRACTERROR;
   float *offsets, *nextoffsets;
   unsigned char *rawdata = NULL, **bytemask = NULL;
   short *srawdata = NULL;
   char *outfilenm, *statsfilenm, *maskfilenm;
   char *bytemaskfilenm, *rfifilenm;
   int numchan = 0, numint = 0, newper = 0, oldper = 0, numfiles, good_padvals = 0;
   int blocksperint, ptsperint = 0, ptsperblock = 0, padding = 0;
   int numcands, candnum, numrfi = 0, numrfivect = NUM_RFI_VECT;
   int ii, jj, kk, slen, numread = 0, insubs = 0;
   int harmsum = RFI_NUMHARMSUM, lobin = RFI_LOBIN, numbetween = RFI_NUMBETWEEN;
   double davg, dvar, freq, dt, T;
   long long N;
   IFs ifs = SUMIFS;
   presto_interptype interptype;
   rfi *rfivect = NULL;
   mask oldmask, newmask;
   fftcand *cands;
   PKMB_tapehdr hdr;
   infodata idata;
   Cmdline *cmd;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      printf("\n");
      usage();
      exit(1);
   }

   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);
   slen = strlen(cmd->outfile) + 20;
   numfiles = cmd->argc;
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
   printf("               Pulsar Data RFI Finder\n");
   printf("                 by Scott M. Ransom\n\n");

   /* The following is the root of all the output files */

   outfilenm = (char *) calloc(slen, sizeof(char));
   sprintf(outfilenm, "%s_rfifind", cmd->outfile);

   /* And here are the output file names */

   maskfilenm = (char *) calloc(slen, sizeof(char));
   sprintf(maskfilenm, "%s.mask", outfilenm);
   bytemaskfilenm = (char *) calloc(slen, sizeof(char));
   sprintf(bytemaskfilenm, "%s.bytemask", outfilenm);
   rfifilenm = (char *) calloc(slen, sizeof(char));
   sprintf(rfifilenm, "%s.rfi", outfilenm);
   statsfilenm = (char *) calloc(slen, sizeof(char));
   sprintf(statsfilenm, "%s.stats", outfilenm);
   sprintf(idata.name, "%s", outfilenm);

   /* Read an input mask if wanted */

   if (cmd->maskfileP) {
      read_mask(cmd->maskfile, &oldmask);
      printf("Read old mask information from '%s'\n\n", cmd->maskfile);
      good_padvals = determine_padvals(cmd->maskfile, &oldmask, &padvals);
   } else {
      oldmask.numchan = oldmask.numint = 0;
   }

   if (!cmd->nocomputeP) {
      char *root, *suffix;
      if (!RAWDATA) {
         /* Split the filename into a rootname and a suffix */
         if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
            printf("\nThe input filename (%s) must have a suffix, or you must\n"
                   "    specify the data type!\n\n", cmd->argv[0]);
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
            } else if (isdigit(suffix[0]) &&
                       isdigit(suffix[1]) && isdigit(suffix[2])) {
               printf("Assuming the data is from the Arecibo WAPP system...\n");
               cmd->wappP = 1;
            } else if (strcmp(suffix, "sub0") == 0 ||
                       strcmp(suffix, "sub00") == 0 ||
                       strcmp(suffix, "sub000") == 0 ||
                       strcmp(suffix, "sub0000") == 0) {
               printf
                   ("Assuming the data is subband data generated by prepsubband...\n");
               insubs = 1;
            }
         }
         if (RAWDATA) {         /* Clean-up a bit */
            free(root);
            free(suffix);
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

      infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
      for (ii = 0; ii < numfiles; ii++) {
         printf("  '%s'\n", cmd->argv[ii]);
         infiles[ii] = chkfopen(cmd->argv[ii], "rb");
      }
      printf("\n");


      if (insubs) {

         /* Set-up values if we are using subbands */

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
         free(root);
         free(suffix);
         ptsperblock = 1;
         /* Compensate for the fact that we have subbands and not channels */
         idata.freq = idata.freq - 0.5 * idata.chan_wid +
             0.5 * idata.chan_wid * (idata.num_chan / numfiles);
         idata.chan_wid = idata.num_chan / numfiles * idata.chan_wid;
         idata.num_chan = numchan = numfiles;
         idata.dm = 0.0;
         sprintf(idata.name, "%s", outfilenm);
         writeinf(&idata);

      } else if (cmd->pkmbP) {

         /* Set-up values if we are using the Parkes multibeam */

         printf("Filterbank input file information:\n");
         get_PKMB_file_info(infiles, numfiles, cmd->clip, &N, 
                            &ptsperblock, &numchan, &dt, &T, 1);

         /* Read the first header file and generate an infofile from it */

         chkfread(&hdr, 1, HDRLEN, infiles[0]);
         rewind(infiles[0]);
         PKMB_hdr_to_inf(&hdr, &idata);
         PKMB_update_infodata(numfiles, &idata);
         idata.dm = 0.0;
         writeinf(&idata);

      } else if (cmd->bcpmP) {

         /* Set-up for the BCPM machines at Green Bank  */

         printf("BCPM input file information:\n");
         get_BPP_file_info(infiles, numfiles, cmd->clip, &N, &ptsperblock,
                           &numchan, &dt, &T, &idata, 1);
         BPP_update_infodata(numfiles, &idata);
         set_BPP_padvals(padvals, good_padvals);
         idata.dm = 0.0;
         writeinf(&idata);

      } else if (cmd->spigotP) {

         /* Set-up values if we are using the NRAO-Caltech Spigot card */

         SPIGOT_INFO *spigots;
         printf("Spigot card input file information:\n");
         spigots = (SPIGOT_INFO *) malloc(sizeof(SPIGOT_INFO) * numfiles);
         for (ii = 0; ii < numfiles; ii++)
            read_SPIGOT_header(cmd->argv[ii], spigots + ii);
         get_SPIGOT_file_info(infiles, spigots, numfiles, cmd->windowP, cmd->clip,
                              &N, &ptsperblock, &numchan, &dt, &T, &idata, 1);
         SPIGOT_update_infodata(numfiles, &idata);
         set_SPIGOT_padvals(padvals, good_padvals);
         idata.dm = 0.0;
         writeinf(&idata);
         free(spigots);

      } else if (cmd->psrfitsP) {

         /* Set-up values if we are using search-mode PSRFITS data */

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
         writeinf(&idata);

      } else if (cmd->wappP) {

         /* Set-up for the WAPP machine at Arecibo */

         printf("WAPP input file information:\n");
         get_WAPP_file_info(infiles, cmd->numwapps, numfiles,
                            cmd->windowP, cmd->clip,
                            &N, &ptsperblock, &numchan, &dt, &T, &idata, 1);
         WAPP_update_infodata(numfiles, &idata);
         set_WAPP_padvals(padvals, good_padvals);
         idata.dm = 0.0;
         writeinf(&idata);

      } else if (cmd->gmrtP) {

         /* Set-up for the GMRT Phase array data */

         printf("GMRT input file information:\n");
         get_GMRT_file_info(infiles, argv + 1, numfiles, cmd->clip,
                            &N, &ptsperblock, &numchan, &dt, &T, 1);
         /* Read the first header file and generate an infofile from it */
         GMRT_hdr_to_inf(argv[1], &idata);
         GMRT_update_infodata(numfiles, &idata);
         set_GMRT_padvals(padvals, good_padvals);
         idata.dm = 0.0;
         writeinf(&idata);

      } else if (cmd->filterbankP) {
         int headerlen;
         sigprocfb fb;

         /* Set-up for SIGPROC filterbank-style data */

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
         idata.dm = 0.0;
         writeinf(&idata);
      }

      /* The number of data points and blocks to work with at a time */

      if (cmd->blocksP) {
         blocksperint = cmd->blocks;
         cmd->time = blocksperint * ptsperblock * idata.dt;
      } else {
         blocksperint = (int) (cmd->time / (ptsperblock * idata.dt) + 0.5);
      }
      ptsperint = blocksperint * ptsperblock;
      numint = (long long) idata.N / ptsperint;
      if ((long long) idata.N % ptsperint)
         numint++;
      inttime = ptsperint * idata.dt;
      printf("Analyzing data sections of length %d points (%.6g sec).\n",
             ptsperint, inttime);
      {
         int *factors, numfactors;

         factors = get_prime_factors(ptsperint, &numfactors);
         printf("  Prime factors are:  ");
         for (ii = 0; ii < numfactors; ii++)
            printf("%d ", factors[ii]);
         printf("\n");
         if (factors[numfactors - 1] > 13) {
            printf("  WARNING:  The largest prime factor is pretty big!  This will\n"
                   "            cause the FFTs to take a long time to compute.  I\n"
                   "            recommend choosing a different -time value.\n");
         }
         printf("\n");
         free(factors);
      }

      /* Allocate our workarrays */

      if (cmd->pkmbP)
         rawdata = gen_bvect(DATLEN * blocksperint);
      else if (cmd->bcpmP || cmd->spigotP || cmd->wappP || cmd->gmrtP
               || cmd->filterbankP || cmd->psrfitsP)
         /* This allocates extra incase both IFs were stored */
         rawdata = gen_bvect(idata.num_chan * ptsperblock * blocksperint);
      else if (insubs)
         srawdata = gen_svect(idata.num_chan * ptsperblock * blocksperint);
      dataavg = gen_fmatrix(numint, numchan);
      datastd = gen_fmatrix(numint, numchan);
      datapow = gen_fmatrix(numint, numchan);
      chandata = gen_fvect(ptsperint);
      bytemask = gen_bmatrix(numint, numchan);
      for (ii = 0; ii < numint; ii++)
         for (jj = 0; jj < numchan; jj++)
            bytemask[ii][jj] = GOODDATA;
      rfivect = rfi_vector(rfivect, numchan, numint, 0, numrfivect);
      if (numbetween == 2)
         interptype = INTERBIN;
      else
         interptype = INTERPOLATE;

      /* Main loop */

      printf("Writing mask data  to '%s'.\n", maskfilenm);
      printf("Writing  RFI data  to '%s'.\n", rfifilenm);
      printf("Writing statistics to '%s'.\n\n", statsfilenm);
      printf("Massaging the data ...\n\n");
      printf("Amount Complete = %3d%%", oldper);
      fflush(stdout);

      /* Prep the offset file if required */
      read_offsets(&offsets, &nextoffsets, ptsperint, numchan);

      for (ii = 0; ii < numint; ii++) { /* Loop over the intervals */
         newper = (int) ((float) ii / numint * 100.0 + 0.5);
         if (newper > oldper) {
            printf("\rAmount Complete = %3d%%", newper);
            fflush(stdout);
            oldper = newper;
         }

         /* Read a chunk of data */

         if (cmd->pkmbP)
            numread = read_PKMB_rawblocks(infiles, numfiles,
                                          rawdata, blocksperint, &padding);
         else if (cmd->bcpmP)
            numread = read_BPP_rawblocks(infiles, numfiles,
                                         rawdata, blocksperint, &padding);
         else if (cmd->spigotP)
            numread = read_SPIGOT_rawblocks(infiles, numfiles,
                                            rawdata, blocksperint, &padding, ifs);
         else if (cmd->wappP)
            numread = read_WAPP_rawblocks(infiles, numfiles,
                                          rawdata, blocksperint, &padding, ifs);
         else if (cmd->psrfitsP)
             numread = read_PSRFITS_rawblocks(rawdata, blocksperint, &padding);
         else if (cmd->gmrtP)
            numread = read_GMRT_rawblocks(infiles, numfiles,
                                          rawdata, blocksperint, &padding);
         else if (cmd->filterbankP)
            numread = read_filterbank_rawblocks(infiles, numfiles,
                                                rawdata, blocksperint, &padding);
         else if (insubs)
            numread = read_subband_rawblocks(infiles, numfiles,
                                             srawdata, blocksperint, &padding);

         if (padding)
            for (jj = 0; jj < numchan; jj++)
               bytemask[ii][jj] |= PADDING;

         /* Read the offset file if required */
         read_offsets(&offsets, &nextoffsets, ptsperint, numchan);

         for (jj = 0; jj < numchan; jj++) {     /* Loop over the channels */

            if (cmd->pkmbP)
               get_PKMB_channel(jj, chandata, rawdata, blocksperint);
            else if (cmd->bcpmP)
               get_BPP_channel(jj, chandata, rawdata, blocksperint, ifs);
            else if (cmd->spigotP)
               get_SPIGOT_channel(jj, chandata, rawdata, blocksperint);
            else if (cmd->wappP)
               get_WAPP_channel(jj, chandata, rawdata, blocksperint);
            else if (cmd->psrfitsP)
               get_PSRFITS_channel(jj, chandata, rawdata, blocksperint);
            else if (cmd->gmrtP)
               get_GMRT_channel(jj, chandata, rawdata, blocksperint);
            else if (cmd->filterbankP)
               get_filterbank_channel(jj, chandata, rawdata, blocksperint);
            else if (insubs)
               get_subband(jj, chandata, srawdata, blocksperint);

            /* Adjust the channels based on the offsets */
            for (kk = 0; kk < ptsperint; kk++) {
                chandata[kk] -= offsets[kk];
            }

            /* Calculate the averages and standard deviations */
            /* for each point in time.                        */

            if (padding) {
               if (cmd->pkmbP) {
                  dataavg[ii][jj] = 0.5;
                  datastd[ii][jj] = 0.5;
                  datapow[ii][jj] = 1.0;
               } else {
                  dataavg[ii][jj] = 0.0;
                  datastd[ii][jj] = 0.0;
                  datapow[ii][jj] = 1.0;
               }
            } else {
               avg_var(chandata, ptsperint, &davg, &dvar);
               dataavg[ii][jj] = davg;
               datastd[ii][jj] = sqrt(dvar);
               realfft(chandata, ptsperint, -1);
               numcands = 0;
               norm = datastd[ii][jj] * datastd[ii][jj] * ptsperint;
               if (norm == 0.0)
                  norm = (chandata[0] == 0.0) ? 1.0 : chandata[0];
               cands = search_fft((fcomplex *) chandata, ptsperint / 2,
                                  lobin, ptsperint / 2, harmsum,
                                  numbetween, interptype, norm, cmd->freqsigma,
                                  &numcands, &powavg, &powstd, &powmax);
               datapow[ii][jj] = powmax;

               /* Record the birdies */

               if (numcands) {
                  for (kk = 0; kk < numcands; kk++) {
                     freq = cands[kk].r / inttime;
                     candnum = find_rfi(rfivect, numrfi, freq, RFI_FRACTERROR);
                     if (candnum >= 0) {
                        update_rfi(rfivect + candnum, freq, cands[kk].sig, jj, ii);
                     } else {
                        update_rfi(rfivect + numrfi, freq, cands[kk].sig, jj, ii);
                        numrfi++;
                        if (numrfi == numrfivect) {
                           numrfivect *= 2;
                           rfivect = rfi_vector(rfivect, numchan, numint,
                                                numrfivect / 2, numrfivect);
                        }
                     }
                  }
                  free(cands);
               }
            }
         }
      }
      printf("\rAmount Complete = 100%%\n");

      /* Write the data to the output files */

      write_rfifile(rfifilenm, rfivect, numrfi, numchan, numint,
                    ptsperint, lobin, numbetween, harmsum,
                    fracterror, cmd->freqsigma);
      write_statsfile(statsfilenm, datapow[0], dataavg[0], datastd[0],
                      numchan, numint, ptsperint, lobin, numbetween);

   } else {                     /* If "-nocompute" */
      float freqsigma;

      /* Read the data from the output files */

      printf("Reading  RFI data  from '%s'.\n", rfifilenm);
      printf("Reading statistics from '%s'.\n", statsfilenm);
      readinf(&idata, outfilenm);
      read_rfifile(rfifilenm, &rfivect, &numrfi, &numchan, &numint,
                   &ptsperint, &lobin, &numbetween, &harmsum,
                   &fracterror, &freqsigma);
      numrfivect = numrfi;
      read_statsfile(statsfilenm, &datapow, &dataavg, &datastd,
                     &numchan, &numint, &ptsperint, &lobin, &numbetween);
      bytemask = gen_bmatrix(numint, numchan);
      printf("Reading  bytemask  from '%s'.\n\n", bytemaskfilenm);
      bytemaskfile = chkfopen(bytemaskfilenm, "rb");
      chkfread(bytemask[0], numint * numchan, 1, bytemaskfile);
      fclose(bytemaskfile);
      for (ii = 0; ii < numint; ii++)
         for (jj = 0; jj < numchan; jj++)
            bytemask[ii][jj] &= PADDING;        /* Clear all but the PADDING bits */
      inttime = ptsperint * idata.dt;
   }

   /* Make the plots and set the mask */

   {
      int *zapints, *zapchan;
      int numzapints = 0, numzapchan = 0;

      if (cmd->zapintsstrP) {
         zapints = ranges_to_ivect(cmd->zapintsstr, 0, numint - 1, &numzapints);
         zapints = (int *) realloc(zapints, (size_t) (sizeof(int) * numint));
      } else {
         zapints = gen_ivect(numint);
      }
      if (cmd->zapchanstrP) {
         zapchan = ranges_to_ivect(cmd->zapchanstr, 0, numchan - 1, &numzapchan);
         zapchan = (int *) realloc(zapchan, (size_t) (sizeof(int) * numchan));
      } else {
         zapchan = gen_ivect(numchan);
      }
      rfifind_plot(numchan, numint, ptsperint, cmd->timesigma, cmd->freqsigma,
                   cmd->inttrigfrac, cmd->chantrigfrac,
                   dataavg, datastd, datapow, zapchan, numzapchan,
                   zapints, numzapints, &idata, bytemask,
                   &oldmask, &newmask, rfivect, numrfi,
                   cmd->rfixwinP, cmd->rfipsP, cmd->xwinP);

      vect_free(zapints);
      vect_free(zapchan);
   }

   /* Write the new mask and bytemask to the file */

   write_mask(maskfilenm, &newmask);
   bytemaskfile = chkfopen(bytemaskfilenm, "wb");
   chkfwrite(bytemask[0], numint * numchan, 1, bytemaskfile);
   fclose(bytemaskfile);

   /* Determine the percent of good and bad data */

   {
      int numpad = 0, numbad = 0, numgood = 0;

      for (ii = 0; ii < numint; ii++) {
         for (jj = 0; jj < numchan; jj++) {
            if (bytemask[ii][jj] == GOODDATA) {
               numgood++;
            } else {
               if (bytemask[ii][jj] & PADDING)
                  numpad++;
               else
                  numbad++;
            }
         }
      }
      printf("\nTotal number of intervals in the data:  %d\n\n", numint * numchan);
      printf("  Number of padded intervals:  %7d  (%6.3f%%)\n",
             numpad, (float) numpad / (float) (numint * numchan) * 100.0);
      printf("  Number of  good  intervals:  %7d  (%6.3f%%)\n",
             numgood, (float) numgood / (float) (numint * numchan) * 100.0);
      printf("  Number of  bad   intervals:  %7d  (%6.3f%%)\n\n",
             numbad, (float) numbad / (float) (numint * numchan) * 100.0);
      qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_sigma);
      printf("  Ten most significant birdies:\n");
      printf("#  Sigma     Period(ms)      Freq(Hz)       Number \n");
      printf("----------------------------------------------------\n");
      for (ii = 0; ii < 10; ii++) {
         double pperr;
         char temp1[40], temp2[40];

         if (rfivect[ii].freq_var == 0.0) {
            pperr = 0.0;
            sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
            sprintf(temp2, " %-14g", 1000.0 / rfivect[ii].freq_avg);
         } else {
            pperr = 1000.0 * sqrt(rfivect[ii].freq_var) /
                (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
            nice_output_2(temp1, rfivect[ii].freq_avg, sqrt(rfivect[ii].freq_var),
                          -15);
            nice_output_2(temp2, 1000.0 / rfivect[ii].freq_avg, pperr, -15);
         }
         printf("%-2d %-8.2f %13s %13s %-8d\n", ii + 1, rfivect[ii].sigma_avg,
                temp2, temp1, rfivect[ii].numobs);
      }
      qsort(rfivect, numrfi, sizeof(rfi), compare_rfi_numobs);
      printf("\n  Ten most numerous birdies:\n");
      printf("#  Number    Period(ms)      Freq(Hz)       Sigma \n");
      printf("----------------------------------------------------\n");
      for (ii = 0; ii < 10; ii++) {
         double pperr;
         char temp1[40], temp2[40];

         if (rfivect[ii].freq_var == 0.0) {
            pperr = 0.0;
            sprintf(temp1, " %-14g", rfivect[ii].freq_avg);
            sprintf(temp2, " %-14g", 1000.0 / rfivect[ii].freq_avg);
         } else {
            pperr = 1000.0 * sqrt(rfivect[ii].freq_var) /
                (rfivect[ii].freq_avg * rfivect[ii].freq_avg);
            nice_output_2(temp1, rfivect[ii].freq_avg, sqrt(rfivect[ii].freq_var),
                          -15);
            nice_output_2(temp2, 1000.0 / rfivect[ii].freq_avg, pperr, -15);
         }
         printf("%-2d %-8d %13s %13s %-8.2f\n", ii + 1, rfivect[ii].numobs,
                temp2, temp1, rfivect[ii].sigma_avg);
      }
      printf("\nDone.\n\n");
   }

   /* Close the files and cleanup */

   free_rfi_vector(rfivect, numrfivect);
   free_mask(newmask);
   if (cmd->maskfileP)
      free_mask(oldmask);
   free(outfilenm);
   free(statsfilenm);
   free(bytemaskfilenm);
   free(maskfilenm);
   free(rfifilenm);
   vect_free(dataavg[0]);
   vect_free(dataavg);
   vect_free(datastd[0]);
   vect_free(datastd);
   vect_free(datapow[0]);
   vect_free(datapow);
   vect_free(bytemask[0]);
   vect_free(bytemask);
   if (!cmd->nocomputeP) {
      for (ii = 0; ii < numfiles; ii++)
         fclose(infiles[ii]);
      free(infiles);
      vect_free(chandata);
      if (insubs)
         vect_free(srawdata);
      else
         vect_free(rawdata);
   }
   return (0);
}

static void write_rfifile(char *rfifilenm, rfi * rfivect, int numrfi,
                          int numchan, int numint, int ptsperint,
                          int lobin, int numbetween, int harmsum,
                          float fracterror, float freqsigma)
{
   FILE *outfile;
   int ii;

   outfile = chkfopen(rfifilenm, "wb");
   chkfwrite(&numchan, sizeof(int), 1, outfile);
   chkfwrite(&numint, sizeof(int), 1, outfile);
   chkfwrite(&ptsperint, sizeof(int), 1, outfile);
   chkfwrite(&lobin, sizeof(int), 1, outfile);
   chkfwrite(&numbetween, sizeof(int), 1, outfile);
   chkfwrite(&harmsum, sizeof(int), 1, outfile);
   chkfwrite(&numrfi, sizeof(int), 1, outfile);
   chkfwrite(&fracterror, sizeof(float), 1, outfile);
   chkfwrite(&freqsigma, sizeof(float), 1, outfile);
   for (ii = 0; ii < numrfi; ii++)
      write_rfi(outfile, rfivect + ii, numchan, numint);
   fclose(outfile);
}

static void write_statsfile(char *statsfilenm, float *datapow,
                            float *dataavg, float *datastd,
                            int numchan, int numint, int ptsperint,
                            int lobin, int numbetween)
{
   FILE *outfile;

   outfile = chkfopen(statsfilenm, "wb");
   chkfwrite(&numchan, sizeof(int), 1, outfile);
   chkfwrite(&numint, sizeof(int), 1, outfile);
   chkfwrite(&ptsperint, sizeof(int), 1, outfile);
   chkfwrite(&lobin, sizeof(int), 1, outfile);
   chkfwrite(&numbetween, sizeof(int), 1, outfile);
   chkfwrite(datapow, sizeof(float), numchan * numint, outfile);
   chkfwrite(dataavg, sizeof(float), numchan * numint, outfile);
   chkfwrite(datastd, sizeof(float), numchan * numint, outfile);
   fclose(outfile);
}

static void read_rfifile(char *rfifilenm, rfi ** rfivect, int *numrfi,
                         int *numchan, int *numint, int *ptsperint,
                         int *lobin, int *numbetween, int *harmsum,
                         float *fracterror, float *freqsigma)
{
   FILE *outfile;
   int ii;

   outfile = chkfopen(rfifilenm, "rb");
   chkfread(numchan, sizeof(int), 1, outfile);
   chkfread(numint, sizeof(int), 1, outfile);
   chkfread(ptsperint, sizeof(int), 1, outfile);
   chkfread(lobin, sizeof(int), 1, outfile);
   chkfread(numbetween, sizeof(int), 1, outfile);
   chkfread(harmsum, sizeof(int), 1, outfile);
   chkfread(numrfi, sizeof(int), 1, outfile);
   chkfread(fracterror, sizeof(float), 1, outfile);
   chkfread(freqsigma, sizeof(float), 1, outfile);
   *rfivect = rfi_vector(*rfivect, *numchan, *numint, 0, *numrfi);
   for (ii = 0; ii < *numrfi; ii++)
      read_rfi(outfile, *rfivect + ii, *numchan, *numint);
   fclose(outfile);
}

static void read_statsfile(char *statsfilenm, float ***datapow,
                           float ***dataavg, float ***datastd,
                           int *numchan, int *numint, int *ptsperint,
                           int *lobin, int *numbetween)
{
   FILE *outfile;

   outfile = chkfopen(statsfilenm, "rb");
   chkfread(numchan, sizeof(int), 1, outfile);
   chkfread(numint, sizeof(int), 1, outfile);
   chkfread(ptsperint, sizeof(int), 1, outfile);
   chkfread(lobin, sizeof(int), 1, outfile);
   chkfread(numbetween, sizeof(int), 1, outfile);
   *dataavg = gen_fmatrix(*numint, *numchan);
   *datastd = gen_fmatrix(*numint, *numchan);
   *datapow = gen_fmatrix(*numint, *numchan);
   chkfread(*(datapow[0]), sizeof(float), *numchan * *numint, outfile);
   chkfread(*(dataavg[0]), sizeof(float), *numchan * *numint, outfile);
   chkfread(*(datastd[0]), sizeof(float), *numchan * *numint, outfile);
   fclose(outfile);
}

int read_subband_rawblocks(FILE * infiles[], int numfiles, short *subbanddata,
                           int numsamples, int *padding)
{
   int ii, jj, index, numread = 0;

   for (ii = 0; ii < numfiles; ii++) {
      index = ii * numsamples;
      numread =
          chkfread(subbanddata + index, sizeof(short), numsamples, infiles[ii]);
      for (jj = numread; jj < numsamples; jj++)
         subbanddata[index + jj] = 0.0;
   }
   /* Subband data cannot currently be padded */
   *padding = 0;
   return numread;
}

void get_subband(int subbandnum, float chandat[], short srawdata[], int numsamples)
{
   int ii, index;

   index = subbandnum * numsamples;
   for (ii = 0; ii < numsamples; ii++, index++)
      chandat[ii] = (float) srawdata[index];
}

static int *get_prime_factors(int nn, int *numfactors)
{
   int ii, factors[30], *cfactors;

   *numfactors = 0;
   while ((nn & 1) == 0) {
      nn >>= 1;
      factors[(*numfactors)++] = 2;
   }
   if (nn == 1)
      goto finish;

   for (ii = 3; ii * ii <= nn; ii += 2) {
      while (!(nn % ii)) {
         nn /= ii;
         factors[(*numfactors)++] = ii;
      }
   }
   if (nn == 1)
      goto finish;

   factors[(*numfactors)++] = nn;

 finish:
   cfactors = (int *) malloc(*numfactors * sizeof(int));
   for (ii = 0; ii < *numfactors; ii++)
      cfactors[ii] = factors[ii];
   return cfactors;
}
