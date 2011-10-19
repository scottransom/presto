#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include "presto.h"
#include "mpiprepsubband_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "mpi.h"
#include "gmrt.h"
#include "spigot.h"
#include "psrfits.h"
#include "sigproc_fb.h"

#define RAWDATA (cmd->pkmbP || cmd->bcpmP || cmd->wappP || cmd->gmrtP || cmd->spigotP || cmd->filterbankP || cmd->psrfitsP)

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 20.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

extern void write_data(FILE * outfiles[], int numfiles, float **outdata,
                       int startpoint, int numtowrite);
extern void write_padding(FILE * outfiles[], int numfiles, float value,
                          int numtowrite);
extern void update_infodata(infodata * idata, int datawrote, int padwrote,
                            int *barybins, int numbarybins, int downsamp);
extern void print_percent_complete(int current, int number);
extern void make_infodata_struct(void);
extern void make_maskbase_struct(void);
extern void broadcast_mask(mask * obsmask, int myid);

extern void get_PKMB_static(int *decreasing_freqs, float *clip_sigma,
                            int *bytesperpt, int *offsetbytes);
extern void set_PKMB_static(int ptsperblk, int bytesperpt,
                            int numchan, int decreasing_freqs, 
                            int offsetbytes, float clip_sigma, double dt);
extern void get_BCPM_static(int *bytesperpt, int *bytesperblk, int *numifs,
                            int *chan_map, float *clip_sigma);
extern void set_BCPM_static(int ptsperblk, int bytesperpt, int bytesperblk,
                            int numchan, int numifs, float clip_sigma,
                            double dt, int *chan_map);
extern void get_SPIGOT_static(int *bytesperpt, int *bytesperblk,
                              int *numifs, float *clip_sigma);
extern void set_SPIGOT_static(int ptsperblk, int bytesperpt,
                              int bytesperblk, int numchan, int numifs,
                              float clip_sigma, double dt);
extern void get_WAPP_static(int *bytesperpt, int *bytesperblk, int *numifs,
                            float *clip_sigma);
extern void set_WAPP_static(int ptsperblk, int bytesperpt, int bytesperblk,
                            int numchan, int numifs, float clip_sigma, double dt);
extern void get_GMRT_static(int *bytesperpt, int *bytesperblk, float *clip_sigma);
extern void set_GMRT_static(int ptsperblk, int bytesperpt, int bytesperblk,
                            int numchan, float clip_sigma, double dt);
extern void get_filterbank_static(int *ptsperbyte, int *bytesperpt,
                                  int *bytesperblk, float *clip_sigma);
extern void set_filterbank_static(int ptsperbyte, int ptsperblk,
                                  int bytesperpt, int bytesperblk,
                                  int numchan, float clip_sigma, double dt);
static int get_data(FILE * infiles[], int numfiles, float **outdata,
                    mask * obsmask, float *padvals, double dt,
                    double *dispdts, int **offsets, int *padding);

MPI_Datatype infodata_type;
MPI_Datatype maskbase_type;

static Cmdline *cmd;
static IFs ifs = SUMIFS;
static int blocklen = 0, blocksperread = 0, bytesperblk = 0, worklen =
    0, numchan = 1, insubs = 0;
static int local_numdms = 1, myid = 0, numprocs = 1;
static PKMB_tapehdr hdr;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   /* Any variable that begins with 't' means topocentric */
   /* Any variable that begins with 'b' means barycentric */
   FILE **infiles = NULL, **outfiles = NULL;
   float **outdata, *padvals;
   double dtmp, *dms, avgdm = 0.0, dsdt = 0, maxdm;
   double *dispdt, tlotoa = 0.0, blotoa = 0.0, BW_ddelay = 0.0;
   double max = -9.9E30, min = 9.9E30, var = 0.0, avg = 0.0;
   double *btoa = NULL, *ttoa = NULL, avgvoverc = 0.0;
   char obs[3], ephem[10], rastring[50], decstring[50];
   int numinfiles, totnumtowrite, **offsets;
   int ii, jj, numadded = 0, numremoved = 0, padding = 0, good_inputs = 1;
   int numbarypts = 0, numread = 0, numtowrite = 0, totwrote = 0, datawrote = 0;
   int padwrote = 0, padtowrite = 0, statnum = 0;
   int numdiffbins = 0, *diffbins = NULL, *diffbinptr = NULL, good_padvals = 0;
   int offsetbytes = 0; // For packed 1-bit Parkes data
   double local_lodm;
   char *datafilenm, *outpath, *outfilenm, *hostname;
   infodata idata;
   mask obsmask;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   {
      FILE *hostfile;
      char tmpname[100];
      int retval;

      hostfile = chkfopen("/etc/hostname", "r");
      retval = fscanf(hostfile, "%s\n", tmpname);
      if (retval==0) {
          printf("Warning:  error reading /etc/hostname on proc %d\n", myid);
      }
      hostname = (char *) calloc(strlen(tmpname) + 1, 1);
      memcpy(hostname, tmpname, strlen(tmpname));
      fclose(hostfile);
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

   make_infodata_struct();
   make_maskbase_struct();

   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);
   numinfiles = cmd->argc;
   if (cmd->noclipP)
      cmd->clip = 0.0;

#ifdef DEBUG
   showOptionValues();
#endif

   if (myid == 0) {             /* Master node only */
      printf("\n\n");
      printf("      Parallel Pulsar Subband De-dispersion Routine\n");
      printf("                 by Scott M. Ransom\n\n");
      if (!RAWDATA) {
         char *root, *suffix;
         /* Split the filename into a rootname and a suffix */
         if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
            printf("\nThe input filename (%s) must have a suffix!\n\n",
                   cmd->argv[0]);
            good_inputs = 0;
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
            } else if (strcmp(suffix, "fits") == 0 &&
                       strstr(root, "spigot_5") != NULL) {
               printf("Assuming the data is from the NRAO/Caltech Spigot card...\n");
               cmd->spigotP = 1;
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
            } else {
               printf
                   ("\nCannot determine the format of the input files '%s'...\n\n",
                    cmd->argv[0]);
               good_inputs = 0;
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
         if ((infiles[ii] = fopen(cmd->argv[ii], "rb")) == NULL) {
            perror("\nError in fopen()");
            printf("   filename = '%s'\n", cmd->argv[ii]);
            good_inputs = 0;
         }
      }
      if (!cmd->numoutP)
         cmd->numout = INT_MAX;
   }
   MPI_Bcast(&cmd->pkmbP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->gmrtP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->bcpmP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->wappP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->spigotP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->psrfitsP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->filterbankP, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&cmd->numout, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&insubs, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (insubs)
      cmd->nsub = cmd->argc;

   MPI_Bcast(&good_inputs, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (!good_inputs) {
      MPI_Finalize();
      exit(1);
   }

   /* Determine the output file names and open them */

   local_numdms = cmd->numdms / (numprocs - 1);
   dms = gen_dvect(local_numdms);
   if (cmd->numdms % (numprocs - 1)) {
      if (myid == 0)
         printf
             ("\nThe number of DMs must be divisible by (the number of processors - 1).\n\n");
      MPI_Finalize();
      exit(1);
   }
   local_lodm = cmd->lodm + (myid - 1) * local_numdms * cmd->dmstep;

   split_path_file(cmd->outfile, &outpath, &outfilenm);
   datafilenm = (char *) calloc(strlen(outfilenm) + 20, 1);
   if (myid > 0) {
      if (chdir(outpath) == -1) {
         printf("\nProcess %d cannot chdir() to '%s'.  Exiting.\n\n", myid, outpath);
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

   /* Read an input mask if wanted */

   if (cmd->maskfileP) {
      if (myid == 0) {
         read_mask(cmd->maskfile, &obsmask);
         printf("Read mask information from '%s'\n\n", cmd->maskfile);
         good_padvals = determine_padvals(cmd->maskfile, &obsmask, &padvals);
      }
      broadcast_mask(&obsmask, myid);
      if (myid != 0)
         padvals = gen_fvect(obsmask.numchan);
      MPI_Bcast(&good_padvals, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(padvals, obsmask.numchan, MPI_FLOAT, 0, MPI_COMM_WORLD);
   } else {
      obsmask.numchan = obsmask.numint = 0;
      MPI_Bcast(&good_padvals, 1, MPI_INT, 0, MPI_COMM_WORLD);
   }

   {
      float clip_sigma = 0.0;
      double dt, T;
      int ptsperbyte, ptsperblk, bytesperpt, numifs = 0, decreasing_freqs = 1;
      int chan_mapping[2 * MAXNUMCHAN];
      long long N;

      if (myid == 0) {          /* Master */

         /* Set-up values if we are reading subband data  */

         if (insubs) {
            char *root, *suffix, *tmpname;
            /* Split the filename into a rootname and a suffix */
            if (split_root_suffix(cmd->argv[0], &root, &suffix) == 0) {
               printf
                   ("\nThe input filename (%s) must have a suffix!\n\n",
                    cmd->argv[0]);
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
            T = idata.N * idata.dt;
            N = (long long) idata.N;
            dt = idata.dt;
            ptsperblk = SUBSBLOCKLEN;
            bytesperpt = sizeof(short); /* Subbands are shorts */
            bytesperblk = cmd->nsub * bytesperpt * ptsperblk;
         }

         /* Set-up values if we are using the Parkes multibeam */

         if (cmd->pkmbP) {
            printf("\nFilterbank input file information:\n");
            get_PKMB_file_info(infiles, numinfiles, cmd->clip, 
                               &N, &ptsperblk, &numchan, &dt, &T, 1);
            get_PKMB_static(&decreasing_freqs, &clip_sigma,
                            &bytesperpt, &offsetbytes);
            bytesperblk = DATLEN;
            chkfread(&hdr, 1, HDRLEN, infiles[0]);
            rewind(infiles[0]);
            PKMB_hdr_to_inf(&hdr, &idata);
            PKMB_update_infodata(numinfiles, &idata);
         }

         /* Set-up values if we are using the GMRT Phased Array system */

         if (cmd->gmrtP) {
            printf("\nGMRT input file information:\n");
            get_GMRT_file_info(infiles, argv + 1, numinfiles,
                               cmd->clip, &N, &ptsperblk, &numchan, &dt, &T, 1);
            get_GMRT_static(&bytesperpt, &bytesperblk, &clip_sigma);
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
                                     &N, &ptsperblk, &numchan, &dt, &T, 1);
            get_filterbank_static(&ptsperbyte, &bytesperpt, &bytesperblk, &clip_sigma);
            filterbank_update_infodata(numinfiles, &idata);
            set_filterbank_padvals(padvals, good_padvals);
         }

         /* Set-up values if we are using the Berkeley-Caltech */
         /* Pulsar Machine (or BPP) format.                    */

         if (cmd->bcpmP) {
            printf("\nBCPM input file information:\n");
            get_BPP_file_info(infiles, numinfiles, cmd->clip, &N,
                              &ptsperblk, &numchan, &dt, &T, &idata, 1);
            get_BCPM_static(&bytesperpt, &bytesperblk, &numifs,
                            chan_mapping, &clip_sigma);
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
            get_SPIGOT_file_info(infiles, spigots, numinfiles,
                                 cmd->windowP, cmd->clip, &N,
                                 &ptsperblk, &numchan, &dt, &T, &idata, 1);
            get_SPIGOT_static(&bytesperpt, &bytesperblk, &numifs, &clip_sigma);
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
            ptsperblk = s.spectra_per_subint;
            numchan = s.num_channels;
            dt = s.dt;
            T = s.T;
            get_PSRFITS_file_info(cmd->argv, cmd->argc, cmd->clip, 
                                  &s, &idata, 1);
            get_PSRFITS_static(&bytesperpt, &bytesperblk, &numifs, &clip_sigma);
            PSRFITS_update_infodata(&idata);
            set_PSRFITS_padvals(padvals, good_padvals);
         }

         /* Set-up values if we are using the Arecobo WAPP */

         if (cmd->wappP) {
            printf("\nWAPP input file information:\n");
            get_WAPP_file_info(infiles, cmd->numwapps, numinfiles,
                               cmd->windowP, cmd->clip,
                               &N, &ptsperblk, &numchan, &dt, &T, &idata, 1);
            get_WAPP_static(&bytesperpt, &bytesperblk, &numifs, &clip_sigma);
            WAPP_update_infodata(numinfiles, &idata);
            set_WAPP_padvals(padvals, good_padvals);
         }

         /* The number of topo to bary time points to generate with TEMPO */
         numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;

         // Identify the TEMPO observatory code
         {
             char *outscope = (char *) calloc(40, sizeof(char));
             telescope_to_tempocode(idata.telescope, outscope, obs);
             free(outscope);
         }
      }

      MPI_Bcast(&ptsperbyte, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ptsperblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&bytesperpt, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&bytesperblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&numchan, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&numifs, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&numbarypts, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&decreasing_freqs, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(chan_mapping, 2 * MAXNUMCHAN, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&clip_sigma, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (cmd->pkmbP) {  // Pass the packed bits info for special Parkes data
          MPI_Bcast(&offsetbytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
      blocklen = ptsperblk;

      if (myid > 0) {           /* Slave */
         if (cmd->pkmbP)
            set_PKMB_static(ptsperblk, bytesperpt, numchan, 
                            decreasing_freqs, offsetbytes, 
                            clip_sigma, dt);
         if (cmd->gmrtP) {
            set_GMRT_static(ptsperblk, bytesperpt, bytesperblk,
                            numchan, clip_sigma, dt);
            set_GMRT_padvals(padvals, good_padvals);
         }
         if (cmd->filterbankP) {
            set_filterbank_static(ptsperbyte, ptsperblk, bytesperpt, bytesperblk,
                                  numchan, clip_sigma, dt);
            set_filterbank_padvals(padvals, good_padvals);
         }
         if (cmd->bcpmP) {
            set_BCPM_static(ptsperblk, bytesperpt, bytesperblk,
                            numchan, numifs, clip_sigma, dt, chan_mapping);
            set_BPP_padvals(padvals, good_padvals);
         }
         if (cmd->spigotP) {
            set_SPIGOT_static(ptsperblk, bytesperpt, bytesperblk,
                              numchan, numifs, clip_sigma, dt);
            set_SPIGOT_padvals(padvals, good_padvals);
         }
         if (cmd->psrfitsP) {
            set_PSRFITS_static(ptsperblk, bytesperpt, bytesperblk,
                               numchan, numifs, clip_sigma, dt);
            set_PSRFITS_padvals(padvals, good_padvals);
         }
         if (cmd->wappP) {
            set_WAPP_static(ptsperblk, bytesperpt, bytesperblk,
                            numchan, numifs, clip_sigma, dt);
            set_WAPP_padvals(padvals, good_padvals);
         }
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
      /* For the WAPP (and others) , the number of bytes returned in         */
      /* get_WAPP_rawblock() is ptsperblk since the correlator lags          */
      /* are converted to 1 byte filterbank channels in read_WAPP_rawblock() */
      if (cmd->wappP || cmd->gmrtP || cmd->filterbankP || cmd->spigotP || cmd->psrfitsP) {
         bytesperblk = ptsperblk * numchan;
         bytesperpt = numchan;
      }
   }

   /* Broadcast or calculate a few extra important values */

   MPI_Bcast(&idata, 1, infodata_type, 0, MPI_COMM_WORLD);
   if (insubs) {
      avgdm = idata.dm;
      numchan = idata.num_chan;
   }
   idata.dm = avgdm;
   dsdt = cmd->downsamp * idata.dt;
   maxdm = cmd->lodm + cmd->numdms * cmd->dmstep;
   BW_ddelay = delay_from_dm(maxdm, idata.freq) - 
      delay_from_dm(maxdm, idata.freq + (idata.num_chan-1) * idata.chan_wid);
   blocksperread = ((int) (BW_ddelay / idata.dt) / blocklen + 1);
   worklen = blocklen * blocksperread;

   if (blocklen % cmd->downsamp) {
      if (myid == 0) {
         printf
             ("Error:  The downsample factor (%d) must be a factor of the\n",
              cmd->downsamp);
         printf("        blocklength (%d).  Exiting.\n\n", blocklen);
      }
      MPI_Finalize();
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

      offsets = gen_imatrix(local_numdms, cmd->nsub);
      for (ii = 0; ii < local_numdms; ii++) {
         double *subdispdt;

         subdispdt = subband_delays(numchan, cmd->nsub, dms[ii],
                                    idata.freq, idata.chan_wid, 0.0);
         dtmp = subdispdt[cmd->nsub - 1];
         for (jj = 0; jj < cmd->nsub; jj++)
            offsets[ii][jj] = (int) ((subdispdt[jj] - dtmp) / dsdt + 0.5);
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

      {                         /* Print the nodes and the DMs they are handling */
         int kk;

         MPI_Barrier(MPI_COMM_WORLD);
         fflush(NULL);
         if (myid == 0)
            printf("Node\t\tDMs\n----\t\t---\n");
         fflush(NULL);

         for (jj = 1; jj < numprocs; jj++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == jj) {
               printf("%s\t\t", hostname);
               for (kk = 0; kk < local_numdms - 1; kk++)
                  printf("%.2f\t", dms[kk]);
               printf("%.2f\n", dms[kk]);
            }
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);
         }
         if (myid == 0)
            printf("\n");
         fflush(NULL);
      }

      outdata = gen_fmatrix(local_numdms, worklen / cmd->downsamp);
      numread = get_data(infiles, numinfiles, outdata,
                         &obsmask, padvals, idata.dt,
                         dispdt, offsets, &padding);

      while (numread == worklen) {

         numread /= cmd->downsamp;
         if (myid == 0)
            print_percent_complete(totwrote, totnumtowrite);

         /* Write the latest chunk of data, but don't   */
         /* write more than cmd->numout points.         */

         numtowrite = numread;
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         if (myid > 0) {
            write_data(outfiles, local_numdms, outdata, 0, numtowrite);
            /* Update the statistics */
            if (!padding) {
               for (ii = 0; ii < numtowrite; ii++)
                  update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg, &var);
               statnum += numtowrite;
            }
         }
         totwrote += numtowrite;

         /* Stop if we have written out all the data we need to */

         if (cmd->numoutP && (totwrote == cmd->numout))
            break;

         numread = get_data(infiles, numinfiles, outdata,
                            &obsmask, padvals, idata.dt,
                            dispdt, offsets, &padding);
      }
      datawrote = totwrote;

   } else {                     /* Main loop if we are barycentering... */

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
      for (ii = 0; ii < numbarypts; ii++)
         ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

      /* Call TEMPO for the barycentering */

      if (myid == 0) {
         double maxvoverc = -1.0, minvoverc = 1.0, *voverc = NULL;

         printf("\nGenerating barycentric corrections...\n");
         voverc = gen_dvect(numbarypts);
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

         printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
         printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
         printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
         printf("De-dispersing using %d subbands.\n", cmd->nsub);
         if (cmd->downsamp > 1)
            printf("Downsampling by a factor of %d (new dt = %.10g)\n",
                   cmd->downsamp, dsdt);
         printf("\n");
      }

      {                         /* Print the nodes and the DMs they are handling */
         int kk;

         MPI_Barrier(MPI_COMM_WORLD);
         fflush(NULL);
         if (myid == 0)
            printf("Node\t\tDMs\n----\t\t---\n");
         fflush(NULL);

         for (jj = 1; jj < numprocs; jj++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == jj) {
               printf("%s\t\t", hostname);
               for (kk = 0; kk < local_numdms - 1; kk++)
                  printf("%.2f\t", dms[kk]);
               printf("%.2f\n", dms[kk]);
            }
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);
         }
         if (myid == 0)
            printf("\n");
         fflush(NULL);
      }

      MPI_Bcast(btoa, numbarypts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&avgvoverc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blotoa = btoa[0];

      /* Dispersion delays (in bins).  The high freq gets no delay   */
      /* All other delays are positive fractions of bin length (dt)  */

      dispdt = subband_search_delays(numchan, cmd->nsub, avgdm,
                                     idata.freq, idata.chan_wid, avgvoverc);
      for (ii = 0; ii < numchan; ii++)
         dispdt[ii] /= idata.dt;

      /* The subband dispersion delays (see note above) */

      offsets = gen_imatrix(local_numdms, cmd->nsub);
      for (ii = 0; ii < local_numdms; ii++) {
         double *subdispdt;

         subdispdt = subband_delays(numchan, cmd->nsub, dms[ii],
                                    idata.freq, idata.chan_wid, avgvoverc);
         dtmp = subdispdt[cmd->nsub - 1];
         for (jj = 0; jj < cmd->nsub; jj++)
            offsets[ii][jj] = (int) ((subdispdt[jj] - dtmp) / dsdt + 0.5);
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
                  *diffbinptr =
                      NEAREST_INT(LININTERP
                                  (calcpt, btoa[ii - 1], btoa[ii], lobin, hibin));
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

      outdata = gen_fmatrix(local_numdms, worklen / cmd->downsamp);
      numread = get_data(infiles, numinfiles, outdata,
                         &obsmask, padvals, idata.dt,
                         dispdt, offsets, &padding);

      while (numread == worklen) {      /* Loop to read and write the data */
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
         if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
            numtowrite = cmd->numout - totwrote;
         if (numtowrite > numread)
            numtowrite = numread;
         if (myid > 0) {
            write_data(outfiles, local_numdms, outdata, 0, numtowrite);
            /* Update the statistics */
            if (!padding) {
               for (ii = 0; ii < numtowrite; ii++)
                  update_stats(statnum + ii, outdata[0][ii], &min, &max, &avg, &var);
               statnum += numtowrite;
            }
         }
         datawrote += numtowrite;
         totwrote += numtowrite;
         numwritten += numtowrite;

         if ((datawrote == abs(*diffbinptr)) && (numwritten != numread) && (totwrote < cmd->numout)) {  /* Add/remove a bin */
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
               if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
                  numtowrite = cmd->numout - totwrote;
               nextdiffbin = abs(*diffbinptr) - datawrote;
               if (numtowrite > nextdiffbin)
                  numtowrite = nextdiffbin;
               if (myid > 0) {
                  write_data(outfiles, local_numdms, outdata, skip, numtowrite);
                  /* Update the statistics and counters */
                  if (!padding) {
                     for (ii = 0; ii < numtowrite; ii++)
                        update_stats(statnum + ii,
                                     outdata[0][skip + ii], &min, &max, &avg, &var);
                     statnum += numtowrite;
                  }
               }
               numwritten += numtowrite;
               datawrote += numtowrite;
               totwrote += numtowrite;
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
                            &obsmask, padvals, idata.dt,
                            dispdt, offsets, &padding);
      }
   }

   if (myid > 0) {

      /* Calculate the amount of padding we need  */

      if (cmd->numoutP && (cmd->numout > totwrote))
         padwrote = padtowrite = cmd->numout - totwrote;

      /* Write the new info file for the output data */

      idata.dt = dsdt;
      update_infodata(&idata, totwrote, padtowrite, diffbins,
                      numdiffbins, cmd->downsamp);
      for (ii = 0; ii < local_numdms; ii++) {
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
   }

   /* Close the files and cleanup */

   if (cmd->maskfileP)
      free_mask(obsmask);
   if (myid == 0) {
      for (ii = 0; ii < numinfiles; ii++)
         fclose(infiles[ii]);
      free(infiles);
   } else {
      for (ii = 0; ii < local_numdms; ii++)
         fclose(outfiles[ii]);
      free(outfiles);
   }
   vect_free(outdata[0]);
   vect_free(outdata);
   vect_free(dms);
   free(hostname);
   vect_free(dispdt);
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


static int read_subbands(FILE * infiles[], int numfiles, short *subbanddata)
/* Read short int subband data written by prepsubband */
{
   int ii, jj, index, numread = 0;
   double run_avg;

   for (ii = 0; ii < numfiles; ii++) {
      index = ii * SUBSBLOCKLEN;
      numread = chkfread(subbanddata + index, sizeof(short),
                         SUBSBLOCKLEN, infiles[ii]);

      if (cmd->runavgP==1) {
          run_avg = 0.0;
          for(jj = 0; jj < numread; jj++)
              run_avg += (float) subbanddata[jj+index];
          run_avg /= numread;
          for(jj = 0; jj < numread; jj++)
              subbanddata[jj+index] = (float) subbanddata[jj+index] - run_avg;
      }

      for (jj = numread; jj < SUBSBLOCKLEN; jj++)
         subbanddata[index + jj] = 0.0;
   }
   return numread;
}

static void convert_subbands(int numfiles, short *shortdata,
                             float *subbanddata, double timeperblk,
                             int *maskchans, int *nummasked, mask * obsmask,
                             float clip_sigma, float *padvals)
/* Convert and transpose the subband data, then mask it*/
{
   int ii, jj, index, shortindex, offset, channum, mask = 0;
   double starttime;
   float subband_sum;
   static int currentblock = 0;

   *nummasked = 0;
   if (obsmask->numchan) mask = 1;

   for (ii = 0; ii < numfiles; ii++) {
      shortindex = ii * SUBSBLOCKLEN;
      for (jj = 0, index = ii; jj < SUBSBLOCKLEN; jj++, index += numfiles)
         subbanddata[index] = (float) shortdata[shortindex + jj];
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
   if (mask) {
       fflush(stdout);
       if (*nummasked == -1) {   /* If all channels are masked */
           for (ii = 0; ii < SUBSBLOCKLEN; ii++)
               memcpy(subbanddata + ii * numfiles, padvals, sizeof(float) * numfiles);
           fflush(stdout);
       } else if (*nummasked > 0) {      /* Only some of the channels are masked */
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
}


static int get_data(FILE * infiles[], int numfiles, float **outdata,
                    mask * obsmask, float *padvals, double dt,
                    double *dispdts, int **offsets, int *padding)
{
   static int firsttime = 1, worklen, *maskchans = NULL, blocksize;
   static int dsworklen;
   static float *tempzz, *data1, *data2, *dsdata1 = NULL, *dsdata2 = NULL;
   static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
   static double blockdt;
   static unsigned char *rawdata = NULL;
   int totnumread = 0, numread = 0, tmpnumread = 0;
   int ii, jj, tmppad = 0, nummasked = 0;

   if (firsttime) {
      /* For rawdata, we need to make two initial reads in order to        */
      /* prepare the prep_*_subbands() functions as well as this routine.  */
      if (RAWDATA)
         firsttime = 2;
      if (cmd->maskfileP)
         maskchans = gen_ivect(numchan);
      worklen = blocklen * blocksperread;
      dsworklen = worklen / cmd->downsamp;
      blocksize = blocklen * cmd->nsub;
      blockdt = blocklen * dt;
      data1 = gen_fvect(cmd->nsub * worklen);
      data2 = gen_fvect(cmd->nsub * worklen);
      rawdata = gen_bvect(bytesperblk);
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
      if (RAWDATA) {
         for (ii = 0; ii < blocksperread; ii++) {
            if (myid == 0) {
               if (cmd->pkmbP)
                  numread =
                      read_PKMB_rawblock(infiles, numfiles, &hdr, rawdata, &tmppad);
               else if (cmd->gmrtP)
                  numread = read_GMRT_rawblock(infiles, numfiles, rawdata, &tmppad);
               else if (cmd->filterbankP)
                  numread =
                      read_filterbank_rawblock(infiles, numfiles, rawdata, &tmppad);
               else if (cmd->bcpmP)
                  numread = read_BPP_rawblock(infiles, numfiles, rawdata, &tmppad);
               else if (cmd->wappP)
                  numread =
                      read_WAPP_rawblock(infiles, numfiles, rawdata, &tmppad, ifs);
               else if (cmd->spigotP)
                  numread =
                      read_SPIGOT_rawblock(infiles, numfiles, rawdata, &tmppad, ifs);
               else if (cmd->psrfitsP)
                  numread =
                      read_PSRFITS_rawblock(rawdata, &tmppad);
               numread *= blocklen;
            }
            MPI_Bcast(&numread, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&tmppad, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(rawdata, bytesperblk, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
            if (myid > 0) {
               if (numread) {
                  if (cmd->pkmbP)
                     tmpnumread =
                         prep_PKMB_subbands(rawdata,
                                            currentdata +
                                            ii * blocksize, dispdts,
                                            cmd->nsub, 0, maskchans,
                                            &nummasked, obsmask);
                  else if (cmd->gmrtP)
                     tmpnumread =
                         prep_GMRT_subbands(rawdata,
                                            currentdata +
                                            ii * blocksize, dispdts,
                                            cmd->nsub, 0, maskchans,
                                            &nummasked, obsmask);
                  else if (cmd->filterbankP)
                     tmpnumread =
                         prep_filterbank_subbands(rawdata,
                                                  currentdata +
                                                  ii * blocksize,
                                                  dispdts,
                                                  cmd->nsub, 0,
                                                  maskchans, &nummasked, obsmask);
                  else if (cmd->bcpmP)
                     tmpnumread =
                         prep_BPP_subbands(rawdata,
                                           currentdata +
                                           ii * blocksize, dispdts,
                                           cmd->nsub, 0, maskchans,
                                           &nummasked, obsmask, ifs);
                  else if (cmd->wappP)
                     tmpnumread =
                         prep_WAPP_subbands(rawdata,
                                            currentdata +
                                            ii * blocksize, dispdts,
                                            cmd->nsub, 0, maskchans,
                                            &nummasked, obsmask);
                  else if (cmd->spigotP)
                     tmpnumread =
                         prep_SPIGOT_subbands(rawdata,
                                              currentdata +
                                              ii * blocksize,
                                              dispdts, cmd->nsub, 0,
                                              maskchans, &nummasked, obsmask);
                  else if (cmd->psrfitsP)
                     tmpnumread =
                         prep_PSRFITS_subbands(rawdata,
                                               currentdata +
                                               ii * blocksize,
                                               dispdts, cmd->nsub, 0,
                                               maskchans, &nummasked, obsmask);
               } else {
                  *padding = 1;
                  for (jj = ii * blocksize; jj < (ii + 1) * blocksize; jj++)
                     currentdata[jj] = 0.0;
               }
               if (tmppad)
                  *padding = 1;
            }
            if (!firsttime)
               totnumread += numread;
         }
      } else if (insubs) {
         short *subsdata = NULL;

         subsdata = gen_svect(SUBSBLOCKLEN * numfiles);
         for (ii = 0; ii < blocksperread; ii++) {
            if (myid == 0)
               numread = read_subbands(infiles, numfiles, subsdata);
            MPI_Bcast(&numread, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(subsdata, SUBSBLOCKLEN * numfiles, MPI_SHORT, 0,
                      MPI_COMM_WORLD);
            convert_subbands(numfiles, subsdata,
                             currentdata + ii * blocksize, blockdt,
                             maskchans, &nummasked, obsmask,
                             cmd->clip, padvals);
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
   /*{
      int jj;
      for (jj=0; jj<numprocs; jj++){
      if (myid==jj)
      printf("%d:  %d  %d  %d\n", myid, numread, totnumread, worklen);
      fflush(NULL);
      MPI_Barrier(MPI_COMM_WORLD);
      }
      } */
   if (totnumread != worklen) {
     if (cmd->maskfileP) {
       vect_free(maskchans);
       vect_free(padvals);
     }
      vect_free(data1);
      vect_free(data2);
      vect_free(rawdata);
      if (cmd->downsamp > 1) {
         vect_free(dsdata1);
         vect_free(dsdata2);
      }
   }
   return totnumread;
}
