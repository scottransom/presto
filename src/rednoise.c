#include "presto.h"
#include "math.h"
#include "rednoise_cmd.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
   FILE *infile, *outfile;
   int ii, bufflen = 10, numread_old, numread_new, numsamp, binnum = 1;
   int old_percent = 0, new_percent = 0;
   float *realbuffer = NULL, samprate, mean_old, mean_new, T, slope;
   fcomplex *inbuffer_old = NULL, *inbuffer_new = NULL, *outbuffer = NULL;
   char *rootfilenm, *outname;
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

#ifdef DEBUG
   showOptionValues();
#endif

   printf("\n\n");
   printf("     Rednoise Removal Routine\n");
   printf("          October, 2002    \n\n");

   {
      int hassuffix = 0;
      char *suffix;

      hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);

      if (hassuffix) {
         if (strcmp(suffix, "fft") != 0) {
            printf("\nInput file ('%s') must be a fourier transform ('.fft')!\n\n",
                   cmd->argv[0]);
            free(suffix);
            exit(0);
         }
         free(suffix);
      } else {
         printf("\nInput file ('%s') must be a fourier transform ('.fft')!\n\n",
                cmd->argv[0]);
         exit(0);
      }
      outname = (char *) calloc(strlen(rootfilenm) + 11, sizeof(char));
      sprintf(outname, "%s_red.fft", rootfilenm);
   }

   /* Read the info file */

   readinf(&idata, rootfilenm);
   samprate = idata.dt;
   numsamp = idata.N;
   T = numsamp * samprate;

   /* Open files and create arrays */

   infile = chkfopen(argv[1], "rb");
   outfile = chkfopen(outname, "wb");

   /* Read and remove rednoise */

   bufflen = cmd->startwidth;
   inbuffer_old = gen_cvect(cmd->endwidth);
   inbuffer_new = gen_cvect(cmd->endwidth);
   realbuffer = gen_fvect(cmd->endwidth);
   outbuffer = gen_cvect(cmd->endwidth);

   /* Takes care of the DC offset */
   chkfread(inbuffer_old, sizeof(fcomplex), 1, infile);
   inbuffer_old[0].r = 1.0;
   inbuffer_old[0].i = 0.0;
   chkfwrite(inbuffer_old, sizeof(fcomplex), 1, outfile);

   /* Calculates the first mean */
   numread_old = chkfread(inbuffer_old, sizeof(fcomplex), bufflen, infile);

   for (ii = 0; ii < numread_old; ii++) {
      realbuffer[ii] = 0;
      realbuffer[ii] = inbuffer_old[ii].r * inbuffer_old[ii].r +
          inbuffer_old[ii].i * inbuffer_old[ii].i;
   }

   mean_old = median(realbuffer, numread_old) / log(2.0);
   binnum += numread_old;
   bufflen = cmd->startwidth * log(binnum);

   while ((numread_new = chkfread(inbuffer_new, sizeof(fcomplex), bufflen, infile))) {
      for (ii = 0; ii < numread_new; ii++) {
         realbuffer[ii] = 0;
         realbuffer[ii] =
             inbuffer_new[ii].r * inbuffer_new[ii].r +
             inbuffer_new[ii].i * inbuffer_new[ii].i;
      }
      mean_new = median(realbuffer, numread_new) / log(2.0);
      slope = (mean_new - mean_old) / (numread_old + numread_new);

      /*
         printf("realbuffer[3] = %f ", realbuffer[3]);
         printf("Mean_old = %f ", mean_old);
         printf("Mean_new = %f ", mean_new);
         printf("Slope = %f\n", slope);
       */

      for (ii = 0; ii < numread_old; ii++) {
         outbuffer[ii].r = 0.0;
         outbuffer[ii].i = 0.0;
         outbuffer[ii].r =
             inbuffer_old[ii].r / sqrt(mean_old +
                                       slope * ((numread_old + numread_new) / 2.0 -
                                                ii));
         outbuffer[ii].i =
             inbuffer_old[ii].i / sqrt(mean_old +
                                       slope * ((numread_old + numread_new) / 2.0 -
                                                ii));
      }

      chkfwrite(outbuffer, sizeof(fcomplex), numread_old, outfile);

      /* printf("Binnum = %d\n", binnum); */
      binnum += numread_new;
      if ((binnum * 1.0) / T < cmd->endfreq)
         bufflen = cmd->startwidth * log(binnum);
      else
         /* exit(0); */
         bufflen = cmd->endwidth;

      numread_old = numread_new;
      mean_old = mean_new;

      for (ii = 0; ii < numread_new; ii++) {
         inbuffer_old[ii].r = 0.0;
         inbuffer_old[ii].i = 0.0;
         inbuffer_old[ii].r = inbuffer_new[ii].r;
         inbuffer_old[ii].i = inbuffer_new[ii].i;
      }

      /* Print percent complete */

      new_percent = (int) 100 *((binnum * 2.0) / numsamp);
      if (new_percent != old_percent) {
         printf("\rAmount Complete = %d%%", new_percent);
         old_percent = new_percent;
         fflush(stdout);
      }
   }

   /* Deal with the last chunk */

   for (ii = 0; ii < numread_old; ii++) {
      outbuffer[ii].r = 0;
      outbuffer[ii].i = 0;
      outbuffer[ii].r = inbuffer_old[ii].r / sqrt(mean_old);
      outbuffer[ii].i = inbuffer_old[ii].i / sqrt(mean_old);
   }
   chkfwrite(outbuffer, sizeof(fcomplex), numread_old, outfile);

   printf("\nDone.  Rednoise removed.\n\n");

   vect_free(inbuffer_old);
   vect_free(inbuffer_new);
   vect_free(realbuffer);
   vect_free(outbuffer);

   fclose(infile);
   fclose(outfile);
   free(rootfilenm);
   free(outname);
   exit(0);
}
