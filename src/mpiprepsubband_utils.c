#include <limits.h>
#include "presto.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "mpi.h"

MPI_Datatype infodata_type;
MPI_Datatype maskbase_type;

typedef struct MASKBASE {
   double timesigma;            /* Cutoff time-domain sigma               */
   double freqsigma;            /* Cutoff freq-domain sigma               */
   double mjd;                  /* MJD of time zero                       */
   double dtint;                /* Time in sec of each interval           */
   double lofreq;               /* Freq (MHz) of low channel              */
   double dfreq;                /* Channel width (MHz)                    */
   int numchan;                 /* Number of channels                     */
   int numint;                  /* Number of intervals                    */
   int ptsperint;               /* Points per interval                    */
   int num_zap_chans;           /* Number of full channels to zap         */
   int num_zap_ints;            /* Number of full intervals to zap        */
} maskbase;


void make_maskbase_struct(void)
{
   int blockcounts[2] = { 6, 5 };
   MPI_Datatype types[2] = { MPI_DOUBLE, MPI_INT };
   MPI_Aint displs[2];
   maskbase mbase;

   MPI_Address(&mbase.timesigma, &displs[0]);
   MPI_Address(&mbase.numchan, &displs[1]);
   displs[1] -= displs[0];
   displs[0] = 0;
   MPI_Type_struct(2, blockcounts, displs, types, &maskbase_type);
   MPI_Type_commit(&maskbase_type);
}


void broadcast_mask(mask * obsmask, int myid)
{
   int ii;
   maskbase mbase;

   if (myid == 0) {
      mbase.timesigma = obsmask->timesigma;
      mbase.freqsigma = obsmask->freqsigma;
      mbase.mjd = obsmask->mjd;
      mbase.dtint = obsmask->dtint;
      mbase.lofreq = obsmask->lofreq;
      mbase.dfreq = obsmask->dfreq;
      mbase.numchan = obsmask->numchan;
      mbase.numint = obsmask->numint;
      mbase.ptsperint = obsmask->ptsperint;
      mbase.num_zap_chans = obsmask->num_zap_chans;
      mbase.num_zap_ints = obsmask->num_zap_ints;
   }
   MPI_Bcast(&mbase, 1, maskbase_type, 0, MPI_COMM_WORLD);
   if (myid > 0) {
      obsmask->timesigma = mbase.timesigma;
      obsmask->freqsigma = mbase.freqsigma;
      obsmask->mjd = mbase.mjd;
      obsmask->dtint = mbase.dtint;
      obsmask->lofreq = mbase.lofreq;
      obsmask->dfreq = mbase.dfreq;
      obsmask->numchan = mbase.numchan;
      obsmask->numint = mbase.numint;
      obsmask->ptsperint = mbase.ptsperint;
      obsmask->num_zap_chans = mbase.num_zap_chans;
      obsmask->num_zap_ints = mbase.num_zap_ints;
      obsmask->zap_chans = gen_ivect(mbase.num_zap_chans);
      obsmask->zap_ints = gen_ivect(mbase.num_zap_ints);
      obsmask->num_chans_per_int = gen_ivect(mbase.numint);
      obsmask->chans = (int **) malloc(mbase.numint * sizeof(int *));
   }
   MPI_Bcast(obsmask->zap_chans, mbase.num_zap_chans, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(obsmask->zap_ints, mbase.num_zap_ints, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(obsmask->num_chans_per_int, mbase.numint, MPI_INT, 0, MPI_COMM_WORLD);
   for (ii = 0; ii < mbase.numint; ii++) {
      if (obsmask->num_chans_per_int[ii] > 0 &&
          obsmask->num_chans_per_int[ii] < obsmask->numchan) {
         if (myid > 0)
            obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
         MPI_Bcast(obsmask->chans[ii], obsmask->num_chans_per_int[ii],
                   MPI_INT, 0, MPI_COMM_WORLD);
         //int jj;
         //printf("%d: ", myid);
         //for (jj=0; jj<obsmask->num_chans_per_int[ii]; jj++)
         //   printf("%d ", obsmask->chans[ii][jj]);
         //printf("\n");
      } else if (obsmask->num_chans_per_int[ii] == obsmask->numchan) {
         int jj;
         //printf("%d: all zapped\n", myid);
         if (myid > 0) {
            obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
            for (jj = 0; jj < obsmask->numchan; jj++)
               obsmask->chans[ii][jj] = jj;
         }
      }
   }
}


void make_infodata_struct(void)
{
   int ii, blockcounts[11] = { MAXNUMONOFF * 2 + 14, 8, 500, 200, 100,
      100, 100, 100, 40, 40, 7
   };
   MPI_Datatype types[11] = { MPI_DOUBLE, MPI_INT, MPI_CHAR, MPI_CHAR,
      MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR,
      MPI_CHAR, MPI_CHAR, MPI_CHAR
   };
   MPI_Aint displs[11];
   infodata idata;

   MPI_Address(&idata.ra_s, &displs[0]);
   MPI_Address(&idata.num_chan, &displs[1]);
   MPI_Address(&idata.notes, &displs[2]);
   MPI_Address(&idata.name, &displs[3]);
   MPI_Address(&idata.object, &displs[4]);
   MPI_Address(&idata.instrument, &displs[5]);
   MPI_Address(&idata.observer, &displs[6]);
   MPI_Address(&idata.analyzer, &displs[7]);
   MPI_Address(&idata.telescope, &displs[8]);
   MPI_Address(&idata.band, &displs[9]);
   MPI_Address(&idata.filt, &displs[10]);
   for (ii = 10; ii >= 0; ii--)
      displs[ii] -= displs[0];
   MPI_Type_struct(11, blockcounts, displs, types, &infodata_type);
   MPI_Type_commit(&infodata_type);
}


void write_data(FILE * outfiles[], int numfiles, float **outdata,
                int startpoint, int numtowrite)
{
   int ii;

   for (ii = 0; ii < numfiles; ii++)
      chkfwrite(outdata[ii] + startpoint, sizeof(float), numtowrite, outfiles[ii]);
}


void write_padding(FILE * outfiles[], int numfiles, float value, int numtowrite)
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


void print_percent_complete(int current, int number)
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


void update_infodata(infodata * idata, int datawrote, int padwrote,
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
