#include "presto.h"
#include "vectors.h"

char bands[NUMBANDS][40] = { "Radio", "IR", "Optical", "UV", "X-ray", "Gamma" };

char scopes[NUMSCOPES][40] =
    { "None (Artificial Data Set)", "Arecibo", "Parkes", "VLA",
   "MMT", "Las Campanas 2.5m", "Mt. Hopkins 48in", "Other"
};


void readinf(infodata * data, char *filenm)
{
   char tmp1[100], tmp2[100], *infofilenm;
   int ii;
   FILE *infofile;

   infofilenm = malloc(strlen(filenm)+5);
   sprintf(infofilenm, "%s.inf", filenm);
   infofile = chkfopen(infofilenm, "r");
   free(infofilenm);

   fscanf(infofile, "%*[^=]= %s", data->name);
   fscanf(infofile, "%*[^=]= %[^\n]\n", data->telescope);

   /* If not using makedata */

   if (strcmp(data->telescope, scopes[0]) != 0) {

      fscanf(infofile, "%*[^=]= %[^\n]\n", data->instrument);
      fscanf(infofile, "%*[^=]= %[^\n]\n", data->object);
      
      fscanf(infofile, "%*[^=]= %s\n", tmp1);
      ra_dec_from_string(tmp1, &data->ra_h, &data->ra_m, &data->ra_s);
      fscanf(infofile, "%*[^=]= %s\n", tmp1);
      ra_dec_from_string(tmp1, &data->dec_d, &data->dec_m, &data->dec_s);
      fscanf(infofile, "%*[^=]= %[^\n]\n", data->observer);
      fscanf(infofile, "%*[^=]= %d.%s", &data->mjd_i, tmp1);
      sprintf(tmp2, "0.%s", tmp1);
      data->mjd_f = strtod(tmp2, (char **) NULL);
/*     fscanf(infofile, "%*[^=]= %d.%lf", &data->mjd_i, \ */
/* 	   &data->mjd_f); */
/*     while (data->mjd_f > 1.0){ */
/*       data->mjd_f /= 10.0; */
/*     } */
      fscanf(infofile, "%*[^)] %*[^=]= %d", &data->bary);

   } else {

      strcpy(data->object, "fake pulsar");

   }

   fscanf(infofile, "%*[^=]= %lf", &data->N);
   fscanf(infofile, "%*[^=]= %lf", &data->dt);
   fscanf(infofile, "%*[^)] %*[^=]= %d\n", &data->numonoff);

   if (data->numonoff) {
      ii = 0;
      do {
         fscanf(infofile, "%*[^=]= %lf %*[ ,] %lf",
                &data->onoff[ii], &data->onoff[ii + 1]);
         ii += 2;
      } while (data->onoff[ii - 1] < data->N - 1 && ii < 2 * MAXNUMONOFF);
      data->numonoff = ii / 2;
      if (data->numonoff == MAXNUMONOFF) {
         printf("Number of onoff pairs (%d) is >= than MAXNUMONOFF (%d).\n",
                data->numonoff, MAXNUMONOFF);
         exit(1);
      }
   } else {
      data->numonoff = 1;
      data->onoff[0] = 0;
      data->onoff[1] = data->N - 1;
   }

   /* If not using makedata */

   if (strcmp(data->telescope, scopes[0]) != 0) {

      fscanf(infofile, "%*[^=]= %s", data->band);

      if (strcmp(data->band, bands[0]) == 0) {

         fscanf(infofile, "%*[^=]= %lf", &data->fov);
         fscanf(infofile, "%*[^=]= %lf", &data->dm);
         fscanf(infofile, "%*[^=]= %lf", &data->freq);
         fscanf(infofile, "%*[^=]= %lf", &data->freqband);
         fscanf(infofile, "%*[^=]= %d", &data->num_chan);
         fscanf(infofile, "%*[^=]= %lf", &data->chan_wid);

      } else if ((strcmp(data->band, bands[4]) == 0) ||
                 (strcmp(data->band, bands[5]) == 0)) {

         fscanf(infofile, "%*[^=]= %lf", &data->fov);
         fscanf(infofile, "%*[^=]= %lf", &data->energy);
         fscanf(infofile, "%*[^=]= %lf", &data->energyband);

      } else {

         fscanf(infofile, "%*[^=]= %s", data->filt);
         fscanf(infofile, "%*[^=]= %lf", &data->fov);
         fscanf(infofile, "%*[^=]= %lf", &data->wavelen);
         fscanf(infofile, "%*[^=]= %lf", &data->waveband);

      }

   }
   fscanf(infofile, "%*[^=]= %[^\n]\n", data->analyzer);
   fscanf(infofile, "%[^\n]\n", tmp1);
   fscanf(infofile, "%[^\n]\n", data->notes);
   fclose(infofile);
}



void writeinf(infodata * data)
{
   char tmp1[100], tmp2[100], *infofilenm;
   int itmp, ii;
   FILE *infofile;

   infofilenm = malloc(strlen(data->name)+5);
   sprintf(infofilenm, "%s.inf", data->name);
   infofile = chkfopen(infofilenm, "w");
   free(infofilenm);

   fprintf(infofile, " Data file name without suffix          =  %s\n", data->name);
   fprintf(infofile,
           " Telescope used                         =  %s\n", data->telescope);

   if (strcmp(data->telescope, scopes[0]) != 0) {       /* If using makedata */

      fprintf(infofile,
              " Instrument used                        =  %s\n", data->instrument);
      fprintf(infofile,
              " Object being observed                  =  %s\n", data->object);
      ra_dec_to_string(tmp1, data->ra_h, data->ra_m, data->ra_s);
      fprintf(infofile, " J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n", tmp1);
      ra_dec_to_string(tmp1, data->dec_d, data->dec_m, data->dec_s);
      fprintf(infofile, " J2000 Declination     (dd:mm:ss.ssss)  =  %s\n", tmp1);
      fprintf(infofile,
              " Data observed by                       =  %s\n", data->observer);
      sprintf(tmp1, "%.15f", data->mjd_f);
      sscanf(tmp1, "%d.%s", &itmp, tmp2);
      fprintf(infofile,
              " Epoch of observation (MJD)             =  %d.%s\n",
              data->mjd_i, tmp2);
      fprintf(infofile,
              " Barycentered?           (1=yes, 0=no)  =  %d\n", data->bary);

   }
   fprintf(infofile,
           " Number of bins in the time series      =  %-11.0f\n", data->N);
   fprintf(infofile, " Width of each time series bin (sec)    =  %.15g\n", data->dt);
   fprintf(infofile,
           " Any breaks in the data? (1=yes, 0=no)  =  %d\n",
           data->numonoff > 1 ? 1 : 0);

   if (data->numonoff > 1) {
      for (ii = 0; ii < data->numonoff; ii++) {
         fprintf(infofile,
                 " On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n",
                 ii + 1, data->onoff[2 * ii], data->onoff[2 * ii + 1]);
      }
   }

   if (strcmp(data->telescope, scopes[0]) != 0) {       /* If using makedata */

      fprintf(infofile,
              " Type of observation (EM band)          =  %s\n", data->band);

      if (strcmp(data->band, bands[0]) == 0) {

         fprintf(infofile,
                 " Beam diameter (arcsec)                 =  %.0f\n", data->fov);
         fprintf(infofile,
                 " Dispersion measure (cm-3 pc)           =  %.12g\n", data->dm);
         fprintf(infofile,
                 " Central freq of low channel (Mhz)      =  %.12g\n", data->freq);
         fprintf(infofile,
                 " Total bandwidth (Mhz)                  =  %.12g\n",
                 data->freqband);
         fprintf(infofile,
                 " Number of channels                     =  %d\n", data->num_chan);
         fprintf(infofile,
                 " Channel bandwidth (Mhz)                =  %.12g\n",
                 data->chan_wid);

      } else if ((strcmp(data->band, bands[4]) == 0) ||
                 (strcmp(data->band, bands[5]) == 0)) {

         fprintf(infofile,
                 " Field-of-view diameter (arcsec)        =  %.2f\n", data->fov);
         fprintf(infofile,
                 " Central energy (kev)                   =  %.1f\n", data->energy);
         fprintf(infofile,
                 " Energy bandpass (kev)                  =  %.1f\n",
                 data->energyband);

      } else {

         fprintf(infofile,
                 " Photometric filter used                =  %s\n", data->filt);
         fprintf(infofile,
                 " Field-of-view diameter (arcsec)        =  %.2f\n", data->fov);
         fprintf(infofile,
                 " Central wavelength (nm)                =  %.1f\n", data->wavelen);
         fprintf(infofile,
                 " Bandpass (nm)                          =  %.1f\n",
                 data->waveband);
      }

   }
   fprintf(infofile,
           " Data analyzed by                       =  %s\n", data->analyzer);

   fprintf(infofile, " Any additional notes:\n    %s\n\n", data->notes);

   fclose(infofile);

}
