#include "presto.h"
#include "plot2d.h"
#include "profile_cmd.h"

/* The number of points to work with at a time      */
/* This must be the same as the WORKLEN in fold.c!  */
#define WORKLEN   16384

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char **argv)
/* profile: gets pulse profile from barycentered data  */
/* Written by Scott Ransom on 4 Jan 98 based on a vers */
/* written by Steve Eikenberry 1/30/95                 */
/* Added binary pulsar capability 30 Jan 98            */
/* Core completely re-written (faster, more accurate)  */
/* in April 1998.                                      */
/* And yet again in November 1999.  Yeesh.             */
{
   FILE *datafile, *proffile, *chifile, *filemarker;
   float *fprof = NULL, *chiarr = NULL, *times = NULL, *errors = NULL;
   double dt, orbdt, endphase, f = 0.0, fd = 0.0, fdd = 0.0;
   double epoch = 0.0, difft = 0.0, p = 0.0, pd = 0.0, pdd = 0.0;
   double *prof = NULL, endtime, N;
   double *Ep = NULL, *tp = NULL, startE = 0.0;
   double chip = 0.0, chiq = 0.0, chidf = 0.0;
   double chixmeas = 0.0, chitmp = 0.0;
   double varph = 0.0, numph = 0.0, avgph = 0.0;
   double normz = 0.0, normmean = 0.0, normstdev = 1.0;
   int chistatus = 0, chiwhich = 1, showerr = 0, flags = 0;
   double dbepoch = 0.0, onoffpairs[40], dtmp;
   int pnum, binary = 0, dochi = 0, numonoffpairs = 1;
   long numreads = 0, numpoints = 0;
   unsigned long filelen;
   long i = 0, proflen = 0;
   char datanm[200], profnm[200], psfilenm[200], chifilenm[200];
   char tmp1[100], tmp2[100], pname[30], *ctmp;
   orbitparams orb;
   psrparams psr;
   fourierprops rzwcand;
   makedata mdata;
   infodata idata;
   Cmdline *cmd;
   foldstats stats;

   /* Call usage() if we have no command line arguments */

   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(1);
   }
   /* Parse the command line using the excellent program Clig */

   cmd = parseCmdline(argc, argv);

#ifdef DEBUG
   showOptionValues();
#endif

   printf("\n\n");
   printf("          Pulse Folding Routine\n");
   printf("            by Scott M. Ransom\n");
   printf("             1 December, 1999\n\n");

   /* Initialize the filenames and some data: */

   sprintf(datanm, "%s.dat", cmd->argv[0]);
   sprintf(profnm, "%s.prof", cmd->argv[0]);
   sprintf(chifilenm, "%s.prof.chi", cmd->argv[0]);
   orb.p = 0.0;
   orb.x = 0.0;
   orb.e = 0.0;
   orb.w = 0.0;
   orb.t = 0.0;
   orb.wd = 0.0;
   orb.pd = 0.0;

   /* If we are going to compute a profile continue, otherwise */
   /* skip way foreward to the display section.....            */

   if (!cmd->dispP) {

      /* Open the data file */

      datafile = chkfopen(datanm, "rb");

      /* The number of points in datafile */

      filelen = chkfilelen(datafile, sizeof(float));
      numreads = filelen / WORKLEN;
      if (numreads % WORKLEN)
         numreads++;
      chiarr = gen_fvect(numreads);
      for (i = 0; i < numreads; i++)
         chiarr[i] = 0.0;

      /* Read the info file */

      readinf(&idata, cmd->argv[0]);
      if (idata.object) {
         printf("Folding a %s candidate from '%s'.\n", idata.object, datanm);
      } else {
         printf("Folding a candidate from '%s'.\n", datanm);
      }

      /* The MJD of the beginning of the observation */

      if (idata.mjd_i && idata.mjd_f) {
         epoch = (double) idata.mjd_i + idata.mjd_f;
      }
      N = idata.N;
      dt = idata.dt;
      endtime = 1.1 * dt * filelen;

      /* Read the pulsar database if needed */

      if (cmd->psrnameP) {
         pnum = get_psr_at_epoch(cmd->psrname, epoch, &psr);
         if (!pnum) {
            printf("The pulsar is not in the database.  Exiting.\n\n");
            exit(1);
         }
         if (psr.orb.p != 0.0) {
            binary = 1;
            orb = psr.orb;
            dbepoch = psr.orb.t / SECPERDAY;
         }
         p = psr.p;
         pd = psr.pd;
         pdd = psr.pdd;
         f = psr.f;
         fd = psr.fd;
         fdd = psr.fdd;
         strcpy(pname, psr.jname);

         /* If the user specifies all of the binaries parameters */

      } else if (cmd->binaryP) {

         /* Assume that the psr characteristics were measured at the time */
         /* of periastron passage (cmd->To)                               */

         difft = SECPERDAY * (epoch - cmd->To);
         orb.p = cmd->pb;
         orb.x = cmd->asinic;
         orb.e = cmd->e;
         orb.t = fmod(difft, orb.p);
         if (orb.t < 0.0)
            orb.t += orb.p;
         orb.w = (cmd->w + difft * cmd->wdot / SECPERJULYR);
         binary = 1;

         /* If the data set was generated using a makefile ('filename.mak')  */
         /* read it to determine the binary parameters.                      */

      } else if (cmd->makefileP) {

         epoch = 0.0;
         read_mak_file(cmd->argv[0], &mdata);
         binary = mdata.binary;
         p = mdata.p;
         pd = mdata.pd;
         pdd = mdata.pdd;
         orb = mdata.orb;
         f = mdata.f;
         fd = mdata.fd;
         fdd = mdata.fdd;
         dbepoch = mdata.orb.t / SECPERDAY;

         /* Determine the pulsar parameters to fold from a _rzw.cand file */

      } else if (cmd->rzwcandP) {

         if (!cmd->rzwfileP) {
            printf("\nYou must enter a name for the rzw candidate ");
            printf("file (-rzwfile filename)\n");
            printf("Exiting.\n\n");
            exit(1);
         }

         get_rzw_cand(cmd->rzwfile, cmd->rzwcand, &rzwcand);
         f = (rzwcand.r - 0.5 * rzwcand.z) / (dt * idata.N);
         p = 1.0 / f;
         fd = rzwcand.z / ((dt * idata.N) * (dt * idata.N));
         pd = -fd / (f * f);
      }

      /* Determine the pulsar parameters to fold if we are not getting   */
      /* the data from a .cand file, the pulsar database, or a makefile. */

      if (!cmd->rzwcandP && !cmd->makefileP && !cmd->psrnameP) {

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
      }

      /* Determine the length of the profile */

      if (cmd->proflenP) {
         proflen = cmd->proflen;
      } else {
         proflen = (long) (p / dt + 0.5);
      }

      /* Determine the frequency shifts caused by the orbit if needed */

      if (binary) {

         /* Save the orbital solution every half second */

         startE = keplers_eqn(orb.t, orb.p, orb.e, 1.0E-15);
         if (endtime > 2048)
            orbdt = 0.5;
         else
            orbdt = endtime / 4096.0;
         numpoints = (long) floor(endtime / orbdt + 0.5) + 1;
         Ep = dorbint(startE, numpoints, orbdt, &orb);
         tp = gen_dvect(numpoints);
         for (i = 0; i < numpoints; i++)
            tp[i] = i * orbdt;

         /* Convert Eccentric anomaly to time delays */

         orb.w *= DEGTORAD;
         E_to_phib(Ep, numpoints, &orb);
      }

      /* Determine if we will take any breaks during the folding */

      for (i = 0; i < 40; i++)
         onoffpairs[i] = 0.0;
      i = 0;
      if (cmd->onoffP) {
         numonoffpairs = 0;
         ctmp = strtok(cmd->onoff, " \t\n");
         while (NULL != ctmp) {
            onoffpairs[i] = strtod(ctmp, NULL);
            if (onoffpairs[i] < 0.0 || onoffpairs[i] > 1.0) {
               printf("\nonoff pairs must be between 0.0 and 1.0 inclusive.\n\n");
               exit(1);
            }
            if (i >= 1 && onoffpairs[i] < onoffpairs[i - 1]) {
               printf("\nonoff values must increase from left to right.\n\n");
               exit(1);
            }
            if (onoffpairs[i] > 1.0)
               ctmp = NULL;
            if (!(i & 1))
               numonoffpairs++;
            i++;
            ctmp = strtok(NULL, " \t\n");
         }
         if (i & 1)
            onoffpairs[i] = 1.0;

         /* Adjust the epoch of the beginning of the first bin */

         epoch += (onoffpairs[0] * N * dt) / SECPERDAY;
      } else {
         onoffpairs[0] = 0.0;
         onoffpairs[1] = 1.0;
      }

      /* Output some informational data on the screen and to the */
      /* output file.                                            */

      proffile = chkfopen(profnm, "w");
      fprintf(stdout, "\n");
      filemarker = stdout;
      for (i = 0; i < 2; i++) {
         if (cmd->psrnameP)
            fprintf(filemarker, "Pulsar                       =  %s\n", pname);
         if (epoch != 0.0)
            fprintf(filemarker, "Folding (Obs) epoch   (MJD)  =  %-17.11f\n", epoch);
         fprintf(filemarker, "Data pt duration (dt)   (s)  =  %-.12f\n", dt);
         fprintf(filemarker, "Total number of data points  =  %-.0f\n", N);
         fprintf(filemarker, "Number of profile bins       =  %ld\n", proflen);
         fprintf(filemarker, "Folding period          (s)  =  %-.15f\n", p);
         if (pd != 0.0)
            fprintf(filemarker, "Folding p-dot         (s/s)  =  %-.10e\n", pd);
         if (pdd != 0.0)
            fprintf(filemarker, "Folding p-dotdot    (s/s^2)  =  %-.10e\n", pdd);
         fprintf(filemarker, "Folding frequency      (hz)  =  %-.12f\n", f);
         if (fd != 0.0)
            fprintf(filemarker, "Folding f-dot        (hz/s)  =  %-.8e\n", fd);
         if (fdd != 0.0)
            fprintf(filemarker, "Folding f-dotdot   (hz/s^2)  =  %-.8e\n", fdd);
         if (binary) {
            fprintf(filemarker, "Orbital period          (s)  =  %-.10f\n", orb.p);
            fprintf(filemarker, "a*sin(i)/c (x)     (lt-sec)  =  %-.10f\n", orb.x);
            fprintf(filemarker, "Eccentricity                 =  %-.10f\n", orb.e);
            fprintf(filemarker,
                    "Longitude of peri (w) (deg)  =  %-.10f\n", orb.w / DEGTORAD);
            if (cmd->psrnameP)
               fprintf(filemarker, "Epoch of periapsis    (MJD)  =  %-17.11f\n",
                       epoch - dbepoch);
            else if (cmd->makefileP || cmd->To)
               fprintf(filemarker,
                       "Epoch of periapsis      (d)  =  %-17.11f\n",
                       -(orb.t - orb.p) / SECPERDAY);
         }
         filemarker = proffile;
      }

      /* Begin the profile generation */
      /* Prep variables for the folding */

      prof = gen_dvect(proflen);
      for (i = 0; i < proflen; i++)
         prof[i] = 0.0;

      /* The heart of the routine: */

      if (binary)
         flags = 3;
      else
         flags = 2;
      endphase = foldfile(datafile, dt, 0.0, prof, proflen, cmd->phs,
                          f, fd, fdd, flags, Ep, tp, numpoints,
                          onoffpairs, &stats, chiarr);
      fclose(datafile);

      /* The total number of "photons"... */

      for (i = 0; i < proflen; i++)
         numph += prof[i];

      /* Average value of a profile bin... */

      avgph = stats.prof_avg;

      /* Variance of a profile bin... */

      varph = stats.prof_var;

      /* Compute the Chi-Squared probability that there is a signal */
      /* See Leahy et al., ApJ, Vol 266, pp. 160-170, 1983 March 1. */

      chixmeas = chisqr(prof, proflen, avgph, varph);
      times = gen_freqs(numreads, 0.0, 1.0 / (float) numreads);
      dochi = 1;

      /* Calculate the values of P and Q since we know X and DF */

      chidf = proflen - 1;
      cdfchi(&chiwhich, &chip, &chiq, &chixmeas, &chidf, &chistatus, &chitmp);
      if (chistatus != 0) {
         if (chistatus < 0)
            printf("\nOne of the input parameters to cdfchi() was out of range.\n");
         else if (chistatus == 3)
            printf("\nP + Q do not equal 1.0 in cdfchi().\n");
         else if (chistatus == 10)
            printf("\nError in cdfgam().\n");
         else
            printf("\nUnknown error in cdfchi().\n");
      }

      /* Write the Chi-Array to a text output file */
      chifile = chkfopen(chifilenm, "w");
      for (i = 0; i < numreads; i++)
         fprintf(chifile, "%9.7f    %9f\n", times[i], chiarr[i]);
      fclose(chifile);

      /* Output our statistics */

      filemarker = stdout;
      for (i = 0; i < 2; i++) {
         fprintf(filemarker, "Number of points folded      =  %-.0f\n",
                 stats.numdata);
         fprintf(filemarker, "Total signal          (cts)  =  %-15.2f\n", numph);
         fprintf(filemarker, "Time series average          =  %-14.6f\n",
                 stats.data_avg);
         fprintf(filemarker, "Time series variance         =  %-14.6f\n",
                 stats.data_var);
         fprintf(filemarker, "Profile average              =  %-14.6f\n", avgph);
         fprintf(filemarker, "Expected profile variance    =  %-14.6f\n", varph);
         fprintf(filemarker, "Ending pulse phase           =  %-.9f\n", endphase);
         fprintf(filemarker, "Degrees of freedom           =  %-.0f\n", chidf);
         fprintf(filemarker, "Measured chi-square   (X^2)  =  %-12.4f\n", chixmeas);
         if (chiq < 0.00001) {
            fprintf(filemarker, "Prob to randomly exceed X^2  =  %-.9g\n", chiq);
         } else {
            fprintf(filemarker, "Prob to randomly exceed X^2  =  %-.9f\n", chiq);
         }
         filemarker = proffile;
      }
      fprintf(stdout, "\n\n");
      fprintf(proffile, "-------------------------------------------------\n");

      /* Write the profile to the output file */

      for (i = 0; i < proflen; i++)
         fprintf(proffile, "%5ld    %9f\n", i + 1, prof[i]);
      fclose(proffile);

      /* Convert the double precision profile to single precision */
      /* in order to plot it.                                     */

      fprof = gen_fvect(proflen);
      for (i = 0; i < proflen; i++)
         fprof[i] = (float) prof[i];
      vect_free(prof);

   } else {

      /* Use an already written profile file... */

      proffile = chkfopen(profnm, "r");

      /* Read and get the number of points */

      do {
         fgets(tmp1, sizeof(tmp1), proffile);
         printf("%s", tmp1);
         sscanf(tmp1, "%[^=]s", tmp2);
         if (!strcmp(tmp2, "Number of profile bins       ")) {
            sscanf(tmp1, "Number of profile bins       =  %ld", &proflen);
            binary = 1;
         }
      } while (!binary);

      /* Generate the profile vector */

      fprof = gen_fvect(proflen);

      /* Read through the header section */

      binary = 0;
      do {
         fgets(tmp1, sizeof(tmp1), proffile);
         if (!strcmp(tmp1, "-------------------------------------------------\n")) {
            binary = 1;
         } else {
            printf("%s", tmp1);
         }
      } while (!binary);

      i = 0;
      while (fscanf(proffile, "%s %s", tmp1, tmp2) != EOF) {
         binary = atoi(tmp1);
         fprof[i] = strtod(tmp2, NULL);
         i++;
      }
      printf("\n");
   }

   /* Calculate the errors for each profile bin */

   if (!cmd->dispP) {
      showerr = 1;
      errors = gen_fvect(proflen);
      for (i = 0; i < proflen; i++)
         errors[i] = sqrt(varph);

      /* Re-calculate chi-squared if needed */

      if (chixmeas == 0.0) {
         avgph = 0.0;
         for (i = 0; i < proflen; i++)
            avgph += fprof[i];
         avgph /= proflen;
         for (i = 0; i < proflen; i++) {
            dtmp = fprof[i];
            chitmp = dtmp - avgph;
            chixmeas += (chitmp * chitmp);
         }
         chixmeas /= varph;
      }

      /* Calculate the values of P and Q since we know X and DF */

      chidf = proflen - 1;
      cdfchi(&chiwhich, &chip, &chiq, &chixmeas, &chidf, &chistatus, &chitmp);
      if (chistatus != 0) {
         if (chistatus < 0)
            printf("\nInput parameter %d to cdfchi() was out of range.\n",
                   chistatus);
         else if (chistatus == 3)
            printf("\nP + Q do not equal 1.0 in cdfchi().\n");
         else if (chistatus == 10)
            printf("\nError in cdfgam().\n");
         else
            printf("\nUnknown error in cdfchi().\n");
      }
      chixmeas /= (proflen - 1.0);

      /* Calculate the equivalent sigma */

      chiwhich = 2;
      cdfnor(&chiwhich, &chip, &chiq, &normz, &normmean, &normstdev,
             &chistatus, &chitmp);
      if (chistatus != 0) {
         if (chistatus < 0)
            printf("\nInput parameter %d to cdfnor() was out of range.\n",
                   chistatus);
         else if (chistatus == 3)
            printf("\nP + Q do not equal 1.0 in cdfnor().\n");
         else
            printf("\nUnknown error in cdfnor().\n");
      }
   }

   /* Initialize PGPLOT using Postscript if requested  */

   if (cmd->noerrP)
      showerr = 0;
   else if (proflen > 100)
      showerr = 0;
   if (cmd->psP || cmd->bothP) {
      sprintf(psfilenm, "%s.prof.ps", cmd->argv[0]);
      cpgstart_ps(psfilenm, "landscape");
      sprintf(tmp1, "p = %16.13f s,  p-dot = %12.7g", p, pd);
      if (!cmd->dispP) {
         sprintf(tmp2, "Reduced Chi = %.3f,  P(Noise) < %.2g,  sigma = %.2f",
                 chixmeas, chiq, normz);
         if (idata.object)
            plot_profile(proflen, fprof, idata.object, tmp1, tmp2,
                         showerr, errors, 1);
         else
            plot_profile(proflen, fprof, "", tmp1, tmp2, showerr, errors, 1);
      } else {
         plot_profile(proflen, fprof, "", tmp1, "", showerr, errors, 1);
      }
      if (dochi)
         xyline(numreads, times, chiarr, "Time", "Reduced Chi-Squared", 1);
      cpgend();
   } else if (cmd->xwinP) {
      cpgstart_x("landscape");
      sprintf(tmp1, "p = %16.13f s,  p-dot = %12.7g", p, pd);
      if (!cmd->dispP) {
         sprintf(tmp2, "Reduced Chi = %.3f,  P(Noise) < %.2g,  sigma = %.2f",
                 chixmeas, chiq, normz);
         if (idata.object)
            plot_profile(proflen, fprof, idata.object, tmp1, tmp2,
                         showerr, errors, 1);
         else
            plot_profile(proflen, fprof, "", tmp1, tmp2, showerr, errors, 1);
      } else {
         plot_profile(proflen, fprof, "", tmp1, "", showerr, errors, 1);
      }
      if (dochi)
         xyline(numreads, times, chiarr, "Time", "Reduced Chi-Squared", 1);
      cpgend();
   }

   /* Send the plot to the screen if necessary */

   if (cmd->bothP) {
      cpgstart_x("landscape");
      sprintf(tmp1, "p = %16.13f s,  p-dot = %12.7g", p, pd);
      if (!cmd->dispP) {
         sprintf(tmp2, "Reduced Chi = %.3f,  P(Noise) < %.2g,  sigma = %.2f",
                 chixmeas, chiq, normz);
         if (idata.object)
            plot_profile(proflen, fprof, idata.object, tmp1, tmp2,
                         showerr, errors, 1);
         else
            plot_profile(proflen, fprof, "", tmp1, tmp2, showerr, errors, 1);
      } else {
         plot_profile(proflen, fprof, "", tmp1, "", showerr, errors, 1);
      }
      if (dochi)
         xyline(numreads, times, chiarr, "Time", "Reduced Chi-Squared", 1);
      cpgend();
   }

   /* Cleanup */

   if (cmd->makefileP)
      if (mdata.onoff)
         free(mdata.onoff);
   if (!cmd->dispP && binary) {
      vect_free(Ep);
      vect_free(tp);
   }
   if (!cmd->dispP)
      vect_free(errors);
   vect_free(times);
   vect_free(fprof);
   vect_free(chiarr);
   return (0);
}
