#include <glib.h>
#include "accel.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *fftfile, *accelcandfile, *fourierpropsfile;
  int nextbin;
  float nph;
  char *rootfilenm, *candnm, *accelcandnm, *accelnm;
  GSList *cands=NULL;
  subharminfo *subharminf;
  accelobs obs;
  fourierprops *props;
  infodata idata;
  struct tms runtimes;
  double ttim, utim, stim, tott;
  Cmdline *cmd;

  /* Prep the timer */

  tott = times(&runtimes) / (double) CLK_TCK;

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
  printf("       Pulsar Acceleration Search Routine II\n");
  printf("              by Scott M. Ransom\n");
  printf("                February, 2001\n\n");

  {
    int hassuffix=0;
    char *suffix;
    
    hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);
    if (hassuffix){
      if (strcmp(suffix, "fft")!=0){
	printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
	       cmd->argv[0]);
	free(suffix);
	exit(0);
      }
      free(suffix);
    } else {
      printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n",
	     cmd->argv[0]);
      exit(0);
    }
  }
  
  /* Read the info file */

  readinf(&idata, cmd->argv[0]);
  if (idata.object) {
    printf("Analyzing %s data from '%s'.\n\n", 
	   remove_whitespace(idata.object), cmd->argv[0]);
  } else {
    printf("Analyzing data from '%s'.\n\n", cmd->argv[0]);
  }

  /* Determine the output filenames */

  candnm = (char *)calloc(strlen(rootfilenm)+25, 1);
  accelcandnm = (char *)calloc(strlen(rootfilenm)+25, 1);
  accelnm = (char *)calloc(strlen(rootfilenm)+25, 1);
  sprintf(candnm, "%s_ACCEL_%d.cand", rootfilenm, cmd->zmax);
  sprintf(accelcandnm, "%s_ACCEL_%d.accelcand", rootfilenm, cmd->zmax);
  sprintf(accelnm, "%s_ACCEL_%d", rootfilenm, cmd->zmax);

  /* Create the accelobs structure */
  
  fftfile = chkfopen(cmd->argv[0], "rb");
  create_accelobs(fftfile, &obs, &idata, cmd);

  /* Generate the correlation kernels */

  printf("Generating fdot kernels for the correlations...");
  fflush(NULL);
  subharminf = create_subharminfo_vect(cmd->numharm, cmd->zmax);
  printf("Done.\n\n");

  /* Start the main search loop */

  nextbin = cmd->rlo;

  do {

    startbin = nextbin;

    /* Get the data from the file */

    filedata = read_fcomplex_file(fftfile, 
				  startbin - kern_half_width,
				  filedatalen);

    /* Get approximate local power statistics */
    
    powlist = gen_fvect(filedatalen);
    
    /* Calculate the powers */
    
    for (ii = 0; ii < filedatalen; ii++) 
      powlist[ii] = POWER(filedata[ii].r, filedata[ii].i);
    
    /* Set the local power level to 1.442695 * median value.    */
    /* The 1.442695 corrects the median to the mean for an      */
    /* exponential distribution.  Then take the reciprocal so   */
    /* that we multiply instead of divide during normalization. */
    
    locpow = 1.0 / (1.442695 * median(powlist, filedatalen));
    free(powlist);

    /*  Do the f-fdot plane correlations: */

    for (zct = 0; zct < nz; zct++) {

      /* Calculate percentage complete */

      newper = (int) (totnumsearched / (numr * 0.01)) + 1;

      if (newper > oldper) {
	newper = (newper > 99) ? 100 : newper;
	printf("\rAmount of search complete = %3d%%", newper);
	fflush(stdout);
	oldper = newper;
      }

      if (zct == 0) datainf = RAW;
      else datainf = SAME;

      /* Perform the correlation */

      nr = corr_complex(filedata, filedatalen, datainf, \
			kernels[zct], corrsize, FFT, \
			corrdata, corrsize, kern_half_width, \
			numbetween, kern_half_width, CORR);
      nextbin = startbin + nr / numbetween;
      worknumbins = (nextbin > highestbin) ? \
	(highestbin - startbin) * numbetween : nr;
      
      /* This loop is the heart of the search */

      for (ii = 0; ii < worknumbins; ii++) {
	chkpow = POWER(corrdata[ii].r, corrdata[ii].i) * locpow;

	/* Check if the measured power is greater than cutoff */

	if (chkpow > minpow) {
	  newpos.pow = chkpow;
	  newpos.p1 = startbin + ii * dr;
	  newpos.p2 = cmd->zlo + zct * dz;
	  newpos.p3 = 0.0;

	  /* If there is a zapfile, check to see if our candidate */
	  /* matches one of the 'birdies'.  If it does, continue. */

	  if (cmd->zapfileP && 
	      check_to_zap(newpos.p1, zapfreqs, zapwidths, numzap))
	    continue;

	  /* Check to see if another candidate with these properties */
	  /* is already in the list.                                 */

	  spot = not_already_there_rzw(&newpos, list, ncand);

	  if (spot >= 0) {
	    list[spot] = newpos;
	    minpow = percolate(list, ncand, spot);
	  }
	}
      }
      totnumsearched += worknumbins;
    }
    free(filedata);
  } while (nextbin <= highestbin);

  /* Free the memory used by the correlation kernels */

  free(kernels[0]);
  free(kernels);

  printf("\rAmount of search complete = %3d%%", 100);
  fflush(stdout);
  printf("\nDone searching.  ");
  printf("Now optimizing each candidate and sorting.\n\n");

  /* Do rough duplicate removal (probably not necessary) */

  newncand = ncand;
  newncand -= remove_dupes(list, ncand);

  /* Save the list of 'rough' candidates to a file */

  poscandfile = chkfopen(poscandnm, "w");
  chkfwrite(list, sizeof(position), newncand, poscandfile);
  fclose(poscandfile);

  /* Now maximize each good candidate */

  newper = 0;
  oldper = 0;

  for (ii = 0; ii < newncand; ii++) {

    /* Calculate percentage complete */

    newper = (int) (ii / (float) (newncand) * 100.0) + 1;
    if (newper > oldper) {
      printf("\rAmount of optimization complete = %3d%%", newper);
      fflush(stdout);
      oldper = newper;
    }
    hipow = max_rz_file(fftfile, list[ii].p1, list[ii].p2, \
			&hir, &hiz, &derivs[ii]);
    calc_props(derivs[ii], hir + cmd->lobin, hiz, 0.0, &props[ii]);
  }
  printf("\rAmount of optimization complete = %3d%%\n\n", 100);

  qsort(props, (unsigned long) newncand, sizeof(fourierprops), \
	compare_fourierprops);

  /* Do fine scale duplicate removal and other cleaning */

  newncand -= remove_dupes2(props, newncand);
  newncand -= remove_other(props, newncand, cmd->rlo, highestbin, lowpowlim, 
			   cmd->zapfileP, zapfreqs, zapwidths, numzap);

  /* Set our candidate notes to all spaces */

  notes = malloc(sizeof(char) * newncand * 20 + 1);
  for (ii = 0; ii < newncand; ii++)
    strncpy(notes + ii * 20, "                         ", 20);

  /* Compare the candidates with the pulsar database */

  if (idata.ra_h && idata.dec_d) {
    for (ii = 0; ii < newncand; ii++) {
      comp_psr_to_cand(&props[ii], &idata, notes + ii * 20, 0);
    }
  }

  /* Compare the candidates with themselves */

  compare_rzw_cands(props, newncand, notes);

  /* Write the binary candidate file */

  candfile = chkfopen(candnm, "wb");
  chkfwrite(props, sizeof(fourierprops), (unsigned long) newncand, \
	    candfile);
  fclose(candfile);

  /* Send the candidates to the text file */

  if (cmd->ncand < newncand) newncand = cmd->ncand;
  file_reg_candidates(props, notes, newncand, dt, \
		      (long) (N + DBLCORRECT), nph, rootfilenm, accelnm);

  /* Finish up */

  printf("Done.\n\n");
  printf("Searched %.0f pts (approximately %.0f were independent).\n\n", \
	 totnumsearched, totnumsearched * 0.5 * dz / 6.95);

  printf("Timing summary:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n", \
	 ttim, utim, stim);
  printf("  Total time: %.3f sec\n\n", tott);

  printf("Candidates in binary format are stored in '%s'.\n", candnm);
  printf("A candidate Postscript table is stored in '%s.ps'.\n\n", accelnm);

/*     readinf(&idata, infonm); */
/*     realpsr = comp_psr_to_cand(&props, &idata, compare); */
/*     printf("%s\n",compare); */

  fclose(fftfile);
  free(candnm);
  free(accelcandnm);
  free(accelnm);
  free(rootfilenm);
  free(list);
  free(derivs);
  free(props);
  free(notes);
  if (cmd->zapfileP){
    free(zapfreqs);
    free(zapwidths);
  }
  return (0);
}


int not_already_there_rzw(position * newpos, position * list, int nlist)
{
  int ii;

  /* Loop through the candidates already in the list */

  for (ii = 0; ii < nlist; ii++) {
    if (list[ii].pow == 0.0)
      break;
    
    /* Check to see if the positions of the candidates match */

    if ((fabs(newpos->p1 - list[ii].p1) < 0.51) &&
	(fabs(newpos->p2 - list[ii].p2) < 2.01) &&
	(fabs(newpos->p3 - list[ii].p3) < 5.0)){
      
      /* If the previous candidate is simply a lower power   */
      /* version of the new candidate, overwrite the old and */
      /* percolate it to its proper position.                */

      if (list[ii].pow < newpos->pow){
	return ii;
	
      /* Otherwise, skip the new candidate.  Its already there. */

      }	else return -1;
    }
  }

  /* Place the new candidate in the last position */

  return nlist - 1;
}


void compare_rzw_cands(fourierprops * list, int nlist, char *notes)
{
  int ii, jj, kk;
  char tmp[30];

  /* Loop through the candidates (reference cands) */

  for (ii = 0; ii < nlist; ii++) {

    /* Loop through the candidates (referenced cands) */

    for (jj = 0; jj < nlist; jj++) {
      if (ii == jj)
	continue;

      /* Look for standard sidelobes */

      if (fabs(list[ii].r - list[jj].r) < 15.0 && \
	  fabs(list[ii].z - list[jj].z) > 1.0 && \
	  list[ii].pow > list[jj].pow) {

	/* Check if the note has already been written */

	sprintf(tmp, "%.20s", notes + jj * 20);
	if (!strcmp("                    ", tmp)) {

	  /* Write the note */

	  sprintf(notes + jj * 20, "SL? of Cand %d", ii + 1);
	}
	continue;
      }
      /* Loop through the possible PSR period harmonics */

      for (kk = 1; kk < 61; kk++) {

	/* Check if the PSR Fourier freqs and z's are close enough */

	if ((fabs(list[ii].r - list[jj].r / kk) < list[jj].rerr * 3) && \
	    (fabs(list[ii].z - list[jj].z / kk) < list[jj].zerr * 2)) {

	  /* Check if the note has already been written */

	  sprintf(tmp, "%.20s", notes + jj * 20);
	  if (!strcmp("                    ", tmp)) {

	    /* Write the note */

	    sprintf(notes + jj * 20, "H %d of Cand %d", kk, ii + 1);

	    break;
	  }
	}
      }
    }
  }
}
#undef SHORTESTFFT
#undef LOSKIP
