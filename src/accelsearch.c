#include "accel.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static void print_percent_complete(int current, int number)
{
  static int newper=0, oldper=-1;

  newper = (int) (current / (float)(number) * 100.0);
  if (newper < 0) newper = 0;
  if (newper > 100) newper = 100;
  if (newper > oldper) {
    printf("\rAmount of optimization complete = %3d%%", newper);
    fflush(stdout);
    oldper = newper;
  }
}


int main(int argc, char *argv[])
{
  double ttim, utim, stim, tott;
  char *rootfilenm, *candnm, *accelcandnm, *accelnm, *workfilenm;
  struct tms runtimes;
  subharminfo *subharminf;
  accelobs obs;
  infodata idata;
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
  workfilenm = (char *)calloc(strlen(rootfilenm)+25, 1);
  sprintf(candnm, "%s_ACCEL_%d.cand", rootfilenm, cmd->zmax);
  sprintf(accelcandnm, "%s_ACCEL_%d.accelcand", rootfilenm, cmd->zmax);
  sprintf(accelnm, "%s_ACCEL_%d", rootfilenm, cmd->zmax);
  sprintf(workfilenm, "%s_ACCEL_%d.txtcand", rootfilenm, cmd->zmax);

  /* Create the accelobs structure */
  
  create_accelobs(chkfopen(cmd->argv[0], "rb"),
		  chkfopen(workfilenm, "w"),
		  &obs, &idata, cmd);

  /* Generate the correlation kernels */

  printf("Generating fdot kernels for the correlations...\n");
  fflush(NULL);
  subharminf = create_subharminfo_vect(cmd->numharm, cmd->zmax);
  printf("Done.\n\n");

  /* Start the main search loop */

  {
    int harm_to_sum=1;
    double startr=obs.rlo, lastr=0, nextr=0;
    GSList *cands=NULL;
    ffdotpows *fundamental;
    
    while (startr < obs.highestbin){  /* Search the fundamental */
      harm_to_sum = 1;
      nextr = startr + ACCEL_USELEN * ACCEL_DR;
      lastr = nextr - ACCEL_DR;
      fundamental = subharm_ffdot_plane(harm_to_sum, 1, startr, lastr, 
					&subharminf[harm_to_sum], &obs);
      search_ffdotpows(fundamental, harm_to_sum, &obs, cands);

      if (obs.numharm > 1){   /* Search the subharmonics */
	ffdotpows *fundamental_copy, *subharmonic;

	harm_to_sum = 2;
	fundamental_copy = copy_ffdotpows(fundamental);

	/* Search the 1/2 subharmonic (i.e. 1/2 fundamental) */

	subharmonic = subharm_ffdot_plane(harm_to_sum, 1, startr, lastr, 
					  &subharminf[harm_to_sum], &obs);
	add_ffdotpows(fundamental, subharmonic, harm_to_sum, 1);
	free_ffdotpows(subharmonic);
	search_ffdotpows(fundamental, harm_to_sum, &obs, cands);

	/* Search the 1/4 subharmonic by building on the 1/2 */

	if (obs.numharm == 4){
	  harm_to_sum = 4;
	  subharmonic = subharm_ffdot_plane(harm_to_sum, 1, startr, lastr, 
					    &subharminf[harm_to_sum], &obs);
	  add_ffdotpows(fundamental, subharmonic, harm_to_sum, 1);
	  free_ffdotpows(subharmonic);
	  subharmonic = subharm_ffdot_plane(harm_to_sum, 3, startr, lastr, 
					    &subharminf[harm_to_sum], &obs);
	  add_ffdotpows(fundamental, subharmonic, harm_to_sum, 3);
	  free_ffdotpows(subharmonic);
	  search_ffdotpows(fundamental, harm_to_sum, &obs, cands);
	}
	
	/* Search the 1/3 subharmonic (work from scratch) */

	if (obs.numharm >= 3){
	  harm_to_sum = 3;
	  free_ffdotpows(fundamental);
	  fundamental = copy_ffdotpows(fundamental_copy);
	  subharmonic = subharm_ffdot_plane(harm_to_sum, 1, startr, lastr, 
					    &subharminf[harm_to_sum], &obs);
	  add_ffdotpows(fundamental, subharmonic, harm_to_sum, 1);
	  free_ffdotpows(subharmonic);
	  subharmonic = subharm_ffdot_plane(harm_to_sum, 2, startr, lastr, 
					    &subharminf[harm_to_sum], &obs);
	  add_ffdotpows(fundamental, subharmonic, harm_to_sum, 2);
	  free_ffdotpows(subharmonic);
	  search_ffdotpows(fundamental, harm_to_sum, &obs, cands);
	}
	free_ffdotpows(fundamental_copy);
      }
      free_ffdotpows(fundamental);
      startr = nextr;
    }
  }

  printf("\nDone searching.  ");
  printf("Now optimizing each candidate.\n\n");

  {
    double hir, hiz;
    int ii, jj, numcands;
    accelcand *cand;
    rderivs

    /* Now maximize each good candidate */
    
    numcands = g_slist_length(cands);
    for (ii=0; ii<numcands; ii++){
      print_percent_complete(ii, numcands);
      cand = (accelcand *)g_slist_nth_data(cands, ii);
      calc_props(derivs[ii], hir + cmd->lobin, hiz, 0.0, &props[ii]);
    }
    print_percent_complete(ii, numcands);
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
  free(workfilenm);
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
