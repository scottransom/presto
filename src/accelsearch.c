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

  /* Read the info file */

  readinf(&idata, cmd->argv[0]);
  if (idata.object) {
    printf("Analyzing %s data from '%s'.\n\n", 
	   remove_whitespace(idata.object), cmd->argv[0]);
  } else {
    printf("Analyzing data from '%s'.\n\n", cmd->argv[0]);
  }

  /* Create the accelobs structure */
  
  create_accelobs(&obs, &idata, cmd);

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

  free_accelobs(obs);
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
