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
  int ii;
  double ttim, utim, stim, tott;
  struct tms runtimes;
  subharminfo *subharminf;
  accelobs obs;
  infodata idata;
  GSList *cands=NULL;
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
    int numcands;
    GSList *listptr;
    accelcand *cand;
    fourierprops *props;

    /* Now optimize each candidate and its harmonics */
    
    numcands = g_slist_length(cands);
    listptr = cands;
    for (ii=0; ii<numcands; ii++){
      print_percent_complete(ii, numcands);
      cand = (accelcand *)(listptr->data);
      optimize_accelcand(cand, &obs);
      listptr = listptr->next;
    }
    print_percent_complete(ii, numcands);
  
    /* Sort the candidates according to the optimized sigmas */

    cands = sort_accelcands(cands);

    /* Calculate the properties of the fundamentals */

    props = (fourierprops *)malloc(sizeof(fourierprops) * numcands);
    listptr = cands;
    for (ii=0; ii<numcands; ii++){
      cand = (accelcand *)(listptr->data);
      calc_props(cand->derivs[0], cand->hirs[0], cand->hizs[0], 
		 0.0, props + ii);
      listptr = listptr->next;
    }

    /* Write the fundamentals to the output text file */

    output_fundamentals(props, cands, &obs, &idata);
 
    /* Write the harmonics to the output text file */

    output_harmonics(props, cands, &obs, &idata);

    /* Write the fundamental fourierprops to the cand file */

    obs.workfile = chkfopen(obs.candnm, "wb");
    chkfwrite(props, sizeof(fourierprops), numcands, obs.workfile);
    fclose(obs.workfile);
    free(props);
  }

  /* Finish up */

  printf("Done.\n\n");
  printf("Searched the following approx numbers of independent points:\n");
  printf("  %d  harmonic:  %9lld\n", 1, obs.numindep[ii]);
  for (ii=1; ii<obs.numharm; ii++)
    printf("  %d harmonics:  %9lld\n", ii+1, obs.numindep[ii]);
  
  printf("\nTiming summary:\n");
  tott = times(&runtimes) / (double) CLK_TCK - tott;
  utim = runtimes.tms_utime / (double) CLK_TCK;
  stim = runtimes.tms_stime / (double) CLK_TCK;
  ttim = utim + stim;
  printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n", \
	 ttim, utim, stim);
  printf("  Total time: %.3f sec\n\n", tott);

  printf("Candidates in binary format are stored in '%s'.\n", obs.candnm);
  printf("Candidates in a text format are stored in '%s'.\n", obs.accelnm);

  free_accelobs(&obs);
  g_slist_foreach(cands, free_accelcand, NULL);
  g_slist_free(cands);
  return (0);
}
