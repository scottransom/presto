#include "accel.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


int calc_required_z(int numharm, int harmnum, double zfull);
double calc_required_r(int numharm, int harmnum, double rfull);
int index_from_r(double r, double lor);
int index_from_z(double z, double loz);


static void print_percent_complete(int current, int number, 
				   char *what, int reset)
{
  static int newper=0, oldper=-1;

  if (reset){
    oldper = -1;
    newper = 0;
  } else {
    newper = (int) (current / (float)(number) * 100.0);
    if (newper < 0) newper = 0;
    if (newper > 100) newper = 100;
    if (newper > oldper) {
      printf("\rAmount of %s complete = %3d%%", what, newper);
      fflush(stdout);
      oldper = newper;
    }
  }
}


int main(int argc, char *argv[])
{
  int ii;
  double ttim, utim, stim, tott;
  struct tms runtimes;
  subharminfo **subharminfs;
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

  /* Create the accelobs structure */
  
  create_accelobs(&obs, &idata, cmd);
  printf("Searching with up to %d harmonics summed:\n", obs.numharm);
  printf("  f = %.1f to %.1f\n", obs.rlo/obs.T, obs.rhi/obs.T);
  printf("  r = %.1f to %.1f\n", obs.rlo, obs.rhi);
  printf("  z = %.1f to %.1f\n\n", obs.zlo, obs.zhi);

  /* Generate the correlation kernels */
  
  printf("Generating correlation kernels:\n");
  subharminfs = create_subharminfos(obs.numharm, (int) obs.zhi);
  printf("Done generating kernels.\n\n");
  printf("Starting the search.\n");
  printf("  Working candidates in a test format are in '%s'.\n\n", 
	 obs.workfilenm);
  
  /* Start the main search loop */
  
  {
    double startr=obs.rlo, lastr=0, nextr=0;
    ffdotpows *fundamental;
    
    while (startr < obs.highestbin){  /* Search the fundamental */
      print_percent_complete(startr-obs.rlo, 
			     obs.highestbin-obs.rlo, "search", 0);
      nextr = startr + ACCEL_USELEN * ACCEL_DR;
      lastr = nextr - ACCEL_DR;
      fundamental = subharm_ffdot_plane(1, 1, startr, lastr, 
					&subharminfs[1][1], &obs);
      cands = search_ffdotpows(fundamental, 1, &obs, cands);
      
      if (obs.numharm > 1){   /* Search the subharmonics */
	ffdotpows *fundamental_copy, *subharmonic;
	
	fundamental_copy = copy_ffdotpows(fundamental);

	/* Search the 1/2 subharmonic (i.e. 1/2 fundamental) */

	subharmonic = subharm_ffdot_plane(2, 1, startr, lastr, 
					  &subharminfs[2][1], &obs);
	add_ffdotpows(fundamental, subharmonic, 2, 1);
	free_ffdotpows(subharmonic);
	cands = search_ffdotpows(fundamental, 2, &obs, cands);

	/* Search the 1/4 subharmonic by building on the 1/2 */

	if (obs.numharm >= 4){
	  subharmonic = subharm_ffdot_plane(4, 1, startr, lastr, 
					    &subharminfs[4][1], &obs);
	  add_ffdotpows(fundamental, subharmonic, 4, 1);
	  free_ffdotpows(subharmonic);
	  subharmonic = subharm_ffdot_plane(4, 3, startr, lastr, 
					    &subharminfs[4][3], &obs);
	  add_ffdotpows(fundamental, subharmonic, 4, 3);
	  free_ffdotpows(subharmonic);
	  cands = search_ffdotpows(fundamental, 4, &obs, cands);
	}
	
	/* Search the 1/3 subharmonic (work from scratch) */

	if (obs.numharm >= 3){
	  free_ffdotpows(fundamental);
	  fundamental = copy_ffdotpows(fundamental_copy);
	  subharmonic = subharm_ffdot_plane(3, 1, startr, lastr, 
					    &subharminfs[3][1], &obs);
	  add_ffdotpows(fundamental, subharmonic, 3, 1);
	  free_ffdotpows(subharmonic);
	  subharmonic = subharm_ffdot_plane(3, 2, startr, lastr, 
					    &subharminfs[3][2], &obs);
	  add_ffdotpows(fundamental, subharmonic, 3, 2);
	  free_ffdotpows(subharmonic);
	  cands = search_ffdotpows(fundamental, 3, &obs, cands);
	}
	free_ffdotpows(fundamental_copy);
      }
      free_ffdotpows(fundamental);
      startr = nextr;
    }
    print_percent_complete(startr-obs.rlo, 
			   obs.highestbin-obs.rlo, "search", 0);
  }

  printf("\n\nDone searching.  Now optimizing each candidate.\n\n");
  free_subharminfos(obs.numharm, subharminfs);

  {
    int numcands;
    GSList *listptr;
    accelcand *cand;
    fourierprops *props;

    /* Now optimize each candidate and its harmonics */
    
    numcands = g_slist_length(cands);

    if (numcands){
      listptr = cands;
      print_percent_complete(0, 0, NULL, 1);
      for (ii=0; ii<numcands; ii++){
	print_percent_complete(ii, numcands, "optimization", 0);
	cand = (accelcand *)(listptr->data);
	optimize_accelcand(cand, &obs);
	listptr = listptr->next;
      }
      print_percent_complete(ii, numcands, "optimization", 0);
  
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
      
      output_harmonics(cands, &obs);
      
      /* Write the fundamental fourierprops to the cand file */
      
      obs.workfile = chkfopen(obs.candnm, "wb");
      chkfwrite(props, sizeof(fourierprops), numcands, obs.workfile);
      fclose(obs.workfile);
      free(props);
      printf("\n\n");
    } else {
      printf("No candidates above sigma = %.2f were found.\n\n", 
	     obs.sigma);
    }
  }

  /* Finish up */

  printf("Searched the following approx numbers of independent points:\n");
  printf("  %d harmonic:   %9lld\n", 1, obs.numindep[0]);
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

  printf("Final candidates in binary format are in '%s'.\n", obs.candnm);
  printf("Final Candidates in a text format are in '%s'.\n\n", obs.accelnm);

  free_accelobs(&obs);
  g_slist_foreach(cands, free_accelcand, NULL);
  g_slist_free(cands);
  return (0);
}
