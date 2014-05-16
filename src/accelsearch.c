#include "accel.h"

/*#undef USEMMAP*/

#ifdef USEMMAP
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

typedef struct accel_cand_gpu{
	float pow;					/*pow of selected candidate*/
	int		nof_cand;			/*number of candidates in sub_array/plane */
	int		z_ind;				/*z_index of the selected candidate*/
	int		r_ind;				/*r_index of the selected candidate*/
}	accel_cand_gpu;

extern float calc_median_powers(fcomplex * amplitudes, int numamps);
extern void zapbirds(double lobin, double hibin, FILE * fftfile, fcomplex * fft);

/**************************** Functions on GPU ****************************/
#ifdef USECUDA
extern fcomplex* prep_data_on_gpu(subharminfo **subharminfs, int numharmstages);
extern fcomplex* prep_result_on_gpu(subharminfo **subharminfs, int numharmstages);

extern fcomplex* cp_kernel_array_to_gpu(subharminfo **subharminfs, int numharmstages, int **offset_array);
extern fcomplex* cp_input_to_gpu(fcomplex * input_vect_on_cpu, long long numbins, long long N);

extern unsigned short* prep_rz_inds_on_gpu(int size_inds);

extern float* prep_float_vect_on_gpu(int size);

extern accel_cand_gpu* prep_cand_array(int size);

extern void cudaFree_kernel_vect(fcomplex *in);

extern void select_cuda_dev(int cuda_inds);

extern int  search_ffdotpows_gpu(float powcut, float *d_fundamental, accel_cand_gpu * cand_array_search_gpu, accel_cand_gpu * cand_array_sort_gpu, int numzs, int numrs, accel_cand_gpu *cand_gpu_cpu);

extern GSList* search_ffdotpows_sort_gpu_result(ffdotpows * ffdot, int numharm,
                         accelobs * obs, GSList * cands, accel_cand_gpu *cand_gpu_cpu, int nof_cand, int numzs, int numrs);

extern void subharm_ffdot_plane_gpu(int numharm, int harmnum,
                               double fullrlo, double fullrhi,
                               subharminfo * shi, accelobs * obs,
                               int stage, int **offset_array,
                               fcomplex *kernel_vect_on_gpu, fcomplex *d_data, fcomplex *d_result,
                               ffdotpows *fundamental,
                               float *d_fundamental,
                               unsigned short *d_zinds, unsigned short *d_rinds
                               );

extern ffdotpows *init_fundamental_ffdot(int numharm, int harmnum,
                               double fullrlo, double fullrhi,
                               subharminfo * shi, accelobs * obs,
                               int stage, int **offset_array,
                               fcomplex *kernel_vect_on_gpu, fcomplex *d_data, fcomplex *d_result
                               );
                               
extern void init_cuFFT_plans(subharminfo **subharminfs, int numharmstages);                               
#endif
/**************************** Functions on GPU, End****************************/

static void print_percent_complete(int current, int number, char *what, int reset)
{
   static int newper = 0, oldper = -1;

   if (reset) {
      oldper = -1;
      newper = 0;
   } else {
      newper = (int) (current / (float) (number) * 100.0);
      if (newper < 0)
         newper = 0;
      if (newper > 100)
         newper = 100;
      if (newper > oldper) {
         printf("\rAmount of %s complete = %3d%%", what, newper);
         fflush(stdout);
         oldper = newper;
      }
   }
}

/* Return x such that 2**x = n */
static inline int twon_to_index(int n)
{
   int x = 0;

   while (n > 1) {
      n >>= 1;
      x++;
   }
   return x;
}

int main(int argc, char *argv[])
{
   int ii, jj;
   double ttim, utim, stim, tott;
   struct tms runtimes;
   subharminfo **subharminfs;
   accelobs obs;
   infodata idata;
   GSList *cands = NULL;
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
   printf("    Fourier-Domain Acceleration Search Routine\n");
   printf("               by Scott M. Ransom\n\n");
   printf("            GPU version by Jintao Luo\n\n");

   /* Create the accelobs structure */
   create_accelobs(&obs, &idata, cmd, 1);

   /* Zap birdies if requested and if in memory */
   if (cmd->zaplistP && !obs.mmap_file && obs.fft) {
      int numbirds;
      double *bird_lobins, *bird_hibins, hibin;

      /* Read the Standard bird list */
      numbirds = get_birdies(cmd->zaplist, obs.T, cmd->baryv,
                             &bird_lobins, &bird_hibins);

      /* Zap the birdies */
      printf("Zapping them using a barycentric velocity of %.5gc.\n\n", cmd->baryv);
      hibin = obs.N / 2;
      for (ii = 0; ii < numbirds; ii++) {
         if (bird_lobins[ii] >= hibin)
            break;
         if (bird_hibins[ii] >= hibin)
            bird_hibins[ii] = hibin - 1;
         zapbirds(bird_lobins[ii], bird_hibins[ii], NULL, obs.fft);
      }

      vect_free(bird_lobins);
      vect_free(bird_hibins);
   }

   printf("Searching with up to %d harmonics summed:\n",
          1 << (obs.numharmstages - 1));
   printf("  f = %.1f to %.1f Hz\n", obs.rlo / obs.T, obs.rhi / obs.T);
   printf("  r = %.1f to %.1f Fourier bins\n", obs.rlo, obs.rhi);
   printf("  z = %.1f to %.1f Fourier bins drifted\n\n", obs.zlo, obs.zhi);

   /* Generate the correlation kernels */

   printf("Generating correlation kernels:\n");
   subharminfs = create_subharminfos(&obs);
   printf("Done generating kernels.\n\n");
   printf("Starting the search.\n");
   /* Don't use the *.txtcand files on short in-memory searches */
   if (!obs.dat_input) {
      printf("  Working candidates in a test format are in '%s'.\n\n",
             obs.workfilenm);
   }

	int use_cuda_flag = 0;
	#ifdef USECUDA
	if(cmd->cudaP == 1)	
	{
		use_cuda_flag = 1 ;
	}
	#endif

#ifdef USECUDA
if(cmd->cudaP == 1){/****************************************call GPU to run the search***************************************************/
//if( use_cuda_flag == 1 ){
	
	printf("\nTry to use a GPU to run the search\n");
	
	#ifdef USECUDA
	printf( "\nUSECUDA\n" );
	#endif
	
	select_cuda_dev(cmd->cuda);
	
	//prepare d_data and d_result on GPU	
	fcomplex *d_data;
	fcomplex *d_result;	
	d_data = prep_data_on_gpu(subharminfs, obs.numharmstages);
	d_result = prep_result_on_gpu(subharminfs, obs.numharmstages);
	
	//prepare the fundamental on GPU
	float *d_fundamental;	
	d_fundamental = prep_float_vect_on_gpu(subharminfs[0][0].numkern * subharminfs[0][0].kern[0].fftlen );
	
	//inds lookup tables for add_ffdot_pow
	unsigned short *d_zinds, *d_rinds;
	d_zinds = prep_rz_inds_on_gpu(subharminfs[0][0].numkern);
	d_rinds = prep_rz_inds_on_gpu(subharminfs[0][0].kern[0].fftlen);//In fact fftlen >= size_d_rinds, but that's ok to have a bit waste
	
	//Generate the kernel offset_array	for stages
	int **offset_array;
	fcomplex *kernel_vect_on_gpu;
	offset_array = (int **)malloc(obs.numharmstages * sizeof(int *));
	offset_array[0] = (int *)malloc( 1 * sizeof(int) );	 
	for(ii=1; ii<obs.numharmstages; ii++){	
			jj = 1 << ii;		
			offset_array[ii] = (int *)malloc( jj * sizeof(int) );		
	}	

	//copy all kernels into the GPU memory
	kernel_vect_on_gpu = cp_kernel_array_to_gpu(subharminfs,                 obs.numharmstages, offset_array);

	//prepare arrays for ffdotpow_search on GPU
	accel_cand_gpu * cand_array_search_gpu;
	accel_cand_gpu * cand_array_sort_gpu;	
	cand_array_search_gpu = prep_cand_array(subharminfs[0][0].numkern * subharminfs[0][0].kern[0].fftlen); 
	cand_array_sort_gpu		= prep_cand_array(subharminfs[0][0].numkern * subharminfs[0][0].kern[0].fftlen); 
	
	//vector on CPU to store the GPU ffdot_search result
	accel_cand_gpu *cand_gpu_cpu ;	
	cand_gpu_cpu = (accel_cand_gpu *) malloc(sizeof(accel_cand_gpu) * subharminfs[0][0].numkern * subharminfs[0][0].kern[0].fftlen );
	
	//initialize cuFFT plans on GPU
	init_cuFFT_plans( subharminfs, obs.numharmstages );

   /* Start the main search loop */
   printf("\n\n/* Start the main search loop */\n\n");

   {
      double startr = obs.rlo, lastr = 0, nextr = 0;	
			ffdotpows *my_fundamental;			
			int nof_cand = 0 ;
			
			my_fundamental= init_fundamental_ffdot(1, 1, startr, lastr,
                               &subharminfs[0][0], &obs,
                               0, offset_array,
                               kernel_vect_on_gpu, d_data, d_result
                               );
			
      while (startr + ACCEL_USELEN * ACCEL_DR < obs.highestbin) {
         /* Search the fundamental */
         print_percent_complete(startr - obs.rlo,
                                obs.highestbin - obs.rlo, "search", 0);
         nextr = startr + ACCEL_USELEN * ACCEL_DR;
         lastr = nextr - ACCEL_DR;

						subharm_ffdot_plane_gpu(1, 1, startr, lastr,
                               &subharminfs[0][0], &obs,
                               0, offset_array,
                               kernel_vect_on_gpu, d_data, d_result, 
                               my_fundamental,
                               d_fundamental,
                               d_zinds, d_rinds);                                               
         nof_cand = search_ffdotpows_gpu(obs.powcut[twon_to_index(1)], d_fundamental, cand_array_search_gpu, cand_array_sort_gpu, my_fundamental->numzs, my_fundamental->numrs, cand_gpu_cpu);
         cands = search_ffdotpows_sort_gpu_result(my_fundamental, 1, &obs, cands, cand_gpu_cpu, nof_cand, my_fundamental->numzs, my_fundamental->numrs);


         if (obs.numharmstages > 1) {   /* Search the subharmonics */
            int stage, harmtosum, harm;

            for (stage = 1; stage < obs.numharmstages; stage++) { 
               harmtosum = 1 << stage;
               for (harm = 1; harm < harmtosum; harm += 2) {                           
                  subharm_ffdot_plane_gpu(harmtosum, harm, startr, lastr,
                               &subharminfs[stage][harm - 1], &obs,
                               stage, offset_array,
                               kernel_vect_on_gpu, d_data, d_result,
                               my_fundamental,
                               d_fundamental,
                               d_zinds, d_rinds);                               
               }
               nof_cand = search_ffdotpows_gpu(obs.powcut[twon_to_index(harmtosum)], d_fundamental, cand_array_search_gpu, cand_array_sort_gpu, my_fundamental->numzs, my_fundamental->numrs, cand_gpu_cpu);
               cands = search_ffdotpows_sort_gpu_result(my_fundamental, harmtosum, &obs, cands, cand_gpu_cpu, nof_cand, my_fundamental->numzs, my_fundamental->numrs);
            }
         }
         startr = nextr;
      }
      free_ffdotpows(my_fundamental);            
      print_percent_complete(obs.highestbin - obs.rlo,
                             obs.highestbin - obs.rlo, "search", 0);
   }

   printf("\n\nDone searching.  Now optimizing each candidate.\n\n");

	 destroy_cuFFT_plans( subharminfs, obs.numharmstages ); //destroy cuFFT plans
	 
	 printf("\ndebug:.......................\n");
	 
   //free_subharminfos(obs.numharmstages, subharminfs); //from presto1
   free_subharminfos(&obs, subharminfs); //from presto2

   {                            /* Candidate list trimming and optimization */
      int numcands;
      GSList *listptr;
      accelcand *cand;
      fourierprops *props;


      numcands = g_slist_length(cands);

      if (numcands) {

         /* Sort the candidates according to the optimized sigmas */

         cands = sort_accelcands(cands);

         /* Eliminate (most of) the harmonically related candidates */
         if ((cmd->numharm > 1) && !(cmd->noharmremoveP))
             eliminate_harmonics(cands, &numcands);

         /* Now optimize each candidate and its harmonics */

         print_percent_complete(0, 0, NULL, 1);
         listptr = cands;
         for (ii = 0; ii < numcands; ii++) {
            print_percent_complete(ii, numcands, "optimization", 0);
            cand = (accelcand *) (listptr->data);
            optimize_accelcand(cand, &obs);
            listptr = listptr->next;
         }
         print_percent_complete(ii, numcands, "optimization", 0);

         /* Calculate the properties of the fundamentals */

         props = (fourierprops *) malloc(sizeof(fourierprops) * numcands);
         listptr = cands;
         for (ii = 0; ii < numcands; ii++) {
            cand = (accelcand *) (listptr->data);
            /* In case the fundamental harmonic is not significant,  */
            /* send the originally determined r and z from the       */
            /* harmonic sum in the search.  Note that the derivs are */
            /* not used for the computations with the fundamental.   */
            calc_props(cand->derivs[0], cand->r, cand->z, 0.0, props + ii);
            /* Override the error estimates based on power */
            props[ii].rerr = (float) (ACCEL_DR) / cand->numharm;
            props[ii].zerr = (float) (ACCEL_DZ) / cand->numharm;
            listptr = listptr->next;
         }

         /* Write the fundamentals to the output text file */

         output_fundamentals(props, cands, &obs, &idata);

         /* Write the harmonics to the output text file */

         output_harmonics(cands, &obs, &idata);

         /* Write the fundamental fourierprops to the cand file */

         obs.workfile = chkfopen(obs.candnm, "wb");
         chkfwrite(props, sizeof(fourierprops), numcands, obs.workfile);
         fclose(obs.workfile);
         free(props);
         printf("\n\n");
      } else {
         printf("No candidates above sigma = %.2f were found.\n\n", obs.sigma);
      }
   }

   /* Finish up */

   printf("Searched the following approx numbers of independent points:\n");
   printf("  %d harmonic:   %9lld\n", 1, obs.numindep[0]);
   for (ii = 1; ii < obs.numharmstages; ii++)
      printf("  %d harmonics:  %9lld\n", 1 << ii, obs.numindep[ii]);

   printf("\nTiming summary:\n");
   tott = times(&runtimes) / (double) CLK_TCK - tott;
   utim = runtimes.tms_utime / (double) CLK_TCK;
   stim = runtimes.tms_stime / (double) CLK_TCK;
   ttim = utim + stim;
   printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n",
          ttim, utim, stim);
   printf("  Total time: %.3f sec\n\n", tott);

   printf("Final candidates in binary format are in '%s'.\n", obs.candnm);
   printf("Final Candidates in a text format are in '%s'.\n\n", obs.accelnm);

   free_accelobs(&obs);
   g_slist_foreach(cands, free_accelcand, NULL);
   g_slist_free(cands);
   
  free(offset_array);    
  cudaFree_kernel_vect(d_data);
  cudaFree_kernel_vect(d_result);
  cudaFree_kernel_vect(d_fundamental);
	cudaFree_kernel_vect(kernel_vect_on_gpu);	
	cudaFree_kernel_vect(d_zinds);
	cudaFree_kernel_vect(d_rinds);	
	cudaFree_kernel_vect(cand_array_search_gpu);
	cudaFree_kernel_vect(cand_array_sort_gpu);
	//cudaFree_kernel_vect(input_vect_on_gpu);
  
  
  free(cand_gpu_cpu);
   
   return (0);
   
	
}
#endif

/****************************************call CPU to run the search***************************************************/
//else{
if( use_cuda_flag == 0 ){
	printf("\nUse CPU to run the search\n");
   /* Start the main search loop */
   {
      double startr = obs.rlo, lastr = 0, nextr = 0;
      ffdotpows *fundamental;

      while (startr + ACCEL_USELEN * ACCEL_DR < obs.highestbin) {
         /* Search the fundamental */
         print_percent_complete(startr - obs.rlo,
                                obs.highestbin - obs.rlo, "search", 0);
         nextr = startr + ACCEL_USELEN * ACCEL_DR;
         lastr = nextr - ACCEL_DR;
         fundamental = subharm_ffdot_plane(1, 1, startr, lastr,
                                           &subharminfs[0][0], &obs);
         cands = search_ffdotpows(fundamental, 1, &obs, cands);

         if (obs.numharmstages > 1) {   /* Search the subharmonics */
            int stage, harmtosum, harm;
            ffdotpows *subharmonic;

            // Copy the fundamental's ffdot plane to the full in-core one
            if (obs.inmem){
                if (cmd->otheroptP)
                    fund_to_ffdotplane_trans(fundamental, &obs);
                else
                    fund_to_ffdotplane(fundamental, &obs);
            }
            for (stage = 1; stage < obs.numharmstages; stage++) {
               harmtosum = 1 << stage;
               for (harm = 1; harm < harmtosum; harm += 2) {
                   if (obs.inmem) {
                       if (cmd->otheroptP)
                           inmem_add_ffdotpows_trans(fundamental, &obs, harmtosum, harm);
                       else
                           inmem_add_ffdotpows(fundamental, &obs, harmtosum, harm);
                   } else {
                       subharmonic = subharm_ffdot_plane(harmtosum, harm, startr, lastr,
                                                         &subharminfs[stage][harm - 1],
                                                         &obs);
                       if (cmd->otheroptP)
                           add_ffdotpows_ptrs(fundamental, subharmonic, harmtosum, harm);
                       else
                           add_ffdotpows(fundamental, subharmonic, harmtosum, harm);
                       free_ffdotpows(subharmonic);
                   }
               }
               cands = search_ffdotpows(fundamental, harmtosum, &obs, cands);
            }
         }
         free_ffdotpows(fundamental);
         startr = nextr;
      }
      print_percent_complete(obs.highestbin - obs.rlo,
                             obs.highestbin - obs.rlo, "search", 0);
   }

   printf("\n\nDone searching.  Now optimizing each candidate.\n\n");
   free_subharminfos(&obs, subharminfs);

   {                            /* Candidate list trimming and optimization */
      int numcands;
      GSList *listptr;
      accelcand *cand;
      fourierprops *props;


      numcands = g_slist_length(cands);

      if (numcands) {

         /* Sort the candidates according to the optimized sigmas */

         cands = sort_accelcands(cands);

         /* Eliminate (most of) the harmonically related candidates */
         if ((cmd->numharm > 1) && !(cmd->noharmremoveP))
             eliminate_harmonics(cands, &numcands);

         /* Now optimize each candidate and its harmonics */

         print_percent_complete(0, 0, NULL, 1);
         listptr = cands;
         for (ii = 0; ii < numcands; ii++) {
            print_percent_complete(ii, numcands, "optimization", 0);
            cand = (accelcand *) (listptr->data);
            optimize_accelcand(cand, &obs);
            listptr = listptr->next;
         }
         print_percent_complete(ii, numcands, "optimization", 0);

         /* Calculate the properties of the fundamentals */

         props = (fourierprops *) malloc(sizeof(fourierprops) * numcands);
         listptr = cands;
         for (ii = 0; ii < numcands; ii++) {
            cand = (accelcand *) (listptr->data);
            /* In case the fundamental harmonic is not significant,  */
            /* send the originally determined r and z from the       */
            /* harmonic sum in the search.  Note that the derivs are */
            /* not used for the computations with the fundamental.   */
            calc_props(cand->derivs[0], cand->r, cand->z, 0.0, props + ii);
            /* Override the error estimates based on power */
            props[ii].rerr = (float) (ACCEL_DR) / cand->numharm;
            props[ii].zerr = (float) (ACCEL_DZ) / cand->numharm;
            listptr = listptr->next;
         }

         /* Write the fundamentals to the output text file */

         output_fundamentals(props, cands, &obs, &idata);

         /* Write the harmonics to the output text file */

         output_harmonics(cands, &obs, &idata);

         /* Write the fundamental fourierprops to the cand file */

         obs.workfile = chkfopen(obs.candnm, "wb");
         chkfwrite(props, sizeof(fourierprops), numcands, obs.workfile);
         fclose(obs.workfile);
         free(props);
         printf("\n\n");
      } else {
         printf("No candidates above sigma = %.2f were found.\n\n", obs.sigma);
      }
   }

   /* Finish up */

   printf("Searched the following approx numbers of independent points:\n");
   printf("  %d harmonic:   %9lld\n", 1, obs.numindep[0]);
   for (ii = 1; ii < obs.numharmstages; ii++)
      printf("  %d harmonics:  %9lld\n", 1 << ii, obs.numindep[ii]);

   printf("\nTiming summary:\n");
   tott = times(&runtimes) / (double) CLK_TCK - tott;
   utim = runtimes.tms_utime / (double) CLK_TCK;
   stim = runtimes.tms_stime / (double) CLK_TCK;
   ttim = utim + stim;
   printf("    CPU time: %.3f sec (User: %.3f sec, System: %.3f sec)\n",
          ttim, utim, stim);
   printf("  Total time: %.3f sec\n\n", tott);

   printf("Final candidates in binary format are in '%s'.\n", obs.candnm);
   printf("Final Candidates in a text format are in '%s'.\n\n", obs.accelnm);

   free_accelobs(&obs);
   g_slist_foreach(cands, free_accelcand, NULL);
   g_slist_free(cands);
   return (0);
}

}
