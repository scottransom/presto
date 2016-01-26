/*
   It:
	1. Accepts a stream of single byte input samples
	   ordered as f_1 f_2 ... f_{chn_per_sample} in
  	   frequency for each time sample.
	2. Analyzes separately a specificed number (n_subbands)
	   of subbands.  I.e. it dedisperses each subband separately 
	   using a single value of DM, a list of DMs, or using the
           Taylor tree algorithm. 
	3. The environment variable dm_dir is used to define
	   the path for the output data (dedispersed time series).
	4. Output data may be written as unsigned bytes, short ints,
           ints, or floats, selectable using the -o option.
           Default is bytes (== -o c).

   JMC 10-26 Sep 1994 
   SCL 1994-1995 changes for SGI and improvements on output file open/close
   JMC 10-11 July 1996 to run on Cornell Sun's
   JMC & MAM: 5 Aug 1996: rework the way 'takeback' samples are
             saved and to correct errors in unflip procedure and the
             number of takeback samples.

	     Also, changed the brute force (-t) option so that it
	     reads in a list of DM values from file dmlist, rather
             than computing the DMs. 

             8 Aug 1996: naming of DM files: now all DM files are of
	     the form DMXXX_Y with Y = subband #, even for the case
	     where there are not multiple subbands.	
	    
             6 Sep 1996: rename lenfft to ntimesamples to signify that
             the number of samples dedispersed is not linked to 
             a 2^n length FFT.   Also, eliminate loopmax as an input variable,
             as this is now calculated. 
*/

#include 	<stdio.h>
#include 	<fcntl.h>
#include 	<malloc.h>
#include        <math.h>
#include        <stdlib.h>

#define  STDOUT 1
#define  FOUT_MAX 256
#define FOPEN_MAX 400
/* #define  ndim   34000 */
#define ndim 1644
#define  mdim  512
#define  ldim   2
#define  ndmmax 512
#define  mod(a,b)  ( (a) - (b)*(int)( (a) / (b) ) )

void processargs();
void dedisperse_tree();
void calc_delays_dmlist();
void dedisperse_dmlist();
void tree_stage();
void dedisperse_single_dm();
void calc_delays();

main(argc, argv)
int argc;
char **argv;
{
/*
 *	dedisperse data stream
 *	
 *       ddisp_p 
 *          -b bins_input	chunk size in time samples <,= ndim 
 *          -c chn_per_sample	# frequency channels per time sample
 *          -d ntimesamples     data span length to dedisperse
 *          -f frequency 	center frequency in MHz of lowest freq channel
 *          -i invert           multiplier of data (to invert) 
 *          -m dm		dedispersion with single dm value
 *          -n n_subbands	# subbands to dedisperse separately
 *          -o outtype          c,s,i,f for byte,short,int, float
 *          -q (no arg)         (-> quiet mode)
 *          -s sample_time  	sample interval in seconds
 *          -u (no arg)         unflip the frequency channels 
 *          -t dedisp_type      1=tree algorithm, > 1: brute force dedisp
 *          -w bandwidth  	total bandwidth in MHz
 *          -z shift            dedisperse shift parameter in bins
 *
 * Further explanations of options:
 * [all are set in ddisp_defs]
 *
 *      -b   data is read in chunks of bins_input time samples <,= ndim
 *      -c   number of channels in input 
 *      -d   data span length for fft (determines length of dedispersed
 *           time series.
 *      -f   frequency in MHz
 *      -i   specify a multiplier of the data (for inversion purposes)
 *      -l   number of chunks of data (bins_input long) to be read
 *      -m   to specify a value of DM to dedisperse separate sidebands with
 *           and  to dedisperse exactly rather than use the tree algorithm
 *      -n   number of subbands to dedisperse
 *      -o   outtype: c,s,i,f for byte,short,int, float
 *      -q   if used, prompts will be shown
 *      -s   time for 1 sample (def. = 1.0 seconds).
 *           A sample consists of -c chn_per_sample numbers of the 
 *           specified data type 
 *           [e.g. 
 *			32 for AO filter bank; 
 *                 	64 for AO correlator;
 *       	    64-256 for Parkes filter banks 
 *           ]
 *      -u   unflips the frequency axis when needed by LO scheme used
 *      -t   method for dedispersion  1=tree algorithm, 
 *                                    >1: using brute force dedisp
 *      -w   total bandwidth in MHz
 *      -z   shift = 1 if largest DM used corresponds to a shift 
 *           of one time sample between adjacent frequency channels;
 *           larger DM's can be searched using shift > 1
 *

 */


#include "dedisp_defs"
    char dm_names[ndmmax][ldim][80];
    int jj, kk, iii;
    int loop_print;
    int nskip, newinit, max_delay;
    float fref, fch1, chbw;
    char *dm_dir;
/* char dm_name_environ[200]; */
    char dm_name_environ[BUFSIZ];
    char tree_dms[80], opt_dms[80], save_dms[80];
    fprintf(stderr, "\n\nstarting dedisp\n\n");
    fflush(stderr);

    ddisp_diag = fopen("ddisp_diag", "w");
    delay_list = fopen("delays", "w");


/* GET THE PARAMETERS */
    unflip = 0;

    processargs(argc, argv, &chn_per_sample, &bins_input,
                &sample_time, &shift, &frequency, &bandwidth,
                &invert, &quiet, &n_subbands, &ntimesamples, &unflip,
                &dm, &do_single_dm, &outtype, &bytesout_datatype, &dedisp_type);


/* check that buffer block length is smaller than declared dimension, ndim */
    if (bins_input > ndim) {
        fprintf(stderr,
                "\n requested buffer length is too large: \n bins_input = %d \n ndim = %d\n",
                bins_input, ndim);
        exit(-1);
    }

    fprintf(stderr, "\n\nstarting dedisp using dedisp_type = %d\n\n", dedisp_type);


/*  INITIALIZE DISPERSION DELAYS */

    if (dedisp_type > 1) {
    } else {
        ndm = chn_per_sample / n_subbands;
        sprintf(tree_dms, "dmlist_tree");
        if ((dmlist_tree = fopen(tree_dms, "w")) == NULL) {
            fprintf(stderr, "cannot open dmlist file %30s\n", tree_dms);
            return 1;
        }
        fprintf(dmlist_tree, "%d\n", ndm);
        for (n = 0; n < ndm; n++) {
            dmvalues_tree[n] =
                n * sample_time * shift * frequency * frequency * frequency / 8.3e3 /
                (bandwidth / n_subbands);
            fprintf(dmlist_tree, "%f\n", dmvalues_tree[n]);
        }
    }

    fclose(dmlist_tree);


    sprintf(opt_dms, "dmlist_opt");
    if ((dmlist_opt = fopen(opt_dms, "r")) == NULL) {
        fprintf(stderr, "cannot open dmlist file %30s\n", opt_dms);
        return 1;
    }
    sprintf(save_dms, "dmlist");
    if ((dmlist = fopen(save_dms, "w")) == NULL) {
        fprintf(stderr, "cannot open dmlist file %30s\n", save_dms);
        return 1;
    }

    ndm_tree = 0;
    fscanf(dmlist_opt, "%d", &ndm_opt);
    fprintf(dmlist, "%d\n", ndm_opt);
    for (n = 0; n < ndm_opt; n++) {
        fscanf(dmlist_opt, "%f", &dm_opt);
        while (dmvalues_tree[ndm_tree] < dm_opt)
            ndm_tree = ndm_tree + 1;
        if (ndm_tree < ndm) {
            fprintf(dmlist, "%f\n", dmvalues_tree[ndm_tree]);
            ndm_save[n] = ndm_tree;
        }
    }

    fclose(dmlist_opt);
    fclose(dmlist);

    if (do_single_dm) {
        /* f_center = cf of entire bandwidth   */
        /* fch1_ns  = cf of ch1 for subband ns */
        /* fch1_0   = cf of ch1 for subband 0  */
        n_dm_channels_per_subband = 1;
        f_center = frequency + (bandwidth / 2.) * (1. - 1. / chn_per_sample);
        fch1_0 = frequency + 0.5 * bandwidth / chn_per_sample;
        fch1_0 = frequency;


/*
       fch1_ns = fch1_0 + (ns * bandwidth) / n_subbands;            
*/
        calc_delays(&delays, chn_per_sample, fch1_0, fch1_0,
                    bandwidth / chn_per_sample, dm, sample_time);

        /* fprintf(ddisp_diag, 
           "fref=%f, fch1=%f, channel bandwidth=%f, sample_time=%f\n",
           fch1_0, fch1_0,  bandwidth/chn_per_sample, sample_time);

           if (quiet == 1) { 
           for(n=0;n<chn_per_sample;n++) {
           fprintf(ddisp_diag, "channel no. %d  delay = %d\n", 
           n,delays[n][0]);
           }
           } */
    } else
        n_dm_channels_per_subband = ndm;

    /* fprintf(ddisp_diag, 
       "n_dm_channels_per_subband = %d\n", n_dm_channels_per_subband); */




/* ooo output setup stuff */

/* GET DIRECTORY FOR DM FILES FROM ENVIRONMENT VARIABLE */

    sprintf(dm_name_environ, "dm_dir");
    dm_dir = getenv(dm_name_environ);

    /*fprintf(ddisp_diag, "environment variable = %s\n", dm_name_environ);
       fprintf(ddisp_diag, "%s = %s\n", dm_name_environ, dm_dir);
       fprintf(ddisp_diag, "outtype = %c   bytesout_datatype = %d\n",
       outtype, bytesout_datatype);

       fflush(ddisp_diag); */

/* ooo end output stuff */

    /*fprintf(ddisp_diag, "do_single_dm = %d,  dm = %f\n", do_single_dm, dm); */


/* NUMBER OF CHANNELS IN EACH SUBBAND: */

    nchans_per_subband = chn_per_sample / n_subbands;

    printf("%d %d %d %d\n", chn_per_sample, bins_input, n_subbands,
           nchans_per_subband);

/* 
   For dedispersion using the tree algorithm 
   (i.e. -t 1 and no use of the -m option),
   the number of output files is chn_per_sample (if dedisp_type=1).

   For dedisp_type > 1, the number of output files is ndm as read in
   from the file ``dmlist''

   With the -m option, a single output channel is produced
   per subband that is dedispersed at the specified dm    

   there are chn_per_sample total output files:

        # output files = n_subbands * nchans_per_subband
        # trial DM's   = nchans_per_subband

   create output files with the following naming convention:
                 DMxxx_y where xxx = dispersion channel number
                                y = subband number

   If do_single_dm is true, the subbands are dedispersed according
   to the exact dispersion law with a single reference frequency, defined
   to be the center of the lowest channel (fch1_0).  A pulsar dedispersed
   with the proper DM will then show pulses that align in the time series
   for each subband
*/

/* ooo output stuff */

/* NAMING & OPENING OUTPUT DM FILES */

    for (j = 0; j < n_subbands; j++) {
        nsave = 0;
        for (k = 0; k < n_dm_channels_per_subband; k++) {
            if (k == ndm_save[nsave]) {
                if (nsave < 10)
                    sprintf(name, "%s/DM00%d_%d", dm_dir, nsave, j);
                if (nsave > 9 && nsave < 100)
                    sprintf(name, "%s/DM0%d_%d", dm_dir, nsave, j);
                if (nsave > 99)
                    sprintf(name, "%s/DM%d_%d", dm_dir, nsave, j);
                nsave = nsave + 1;
                fout[k][j] = fopen(name, "wb");
                sprintf(dm_names[k][j], "%s", name);
                fclose(fout[k][j]);
            }
        }
    }

    fflush(stderr);
/* ooo end output stuff */

/* tkbk_samples  = number of samples needed to allow for
                   dedispersing lags; this is determined
                   by the number of channels used in the
                   dedispersion, which is the number of
                   channels in a subband.                       

   For the dispersion tree mode tkbk_samples is determined by 
   the number of channels in a subband and the shift parameter.
 
   For -m mode with dedispersion at a specified dm, tkbk_samples
   is determined by the dispersion delay across the entire 
   bandpass (not just the subband bandpass).

   For the brute force dedispersion using a list of DMs, tkbk_samples
   is determined by the maximum dispersion delay across the entire
   passband for the maximum DM.
*/

    delta_nu = bandwidth * MHz;
    delta_nu_subband = delta_nu / n_subbands;
    f = frequency * MHz;
    delta_dm = 0.5 * dc * f * f * f * sample_time * shift / delta_nu;

    if (do_single_dm) {
/*   f_low = f - delta_nu/2.;
     f_up = f_low + delta_nu_subband;
     delta_t_dm = (dm/dc) * (1./(f_low*f_low) - 1./(f_up*f_up));
     delta_t_dm_samples = delta_t_dm / sample_time;
     tkbk_samples = delta_t_dm_samples + 0.5; }
*/
        tkbk_samples = -delays[chn_per_sample - 1][0];
    } else {
        if (dedisp_type > 1) {
            tkbk_samples = -delays[chn_per_sample - 1][ndm - 1];
        } else {
            tkbk_samples = nchans_per_subband * shift;
        }
    }

    fprintf(stderr, " %i %i \n", shift, tkbk_samples);

    if (tkbk_samples > bins_input) {
        fprintf(stderr, "TKBK_SAMPLES (%d) LARGER THAN BINS_INPUT (%d)!\n",
                tkbk_samples, bins_input);
        exit(-1);
    }

    if (dedisp_type > 1) {
        for (idm = 0; idm < ndm; idm++) {
            for (ns = 0; ns < n_subbands; ns++) {
                max_delay = -delays[(ns + 1) * nchans_per_subband - 1][idm];
                fprintf(delay_list, "%d %d %d\n", idm, ns, max_delay);
            }
        }
    }

/* ALLOCATE SPACE FOR INPUT BUFFER
   INPUT BUFFER: READ IN ONE FULL INTERVAL = bins_input TIME SAMPLES   
   tkbk_bytes = NUMBER OF TAKE BACK BYTES IN ALL THE CHANNELS ON INPUT 
*/
    tkbk_bytes_in = tkbk_samples * chn_per_sample;
    samples_inp_req = bins_input;
    bytes_per_sample = chn_per_sample * bytes_datatype;
    bytes_inp_req = samples_inp_req * bytes_per_sample;
    inp_buf = (unsigned char *) calloc((unsigned) bytes_inp_req, sizeof(char));
    if (inp_buf == NULL) {
        perror("ddisp: Allocating input buffer");
        exit(-1);
    }
    samples_all = bins_input;

/* Before starting loop. read in the take back samples that can be
   used for the first time through the loop: */

    bytes_inp = read_pipe(inp_buf, tkbk_bytes_in);

    if (bytes_inp < 0) {
        perror("dedisp: failed to get take back samples");
        exit(-1);
    }

    if (unflip == 1) {
        kstart = chn_per_sample - 1;    /* mam changed 10-22-96 */
        kinc = -1;
    } else {
        kstart = 0;
        kinc = 1;
    }

    tkbk_data = (float *) malloc(tkbk_samples * chn_per_sample * sizeof(float));

    for (i = 0; i < tkbk_samples; i++) {
        k = kstart + i * chn_per_sample;        /* mam changed 10-22-96 */
        for (j = 0; j < chn_per_sample; j++) {
            tkbk_data[i * chn_per_sample + j] = inp_buf[k];
            k += kinc;
        }
    }

    fprintf(stderr, " ok through loop!\n");
    fflush(stderr);

    bts_inp_req = bytes_inp_req - tkbk_bytes_in;
    bytes_inp = read_pipe(inp_buf, bts_inp_req);

    fprintf(stderr, "bytes_inp = %d\n", bytes_inp);
    fflush(stderr);

    for (ns = 0; ns < n_subbands; ns++) {       /* LOOP OVER SUBBANDS */
        nskip = ns * nchans_per_subband;
        mean[ns] = 0;
        for (i = samps_skip; i < samples_all - tkbk_samples; i++) {
            k = i * chn_per_sample;
            testddata = 0;
            for (j = 0; j < nchans_per_subband; j++) {
                testddata = testddata + inp_buf[k];
                k += 1;
            }
            mean[ns] = mean[ns] + testddata;
        }
        mean[ns] = mean[ns] / (samples_all - tkbk_samples - samps_skip);
    }

    for (ns = 0; ns < n_subbands; ns++) {       /* LOOP OVER SUBBANDS */
        nskip = ns * nchans_per_subband;
        rms[ns] = 0;
        for (i = samps_skip; i < samples_all - tkbk_samples; i++) {
            k = i * chn_per_sample;
            testddata = 0;
            for (j = 0; j < nchans_per_subband; j++) {
                testddata = testddata + inp_buf[k];
                k += 1;
            }
            rms[ns] = rms[ns] + (testddata - mean[ns]) * (testddata - mean[ns]);
        }
        rms[ns] = sqrt(rms[ns] / (samples_all - tkbk_samples - samps_skip - 1));
        min[ns] = mean[ns] - rms[ns] * rms_cutoff;
        max[ns] = mean[ns] + rms[ns] * rms_cutoff;
    }



/**********************  LOOP OVER BLOCK READS ************************/

    loop_max = ntimesamples / (bins_input - tkbk_samples) + 1;
    fprintf(stderr, "ntimesamples, loop_max = %d %d \n", ntimesamples, loop_max);
    fprintf(stderr, "bins_input, tkbk_samples = %d %d \n", bins_input, tkbk_samples);
    fflush(stderr);
/*
   if (loops_asked < loop_max)
     loop_max = loops_asked;
*/
    loop_print = 0.1 * loop_max;

    /*fprintf(ddisp_diag, "chn_per_sample = %d\n", chn_per_sample);
       fprintf(ddisp_diag, "nchans_per_subband = %d\n", nchans_per_subband);
       fprintf(ddisp_diag, "n_subbands = %d\n", n_subbands);
       fprintf(ddisp_diag, "bins_input = %d\n", bins_input);
       fprintf(ddisp_diag, "tkbk_samples = %d\n", tkbk_samples);
       fprintf(ddisp_diag, "tkbk_bytes_in = %d\n", tkbk_bytes_in);
       fprintf(ddisp_diag, "bytes_per_sample = %d\n", bytes_per_sample);
       fprintf(ddisp_diag, "bytes_inp_req = %d\n", bytes_inp_req);
       fprintf(ddisp_diag, "samples_all = %d\n", samples_all);
       fprintf(ddisp_diag, "ntimesamples= %d\n", ntimesamples);
       fprintf(ddisp_diag, "loop_max = %d\n", loop_max); */

    for (loop = 0; loop <= loop_max; loop++) {

/* iii data input section */
/* READ DATA INTO inp_buf: */

        bts_inp_req = bytes_inp_req - tkbk_bytes_in;
        /*if(quiet == 1) {
           fprintf(ddisp_diag, "bts_inp_req = %d\n", bts_inp_req); 
           fflush(ddisp_diag);
           } */

        if (loop != 0)
            bytes_inp = read_pipe(inp_buf, bts_inp_req);
        /*if(quiet == 1) {
           fprintf(ddisp_diag, "bytes_inp = %d\n", bytes_inp);
           fflush(ddisp_diag);
           } */

        if (bytes_inp == 0)
            goto done;          /* HIT EOF */
        if (bytes_inp < 0) {
            perror("ddisp: inputting data");
            exit(-1);
        }

/*
     if (quiet == 1) { 
       for (i=0; i<bytes_inp/chn_per_sample; i++) { 
	 for (j=0; j<chn_per_sample;j++) {
	   fprintf(ddisp_diag, "%c", inp_buf[i*chn_per_sample+j]);
	 }
	 fprintf(ddisp_diag, "\n");
       }
       fprintf(ddisp_diag, "done printing to ddisp_diag\n ");
     }
*/

/* iii end data input section */

/* samples_inp = NUMBER OF TIME SAMPLES INPUTTED 
   bytes_out_req = NUMBER OF BYTES TO BE WRITTEN FOR EACH DM CHANNEL 
*/

/*   if (loop == 0) 
       bytes_out_req = (bytesout_datatype*(bytes_inp - tkbk_bytes_in)) 
	 / chn_per_sample;
     else */
        bytes_out_req = (bytesout_datatype * bytes_inp) / chn_per_sample;

        /* if(quiet == 1) {
           fprintf(ddisp_diag, "bytes_out_req = %d\n", bytes_out_req);
           fflush(ddisp_diag);
           } */

        samples_inp = bytes_inp / bytes_per_sample;
        samples_out = bytes_out_req / bytesout_datatype;
        samples_done += samples_out;

        /*  if(quiet == 1) {
           fprintf(ddisp_diag, 
           "samples_inp, samples_out, samples_done, %d %d %d\n",
           samples_inp, samples_out, samples_done);

           fflush(ddisp_diag);
           } */

/*     printf("\t loop = %d \t samples_out = %d \t samples_done = %d\r", 
	     loop, samples_out, samples_done);
*/




/* LOOP OVER SUBBANDS; MOVE DATA FROM EACH SUBBAND INTO THE DEDISPERSING
   ARRAY (ddata), DEDISPERSE, AND WRITE OUT TO DISK.

   DATA IN INPUT ARRAY ARE ARRANGED AS f1,f2...fchn_per_sample
   (OR POSSIBLY IN FLIPPED ORDER)   
*/


        for (ns = 0; ns < n_subbands; ns++) {   /* LOOP OVER SUBBANDS */
            nskip = ns * nchans_per_subband;

            /* if(quiet == 1) {
               fprintf(ddisp_diag, "working on subband %d\n", ns);
               fflush(ddisp_diag);
               } */

            /* RESTORE RAW DATA TO FIRST tkbk_samples IN DDATA:
               SUBSTITUTE THOSE WITH LAST tkbk_samples IN DATA */

/*if(quiet == 1) {
	   fprintf(ddisp_diag, 
		   "restoring from tkbk %d samples  init = %d samples\n", 
		   tkbk_samples, init);
	   fflush(ddisp_diag);
         } */

            for (i = 0; i < tkbk_samples; i++) {
                for (j = 0; j < nchans_per_subband; j++) {
                    ddata[i][j] = tkbk_data[i * chn_per_sample + j + nskip];
                }
            }

            if (unflip == 1) {
                kstart = chn_per_sample - nskip - 1;    /* mam changed 10-22-96 */
                kinc = -1;
            } else {
                kstart = nskip;
                kinc = 1;
            }

            for (i = tkbk_samples; i < samples_all; i++) {
                k = kstart + (i - tkbk_samples) * chn_per_sample;       /* mam changed 10-22-96 */
                for (j = 0; j < nchans_per_subband; j++) {
                    ddata[i][j] = invert * inp_buf[k];
                    k += kinc;
                }
            }


            if (quiet == 1)
                printf("done with filling subband array");



/* SAVE LAST tkbk_samples INTO tkbk_data  for this subband */

            for (i = samples_all - tkbk_samples; i < samples_all; i++) {
                ii = i - (samples_all - tkbk_samples);
                for (j = 0; j < nchans_per_subband; j++)
                    tkbk_data[ii * chn_per_sample + j + nskip] = ddata[i][j];
            }

/* DEDISPERSE THIS SUBBAND 
   DEFAULT OPTION IS TO USE THE TREE ALGORITHM
   -m DM OPTION SETS do_single_dm == true --> DEDISPERSE WITH DM
*/

            /* if(quiet == 1) {
               fprintf(ddisp_diag, "calling dedisperse %d %d %d %d %d\n",
               loop, loop_max,  shift, samples_all, 
               nchans_per_subband); 
               fflush(ddisp_diag);
               } */

            if (do_single_dm) {
                if (quiet == 1)
                    fprintf(stderr, "dedispersing with single dm value\n");
                dedisperse_single_dm(ddata, samples_all, nchans_per_subband,
                                     chn_per_sample, &delays, ns, output);
            } else {
                if (dedisp_type > 1) {
                    if (quiet == 1)
                        fprintf(stderr, "dedispersing with dmlist\n");
                    dedisperse_dmlist(ddata, samples_all, nchans_per_subband,
                                      chn_per_sample, &delays, ns, ndm, output);
                } else {
                    if (quiet == 1)
                        fprintf(stderr, "dedispersing with tree\n");
                    dedisperse_tree(loop, shift, samples_all, nchans_per_subband,
                                    nchans_per_subband, output, ddata);
                }
            }

            /* if(quiet == 1) {
               fprintf(ddisp_diag, "finished dedispersion\n");
               fflush(ddisp_diag);
               } */

/* ooo output section */
/* ***rewrite to open files once at the beginning instead of open/close
      every time */
/* WRITE OUT DATA: n_dm_channels_per_subband FILES THIS TIME */

            /* if(quiet == 1) {
               fprintf(ddisp_diag, "begin writing out data\n");
               fflush(ddisp_diag);
               } */

            nsave = 0;

            for (j = 0; j < n_dm_channels_per_subband; j++) {
                for (i = 0; i < samples_out; i++) {
                    ii = i + tkbk_samples;
                    ddata[ii][j] =
                        (ddata[ii][j] - min[ns]) / (max[ns] - min[ns]) * 255.;
                    if (ddata[ii][j] < 0)
                        ddata[ii][j] = 0;
                    if (ddata[ii][j] > 255)
                        ddata[ii][j] = 255;
                    ii = i + tkbk_samples;
                    if (outtype == 'c') {
                        outdata.c[i] = ddata[ii][j];
                    } else if (outtype == 's') {
                        outdata.s[i] = ddata[ii][j];
                    } else if (outtype == 'i') {
                        outdata.i[i] = ddata[ii][j];
                    } else if (outtype == 'f') {
                        outdata.f[i] = ddata[ii][j];
                    }           /* if outtype */
                }               /* for i */

/*            if(j==0 || j==n_dm_channels_per_subband-1)
              fprintf(ddisp_diag, "opening file %s   ", dm_names[j][ns]);
              fflush(ddisp_diag); */

                if (j == ndm_save[nsave]) {

                    fdd = fopen(dm_names[j][ns], "ab");
                    /*       if(fdd == NULL) fprintf(ddisp_diag, "error opening file %s\n",
                       dm_names[j][ns]); */
                    fdd_fileno[j][ns] = fileno(fdd);


                    bytes_out = write(fdd_fileno[j][ns], outdata.c, bytes_out_req);

                    fclose(fdd);
                    nsave += 1;
                }


            }                   /* for j */
/* ooo end output section */

        }                       /* END OF LOOP OVER SUBBANDS */

    } /*********************************** END LOOP OVER ALL CHUNKS*/
/* close dm files */

    for (j = 0; j < n_subbands; j++) {
        nsave = 0;
        for (k = 0; k < n_dm_channels_per_subband; k++) {
            if (k == ndm_save[nsave]) {
                fclose(fout[k][j]);
                nsave += 1;
            }
/*if(fout[k][j] == NULL) fprintf(ddisp_diag, "error closing file %s\n",
				 dm_names[k][j]); */
        }
    }
    free(inp_buf);


  done:;
    if (bytes_inp <= 0)
        fprintf(stderr, "EOF!\n");
    fclose(ddisp_diag);
    fprintf(stderr, "\n");
    exit(0);
}



/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void processargs(argc, argv, chn_per_sample, bins_input,
                 sample_time, shift, frequency, bandwidth,
                 invert, quiet, n_subbands, ntimesamples, unflip, dm, do_single_dm,
                 outtype, bytesout_datatype, dedisp_type)

int argc;
char **argv;
char *outtype;
int *chn_per_sample, *quiet, *unflip;
int *bins_input, *shift, *invert, *n_subbands;
int *ntimesamples, *do_single_dm, *bytesout_datatype;
int *dedisp_type;
float *frequency, *bandwidth, *dm;
double *sample_time;

{
/*
	function to process a program's input command line.
	This is a template that can be customized for individual programs
	To use it you should:

	- pass in the parameters that may be changed.
	- edit the case statement below to correspond to what you want.
	- stdio.h must be added for this routine to work

	Don't forget the ** on the arguments coming in (since you want to 
	pass back data.
*/
    int getopt();               /* c lib function returns next opt */
    extern char *optarg;        /* if arg with option, this pts to it */
    extern int optind;          /* after call, ind into argv for next */
    extern int opterr;          /* if 0, getopt won't output err mesg */

    int c;                      /* Option letter returned by getopt */

    /* options to search for. :--> needs an argument */
    char *myoptions = "b:c:d:f:i:l:m:n:o:s:t:uw:z:q";


    char *USAGE =
        "Usage:  dedisp  -b bins_input  -c chn_per_sample  -d ntimesamples -f frequency  -i invert   -m dm  -n n_subbands  -q  -s sec/sample  -t dedisp_type -u unflip  -w bandwidth  -z shift";

    opterr = 0;                 /* turn off there message */


    /* LOOP OVER ALL THE OPTIONS IN LIST */

    while ((c = getopt(argc, argv, myoptions)) != -1) {
        switch (c) {
        case 'b':
            sscanf(optarg, "%d", bins_input);
            break;
        case 'c':
            sscanf(optarg, "%d", chn_per_sample);
            break;
        case 'd':
            sscanf(optarg, "%d", ntimesamples);
            break;
        case 'f':
            sscanf(optarg, "%f", frequency);
            break;
        case 'i':
            sscanf(optarg, "%d", invert);
            break;
/*
	  case 'l':
	           sscanf(optarg,"%d",loopmax);
		   break; 
*/
        case 'm':
            sscanf(optarg, "%f", dm);
            *do_single_dm = 1;
            break;
        case 'n':
            sscanf(optarg, "%d", n_subbands);
            break;
        case 'o':
            sscanf(optarg, "%c", outtype);
            if (*outtype == 'c')
                *bytesout_datatype = 1;
            if (*outtype == 's')
                *bytesout_datatype = 2;
            if (*outtype == 'i')
                *bytesout_datatype = 4;
            if (*outtype == 'f')
                *bytesout_datatype = 4;
            break;
        case 'q':
            *quiet = 1;
            break;
        case 's':
            sscanf(optarg, "%lf", sample_time);
            break;
        case 't':
            sscanf(optarg, "%d", dedisp_type);
            break;
        case 'u':
            *unflip = 1;
            break;
        case 'w':
            sscanf(optarg, "%f", bandwidth);
            break;
        case 'z':
            sscanf(optarg, "%d", shift);
            break;
        case '?':              /* if c not in myoptions, getopt rets ? */
            goto errout;
            break;
        }
    }
    fprintf(stderr, " %i %i %f %i  %f %f %i %i %i %i %i %i %i \n",
            *chn_per_sample, *bins_input, *sample_time, *shift,
            *frequency, *bandwidth, *invert, *quiet, *n_subbands, *ntimesamples,
            *unflip, *dm, *do_single_dm);
    return;

/*	HERE IF ILLEGAL OPTION OR ARGUMENT     */
  errout:fprintf(stderr, "%s\n", USAGE);
    exit(1);
}

#include        "dedisperse_tree.c"
#include        "dedisperse_dmlist.c"
#include        "dedisperse_single_dm96.c"
