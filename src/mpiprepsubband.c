#include <limits.h>
#include "presto.h"
#include "mpiprepsubband_cmd.h"
#include "mask.h"
#include "multibeam.h"
#include "bpp.h"
#include "wapp.h"
#include "mpi.h"

/* This causes the barycentric motion to be calculated once per TDT sec */
#define TDT 10.0

/* Simple linear interpolation macro */
#define LININTERP(X, xlo, xhi, ylo, yhi) ((ylo)+((X)-(xlo))*((yhi)-(ylo))/((xhi)-(xlo)))

/* Round a double or float to the nearest integer. */
/* x.5s get rounded away from zero.                */
#define NEAREST_INT(x) (int) (x < 0 ? ceil(x - 0.5) : floor(x + 0.5))

extern void write_data(FILE *outfiles[], int numfiles, float **outdata, 
		       int startpoint, int numtowrite);
extern void write_padding(FILE *outfiles[], int numfiles, float value, 
			  int numtowrite);
extern void update_stats(int N, double x, double *min, double *max,
			 double *avg, double *var);
extern void update_infodata(infodata *idata, int datawrote, int padwrote, 
			    int *barybins, int numbarybins, int downsamp);
extern void print_percent_complete(int current, int number);
extern void make_infodata_struct(void);
extern void make_maskbase_struct(void);
extern void broadcast_mask(mask *obsmask, int myid);

extern void set_PKMB_static(int ptsperblk, int bytesperpt, 
			    int numchan, double dt);
extern void get_BCPM_static(int *bytesperpt, int *bytesperblk, int *numifs, int *chan_map);
extern void set_BCPM_static(int ptsperblk, int bytesperpt, int bytesperblk, 
			    int numchan, int numifs, double dt, int *chan_map);
extern void get_WAPP_static(int *bytesperpt, int *bytesperblk, float *clip_sigma);
extern void set_WAPP_static(int ptsperblk, int bytesperpt, int bytesperblk, 
			    int numchan, float clip_sigma, double dt);

static int get_data(FILE *infiles[], int numfiles, float **outdata, 
		    mask *obsmask, double *dispdts, int **offsets, 
		    int *padding);

MPI_Datatype infodata_type;
MPI_Datatype maskbase_type;

static Cmdline *cmd;
static BPP_ifs bppifs=SUMIFS;
static int blocklen=0, blocksperread=0, bytesperblk=0, worklen=0, numchan=1;
static int local_numdms=1, myid=0, numprocs=1;
static PKMB_tapehdr hdr;

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  /* Any variable that begins with 't' means topocentric */
  /* Any variable that begins with 'b' means barycentric */
  FILE **infiles=NULL, **outfiles=NULL;
  float **outdata;
  double dtmp, *dms, avgdm=0.0, dsdt=0;
  double *dispdt, tlotoa=0.0, blotoa=0.0;
  double max=-9.9E30, min=9.9E30, var=0.0, avg=0.0;
  double *btoa=NULL, *ttoa=NULL, avgvoverc=0.0;
  char obs[3], ephem[10], rastring[50], decstring[50];
  int numinfiles, totnumtowrite, **offsets;
  int ii, jj, numadded=0, numremoved=0, padding=0;
  int numbarypts=0, numread=0, numtowrite=0, totwrote=0, datawrote=0;
  int padwrote=0, padtowrite=0, statnum=0;
  int numdiffbins=0, *diffbins=NULL, *diffbinptr=NULL;
  double local_lodm;
  char *datafilenm;
  infodata idata;
  mask obsmask;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /* Call usage() if we have no command line arguments */

  if (argc==1) {
    if (myid==0){
      Program = argv[0];
      usage();
    }
    MPI_Finalize();
    exit(1);
  }

  make_infodata_struct();
  make_maskbase_struct();

  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);
  numinfiles = cmd->argc;

#ifdef DEBUG
  showOptionValues();
#endif

  if (myid==0){
    printf("\n\n");
    printf("      Parallel Pulsar Subband De-dispersion Routine\n");
    printf("                 by Scott M. Ransom\n");
    printf("           Last Modification:  12 Feb, 2002\n\n");
    if (cmd->pkmbP){
      if (numinfiles > 1)
	printf("Reading Parkes PKMB data from %d files:\n", numinfiles);
      else
	printf("Reading Parkes PKMB data from 1 file:\n");
    } else if (cmd->bcpmP){
      if (numinfiles > 1)
	printf("Reading Green Bank BCPM data from %d files:\n", numinfiles);
      else
	printf("Reading Green Bank BCPM data from 1 file:\n");
    } else if (cmd->wappP){
      if (numinfiles > 1)
	printf("Reading Arecibo WAPP data from %d files:\n", numinfiles);
      else
	printf("Reading Arecibo WAPP data from 1 file:\n");
    }

    /* Open the raw data files */

    infiles = (FILE **)malloc(numinfiles * sizeof(FILE *));
    for (ii=0; ii<numinfiles; ii++){
      printf("  '%s'\n", cmd->argv[ii]);
      infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }
  }
  if (!cmd->numoutP)
    cmd->numout = INT_MAX;

  /* Determine the output file names and open them */

  datafilenm = (char *)calloc(strlen(cmd->outfile)+20, 1);
  local_numdms = cmd->numdms / (numprocs-1);
  if (cmd->numdms % (numprocs-1)){
    if (myid==0)
      printf("\nThe number of DMs must be divisible by (the number of processors - 1).\n\n");
    MPI_Finalize();
    exit(1);
  }
  local_lodm = cmd->lodm + (myid-1) * local_numdms * cmd->dmstep;
  dms = gen_dvect(local_numdms);
  if (myid>0){
    outfiles = (FILE **)malloc(local_numdms * sizeof(FILE *));
    for (ii=0; ii<local_numdms; ii++){
      dms[ii] = local_lodm + ii * cmd->dmstep;
      avgdm += dms[ii];
      sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
      outfiles[ii] = chkfopen(datafilenm, "wb");
    }
  }
  avgdm /= local_numdms;

  /* Read an input mask if wanted */

  if (cmd->maskfileP){
    if (myid==0)
      read_mask(cmd->maskfile, &obsmask);
    broadcast_mask(&obsmask, myid);
  } else {
    obsmask.numchan = obsmask.numint = 0;
  }  

  {
    float clip_sigma=0.0;
    double dt, T;
    int ptsperblk, bytesperpt, numifs=0;
    int chan_mapping[2*MAXNUMCHAN];
    long long N;

    if (myid==0){  /* Master */

      /* Set-up values if we are using the Parkes multibeam */

      if (cmd->pkmbP) {
	printf("\nPKMB input file information:\n");
	get_PKMB_file_info(infiles, numinfiles, &N, &ptsperblk, &numchan, 
			   &dt, &T, 1);
	bytesperpt = numchan / 8;
	chkfread(&hdr, 1, HDRLEN, infiles[0]);
	rewind(infiles[0]);
	PKMB_hdr_to_inf(&hdr, &idata);
	PKMB_update_infodata(numinfiles, &idata);
	strcpy(obs, "PK");  /* OBS code for TEMPO */
      }

      /* Set-up values if we are using the Berkeley-Caltech */
      /* Pulsar Machine (or BPP) format.                    */

      if (cmd->bcpmP) {
	printf("\nBCPM input file information:\n");
	get_BPP_file_info(infiles, numinfiles, &N, &ptsperblk, &numchan, 
			  &dt, &T, &idata, 1);
	get_BCPM_static(&bytesperpt, &bytesperblk, &numifs, chan_mapping);
	BPP_update_infodata(numinfiles, &idata);
	strcpy(obs, "GB");  /* OBS code for TEMPO */
      }

      /* Set-up values if we are using the Arecobo WAPP */
      
      if (cmd->wappP){
	printf("\nWAPP input file information:\n");
	get_WAPP_file_info(infiles, numinfiles, cmd->clip,
			   &N, &ptsperblk, &numchan, 
			   &dt, &T, &idata, 1);
	get_WAPP_static(&bytesperpt, &bytesperblk, &clip_sigma);
	WAPP_update_infodata(numinfiles, &idata);
	strcpy(obs, "AO");  /* OBS code for TEMPO */
      }
      
      /* The number of topo to bary time points to generate with TEMPO */
      numbarypts = (int) (T * 1.1 / TDT + 5.5) + 1;
    }

    MPI_Bcast(&ptsperblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bytesperpt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bytesperblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numchan, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numifs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numbarypts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(chan_mapping, 2*MAXNUMCHAN, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&clip_sigma, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    blocklen = ptsperblk;
    
    if (myid>0){ /* Slave */
      if (cmd->pkmbP)
	set_PKMB_static(ptsperblk, bytesperpt, numchan, dt);
      if (cmd->bcpmP)
	set_BCPM_static(ptsperblk, bytesperpt, bytesperblk, 
			numchan, numifs, dt, chan_mapping);
      if (cmd->wappP)
	set_WAPP_static(ptsperblk, bytesperpt, bytesperblk, 
			numchan, clip_sigma, dt);
    }
    if (cmd->bcpmP) {
      /* Which IFs will we use? */
      if (cmd->ifsP){
	if (cmd->ifs==0)
	  bppifs = IF0;
	else if (cmd->ifs==1)
	  bppifs = IF1;
	else
	  bppifs = SUMIFS;
      }
    }
  }

  /* Broadcast or calculate a few extra important values */

  MPI_Bcast(&idata, 1, infodata_type, 0, MPI_COMM_WORLD);
  dsdt = cmd->downsamp * idata.dt;
  idata.dm = avgdm;
  blocksperread = ((int)(delay_from_dm(cmd->lodm+cmd->numdms*cmd->dmstep, 
				       idata.freq)/dsdt) / blocklen + 1);
  worklen = blocklen * blocksperread;

  if (blocklen % cmd->downsamp){
    if (myid==0){
      printf("Error:  The downsample factor (%d) must be a factor of the\n",
	     cmd->downsamp);
      printf("        blocklength (%d).  Exiting.\n\n", 
	     blocklen);
    }
    MPI_Finalize();
    exit(1);
  }

  tlotoa = idata.mjd_i + idata.mjd_f;  /* Topocentric epoch */

  if (cmd->numoutP)
    totnumtowrite = cmd->numout;
  else
    totnumtowrite = (int) idata.N/cmd->downsamp;

  if (cmd->nobaryP) { /* Main loop if we are not barycentering... */
    
    /* Dispersion delays (in bins).  The high freq gets no delay   */
    /* All other delays are positive fractions of bin length (dt)  */
    
    dispdt = subband_search_delays(numchan, cmd->numsub, avgdm, 
				   idata.freq, idata.chan_wid, 0.0);
    for (ii=0; ii<numchan; ii++)
      dispdt[ii] /= idata.dt;
    
    /* The subband dispersion delays (see note above) */

    offsets = gen_imatrix(local_numdms, cmd->numsub);
    for (ii=0; ii<local_numdms; ii++){
      double *subdispdt;

      subdispdt = subband_delays(numchan, cmd->numsub, dms[ii], 
				 idata.freq, idata.chan_wid, 0.0);
      dtmp = subdispdt[cmd->numsub-1];
      for (jj=0; jj<cmd->numsub; jj++)
	offsets[ii][jj] = (int)((subdispdt[jj] - dtmp) / dsdt + 0.5);
      free(subdispdt);
    }

    /* Allocate our data array and start getting data */
    
    if (myid==0){
      printf("De-dispersing using %d subbands.\n", cmd->numsub);
      if (cmd->downsamp > 1)
	printf("Downsampling by a factor of %d (new dt = %.10g)\n", cmd->downsamp, dsdt);
      printf("\n");
    }

    outdata = gen_fmatrix(cmd->numsub, worklen/cmd->downsamp);
    numread = get_data(infiles, numinfiles, outdata, 
		       &obsmask, dispdt, offsets, &padding);

    while (numread==worklen){

      numread /= cmd->downsamp;
      if (myid==0)
	print_percent_complete(totwrote, totnumtowrite);
      
      /* Write the latest chunk of data, but don't   */
      /* write more than cmd->numout points.         */
      
      numtowrite = numread;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      if (myid>0){
	write_data(outfiles, local_numdms, outdata, 0, numtowrite);
	/* Update the statistics */
	if (!padding){
	  for (ii=0; ii<numtowrite; ii++)
	    update_stats(statnum+ii, outdata[0][ii], &min, &max, &avg, &var);
	  statnum += numtowrite;
	}
      }
      totwrote += numtowrite;
      
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote == cmd->numout))
	break;
      
      numread = get_data(infiles, numinfiles, outdata, 
			 &obsmask, dispdt, offsets, &padding);
    }
    datawrote = totwrote;

  } else { /* Main loop if we are barycentering... */

    /* What ephemeris will we use?  (Default is DE200) */

    if (cmd->de405P)
      strcpy(ephem, "DE405");
    else
      strcpy(ephem, "DE200");

    /* Define the RA and DEC of the observation */
    
    ra_dec_to_string(rastring, idata.ra_h, idata.ra_m, idata.ra_s);
    ra_dec_to_string(decstring, idata.dec_d, idata.dec_m, idata.dec_s);

    /* Allocate some arrays */

    btoa = gen_dvect(numbarypts);
    ttoa = gen_dvect(numbarypts);
    for (ii=0; ii<numbarypts; ii++)
      ttoa[ii] = tlotoa + TDT * ii / SECPERDAY;

    /* Call TEMPO for the barycentering */

    if (myid==0){
      double maxvoverc=-1.0, minvoverc=1.0, *voverc=NULL;

      printf("Generating barycentric corrections...\n");
      voverc = gen_dvect(numbarypts);
      barycenter(ttoa, btoa, voverc, numbarypts, \
		 rastring, decstring, obs, ephem);
      for (ii=0; ii<numbarypts; ii++){
	if (voverc[ii] > maxvoverc) maxvoverc = voverc[ii];
	if (voverc[ii] < minvoverc) minvoverc = voverc[ii];
	avgvoverc += voverc[ii];
      }
      avgvoverc /= numbarypts;
      free(voverc);
      
      printf("   Insure you check the files tempoout_times.tmp and\n");
      printf("   tempoout_vels.tmp for errors from TEMPO when complete.\n");
      printf("   Average topocentric velocity (c) = %.7g\n", avgvoverc);
      printf("   Maximum topocentric velocity (c) = %.7g\n", maxvoverc);
      printf("   Minimum topocentric velocity (c) = %.7g\n\n", minvoverc);
      printf("De-dispersing using %d subbands.\n", cmd->numsub);
      if (cmd->downsamp > 1)
	printf("Downsampling by a factor of %d (new dt = %.10g)\n", 
	       cmd->downsamp, dsdt);
      printf("\n");
    }
    MPI_Bcast(btoa, numbarypts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&avgvoverc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    blotoa = btoa[0];

    /* Dispersion delays (in bins).  The high freq gets no delay   */
    /* All other delays are positive fractions of bin length (dt)  */
    
    dispdt = subband_search_delays(numchan, cmd->numsub, avgdm, 
				   idata.freq, idata.chan_wid, avgvoverc);
    for (ii=0; ii<numchan; ii++)
      dispdt[ii] /= idata.dt;
    
    /* The subband dispersion delays (see note above) */

    offsets = gen_imatrix(local_numdms, cmd->numsub);
    for (ii=0; ii<local_numdms; ii++){
      double *subdispdt;

      subdispdt = subband_delays(numchan, cmd->numsub, dms[ii], 
				 idata.freq, idata.chan_wid, avgvoverc);
      dtmp = subdispdt[cmd->numsub-1];
      for (jj=0; jj<cmd->numsub; jj++)
	offsets[ii][jj] = (int)((subdispdt[jj] - dtmp) / dsdt + 0.5);
      free(subdispdt);
    }

    /* Convert the bary TOAs to differences from the topo TOAs in */
    /* units of bin length (dt) rounded to the nearest integer.   */

    dtmp = (btoa[0] - ttoa[0]);
    for (ii=0; ii<numbarypts; ii++)
      btoa[ii] = ((btoa[ii] - ttoa[ii]) - dtmp) * SECPERDAY / dsdt;

    { /* Find the points where we need to add or remove bins */

      int oldbin=0, currentbin;
      double lobin, hibin, calcpt;

      numdiffbins = abs(NEAREST_INT(btoa[numbarypts-1])) + 1;
      diffbins = gen_ivect(numdiffbins);
      diffbinptr = diffbins;
      for (ii=1; ii<numbarypts; ii++){
	currentbin = NEAREST_INT(btoa[ii]);
	if (currentbin != oldbin){
	  if (currentbin > 0){
	    calcpt = oldbin + 0.5;
	    lobin = (ii-1) * TDT / dsdt;
	    hibin = ii * TDT / dsdt;
	  } else {
	    calcpt = oldbin - 0.5;
	    lobin = -((ii-1) * TDT / dsdt);
	    hibin = -(ii * TDT / dsdt);
	  }
	  while(fabs(calcpt) < fabs(btoa[ii])){
	    /* Negative bin number means remove that bin */
	    /* Positive bin number means add a bin there */
	    *diffbinptr = NEAREST_INT(LININTERP(calcpt, btoa[ii-1],
						btoa[ii], lobin, hibin));
	    diffbinptr++;
	    calcpt = (currentbin > 0) ? calcpt + 1.0 : calcpt - 1.0;
	  }
	  oldbin = currentbin;
	}
      }
      *diffbinptr = INT_MAX;  /* Used as a marker */
    }
    diffbinptr = diffbins;

    /* Now perform the barycentering */

    outdata = gen_fmatrix(cmd->numsub, worklen/cmd->downsamp);
    numread = get_data(infiles, numinfiles, outdata, 
		       &obsmask, dispdt, offsets, &padding);
    
    while (numread==worklen){ /* Loop to read and write the data */
      int numwritten=0;

      numread /= cmd->downsamp;
      if (myid==0)
	print_percent_complete(totwrote, totnumtowrite);

      /* Simply write the data if we don't have to add or */
      /* remove any bins from this batch.                 */
      /* OR write the amount of data up to cmd->numout or */
      /* the next bin that will be added or removed.      */
      
      numtowrite = abs(*diffbinptr) - datawrote;
      if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	numtowrite = cmd->numout - totwrote;
      if (numtowrite > numread)
	numtowrite = numread;
      if (myid>0){
	write_data(outfiles, local_numdms, outdata, 0, numtowrite);
	/* Update the statistics */
	if (!padding){
	  for (ii = 0; ii < numtowrite; ii++)
	    update_stats(statnum + ii, outdata[0][ii], &min, &max, 
			 &avg, &var);
	  statnum += numtowrite;
	}
      }
      datawrote += numtowrite;
      totwrote += numtowrite;
      numwritten += numtowrite;
      
      if ((datawrote == abs(*diffbinptr)) &&
	  (numwritten != numread) &&
	  (totwrote < cmd->numout)){ /* Add/remove a bin */
	int skip, nextdiffbin;
	
	skip = numtowrite;
	
	do { /* Write the rest of the data after adding/removing a bin  */
	  
	  if (*diffbinptr > 0){
	    /* Add a bin */
	    if (myid>0)
	      write_padding(outfiles, local_numdms, avg, 1);
	    numadded++;
	    totwrote++;
	  } else {
	    /* Remove a bin */
	    numremoved++;
	    datawrote++;
	    numwritten++;
	    skip++;
	  }
	  diffbinptr++;
	  
	  /* Write the part after the diffbin */
	  
	  numtowrite = numread - numwritten;
	  if (cmd->numoutP && (totwrote + numtowrite) > cmd->numout)
	    numtowrite = cmd->numout - totwrote;
	  nextdiffbin = abs(*diffbinptr) - datawrote;
	  if (numtowrite > nextdiffbin)
	    numtowrite = nextdiffbin;
	  if (myid>0){
	    write_data(outfiles, local_numdms, outdata, skip, numtowrite);
	    /* Update the statistics and counters */
	    if (!padding){
	      for (ii=0; ii<numtowrite; ii++)
		update_stats(statnum+ii, outdata[0][skip+ii], 
			     &min, &max, &avg, &var);
	      statnum += numtowrite;
	    }
	  }
	  numwritten += numtowrite;
	  datawrote += numtowrite;
	  totwrote += numtowrite;
	  skip += numtowrite;

	   /* Stop if we have written out all the data we need to */
      
	  if (cmd->numoutP && (totwrote == cmd->numout)) 
	    break;
	} while (numwritten < numread);
      }
      /* Stop if we have written out all the data we need to */
      
      if (cmd->numoutP && (totwrote == cmd->numout))
	break;

      numread = get_data(infiles, numinfiles, outdata, 
			 &obsmask, dispdt, offsets, &padding);
    }
  }

  if (myid>0){

    /* Calculate the amount of padding we need  */

    if (cmd->numoutP && (cmd->numout > totwrote))
      padwrote = padtowrite = cmd->numout - totwrote;
    
    /* Write the new info file for the output data */
    
    idata.dt = dsdt;
    update_infodata(&idata, totwrote, padtowrite, diffbins, 
		    numdiffbins, cmd->downsamp);
    for (ii=0; ii<local_numdms; ii++){
      idata.dm = dms[ii];
      if (!cmd->nobaryP) {
	double baryepoch, barydispdt, baryhifreq;
	
	baryhifreq = idata.freq+(numchan-1)*idata.chan_wid;
	barydispdt = delay_from_dm(dms[ii], doppler(baryhifreq, avgvoverc));
	baryepoch = blotoa - (barydispdt / SECPERDAY);
	idata.bary = 1;
	idata.mjd_i = (int) floor(baryepoch);
	idata.mjd_f = baryepoch - idata.mjd_i;
      }
      sprintf(idata.name, "%s_DM%.2f", cmd->outfile, dms[ii]);
      writeinf(&idata);
    }
    
    /* Set the padded points equal to the average data point */
    
    if (idata.numonoff > 1){
      int index, startpad, endpad;
      
      for (ii=0; ii<local_numdms; ii++){
	fclose(outfiles[ii]);
	sprintf(datafilenm, "%s_DM%.2f.dat", cmd->outfile, dms[ii]);
	outfiles[ii] = chkfopen(datafilenm, "rb+");
      }
      for (ii=0; ii<idata.numonoff; ii++){
	index = 2 * ii;
	startpad = idata.onoff[index+1];
	if (ii==idata.numonoff-1)
	  endpad = idata.N - 1;
	else
	  endpad = idata.onoff[index+2];
	for (jj=0; jj<local_numdms; jj++)
	  chkfseek(outfiles[jj], (startpad+1)*sizeof(float), SEEK_SET);
	padtowrite = endpad - startpad;
	write_padding(outfiles, local_numdms, avg, padtowrite);
      }
    }
  }

  /* Print simple stats and results */

  var /= (datawrote - 1);
  if (myid==0)
    print_percent_complete(1, 1);
  if (myid==1){
    printf("\n\nDone.\n\nSimple statistics of the output data:\n");
    printf("             Data points written:  %d\n", totwrote);
    if (padwrote)
      printf("          Padding points written:  %d\n", padwrote);
    if (!cmd->nobaryP) {
      if (numadded)
	printf("    Bins added for barycentering:  %d\n", numadded);
      if (numremoved)
	printf("  Bins removed for barycentering:  %d\n", numremoved);
    }
    printf("           Maximum value of data:  %.2f\n", max);
    printf("           Minimum value of data:  %.2f\n", min);
    printf("              Data average value:  %.2f\n", avg);
    printf("         Data standard deviation:  %.2f\n", sqrt(var));
    printf("\n");
  }

  /* Close the files and cleanup */

  if (cmd->maskfileP)
    free_mask(obsmask);
  if (myid==0){
    for (ii=0; ii<numinfiles; ii++)
      fclose(infiles[ii]);
    free(infiles);
  } else {
    for (ii=0; ii<local_numdms; ii++)
      fclose(outfiles[ii]);
    free(outfiles);
  }
  free(outdata[0]);
  free(outdata);
  free(dms);
  free(dispdt);
  free(offsets[0]);
  free(offsets);
  free(datafilenm);
  if (!cmd->nobaryP){
    free(btoa);
    free(ttoa);
    free(diffbins);
  }
  MPI_Finalize();
  return (0);
}


static int get_data(FILE *infiles[], int numfiles, float **outdata, 
		    mask *obsmask, double *dispdts, int **offsets, 
		    int *padding)
{
  static int firsttime=1, worklen, *maskchans=NULL, blocksize;
  static int dsworklen;
  static float *tempzz, *data1, *data2, *dsdata1=NULL, *dsdata2=NULL; 
  static float *currentdata, *lastdata, *currentdsdata, *lastdsdata;
  static unsigned char *rawdata=NULL;
  int totnumread=0, numread=0, tmpnumread=0, ii, jj, tmppad=0, nummasked=0;
  
  if (firsttime){
    if (cmd->maskfileP)
      maskchans = gen_ivect(numchan);
    worklen = blocklen * blocksperread;
    dsworklen = worklen / cmd->downsamp;
    blocksize = blocklen * cmd->numsub;
    data1 = gen_fvect(cmd->numsub * worklen);
    data2 = gen_fvect(cmd->numsub * worklen);
    rawdata = gen_bvect(bytesperblk);
    currentdata = data1;
    lastdata = data2;
    if (cmd->downsamp > 1){
      dsdata1 = gen_fvect(cmd->numsub * dsworklen);
      dsdata2 = gen_fvect(cmd->numsub * dsworklen);
      currentdsdata = dsdata1;
      lastdsdata = dsdata2;
    } else {
      currentdsdata = data1;
      lastdsdata = data2;
    }
  }
  while (firsttime >= 0){
    if (cmd->pkmbP || cmd->bcpmP || cmd->wappP){
      for (ii=0; ii<blocksperread; ii++){
	if (myid==0){
	  if (cmd->pkmbP)
	    numread = read_PKMB_rawblock(infiles, numfiles, &hdr, rawdata, &tmppad);
	  if (cmd->bcpmP)
	    numread = read_BPP_rawblock(infiles, numfiles, rawdata, &tmppad);
	  if (cmd->wappP)
	    numread = read_WAPP_rawblock(infiles, numfiles, rawdata, &tmppad);
	  numread *= blocklen;
	}
	MPI_Bcast(&tmppad, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numread, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(rawdata, bytesperblk, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	if (myid>0){
	  if (numread){
	    if (cmd->pkmbP)
	      tmpnumread = prep_PKMB_subbands(rawdata, currentdata+ii*blocksize, 
					      dispdts, cmd->numsub, 0, 
					      maskchans, &nummasked, obsmask);
	    if (cmd->bcpmP)
	      tmpnumread = prep_BPP_subbands(rawdata, currentdata+ii*blocksize, 
					     dispdts, cmd->numsub, 0, 
					     maskchans, &nummasked, obsmask, bppifs);
	    if (cmd->wappP)
	      tmpnumread = prep_WAPP_subbands(rawdata, currentdata+ii*blocksize, 
					      dispdts, cmd->numsub, 0, 
					      maskchans, &nummasked, obsmask);
	  } else {
	    *padding = 1;
	    for (jj=ii*blocksize; jj<(ii+1)*blocksize; jj++)
	      currentdata[jj] = 0.0;
	  }
	  if (tmppad) 
	    *padding = 1;
	}
	if (firsttime==0) totnumread += numread;
      }
    }
    /* Downsample the subband data if needed */
    if (myid>0){
      if (cmd->downsamp > 1){
	int kk, offset, dsoffset, index, dsindex;
	for (ii=0; ii<dsworklen; ii++){
	  dsoffset = ii * cmd->numsub;
	  offset = dsoffset * cmd->downsamp;
	  for (jj=0; jj<cmd->numsub; jj++){
	    dsindex = dsoffset + jj;
	    index = offset + jj;
	    currentdsdata[dsindex] = 0.0;
	    for (kk=0; kk<cmd->downsamp; kk++){
	      currentdsdata[dsindex] += currentdata[index];
	      index += cmd->numsub;
	    }
	  }
	}
      }
      for (ii=0; ii<local_numdms; ii++)
	float_dedisp(currentdsdata, lastdsdata, dsworklen, 
		     cmd->numsub, offsets[ii], 0.0, outdata[ii]);
      SWAP(currentdata, lastdata);
      SWAP(currentdsdata, lastdsdata);
    }
    firsttime--;
  }
  firsttime = 0;
/*
{
  int jj;
  for (jj=0; jj<numprocs; jj++){
    if (myid==jj)
      printf("%d:  %d  %d  %d\n", myid, numread, totnumread, worklen);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
*/
  if (totnumread != worklen){
    if (cmd->maskfileP)
      free(maskchans);
    free(data1);
    free(data2);
    free(rawdata);
    if (cmd->downsamp > 1){
      free(dsdata1);
      free(dsdata2);
    }
  }
  return totnumread;
}


