#include "presto.h"
#include "cpgplot.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define DEBUGOUT 0
 
/* Note:  zoomlevel is simply (LOGDISPLAYNUM-Log_2(numsamps)) */
#define LOGDISPLAYNUM      11 /* 2048: Maximum number of points to display at once */
#define LOGMAXBINS         23 /* 8M points */
#define LOGINITIALNUMSAMPS 16 /* 65536: The initial number of samples to plot */
#define DISPLAYNUM      (1<<LOGDISPLAYNUM)
#define MAXBINS         (1<<LOGMAXBINS)
#define INITIALNUMSAMPS (1<<LOGINITIALNUMSAMPS)

static long long N;    /* Number of points in the time series */
static double T;       /* The time duration of data */
static float mjd0;     /* The MJD of the beginning of the first point */
static infodata idata;
static FILE *datfile;

typedef struct datapart {
  double tlo;  /* Elapsed time (s) of the first point */
  double avg;  /* The average of this chunk of data */
  double med;  /* The median of this chunk of data */
  double std;  /* The standard deviation of this chunk of data */
  int nlo      /* The sample number of the first point */
  int nn;      /* The total number samples in *data */
  float *data; /* Raw data  */
} datapart;

typedef struct dataview {
  double vdt;       /* Data view time step (2.0**(-zoomlevel))*dt */
  int centern;      /* The center sample to plot */
  int lon;          /* The lowest sample to plot */
  int zoomlevel;    /* Positive = zoomed in, Negative = zoomed out */
  int numsamps;     /* The number of samples from low to high to display */
  int chunklen;     /* The length of the chunk of samples used to calculate stats */
  float avgs[DISPLAYNUM];  /* The average samples for each chunk */
  float meds[DISPLAYNUM];  /* The median samples for each chunk */
  float stds[DISPLAYNUM];  /* The atandard deviation of the samples for each chunk */
  float maxs[DISPLAYNUM];  /* The maximum samples for each chunk */
  float mins[DISPLAYNUM];  /* The minimum samples for each chunk */
} dataview;

typedef enum {
  MINMAX, STATS
} plottype;


static double plot_dataview(dataview *dv, float minval, float maxval, 
			    float charhgt, plottype type)
/* The return value is offsetf */
{
  int ii;
  double lof, hif, offsetf=0.0;
  float *freqs;

  cpgsave();
  cpgbbuf();

  /* Set the "Normal" plotting attributes */

  cpgsls(1);
  cpgslw(1);
  cpgsch(charhgt);
  cpgsci(1);
  cpgvstd();

  if (maxpow==0.0) /* Autoscale for the maximum value */
    maxpow = 1.1 * dv->maxpow;

  lof = dv->lor / T;
  hif = (dv->lor + dv->dr * DISPLAYNUM) / T;
  offsetf = 0.0;

  /* Period Labels */

  if (dv->zoomlevel >= 0 && lof > 1.0){
    double lop, hip, offsetp=0.0;
    lop = 1.0/lof;
    hip = 1.0/hif;
    offsetp = 0.0;
    
    if ((lop-hip)/hip < 0.001){
      int numchar;
      char label[50];
      
      offsetp = 0.5*(hip+lop);
      numchar = snprintf(label, 50, "Period - %.15g (s)", offsetp);
      cpgmtxt("T", 2.5, 0.5, 0.5, label);
    } else {
      cpgmtxt("T", 2.5, 0.5, 0.5, "Period (s)");
    }
    cpgswin(lop-offsetp, hip-offsetp, 0.0, maxpow);
    cpgbox("CIMST", 0.0, 0, "", 0.0, 0);
  }

  /* Frequency Labels */

  if ((hif-lof)/hif < 0.001){
    int numchar;
    char label[50];

    offsetf = 0.5*(hif+lof);
    numchar = snprintf(label, 50, "Frequency - %.15g (Hz)", offsetf);
    cpgmtxt("B", 2.8, 0.5, 0.5, label);
  } else {
    cpgmtxt("B", 2.8, 0.5, 0.5, "Frequency (Hz)");
  }
  cpgswin(lof-offsetf, hif-offsetf, 0.0, maxpow);
  if (dv->zoomlevel >= 0 && lof > 1.0)
    cpgbox("BINST", 0.0, 0, "BCNST", 0.0, 0);
  else
    cpgbox("BCINST", 0.0, 0, "BCNST", 0.0, 0);

  /* Plot the spectrum */

  freqs = gen_fvect(DISPLAYNUM);
  for (ii=0; ii<DISPLAYNUM; ii++)
    freqs[ii] = dv->rs[ii]/T - offsetf;
  if (dv->zoomlevel > 0){ /* Magnified power spectrum */
    cpgline(DISPLAYNUM, freqs, dv->powers);
  } else { /* Down-sampled power spectrum */
    for (ii=0; ii<DISPLAYNUM; ii++){
      cpgmove(freqs[ii], 0.0);
      cpgdraw(freqs[ii], dv->powers[ii]);
    }
  }
  free(freqs);
  cpgmtxt("L", 2.5, 0.5, 0.5, "Normalized Power");
  cpgebuf();
  cpgunsa();
  return offsetf;
}


static dataview *get_dataview(int centern, int zoomlevel, datapart *fp)
{
  int ii, jj, powindex, normindex;
  float tmpchunk;
  double split=1.0/(double)LOCALCHUNK;
  dataview *dv;

  fv = (dataview *)malloc(sizeof(dataview));
  dv->zoomlevel = zoomlevel;
  dv->chunklen = (1 << abs(zoomlevel));
  dv->vdt = dv->chunklen*dt;
  dv->centern = centern;
  dv->numsamps = DISPLAYNUM * dv->chunklen;
  dv->dr = (double) dv->chunklen;
  dv->lon = (int) floor(centern - 0.5 * dv->numsamps);
  if (dv->lon < 0) dv->lon = 0;
  tmprawpwrs = gen_fvect(dv->numsamps);
  if (norm_const==0.0){
    for (ii=0; ii<dv->numsamps; ii++){
      powindex = (int) (dv->lor - dp->rlo + ii + 0.5);
      normindex = (int) (powindex * split);
      tmprawpwrs[ii] = dp->rawpowers[powindex] * dp->normvals[normindex];
    }
  } else {
    for (ii=0; ii<dv->numsamps; ii++){
      powindex = (int) (dv->lor - dp->rlo + ii + 0.5);
      tmprawpwrs[ii] = dp->rawpowers[powindex] * norm_const;
    }
  }
  powindex = 0;
  for (ii=0; ii<DISPLAYNUM; ii++){
    maxpow = 0.0;
    for (jj=0; jj<dv->chunklen; jj++, powindex++)
      if (tmprawpwrs[powindex] > maxpow) maxpow = tmprawpwrs[powindex];
    dv->rs[ii] = dv->lor + ii * dv->dr;
    dv->powers[ii] = maxpow;
  }
  dv->maxpow = 0.0;
  for (ii=0; ii<DISPLAYNUM; ii++)
    if (dv->powers[ii] > dv->maxpow) dv->maxpow = dv->powers[ii];
  return fv;
}


static datapart *get_datapart(int nlo, int numn)
{
  datapart *dp;

  if (nlo+numn > N)
    return NULL;
  else {
    float *tmpdata;

    dp = (datapart *)malloc(sizeof(datapart));
    dp->nlo = nlo;
    dp->numamps = numn;
    dp->tlo = idata.dt * nlo;
    dp->data = read_float_file(datfile, nlo, numn);
    tmpdata = gen_fvect(numn);
    memcpy(tmpdata, dp->data, sizeof(float)*numn);
    dp->median = median(tmpdata, numn);
    free(tmpdata);
    avg_var(dp->data, numn, &(dp->avg), &(dp->std));
    dp->std = sqrt(dp->std);
    return dp;
  }
}

static void free_datapart(datapart *dp)
{
  free(dp->data);
  free(dp);
}


static void print_help(void)
{
  printf("\n"
	 " Button or Key            Effect\n"
	 " -------------            ------\n"
	 " Left Mouse or I or A     Zoom in  by a factor of 2\n"
	 " Right Mouse or O or X    Zoom out by a factor of 2\n"
	 " Middle Mouse or D        Show details about a selected frequency\n"
	 " < or ,                   Shift left  by 15%% of the screen width\n"
	 " > or .                   Shift right by 15%% of the screen width\n"
	 " + or =                   Increase the power scale (make them taller)\n"
	 " - or _                   Decrease the power scale (make them shorter)\n"
	 " S                        Scale the powers automatically\n"
	 " N                        Re-normalize the nowers by one of several methods\n"
	 " P                        Print the current plot to a file\n"
	 " G                        Go to a specified frequency\n"
	 " H                        Show the harmonics of the center frequency\n"
	 " ?                        Show this help screen\n"
	 " Q                        Quit\n"
	 "\n");
}


int main(int argc, char *argv[])
{
  float maxpow=0.0;
  double centern, offsetf;
  int numamps, zoomlevel, maxzoom, minzoom, xid, psid;
  char *rootfilenm, inchar;
  datapart *lofp;
  dataview *fv;
 
  printf("\n\n");
  printf("      Interactive FFT Explorer\n");
  printf("         by Scott M. Ransom\n");
  printf("            October, 2001\n");
  print_help();

  {
    int hassuffix=0;
    char *suffix;
    
    hassuffix = split_root_suffix(argv[1], &rootfilenm, &suffix);
    if (hassuffix){
      if (strcmp(suffix, "fft")!=0){
        printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n", argv[1]);
        free(suffix);
        exit(0);
      }
      free(suffix);
    } else {
      printf("\nInput file ('%s') must be a FFT file ('.fft')!\n\n", argv[1]);
      exit(0);
    }
  }

  /* Read the info file */
 
  readinf(&idata, rootfilenm);
  if (idata.object) {
    printf("Examining %s data from '%s'.\n\n",
           remove_whitespace(idata.object), argv[1]);
  } else {
    printf("Examining data from '%s'.\n\n", argv[1]);
  }
  N = idata.N;
  T = idata.dt * idata.N;
  datfile = chkfopen(argv[1], "rb");
  Nfft = chkfilelen(datfile, sizeof(fcomplex));

  /* Get and plot the initial data */
  
  numamps = (Nfft > MAXBINS) ? (int) MAXBINS : (int) Nfft;
  lofp = get_datapart(0, numamps);
  r0 = lodp->amps[0].r;
  centern = 0.5 * INITIALNUMSAMPS;
  zoomlevel = LOGDISPLAYNUM - LOGINITIALNUMSAMPS;
  minzoom = LOGDISPLAYNUM - LOGMAXBINS;
  maxzoom = LOGDISPLAYNUM - LOGMINBINS;
  fv = get_dataview(centern, zoomlevel, lofp);

  /* Prep the XWIN device for PGPLOT */

  xid = cpgopen("/XWIN");
  if (xid <= 0){
    fclose(datfile);
    free(fv);
    free_datapart(lofp);
    exit(EXIT_FAILURE);
  }
  cpgask(0);
  cpgpage();
  offsetf = plot_dataview(fv, maxpow, 1.0);

  do {
    float inx, iny;
    
    cpgcurs(&inx, &iny, &inchar);
    if (DEBUGOUT) printf("You pressed '%c'\n", inchar);

    switch (inchar){
    case 'A': /* Zoom in */
    case 'a':
      centern = (inx + offsetf) * T;
    case 'I':
    case 'i':
      if (DEBUGOUT) printf("  Zooming in  (zoomlevel = %d)...\n", zoomlevel);
      if (zoomlevel < maxzoom){
	zoomlevel++;
	free(fv);
	fv = get_dataview(centern, zoomlevel, lofp);
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
      } else 
	printf("  Already at maximum zoom level (%d).\n", zoomlevel);
      break;
    case 'X': /* Zoom out */
    case 'x':
    case 'O':
    case 'o':
      if (DEBUGOUT) printf("  Zooming out  (zoomlevel = %d)...\n", zoomlevel);
      if (zoomlevel > minzoom){
	zoomlevel--;
	free(fv);
	fv = get_dataview(centern, zoomlevel, lofp);
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
      } else 
	printf("  Already at minimum zoom level (%d).\n", zoomlevel);
      break;
    case '<': /* Shift left */
    case ',':
      if (DEBUGOUT) printf("  Shifting left...\n");
      centern -= 0.15 * dv->numsamps;
      { /* Should probably get the previous chunk from the datfile... */
	double lowestr;

	lowestr = 0.5 * dv->numsamps;
	if (centern < lowestr)
	  centern = lowestr;
      }
      free(fv);
      fv = get_dataview(centern, zoomlevel, lofp);
      cpgpage();
      offsetf = plot_dataview(fv, maxpow, 1.0);
      break;
    case '>': /* Shift right */
    case '.':
      if (DEBUGOUT) printf("  Shifting right...\n");
      centern += 0.15 * dv->numsamps;
      { /* Should probably get the next chunk from the datfile... */
	double highestr;

	highestr = lodp->nlo + lodp->numamps - 0.5 * dv->numsamps;
	if (centern > highestr)
	  centern = highestr;
      }
      free(fv);
      fv = get_dataview(centern, zoomlevel, lofp);
      cpgpage();
      offsetf = plot_dataview(fv, maxpow, 1.0);
      break;
    case '+': /* Increase height of powers */
    case '=':
      if (maxpow==0.0){
	printf("  Auto-scaling is off.\n");
	maxpow = 1.1 * dv->maxpow;
      }
      maxpow = 3.0/4.0 * maxpow;
      cpgpage();
      offsetf = plot_dataview(fv, maxpow, 1.0);
      break;
    case '-': /* Decrease height of powers */
    case '_':
      if (maxpow==0.0){ 
	printf("  Auto-scaling is off.\n");
	maxpow = 1.1 * dv->maxpow;
      }
      maxpow = 4.0/3.0 * maxpow;
      cpgpage();
      offsetf = plot_dataview(fv, maxpow, 1.0);
      break;
    case 'S': /* Auto-scale */
    case 's':
      if (maxpow==0.0)
	break;
      else {
	printf("  Auto-scaling is on.\n");
	maxpow = 0.0;
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
	break;
      }
    case 'G': /* Goto a frequency */
    case 'g':
      {
	char freqstr[50];
	double freq=-1.0;

	while (freq < 0.0){
	  printf("  Enter the frequency (Hz) to go to:\n");
	  fgets(freqstr, 50, stdin);
	  freqstr[strlen(freqstr)-1]= '\0';
	  freq = atof(freqstr);
	}
	offsetf = 0.0;
	centern = freq * T;
	printf("  Moving to frequency %.15g.\n", freq);
	free(fv);
	fv = get_dataview(centern, zoomlevel, lofp);
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
      }
      break;
    case 'H': /* Show harmonics */
    case 'h':
      {
	double retval;
	
	retval = harmonic_loop(xid, centern);
	if (retval > 0.0){
	  offsetf = 0.0;
	  centern = retval;
	  zoomlevel = LOGDISPLAYNUM - LOGNUMHARMBINS;
	  free(fv);
	  fv = get_dataview(centern, zoomlevel, lofp);
	  cpgpage();
	  offsetf = plot_dataview(fv, maxpow, 1.0);
	}
      }
      break;
    case '?': /* Print help screen */
      print_help();
      break;
    case 'D': /* Show details about a selected point  */
    case 'd':
      {
	double newr;
	
	printf("  Searching for peak near freq = %.7g Hz...\n", (inx + offsetf));
	newr = find_peak(inx+offsetf, fv, lofp);
	centern = newr;
	if (zoomlevel < maxzoom)
	  zoomlevel++;
	free(fv);
	fv = get_dataview(centern, zoomlevel, lofp);
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
	/*
	  cpgsci(2);
	  cpgmove(newr/T-offsetf, 0.0);
	  cpgdraw(newr/T-offsetf, 10.0*dv->maxpow);
	*/
      }
      break;
    case 'P': /* Print the current plot */
    case 'p':
      {
	int len;
	char filename[200];

	printf("  Enter the filename to save the plot as:\n");
	fgets(filename, 196, stdin);
	len = strlen(filename)-1;
	filename[len+0] = '/';
	filename[len+1] = 'P';
	filename[len+2] = 'S';
	filename[len+3] = '\0';
	psid = cpgopen(filename);
	cpgslct(psid);
	cpgpap(10.25, 8.5/11.0);
	cpgiden();
	offsetf = plot_dataview(fv, maxpow, 1.0);
	cpgclos();
	cpgslct(xid);
	filename[len] = '\0';
	printf("  Wrote the plot to the file '%s'.\n", filename);
      }
      break;
    case 'N': /* Changing power normalization */
    case 'n':
      {
	float inx2, iny2;
	char choice;
	unsigned char badchoice=1;

	printf("  Specify the type of power normalization:\n"
	       "       m,M  :  Median values determined locally\n"
	       "       d,D  :  DC frequency amplitude\n"
	       "       r,R  :  Raw powers (i.e. no normalization)\n"
	       "       u,U  :  User specified interval (the average powers)\n");
	while (badchoice){
	  cpgcurs(&inx2, &iny2, &choice);
	  switch (choice){
	  case 'M':
	  case 'm':
	    norm_const = 0.0;
	    maxpow = 0.0;
	    badchoice = 0;
	    printf("  Using local median normalization.  Autoscaling is on.\n");
	    break;
	  case 'D':
	  case 'd':
	    norm_const = 1.0/r0;
	    maxpow = 0.0;
	    badchoice = 0;
	    printf("  Using DC frequency (%f) normalization.  Autoscaling is on.\n", r0);
	    break;
	  case 'R':
	  case 'r':
	    norm_const = 1.0;
	    maxpow = 0.0;
	    badchoice = 0;
	    printf("  Using raw powers (i.e. no normalization).  Autoscaling is on.\n");
	    break;
	  case 'U':
	  case 'u':
	    {
	      char choice2;
	      float xx, yy;
	      int lor, hir, numn;
	      double avg, var;

	      printf("  Use the left mouse button to select a left and right boundary\n"
		     "  of a region to calculate the average power.\n");
	      do {
		cpgcurs(&xx, &yy, &choice2);
	      } while (choice2 != 'A' && choice2 != 'a');
	      lor = (int)((xx + offsetf) * T);
	      cpgsci(7);
	      cpgmove(xx, 0.0);
	      cpgdraw(xx, 10.0*dv->maxpow);
	      do {
		cpgcurs(&xx, &yy, &choice2);
	      } while (choice2 != 'A' && choice2 != 'a');
	      hir = (int)((xx + offsetf) * T);
	      cpgmove(xx, 0.0);
	      cpgdraw(xx, 10.0*dv->maxpow);
	      cpgsci(1);
	      if (lor > hir){
		int tempr;
		tempr = hir;
		hir = lor;
		lor = tempr;
	      }
	      numn = hir - lor + 1;
	      avg_var(lodp->rawpowers+lor-lodp->nlo, numn, &avg, &var);
	      printf("  Selection has:  average = %.5g\n"
		     "                  std dev = %.5g\n", avg, sqrt(var));
	      norm_const = 1.0/avg;
	      maxpow = 0.0;
	      badchoice = 0;
	      printf("  Using %.5g as the normalization constant.  Autoscaling is on.\n", avg);
	      break;
	    }
	  default:
	    printf("  Unrecognized choice '%c'.\n", choice);
	    break;
	  }
	}
	free(fv);
	fv = get_dataview(centern, zoomlevel, lofp);
	cpgpage();
	offsetf = plot_dataview(fv, maxpow, 1.0);
      }
      break;
    case 'Q': /* Quit */
    case 'q':
      printf("  Quitting...\n");
      free(fv);
      cpgclos();
      break;
    default:
      printf("  Unrecognized option '%c'.\n", inchar);
      break;
    }
  } while (inchar != 'Q' && inchar != 'q');

  free_datapart(lofp);
  fclose(datfile);
  printf("Done\n\n");
  return 0;
}
