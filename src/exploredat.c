#include "presto.h"
#include "cpgplot.h"
#include "float.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define DEBUGOUT 0
 
/* Note:  zoomlevel is simply (LOGMAXDISPNUM-Log_2(numsamps)) */
#define LOGMAXDISPNUM      11 /* 2048: Maximum number of points to display at once */
#define LOGMINDISPNUM      3  /* 8: Minimum number of points to display at once */
#define LOGMINCHUNKLEN     3  /* 8: The minimum number of real points in a stats chunk */
#define LOGMAXPTS          23 /* 8M points */
#define LOGINITIALNUMPTS   16 /* 65536: The initial number of samples to plot */
#define MAXDISPNUM      (1<<LOGMAXDISPNUM)
#define MINDISPNUM      (1<<LOGMINDISPNUM)
#define MINCHUNKLEN     (1<<LOGMINCHUNKLEN)
#define MAXPTS          (1<<LOGMAXPTS)
#define INITIALNUMPTS   (1<<LOGINITIALNUMPTS)

static long long N;    /* Number of points in the time series */
static double T;       /* The time duration of data */
static float mjd0;     /* The MJD of the first point in the file */
static infodata idata;
static FILE *datfile;
static int plotstats=1;

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
  int numchunks;    /* The number of chunks that are being displayed */
  float avgs[MAXDISPNUM];  /* The average samples for each chunk */
  float meds[MAXDISPNUM];  /* The median samples for each chunk */
  float stds[MAXDISPNUM];  /* The atandard deviation of the samples for each chunk */
  float maxs[MAXDISPNUM];  /* The maximum samples for each chunk */
  float mins[MAXDISPNUM];  /* The minimum samples for each chunk */
  float vals[MAXDISPNUM];  /* The raw data values when zoomlevel > 0 */
} dataview;


static int plot_dataview(dataview *dv, float minval, float maxval, 
			 float charhgt)
/* The return value is offsetn */
{
  int ii, lon, hin, offsetn=0, tmpn;
  double lot, hit, offsett;
  float ns[MAXDISPNUM], scalemin, scalemax, dscale;
  
  cpgsave();
  cpgbbuf();

  /* Set the "Normal" plotting attributes */

  cpgsls(1);
  cpgslw(1);
  cpgsch(charhgt);
  cpgsci(1);
  cpgvstd();
  
  /* Autoscale for the maximum value */
  if (maxval==FLT_MAX)
    scalemax = dv->maxval;
  else
    scalemax = maxval;
  /* Autoscale for the minimum value */
  if (minval==FLT_MIN)
    scalemin = dv->minval;
  else
    scalemax = minval;
  dscale = 0.1 * (scalemax - scalemin);
  if (maxval==FLT_MAX)
    maxval = scalemax + dscale;
  if (minval==FLT_MIN)
    minval = scalemin - dscale;
  
  lon = dv->lon;
  lot = lon * idata.dt;
  hin = lon + dv->numsamps;
  hit = hin * idata.dt;

  /* Time Labels (top of box) */

  if ((lot-hit)/hit < 0.001){
    int numchar;
    char label[50];
    
    offsett = 0.5*(hit+lot);
    numchar = snprintf(label, 50, "Time - %.15g (s)", offsett);
    cpgmtxt("T", 2.5, 0.5, 0.5, label);
  } else {
    cpgmtxt("T", 2.5, 0.5, 0.5, "Time (s)");
  }
  cpgswin(lot-offsett,  hit-offsett, minval, maxval);
  cpgbox("CMST", 0.0, 0, "", 0.0, 0);
  
  /* Sample number labels */

  if (lon > 10000000 ||
      (double)(hin-lon)/(double)hin < 0.001){
    int numchar;
    char label[50];
    
    offsetn = (lon / 1000) * 1000;
    numchar = snprintf(label, 50, "Sample - %d", offsetn);
    cpgmtxt("B", 2.8, 0.5, 0.5, label);
  } else {
    cpgmtxt("B", 2.8, 0.5, 0.5, "Sample");
  }
  cpgswin(lon-offsetn, hin-offsetn, min, maxval);
  cpgbox("BNST", 0.0, 0, "BCNST", 0.0, 0);

  /* Plot the rawdata if required */

  tmpn = lon - offsetn;
  if (zoomlevel > 0){
    for (ii=0; ii<dv->dispnum; ii++)
      ns[ii] = tmpn + ii;
    cpgbin(dv->dispnum, ns, dv->vals, 0);
  } else {  /* Plot the min/max values */
    for (ii=0; ii<dv->numchunks; ii++, tmpn += dv->chunklen){
      cpgmove((float) tmpn, dv->mins[ii]);
      cpgdraw((float) tmpn, dv->maxs[ii]);
    }
  }

  /* Plot the other statistics if requested */

  if (plotstats){
    tmpn = lon - offsetn;
    for (ii=0; ii<dv->numchunks; ii++, tmpn += dv->chunklen)
      ns[ii] = tmpn;
    if (dv->numchunks > 512){
      
      cpgline(dv->numchunks, ns, dv->avgs);
    } else {
      cpgbin(dv->numchunks, ns, dv->avgs, 0);
    }
  }
  cpgmtxt("L", 2.5, 0.5, 0.5, "Sample Value");
  cpgebuf();
  cpgunsa();
  return offsetn;
}


static dataview *get_dataview(int centern, int zoomlevel, datapart *dp)
{
  int ii, jj, offset;
  float *tmpchunk;
  dataview *dv;

  dv = (dataview *)malloc(sizeof(dataview));
  dv->zoomlevel = zoomlevel;
  dv->numsamps = (1 << (LOGMAXDISPNUM - zoomlevel));
  dv->chunklen = (zoomlevel < -LOGMINCHUNKLEN) ?
    (1 << abs(zoomlevel)) : dv->chunklen = (1 << LOGMINCHUNKLEN);
  dv->dispnum = (dv->numsamps > MAXDISPNUM) ? MAXDISPNUM : dv->numsamps;
  dv->numchunks = dv->numsamps / dv->chunklen;
  dv->vdt = dv->chunklen * idata.dt;
  dv->centern = centern;
  dv->lon = (int) floor(centern - 0.5 * dv->numsamps);
  if (dv->lon < 0) dv->lon = 0;
  tmpchunk = gen_fvect(dv->chunklen);
  for (ii=0; ii<dv->dispnum; ii++){
    float tmpmin=FLT_MAX, tmpmax=FLT_MIN, tmpval;
    offset = dv->lon + ii * dv->chunklen;
    memcpy(tmpchunk, dp->data+offset, sizeof(float)*dv->chunklen);
    dv->meds[ii] = median(tmpchunk, dv->chunklen);
    avg_var(dp->data+offset, dv->chunklen, dv->avgs+ii, dv->stds+ii);
    dv->stds[ii] = sqrt(dv->stds[ii]);
    for (jj=0; jj<dv->chunklen; jj++, offset++){
      tmpval = dp->data[offset];
      if (tmpval > tmpmax) tmpmax = tmpval;
      if (tmpval < tmpmin) tmpmin = tmpval;
    }
    dv->maxs[ii] = tmpmax;
    dv->mins[ii] = tmpmin;
  }
  free(tmpchunk);
  offset = dv->lon;
  if (zoomlevel > 0){
    for (ii=0, offset=dv->lon; ii<dv->numsamps; ii++, offset++)
      *(dv->vals+ii) = *(dp->data+offset);
  }
  return dv;
}
  

static datapart *get_datapart(int nlo, int numn)
{
  datapart *dp;

  if (nlo+numn > N)
    return NULL;
  else {
    float *tmpdata;

    dp = (datapart *)malloc(sizeof(datapart));
    dp->nn = numn;
    dp->nlo = nlo;
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
  float maxval=0.0;
  double centern, offsetn;
  int numamps, zoomlevel, maxzoom=0, minzoom, xid, psid;
  char *rootfilenm, inchar;
  datapart *lodp;
  dataview *dv;
 
  printf("\n\n");
  printf("      Interactive Data Explorer\n");
  printf("         by Scott M. Ransom\n");
  printf("            November, 2001\n");
  print_help();

  {
    int hassuffix=0;
    char *suffix;
    
    hassuffix = split_root_suffix(argv[1], &rootfilenm, &suffix);
    if (hassuffix){
      if (strcmp(suffix, "dat")!=0){
        printf("\nInput file ('%s') must be a single PRESTO data file ('.dat')!\n\n", argv[1]);
        free(suffix);
        exit(0);
      }
      free(suffix);
    } else {
      printf("\nInput file ('%s') must be a PRESTO data file ('.dat')!\n\n", argv[1]);
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
  Ndat = chkfilelen(datfile, sizeof(fcomplex));

  /* Get and plot the initial data */
  
  numamps = (Ndat > MAXPTS) ? (int) MAXPTS : (int) Ndat;
  lodp = get_datapart(0, numamps);
  centern = 0.5 * INITIALNUMPTS;
  zoomlevel = LOGMAXDISPNUM - LOGINITIALNUMPTS;
  minzoom = LOGMAXDISPNUM - LOGMAXPTS;
  maxzoom = LOGMINDISPNUM - LOGMINPTS;
  dv = get_dataview(centern, zoomlevel, lodp);

  /* Prep the XWIN device for PGPLOT */

  xid = cpgopen("/XWIN");
  if (xid <= 0){
    fclose(datfile);
    free(dv);
    free_datapart(lodp);
    exit(EXIT_FAILURE);
  }
  cpgask(0);
  cpgpage();
  offsetn = plot_dataview(dv, maxval, 1.0);

  do {
    float inx, iny;
    
    cpgcurs(&inx, &iny, &inchar);
    if (DEBUGOUT) printf("You pressed '%c'\n", inchar);

    switch (inchar){
    case 'A': /* Zoom in */
    case 'a':
      centern = (inx + offsetn) * T;
    case 'I':
    case 'i':
      if (DEBUGOUT) printf("  Zooming in  (zoomlevel = %d)...\n", zoomlevel);
      if (zoomlevel < maxzoom){
	zoomlevel++;
	free(dv);
	dv = get_dataview(centern, zoomlevel, lodp);
	cpgpage();
	offsetn = plot_dataview(dv, maxval, 1.0);
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
	free(dv);
	dv = get_dataview(centern, zoomlevel, lodp);
	cpgpage();
	offsetn = plot_dataview(dv, maxval, 1.0);
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
      free(dv);
      dv = get_dataview(centern, zoomlevel, lodp);
      cpgpage();
      offsetn = plot_dataview(dv, maxval, 1.0);
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
      free(dv);
      dv = get_dataview(centern, zoomlevel, lodp);
      cpgpage();
      offsetn = plot_dataview(dv, maxval, 1.0);
      break;
    case '+': /* Increase height of powers */
    case '=':
      if (maxval==0.0){
	printf("  Auto-scaling is off.\n");
	maxval = 1.1 * dv->maxval;
      }
      maxval = 3.0/4.0 * maxval;
      cpgpage();
      offsetn = plot_dataview(dv, maxval, 1.0);
      break;
    case '-': /* Decrease height of powers */
    case '_':
      if (maxval==0.0){ 
	printf("  Auto-scaling is off.\n");
	maxval = 1.1 * dv->maxval;
      }
      maxval = 4.0/3.0 * maxval;
      cpgpage();
      offsetn = plot_dataview(dv, maxval, 1.0);
      break;
    case 'S': /* Auto-scale */
    case 's':
      if (maxval==0.0)
	break;
      else {
	printf("  Auto-scaling is on.\n");
	maxval = 0.0;
	cpgpage();
	offsetn = plot_dataview(dv, maxval, 1.0);
	break;
      }
    case 'G': /* Goto a time */
    case 'g':
      {
	char timestr[50];
	double time=-1.0;

	while (time < 0.0){
	  printf("  Enter the time (s) from the beginning of the file to go to:\n");
	  fgets(timestr, 50, stdin);
	  timestr[strlen(timestr)-1]= '\0';
	  time = atof(timestr);
	}
	offsetn = 0.0;
	centern = time / idata.dt;
	printf("  Moving to time %.15g (data point %d).\n", time, centern);
	free(dv);
	dv = get_dataview(centern, zoomlevel, lodp);
	cpgpage();
	offsetn = plot_dataview(dv, maxval, 1.0);
      }
      break;
    case '?': /* Print help screen */
      print_help();
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
	offsetn = plot_dataview(dv, maxval, 1.0);
	cpgclos();
	cpgslct(xid);
	filename[len] = '\0';
	printf("  Wrote the plot to the file '%s'.\n", filename);
      }
      break;
    case 'Q': /* Quit */
    case 'q':
      printf("  Quitting...\n");
      free(dv);
      cpgclos();
      break;
    default:
      printf("  Unrecognized option '%c'.\n", inchar);
      break;
    }
  } while (inchar != 'Q' && inchar != 'q');

  free_datapart(lodp);
  fclose(datfile);
  printf("Done\n\n");
  return 0;
}
