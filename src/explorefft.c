#include "presto.h"
#include "cpgplot.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
 
/* Note:  zoomlevel is simply (LOGMAXDISPLAYNUM-Log_2(numbins)) */
#define LOGNUMLOCALBINS   7  /* 128: Number of bins to show for each harmonic */
#define LOGMAXDISPLAYNUM  10 /* 1024: Maximum number of points to display at once */
#define LOGLOCALCHUNK     4  /* 16: Chunk size for polynomial fit */
#define LOGBIGGESTCHUNK   22 /* 4M points */
#define LOGINITIALNUMBINS 16 /* 65536: The initial number of bins to plot */
#define NUMLOCALBINS   (1<<LOGNUMLOCALBINS)
#define MAXDISPLAYNUM  (1<<LOGMAXDISPLAYNUM)
#define LOCALCHUNK     (1<<LOGLOCALCHUNK)
#define BIGGESTCHUNK   (1<<LOGBIGGESTCHUNK)
#define INITIALNUMBINS (1<<LOGINITIALNUMBINS)

static long long N;    /* Number of points in the original time series */
static long long Nfft; /* Number of bins in the FFT */
static double T;    /* The time duration of FFT */
static float r0;    /* The value of the zeroth Fourier freq */
static infodata idata;
static FILE *fftfile;

typedef struct fftpart {
  int rlo;          /* Lowest Fourier freq displayed */
  int numamps;      /* Number of raw amplitudes */
  float maxrawpow;  /* The highest raw power present */
  float *rawpowers; /* The raw powers */
  float *medians;   /* The local median values (chunks of size LOCALCHUNK bins) */
  float *normvals;  /* The values to use for normalization (default is median/-log(0.5)) */
  fcomplex *amps;   /* Raw FFT amplitudes    */
} fftpart;

typedef struct fftview {
  double dr;        /* Fourier frequency stepsize (2.0**(-zoomlevel)) */
  double centerr;   /* The center Fourier freq to plot */
  int lor;          /* The lowest Fourier freq to plot */
  int zoomlevel;    /* Positive = zoomed in, Negative = zoomed out */
  int numbins;      /* The number of full bins from low to high to display */
  float maxpow;     /* The maximum normalized power in the view */
  float powers[MAXDISPLAYNUM];  /* The normalized powers to display */
  double rs[MAXDISPLAYNUM];     /* The corresponding Fourier freqs */
} fftview;


static int floor_log2(int nn)
{
  int ii=0;
  
  if (nn <= 0)
    return -1;
  while (nn >> 1){
    ii++;
    nn >>= 1;
  }
  return ii;
}


static double plot_fftview(fftview *fv, float maxpow)
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
  cpgsch(1);
  cpgsci(1);

  if (maxpow==0.0) /* Autoscale for the maximum value */
    maxpow = 1.1 * fv->maxpow;

  lof = fv->lor / T;
  hif = (fv->lor + fv->dr * MAXDISPLAYNUM) / T;
  offsetf = 0.0;

  /* Period Labels */
  if (fv->zoomlevel >= 0 && fv->lor > 0.1 * Nfft){
    double lop, hip, offsetp=0.0;
    lop = 1.0/lof;
    hip = 1.0/hif;
    offsetp = 0.0;

    
    if ((lop-hip)/hip < 0.001){
      int numchar;
      char label[50];
      
      offsetp = 0.5*(hip+lop);
      numchar = snprintf(label, 50, "Period - %.15g (s)", offsetp);
      cpgmtxt("T", 3, 0.5, 0.5, label);
    } else {
      cpgmtxt("T", 3, 0.5, 0.5, "Period (s)");
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
    cpgmtxt("B", 3, 0.5, 0.5, label);
  } else {
    cpgmtxt("B", 3, 0.5, 0.5, "Frequency (Hz)");
  }
  cpgswin(lof-offsetf, hif-offsetf, 0.0, maxpow);
  if (fv->zoomlevel >= 0 && fv->lor > 0.1 * Nfft)
    cpgbox("BINST", 0.0, 0, "BCNST", 0.0, 0);
  else
    cpgbox("BCINST", 0.0, 0, "BCNST", 0.0, 0);

  /* Plot the spectrum */

  freqs = gen_fvect(MAXDISPLAYNUM);
  for (ii=0; ii<MAXDISPLAYNUM; ii++)
    freqs[ii] = fv->rs[ii]/T - offsetf;
  if (fv->zoomlevel > 0){ /* Magnified power spectrum */
    cpgline(MAXDISPLAYNUM, freqs, fv->powers);
  } else { /* Down-sampled power spectrum */
    for (ii=0; ii<MAXDISPLAYNUM; ii++){
      cpgmove(freqs[ii], 0.0);
      cpgdraw(freqs[ii], fv->powers[ii]);
    }
  }
  free(freqs);
  cpgmtxt("L", 3, 0.5, 0.5, "Normalized Power");
  cpgebuf();
  cpgunsa();
  return offsetf;
}


static fftview *get_fftview(double centerr, int zoomlevel, fftpart *fp)
{
  int ii;
  double split=1.0/(double)LOCALCHUNK;
  float powargr, powargi;
  fftview *fv;

  fv = (fftview *)malloc(sizeof(fftview));
  fv->zoomlevel = zoomlevel;
  fv->centerr = centerr;
  if (zoomlevel > 0){ /* Magnified power spectrum */
    int numbetween, nextbin, index;
    fcomplex *interp;

    numbetween = (1 << zoomlevel);
    fv->numbins = MAXDISPLAYNUM / numbetween;
    fv->dr = 1.0 / (double) numbetween;
    fv->lor = (int) floor(centerr - 0.5 * fv->numbins);
    /* fv->lor = floor(numbewteen * fv->lor + 0.5) / numbetween; */
    interp = corr_rz_interp(fp->amps, fp->numamps, numbetween,
			    fv->lor - fp->rlo, 0.0, MAXDISPLAYNUM*2,
			    LOWACC, &nextbin);
    for (ii=0; ii<MAXDISPLAYNUM; ii++){
      fv->rs[ii] = fv->lor + ii * fv->dr;
      index = (int) ((fv->rs[ii] - fp->rlo) * split + 0.5);
      fv->powers[ii] = POWER(interp[ii].r, interp[ii].i) * fp->normvals[index];
    }
    free(interp);
  } else { /* Down-sampled power spectrum */
    int jj, powindex, normindex, binstocombine;
    float *tmprawpwrs, maxpow;

    binstocombine = (1 << abs(zoomlevel));
    fv->numbins = MAXDISPLAYNUM * binstocombine;
    fv->dr = (double) binstocombine;
    fv->lor = (int) floor(centerr - 0.5 * fv->numbins);
    tmprawpwrs = gen_fvect(fv->numbins);
    for (ii=0; ii<fv->numbins; ii++){
      powindex = (int) (fv->lor - fp->rlo + ii + 0.5);
      normindex = (int) (powindex * split);
      tmprawpwrs[ii] = fp->rawpowers[powindex] * fp->normvals[normindex];
    }
    powindex = 0;
    for (ii=0; ii<MAXDISPLAYNUM; ii++){
      maxpow = 0.0;
      for (jj=0; jj<binstocombine; jj++, powindex++)
	if (tmprawpwrs[powindex] > maxpow) maxpow = tmprawpwrs[powindex];
      fv->rs[ii] = fv->lor + ii * fv->dr;
      fv->powers[ii] = maxpow;
    }
  }
  fv->maxpow = 0.0;
  for (ii=0; ii<MAXDISPLAYNUM; ii++)
    if (fv->powers[ii] > fv->maxpow) fv->maxpow = fv->powers[ii];
  return fv;
}


static fftpart *get_fftpart(int rlo, int numr)
{
  int ii, jj, index;
  float powargr, powargi, tmppwr, chunk[LOCALCHUNK];
  fftpart *fp;

  fp = (fftpart *)malloc(sizeof(fftpart));
  fp->rlo = rlo;
  fp->numamps = numr;
  fp->amps = read_fcomplex_file(fftfile, rlo, fp->numamps);
  if (rlo==0){
    fp->amps[0].r = 1.0;
    fp->amps[0].i = 0.0;
  }
  fp->rawpowers = gen_fvect(fp->numamps);
  fp->medians = gen_fvect(fp->numamps/LOCALCHUNK);
  fp->normvals = gen_fvect(fp->numamps/LOCALCHUNK);
  fp->maxrawpow = 0.0;
  for (ii=0; ii<fp->numamps/LOCALCHUNK; ii++){
    index = ii * LOCALCHUNK;
    for (jj=0; jj<LOCALCHUNK; jj++, index++){
      tmppwr = POWER(fp->amps[index].r, fp->amps[index].i);
      if (tmppwr > fp->maxrawpow)
	fp->maxrawpow = tmppwr;
      fp->rawpowers[index] = tmppwr;
      chunk[jj] = tmppwr;
    }
    fp->medians[ii] = median(chunk, LOCALCHUNK);
    fp->normvals[ii] = 1.0 / (1.4426950408889634 * fp->medians[ii]);
  }
  return fp;
}

static void free_fftpart(fftpart *fp)
{
  free(fp->normvals);
  free(fp->medians);
  free(fp->rawpowers);
  free(fp->amps);
  free(fp);
}


int main(int argc, char *argv[])
{
  double centerr;
  int numamps, zoomlevel;
  char *rootfilenm, inchar;
  fftpart *lofp;
 
  printf("\n\n");
  printf("      Interactive FFT Explorer\n");
  printf("         by Scott M. Ransom\n");
  printf("            October, 2001\n\n");

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
  fftfile = chkfopen(argv[1], "rb");
  Nfft = chkfilelen(fftfile, sizeof(fcomplex));

  /* Get and plot the initial data */
  
  numamps = (Nfft > BIGGESTCHUNK) ? (int) BIGGESTCHUNK : (int) Nfft;
  lofp = get_fftpart(0, numamps);
  r0 = lofp->amps[0].r;
  centerr = 0.5 * INITIALNUMBINS;
  zoomlevel = LOGMAXDISPLAYNUM - LOGINITIALNUMBINS;

  /* Prep the XWIN device for PGPLOT */

  if (cpgopen("/XWIN") <= 0)
    exit(EXIT_FAILURE);
  /* cpgpap(11.0, 8.5/11.0); */
  cpgpap(6.0, 8.5/11.0);
  cpgask(0);

  do {
    double offsetf=0.0;
    float inx, iny, maxpow=0.0;
    fftview *fv;
    
    fv = get_fftview(centerr, zoomlevel, lofp);
    cpgpage();
    offsetf = plot_fftview(fv, maxpow);

    cpgcurs(&inx, &iny, &inchar);
    printf("You pressed '%c'\n", inchar);

    switch (inchar){
    case 'A': /* Zoom in */
    case 'a':
    case 'I':
    case 'i':
      printf("  Zooming in  (zoomlevel = %d)...\n", zoomlevel);
      zoomlevel++;
      free(fv);
      fv = get_fftview(centerr, zoomlevel, lofp);
      cpgpage();
      offsetf = plot_fftview(fv, maxpow);
      break;
    case 'X': /* Zoom out */
    case 'x':
    case 'O':
    case 'o':
      printf("  Zooming out  (zoomlevel = %d)...\n", zoomlevel);
      zoomlevel--;
      free(fv);
      fv = get_fftview(centerr, zoomlevel, lofp);
      cpgpage();
      offsetf = plot_fftview(fv, maxpow);
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

  free_fftpart(lofp);
  fclose(fftfile);
  printf("Done\n\n");
  return 0;
}
