#include <glib.h>
#include "plot2d.h"
#include "presto.h"
#include "zapbirds_cmd.h"

#define NUMBETWEEN 4
#define FFTLEN     65536
#define NUMTOGET   10000
#define NUMTOSHOW  NUMBETWEEN*NUMTOGET

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
 
static FILE *fftfile;
static int khw;
static fcomplex *kernel;
static double T, dr;

typedef struct birdie {
  double freq;
  double width;
} birdie;


static birdie *birdie_create(double lofreq, double hifreq)
{
  birdie *obj;
 
  obj = (birdie *)malloc(sizeof(birdie));
  obj->freq = 0.5 * (hifreq + lofreq);
  obj->width = hifreq - lofreq;
  return obj;
}


static void birdie_free(gpointer data, gpointer user_data)
{
  free((birdie *)data);
}


static int birdie_compare(gconstpointer ca, gconstpointer cb)
/* Sorts from high to low */
{
  birdie *a, *b;
 
  a = (birdie *) ca;
  b = (birdie *) cb;
  return (a->freq > b->freq) - (a->freq < b->freq);
}


static void birdie_print(gpointer data, gpointer user_data)
{
  FILE *file=(FILE *)user_data;
  birdie *obj=(birdie *)data;
  
  fprintf(file, "%17.14g  %17.14g\n", 
	  obj->freq, obj->width);
}


static fcomplex *get_rawbins(FILE *fftfile, double bin, 
			     int numtoget, float *median, int *lobin)
{
  int ii;
  float *powers;
  fcomplex *result;

  *lobin = (int) bin - numtoget / 2;
  result = read_fcomplex_file(fftfile, *lobin, numtoget);

  /* Calculate the median power */
  powers = gen_fvect(numtoget);
  for (ii=0; ii<numtoget; ii++)
    powers[ii] = POWER(result[ii].r, result[ii].i);
  *median = selectkth(numtoget/2, numtoget, powers);  
  free(powers);

  return result;
}


static void process_bird(double basebin, int harm,
			 double *lofreq, double *hifreq)
{
  int lobin, plotlobin, plotnumbins=1024, not_done_yet=1;
  int plotoffset, inchar;
  float median, xx[2], yy[2], inx, iny;
  double bin, firstbin, pred_freq, average, tmpfreq;
  double firstfreq=0.0, lastfreq=0.0;
  fcomplex *data, *result;

  bin = basebin * harm;
  pred_freq = bin / T;
  xx[0] = xx[1] = pred_freq;
  data *get_rawbins(fftfile, bin, NUMTOGET, &median, &lobin);
  average = median / -log(0.5);
  result = gen_fcomplex(FFTLEN);
  corr_complex(data, NUMTOGET, RAW, kernel, FFTLEN, FFT, \
	       result, FFTLEN, khw, NUMBETWEEN, khw, CORR);
  firstbin = lobin + 2 * NUMBETWEEN * khw;
  do {
    plotoffset = (((int) bin - plotnumbins/2) - firstbin) * NUMBETWEEN;
    if (plotoffset < 0)
      printf("plotoffset < 0!!!\n");
    firstfreq = (firstbin+dr*plotoffset)/T;
    lastfreq = (firstbin+dr*(plotoffset+plotnumbins))/T;
    plot_spectrum(result+plotoffset, plotnumbins,
		  firstbin+dr*plotoffset, dr, T, average);
    /* Plot the marker lines */
    cpgsci(4); /* blue */
    xx[0] = firstfreq; xx[1] = lastfreq;
    yy[0] = yy[1] = 1.0;
    cpgline(2, xx, yy); /* Average power level */
    xx[0] = xx[1] = pred_freq;
    yy[0] = 0.0; yy[1] = 1.0e8;
    cpgline(2, xx, yy); /* Predicted freq */
    cpgsci(2); /* red */
    cpgswin(0.0, 1.0, 0.0, 1.0);
    yy[0] = 0.0; yy[1] = 1.0;
    cpgcurs(&inx, &iny, &inchar);
    switch (inchar){
    case ' ':
    case 'A':
    case 'a':
      printf("Add a point\n");
      xx[0] = xx[1] = inx;
      if (firstfreq)
	*lofreq = inx * (lastfreq - firstfreq) + firstfreq;
      else
	*hifreq = inx * (lastfreq - firstfreq) + firstfreq;
      cpgline(2, xx, yy); /* Predicted freq */
      break;
    case 'I': /* Zoom in */
    case 'i':
      plotnumbins /= 2;
      if (plotnumbins <= 8)
	plotnumbins = 8;
      break;
    case 'O': /* Zoom out */
    case 'o':
      plotnumbins *= 2;
      if (plotnumbins > NUMTOSHOW * NUMBETWEEN)
	plotnumbins = NUMTOSHOW * NUMBETWEEN;
      break;
    case 'C': /* Clear/Delete the points */
    case 'c':
    case 'D':
    case 'd':
      *lofreq = *hifreq = 0.0;
      continue;
    case 'Q': /* Quit/Next birdie */
    case 'q':
    case 'N':
    case 'n':
      *lofreq = *hifreq = 0.0;
      free(data);
      free(result);
      return;
    }
    cpgsci(1); /* foreground */
    if (*lofreq && *hifreq)
      not_done_yet = 0;
  } while (not_done_yet);
  if (*hifreq < *lofreq){
    tmpfreq = *lofreq;
    *lofreq = *hifreq;
    *hifreq = tmpfreq;
  }
  free(data);
  free(result);
}


int main(int argc, char *argv[])
{
  FILE *infile, *outfile;
  int ii, jj, *bird_basebins, *bird_numharms, numbirds;
  char *rootfilenm;
  GSList *zapped=NULL;
  infodata idata;
  Cmdline *cmd;
 
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
  printf("     Interactive/Automatic Bug Zapping Program\n");
  printf("              by Scott M. Ransom\n");
  printf("                 January, 2001\n\n");

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
    printf("Examining %s data from '%s'.\n\n",
           remove_whitespace(idata.object), cmd->argv[0]);
  } else {
    printf("Examining data from '%s'.\n\n", cmd->argv[0]);
  }
  T = idata.dt * idata.N;
  dr = 1.0 / numbewteen;
  fftfile = chkfopen(cmd->argv[0], "rb");

  /* Read the Standard bird list */

  numbirds = get_std_birds(cmd->inzapfile, T, cmd->baryv,
			   &bird_basebins, &bird_numharms);

  /* Create our correlation kernel */

  {
    int numkern;
    fcomplex *resp;

    khw = r_resp_halfwidth(LOWACC);
    numkern = 2 * NUMBETWEEN * khw;
    resp = gen_r_response(0.0, NUMBETWEEN, numkern);
    kernel = gen_cvect(FFTLEN);
    place_complex_kernel(resp, numkern, kernel, FFTLEN);
    COMPLEXFFT(kernel, FFTLEN, -1);
    free(resp);
  }

  /* Loop over the birdies */

  for (ii=0; ii<numbirds; ii++);
  process_bird(double basebin, int harm,
			 double *lofreq, double *hifreq)

  free(rootfilenm);
  free(kernel);
}
