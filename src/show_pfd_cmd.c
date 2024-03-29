/*****
  command line parser -- generated by clig
  (http://wsd.iitb.fhg.de/~kir/clighome/)

  The command line parser `clig':
  (C) 1995-2004 Harald Kirsch (clig@geggus.net)
*****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "show_pfd_cmd.h"

char *Program;

/*@-null*/

static Cmdline cmd = {
  /***** -noxwin: Do not show the result plots on-screen, only make postscript files */
  /* noxwinP = */ 0,
  /***** -showfold: Use the input fold paramters (i.e. not the optimized values) when showing the plot */
  /* showfoldP = */ 0,
  /***** -scaleparts: Scale the part profiles independently */
  /* scalepartsP = */ 0,
  /***** -allgrey: Make all the images greyscale instead of color */
  /* allgreyP = */ 0,
  /***** -justprofs: Only output the profile portions of the plot */
  /* justprofsP = */ 0,
  /***** -portrait: Orient the output in portrait mode (for -justprofs) */
  /* portraitP = */ 0,
  /***** -events: The folded data were events instead of samples or bins */
  /* eventsP = */ 0,
  /***** -infoonly: Display the pfd info and exit without generating plots. */
  /* infoonlyP = */ 0,
  /***** -fixchi: Adjust the reduced chi^2 values so that off-pulse reduced chi^2 = 1 */
  /* fixchiP = */ 0,
  /***** -samples: Treat the data as samples and not as finite-duration integrated data */
  /* samplesP = */ 0,
  /***** -normalize: Normalize stats for each fold (i.e. to bandpass flatten subbands) */
  /* normalizeP = */ 0,
  /***** -killsubs: Comma separated string (no spaces!) of subbands to explicitly remove from analysis (i.e. zero out).  Ranges are specified by min:max[:step] */
  /* killsubsstrP = */ 0,
  /* killsubsstr = */ (char*)0,
  /* killsubsstrC = */ 0,
  /***** -killparts: Comma separated string (no spaces!) of intervals to explicitly remove from analysis (i.e. zero-out).  Ranges are specified by min:max[:step] */
  /* killpartsstrP = */ 0,
  /* killpartsstr = */ (char*)0,
  /* killpartsstrC = */ 0,
  /***** uninterpreted rest of command line */
  /* argc = */ 0,
  /* argv = */ (char**)0,
  /***** the original command line concatenated */
  /* full_cmd_line = */ NULL
};

/*@=null*/

/***** let LCLint run more smoothly */
/*@-predboolothers*/
/*@-boolops*/


/******************************************************************/
/*****
 This is a bit tricky. We want to make a difference between overflow
 and underflow and we want to allow v==Inf or v==-Inf but not
 v>FLT_MAX. 

 We don't use fabs to avoid linkage with -lm.
*****/
static void
checkFloatConversion(double v, char *option, char *arg)
{
  char *err = NULL;

  if( (errno==ERANGE && v!=0.0) /* even double overflowed */
      || (v<HUGE_VAL && v>-HUGE_VAL && (v<0.0?-v:v)>(double)FLT_MAX) ) {
    err = "large";
  } else if( (errno==ERANGE && v==0.0) 
	     || (v!=0.0 && (v<0.0?-v:v)<(double)FLT_MIN) ) {
    err = "small";
  }
  if( err ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of option `%s' to %s to represent\n",
	    Program, arg, option, err);
    exit(EXIT_FAILURE);
  }
}

int
getIntOpt(int argc, char **argv, int i, int *value, int force)
{
  char *end;
  long v;

  if( ++i>=argc ) goto nothingFound;

  errno = 0;
  v = strtol(argv[i], &end, 0);

  /***** check for conversion error */
  if( end==argv[i] ) goto nothingFound;

  /***** check for surplus non-whitespace */
  while( isspace((int) *end) ) end+=1;
  if( *end ) goto nothingFound;

  /***** check if it fits into an int */
  if( errno==ERANGE || v>(long)INT_MAX || v<(long)INT_MIN ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of option `%s' to large to represent\n",
	    Program, argv[i], argv[i-1]);
    exit(EXIT_FAILURE);
  }
  *value = (int)v;

  return i;

nothingFound:
  if( !force ) return i-1;

  fprintf(stderr, 
	  "%s: missing or malformed integer value after option `%s'\n",
	  Program, argv[i-1]);
    exit(EXIT_FAILURE);
}
/**********************************************************************/

int
getIntOpts(int argc, char **argv, int i, 
	   int **values,
	   int cmin, int cmax)
/*****
  We want to find at least cmin values and at most cmax values.
  cmax==-1 then means infinitely many are allowed.
*****/
{
  int alloced, used;
  char *end;
  long v;
  if( i+cmin >= argc ) {
    fprintf(stderr, 
	    "%s: option `%s' wants at least %d parameters\n",
	    Program, argv[i], cmin);
    exit(EXIT_FAILURE);
  }

  /***** 
    alloc a bit more than cmin values. It does not hurt to have room
    for a bit more values than cmax.
  *****/
  alloced = cmin + 4;
  *values = (int*)calloc((size_t)alloced, sizeof(int));
  if( ! *values ) {
outMem:
    fprintf(stderr, 
	    "%s: out of memory while parsing option `%s'\n",
	    Program, argv[i]);
    exit(EXIT_FAILURE);
  }

  for(used=0; (cmax==-1 || used<cmax) && used+i+1<argc; used++) {
    if( used==alloced ) {
      alloced += 8;
      *values = (int *) realloc(*values, alloced*sizeof(int));
      if( !*values ) goto outMem;
    }

    errno = 0;
    v = strtol(argv[used+i+1], &end, 0);

    /***** check for conversion error */
    if( end==argv[used+i+1] ) break;

    /***** check for surplus non-whitespace */
    while( isspace((int) *end) ) end+=1;
    if( *end ) break;

    /***** check for overflow */
    if( errno==ERANGE || v>(long)INT_MAX || v<(long)INT_MIN ) {
      fprintf(stderr, 
	      "%s: parameter `%s' of option `%s' to large to represent\n",
	      Program, argv[i+used+1], argv[i]);
      exit(EXIT_FAILURE);
    }

    (*values)[used] = (int)v;

  }
    
  if( used<cmin ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of `%s' should be an "
	    "integer value\n",
	    Program, argv[i+used+1], argv[i]);
    exit(EXIT_FAILURE);
  }

  return i+used;
}
/**********************************************************************/

int
getLongOpt(int argc, char **argv, int i, long *value, int force)
{
  char *end;

  if( ++i>=argc ) goto nothingFound;

  errno = 0;
  *value = strtol(argv[i], &end, 0);

  /***** check for conversion error */
  if( end==argv[i] ) goto nothingFound;

  /***** check for surplus non-whitespace */
  while( isspace((int) *end) ) end+=1;
  if( *end ) goto nothingFound;

  /***** check for overflow */
  if( errno==ERANGE ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of option `%s' to large to represent\n",
	    Program, argv[i], argv[i-1]);
    exit(EXIT_FAILURE);
  }
  return i;

nothingFound:
  /***** !force means: this parameter may be missing.*/
  if( !force ) return i-1;

  fprintf(stderr, 
	  "%s: missing or malformed value after option `%s'\n",
	  Program, argv[i-1]);
    exit(EXIT_FAILURE);
}
/**********************************************************************/

int
getLongOpts(int argc, char **argv, int i, 
	    long **values,
	    int cmin, int cmax)
/*****
  We want to find at least cmin values and at most cmax values.
  cmax==-1 then means infinitely many are allowed.
*****/
{
  int alloced, used;
  char *end;

  if( i+cmin >= argc ) {
    fprintf(stderr, 
	    "%s: option `%s' wants at least %d parameters\n",
	    Program, argv[i], cmin);
    exit(EXIT_FAILURE);
  }

  /***** 
    alloc a bit more than cmin values. It does not hurt to have room
    for a bit more values than cmax.
  *****/
  alloced = cmin + 4;
  *values = (long int *)calloc((size_t)alloced, sizeof(long));
  if( ! *values ) {
outMem:
    fprintf(stderr, 
	    "%s: out of memory while parsing option `%s'\n",
	    Program, argv[i]);
    exit(EXIT_FAILURE);
  }

  for(used=0; (cmax==-1 || used<cmax) && used+i+1<argc; used++) {
    if( used==alloced ) {
      alloced += 8;
      *values = (long int*) realloc(*values, alloced*sizeof(long));
      if( !*values ) goto outMem;
    }

    errno = 0;
    (*values)[used] = strtol(argv[used+i+1], &end, 0);

    /***** check for conversion error */
    if( end==argv[used+i+1] ) break;

    /***** check for surplus non-whitespace */
    while( isspace((int) *end) ) end+=1; 
    if( *end ) break;

    /***** check for overflow */
    if( errno==ERANGE ) {
      fprintf(stderr, 
	      "%s: parameter `%s' of option `%s' to large to represent\n",
	      Program, argv[i+used+1], argv[i]);
      exit(EXIT_FAILURE);
    }

  }
    
  if( used<cmin ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of `%s' should be an "
	    "integer value\n",
	    Program, argv[i+used+1], argv[i]);
    exit(EXIT_FAILURE);
  }

  return i+used;
}
/**********************************************************************/

int
getFloatOpt(int argc, char **argv, int i, float *value, int force)
{
  char *end;
  double v;

  if( ++i>=argc ) goto nothingFound;

  errno = 0;
  v = strtod(argv[i], &end);

  /***** check for conversion error */
  if( end==argv[i] ) goto nothingFound;

  /***** check for surplus non-whitespace */
  while( isspace((int) *end) ) end+=1;
  if( *end ) goto nothingFound;

  /***** check for overflow */
  checkFloatConversion(v, argv[i-1], argv[i]);

  *value = (float)v;

  return i;

nothingFound:
  if( !force ) return i-1;

  fprintf(stderr,
	  "%s: missing or malformed float value after option `%s'\n",
	  Program, argv[i-1]);
  exit(EXIT_FAILURE);
 
}
/**********************************************************************/

int
getFloatOpts(int argc, char **argv, int i, 
	   float **values,
	   int cmin, int cmax)
/*****
  We want to find at least cmin values and at most cmax values.
  cmax==-1 then means infinitely many are allowed.
*****/
{
  int alloced, used;
  char *end;
  double v;

  if( i+cmin >= argc ) {
    fprintf(stderr, 
	    "%s: option `%s' wants at least %d parameters\n",
	    Program, argv[i], cmin);
    exit(EXIT_FAILURE);
  }

  /***** 
    alloc a bit more than cmin values.
  *****/
  alloced = cmin + 4;
  *values = (float*)calloc((size_t)alloced, sizeof(float));
  if( ! *values ) {
outMem:
    fprintf(stderr, 
	    "%s: out of memory while parsing option `%s'\n",
	    Program, argv[i]);
    exit(EXIT_FAILURE);
  }

  for(used=0; (cmax==-1 || used<cmax) && used+i+1<argc; used++) {
    if( used==alloced ) {
      alloced += 8;
      *values = (float *) realloc(*values, alloced*sizeof(float));
      if( !*values ) goto outMem;
    }

    errno = 0;
    v = strtod(argv[used+i+1], &end);

    /***** check for conversion error */
    if( end==argv[used+i+1] ) break;

    /***** check for surplus non-whitespace */
    while( isspace((int) *end) ) end+=1;
    if( *end ) break;

    /***** check for overflow */
    checkFloatConversion(v, argv[i], argv[i+used+1]);
    
    (*values)[used] = (float)v;
  }
    
  if( used<cmin ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of `%s' should be a "
	    "floating-point value\n",
	    Program, argv[i+used+1], argv[i]);
    exit(EXIT_FAILURE);
  }

  return i+used;
}
/**********************************************************************/

int
getDoubleOpt(int argc, char **argv, int i, double *value, int force)
{
  char *end;

  if( ++i>=argc ) goto nothingFound;

  errno = 0;
  *value = strtod(argv[i], &end);

  /***** check for conversion error */
  if( end==argv[i] ) goto nothingFound;

  /***** check for surplus non-whitespace */
  while( isspace((int) *end) ) end+=1;
  if( *end ) goto nothingFound;

  /***** check for overflow */
  if( errno==ERANGE ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of option `%s' to %s to represent\n",
	    Program, argv[i], argv[i-1],
	    (*value==0.0 ? "small" : "large"));
    exit(EXIT_FAILURE);
  }

  return i;

nothingFound:
  if( !force ) return i-1;

  fprintf(stderr,
	  "%s: missing or malformed value after option `%s'\n",
	  Program, argv[i-1]);
  exit(EXIT_FAILURE);
 
}
/**********************************************************************/

int
getDoubleOpts(int argc, char **argv, int i, 
	   double **values,
	   int cmin, int cmax)
/*****
  We want to find at least cmin values and at most cmax values.
  cmax==-1 then means infinitely many are allowed.
*****/
{
  int alloced, used;
  char *end;

  if( i+cmin >= argc ) {
    fprintf(stderr, 
	    "%s: option `%s' wants at least %d parameters\n",
	    Program, argv[i], cmin);
    exit(EXIT_FAILURE);
  }

  /***** 
    alloc a bit more than cmin values.
  *****/
  alloced = cmin + 4;
  *values = (double*)calloc((size_t)alloced, sizeof(double));
  if( ! *values ) {
outMem:
    fprintf(stderr, 
	    "%s: out of memory while parsing option `%s'\n",
	    Program, argv[i]);
    exit(EXIT_FAILURE);
  }

  for(used=0; (cmax==-1 || used<cmax) && used+i+1<argc; used++) {
    if( used==alloced ) {
      alloced += 8;
      *values = (double *) realloc(*values, alloced*sizeof(double));
      if( !*values ) goto outMem;
    }

    errno = 0;
    (*values)[used] = strtod(argv[used+i+1], &end);

    /***** check for conversion error */
    if( end==argv[used+i+1] ) break;

    /***** check for surplus non-whitespace */
    while( isspace((int) *end) ) end+=1;
    if( *end ) break;

    /***** check for overflow */
    if( errno==ERANGE ) {
      fprintf(stderr, 
	      "%s: parameter `%s' of option `%s' to %s to represent\n",
	      Program, argv[i+used+1], argv[i],
	      ((*values)[used]==0.0 ? "small" : "large"));
      exit(EXIT_FAILURE);
    }

  }
    
  if( used<cmin ) {
    fprintf(stderr, 
	    "%s: parameter `%s' of `%s' should be a "
	    "double value\n",
	    Program, argv[i+used+1], argv[i]);
    exit(EXIT_FAILURE);
  }

  return i+used;
}
/**********************************************************************/

/**
  force will be set if we need at least one argument for the option.
*****/
int
getStringOpt(int argc, char **argv, int i, char **value, int force)
{
  i += 1;
  if( i>=argc ) {
    if( force ) {
      fprintf(stderr, "%s: missing string after option `%s'\n",
	      Program, argv[i-1]);
      exit(EXIT_FAILURE);
    } 
    return i-1;
  }
  
  if( !force && argv[i][0] == '-' ) return i-1;
  *value = argv[i];
  return i;
}
/**********************************************************************/

int
getStringOpts(int argc, char **argv, int i, 
	   char*  **values,
	   int cmin, int cmax)
/*****
  We want to find at least cmin values and at most cmax values.
  cmax==-1 then means infinitely many are allowed.
*****/
{
  int alloced, used;

  if( i+cmin >= argc ) {
    fprintf(stderr, 
	    "%s: option `%s' wants at least %d parameters\n",
	    Program, argv[i], cmin);
    exit(EXIT_FAILURE);
  }

  alloced = cmin + 4;
    
  *values = (char**)calloc((size_t)alloced, sizeof(char*));
  if( ! *values ) {
outMem:
    fprintf(stderr, 
	    "%s: out of memory during parsing of option `%s'\n",
	    Program, argv[i]);
    exit(EXIT_FAILURE);
  }

  for(used=0; (cmax==-1 || used<cmax) && used+i+1<argc; used++) {
    if( used==alloced ) {
      alloced += 8;
      *values = (char **)realloc(*values, alloced*sizeof(char*));
      if( !*values ) goto outMem;
    }

    if( used>=cmin && argv[used+i+1][0]=='-' ) break;
    (*values)[used] = argv[used+i+1];
  }
    
  if( used<cmin ) {
    fprintf(stderr, 
    "%s: less than %d parameters for option `%s', only %d found\n",
	    Program, cmin, argv[i], used);
    exit(EXIT_FAILURE);
  }

  return i+used;
}
/**********************************************************************/

void
checkIntLower(char *opt, int *values, int count, int max)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]<=max ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' greater than max=%d\n",
	    Program, i+1, opt, max);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkIntHigher(char *opt, int *values, int count, int min)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]>=min ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' smaller than min=%d\n",
	    Program, i+1, opt, min);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkLongLower(char *opt, long *values, int count, long max)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]<=max ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' greater than max=%ld\n",
	    Program, i+1, opt, max);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkLongHigher(char *opt, long *values, int count, long min)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]>=min ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' smaller than min=%ld\n",
	    Program, i+1, opt, min);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkFloatLower(char *opt, float *values, int count, float max)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]<=max ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' greater than max=%f\n",
	    Program, i+1, opt, max);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkFloatHigher(char *opt, float *values, int count, float min)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]>=min ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' smaller than min=%f\n",
	    Program, i+1, opt, min);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkDoubleLower(char *opt, double *values, int count, double max)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]<=max ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' greater than max=%f\n",
	    Program, i+1, opt, max);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

void
checkDoubleHigher(char *opt, double *values, int count, double min)
{
  int i;

  for(i=0; i<count; i++) {
    if( values[i]>=min ) continue;
    fprintf(stderr, 
	    "%s: parameter %d of option `%s' smaller than min=%f\n",
	    Program, i+1, opt, min);
    exit(EXIT_FAILURE);
  }
}
/**********************************************************************/

static char *
catArgv(int argc, char **argv)
{
  int i;
  size_t l;
  char *s, *t;

  for(i=0, l=0; i<argc; i++) l += (1+strlen(argv[i]));
  s = (char *)malloc(l);
  if( !s ) {
    fprintf(stderr, "%s: out of memory\n", Program);
    exit(EXIT_FAILURE);
  }
  strcpy(s, argv[0]);
  t = s;
  for(i=1; i<argc; i++) {
    t = t+strlen(t);
    *t++ = ' ';
    strcpy(t, argv[i]);
  }
  return s;
}
/**********************************************************************/

void
showOptionValues(void)
{
  int i;

  printf("Full command line is:\n`%s'\n", cmd.full_cmd_line);

  /***** -noxwin: Do not show the result plots on-screen, only make postscript files */
  if( !cmd.noxwinP ) {
    printf("-noxwin not found.\n");
  } else {
    printf("-noxwin found:\n");
  }

  /***** -showfold: Use the input fold paramters (i.e. not the optimized values) when showing the plot */
  if( !cmd.showfoldP ) {
    printf("-showfold not found.\n");
  } else {
    printf("-showfold found:\n");
  }

  /***** -scaleparts: Scale the part profiles independently */
  if( !cmd.scalepartsP ) {
    printf("-scaleparts not found.\n");
  } else {
    printf("-scaleparts found:\n");
  }

  /***** -allgrey: Make all the images greyscale instead of color */
  if( !cmd.allgreyP ) {
    printf("-allgrey not found.\n");
  } else {
    printf("-allgrey found:\n");
  }

  /***** -justprofs: Only output the profile portions of the plot */
  if( !cmd.justprofsP ) {
    printf("-justprofs not found.\n");
  } else {
    printf("-justprofs found:\n");
  }

  /***** -portrait: Orient the output in portrait mode (for -justprofs) */
  if( !cmd.portraitP ) {
    printf("-portrait not found.\n");
  } else {
    printf("-portrait found:\n");
  }

  /***** -events: The folded data were events instead of samples or bins */
  if( !cmd.eventsP ) {
    printf("-events not found.\n");
  } else {
    printf("-events found:\n");
  }

  /***** -infoonly: Display the pfd info and exit without generating plots. */
  if( !cmd.infoonlyP ) {
    printf("-infoonly not found.\n");
  } else {
    printf("-infoonly found:\n");
  }

  /***** -fixchi: Adjust the reduced chi^2 values so that off-pulse reduced chi^2 = 1 */
  if( !cmd.fixchiP ) {
    printf("-fixchi not found.\n");
  } else {
    printf("-fixchi found:\n");
  }

  /***** -samples: Treat the data as samples and not as finite-duration integrated data */
  if( !cmd.samplesP ) {
    printf("-samples not found.\n");
  } else {
    printf("-samples found:\n");
  }

  /***** -normalize: Normalize stats for each fold (i.e. to bandpass flatten subbands) */
  if( !cmd.normalizeP ) {
    printf("-normalize not found.\n");
  } else {
    printf("-normalize found:\n");
  }

  /***** -killsubs: Comma separated string (no spaces!) of subbands to explicitly remove from analysis (i.e. zero out).  Ranges are specified by min:max[:step] */
  if( !cmd.killsubsstrP ) {
    printf("-killsubs not found.\n");
  } else {
    printf("-killsubs found:\n");
    if( !cmd.killsubsstrC ) {
      printf("  no values\n");
    } else {
      printf("  value = `%s'\n", cmd.killsubsstr);
    }
  }

  /***** -killparts: Comma separated string (no spaces!) of intervals to explicitly remove from analysis (i.e. zero-out).  Ranges are specified by min:max[:step] */
  if( !cmd.killpartsstrP ) {
    printf("-killparts not found.\n");
  } else {
    printf("-killparts found:\n");
    if( !cmd.killpartsstrC ) {
      printf("  no values\n");
    } else {
      printf("  value = `%s'\n", cmd.killpartsstr);
    }
  }
  if( !cmd.argc ) {
    printf("no remaining parameters in argv\n");
  } else {
    printf("argv =");
    for(i=0; i<cmd.argc; i++) {
      printf(" `%s'", cmd.argv[i]);
    }
    printf("\n");
  }
}
/**********************************************************************/

void
usage(void)
{
  fprintf(stderr,"%s","   [-noxwin] [-showfold] [-scaleparts] [-allgrey] [-justprofs] [-portrait] [-events] [-infoonly] [-fixchi] [-samples] [-normalize] [-killsubs killsubsstr] [-killparts killpartsstr] [--] infile ...\n");
  fprintf(stderr,"%s","      Displays or regenerates the Postscript for a 'pfd' file created by prepfold.\n");
  fprintf(stderr,"%s","        -noxwin: Do not show the result plots on-screen, only make postscript files\n");
  fprintf(stderr,"%s","      -showfold: Use the input fold paramters (i.e. not the optimized values) when showing the plot\n");
  fprintf(stderr,"%s","    -scaleparts: Scale the part profiles independently\n");
  fprintf(stderr,"%s","       -allgrey: Make all the images greyscale instead of color\n");
  fprintf(stderr,"%s","     -justprofs: Only output the profile portions of the plot\n");
  fprintf(stderr,"%s","      -portrait: Orient the output in portrait mode (for -justprofs)\n");
  fprintf(stderr,"%s","        -events: The folded data were events instead of samples or bins\n");
  fprintf(stderr,"%s","      -infoonly: Display the pfd info and exit without generating plots.\n");
  fprintf(stderr,"%s","        -fixchi: Adjust the reduced chi^2 values so that off-pulse reduced chi^2 = 1\n");
  fprintf(stderr,"%s","       -samples: Treat the data as samples and not as finite-duration integrated data\n");
  fprintf(stderr,"%s","     -normalize: Normalize stats for each fold (i.e. to bandpass flatten subbands)\n");
  fprintf(stderr,"%s","      -killsubs: Comma separated string (no spaces!) of subbands to explicitly remove from analysis (i.e. zero out).  Ranges are specified by min:max[:step]\n");
  fprintf(stderr,"%s","                 1 char* value\n");
  fprintf(stderr,"%s","     -killparts: Comma separated string (no spaces!) of intervals to explicitly remove from analysis (i.e. zero-out).  Ranges are specified by min:max[:step]\n");
  fprintf(stderr,"%s","                 1 char* value\n");
  fprintf(stderr,"%s","         infile: The input 'pfd' file name.\n");
  fprintf(stderr,"%s","                 1...100 values\n");
  fprintf(stderr,"%s","  version: 04Feb23\n");
  fprintf(stderr,"%s","  ");
  exit(EXIT_FAILURE);
}
/**********************************************************************/
Cmdline *
parseCmdline(int argc, char **argv)
{
  int i;

  Program = argv[0];
  cmd.full_cmd_line = catArgv(argc, argv);
  for(i=1, cmd.argc=1; i<argc; i++) {
    if( 0==strcmp("--", argv[i]) ) {
      while( ++i<argc ) argv[cmd.argc++] = argv[i];
      continue;
    }

    if( 0==strcmp("-noxwin", argv[i]) ) {
      cmd.noxwinP = 1;
      continue;
    }

    if( 0==strcmp("-showfold", argv[i]) ) {
      cmd.showfoldP = 1;
      continue;
    }

    if( 0==strcmp("-scaleparts", argv[i]) ) {
      cmd.scalepartsP = 1;
      continue;
    }

    if( 0==strcmp("-allgrey", argv[i]) ) {
      cmd.allgreyP = 1;
      continue;
    }

    if( 0==strcmp("-justprofs", argv[i]) ) {
      cmd.justprofsP = 1;
      continue;
    }

    if( 0==strcmp("-portrait", argv[i]) ) {
      cmd.portraitP = 1;
      continue;
    }

    if( 0==strcmp("-events", argv[i]) ) {
      cmd.eventsP = 1;
      continue;
    }

    if( 0==strcmp("-infoonly", argv[i]) ) {
      cmd.infoonlyP = 1;
      continue;
    }

    if( 0==strcmp("-fixchi", argv[i]) ) {
      cmd.fixchiP = 1;
      continue;
    }

    if( 0==strcmp("-samples", argv[i]) ) {
      cmd.samplesP = 1;
      continue;
    }

    if( 0==strcmp("-normalize", argv[i]) ) {
      cmd.normalizeP = 1;
      continue;
    }

    if( 0==strcmp("-killsubs", argv[i]) ) {
      int keep = i;
      cmd.killsubsstrP = 1;
      i = getStringOpt(argc, argv, i, &cmd.killsubsstr, 1);
      cmd.killsubsstrC = i-keep;
      continue;
    }

    if( 0==strcmp("-killparts", argv[i]) ) {
      int keep = i;
      cmd.killpartsstrP = 1;
      i = getStringOpt(argc, argv, i, &cmd.killpartsstr, 1);
      cmd.killpartsstrC = i-keep;
      continue;
    }

    if( argv[i][0]=='-' ) {
      fprintf(stderr, "\n%s: unknown option `%s'\n\n",
              Program, argv[i]);
      usage();
    }
    argv[cmd.argc++] = argv[i];
  }/* for i */


  /*@-mustfree*/
  cmd.argv = argv+1;
  /*@=mustfree*/
  cmd.argc -= 1;

  if( 1>cmd.argc ) {
    fprintf(stderr, "%s: there should be at least 1 non-option argument(s)\n",
            Program);
    exit(EXIT_FAILURE);
  }
  if( 100<cmd.argc ) {
    fprintf(stderr, "%s: there should be at most 100 non-option argument(s)\n",
            Program);
    exit(EXIT_FAILURE);
  }
  /*@-compmempass*/  return &cmd;
}

