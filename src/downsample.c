#include "presto.h"
#include "downsample_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(int argc, char *argv[])
{
  FILE *infile, *outfile;
  int ii, jj, bufflen=10000, numread;
  long long N=0;
  float *inbuffer, *outbuffer;
  char *rootfilenm, *outname;
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
  printf("     Time Series Downsampling Routine\n");
  printf("               Sept, 2002\n\n");

  {
    int hassuffix=0;
    char *suffix;
    
    hassuffix = split_root_suffix(cmd->argv[0], &rootfilenm, &suffix);
    if (hassuffix){
      if (strcmp(suffix, "dat")!=0){
        printf("\nInput file ('%s') must be a time series ('.dat')!\n\n",
               cmd->argv[0]);
        free(suffix);
        exit(0);
      }
      free(suffix);
    } else {
      printf("\nInput file ('%s') must be a time series ('.dat')!\n\n",
             cmd->argv[0]);
      exit(0);
    }
    if (cmd->outfileP){
      outname = cmd->outfile;
    } else {
      outname = (char *)calloc(strlen(rootfilenm)+10, sizeof(char));
      sprintf(outname, "%s_D%d.dat", rootfilenm, cmd->factor);
    }
  }

  /* Read the info file */
 
  readinf(&idata, rootfilenm);
  if (idata.object) {
    printf("Downsampling %s data from '%s'.\n\n",
           remove_whitespace(idata.object), cmd->argv[0]);
  } else {
    printf("Downsampling  data from '%s'.\n\n", cmd->argv[0]);
  }

  /* Open files and create arrays */

  infile = chkfopen(argv[1], "rb");
  outfile = chkfopen(outname, "wb");
  inbuffer = gen_fvect(bufflen*cmd->factor);
  outbuffer = gen_fvect(bufflen);

  /* Read and downsample */

  while ((numread=chkfread(inbuffer, sizeof(float), bufflen*cmd->factor, infile))){
    for (ii=0; ii<numread/cmd->factor; ii++){
      outbuffer[ii] = 0; 
      for (jj=0; jj<cmd->factor; jj++)
	outbuffer[ii] += inbuffer[cmd->factor*ii+jj];
    }
    chkfwrite(outbuffer, sizeof(float), numread/cmd->factor, outfile);
    N += numread/cmd->factor;
  }
  printf("done.  Wrote %lld points.\n\n", N);

  /* Write the new info file */

  idata.dt = idata.dt * cmd->factor;
  idata.numonoff = 0;
  idata.N = (double) N;
  strncpy(idata.name, outname, strlen(outname)-4);
  idata.name[strlen(outname)-4] = '\0';
  writeinf(&idata);

  fclose(infile);
  fclose(outfile);
  free(inbuffer);
  free(outbuffer);
  free(rootfilenm);
  if (!cmd->outfileP)
    free(outname);
  exit(0);
}

