#include "prepfold.h"
#include "show_pfd_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* 
 * The main program 
 */

int main(int argc, char *argv[])
{
  prepfoldinfo search;
  Cmdline *cmd;
  plotflags flags;

  /* Call usage() if we have no command line arguments */

  if (argc == 1) {
    Program = argv[0];
    usage();
    exit(0);
  }
  /* Parse the command line using the excellent program Clig */

  cmd = parseCmdline(argc, argv);
  flags.events = cmd->eventsP;
  flags.scaleparts = cmd->scalepartsP;
  flags.justprofs = cmd->justprofsP;
  flags.allgrey = cmd->allgreyP;

  /*
   *   Read the raw prepfoldinfo structure
   */

  read_prepfoldinfo(&search, argv[1]);

  /*
   *   Print the main prepfoldinfo structure values
   */

  print_prepfoldinfo(&search);

  /* Switch to portrait mode */

  if (cmd->portraitP){
    int goodlen;
    char *substr, *tmpdev;

    substr = strstr(search.pgdev, "/CPS");
    goodlen = substr - search.pgdev;
    *substr = '\0';
    tmpdev = calloc(goodlen+6, sizeof(char));
    sprintf(tmpdev, "%s/VCPS", search.pgdev);
    free(search.pgdev);
    search.pgdev = calloc(goodlen+6, sizeof(char));
    strncpy(search.pgdev, tmpdev, strlen(tmpdev));
    free(tmpdev);
    printf("New device is '%s'\n", search.pgdev);
  }

  if (0){
    int goodlen;
    char *substr, *tmpdev;

    substr = strstr(search.pgdev, "ps/CPS");
    goodlen = substr - search.pgdev;
    *substr = '\0';
    tmpdev = calloc(goodlen+9, sizeof(char));
    strncpy(tmpdev, search.pgdev, goodlen);   
    strcpy(tmpdev+goodlen, "png/TPNG");   
    free(search.pgdev);
    search.pgdev = calloc(goodlen+9, sizeof(char));
    strncpy(search.pgdev, tmpdev, strlen(tmpdev));
    free(tmpdev);
    printf("New device is '%s'\n", search.pgdev);
  }

  /*
   *   Plot our results
   */

  if (cmd->xwinP)
    prepfold_plot(&search, &flags, 1);
  else
    prepfold_plot(&search, &flags, 0);
  

  /* Free our memory  */

  delete_prepfoldinfo(&search);
  return (0);
}

