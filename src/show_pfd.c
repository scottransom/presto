#include "prepfold.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* 
 * The main program 
 */

int main(int argc, char *argv[])
{
  prepfoldinfo search;

  if (argc==1 || argc > 2){
    printf("\nusage:  show_pfd pfdfilename\n\n");
    exit(0);
  }

  /*
   *   Read the raw prepfoldinfo structure
   */

  read_prepfoldinfo(&search, argv[1]);

  /*
   *   Print the main prepfoldinfo structure values
   */

  print_prepfoldinfo(&search);

  /* Switch to portrait mode */

  if (0){
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

  /*
   *   Plot our results
   */

  prepfold_plot(&search, 1);

  /* Free our memory  */

  delete_prepfoldinfo(&search);
  return (0);
}

