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
   *   Plot our results
   */

  prepfold_plot(&search, 1);

  /* Free our memory  */

  delete_prepfoldinfo(&search);
  return (0);
}
