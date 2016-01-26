#include "prepfold.h"
#include "show_pfd_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

extern int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals);

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
    flags.fixchi = cmd->fixchiP;
    flags.showfold = cmd->showfoldP;
    flags.nosearch = 1;

    /*
     *   Read the raw prepfoldinfo structure
     */

    read_prepfoldinfo(&search, argv[1]);

    /*
     *   Print the main prepfoldinfo structure values
     */

    print_prepfoldinfo(&search);
    if (cmd->infoonlyP)
        return 0;
    /*
     *   Zap requested subbands or intervals
     */

    {
        int *killparts, *killsubs, ii, jj, kk, index;
        int numkillparts = 0, numkillsubs = 0;

        if (cmd->killpartsstrP) {
            killparts = ranges_to_ivect(cmd->killpartsstr, 0,
                                        search.npart - 1, &numkillparts);
            for (ii = 0; ii < numkillparts; ii++) {
                if ((killparts[ii] >= 0) && (killparts[ii] < search.npart)) {
                    index = killparts[ii] * search.proflen * search.nsub;
                    for (jj = 0; jj < search.nsub; jj++) {
                        search.stats[killparts[ii] * search.nsub + jj].prof_var =
                            0.0;
                        search.stats[killparts[ii] * search.nsub + jj].prof_avg =
                            0.0;
                        for (kk = 0; kk < search.proflen; kk++)
                            search.rawfolds[index + kk] = 0.0;
                        index += search.proflen;
                    }
                }
            }
            free(killparts);
        }
        if (cmd->killsubsstrP) {
            killsubs = ranges_to_ivect(cmd->killsubsstr, 0,
                                       search.nsub - 1, &numkillsubs);
            for (ii = 0; ii < numkillsubs; ii++) {
                if ((killsubs[ii] >= 0) && (killsubs[ii] < search.nsub)) {
                    for (jj = 0; jj < search.npart; jj++) {
                        index = search.proflen * (jj * search.nsub + killsubs[ii]);
                        search.stats[jj * search.nsub + killsubs[ii]].prof_var = 0.0;
                        search.stats[jj * search.nsub + killsubs[ii]].prof_avg = 0.0;
                        for (kk = 0; kk < search.proflen; kk++)
                            search.rawfolds[index + kk] = 0.0;
                    }
                }
            }
            free(killsubs);
        }
    }

    /* Switch to portrait mode */

    if (cmd->portraitP) {
        int goodlen;
        char *substr, *tmpdev;

        substr = strstr(search.pgdev, "/CPS");
        goodlen = substr - search.pgdev;
        *substr = '\0';
        tmpdev = calloc(goodlen + 6, sizeof(char));
        sprintf(tmpdev, "%s/VCPS", search.pgdev);
        free(search.pgdev);
        search.pgdev = calloc(goodlen + 6, sizeof(char));
        strncpy(search.pgdev, tmpdev, strlen(tmpdev));
        free(tmpdev);
        printf("New device is '%s'\n", search.pgdev);
    }

    if (0) {
        int goodlen;
        char *substr, *tmpdev;

        substr = strstr(search.pgdev, "ps/CPS");
        goodlen = substr - search.pgdev;
        *substr = '\0';
        tmpdev = calloc(goodlen + 9, sizeof(char));
        strncpy(tmpdev, search.pgdev, goodlen);
        strcpy(tmpdev + goodlen, "png/TPNG");
        free(search.pgdev);
        search.pgdev = calloc(goodlen + 9, sizeof(char));
        strncpy(search.pgdev, tmpdev, strlen(tmpdev));
        free(tmpdev);
        printf("New device is '%s'\n", search.pgdev);
    }

    /*
     *   Plot our results
     */

    prepfold_plot(&search, &flags, !cmd->noxwinP, NULL);

    /* Free our memory  */

    delete_prepfoldinfo(&search);
    return (0);
}
