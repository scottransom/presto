#include <string.h>
#include "chkio.h"

#define DEFAULT_MAXLEN 2000000000

typedef struct MULTIFILE {
  int numfiles;         /* Number of files to use together           */
  int currentfile;      /* Index of the current file (0, numfiles-1) */
  long long currentpos; /* Position in the current file              */
  long long length;     /* Length of total unified file (bytes)      */
  long long position;   /* Position in the full multifile (bytes)    */
  long long maxfilelen; /* Maximum length (bytes) of each file       */
  long long *filelens;  /* Array of lengths in bytes of each file    */
  char **filenames;     /* Array of the full paths of the files      */
  char mode[4];         /* Mode that the current file was opened     */
  FILE **fileptrs;      /* Array of FILE pointers for the files      */
} multifile;

void print_multifile(multifile *mfile, int full);
/* Print a multifile structure.  If 'full' is true, */
/* print everything.  Used for debugging.           */     

multifile *fopen_multifile(int numfiles, char **filenames, char *mode, 
			   long long maxlen);
/* Open a multifile for use and return the multifile structure.  */
/*  'numfiles' is the number of files in the multifile.          */
/*  'filenames' is an array of the names of the component files  */
/*  'mode' is the method of opening the file (binary is assumed) */
/* 	"r" : read only, do not create (truncate) the files      */
/*      "r+": read+write, do not create (truncate) the files     */
/*      "w" : read+write, create (truncate) the files            */
/*      "a" : read+write, create files or open at end of files   */
/*  'maxlen' is the maximum length in bytes of each file.  This  */
/*      number is only used if a file is opened fo writing.  The */
/*      default value is DEFAULT_MAXLEN.  If you want to use     */
/*      the default, simply set 'maxlen' to 0.                   */

int fclose_multifile(multifile *mfile);
/* Close a multifile and free its resources. */

int fread_multifile(void *data, size_t type, size_t number, 
		    multifile *mfile);
/* Read binary data from a multifile.                        */
/*   'data' is an array of the proper type to store the data */
/*   'type' is the size of each nugget of data to read       */
/*   'number' is the number of nuggets to read               */
/*   'mfile' is a pointer to a valid multifile structure     */

int fwrite_multifile(void *data, size_t type, size_t number, 
		     multifile *mfile);
/* Write binary data to a multifile.                         */
/*   'data' is an array of data to write                     */
/*   'type' is the size of each nugget of data to write      */
/*   'number' is the number of nuggets to write              */
/*   'mfile' is a pointer to a valid multifile structure     */

int fseek_multifile(multifile *mfile, long long offset, int whence);
/* Set the multifile position indicator as in fseek().   */
/*   'offset' is the file offset in bytes.               */
/*   'whence' is either SEEK_SET, SEEK_CUR, or SEEK_END. */
/*   Note:  Return is 0 for success, -1 for failure      */

void rewind_multifile(multifile *mfile);
/* Rewind a multifile. */

long long ftell_multifile(multifile *mfile);
/* Report the current position of a muiltifile. */
