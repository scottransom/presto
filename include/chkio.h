#include <sys/types.h>
#ifdef USE_PIOFS
#include <sys/limits.h>
#include <piofs/piofs_ioctl.h>
#else
#include <sys/stat.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

FILE *chkfopen(char *path, const char *mode);
/* Preform a file open with error checking.  */

size_t chkfread(void *data, size_t type, size_t number, FILE * stream);
/* Preform a file read with error checking.  */

size_t chkfwrite(void *data, size_t type, size_t number, FILE * stream);
/* Preform a file write with error checking. */

size_t chkfseek(FILE * stream, long offset, int whence);
/* NOTE:  This is meant only for backwards compatibility.  */
/* You should probably be calling chkfileseek() directly.  */

size_t chkfileseek(FILE * stream, off_t offset, size_t size, int whence);
/* Preform a file seek with error checking.  */

long long chkfilelen(FILE *file, size_t size);
/* Return the length of a file (in blocks of 'size').  */

int read_int(FILE *infile, int byteswap);
/* Reads a binary integer value from the file 'infile' */

float read_float(FILE *infile, int byteswap);
/* Reads a binary float value from the file 'infile' */

double read_double(FILE *infile, int byteswap);
/* Reads a double precision value from the file 'infile' */

