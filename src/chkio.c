#include "chkio.h"

FILE *chkfopen(char *path, const char *mode)
{
  FILE *file;

#ifdef USE_PIOFS
  if ((file = fopen64(path, mode)) == NULL) {
    perror("\nError in chkfopen()");
    printf("\n");
    exit(-1);
  }
#else
  if ((file = fopen(path, mode)) == NULL) {
    perror("\nError in chkfopen()");
    printf("\n");
    exit(-1);
  }
#endif
  return (file);
}


int chkfread(void *data, size_t type, size_t number, FILE * stream)
{
  unsigned int num;

  num = fread(data, type, number, stream);
  if (num != number && ferror(stream)) {
    perror("\nError in chkfread()");
    printf("\n");
    exit(-1);
  }
  return num;
}


int chkfwrite(void *data, size_t type, size_t number, FILE * stream)
{
  unsigned int num;
 
  num = fwrite(data, type, number, stream);
  if (num != number && ferror(stream)) {
    perror("\nError in chkfwrite()");
    printf("\n");
    exit(-1);
  }
  return num;
}


int chkfseek(FILE * stream, long offset, int whence)
/* NOTE:  This is meant only for backwards compatibility.  */
/* You should probably be calling chkfileseek() directly.  */
{
  return chkfileseek(stream, offset, 1, whence);
}


int chkfileseek(FILE * stream, long offset, size_t size, int whence)
{
  int rt;

#ifdef USE_PIOFS
  {
    off64_t pos;

    pos = (long long) offset * (long long) size;
    if ((rt = fseeko64(stream, pos, whence)) == -1){
      perror("\nError in chkfileseek()");
      printf("\n");
      exit(-1);
    }
  }
#else
  if ((rt = fseek(stream, offset * size, whence)) == -1) {
    perror("\nError in chkfileseek()");
    printf("\n");
    exit(-1);
  }
#endif
  return (rt);
}


unsigned long chkfilelen(FILE *file, size_t size)
{
  int filenum, rt;
#ifdef USE_PIOFS
  piofs_fstat_t buf;

  filenum = fileno(file);
  rt = piofsioctl(filenum, PIOFS_FSTAT, &buf);
  if (rt == -1){
    perror("\nError in chkfilelen()");
    printf("\n");
    exit(-1);
  }
#else
  struct stat buf;

  filenum = fileno(file);
  rt = fstat(filenum, &buf);
  if (rt == -1){
    perror("\nError in chkfilelen()");
    printf("\n");
    exit(-1);
  }
#endif
  return (unsigned long) (buf.st_size / size);
}



