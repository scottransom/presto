#include "chkio.h"

#ifndef __USE_FILE_OFFSET64
# ifndef __USE_LARGEFILE
#  define fseeko  fseek
# endif
#endif

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tmpswap=(a);(a)=(b);(b)=tmpswap;
#endif

static unsigned char tmpswap;

FILE *chkfopen(char *path, const char *mode)
{
   FILE *file;

#ifdef USE_PIOFS
   if ((file = fopen64(path, mode)) == NULL) {
      perror("\nError in chkfopen()");
      printf("   path = '%s'\n", path);
      exit(-1);
   }
#else
   if ((file = fopen(path, mode)) == NULL) {
      perror("\nError in chkfopen()");
      printf("   path = '%s'\n", path);
      exit(-1);
   }
#endif
   return (file);
}


size_t chkfread(void *data, size_t type, size_t number, FILE * stream)
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


size_t chkfwrite(void *data, size_t type, size_t number, FILE * stream)
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


size_t chkfseek(FILE * stream, long offset, int whence)
/* NOTE:  This is meant only for backwards compatibility.  */
/* You should probably be calling chkfileseek() directly.  */
{
   return chkfileseek(stream, offset, 1, whence);
}


size_t chkfileseek(FILE * stream, long offset, size_t size, int whence)
{
   int rt;

#ifdef USE_PIOFS
   {
      off64_t pos;

      pos = (long long) offset *(long long) size;
      if ((rt = fseeko64(stream, pos, whence)) == -1) {
         perror("\nError in chkfileseek()");
         printf("\n");
         exit(-1);
      }
   }
#else
   if ((rt = fseeko(stream, offset * size, whence)) == -1) {
      perror("\nError in chkfileseek()");
      printf("\n");
      exit(-1);
   }
#endif
   return (rt);
}


long long chkfilelen(FILE * file, size_t size)
{
   int filenum, rt;
#ifdef USE_PIOFS
   piofs_fstat_t buf;

   filenum = fileno(file);
   rt = piofsioctl(filenum, PIOFS_FSTAT, &buf);
   if (rt == -1) {
      perror("\nError in chkfilelen()");
      printf("\n");
      exit(-1);
   }
#else
   struct stat buf;

   filenum = fileno(file);
   rt = fstat(filenum, &buf);
   if (rt == -1) {
      perror("\nError in chkfilelen()");
      printf("\n");
      exit(-1);
   }
#endif
   return (long long) (buf.st_size / size);
}

int read_int(FILE * infile, int byteswap)
/* Reads a binary integer value from the file 'infile' */
{
   int itmp;

   chkfread(&itmp, sizeof(int), 1, infile);
   if (byteswap) {
      unsigned char *buffer = (unsigned char *) (&itmp);
      SWAP(buffer[0], buffer[3]);
      SWAP(buffer[1], buffer[2]);
   }
   return itmp;
}

float read_float(FILE * infile, int byteswap)
/* Reads a binary float value from the file 'infile' */
{
   float ftmp;

   chkfread(&ftmp, sizeof(float), 1, infile);
   if (byteswap) {
      unsigned char *buffer = (unsigned char *) (&ftmp);
      SWAP(buffer[0], buffer[3]);
      SWAP(buffer[1], buffer[2]);
   }
   return ftmp;
}

double read_double(FILE * infile, int byteswap)
/* Reads a double precision value from the file 'infile' */
{
   double dtmp;

   chkfread(&dtmp, sizeof(double), 1, infile);
   if (byteswap) {
      unsigned char *buffer = (unsigned char *) (&dtmp);
      SWAP(buffer[0], buffer[7]);
      SWAP(buffer[1], buffer[6]);
      SWAP(buffer[2], buffer[5]);
      SWAP(buffer[3], buffer[4]);
   }
   return dtmp;
}
