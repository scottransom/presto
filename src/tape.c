#include <stdlib.h>
#include <stdio.h>

#define BLOCKLEN 49792

int main(int argc, char *argv[])
{
  FILE *rawfile;
  unsigned char block[BLOCKLEN];
  int ii, numwrite;

  if (argc != 4){
    printf("\nusage:  tape blockfile tapefile number\n\n");
    printf("    'blockfile' is the file where we will steal the 1st block to tape\n");
    printf("    'tapefile' is the file where we will append the block\n");
    printf("    'number' is the number of times we will tape the block\n\n");
    exit(0);
  }

  numwrite = atoi(argv[3]);
  printf("\nReading block from '%s' and appending it %d times to '%s'...\n\n", 
	 argv[1], numwrite, argv[2]);

  rawfile = fopen(argv[1], "rb");
  fread(block, sizeof(unsigned char), BLOCKLEN, rawfile);
  fclose(rawfile);

  rawfile = fopen(argv[2], "ab");
  for (ii=0; ii<numwrite; ii++)
    fwrite(block, sizeof(unsigned char), BLOCKLEN, rawfile);
  
  printf("Done.\n\n");
  return 0;
}
