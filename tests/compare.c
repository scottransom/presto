#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  FILE *a, *b;
  double perc;
  float adat[1000], bdat[1000];
  int ard, brd, i, j=0;

  if (argc < 4){
    printf("\nUsage:  compare file1 file2 error\n\n");
    exit(1);
  }

  a = fopen(argv[1], "rb");
  b = fopen(argv[2], "rb");
  perc = atof(argv[3])/100.0;
  while ((ard = fread(adat, sizeof(float), 1000, a)) == 1000 &&
	 (brd = fread(bdat, sizeof(float), 1000, b)) == 1000){
    for(i = 0 ; i<1000 ; i++){
      if(adat[i]==0.0){
	if(bdat[i]==0){
	  continue;
	} else {
	printf("Problem at %7d:  %s=%17.12e  %s=%17.12e\n",\
	       1000*j+i,argv[1],adat[i],argv[2],bdat[i]);
	}
      }	  
      if(fabs((adat[i]-bdat[i])/adat[i]) > perc)
	printf("Problem at %7d:  %s=%17.12e  %s=%17.12e\n",\
	       1000*j+i,argv[1],adat[i],argv[2],bdat[i]);
    }
    j++;
  }
  fclose(a);
  fclose(b);
  return(0);
}
  
