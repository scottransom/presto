#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct FCOMPLEX {
  float r, i;
} fcomplex;

#define TOLERANCE (1e-2)
#define SQR(x) ((x) * (x))

float compute_error_complex(fcomplex A, fcomplex B)
{
  /* compute the relative error */
  float error = 0.0, a, mag;
  
  a = sqrt(SQR(A.r - B.r) + SQR(A.i - B.i));
  mag = 0.5 * (sqrt(SQR(A.r) + SQR(A.i)) +
	       sqrt(SQR(B.r) + SQR(B.i))) + TOLERANCE;
  a /= mag;
  if (a > error)
    error = a;
  return error;
}

int main(int argc, char *argv[])
{
  FILE *a, *b;
  fcomplex adat[1000], bdat[1000];
  float perc, err;
  int ard, brd, i, j=0;

  if (argc < 4){
    printf("\nUsage:  compare file1 file2 error(%%)\n\n");
    exit(1);
  }

  a = fopen(argv[1], "rb");
  b = fopen(argv[2], "rb");
  perc = atof(argv[3])/100.0;
  while ((ard = fread(adat, sizeof(fcomplex), 1000, a)) &&
	 (brd = fread(bdat, sizeof(fcomplex), 1000, b))){
    for(i = 0 ; i<ard ; i++){
      if((err=compute_error_complex(adat[i], bdat[i])) > perc)
	printf("Problem at %d:  %s=%g, %g  %s=%g, %g  err=%g\n",\
	       1000*j+i,argv[1],adat[i].r,adat[i].i
	       ,argv[2],bdat[i].r,bdat[i].i,err);
    }
    j++;
  }
  fclose(a);
  fclose(b);
  return(0);
}
  
