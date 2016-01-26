#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <math.h>

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

int main(void)
{
    FILE *file;
    double tott;
    struct tms runtimes;
    int clk_tck = 1;

    tott = times(&runtimes);
    sleep(1);
    tott = times(&runtimes) - tott;
    clk_tck = (int) (tott / 10.0) * 10;
    printf("\nSlept for 1 sec.\n");
    printf("Measured time in clock_t was: %f\n\n", tott);
    printf("Therefore, CLK_TCK is probably = %d\n\n", clk_tck);
    file = fopen("../include/clk_tck.h", "w");
    fprintf(file, "#ifndef CLK_TCK\n");
    fprintf(file, "#define CLK_TCK %d\n", clk_tck);
    fprintf(file, "#endif\n");
    fclose(file);
    exit(0);
}
