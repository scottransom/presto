#include <time.h>
#include <sys/times.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "clk_tck.h"
#include "vectors.h"
#include "fftw3.h"
#include "assert.h"

extern short transpose_float(float *a, int nx, int ny, unsigned char *move, 
                             int move_size);

fftwf_plan plan_transpose(int rows, int cols, float *in, float *out) {
    const unsigned flags = FFTW_MEASURE; /* other flags are possible */
    fftwf_iodim howmany_dims[2];
    
    howmany_dims[0].n  = rows;
    howmany_dims[0].is = cols;
    howmany_dims[0].os = 1;
    howmany_dims[1].n  = cols;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = rows;
    return fftwf_plan_guru_r2r(/*rank=*/ 0, /*dims=*/ NULL,
                               /*howmany_rank=*/ 2, howmany_dims,
                               in, out, /*kind=*/ NULL, flags);
}

void print_array(float *arr, int N, int M) {
    int ii, jj;
    for (ii=0; ii<N; ii++){
        for (jj=0; jj<M; jj++){
            printf("%6.4f ", arr[ii*M+jj]);
        }
        printf("\n");
    }
    printf("\n");
}


int main(int argc, char *argv[]) {
    float *array1, *array2;
    fftwf_plan tplan1, tplan2;
    struct tms runtimes;
    double ttim, stim, utim, tott;
    int ii, N, M, numtimes, move_size;
    unsigned char *tmpspace;

    if (argc <= 1 || argc > 4) {
        printf("\nUsage:  test_transpose N M #times (array[N][M])\n\n");
        exit(0);
    } else {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
        numtimes = atoi(argv[3]);
    }
    
    // Setup arrays
    array1 = gen_fvect(N * M);
    array2 = gen_fvect(N * M);
    for (ii = 0; ii < N * M; ii++)
        array1[ii] = array2[ii] = (float) rand() / (float) RAND_MAX;

    // Setup things for TOMs transpose
    move_size = N * M / 2;
    tmpspace = gen_bvect(move_size);

    // Start the timing for TOMs transpose
    tott = times(&runtimes) / (double) CLK_TCK;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    
    if (N * M <= 20)
        print_array(array1, N, M);
    for (ii = 0; ii < numtimes; ii++) {
        if (ii % 2)
            transpose_float(array1, M, N, tmpspace, move_size);
        else
            transpose_float(array1, N, M, tmpspace, move_size);
    }
    if (N * M <= 20) {
        if (numtimes % 2)
            print_array(array1, M, N);
        else
            print_array(array1, N, M);
    }

    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK - utim;
    stim = runtimes.tms_stime / (double) CLK_TCK - stim;
    ttim = utim + stim;

    // Check for correctness
    if (numtimes % 2 == 0) {
        for (ii = 0; ii < N * M; ii++)
            assert(fabs(array1[ii] - array2[ii]) < 1e-6);
    }

    printf("Timing summary (TOMS) NxM = %dx%d:\n", N, M);
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);
    
    // Now do the FFTW transpose

    // The plan messes up the input array!
    tplan1 = plan_transpose(N, M, array1, array1);
    tplan2 = plan_transpose(M, N, array1, array1);
    memcpy(array1, array2, sizeof(float) * N * M);

    tott = times(&runtimes) / (double) CLK_TCK;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;

    if (N * M <= 20)
        print_array(array1, N, M);
    for (ii = 0; ii < numtimes; ii++) {
        if (ii % 2)
            fftwf_execute_r2r(tplan2, array1, array1);
        else
            fftwf_execute_r2r(tplan1, array1, array1);
    }
    if (N * M <= 20) {
        if (ii % 2)
            print_array(array1, M, N);
        else
            print_array(array1, N, M);
    }

    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK - utim;
    stim = runtimes.tms_stime / (double) CLK_TCK - stim;
    ttim = utim + stim;

    // Check for correctness
    if (numtimes % 2 == 0) {
        for (ii = 0; ii < N * M; ii++)
            assert(fabs(array1[ii] - array2[ii]) < 1e-6);
    }

    printf("Timing summary (FFTW) NxM = %dx%d:\n", N, M);
    printf("CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n", \
	   ttim, utim, stim);
    printf("Total time elapsed:  %.3f sec.\n\n", tott);

    vect_free(array1);
    vect_free(array2);
    vect_free(tmpspace);
    return 0;
}

