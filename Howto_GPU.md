#GPUed PRESTO

The program of accelsearch in PRESTO was accelerated with CUDA by Jintao Luo. Any questions on how to compile and use the GPU routine may go to jluo@nrao.edu

##How to compile

###compile for CPU-only
By default PRESTO will be compiled for CPU-only. This is controlled by line 22 in $PRESTO/src/Makefile. For CPU-only, this line should be as:
>use_cuda = no# yes/no

###compile with GPU routines
Open the $PRESTO/src/Makefile

1.	Make sure line 22 is:
>use_cuda = yes# yes/no

2.	Check the CUDA-related variables and flags, and modify them if necessary.

##How to use
To run accelsearch on GPU, use the -cuda option. Or the program will run on CPU. For example: 
>accelsearch -numharm 16 -zmax 256 ur_data.dat -cuda 0

0 means using the 1st GPU in your machine.