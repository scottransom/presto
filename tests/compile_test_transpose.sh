gcc -g -O3 -Wall -W -ffast-math -I../include/ -o test_transpose test_transpose.c ../src/transpose.o ../src/vectors.o -lfftw3f -lm
