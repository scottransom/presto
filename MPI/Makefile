#  Makefile for MPI_PRESTO
#    for Pentium/Linux
#    by Scott M. Ransom
#     V 0.9  1 Dec 97

# Path to pulsar database files
DATABASE = "$(PRESTO)/lib"

# Path to MPE headers

MPEINC = /incoming/mpich/mpe

# The names of the libraries to produce
FFTAPPS  = libfftapps.a
FASTFFTS = libfastffts.a

# For Pentium with GCC
CC = mpicc
F77 = g77
CFLAGS = -I$(PRESTO)/include -I$(MPEINC) -DDATABASE='$(DATABASE)' \
	-Wall -ansi -pedantic \
	-O6 -ffast-math -m486 -malign-loops=2 -malign-jumps=2 \
	-malign-functions=2 -fno-strength-reduce

# For DEC Alpha
#CFLAGS =  -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-std1 -O3 -D_FASTMATH

# For Solaris
#CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-fast

# For SGI 64bit
#CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-64 -O2

# For SGI 32bit
#CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-32 -O2

# For Cray T3E
#CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-O3 -h aggress -X m

# For Cray T90
#CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) -DDATABASE='$(DATABASE)'\
#	-O3

MPIRZWFILES = mpisearch_rzw.o mpisearch_rzw_m.o mpisearch_rzw_s.o 


mpisearch_rzw: $(MPIRZWFILES)
	$(CC) $(CFLAGS) -o $@ $(MPIFILES) -L$(PRESTO)/lib -lfftapps

clean:
	rm -f *.o *~ *#

squeaky:
	rm -f core *.dat *.fft *.inf *.o *~ *.win* *.ps *_rzw *# 
	rm -f mpisearch_rzw












