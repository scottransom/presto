#  Makefile for PRESTO Tests
#   for Pentium/Linux
#   by Scott M. Ransom
#    V 0.9  23 Oct 97

# Path to pulsar database files
DATABASE = "$(PRESTO)/lib"

# Path to PGPLOT includes
PGPLOTINC = "/usr/local/pgplot"
# Path to PGPLOT libraries
PGPLOTLIB = "/usr/local/pgplot"
# How to link with the PGPLOT libs
PGPLOTLINK = -L$(PGPLOTLIB) -lcpgplot -lpgplot -L/usr/X11/lib -lX11

# You must define one of the following FFTTYPES

# If you want to use FFTW for the FFTs, place fftw.h and libfftw.a
# in $(PRESTO)/LIB and then define the following:
# You must also run make makewisdom after you make prep.
FFTTYPE = USEFFTW
LINKPRESTO = -L$(PRESTO)/lib -lpresto -lsfftw
# If you want to use SGI's SGIMATH FFT, define the following
#FFTTYPE = USESGIFFT
#FFTLIBS = -L$(PRESTO)/lib -lpresto -lcomplib.sgimath
# The default FFT is a six-step FFT by Scott Ransom
#FFTTYPE = REGFFT
#FFTLIBS = -L$(PRESTO)/lib -lpresto
 
# The names of the libraries to produce
LIBPRESTO  = libpresto.a
FASTFFTS = libransomffts.a

# For Pentium with GCC
CC = egcc
F77 = g77
CFLAGS = -I$(PRESTO)/include -I$(PGPLOTINC) \
	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
	-Wall -ansi -pedantic -pipe \
	-O3 -ffast-math -mpentium -funroll-loops -fforce-addr

# For DEC Alpha
#CC = cc
#F77 = f77
#CFLAGS += -std1 -O3 -D_FASTMATH

# For Solaris
#CC = cc
#F77 = f77
#CFLAGS += -O3 -msupersparc

# For SGI 64bit
#CC = cc
#F77 = f77
#CFLAGS += -64 -O2

# For SGI 32bit
#CC = cc
#F77 = f77
#CFLAGS += -32 -O2

# For Cray T3E
#CC = cc
#F77 = f90 -dp
#CFLAGS += -O3 -h aggress -X m

# For Cray Y-MP
#CC = cc
#F77 = cf77
#CFLAGS += -O3

TEST_FFTS = test_ffts.o ../src/randlib.o ../src/com.o

all: test_ffts test_3d test_2d test_binsearch test_rzwsearch

prep:
	@../bin/makedata test_bin
	@../bin/realfft_sing -1 test_bin
	@../bin/makedata test_rzw
	@../bin/realfft_sing -1 test_rzw

test_ffts: $(TEST_FFTS)
	$(CC) $(CFLAGS) -o $@ $(TEST_FFTS) -lm $(LINKPRESTO)

tester: tester.o ../src/xyline.o ../src/powerplot.o
	$(F77) -o $@ ../src/xyline.o ../src/powerplot.o \
	$(FFTLIBS) $(PGPLOTLINK) -lm

test_3d: test_3d.o
	$(CC) $(CFLAGS) -o $@ test_3d.o -lm $(LINKPRESTO)

test_2d: test_2d.o
	$(CC) $(CFLAGS) -o $@ test_2d.o -lm $(LINKPRESTO) $(LINKPGPLOT)

test_binsearch: test_binsearch.o
	$(CC) $(CFLAGS) -o $@ test_binsearch.o -lm $(LINKPRESTO)

test_rzwsearch: test_rzwsearch.o
	$(CC) $(CFLAGS) -o $@ test_rzwsearch.o -lm $(LINKPRESTO)

clean:
	rm -f *.o *~ *# test_ffts test_3d test_2d
	rm -f test_binsearch test_rzwsearch














