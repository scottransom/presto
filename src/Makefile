#  Makefile for PRESTO
#   for Pentium/Linux
#   by Scott M. Ransom

# Path to pulsar database files
DATABASE = "$(PRESTO)/lib"

# Path to PGPLOT includes
PGPLOTINCDIR = -I$(PGPLOT_DIR)
# Path to PGPLOT libraries
PGPLOTLIBDIR = -L$(PGPLOT_DIR)
# How to link with the PGPLOT libs
PGPLOTLINK = $(PGPLOTLIBDIR) -lcpgplot -lpgplot -L/usr/X11/lib -lX11

# You must define one of the following FFTTYPES (USEFFTW is default).
# If you want to use FFTW for the FFTs, place sfftw.h and libsfftw.a
# (or the shared lib version) in a place where we can find it,
# and then define the FFTTYPE and FFTLIBS variables appropriately.
# Notes:  PRESTO uses the single precision version of FFTW.
#         You must also run make makewisdom after you make prep.
# Other options are USESGIFFT or REGFFT
FFTTYPE = USEFFTW
# Path to FFT includes
FFTINCDIR = -I/usr/local/include
# Path to FFT libraries
FFTLIBDIR = -L/usr/local/lib
# How to link with the FFTW libs
FFTLINK = $(FFTLIBDIR) -lsfftw

# The name of the library to produce
PRESTOLIB = presto

# Use the following to define the propoer library suffix for your
# computer (usually 'so' or 'sl' for shared libraries or 'a' for
# static libraries.
LIBSUFFIX = so

# The true library names
PRESTOLIBNM = lib$(PRESTOLIB).$(LIBSUFFIX)

# Use the following if we are using DMALLOC
#DMALLOCLINK = -ldmalloc

# The standard PRESTO libraries to link into executables
PRESTOLINK = -L$(PRESTO)/lib -l$(PRESTOLIB) $(FFTLINK) $(DMALLOCLINK)

# For Pentium with GCC
CC = gcc
FC = g77
CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
	-g -Wall -W -O2 -mpentiumpro -fpic \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
#	-DUSEDMALLOC
#	-g -ansi -pedantic -Wall -W -O2 -mpentiumpro -fpic \
#	-g -ansi -pedantic -Wall -W -fpic \
#	-Dinline= -Wmissing-prototypes \
#	-Wstrict-prototypes -Wtraditional -Wconversion -Wshadow \
#	-Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings \
#	-Waggregate-return -fno-common -Wnested-externs
#	-O4 -mpentiumpro -funroll-loops -fforce-addr -fno-strict-aliasing 
LIBLINKFLAGS = -shared
LINKCOMMAND = $(CC) $(LIBLINKFLAGS) -o
FFLAGS = -g -O -fpic
FLINK = $(FC)
FLINKFLAGS = $(FFLAGS)

# For DEC Alpha
#CC = cc
#FC = f77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-newc -w0 -O4 -fast -ansi_alias -ansi_args -tune host -std
#	-g -newc -w0 -ansi_alias -ansi_args -std
#LINKCOMMAND = ar rcs
#FFLAGS = -O
#FLINK = $(FC) -nofor_main
#FLINKFLAGS = $(CFLAGS)

# For Solaris
#CC = cc
#FC = f77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-g -O -Xc -v -KPIC
#	-fast
#LIBLINKFLAGS = -G
#LINKCOMMAND = ld $(LIBLINKFLAGS) -o
#FFLAGS = -O
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

# For IBM SP2
#CC = cc
#FC = f77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-DUSE_PIOFS -g -O2 -qarch=pwr2 -qtune=pwr2
#EXTRALINK = -bI:/usr/include/piofs/piofs.exp
#LINKCOMMAND = ar rcs
#FFLAGS = -O
#FLINK = xlf
#FLINKFLAGS = 

# For SGI 64bit
#CC = cc
#FC = f77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-64 -Ofast=ip27 -r10000
#LINKCOMMAND = ar rcs
#FFLAGS = -64 -Ofast=ip27 -r10000
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

# For SGI 32bit
#CC = cc
#FC = f77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-32 -O2
#LINKCOMMAND = ar rcs
#FFLAGS = -O
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

# For Cray T3E
#CC = cc
#FC = f90 -dp
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-O3 -h aggress -X m
#LINKCOMMAND = ar rcs
#FFLAGS = -O
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

# For Cray T90
#CC = cc
#FC = cf77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-O3
#LINKCOMMAND = ar rcs
#FFLAGS = -O
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

PRESTOOBJS = apprgrdn.o barycenter.o birdzap.o cand_output.o\
	characteristics.o chkio.o corr_prep.o corr_routines.o\
	correlations.o database.o dcdflib.o dispersion.o\
	fastffts.o fftcalls.o fminbr.o fold.o fresnl.o ioinf.o\
	get_candidates.o iomak.o ipmpar.o maximize_r.o maximize_rz.o\
	maximize_rzw.o minifft.o multibeam.o misc_utils.o orbint.o\
	output.o PMsurv_phasemod.o read_fft.o responses.o\
	rzinterp.o rzwinterp.o select.o solvopt.o sorter.o\
	transpose.o twopass.o twopass_real_fwd.o\
	twopass_real_inv.o vectors.o multifiles.o

PLOT2DOBJS = powerplot.o xyline.o

all: $(PRESTOLIBNM) text xwin

prep: meminit timetest bincat

text: makedata makeinf mjd2cal realfft quicklook\
	search_bin search_rzw swap_endian prepdata\
	check_parkes_raw bary shiftdata dftfold\
	patchdata readfile toas2dat bincand taperaw

xwin: psrorbit window plotbincand profile showmulti\
	prepfold showmulti_dm searchmultidms show_pfd

geomview: togeomview

cligs:
	rm -f *_cmd.*
	cd ../clig ; clig -o profile_cmd -d profile_cmd.cli
	cd ../clig ; clig -o prepdata_cmd -d prepdata_cmd.cli
	cd ../clig ; clig -o prepfold_cmd -d prepfold_cmd.cli
	cd ../clig ; clig -o bincand_cmd -d bincand_cmd.cli
	cd ../clig ; clig -o search_bin_cmd -d search_bin_cmd.cli
	cd ../clig ; clig -o search_rzw_cmd -d search_rzw_cmd.cli
	cd ../clig ; clig -o readfile_cmd -d readfile_cmd.cli
	cd ../clig ; clig -o toas2dat_cmd -d toas2dat_cmd.cli
	cd ../clig ; clig -o dftfold_cmd -d dftfold_cmd.cli
	cd ../clig ; clig -o realfft_cmd -d realfft_cmd.cli
	mv ../clig/*_cmd.h ../include/
	mv ../clig/*_cmd.c .
	cp ../clig/*.1 ../docs/

$(PRESTOLIBNM): $(PRESTOOBJS)
	$(LINKCOMMAND) $@ $(PRESTOOBJS)
	cp $@ $(PRESTO)/lib

meminit:
	$(CC) $(CFLAGS) -o $@ meminit.c $(DMALLOCLINK)
	./meminit
	rm -f meminit

bincat:
	$(FC) -o $@ bincat.f rdascat.f
	./bincat
	rm -f bincat

timetest: 
	$(CC) -o $@ timetest.c $(DMALLOCLINK)
	./timetest
	rm -f timetest

bary: bary.o 
	$(CC) $(CFLAGS) -o $@ bary.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

bincand: bincand.o bincand_cmd.o
	$(CC) $(CFLAGS) -o $@ bincand.o bincand_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

dftfold: dftfold.o dftfold_cmd.o
	$(CC) $(CFLAGS) -o $@ dftfold.o dftfold_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

shiftdata: shiftdata.o 
	$(CC) $(CFLAGS) -o $@ shiftdata.o $(DMALLOCLINK) -lm
	cp $@ ../bin

patchdata: patchdata.o 
	$(CC) $(CFLAGS) -o $@ patchdata.o $(DMALLOCLINK)
	cp $@ ../bin

check_parkes_raw: check_parkes_raw.o
	$(CC) $(CFLAGS) -o $@ check_parkes_raw.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

test_multifiles: test_multifiles.o
	$(CC) $(CFLAGS) -o $@ test_multifiles.o $(EXTRALINK) $(PRESTOLINK) -lm

showmulti: showmulti.o multi_prof_plot.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti.o multi_prof_plot.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

showmulti_dm: showmulti_dm.o multi_prof_plot.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti_dm.o multi_prof_plot.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

searchmultidms: searchmultidms.o com.o randlib.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) searchmultidms.o com.o randlib.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

test_fold: test_fold.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) test_fold.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm

prepdata: prepdata.o prepdata_cmd.o
	$(CC) $(CFLAGS) -o $@ prepdata.o prepdata_cmd.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

prepfold: prepfold.o prepfold_utils.o prepfold_plot.o prepfold_cmd.o least_squares.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ prepfold.o prepfold_utils.o prepfold_plot.o prepfold_cmd.o least_squares.o $(PLOT2DOBJS) $(LAPACKLINK) $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

show_pfd: show_pfd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ show_pfd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS) $(LAPACKLINK) $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

testresp: ../tests/testresp.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) testresp.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm

makedata: com.o randlib.o makedata.o
	$(CC) $(CFLAGS) -o $@ com.o randlib.o makedata.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

makeinf: makeinf.o ioinf.o
	$(CC) $(CFLAGS) -o $@ makeinf.o ioinf.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

makewisdom:
	$(CC) $(CFLAGS) -o $@ makewisdom.c $(FFTLINK) $(DMALLOCLINK) -lm
	./makewisdom
	cp fftw_wisdom.txt $(PRESTO)/lib

mjd2cal: djcl.o mjd2cal.o
	$(CC) $(CFLAGS) -o $@ djcl.o mjd2cal.o $(DMALLOCLINK) -lm
	cp $@ ../bin

plotbincand: plotbincand.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) plotbincand.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

profile: profile.o profile_cmd.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) profile.o profile_cmd.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

psrorbit: psrorbit.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) psrorbit.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

testbinresp: testbinresp.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ testbinresp.o $(PLOT2DOBJS) $(PGPLOTLINK) $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

quicklook: quicklook.o
	$(CC) $(CFLAGS) -o $@ quicklook.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

readfile:  readfile.o readfile_cmd.o
	$(CC) $(CFLAGS) -o $@ readfile.o readfile_cmd.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

realfft: realfft.o realfft_cmd.o
	$(CC) $(CFLAGS) -o $@ realfft.o realfft_cmd.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

#realfft_4step: realfft_4step.o
#	$(CC) $(CFLAGS) -o $@ realfft_4step.o $(EXTRALINK) $(PRESTOLINK) -lm
#	cp $@ ../bin

search_bin: search_bin.o search_bin_cmd.o
	$(CC) $(CFLAGS) -o $@ search_bin.o search_bin_cmd.o $(EXTRALINK) $(PRESTOLINK)  -lm
	cp $@ ../bin

search_rzw: search_rzw.o search_rzw_cmd.o
	$(CC) $(CFLAGS) -o $@ search_rzw.o search_rzw_cmd.o $(EXTRALINK) $(PRESTOLINK) -lm
	cp $@ ../bin

taperaw: taperaw.o
	$(CC) $(CFLAGS) -o $@ taperaw.o $(DMALLOCLINK) 
	cp $@ ../bin

toas2dat: toas2dat.o toas2dat_cmd.o
	$(CC) $(CFLAGS) -o $@ toas2dat.o toas2dat_cmd.o $(DMALLOCLINK) 
	cp $@ ../bin

swap_endian:
	$(CC) $(CFLAGS) -o $@ swap_endian.c $(DMALLOCLINK) 
	cp $@ ../bin

togeomview:
	$(CC) -O -o $@ togeomview.c $(DMALLOCLINK)
	cp $@ ../bin

window: window.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) window.o $(EXTRALINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

clean:
	rm -f *.o *~ *#

cleaner:  clean
	rm -f makedata makeinf mjd2cal realfft timetest
	rm -f quicklook search_bin search_rzw swap_endian        
	rm -f prepdata prepfold check_parkes_raw bary toas2dat
	rm -f shiftdata patchdata readfile psrorbit window       
	rm -f plotbincand profile showmulti showmulti_dm show_pfd         
	rm -f searchmultidms togeomview meminit bincat makewisdom
	rm -f taperaw bincand $(PRESTOLIBNM)
	cd $(PRESTO)/lib ; rm -f $(PRESTOLIBNM)

squeaky:  cleaner
	rm -f *.dat *.fft *.inf fftw_wisdom.txt
	rm -f core *.win* *.ps *_rzw *.tmp
	cd $(PRESTO)/include ; rm -f *# *~
	cd $(PRESTO)/lib ; rm -f bincat.dat
	cd $(PRESTO)/bin ; rm -f makedata makeinf mjd2cal realfft
	cd $(PRESTO)/bin ; rm -f quicklook search_bin search_rzw swap_endian
	cd $(PRESTO)/bin ; rm -f prepdata prepfold check_parkes_raw bary
	cd $(PRESTO)/bin ; rm -f shiftdata patchdata readfile psrorbit
	cd $(PRESTO)/bin ; rm -f plotbincand profile showmulti showmulti_dm
	cd $(PRESTO)/bin ; rm -f searchmultidms togeomview window toas2dat
	cd $(PRESTO)/bin ; rm -f show_pfd taperaw bincand
