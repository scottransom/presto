#  Makefile for PRESTO
#   for Pentium/Linux
#   by Scott M. Ransom

# Path to pulsar database files
DATABASE = "$(PRESTO)/lib"

# Path to PGPLOT includes
PGPLOTINCDIR = -I/usr/local/include
# Path to PGPLOT libraries
PGPLOTLIBDIR = -L$(PGPLOT_DIR)
# How to link with the PGPLOT libs
PGPLOTLINK = $(PGPLOTLIBDIR) -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11

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
RFFTLINK = -lsrfftw

# Other include directoried (i.e. glib.h)
OTHERINCDIR = -I/usr/include/glib-1.2

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
CC = gcc-3.0
FC = g77-3.0
CFLAGS = -I$(PRESTO)/include $(OTHERINCDIR) $(PGPLOTINCDIR) $(FFTINCDIR) \
	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
	-g -Wall -W -O2 -mcpu=pentiumpro -march=pentiumpro -fPIC
#	-DUSEDMALLOC \
#	-g -Wall -W -O -fPIC \
#	-g -ansi -pedantic -Wall -W -fpic \
#	-Dinline= -Wmissing-prototypes \
#	-Wstrict-prototypes -Wtraditional -Wconversion -Wshadow \
#	-Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings \
#	-Waggregate-return -fno-common -Wnested-externs
#	-O4 -mcpu=pentiumpro -march=pentiumpro -funroll-loops -fforce-addr -fno-strict-aliasing 
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
#CC = gcc
#FC = g77
#CFLAGS = -I$(PRESTO)/include $(PGPLOTINCDIR) $(FFTINCDIR) \
#	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) \
#	-g -O2 -fPIC -Wall -W 
#	-g -O -v -KPIC
#	-fast
#LIBLINKFLAGS = -shared
#LINKCOMMAND = $(CC) $(LIBLINKFLAGS) -o
#FFLAGS = -g -O -fPIC
#FLINK = f77
#FLINKFLAGS = 
#LIBLINKFLAGS = -G
#LINKCOMMAND = ld $(LIBLINKFLAGS) -o
#FFLAGS = -O
#FLINK = $(FC)
#FLINKFLAGS = $(FFLAGS)

# Rules for CLIG generated files

%_cmd.c : ../clig/%_cmd.cli
	cd ../clig ; clig -o $*_cmd -d $<
	mv ../clig/$*_cmd.h ../include/
	mv ../clig/$*_cmd.c .
	cp ../clig/$*.1 ../docs/

PRESTOOBJS = amoeba.o barycenter.o birdzap.o cand_output.o\
	characteristics.o cldj.o chkio.o corr_prep.o corr_routines.o\
	correlations.o database.o dcdflib.o dispersion.o\
	fastffts.o fftcalls.o fminbr.o fold.o fresnl.o ioinf.o\
	get_candidates.o iomak.o ipmpar.o maximize_r.o maximize_rz.o\
	median.o minifft.o misc_utils.o\
	orbint.o output.o read_fft.o responses.o\
	rzinterp.o rzwinterp.o select.o sorter.o swapendian.o\
	transpose.o twopass.o twopass_real_fwd.o\
	twopass_real_inv.o vectors.o multifiles.o mask.o

INSTRUMENTOBJS = multibeam.o bpp.o wapp.o

PLOT2DOBJS = powerplot.o xyline.o

all: $(PRESTOLIBNM) text xwin

prep: meminit timetest bincat

text: makedata makeinf mjd2cal realfft quicklook\
	search_bin search_rzw swap_endian prepdata\
	check_parkes_raw bary shiftdata dftfold\
	patchdata readfile toas2dat bincand taperaw\
	accelsearch prepsubband cal2mjd

xwin: psrorbit window plotbincand profile showmulti\
	prepfold showmulti_dm searchmultidms show_pfd\
	rfifind zapbirds

geomview: togeomview

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

accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o 
	$(CC) $(CFLAGS) -o $@ accelsearch_cmd.o accel_utils.o accelsearch.o $(PRESTOLINK) -lglib -lm
	cp $@ ../bin

bary: bary.o 
	$(CC) $(CFLAGS) -o $@ bary.o $(PRESTOLINK) -lm
	cp $@ ../bin

bincand: bincand_cmd.c bincand_cmd.o bincand.o
	$(CC) $(CFLAGS) -o $@ bincand.o bincand_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

dftfold: dftfold_cmd.c dftfold_cmd.o dftfold.o 
	$(CC) $(CFLAGS) -o $@ dftfold.o dftfold_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

shiftdata: shiftdata.o 
	$(CC) $(CFLAGS) -o $@ shiftdata.o $(DMALLOCLINK) -lm
	cp $@ ../bin

patchdata: patchdata.o 
	$(CC) $(CFLAGS) -o $@ patchdata.o $(DMALLOCLINK)
	cp $@ ../bin

check_parkes_raw: check_parkes_raw.o multibeam.o
	$(CC) $(CFLAGS) -o $@ check_parkes_raw.o multibeam.o $(PRESTOLINK) -lm
	cp $@ ../bin

test_multifiles: test_multifiles.o
	$(CC) $(CFLAGS) -o $@ test_multifiles.o $(PRESTOLINK) -lm

rfifind: rfifind_cmd.c rfifind_cmd.o rfifind.o rfi_utils.o rfifind_plot.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@  $(INSTRUMENTOBJS) $(PLOT2DOBJS) rfifind.o rfi_utils.o rfifind_cmd.o rfifind_plot.o $(PRESTOLINK) $(RFFTLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

showmulti: showmulti.o multi_prof_plot.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti.o multi_prof_plot.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

showmulti_dm: showmulti_dm.o multi_prof_plot.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti_dm.o multi_prof_plot.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

searchmultidms: searchmultidms.o com.o randlib.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) searchmultidms.o com.o randlib.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

test_fold: test_fold.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) test_fold.o $(PRESTOLINK) $(PGPLOTLINK) -lm

prepdata: prepdata_cmd.c prepdata_cmd.o prepdata.o $(INSTRUMENTOBJS)
	$(CC) $(CFLAGS) -o $@ prepdata.o prepdata_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) $(RFFTLINK) -lm
	cp $@ ../bin

prepsubband: prepsubband_cmd.c prepsubband_cmd.o prepsubband.o $(INSTRUMENTOBJS)
	$(CC) $(CFLAGS) -o $@ prepsubband.o prepsubband_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) $(RFFTLINK) -lm
	cp $@ ../bin

prepfold: prepfold_cmd.c prepfold_cmd.o prepfold.o prepfold_utils.o prepfold_plot.o least_squares.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ prepfold.o prepfold_utils.o prepfold_plot.o prepfold_cmd.o least_squares.o $(PLOT2DOBJS) $(INSTRUMENTOBJS) $(LAPACKLINK) $(PRESTOLINK) $(RFFTLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

show_pfd: show_pfd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ show_pfd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS) $(LAPACKLINK) $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

testresp: ../tests/testresp.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) testresp.o $(PRESTOLINK) $(PGPLOTLINK) -lm

makedata: com.o randlib.o makedata.o
	$(CC) $(CFLAGS) -o $@ com.o randlib.o makedata.o $(PRESTOLINK) -lm
	cp $@ ../bin

makeinf: makeinf.o ioinf.o
	$(CC) $(CFLAGS) -o $@ makeinf.o ioinf.o $(PRESTOLINK) -lm
	cp $@ ../bin

makewisdom:
	$(CC) $(CFLAGS) -o $@ makewisdom.c $(FFTLINK) $(DMALLOCLINK) -lm
	./makewisdom
	cp fftw_wisdom.txt $(PRESTO)/lib

mjd2cal: djcl.o mjd2cal.o
	$(CC) $(CFLAGS) -o $@ djcl.o mjd2cal.o $(DMALLOCLINK) -lm
	cp $@ ../bin

cal2mjd: cldj.o cal2mjd.o
	$(CC) $(CFLAGS) -o $@ cldj.o cal2mjd.o $(DMALLOCLINK) -lm
	cp $@ ../bin

plotbincand: plotbincand.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) plotbincand.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

profile: profile_cmd.c profile_cmd.o profile.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) profile.o profile_cmd.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

psrorbit: psrorbit.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) psrorbit.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

testbinresp: testbinresp.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ testbinresp.o $(PLOT2DOBJS) $(PGPLOTLINK) $(PRESTOLINK) -lm
	cp $@ ../bin

quicklook: quicklook.o
	$(CC) $(CFLAGS) -o $@ quicklook.o $(PRESTOLINK) -lm
	cp $@ ../bin

readfile: readfile_cmd.c readfile_cmd.o readfile.o $(INSTRUMENTOBJS)
	$(CC) $(CFLAGS) -o $@ readfile.o readfile_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) $(RFFTLINK) -lm
	cp $@ ../bin

realfft: realfft_cmd.c realfft_cmd.o realfft.o 
	$(CC) $(CFLAGS) -o $@ realfft.o realfft_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

#realfft_4step: realfft_4step.o
#	$(CC) $(CFLAGS) -o $@ realfft_4step.o $(PRESTOLINK) -lm
#	cp $@ ../bin

search_bin: search_bin_cmd.c search_bin_cmd.o search_bin.o 
	$(CC) $(CFLAGS) -o $@ search_bin.o search_bin_cmd.o $(PRESTOLINK)  -lm
	cp $@ ../bin

search_rzw: search_rzw_cmd.c search_rzw_cmd.o search_rzw.o 
	$(CC) $(CFLAGS) -o $@ search_rzw.o search_rzw_cmd.o $(PRESTOLINK) -lm
	cp $@ ../bin

taperaw: taperaw.o
	$(CC) $(CFLAGS) -o $@ taperaw.o $(DMALLOCLINK) 
	cp $@ ../bin

toas2dat: toas2dat_cmd.c toas2dat_cmd.o toas2dat.o 
	$(CC) $(CFLAGS) -o $@ toas2dat.o toas2dat_cmd.o $(DMALLOCLINK) 
	cp $@ ../bin

swap_endian:
	$(CC) $(CFLAGS) -o $@ swap_endian.c $(DMALLOCLINK) 
	cp $@ ../bin

togeomview:
	$(CC) -O -o $@ togeomview.c $(DMALLOCLINK)
	cp $@ ../bin

window: window.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) window.o $(PRESTOLINK) $(PGPLOTLINK) -lm
	cp $@ ../bin

zapbirds: zapbirds_cmd.c zapbirds_cmd.o zapbirds.o $(PLOT2DOBJS)
	$(FLINK) $(FLINKFLAGS) -o $@ zapbirds_cmd.o zapbirds.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lglib -lm
	cp $@ ../bin

clean:
	rm -f *.o *~ *#

cleaner:  clean
	rm -f makedata makeinf mjd2cal realfft timetest cal2mjd
	rm -f quicklook search_bin search_rzw swap_endian        
	rm -f prepdata prepfold check_parkes_raw bary toas2dat
	rm -f shiftdata patchdata readfile psrorbit window       
	rm -f plotbincand profile showmulti showmulti_dm show_pfd         
	rm -f searchmultidms togeomview meminit bincat makewisdom
	rm -f taperaw bincand dftfold rfifind zapbirds test_multifiles 
	rm -f accelsearch prepsubband $(PRESTOLIBNM)
	cd $(PRESTO)/lib ; rm -f $(PRESTOLIBNM)

squeaky:  cleaner
	rm -f *.dat *.fft *.inf fftw_wisdom.txt
	rm -f core *.win* *.ps *_rzw *.tmp
	cd $(PRESTO)/clig ; rm -f *# *~
	cd $(PRESTO)/docs ; rm -f *# *~
	cd $(PRESTO)/python ; rm -f *# *~ *.o *.pyc *.pyo
	cd $(PRESTO)/include ; rm -f *# *~
	cd $(PRESTO)/lib ; rm -f bincat.dat
	cd $(PRESTO)/bin ; rm -f makedata makeinf mjd2cal realfft cal2mjd
	cd $(PRESTO)/bin ; rm -f quicklook search_bin search_rzw swap_endian
	cd $(PRESTO)/bin ; rm -f prepdata prepfold check_parkes_raw bary
	cd $(PRESTO)/bin ; rm -f shiftdata patchdata readfile psrorbit
	cd $(PRESTO)/bin ; rm -f plotbincand profile showmulti showmulti_dm
	cd $(PRESTO)/bin ; rm -f searchmultidms togeomview window toas2dat
	cd $(PRESTO)/bin ; rm -f show_pfd taperaw bincand dftfold rfifind
	cd $(PRESTO)/bin ; rm -f zapbirds test_multifiles accelsearch
	cd $(PRESTO)/bin ; rm -f prepsubband
