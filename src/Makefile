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
PGPLOTLINK = $(PGPLOTLIBDIR) -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11 -lpthread -ldl

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
FFTLINK = $(FFTLIBDIR) -lfftw3f

# Other include directoried (i.e. glib.h)
OTHERINCDIR = -I/usr/include/glib-1.2

# Use the following if we are using DMALLOC
#DMALLOCLINK = -ldmalloc

# The standard PRESTO libraries to link into executables
PRESTOLINK = -L$(PRESTO)/lib -lpresto $(FFTLINK) $(DMALLOCLINK)

# When modifying the CLIG files, the is the location of the clig binary
CLIG = /usr/local/clig/bin/clig

# For Pentium or Athlon with GCC
# For Pentium3 or Pentium4s, use -march=pentiumpro
# For Athlons use -march=k6 or -march=k7
CC = gcc
FC = g77
CFLAGS = -I$(PRESTO)/include $(OTHERINCDIR) $(PGPLOTINCDIR) $(FFTINCDIR) \
	-DDATABASE='$(DATABASE)' -D$(FFTTYPE) -DUSEMMAP\
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DMACHINE_IS_LITTLE_ENDIAN=1\
	-g -O2 -march=athlon -Wall -W
#	-g -pg -O2 -march=athlon -Wall -W
#  
CLINKFLAGS = $(CFLAGS)
FFLAGS = -g -O2 -march=athlon
FLINKFLAGS = $(FFLAGS)
LINKCOMMAND1 = $(CC) -shared -o
LINKCOMMAND2 = ar rcs
#  The following are other possibilities and/or options for the CCFLAGS
#	-DUSEDMALLOC \

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
	cd ../clig ; $(CLIG) -o $*_cmd -d $<
	mv ../clig/$*_cmd.h ../include/
	mv ../clig/$*_cmd.c .
	cp ../clig/$*.1 ../docs/

PRESTOOBJS = amoeba.o barycenter.o birdzap.o cand_output.o\
	characteristics.o cldj.o chkio.o corr_prep.o corr_routines.o\
	correlations.o database.o dcdflib.o dispersion.o\
	fastffts.o fftcalls.o fminbr.o fold.o fresnl.o ioinf.o\
	get_candidates.o iomak.o ipmpar.o maximize_r.o maximize_rz.o\
	median.o minifft.o misc_utils.o clipping.o\
	orbint.o output.o read_fft.o responses.o\
	rzinterp.o rzwinterp.o select.o sorter.o swapendian.o\
	transpose.o twopass.o twopass_real_fwd.o\
	twopass_real_inv.o vectors.o multifiles.o mask.o\
	fitsfile.o hget.o hput.o imio.o

INSTRUMENTOBJS = multibeam.o bpp.o wapp.o gmrt.o spigot.o

PLOT2DOBJS = powerplot.o xyline.o

all: libpresto text xwin

prep: meminit timetest bincat

makewisdom:
	$(CC) $(CLINKFLAGS) -o $@ makewisdom.c $(FFTLINK) $(DMALLOCLINK) -lm
	./makewisdom
	cp fftw_wisdom.txt $(PRESTO)/lib

meminit:
	$(CC) $(CLINKFLAGS) -o $@ meminit.c $(DMALLOCLINK)
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

libpresto: libpresto.a libpresto.so
	@cd ../lib ; cp -sf ../src/libpresto* .

libpresto.so: $(PRESTOOBJS)
	$(LINKCOMMAND1) $@ $(PRESTOOBJS)

libpresto.a: $(PRESTOOBJS)
	$(LINKCOMMAND2) $@ $(PRESTOOBJS)

text: makedata makeinf mjd2cal realfft quicklook\
	search_bin search_rzw swap_endian prepdata\
	check_parkes_raw bary shiftdata dftfold\
	patchdata readfile toas2dat taperaw\
	accelsearch prepsubband cal2mjd split_parkes_beams\
	dat2sdat sdat2dat downsample rednoise un_sc_td

xwin: psrorbit window plotbincand profile \
	prepfold show_pfd\
	rfifind zapbirds plotdat plotfft\
	explorefft exploredat

links:
	cd ../bin ; find ../src -type f -perm 755 | xargs -n 1 -iXXX cp -s XXX .
	rm -f ../bin/libpresto.so

geomview: togeomview

mpi: mpiprepsubband

mpiprepsubband_utils.o: mpiprepsubband_utils.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband_utils.c

mpiprepsubband.o: mpiprepsubband.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband.c

mpiprepsubband: mpiprepsubband_cmd.c mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS)
	mpicc $(CLINKFLAGS) -o $@ mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lm

accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o
	$(CC) $(CLINKFLAGS) -o $@ accelsearch_cmd.o accel_utils.o accelsearch.o $(PRESTOLINK) -L/home/ransom/lib -lglib -lm

bary: bary.o
	$(CC) $(CLINKFLAGS) -o $@ bary.o $(PRESTOLINK) -lm

bincand: bincand_cmd.c bincand_cmd.o bincand.o
	$(CC) $(CLINKFLAGS) -o $@ bincand.o bincand_cmd.o $(PRESTOLINK) -lm

dftfold: dftfold_cmd.c dftfold_cmd.o dftfold.o
	$(CC) $(CLINKFLAGS) -o $@ dftfold.o dftfold_cmd.o $(PRESTOLINK) -lm

shiftdata: shiftdata.o
	$(CC) $(CLINKFLAGS) -o $@ shiftdata.o $(DMALLOCLINK) -lm

patchdata: patchdata.o
	$(CC) $(CLINKFLAGS) -o $@ patchdata.o $(DMALLOCLINK)

dat2sdat: dat2sdat.o
	$(CC) $(CLINKFLAGS) -o $@ dat2sdat.o $(PRESTOLINK) -lm

sdat2dat: sdat2dat.o
	$(CC) $(CLINKFLAGS) -o $@ sdat2dat.o $(PRESTOLINK) -lm

check_parkes_raw: check_parkes_raw.o multibeam.o
	$(CC) $(CLINKFLAGS) -o $@ check_parkes_raw.o multibeam.o $(PRESTOLINK) -lm

downsample: downsample_cmd.c downsample.o downsample_cmd.o
	$(CC) $(CLINKFLAGS) -o $@ downsample.o downsample_cmd.o $(PRESTOLINK) -lm

split_parkes_beams: split_parkes_beams.o
	$(CC) $(CLINKFLAGS) -o $@ split_parkes_beams.o

test_multifiles: test_multifiles.o
	$(CC) $(CLINKFLAGS) -o $@ test_multifiles.o $(PRESTOLINK) -lm

rfifind: rfifind_cmd.c rfifind_cmd.o rfifind.o rfi_utils.o rfifind_plot.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@  $(INSTRUMENTOBJS) $(PLOT2DOBJS) rfifind.o rfi_utils.o rfifind_cmd.o rfifind_plot.o $(PRESTOLINK) $(PGPLOTLINK) -lm

#showmulti: showmulti.o multi_prof_plot.o $(PLOT2DOBJS)
#	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti.o multi_prof_plot.o $(PRESTOLINK) $(PGPLOTLINK) -lm

#showmulti_dm: showmulti_dm.o multi_prof_plot.o $(PLOT2DOBJS)
#	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) showmulti_dm.o multi_prof_plot.o $(PRESTOLINK) $(PGPLOTLINK) -lm

#searchmultidms: searchmultidms.o com.o randlib.o $(PLOT2DOBJS)
#	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) searchmultidms.o com.o randlib.o $(PRESTOLINK) $(PGPLOTLINK) -lm

test_fold: test_fold.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) test_fold.o $(PRESTOLINK) $(PGPLOTLINK) -lm

plotdat: plotdat.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) plotdat.o $(PRESTOLINK) $(PGPLOTLINK) -lm

plotfft: plotfft.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) plotfft.o $(PRESTOLINK) $(PGPLOTLINK) -lm

prepdata: prepdata_cmd.c prepdata_cmd.o prepdata.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $@ prepdata.o prepdata_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lm

prepsubband: prepsubband_cmd.c prepsubband_cmd.o prepsubband.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $@ prepsubband.o prepsubband_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lm

prepfold: prepfold_cmd.c prepfold_cmd.o prepfold.o prepfold_utils.o prepfold_plot.o least_squares.o polycos.o readpar.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ prepfold.o prepfold_utils.o prepfold_plot.o prepfold_cmd.o least_squares.o polycos.o readpar.o $(PLOT2DOBJS) $(INSTRUMENTOBJS) $(LAPACKLINK) $(PRESTOLINK) $(PGPLOTLINK) -lm

show_pfd: show_pfd_cmd.c show_pfd.o show_pfd_cmd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ show_pfd.o show_pfd_cmd.o prepfold_utils.o prepfold_plot.o least_squares.o $(PLOT2DOBJS) $(LAPACKLINK) $(PRESTOLINK) $(PGPLOTLINK) -lm

testresp: ../tests/testresp.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) testresp.o $(PRESTOLINK) $(PGPLOTLINK) -lm

makedata: com.o randlib.o makedata.o
	$(CC) $(CLINKFLAGS) -o $@ com.o randlib.o makedata.o $(PRESTOLINK) -lm

makeinf: makeinf.o ioinf.o
	$(CC) $(CLINKFLAGS) -o $@ makeinf.o ioinf.o $(PRESTOLINK) -lm

mjd2cal: djcl.o mjd2cal.o
	$(CC) $(CLINKFLAGS) -o $@ djcl.o mjd2cal.o $(DMALLOCLINK) -lm

cal2mjd: cldj.o cal2mjd.o
	$(CC) $(CLINKFLAGS) -o $@ cldj.o cal2mjd.o $(DMALLOCLINK) -lm

plotbincand: plotbincand.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) plotbincand.o $(PRESTOLINK) $(PGPLOTLINK) -lm

profile: profile_cmd.c profile_cmd.o profile.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) profile.o profile_cmd.o $(PRESTOLINK) $(PGPLOTLINK) -lm

simplefold: simplefold.o polycos.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) simplefold.o polycos.o  $(PRESTOLINK) $(PGPLOTLINK) -lm

psrorbit: psrorbit.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) psrorbit.o $(PRESTOLINK) $(PGPLOTLINK) -lm

testbinresp: testbinresp.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ testbinresp.o $(PLOT2DOBJS) $(PGPLOTLINK) $(PRESTOLINK) -lm

quicklook: quicklook.o
	$(CC) $(CLINKFLAGS) -o $@ quicklook.o $(PRESTOLINK) -lm

readfile: readfile_cmd.c readfile_cmd.o readfile.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $@ readfile.o readfile_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lm

realfft: realfft_cmd.c realfft_cmd.o realfft.o
	$(CC) $(CLINKFLAGS) -o $@ realfft.o realfft_cmd.o $(PRESTOLINK) -lm

rednoise: rednoise_cmd.c rednoise.o rednoise_cmd.o
	$(CC) $(CLINKFLAGS) -o $@ rednoise.o rednoise_cmd.o $(PRESTOLINK) -lm

search_bin: search_bin_cmd.c search_bin_cmd.o search_bin.o
	$(CC) $(CLINKFLAGS) -o $@ search_bin.o search_bin_cmd.o $(PRESTOLINK) -lm

search_rzw: search_rzw_cmd.c search_rzw_cmd.o search_rzw.o
	$(CC) $(CLINKFLAGS) -o $@ search_rzw.o search_rzw_cmd.o $(PRESTOLINK) -lm

taperaw: taperaw.o
	$(CC) $(CLINKFLAGS) -o $@ taperaw.o $(DMALLOCLINK)

toas2dat: toas2dat_cmd.c toas2dat_cmd.o toas2dat.o
	$(CC) $(CLINKFLAGS) -o $@ toas2dat.o toas2dat_cmd.o $(DMALLOCLINK)

swap_endian:
	$(CC) $(CLINKFLAGS) -o $@ swap_endian.c $(DMALLOCLINK)

togeomview:
	$(CC) -O -o $@ togeomview.c $(DMALLOCLINK)

window: window.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ $(PLOT2DOBJS) window.o $(PRESTOLINK) $(PGPLOTLINK) -lm

zapbirds: zapbirds_cmd.c zapbirds_cmd.o zapbirds.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ zapbirds_cmd.o zapbirds.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -L/home/ransom/lib -lglib -lm

explorefft: explorefft.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ explorefft.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

exploredat: exploredat.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $@ exploredat.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

clean:
	rm -f *.o *~ *#

cleaner: clean
	find . -type f -perm 755 | xargs rm -f
	find . -type f -perm 775 | xargs rm -f
	rm -f libpresto.*

squeaky:  cleaner
	rm -f *.dat *.fft *.inf fftw_wisdom.txt
	rm -f core *.win* *.ps *_rzw *.tmp
	cd $(PRESTO)/lib ; rm -f bincat.dat
	cd $(PRESTO)/clig ; rm -f *# *~
	cd $(PRESTO)/docs ; rm -f *# *~
	cd $(PRESTO)/python ; rm -f *# *~ *.o *.pyc *.pyo
	cd $(PRESTO)/include ; rm -f *# *~
