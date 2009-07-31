#  Makefile for PRESTO
#   for Pentium/Linux
#   by Scott M. Ransom

# OS type
OS = Linux
#OS = OSX

# Linux is the first choice
ifeq ($(OS),Linux)
	LIBSUFFIX = .so
	LIBCMD = -shared
	SYSDIR = /usr
	LOCDIR = /usr/local
# else assume Darwin (i.e. OSX)
else
	LIBSUFFIX = .dylib
	LIBCMD = -dynamiclib
	SYSDIR = /sw
	LOCDIR = /sw
endif

# Path to PGPLOT includes
PGPLOTINCDIR = -I$(PGPLOT_DIR)
# Path to PGPLOT libraries
PGPLOTLIBDIR = -L$(PGPLOT_DIR)
# How to link with the PGPLOT libs
PGPLOTLINK = $(PGPLOTLIBDIR) -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11 -lpthread -lpng -lz -ldl

# You must define one of the following FFTTYPES (USEFFTW is default).
# If you want to use FFTW for the FFTs, place sfftw.h and libsfftw.a
# (or the shared lib version) in a place where we can find it,
# and then define the FFTTYPE and FFTLIBS variables appropriately.
# Notes:  PRESTO uses the single precision version of FFTW.
#         You must also run make makewisdom after you make prep.
# Other options are USESGIFFT or REGFFT
FFTTYPE = USEFFTW
# Path to FFT includes
FFTINCDIR = -I$(LOCDIR)/include
# Path to FFT libraries
FFTLIBDIR = -L$(LOCDIR)/lib
# How to link with the FFTW libs
FFTLINK = $(FFTLIBDIR) -lfftw3f

# Other include directories (i.e. glib.h)
OTHERINCDIR = -I$(SYSDIR)/include/glib-2.0 -I$(SYSDIR)/lib/glib-2.0/include
# Other link directories (i.e. libcfitsio)
OTHERLIBDIR = -L.

# The standard PRESTO libraries to link into executables
PRESTOLINK = -L$(PRESTO)/lib -lpresto $(FFTLINK)

# When modifying the CLIG files, the is the location of the clig binary
CLIG = /usr/bin/clig

# For Pentium IVs, use -march=pentium4
# For Opterons, use -march=opteron
# For Intel Core2 chips, use -march=nocona
# For gcc-4.X, use -mtune=native
CC = gcc
FC = gfortran
#FC = g77
CFLAGS = -I$(PRESTO)/include $(OTHERINCDIR) $(PGPLOTINCDIR) $(FFTINCDIR) \
	-D$(FFTTYPE) -DUSEMMAP -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
	-O3 -ffast-math -Wall -W -fPIC
#	-g -Wall -W
CLINKFLAGS = $(CFLAGS)
# NOTE:  Be careful of upping the optimization on the
#        FFLAGs.  Certain compilers (i.e. on Intel Macs) will
#        cause errors for the code in least_squares.f
FFLAGS = -g -fPIC
FLINKFLAGS = $(FFLAGS)
LINKCOMMAND1 = $(CC) $(LIBCMD) -o
LINKCOMMAND2 = ar rcs

# Add to the search path for the executables
VPATH = ../lib:../bin

# Rules for CLIG generated files
%_cmd.c : ../clig/%_cmd.cli
	cd ../clig ; $(CLIG) -o $*_cmd -d $<
	mv ../clig/$*_cmd.h ../include/
	mv ../clig/$*_cmd.c .
	cp ../clig/$*.1 ../docs/

PRESTOOBJS = amoeba.o atwood.o barycenter.o birdzap.o cand_output.o\
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

INSTRUMENTOBJS = multibeam.o bpp.o gmrt.o spigot.o sigproc_fb.o\
	wapp.o wapp_head_parse.o wapp_y.tab.o psrfits.o

PLOT2DOBJS = powerplot.o xyline.o

BINARIES = makedata makeinf mjd2cal realfft quicklook\
	search_bin search_rzw swap_endian prepdata\
	check_parkes_raw bary shiftdata dftfold\
	patchdata readfile toas2dat taperaw\
	accelsearch prepsubband cal2mjd split_parkes_beams\
	dat2sdat sdat2dat downsample rednoise un_sc_td bincand\
	dump_spigot_zerolag spigot2filterbank\
	spigotSband2filterbank GBT350filterbank\
	psrorbit window plotbincand prepfold show_pfd\
	rfifind zapbirds explorefft exploredat weight_psrfits

all: libpresto slalib binaries

# Default indentation is K&R style with no-tabs,
# an indentation level of 3, and a line-length of 85
indent:
	indent -kr -nut -i3 -l85 *.c
	rm *.c~

prep:
	touch *_cmd.c

makewisdom:
	$(CC) $(CLINKFLAGS) -o $@ makewisdom.c $(FFTLINK) $(DMALLOCLINK) -lm
	./makewisdom
	cp fftw_wisdom.txt $(PRESTO)/lib

timetest: 
	$(CC) -o $@ timetest.c $(DMALLOCLINK)
	./timetest
	rm -f timetest

libpresto: libpresto$(LIBSUFFIX)

libpresto$(LIBSUFFIX): $(PRESTOOBJS)
	$(LINKCOMMAND1) $(PRESTO)/lib/$@ $(PRESTOOBJS) $(FFTLINK)

slalib: libsla$(LIBSUFFIX)
	cd slalib ; $(FC) -o sla_test sla_test.f -fno-second-underscore -L$(PRESTO)/lib -lsla
	slalib/sla_test

libsla$(LIBSUFFIX):
	cd slalib ; $(FC) $(FFLAGS) -fno-second-underscore -c -I. *.f *.F
	rm slalib/sla_test.o
	cd slalib ; $(FC) $(LIBCMD) -o $(PRESTO)/lib/libsla$(LIBSUFFIX) -fno-second-underscore *.o

binaries: $(BINARIES)

mpi: mpiprepsubband

mpiprepsubband_utils.o: mpiprepsubband_utils.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband_utils.c

mpiprepsubband.o: mpiprepsubband.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband.c

mpiprepsubband: mpiprepsubband_cmd.c mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS)
	mpicc $(CLINKFLAGS) -o $(PRESTO)/bin/$@ mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lcfitsio -lm

accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o $(PRESTOLINK) $(OTHERLIB) -lglib-2.0 -lm

bary: bary.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ bary.o $(PRESTOLINK) -lm

bincand: bincand_cmd.c bincand_cmd.o bincand.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ bincand.o bincand_cmd.o $(PRESTOLINK) -lm

dftfold: dftfold_cmd.c dftfold_cmd.o dftfold.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ dftfold.o dftfold_cmd.o $(PRESTOLINK) -lm

shiftdata: shiftdata.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ shiftdata.o $(DMALLOCLINK) -lm

patchdata: patchdata.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ patchdata.o $(DMALLOCLINK)

dat2sdat: dat2sdat.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ dat2sdat.o $(PRESTOLINK) -lm

sdat2dat: sdat2dat.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ sdat2dat.o $(PRESTOLINK) -lm

check_parkes_raw: check_parkes_raw.o multibeam.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ check_parkes_raw.o multibeam.o $(PRESTOLINK) -lm

downsample: downsample_cmd.c downsample.o downsample_cmd.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ downsample.o downsample_cmd.o $(PRESTOLINK) -lm

split_parkes_beams: split_parkes_beams.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ split_parkes_beams.o

test_multifiles: test_multifiles.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ test_multifiles.o $(PRESTOLINK) -lm

rfifind: rfifind_cmd.c rfifind_cmd.o rfifind.o rfi_utils.o rfifind_plot.o range_parse.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@  $(INSTRUMENTOBJS) $(PLOT2DOBJS) rfifind.o rfi_utils.o rfifind_cmd.o rfifind_plot.o range_parse.o $(PRESTOLINK) $(PGPLOTLINK) -lcfitsio -lm

prepdata: prepdata_cmd.c prepdata_cmd.o prepdata.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ prepdata.o prepdata_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lcfitsio -lm

prepsubband: prepsubband_cmd.c prepsubband_cmd.o prepsubband.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ prepsubband.o prepsubband_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lcfitsio -lm

prepfold: prepfold_cmd.c prepfold_cmd.o prepfold.o prepfold_utils.o prepfold_plot.o least_squares.o polycos.o readpar.o $(INSTRUMENTOBJS) $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ prepfold.o prepfold_utils.o prepfold_plot.o prepfold_cmd.o least_squares.o polycos.o readpar.o $(PLOT2DOBJS) $(INSTRUMENTOBJS) $(LAPACKLINK) $(PRESTOLINK) $(PGPLOTLINK) -lcfitsio -lm

dump_spigot_zerolag: dump_spigot_zerolag.o spigot.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ dump_spigot_zerolag.o spigot.o $(PRESTOLINK) -lm

spigot2filterbank: spigot2filterbank_cmd.c spigot2filterbank_cmd.o spigot2filterbank.o spigot.o sigproc_fb.o sla.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ spigot2filterbank.o spigot.o sigproc_fb.o spigot2filterbank_cmd.o sla.o $(PRESTOLINK) -lsla -lm

GBT350filterbank: GBT350filterbank.o spigot.o sigproc_fb.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ GBT350filterbank.o spigot.o sigproc_fb.o $(PRESTOLINK) -lm

spigotSband2filterbank: spigotSband2filterbank.o spigot.o sigproc_fb.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ spigotSband2filterbank.o spigot.o sigproc_fb.o $(PRESTOLINK) -lm

show_pfd: show_pfd_cmd.c show_pfd.o show_pfd_cmd.o prepfold_utils.o prepfold_plot.o least_squares.o range_parse.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ show_pfd.o show_pfd_cmd.o prepfold_utils.o prepfold_plot.o least_squares.o range_parse.o $(PLOT2DOBJS) $(LAPACKLINK) $(PRESTOLINK) $(PGPLOTLINK) -lm

makedata: com.o randlib.o makedata.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ com.o randlib.o makedata.o $(PRESTOLINK) -lm

makeinf: makeinf.o ioinf.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ makeinf.o ioinf.o $(PRESTOLINK) -lm

mjd2cal: djcl.o mjd2cal.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ djcl.o mjd2cal.o $(DMALLOCLINK) -lm

cal2mjd: cldj.o cal2mjd.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ cldj.o cal2mjd.o $(DMALLOCLINK) -lm

plotbincand: plotbincand.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ $(PLOT2DOBJS) plotbincand.o $(PRESTOLINK) $(PGPLOTLINK) -lm

profile: profile_cmd.c profile_cmd.o profile.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ $(PLOT2DOBJS) profile.o profile_cmd.o $(PRESTOLINK) $(PGPLOTLINK) -lm

psrorbit: psrorbit.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ $(PLOT2DOBJS) psrorbit.o $(PRESTOLINK) $(PGPLOTLINK) -lm

testbinresp: testbinresp.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ testbinresp.o $(PLOT2DOBJS) $(PGPLOTLINK) $(PRESTOLINK) -lm

quicklook: quicklook.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ quicklook.o $(PRESTOLINK) -lm

readfile: readfile_cmd.c readfile_cmd.o readfile.o $(INSTRUMENTOBJS)
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ readfile.o readfile_cmd.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lcfitsio -lm

realfft: realfft_cmd.c realfft_cmd.o realfft.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ realfft.o realfft_cmd.o $(PRESTOLINK) -lm

rednoise: rednoise_cmd.c rednoise.o rednoise_cmd.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ rednoise.o rednoise_cmd.o $(PRESTOLINK) -lm

search_bin: search_bin_cmd.c search_bin_cmd.o search_bin.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ search_bin.o search_bin_cmd.o $(PRESTOLINK) -lm

search_rzw: search_rzw_cmd.c search_rzw_cmd.o search_rzw.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ search_rzw.o search_rzw_cmd.o $(PRESTOLINK) -lm

taperaw: taperaw.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ taperaw.o $(DMALLOCLINK)

toas2dat: toas2dat_cmd.c toas2dat_cmd.o toas2dat.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ toas2dat.o toas2dat_cmd.o $(DMALLOCLINK)

un_sc_td:
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ un_sc_td.c $(DMALLOCLINK)

swap_endian:
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ swap_endian.c $(DMALLOCLINK)

window: window.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ $(PLOT2DOBJS) window.o $(PRESTOLINK) $(PGPLOTLINK) -lm

zapbirds: zapbirds_cmd.c zapbirds_cmd.o zapbirds.o zapping.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ zapbirds_cmd.o zapbirds.o zapping.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) $(OTHERLIB) -lglib-2.0 -lm

explorefft: explorefft.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ explorefft.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

exploredat: exploredat.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ exploredat.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

weight_psrfits: weight_psrfits.o psrfits.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ weight_psrfits.o psrfits.o $(PRESTOLINK) -lcfitsio -lm

clean:
	rm -f *.o *~ *#
	rm -f slalib/*.o slalib/sla_test

cleaner: clean
	cd ../bin ; rm -f $(BINARIES)
	rm -f $(PRESTO)/lib/libpresto.* $(PRESTO)/lib/libsla.*

squeaky:  cleaner
	rm -f *.dat *.fft *.inf fftw_wisdom.txt
	rm -f core *.win* *.ps *_rzw *.tmp
	cd $(PRESTO)/clig ; rm -f *# *~
	cd $(PRESTO)/docs ; rm -f *# *~
	cd $(PRESTO)/python ; rm -f *# *~ *.o *.pyc *.pyo
	cd $(PRESTO)/include ; rm -f *# *~
