#  Makefile for PRESTO
#   by Scott M. Ransom
#		GPU-related sections by Jintao Luo

# OS type
OS = Linux
#OS = OSX

# Linux is the first choice
ifeq ($(OS),Linux)
	LIBSUFFIX = .so
	LIBCMD = -shared	
# else assume Darwin (i.e. OSX)
else
	LIBSUFFIX = .dylib
	LIBCMD = -dynamiclib
endif


#------------------------------------------------------------------------------------------
#determine if compile with CUDA
use_cuda = no# yes/no

ifeq ($(use_cuda),yes)
$(info Compile with CUDA)
#******************************************************************
#                     for GUP version, CUDA 
# CUDA compile option
# CUDA code generation flags
	GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
	GENCODE_SM10    := -gencode arch=compute_11,code=sm_11
	GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
	GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35
	GENCODE_FLAGS   := $(GENCODE_SM10) $(GENCODE_SM20) $(GENCODE_SM30)
# Location of the CUDA Toolkit binaries and libraries
#modify it if necessary
#CUDA_PATH       ?= /usr/local/cuda-5.0
	CUDA_PATH       ?= /usr/local/cuda

#modify it if necessary
	CUDA_INC_PATH   ?= -I $(CUDA_PATH)/include -I. -I.. -I $(CUDA_PATH)/samples/common/inc/ -I$(PRESTO)/include

	CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
	NVCC            ?= $(CUDA_BIN_PATH)/nvcc

#modify it if necessary, for example, if 32-bit machine, -m32
	NVCCFLAGS 			?= -m64 -ccbin g++

	EXTRA_NVCCFLAGS ?= -Xcompiler -fPIC -O2

#modify it if necessary, -m64 to -m32 for 32-bit machine
	LINKCOMMAND_CUDA = g++ -m64 -O2 -o	

#modify it if necessary, maybe lib64 should be lib32?
	CUDA_LINK_EXTRA_FLAGS = -L $(CUDA_PATH)/lib64 -lcufft  -lcudart
#******************************************************************
endif

ifeq ($(use_cuda),no)
$(info Compile without CUDA)
endif
#------------------------------------------------------------------------------------------



# How to link with some needed libraries of PGPLOT
X11LINK := $(shell pkg-config --libs x11)
PNGLINK := $(shell pkg-config --libs libpng)

# Include and link information for PGPLOT v5.X (including shared libs!)
# Typically you need to have your PGPLOT_DIR environment variable set
PGPLOTINC = -I$(PGPLOT_DIR)
PGPLOTLINK = -L$(PGPLOT_DIR) -lcpgplot -lpgplot $(X11LINK) $(PNGLINK)

# Include and link information for the FFTW 3.X single-precision library
FFTINC := $(shell pkg-config --cflags fftw3f)
FFTLINK := $(shell pkg-config --libs fftw3f)

# Include and link information for the GLIB 2.0 library
GLIBINC := $(shell pkg-config --cflags glib-2.0)
GLIBLINK := $(shell pkg-config --libs glib-2.0)

# Include and link information for CFITSIO
CFITSIOINC := $(shell pkg-config --cflags cfitsio)
CFITSIOLINK := $(shell pkg-config --libs cfitsio)

# The standard PRESTO libraries to link into executables
#PRESTOLINK = $(CFITSIOLINK) -L$(PRESTO)/lib -lpresto $(FFTLINK)
ifeq ($(use_cuda),no)
PRESTOLINK = $(CFITSIOLINK) -L$(PRESTO)/lib -lpresto $(FFTLINK)
endif

ifeq ($(use_cuda),yes)
PRESTOLINK = $(CFITSIOLINK) -L$(PRESTO)/lib -lpresto $(FFTLINK) $(CUDA_LINK_EXTRA_FLAGS) -L $(SYSDIR)/lib/x86_64-linux-gnu -lstdc++
endif

CC = gcc
FC = gfortran
#FC = g77
# Very recent Intel CPUs might see a few percent speedup using -mavx
CFLAGS = -I$(PRESTO)/include $(GLIBINC) $(CFITSIOINC) $(PGPLOTINC) $(FFTINC) \
	-DUSEFFTW -DUSEMMAP -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
	-g -O3 -ffast-math -Wall -W -fPIC

ifeq ($(use_cuda),yes)
	CFLAGS += -DUSECUDA
endif
	
#	-g -Wall -W -fPIC
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

# When modifying the CLIG files, the is the location of the clig binary
CLIG = clig

# Rules for CLIG generated files
%_cmd.c : ../clig/%_cmd.cli
	cd ../clig ; $(CLIG) -o $*_cmd -d $<
	mv ../clig/$*_cmd.h ../include/
	mv ../clig/$*_cmd.c .
	cp ../clig/$*.1 ../docs/

ifeq ($(use_cuda),yes)
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
	fitsfile.o hget.o hput.o imio.o djcl.o accel_utils_gpu.o
endif
ifeq ($(use_cuda),no)
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
	fitsfile.o hget.o hput.o imio.o djcl.o 
endif

INSTRUMENTOBJS = backend_common.o zerodm.o sigproc_fb.o psrfits.o

# Use old header reading stuff for readfile
READFILEOBJS = $(INSTRUMENTOBJS) multibeam.o bpp.o spigot.o \
	wapp.o wapp_head_parse.o wapp_y.tab.o 

PLOT2DOBJS = powerplot.o xyline.o

BINARIES = makedata makeinf mjd2cal realfft quicklook\
	search_bin search_rzw swap_endian prepdata\
	check_parkes_raw bary shiftdata dftfold\
	patchdata readfile toas2dat taperaw\
	accelsearch prepsubband cal2mjd split_parkes_beams\
	dat2sdat sdat2dat downsample rednoise un_sc_td bincand\
	psrorbit window plotbincand prepfold show_pfd\
	rfifind zapbirds explorefft exploredat\
	weight_psrfits fitsdelrow fitsdelcol psrfits_dumparrays
#	dump_spigot_zerolag spigot2filterbank\
#	spigotSband2filterbank GBT350filterbank\

all: libpresto slalib binaries

# Default indentation is K&R style with no-tabs,
# an indentation level of 4 (default), and a line-length of 85
indent:
	indent -kr -nut -l85 *.c
	rm *.c~

prep:
	touch *_cmd.c

ifeq ($(use_cuda),yes)
accel_utils_gpu.o:	accel_utils_gpu.cu	  
	$(NVCC) $(EXTRA_NVCCFLAGS) $(NVCCFLAGS) $(CUDA_INC_PATH)  $(GENCODE_FLAGS) -o $@ -c $<
endif
ifeq ($(use_cuda),no)
$(info accel_utils_gpu.o will not be compiled)
endif

makewisdom:
	$(CC) $(CLINKFLAGS) -o $@ makewisdom.c $(FFTLINK)
	./makewisdom
	cp fftw_wisdom.txt $(PRESTO)/lib

timetest: 
	$(CC) -o $@ timetest.c
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

ifeq ($(use_cuda),yes)
$(info link accelsearch with CUDA)
accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o accel_utils_gpu.o 
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o accel_utils_gpu.o $(PRESTOLINK) -lglib-2.0 -lm $(CUDA_LINK_EXTRA_FLAGS)
endif
ifeq ($(use_cuda),no)
$(info link accelsearch without CUDA)
accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o $(PRESTOLINK) $(GLIBLINK) -lm
endif

mpi: mpiprepsubband

mpiprepsubband_utils.o: mpiprepsubband_utils.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband_utils.c

mpiprepsubband.o: mpiprepsubband.c
	mpicc $(CLINKFLAGS) -c mpiprepsubband.c

mpiprepsubband: mpiprepsubband_cmd.c mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS)
	mpicc $(CLINKFLAGS) -o $(PRESTO)/bin/$@ mpiprepsubband_cmd.o mpiprepsubband_utils.o mpiprepsubband.o $(INSTRUMENTOBJS) $(PRESTOLINK) -lcfitsio -lm

#accelsearch: accelsearch_cmd.c accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o
#	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ accelsearch_cmd.o accel_utils.o accelsearch.o zapping.o $(PRESTOLINK) $(GLIBLINK) -lm

bary: bary.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ bary.o $(PRESTOLINK) -lm

bincand: bincand_cmd.c bincand_cmd.o bincand.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ bincand.o bincand_cmd.o $(PRESTOLINK) -lm

dftfold: dftfold_cmd.c dftfold_cmd.o dftfold.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ dftfold.o dftfold_cmd.o $(PRESTOLINK) -lm

shiftdata: shiftdata.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ shiftdata.o -lm

patchdata: patchdata.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ patchdata.o

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
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ djcl.o mjd2cal.o -lm

cal2mjd: cldj.o cal2mjd.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ cldj.o cal2mjd.o -lm

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

readfile: readfile_cmd.c readfile_cmd.o readfile.o $(READFILEOBJS) 
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ readfile.o readfile_cmd.o $(READFILEOBJS) $(PRESTOLINK) -lcfitsio -lm

realfft: realfft_cmd.c realfft_cmd.o realfft.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ realfft.o realfft_cmd.o $(PRESTOLINK) -lm

rednoise: rednoise_cmd.c rednoise.o rednoise_cmd.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ rednoise.o rednoise_cmd.o $(PRESTOLINK) -lm

search_bin: search_bin_cmd.c search_bin_cmd.o search_bin.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ search_bin.o search_bin_cmd.o $(PRESTOLINK) -lm

search_rzw: search_rzw_cmd.c search_rzw_cmd.o search_rzw.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ search_rzw.o search_rzw_cmd.o $(PRESTOLINK) -lm

taperaw: taperaw.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ taperaw.o

toas2dat: toas2dat_cmd.c toas2dat_cmd.o toas2dat.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ toas2dat.o toas2dat_cmd.o

un_sc_td:
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ un_sc_td.c

swap_endian:
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ swap_endian.c

window: window.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ $(PLOT2DOBJS) window.o $(PRESTOLINK) $(PGPLOTLINK) -lm

zapbirds: zapbirds_cmd.c zapbirds_cmd.o zapbirds.o zapping.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ zapbirds_cmd.o zapbirds.o zapping.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) $(GLIBLINK) -lm

explorefft: explorefft.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ explorefft.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

exploredat: exploredat.o $(PLOT2DOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ exploredat.o $(PLOT2DOBJS) $(PRESTOLINK) $(PGPLOTLINK) -lm

weight_psrfits: weight_psrfits.o $(INSTRUMENTOBJS)
	$(FC) $(FLINKFLAGS) -o $(PRESTO)/bin/$@ weight_psrfits.o $(INSTRUMENTOBJS) $(PRESTOLINK)

psrfits_dumparrays: psrfits_dumparrays.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ psrfits_dumparrays.o $(CFITSIOLINK)

fitsdelrow: fitsdelrow.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ fitsdelrow.o $(CFITSIOLINK)

fitsdelcol: fitsdelcol.o
	$(CC) $(CLINKFLAGS) -o $(PRESTO)/bin/$@ fitsdelcol.o $(CFITSIOLINK)

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

