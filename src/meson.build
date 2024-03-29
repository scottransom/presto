libpresto = library(
    'presto', 'amoeba.c', 'atwood.c', 'barycenter.c', 'birdzap.c',
    'cand_output.c', 'characteristics.c', 'chkio.c', 'cldj.c',
    'clipping.c', 'corr_prep.c', 'corr_routines.c', 'correlations.c',
    'database.c', 'dcdflib.c', 'dispersion.c', 'djcl.c', 'fastffts.c',
    'fftcalls.c', 'fitsfile.c', 'fminbr.c', 'fold.c', 'fresnl.c',
    'get_candidates.c', 'hget.c', 'hput.c', 'imio.c', 'ioinf.c',
    'iomak.c', 'ipmpar.c', 'mask.c', 'maximize_r.c', 'maximize_rz.c',
    'maximize_rzw.c', 'median.c', 'minifft.c', 'misc_utils.c', 'orbint.c',
    'output.c', 'range_parse.c', 'read_fft.c', 'readpar.c', 'responses.c',
    'rzinterp.c', 'rzwinterp.c', 'select.c', 'sorter.c', 'swapendian.c',
    'transpose.c', 'twopass.c', 'twopass_real_fwd.c',
    'twopass_real_inv.c', 'vectors.c',
    dependencies: [glib, fftw, libm, omp],
    include_directories: inc,
    install: true,
)

INSTRUMENTOBJS= ['backend_common.c', 'zerodm.c', 'sigproc_fb.c', 'psrfits.c']
PLOT2DOBJS = ['powerplot.c', 'xyline.c']

accelsearch = executable('accelsearch', 'accelsearch.c', 'accelsearch_cmd.c', 'accel_utils.c', 'zapping.c',
    dependencies: [glib, fftw, libm, omp], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

bary = executable('bary', 'bary.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

bincand = executable('bincand', 'bincand.c', 'bincand_cmd.c',
    dependencies: [glib, fftw, libm, omp],
    include_directories: inc, link_with: libpresto, install: true)

check_parkes_raw = executable('check_parkes_raw', 'check_parkes_raw.c', 'multibeam.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

dat2sdat = executable('dat2sdat', 'dat2sdat.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

dftfold = executable('dftfold', 'dftfold_cmd.c', 'dftfold.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

downsample = executable('downsample', 'downsample.c', 'downsample_cmd.c', 
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

exploredat = executable('exploredat',
    sources: ['exploredat.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

explorefft = executable('explorefft',
    sources: ['explorefft.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

makedata = executable('makedata', 'makedata.c', 'com.c', 'randlib.c', 
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

makeinf = executable('makeinf', 'makeinf.c', 'ioinf.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

plotbincand = executable('plotbincand',
    sources: ['plotbincand.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

prepdata = executable('prepdata',
    sources: ['prepdata.c', 'prepdata_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits],
    include_directories: inc, link_with: libpresto, install: true)

prepfold = executable('prepfold',
    sources: ['prepfold.c', 'prepfold_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c', 'polycos.c', 'least_squares.f'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

prepsubband = executable('prepsubband',
    sources: ['prepsubband.c', 'prepsubband_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits],
    include_directories: inc, link_with: libpresto, install: true)

psrorbit = executable('psrorbit',
    sources: ['psrorbit.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

quicklook = executable('quicklook', 'quicklook.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

readfile = executable('readfile',
    sources: ['readfile.c', 'readfile_cmd.c', 'multibeam.c', 'bpp.c', 'spigot.c', 'wapp.c', 'wapp_head_parse.c', 'wapp_y.tab.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits],
    include_directories: inc, link_with: libpresto, install: true)

realfft = executable('realfft', 'realfft.c', 'realfft_cmd.c',
    dependencies: [fftw, libm],
    include_directories: inc, link_with: libpresto, install: true)

rednoise = executable('rednoise', 'rednoise.c', 'rednoise_cmd.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

rfifind = executable('rfifind',
    sources: ['rfifind.c', 'rfifind_cmd.c', 'rfi_utils.c', 'rfifind_plot.c'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

sdat2dat = executable('sdat2dat', 'sdat2dat.c',
    dependencies: libm, include_directories: inc, link_with: libpresto, install: true)

search_bin = executable('search_bin', 'search_bin.c', 'search_bin_cmd.c',
    dependencies: [glib, fftw, libm],
    include_directories: inc, link_with: libpresto, install: true)

show_pfd = executable('show_pfd',
    sources: ['show_pfd.c', 'show_pfd_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c', 'least_squares.f'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

weight_psrfits = executable('weight_psrfits',
    sources: ['weight_psrfits.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

window = executable('window',
    sources: ['window.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

zapbirds = executable('zapbirds',
    sources: ['zapbirds.c', 'zapbirds_cmd.c', 'zapping.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)
