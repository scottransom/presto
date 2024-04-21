executable('cal2mjd', 'cal2mjd.c', 'cldj.c',
    dependencies: libm, include_directories: inc, install: true)

executable('fitsdelcol', 'fitsdelcol.c',
    dependencies: [fits, libm], include_directories: inc, install: true)

executable('fitsdelrow', 'fitsdelrow.c',
    dependencies: [fits, libm], include_directories: inc, install: true)

executable('makewisdom', 'makewisdom.c',
    dependencies: [fftw, libm], include_directories: inc, install: false)

executable('mjd2cal', 'mjd2cal.c', 'djcl.c',
    dependencies: libm, include_directories: inc, install: true)

executable('patchdata', 'patchdata.c', install: false)

executable('psrfits_dumparrays', 'psrfits_dumparrays.c',
    dependencies: [fits, libm], include_directories: inc, install: true)

executable('shiftdata', 'shiftdata.c', install: false)

executable('split_parkes_beams', 'split_parkes_beams.c', install: false)

executable('swap_endian', 'swap_endian.c', install: false)

executable('taperaw', 'taperaw.c', install: false)

executable('toas2dat', 'toas2dat.c', 'toas2dat_cmd.c',
    dependencies: libm, include_directories: inc, install: true)

executable('un_sc_td', 'un_sc_td.c', install: false)

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
    install: true
)

INSTRUMENTOBJS= ['backend_common.c', 'psrfits.c', 'sigproc_fb.c', 'zerodm.c']
PLOT2DOBJS = ['powerplot.c', 'xyline.c']

executable('accelsearch', 'accelsearch.c', 'accelsearch_cmd.c', 'accel_utils.c', 'zapping.c',
    dependencies: [glib, fftw, libm, omp], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

executable('bary', 'bary.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('bincand', 'bincand.c', 'bincand_cmd.c',
    dependencies: [glib, fftw, libm, omp],
    include_directories: inc, link_with: libpresto, install: true)

executable('check_parkes_raw', 'check_parkes_raw.c', 'multibeam.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('dat2sdat', 'dat2sdat.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('dftfold', 'dftfold_cmd.c', 'dftfold.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('downsample', 'downsample.c', 'downsample_cmd.c', 
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('exploredat',
    sources: ['exploredat.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

executable('explorefft',
    sources: ['explorefft.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true)

executable('makedata', 'makedata.c', 'com.c', 'randlib.c', 
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('makeinf', 'makeinf.c', 'ioinf.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

if mpi.found()
    executable('mpiprepsubband',
        sources: ['mpiprepsubband.c', 'mpiprepsubband_cmd.c', 'mpiprepsubband_utils.c'] + INSTRUMENTOBJS,
        dependencies: [glib, fftw, libm, fits, omp, mpi],
        include_directories: inc, link_with: libpresto, install: true)
endif

executable('plotbincand',
    sources: ['plotbincand.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('prepdata',
    sources: ['prepdata.c', 'prepdata_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits, omp],
    include_directories: inc, link_with: libpresto, install: true)

executable('prepfold',
    sources: ['prepfold.c', 'prepfold_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c', 'polycos.c', 'least_squares.f'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('prepsubband',
    sources: ['prepsubband.c', 'prepsubband_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits, omp],
    include_directories: inc, link_with: libpresto, install: true)

executable('psrorbit',
    sources: ['psrorbit.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('quicklook', 'quicklook.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('readfile',
    sources: ['readfile.c', 'readfile_cmd.c', 'multibeam.c', 'bpp.c', 'spigot.c', 'wapp.c', 'wapp_head_parse.c', 'wapp_y.tab.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits],
    include_directories: inc, link_with: libpresto, install: true)

executable('realfft', 'realfft.c', 'realfft_cmd.c',
    dependencies: [fftw, libm],
    include_directories: inc, link_with: libpresto, install: true)

executable('rednoise', 'rednoise.c', 'rednoise_cmd.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('rfifind',
    sources: ['rfifind.c', 'rfifind_cmd.c', 'rfi_utils.c', 'rfifind_plot.c'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('sdat2dat', 'sdat2dat.c',
    dependencies: [fftw, libm], include_directories: inc, link_with: libpresto, install: true)

executable('search_bin', 'search_bin.c', 'search_bin_cmd.c',
    dependencies: [glib, fftw, libm],
    include_directories: inc, link_with: libpresto, install: true)

executable('show_pfd',
    sources: ['show_pfd.c', 'show_pfd_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c', 'least_squares.f'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('weight_psrfits',
    sources: ['weight_psrfits.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('window',
    sources: ['window.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)

executable('zapbirds',
    sources: ['zapbirds.c', 'zapbirds_cmd.c', 'zapping.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true)
