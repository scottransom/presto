executable('cal2mjd', 'cal2mjd.c',
    dependencies: [erfa, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('fitsdelcol', 'fitsdelcol.c',
    dependencies: [fits, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('fitsdelrow', 'fitsdelrow.c',
    dependencies: [fits, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('makewisdom', 'makewisdom.c', 'datadir.c',
    dependencies: [fftw, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('mjd2cal', 'mjd2cal.c',
    dependencies: [erfa, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('patchdata', 'patchdata.c', install: false)

executable('psrfits_dumparrays', 'psrfits_dumparrays.c',
    dependencies: [fits, libm], include_directories: inc, install: true, install_rpath: presto_rpath)

executable('shiftdata', 'shiftdata.c', install: false)

executable('split_parkes_beams', 'split_parkes_beams.c', install: false)

executable('swap_endian', 'swap_endian.c', install: false)

executable('taperaw', 'taperaw.c', install: false)

executable('toas2dat', 'toas2dat.c', 'toas2dat_cmd.c',
    dependencies: libm, include_directories: inc, install: true, install_rpath: presto_rpath)

executable('un_sc_td', 'un_sc_td.c', install: false)

libpresto = library(
    'presto', 'atwood.c', 'barycenter.c', 'birdzap.c',
    'characteristics.c', 'chkio.c',
    'clipping.c', 'corr_prep.c', 'corr_routines.c', 'correlations.c',
    'database.c', 'datadir.c', 'dispersion.c', 'fastffts.c',
    'fftcalls.c', 'fitsfile.c', 'fold.c', 'fresnl.c',
    'get_candidates.c', 'hget.c', 'hput.c', 'imio.c', 'ioinf.c',
    'iomak.c', 'mask.c', 'maximize_r.c', 'maximize_rz.c',
    'maximize_rzw.c', 'median.c', 'minifft.c', 'misc_utils.c',
    'observatories.c', 'orbint.c',
    'output.c', 'range_parse.c', 'read_fft.c', 'readpar.c', 'responses.c',
    'rzinterp.c', 'rzwinterp.c', 'select.c', 'sorter.c', 'swapendian.c',
    'transpose.c', 'twopass.c', 'twopass_real_fwd.c',
    'twopass_real_inv.c', 'vectors.c',
    dependencies: [glib, fftw, gsl, erfa, libm],
    include_directories: inc,
    install: true,
    install_rpath: presto_rpath
)

INSTRUMENTOBJS= ['backend_common.c', 'psrfits.c', 'sigproc_fb.c', 'zerodm.c']
PLOT2DOBJS = ['powerplot.c', 'xyline.c']

executable('accelsearch', 'accelsearch.c', 'accelsearch_cmd.c', 'accel_utils.c', 'zapping.c',
    dependencies: [glib, fftw, gsl, libm, omp], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('bary', 'bary.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('bincand', 'bincand.c', 'bincand_cmd.c',
    dependencies: [glib, fftw, gsl, libm, omp],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('check_bary_erfa', 'check_bary_erfa.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: false)

executable('check_polycos_t2', 'check_polycos_t2.c', 'polycos.c',
    dependencies: [glib, fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: false)

executable('check_parkes_raw', 'check_parkes_raw.c', 'multibeam.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('dat2sdat', 'dat2sdat.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('downsample', 'downsample.c', 'downsample_cmd.c', 
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('exploredat',
    sources: ['exploredat.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('explorefft',
    sources: ['explorefft.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png], c_args: '-DUSEMMAP',
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('makedata', 'makedata.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('makeinf', 'makeinf.c', 'ioinf.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

if mpi.found()
    executable('mpiprepsubband',
        sources: ['mpiprepsubband.c', 'mpiprepsubband_cmd.c', 'mpiprepsubband_utils.c'] + INSTRUMENTOBJS,
        dependencies: [glib, fftw, gsl, libm, fits, mpi],
        include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)
endif

executable('plotbincand',
    sources: ['plotbincand.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('prepdata',
    sources: ['prepdata.c', 'prepdata_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, gsl, libm, fits],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('prepfold',
    sources: ['prepfold.c', 'prepfold_pipeline.c', 'prepfold_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c', 'polycos.c'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, fits, pgplot, cpgplot, x11, png, omp],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('prepsubband',
    sources: ['prepsubband.c', 'prepsubband_cmd.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, gsl, libm, fits],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('prepfold_multi',
    sources: ['prepfold_multi.c', 'prepfold_multi_cmd.c', 'prepfold_pipeline.c', 'prepfold_utils.c', 'prepfold_plot.c', 'polycos.c'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, fits, pgplot, cpgplot, x11, png, omp],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

# Stage 3 unit check: factored reader vs original read_subbands (not installed)
executable('reader_equiv_check',
    sources: [files('../tests/prepfold/reader_equiv_check.c')] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, gsl, libm, fits],
    include_directories: inc, link_with: libpresto,
    build_by_default: false, install: false)

# coherentsearch (C port of CoherentSearch.jl; see ROADMAP.md).  Linked into the
# executables that use it rather than into libpresto, for now.
CSOBJS = ['cs_interp.c', 'cs_metric.c', 'cs_chunk.c', 'cs_search.c']

executable('coherentsearch', ['coherentsearch.c', 'coherentsearch_cmd.c', 'cs_plot.c'] + CSOBJS,
    dependencies: [glib, fftw, gsl, libm, omp, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

# Dumps one chunk's stages for tests/test_coherentsearch.py (not installed).
executable('cs_dump',
    sources: [files('../tests/coherentsearch/cs_dump.c')] + CSOBJS,
    dependencies: [glib, fftw, gsl, libm, omp],
    include_directories: inc, link_with: libpresto,
    install: false)

executable('psrorbit',
    sources: ['psrorbit.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('quicklook', 'quicklook.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('readfile',
    sources: ['readfile.c', 'readfile_cmd.c', 'multibeam.c', 'bpp.c', 'spigot.c', 'wapp.c', 'wapp_head_parse.c', 'wapp_y.tab.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, gsl, libm, fits],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('realfft', 'realfft.c', 'realfft_cmd.c',
    dependencies: [fftw, gsl, libm],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('rednoise', 'rednoise.c', 'rednoise_cmd.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('rfifind',
    sources: ['rfifind.c', 'rfifind_cmd.c', 'rfi_utils.c', 'rfifind_plot.c'] + INSTRUMENTOBJS + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, fits, pgplot, cpgplot, x11, png, omp],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('sdat2dat', 'sdat2dat.c',
    dependencies: [fftw, gsl, libm], include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('search_bin', 'search_bin.c', 'search_bin_cmd.c',
    dependencies: [glib, fftw, gsl, libm],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('show_pfd',
    sources: ['show_pfd.c', 'show_pfd_cmd.c', 'prepfold_utils.c', 'prepfold_plot.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('weight_psrfits',
    sources: ['weight_psrfits.c'] + INSTRUMENTOBJS,
    dependencies: [glib, fftw, gsl, libm, fits, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('window',
    sources: ['window.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)

executable('zapbirds',
    sources: ['zapbirds.c', 'zapbirds_cmd.c', 'zapping.c'] + PLOT2DOBJS,
    dependencies: [glib, fftw, gsl, libm, pgplot, cpgplot, x11, png],
    include_directories: inc, link_with: libpresto, install: true, install_rpath: presto_rpath)
