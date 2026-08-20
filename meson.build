project('PRESTO', 'c',
  version: '6.0.1',
  license: 'GPL-2.0',
  meson_version: '>= 1.3.0',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
# Other buildtype possibilities:
#    'buildtype=release',        # -O3 with no -g.  Measurably no faster than
#                                #   the default, so we keep the debug symbols.
#    'buildtype=debugoptimized', # -O2 plus -g  (the default, above)
#    'buildtype=debug',
# The optimization level is not the interesting knob:  what actually speeds up
# the compute-bound tools is letting the compiler use the vector instructions
# of the machine you are on.  Adding '-Dc_args=-march=native' to 'meson setup'
# makes accelsearch roughly 1.5-2x faster on a modern x86-64 CPU (AVX2+FMA).
# Only do that for builds you will run on that same machine, though -- see the
# notes in INSTALL.md.
# Can also set install prefix as option via e.g. 'prefix=/home/sransom'
#    'prefix=/home/sransom'

cc = meson.get_compiler('c')

glib = dependency('glib-2.0')
fftw = dependency('fftw3f') # Note the 'f' for single-precision version!
fits = dependency('cfitsio')
gsl = dependency('gsl')
# If ERFA is not installed on the system, meson automatically downloads
# and builds it via subprojects/erfa.wrap (see INSTALL.md)
erfa = dependency('erfa')
x11 = dependency('x11')
png = dependency('libpng')
omp = dependency('openmp', required: false)
# omp = dependency('NO_openmp', required: false)
mpi = dependency('mpi', language: 'c', required: false,
  not_found_message: 'MPI not found. Skipping mpiprepsubband build.')

libm = cc.find_library('m', required: false)
# PGPLOT: prefer pkg-config (e.g. conda-forge's pgplot ships .pc files);
# fall back to a bare library search for distros without them (e.g.
# Ubuntu's pgplot5)
pgplot = dependency('pgplot', required: false)
if not pgplot.found()
  pgplot = cc.find_library('pgplot', required: true)
endif
cpgplot = dependency('cpgplot', required: false)
if not cpgplot.found()
  cpgplot = cc.find_library('cpgplot', required: true)
endif
# Make sure cpgplot.h is actually reachable: it is in the default include
# path on e.g. Ubuntu, but in <includedir>/pgplot for conda-forge's pgplot
# (whose .pc file does not say so), and in $PGPLOT_DIR for classic source
# installs
if not cc.has_header('cpgplot.h', dependencies: cpgplot)
  fs = import('fs')
  pg_hdr_cands = []
  if cpgplot.type_name() == 'pkgconfig'
    pg_hdr_cands += cpgplot.get_variable(pkgconfig: 'includedir') / 'pgplot'
  endif
  pg_env = run_command('sh', '-c', 'printf %s "$PGPLOT_DIR"', check: false).stdout()
  if pg_env != ''
    pg_hdr_cands += pg_env
  endif
  pg_hdr_found = false
  foreach dir : pg_hdr_cands
    # A raw -I (not include_directories) since e.g. a pixi env lives
    # inside the source tree, which include_directories() refuses
    if not pg_hdr_found and fs.exists(dir / 'cpgplot.h')
      cpgplot = declare_dependency(dependencies: cpgplot,
                                   compile_args: '-I' + dir)
      pg_hdr_found = true
    endif
  endforeach
  if not pg_hdr_found
    error('Cannot find cpgplot.h -- install PGPLOT dev files or set PGPLOT_DIR.')
  endif
endif

inc = include_directories('include')

_global_c_args = cc.get_supported_arguments(
  '-DUSE_FFTW_MALLOC',
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
  '-Wno-unused-result',
  '-Wno-unused-but-set-parameter'
)
add_project_arguments(_global_c_args, language : 'c')

# Installed runtime data (pulsar catalog, optional FFTW wisdom) lives here.
# The absolute path is compiled in so the tools locate it without the $PRESTO
# environment variable being set (see src/datadir.c).  $PRESTO/lib still works
# as an override for running from a source tree.
presto_datadir = get_option('prefix') / get_option('datadir') / 'presto'
add_project_arguments('-DPRESTO_DATADIR="' + presto_datadir + '"', language : 'c')

# Bake an RPATH to the install libdir into libpresto and every installed
# executable, so at runtime they load the FFTW (and other deps) that ship in
# the install prefix rather than falling through the loader cache to a
# possibly-different system copy.  This matters in conda/pixi environments,
# where <prefix>/lib holds a libfftw3f of a different version than the system
# one: makewisdom and the tools would otherwise write/read FFTW wisdom with
# mismatched libraries, producing spurious "wisdom not up-to-date" warnings
# (the wisdom version cookie is per-FFTW-build).  The absolute <prefix>/libdir
# path is exactly right for a dev/pixi tree and, on conda-forge, is rewritten
# to a relocatable $ORIGIN form by conda-build's binary-relocation step.  This
# mirrors the absolute PRESTO_DATADIR baked in just above.
presto_rpath = get_option('prefix') / get_option('libdir')

# These next couple things are from the Scipy top-level meson.build
if host_machine.system() == 'darwin'
  if cc.has_link_argument('-Wl,-ld_classic')
    # New linker introduced in macOS 14 not working yet, see gh-19357 and gh-19387
    add_project_link_arguments('-Wl,-ld_classic', language : ['c', 'cpp'])
  endif
  # Work around duplicate RPATH issue on macOS with gcc15
  # See: https://github.com/mesonbuild/meson/issues/12314
  if cc.has_link_argument('-Wl,-no_warn_duplicate_libraries')
    add_project_link_arguments('-Wl,-no_warn_duplicate_libraries', language : ['c'])
  endif
endif

subdir('src')
subdir('bin')

# Runtime data files, installed to <prefix>/share/presto so they are found
# without $PRESTO.  The pulsar catalog is the binary form read by libpresto's
# database.c; the text sources it is built from are installed alongside it for
# the Python pypsrcat module.
install_data(
  'lib/pulsars.cat', 'lib/psr_catalog.txt', 'lib/aliases.txt',
  install_dir: presto_datadir,
)

# Man pages (generated from the clig/*.cli specs) install to the standard
# <prefix>/share/man/man1, so `man rfifind` etc. work with no $PRESTO or
# extra MANPATH setup (conda activation already adds <prefix>/share/man).
install_man(
  'docs/accelsearch.1', 'docs/bincand.1', 'docs/downsample.1',
  'docs/mpiprepsubband.1', 'docs/prepdata.1', 'docs/prepfold.1',
  'docs/prepfold_multi.1', 'docs/prepsubband.1', 'docs/readfile.1',
  'docs/realfft.1', 'docs/rednoise.1', 'docs/rfifind.1',
  'docs/search_bin.1', 'docs/show_pfd.1', 'docs/toas2dat.1',
  'docs/zapbirds.1',
)

import('pkgconfig').generate(
  libpresto,
  filebase: 'presto',
  name: 'PRESTO',
  description: 'PulsaR Exploration and Search TOolkit by Scott Ransom',
)

# Post-install script to fix duplicate RPATHs on macOS
if host_machine.system() == 'darwin'
  meson.add_install_script('fix_rpath_duplicates.py', get_option('prefix'))
endif
